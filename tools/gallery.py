#!/usr/bin/env python3
"""
gallery.py -- turn the hyperchoreography catalogue into ONE self-contained HTML page.

Re-runnable:  python3 tools/gallery.py                       # every catalog/*.bin
              python3 tools/gallery.py catalog/d7_n10.bin --out d7.html --samples 120

Stdlib only.  No numpy, no CDN, no external assets: everything (CSS, JS, geometry)
is inlined into a single .html file.

What it does per record
-----------------------
  list  -> `rigid` (that column exists nowhere else) and the real record ids
  show  -> scalar metadata + the rotating-frame generator Omega
  export-> S inertial samples of all N bodies

  then, in pure Python:
    * pooled second moment over ALL N bodies and ALL S samples  -> d x d symmetric
    * cyclic Jacobi eigendecomposition                          -> principal frame V
    * body 0's curve rotated into V, quantised to Int16, base64
    * G = V^T exp(-2 pi Omega / N) V   so the browser can rebuild the other N-1
      bodies exactly:  body_j(f) = G^j . body_0((f + j S/N) mod S)

Storing one curve + G instead of N curves is what keeps the page ~1.3 MB
instead of ~11 MB.
"""

import argparse
import base64
import glob
import json
import math
import os
import struct
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

# --------------------------------------------------------------------------
# linear algebra (hand written -- numpy is not available)
# --------------------------------------------------------------------------

def jacobi_eig(a, n, max_sweeps=100):
    """Cyclic Jacobi on a symmetric n x n matrix.

    Returns (evals, V) with evals descending and V's COLUMNS the eigenvectors,
    i.e. A = V diag(evals) V^T and  p = V^T q  is q in the principal frame.

    Deterministic: fixed p<q row-major sweep order, ties broken by index.
    """
    A = [row[:] for row in a]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = 0.0
        for p in range(n - 1):
            Ap = A[p]
            for q in range(p + 1, n):
                off += Ap[q] * Ap[q]
        if off < 1e-30:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = A[p][q]
                if abs(apq) < 1e-300:
                    continue
                theta = (A[q][q] - A[p][p]) / (2.0 * apq)
                if theta >= 0.0:
                    t = 1.0 / (theta + math.sqrt(theta * theta + 1.0))
                else:
                    t = -1.0 / (-theta + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):                      # columns
                    akp = A[k][p]; akq = A[k][q]
                    A[k][p] = c * akp - s * akq
                    A[k][q] = s * akp + c * akq
                for k in range(n):                      # rows
                    apk = A[p][k]; aqk = A[q][k]
                    A[p][k] = c * apk - s * aqk
                    A[q][k] = s * apk + c * aqk
                for k in range(n):                      # accumulate V
                    vkp = V[k][p]; vkq = V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    ev = [A[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: (-ev[i], i))
    ev = [ev[i] for i in order]
    V = [[V[r][c] for c in order] for r in range(n)]
    # sign fix: largest-magnitude component of each eigenvector is positive
    for c in range(n):
        best, bi = -1.0, 0
        for r in range(n):
            m = abs(V[r][c])
            if m > best + 1e-15:
                best, bi = m, r
        if V[bi][c] < 0.0:
            for r in range(n):
                V[r][c] = -V[r][c]
    # handedness fix: keep det(V) = +1 so re-runs never mirror-flip
    if det(V, n) < 0.0:
        for r in range(n):
            V[r][n - 1] = -V[r][n - 1]
    return ev, V


def det(a, n):
    """Determinant by Gaussian elimination with partial pivoting."""
    M = [row[:] for row in a]
    d = 1.0
    for c in range(n):
        piv, pv = c, abs(M[c][c])
        for r in range(c + 1, n):
            if abs(M[r][c]) > pv:
                piv, pv = r, abs(M[r][c])
        if pv < 1e-300:
            return 0.0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        inv = 1.0 / M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] * inv
            if f:
                Mr, Mc = M[r], M[c]
                for k in range(c, n):
                    Mr[k] -= f * Mc[k]
    return d


def matmul(A, B, n):
    C = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Ai, Ci = A[i], C[i]
        for k in range(n):
            aik = Ai[k]
            if aik == 0.0:
                continue
            Bk = B[k]
            for j in range(n):
                Ci[j] += aik * Bk[j]
    return C


def matvec(A, v, n):
    return [sum(A[i][k] * v[k] for k in range(n)) for i in range(n)]


def expm(A, n, terms=25):
    """exp(A) by scaling-and-squaring.  A naive Taylor series diverges here:
    ||Omega|| reaches 13 in the catalogue."""
    nrm = 0.0
    for row in A:
        s = 0.0
        for x in row:
            s += abs(x)
        if s > nrm:
            nrm = s
    sq = 0
    while nrm > 0.5:
        nrm *= 0.5
        sq += 1
    f = 2.0 ** (-sq)
    B = [[A[i][j] * f for j in range(n)] for i in range(n)]
    R = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    T = [row[:] for row in R]
    for k in range(1, terms + 1):
        T = matmul(T, B, n)
        ik = 1.0 / k
        for i in range(n):
            Ti, Ri = T[i], R[i]
            for j in range(n):
                Ti[j] *= ik
                Ri[j] += Ti[j]
    for _ in range(sq):
        R = matmul(R, R, n)
    return R


def transpose(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]


def maxabs_diff_identity(A, n):
    m = 0.0
    for i in range(n):
        for j in range(n):
            v = abs(A[i][j] - (1.0 if i == j else 0.0))
            if v > m:
                m = v
    return m


# --------------------------------------------------------------------------
# talking to ./hyperchoreography
# --------------------------------------------------------------------------

class CliError(Exception):
    pass


def _run(binpath, args, timeout=180):
    p = subprocess.run([binpath] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, timeout=timeout)
    return p.returncode, p.stdout.decode("utf-8", "replace"), p.stderr.decode("utf-8", "replace")


def cat_list(binpath, path):
    """Parse `list`.  Rows are 15 fixed whitespace fields then `sym`, which may be
    empty or contain spaces.  A header line precedes and an `<n> records` footer
    follows; both are rejected by requiring 15 fields with an integer id."""
    rc, out, err = _run(binpath, ["list", path])
    rows = []
    for line in out.splitlines():
        p = line.split(None, 15)
        if len(p) < 15:
            continue
        try:
            rid = int(p[0])
        except ValueError:
            continue
        try:
            rigid = float(p[12])
        except ValueError:
            rigid = float("nan")
        rows.append({"id": rid, "rigid": rigid,
                     "sym": p[15].strip() if len(p) > 15 else ""})
    if not rows and rc != 0:
        raise CliError("list %s failed: %s" % (path, err.strip()[:200]))
    return rows


def fingerprint(path):
    """(size, mtime_ns) of a catalogue file, or None if it vanished.  Compared
    before and after the build so a concurrent harvest can never be committed
    to HTML silently."""
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_size, st.st_mtime_ns)


def cat_show(binpath, path, rid):
    """`show` prints one line of JSON.  On a bad id it prints `error: ...` to
    stderr and can still exit 0, so validate by parsing stdout, not by rc."""
    rc, out, err = _run(binpath, ["show", path, "--id", str(rid)])
    out = out.strip()
    if not out:
        raise CliError("show %s#%d: no output (%s)" % (path, rid, err.strip()[:160]))
    try:
        js = json.loads(out.splitlines()[0])
    except Exception as e:
        raise CliError("show %s#%d: unparseable JSON (%s)" % (path, rid, e))
    # Ids are positions in the file, so a catalogue rewritten by a concurrent
    # search renumbers them.  Never accept a record that is not the one asked
    # for: silently pairing record A's scalars with record B's geometry is the
    # one failure this tool must not have.
    try:
        got = int(js.get("id", -1))
    except (TypeError, ValueError):
        got = -1
    if got != rid:
        raise CliError("show %s#%d: catalogue returned id %s (file rewritten mid-run?)"
                       % (path, rid, js.get("id")))
    return js


def cat_export(binpath, path, rid, samples, outfile, N, d):
    """`export` writes a CSV: header `t,q0_x,...` then `samples` rows of
    1 + N*d fixed-point numbers, t = 2 pi k / samples over the HALF-OPEN [0,2pi).
    Parsed by POSITION, never by header name (the d>6 header names differ
    between HEAD and the working tree of src/main.cpp)."""
    rc, out, err = _run(binpath, ["export", path, "--id", str(rid),
                                  "--samples", str(samples), "--out", outfile])
    if not os.path.exists(outfile):
        raise CliError("export %s#%d: no file (%s)" % (path, rid, err.strip()[:160]))
    want = N * d
    rows = []
    with open(outfile, "r") as fh:
        fh.readline()                                   # header
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) != want + 1:
                raise CliError("export %s#%d: %d columns, expected %d"
                               % (path, rid, len(parts), want + 1))
            rows.append([float(x) for x in parts[1:]])
    if len(rows) != samples:
        raise CliError("export %s#%d: %d rows, expected %d" % (path, rid, len(rows), samples))
    return rows


# --------------------------------------------------------------------------
# per-record geometry
# --------------------------------------------------------------------------

MULT = __import__("operator").mul

# Relative-equilibrium threshold.  The `rigid` column is trimodal over the
# catalogue: 43 records at <= 1.6e-11 (exactly rigid), 8 more in 1.3e-6..6.0e-6,
# then a 2850x gap to the first genuine choreography at 1.7e-2.  Those 8 have
# mutual distances constant to within the Int16 quantisation floor (measured
# ribbon variation <= 2.7e-5), i.e. they are relative equilibria too, just less
# tightly converged.  So the cut goes in the gap, not at 1e-6.
RIGID_TOL = 1e-4


def pooled_moment(rows, N, d):
    """Second moment pooled over ALL N bodies and ALL samples.  Its rank (at
    1e-5 relative) is exactly the stored `deff` -- pooled PCA *is* the object
    deff measures, and it is the only correct one for rotating records where
    the N bodies trace N different (rotated) curves."""
    cols = [[] for _ in range(d)]
    for j in range(N):
        lo, hi = j * d, (j + 1) * d
        sub = [r[lo:hi] for r in rows]
        for a, colvals in enumerate(zip(*sub)):
            cols[a].extend(colvals)
    inv = 1.0 / float(len(cols[0]))
    M = [[0.0] * d for _ in range(d)]
    for a in range(d):
        ca = cols[a]
        for b in range(a, d):
            v = sum(map(MULT, ca, cols[b])) * inv
            M[a][b] = v
            M[b][a] = v
    return M


def quantise(curve, S, d):
    """Per-axis Int16 scaling.  A single global scale wastes ~4 bits on the
    smallest records; per-axis costs d floats of metadata and is exact to
    ~3e-5 of each axis' extent."""
    sc = []
    for k in range(d):
        m = 0.0
        for f in range(S):
            v = abs(curve[f][k])
            if v > m:
                m = v
        sc.append((m / 32767.0) if m > 0.0 else 1.0)
    buf = bytearray(S * d * 2)
    off = 0
    for f in range(S):
        cf = curve[f]
        for k in range(d):
            q = int(round(cf[k] / sc[k]))
            if q > 32767:
                q = 32767
            elif q < -32767:
                q = -32767
            struct.pack_into("<h", buf, off, q)
            off += 2
    return sc, base64.b64encode(bytes(buf)).decode("ascii")


def g6(x):
    """Compact float for the JSON payload."""
    if x == 0.0:
        return 0
    if x != x or x in (float("inf"), float("-inf")):
        return 0
    return float("%.6g" % x)


def path_labels(paths):
    """Display key per catalogue file: the basename when it is unique, else the
    shortest trailing path fragment that separates the collision (full path as
    a last resort).  Records are keyed by this label, so two catalogues with
    the same filename in different directories never share a record key."""
    labels = {}
    by_base = {}
    for p in paths:
        by_base.setdefault(os.path.basename(p), []).append(p)
    for base, group in by_base.items():
        if len(group) == 1:
            labels[group[0]] = base
            continue
        part = {}
        for p in group:
            part[p] = os.path.normpath(os.path.abspath(p)).split(os.sep)
        deep = max(len(v) for v in part.values())
        for p in group:
            lab = base
            for k in range(2, deep + 1):
                cand = "/".join(part[p][-k:])
                lab = cand
                if sum(1 for q in group if "/".join(part[q][-k:]) == cand) == 1:
                    break
            labels[p] = lab
    return labels


def build_record(binpath, path, rid, rigid, samples_max, tmpdir, warn, label=None):
    rec = cat_show(binpath, path, rid)
    N = int(rec["N"]); d = int(rec["d"]); deff = int(rec["deff"])
    S = N * (samples_max // N)
    if S < N:
        S = N
    tmp = os.path.join(tmpdir, "e%d.csv" % (threading_ident(),))
    try:
        rows = cat_export(binpath, path, rid, S, tmp, N, d)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    om = rec.get("omega")
    rot = bool(om) and any(abs(v) > 1e-10 for v in om)
    in_frame = rot          # was this record found in a rotating frame at all?
    W = None
    rates = []
    closed = 1
    if rot:
        W = [[float(om[i * d + j]) for j in range(d)] for i in range(d)]
        # rotation rates: sqrt of the eigenvalues of -Omega^2
        W2 = matmul(W, W, d)
        nW2 = [[-W2[i][j] for j in range(d)] for i in range(d)]
        evw, _ = jacobi_eig(nW2, d)
        seen = []
        for v in evw:
            r = math.sqrt(v) if v > 0.0 else 0.0
            r = round(r, 4)
            if r not in seen:
                seen.append(r)
        rates = sorted(seen)
        TP = 2.0 * math.pi
        closed = 1 if maxabs_diff_identity(
            expm([[W[i][j] * TP for j in range(d)] for i in range(d)], d), d) < 1e-8 else 0
        if not closed:
            # Relative periodic orbit whose inertial motion does NOT close.
            # Leaving the polygon open would misrepresent the orbit, so draw it
            # in the co-rotating frame instead: exp(-Omega t) q_j(t) = q(t+2 pi j/N),
            # one shared closed loop.  Badged "open r.p.o." on the tile.
            # (No record in today's catalogue takes this path -- every Omega has
            # integer rates -- but nothing in the format forbids it.)
            E1 = expm([[-W[i][j] * (TP / S) for j in range(d)] for i in range(d)], d)
            E = [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
            for f in range(S):
                r = rows[f]
                for j in range(N):
                    r[j * d:(j + 1) * d] = matvec(E, r[j * d:(j + 1) * d], d)
                E = matmul(E, E1, d)     # E1 is orthogonal, so powers stay stable
            rot = False                  # one curve now: draw it like an inertial record
            W = None

    # ---- pooled PCA -------------------------------------------------------
    M = pooled_moment(rows, N, d)
    ev, V = jacobi_eig(M, d)
    lam0 = ev[0] if ev[0] > 0.0 else 1.0
    rank = sum(1 for v in ev if v > 1e-5 * lam0)
    if rank != deff:
        warn("%s#%d: pooled PCA rank %d != deff %d" % (os.path.basename(path), rid, rank, deff))

    # ---- body 0 in the principal frame ------------------------------------
    curve0 = []
    for f in range(S):
        r = rows[f]
        curve0.append([sum(V[a][k] * r[a] for a in range(d)) for k in range(d)])

    # ---- G, and the reconstruction check ----------------------------------
    Gflat = None
    store_all = 0
    if rot:
        E = expm([[-W[i][j] * (2.0 * math.pi / N) for j in range(d)] for i in range(d)], d)
        G = matmul(matmul(transpose(V, d), E, d), V, d)
        # verify  body_j(f) = G^j . body_0((f + j S/N) mod S)   in the principal frame
        rms = float(rec.get("rms", 1.0)) or 1.0
        tol = 1e-6 * rms
        worst = 0.0
        R = [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
        step = S // N
        for j in range(1, N):
            R = matmul(R, G, d)
            for f in range(0, S, max(1, S // 24)):
                src = curve0[(f + j * step) % S]
                pred = matvec(R, src, d)
                r = rows[f]
                for k in range(d):
                    act = sum(V[a][k] * r[j * d + a] for a in range(d))
                    e = abs(pred[k] - act)
                    if e > worst:
                        worst = e
        if worst > tol:
            warn("%s#%d: reconstruction err %.2e > %.2e -- storing all %d curves"
                 % (os.path.basename(path), rid, worst, tol, N))
            store_all = 1
        else:
            Gflat = [g6(G[i][j]) for i in range(d) for j in range(d)]

    if store_all:
        allc = []
        for j in range(N):
            for f in range(S):
                r = rows[f]
                allc.append([sum(V[a][k] * r[j * d + a] for a in range(d)) for k in range(d)])
        sc, q = quantise(allc, S * N, d)
    else:
        sc, q = quantise(curve0, S, d)

    out = {
        "f": label if label is not None else os.path.basename(path), "i": rid,
        "N": N, "d": d, "de": deff, "K": int(rec["K"]), "M": int(rec["M"]),
        "cov": int(rec.get("cover", 1)),
        "mo": int(rec["morse"]), "nu": int(rec["nullity"]),
        "A": float("%.10g" % rec["action"]), "_a": float(rec["action"]),
        "E": g6(rec["energy"]),
        "rms": g6(rec["rms"]), "mxr": g6(rec["maxr"]), "ms": g6(rec["minsep"]),
        "Ln": g6(rec["Lnorm"]), "tw": g6(rec["twist"]), "twr": g6(rec["twist_rel"]),
        "ck": int(rec.get("calib_k", 0)), "jr": g6(rec.get("jet_rel", 0.0)),
        "re": g6(rec["ret_err"]), "ce": g6(rec.get("coef_err", -1.0)),
        "hits": int(rec.get("hits", 0)), "seed": int(rec.get("seed", 0)),
        "trial": int(rec.get("trial", 0)), "secs": g6(rec.get("secs", 0.0)),
        "sym": rec.get("sym", ""), "rig": g6(rigid),
        "rot": 1 if rot else 0, "rotf": 1 if in_frame else 0,
        "rates": rates, "closed": closed,
        "ev": [g6(v) for v in ev],
        "Lsv": [g6(v) for v in rec.get("Lsv", [])],
        "pca": [g6(v) for v in rec.get("pca", [])],
        "S": S, "sc": [float("%.9g" % v) for v in sc], "q": q,
    }
    if Gflat is not None:
        out["G"] = Gflat
    if store_all:
        out["all"] = 1
    return out


def threading_ident():
    import threading
    return threading.get_ident()


# ==========================================================================
# the page
# ==========================================================================

CSS = r"""
:root{
  color-scheme: light dark;
  --bg:#f2f3f6; --tile:#ffffff; --canvas:#fafbfd; --panel:#eceef4;
  --rule:#a6adbd; --ink:#2a3140; --fg:#171b23; --fg2:#5d6675; --fg3:#89909f;
  --line:#dcdfe8; --ribbon:#4d5a72; --edge:#d5d9e3; --hdr:#ffffffee;
  --sel:#2f6fd0; --shadow:0 1px 2px rgba(20,26,40,.08),0 4px 14px rgba(20,26,40,.05);
  --b-rig:#c2492a; --b-rot:#8f6a10; --b-cov:#6252bd; --b-fam:#22718f; --b-warn:#b03a3a;
  --c2:#34619e; --c3:#1c7ba0; --c4:#12878c; --c5:#178a68; --c6:#2f8b3c;
  --c7:#5c8c1f; --c8:#7d8410; --c9:#8f760d; --c10:#9c6a0f; --c11:#a75513;
  --trail-a0:.34; --trail-a1:.52; --rtrail-a0:.16; --rtrail-a1:.40;
  --tw:168px;
}
/* Three theme states, not two: the un-stamped default resolves through the media query, while a
   host that stamps data-theme on the root must win in both directions. Tokens only, never rules. */
@media (prefers-color-scheme: dark){
:root:not([data-theme="light"]){
  --bg:#0d0e12; --tile:#16181d; --canvas:#171a20; --panel:#14161b;
  --rule:#34373f; --ink:#eef2f8; --fg:#e7eaf1; --fg2:#8d95a5; --fg3:#666e7e;
  --line:#262931; --ribbon:#93a7c4; --edge:#23262e; --hdr:#0d0e12ee;
  --sel:#5c9bff; --shadow:0 1px 2px rgba(0,0,0,.4);
  --b-rig:#e8785a; --b-rot:#d3a14a; --b-cov:#9b8cd6; --b-fam:#6fa8c9; --b-warn:#e06666;
  --c2:#4f7fbf; --c3:#4ba0cc; --c4:#4dbfc4; --c5:#55d2a4; --c6:#63dd7d;
  --c7:#8ce56a; --c8:#b9ec6d; --c9:#ddf278; --c10:#f6f59a; --c11:#ffe89c;
  --trail-a0:.30; --trail-a1:.45; --rtrail-a0:.13; --rtrail-a1:.32;
}
}
:root[data-theme="dark"]{
  --bg:#0d0e12; --tile:#16181d; --canvas:#171a20; --panel:#14161b;
  --rule:#34373f; --ink:#eef2f8; --fg:#e7eaf1; --fg2:#8d95a5; --fg3:#666e7e;
  --line:#262931; --ribbon:#93a7c4; --edge:#23262e; --hdr:#0d0e12ee;
  --sel:#5c9bff; --shadow:0 1px 2px rgba(0,0,0,.4);
  --b-rig:#e8785a; --b-rot:#d3a14a; --b-cov:#9b8cd6; --b-fam:#6fa8c9; --b-warn:#e06666;
  --c2:#4f7fbf; --c3:#4ba0cc; --c4:#4dbfc4; --c5:#55d2a4; --c6:#63dd7d;
  --c7:#8ce56a; --c8:#b9ec6d; --c9:#ddf278; --c10:#f6f59a; --c11:#ffe89c;
  --trail-a0:.30; --trail-a1:.45; --rtrail-a0:.13; --rtrail-a1:.32;
}
*{box-sizing:border-box}
html,body{margin:0;padding:0}
body{background:var(--bg);color:var(--fg);
  font:13px/1.45 ui-sans-serif,-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;
  -webkit-font-smoothing:antialiased}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace}
a{color:var(--sel)}
header{position:sticky;top:0;z-index:40;background:var(--hdr);backdrop-filter:blur(10px);
  border-bottom:1px solid var(--edge);padding:10px 16px 8px}
h1{font-size:15px;margin:0;font-weight:650;letter-spacing:-.01em;display:inline}
.sub{color:var(--fg2);font-size:11.5px;margin-left:10px}
.prov{color:var(--fg3);font-size:10.5px;margin-top:3px;word-break:break-all}
.ctl{display:flex;flex-wrap:wrap;gap:6px 12px;align-items:center;margin-top:8px}
.ctl label{color:var(--fg2);font-size:11px;display:flex;align-items:center;gap:5px}
input[type=search],select,input[type=number]{background:var(--tile);color:var(--fg);
  border:1px solid var(--edge);border-radius:5px;padding:3px 6px;font-size:11.5px;font-family:inherit}
input[type=search]{width:200px}
button{background:var(--tile);color:var(--fg);border:1px solid var(--edge);border-radius:5px;
  padding:3px 9px;font-size:11.5px;cursor:pointer;font-family:inherit}
button:hover{border-color:var(--fg3)}
button.on{background:var(--sel);border-color:var(--sel);color:#fff}
.chips{display:flex;flex-wrap:wrap;gap:4px}
.chip{padding:2px 7px;border:1px solid var(--edge);border-radius:11px;font-size:10.5px;
  cursor:pointer;background:var(--tile);color:var(--fg2);user-select:none;white-space:nowrap}
.chip.on{background:var(--sel);border-color:var(--sel);color:#fff}
.count{color:var(--fg2);font-size:11px;margin-left:auto}
details.legend{margin-top:7px;font-size:11px;color:var(--fg2)}
details.legend summary{cursor:pointer;color:var(--sel);font-size:11px;outline:none}
.legend .lg{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:8px 20px;
  margin-top:7px;padding:9px 11px;background:var(--tile);border:1px solid var(--edge);border-radius:7px}
.legend b{color:var(--fg);font-weight:600}
.legend .ramp{display:flex;gap:2px;align-items:center;flex-wrap:wrap;margin-top:3px}
.legend .ramp i{width:20px;height:9px;border-radius:2px;display:inline-block;font-style:normal}
main{padding:14px 16px 60px}
.sect{margin:0 0 6px;font-size:12px;font-weight:650;color:var(--fg2);
  border-bottom:1px solid var(--edge);padding-bottom:3px;margin-top:18px}
.sect:first-child{margin-top:0}
.sect em{font-style:normal;color:var(--fg3);font-weight:400}
/* --tw is the CANVAS width; the card adds its own 6px padding + 1px border on
   each side (box-sizing:border-box), so the card is --tw + 14. */
.grid{display:grid;grid-template-columns:repeat(auto-fill,calc(var(--tw) + 14px));
  gap:14px;justify-content:center}
.tile{width:calc(var(--tw) + 14px);background:var(--tile);border:1px solid var(--edge);
  border-radius:8px;
  padding:6px;box-shadow:var(--shadow);cursor:pointer;position:relative}
.tile:hover{border-color:var(--sel)}
.cw{position:relative;line-height:0;border-radius:4px;overflow:hidden;background:var(--canvas)}
.cw canvas{position:absolute;left:0;top:0}
.cw canvas.b{position:relative}
.tile .rb{margin-top:4px;border-radius:3px;overflow:hidden;line-height:0;background:var(--canvas)}
.mrow{margin-top:4px;font-size:10px;color:var(--fg2);line-height:1.35;
  display:flex;justify-content:space-between;gap:4px}
.mrow .id{color:var(--fg);font-weight:600}
.m2{font-size:10px;color:var(--fg2);line-height:1.35;margin-top:1px}
.m2 b{color:var(--fg);font-weight:600}
.bd{display:inline-block;padding:0 4px;border-radius:3px;font-size:9px;line-height:13px;
  border:1px solid currentColor;margin-left:3px;white-space:nowrap}
.deffdot{display:inline-block;width:7px;height:7px;border-radius:2px;vertical-align:-1px;margin-right:3px}
/* ---- detail overlay ---- */
#ov{position:fixed;inset:0;z-index:100;background:rgba(10,12,18,.6);backdrop-filter:blur(3px);
  display:none;overflow:auto;padding:24px}
#ov.on{display:block}
.card{max-width:1000px;margin:0 auto;background:var(--tile);border:1px solid var(--edge);
  border-radius:10px;padding:16px 18px;box-shadow:0 20px 60px rgba(0,0,0,.35)}
.card h2{margin:0 0 2px;font-size:15px;font-weight:650}
.card .x{position:absolute;right:26px;top:26px}
.dcols{display:grid;grid-template-columns:440px 1fr;gap:18px;margin-top:12px}
@media(max-width:900px){.dcols{grid-template-columns:1fr}}
.dcols canvas{border-radius:6px;background:var(--canvas);display:block;max-width:100%}
table.kv{width:100%;border-collapse:collapse;font-size:11px}
table.kv td{padding:1.5px 6px 1.5px 0;border-bottom:1px solid var(--line);vertical-align:top}
table.kv td:first-child{color:var(--fg2);white-space:nowrap;width:96px}
.spec{margin-top:8px}
.specrow{display:flex;align-items:center;gap:6px;font-size:10px;color:var(--fg2);margin:2px 0}
.specrow .bar{height:7px;border-radius:2px;min-width:1px}
.repro{margin-top:10px;font-size:10.5px;background:var(--panel);border:1px solid var(--edge);
  border-radius:6px;padding:7px 9px;color:var(--fg2);overflow-x:auto;white-space:pre}
.fam a{margin-right:8px;font-size:11px}
.hint{color:var(--fg3);font-size:10.5px;margin-top:4px}
.warnbox{margin-top:8px;padding:7px 10px;border:1px solid var(--b-warn);border-radius:6px;
  color:var(--b-warn);font-size:11px;background:color-mix(in srgb,var(--b-warn) 8%,transparent)}
"""


JS = r"""
'use strict';
var D    = JSON.parse(document.getElementById('DATA').textContent);
var RECS = D.recs, META = D.meta;
var DMAX = META.dmax;
var P    = Math.ceil(DMAX/2);            // strip panels
var MH   = 8, RH = 22;                   // meter height, ribbon height
var TW   = META.tile||168, PW = 32;      // tile width, strip panel width
function dprNow(){ return Math.min(window.devicePixelRatio||1, 2); }
var DPR  = dprNow();          /* refreshed on zoom / monitor change, see armDPR */
var PERIOD = 5.0;                        // seconds per orbit at speed 1
var MAXACT = 96;
var TAU  = Math.PI*2;
var TH   = {};

/* ---------------------------------------------------------------- theme */
function css(n){ return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
function hex2rgb(h){
  h=(h||'#888').trim(); if(h.charAt(0)==='#') h=h.slice(1);
  if(h.length===3) h=h[0]+h[0]+h[1]+h[1]+h[2]+h[2];
  return [parseInt(h.slice(0,2),16)||0,parseInt(h.slice(2,4),16)||0,parseInt(h.slice(4,6),16)||0];
}
function rgba(c,a){ return 'rgba('+c[0]+','+c[1]+','+c[2]+','+(a<0?0:a>1?1:a).toFixed(3)+')'; }
function readTheme(){
  TH = { canvas:css('--canvas'), panel:css('--panel'), rule:css('--rule'),
         ink:hex2rgb(css('--ink')), ribbon:hex2rgb(css('--ribbon')),
         ta0:parseFloat(css('--trail-a0')), ta1:parseFloat(css('--trail-a1')),
         ra0:parseFloat(css('--rtrail-a0')), ra1:parseFloat(css('--rtrail-a1')),
         ramp:{} };
  for(var k=2;k<=11;k++) TH.ramp[k]=hex2rgb(css('--c'+k));
}
function dcol(de){ var k=de<2?2:(de>11?11:de); return TH.ramp[k]; }
function dcolS(de){ return 'var(--c'+(de<2?2:(de>11?11:de))+')'; }

/* ------------------------------------------------------------- geometry */
var LE=(function(){var b=new ArrayBuffer(2);new DataView(b).setInt16(0,1,true);
                   return new Int16Array(b)[0]===1;})();
function b64i16(s){
  var bin=atob(s), n=bin.length, u8=new Uint8Array(n), i;
  for(i=0;i<n;i++) u8[i]=bin.charCodeAt(i);
  if(!LE){ for(i=0;i<n;i+=2){ var t=u8[i]; u8[i]=u8[i+1]; u8[i+1]=t; } }
  return new Int16Array(u8.buffer);
}
var YAW=34*Math.PI/180, PIT=24*Math.PI/180;

/* body j at frame f lives on curve ci(g,j) at frame fi(g,j,f) */
function ci(g,j){ return g.multi?j:0; }
function fi(g,j,f){ return g.multi?f:((f+j*g.shift)%g.S); }

function geom(r){
  if(r._g) return r._g;
  var S=r.S, d=r.d, N=r.N, i,j,f,k,a,b,c;
  var raw=b64i16(r.q), nc=r.all?N:1, tot=S*nc, sc=r.sc;
  var base=new Float32Array(tot*d);
  for(i=0;i<tot;i++){ var o=i*d; for(k=0;k<d;k++) base[o+k]=raw[o+k]*sc[k]; }

  var arrs=[], shift=0, multi=false;
  if(r.all){
    multi=true;
    for(j=0;j<N;j++) arrs.push(base.subarray(j*S*d,(j+1)*S*d));
  } else if(r.rot && r.G){
    /* body_j = G^j . body_0(f + j S/N):  the other N-1 curves are rebuilt here,
       which is why the payload only ships one curve. */
    multi=true; arrs.push(base);
    var Rm=new Float64Array(d*d), tmp=new Float64Array(d*d), G=r.G, step=S/N;
    for(i=0;i<d;i++) Rm[i*d+i]=1;
    for(j=1;j<N;j++){
      for(a=0;a<d;a++) for(b=0;b<d;b++){ var s=0;
        for(k=0;k<d;k++) s+=Rm[a*d+k]*G[k*d+b]; tmp[a*d+b]=s; }
      Rm.set(tmp);
      var pj=new Float32Array(S*d);
      for(f=0;f<S;f++){
        var src=((f+j*step)%S)*d, dst=f*d;
        for(a=0;a<d;a++){ var sv=0, ro=a*d;
          for(b=0;b<d;b++) sv+=Rm[ro+b]*base[src+b]; pj[dst+a]=sv; }
      }
      arrs.push(pj);
    }
  } else { arrs.push(base); shift=S/N; }
  nc = arrs.length;

  var g={S:S,d:d,N:N,arrs:arrs,multi:multi,shift:shift,nc:nc};

  /* main orthographic projection of principal axes 1,2,3 */
  var cy=Math.cos(YAW),sy=Math.sin(YAW),cp=Math.cos(PIT),sp=Math.sin(PIT);
  var mx=[],my=[],md=[], ext=1e-12, dmin=1e30, dmax=-1e30;
  for(c=0;c<nc;c++){
    var A=arrs[c], X=new Float32Array(S), Y=new Float32Array(S), Dp=new Float32Array(S);
    for(f=0;f<S;f++){
      var o=f*d, p1=A[o], p2=d>1?A[o+1]:0, p3=d>2?A[o+2]:0;
      var xx=p1*cy-p2*sy, yt=p1*sy+p2*cy;
      var ys=yt*sp+p3*cp, dd=yt*cp-p3*sp;
      X[f]=xx; Y[f]=ys; Dp[f]=dd;
      var m=xx<0?-xx:xx; if(m>ext)ext=m; m=ys<0?-ys:ys; if(m>ext)ext=m;
      if(dd<dmin)dmin=dd; if(dd>dmax)dmax=dd;
    }
    mx.push(X); my.push(Y); md.push(Dp);
  }
  g.mx=mx; g.my=my; g.md=md; g.ext=ext; g.dmin=dmin; g.dmax=dmax;

  /* per-axis extents, for the anisotropically normalised strip */
  var axe=new Float64Array(d);
  for(c=0;c<nc;c++){ var Ac=arrs[c];
    for(f=0;f<S;f++){ var oo=f*d;
      for(k=0;k<d;k++){ var v=Ac[oo+k]; if(v<0)v=-v; if(v>axe[k])axe[k]=v; } } }
  for(k=0;k<d;k++) if(axe[k]<1e-12) axe[k]=1e-12;
  g.axe=axe;

  /* live axes: decided from the eigenvalues, never from the drawn extent
     (the weakest live principal value gets down to 7.7e-5 of the first) */
  var lam0=r.ev[0]||1, live=[];
  for(k=0;k<2*P;k++) live.push(k<d && r.ev[k]>1e-5*lam0);
  g.live=live;

  /* mutual-distance ribbon: r_k(t)=|q_0-q_k|, k=1..floor(N/2) */
  var K=Math.floor(N/2)||1, rib=new Float32Array(K*S), rmax=1e-12;
  for(k=1;k<=K;k++){
    var Ak=arrs[ci(g,k)], A0=arrs[ci(g,0)];
    for(f=0;f<S;f++){
      var o0=fi(g,0,f)*d, ok=fi(g,k,f)*d, s2=0;
      for(a=0;a<d;a++){ var dv=A0[o0+a]-Ak[ok+a]; s2+=dv*dv; }
      var rv=Math.sqrt(s2); rib[(k-1)*S+f]=rv; if(rv>rmax)rmax=rv;
    }
  }
  g.rib=rib; g.ribK=K; g.ribMax=rmax;
  r._g=g; return g;
}

/* ------------------------------------------------------------- painting */
function mkCanvas(w,h,cls){
  var c=document.createElement('canvas');
  c.width=Math.max(1,Math.round(w*DPR)); c.height=Math.max(1,Math.round(h*DPR));
  c.style.width=w+'px'; c.style.height=h+'px';
  if(cls) c.className=cls;
  var x=c.getContext('2d'); x.setTransform(DPR,0,0,DPR,0,0);
  return c;
}
function killCanvas(c){ if(!c) return; c.width=0; c.height=0; if(c.parentNode) c.parentNode.removeChild(c); }

function paintMain(r,g,ctx,W){
  ctx.setTransform(DPR,0,0,DPR,0,0);
  ctx.fillStyle=TH.canvas; ctx.fillRect(0,0,W,W);
  var col=dcol(r.de), a0=r.rot?TH.ra0:TH.ta0, a1=r.rot?TH.ra1:TH.ta1;
  var NB=6, span=(g.dmax-g.dmin)||1, s=0.44*W/g.ext, cx=W/2, S=g.S, c,f,t;
  ctx.lineWidth=1; ctx.lineCap='round';
  for(c=0;c<g.nc;c++){
    var X=g.mx[c],Y=g.my[c],Dp=g.md[c], paths=[];
    for(t=0;t<NB;t++) paths.push(new Path2D());
    var prev=-1;
    for(f=0;f<S;f++){
      var f2=(f+1)%S;
      var bk=(((Dp[f]+Dp[f2])*0.5-g.dmin)/span*NB)|0;
      if(bk<0)bk=0; else if(bk>=NB)bk=NB-1;
      var pt=paths[bk];
      if(bk!==prev){ pt.moveTo(cx+X[f]*s, cx-Y[f]*s); prev=bk; }
      pt.lineTo(cx+X[f2]*s, cx-Y[f2]*s);
    }
    for(t=0;t<NB;t++){ ctx.strokeStyle=rgba(col,a0+(a1-a0)*(t/(NB-1))); ctx.stroke(paths[t]); }
  }
  return s;
}

function paintStrip(r,g,ctx,W,pw){
  ctx.setTransform(DPR,0,0,DPR,0,0);
  ctx.clearRect(0,0,W,pw+2+MH);
  var col=dcol(r.de), lam0=r.ev[0]||1, S=g.S, p,c,f,k;
  var al=r.rot?0.33:0.55;
  for(p=0;p<P;p++){
    var x0=p*(pw+2), a=2*p, b=2*p+1, la=g.live[a], lb=g.live[b], cxp=x0+pw/2, cyp=pw/2;
    ctx.fillStyle=TH.panel; ctx.fillRect(x0,0,pw,pw);
    if(la&&lb){
      var sa=0.42*pw/g.axe[a], sb=0.42*pw/g.axe[b];
      ctx.lineWidth=0.8; ctx.strokeStyle=rgba(col,al);
      for(c=0;c<g.nc;c++){
        var A=g.arrs[c]; ctx.beginPath();
        for(f=0;f<S;f++){
          var o=f*g.d, X=cxp+A[o+a]*sa, Y=cyp-A[o+b]*sb;
          if(f===0) ctx.moveTo(X,Y); else ctx.lineTo(X,Y);
        }
        ctx.closePath(); ctx.stroke();
      }
    } else if(la||lb){
      ctx.strokeStyle=rgba(col,0.30); ctx.lineWidth=1;
      ctx.beginPath(); ctx.moveTo(x0+pw*0.08,cyp); ctx.lineTo(x0+pw*0.92,cyp); ctx.stroke();
    } else {
      ctx.fillStyle=TH.rule; ctx.fillRect(cxp-pw*0.2, cyp-0.5, pw*0.4, 1);
    }
    /* meter: two bars, one per axis of this panel */
    var bw=Math.floor((pw-3)/2), y0=pw+2+MH;
    for(var q=0;q<2;q++){
      k=a+q; var bx=x0+q*(bw+3);
      if(k>=r.d){ ctx.fillStyle=TH.rule; ctx.fillRect(bx,y0-1,bw,1); continue; }
      var h=Math.max(1, MH*Math.sqrt(Math.max(0,r.ev[k])/lam0));
      ctx.fillStyle = g.live[k] ? rgba(col,0.85) : TH.rule;
      ctx.fillRect(bx, y0-h, bw, h);
    }
  }
}

function paintRibbon(r,g,ctx,W,H){
  ctx.setTransform(DPR,0,0,DPR,0,0);
  ctx.fillStyle=TH.canvas; ctx.fillRect(0,0,W,H);
  var S=g.S, K=g.ribK, mx=g.ribMax, k,f;
  ctx.lineWidth=1; ctx.strokeStyle=rgba(TH.ribbon,0.8);
  for(k=0;k<K;k++){
    ctx.beginPath();
    for(f=0;f<S;f++){
      var X=f/(S-1)*W, Y=H-1.5-g.rib[k*S+f]/mx*(H-3);
      if(f===0) ctx.moveTo(X,Y); else ctx.lineTo(X,Y);
    }
    ctx.stroke();
  }
}

function paintDots(t,f){
  var g=t.g, r=t.r, ctx=t.dctx, W=t.W, s=t.mscale, cx=W/2, j;
  ctx.clearRect(0,0,W,W);
  var col=dcol(r.de), k=W/168;
  for(j=0;j<r.N;j++){
    var c=ci(g,j), ff=fi(g,j,f);
    var x=cx+g.mx[c][ff]*s, y=cx-g.my[c][ff]*s;
    ctx.fillStyle=rgba(TH.ink,0.55+0.45*(1-j/r.N));
    ctx.beginPath(); ctx.arc(x,y,2.4*k,0,TAU); ctx.fill();
    if(j===0){ ctx.lineWidth=1; ctx.strokeStyle=rgba(col,0.95);
      ctx.beginPath(); ctx.arc(x,y,4.2*k,0,TAU); ctx.stroke(); }
  }
}

function paintSDots(t,f){
  var g=t.g, r=t.r, ctx=t.sdctx, pw=t.pw, W=t.W, p,j;
  ctx.clearRect(0,0,W,pw+2+MH);
  ctx.fillStyle=rgba(TH.ink,0.85);
  var rr=1.4, dd=rr*2;
  for(p=0;p<P;p++){
    var x0=p*(pw+2), a=2*p, b=2*p+1, la=g.live[a], lb=g.live[b];
    if(!la&&!lb) continue;
    var cxp=x0+pw/2, cyp=pw/2;
    for(j=0;j<r.N;j++){
      var A=g.arrs[ci(g,j)], o=fi(g,j,f)*g.d, X,Y;
      if(la&&lb){ X=cxp+A[o+a]*(0.42*pw/g.axe[a]); Y=cyp-A[o+b]*(0.42*pw/g.axe[b]); }
      else { var kk=la?a:b; X=cxp+A[o+kk]*(0.42*pw/g.axe[kk]); Y=cyp; }
      ctx.fillRect(X-rr,Y-rr,dd,dd);
    }
  }
}
"""


JS2 = r"""
/* ------------------------------------------------------------ tile DOM */
var TILES=[], ACT=[], VIS=new Set();
var REDUCED = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

function fmt(x,n){ if(x===0) return '0';
  var a=Math.abs(x);
  if(a<1e-4||a>=1e5) return x.toExponential(1);
  return x.toFixed(n===undefined?3:n); }
function esc(s){ return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }
function isRigid(r){ return r.rig<META.rigtol; }
function budget(r){ return 2*Math.floor(r.N/2); }

function badges(r){
  var h='';
  if(isRigid(r)) h+='<span class="bd" style="color:var(--b-rig)">rigid</span>';
  if(r.rotf) h+='<span class="bd" style="color:var(--b-rot)">frame {'+r.rates.join(',')+'}</span>';
  if(!r.closed) h+='<span class="bd" style="color:var(--b-rig)">open r.p.o.</span>';
  if(r.cov>1) h+='<span class="bd" style="color:var(--b-cov)">cover '+r.cov+'</span>';
  if(r.fam>0) h+='<span class="bd" style="color:var(--b-fam)">family '+r.fam+'</span>';
  if(r.re>1e-9||r.ce>1e-6) h+='<span class="bd" style="color:var(--b-warn)">residual '+fmt(Math.max(r.re,r.ce),2)+'</span>';
  return h;
}
function metaHTML(r){
  return '<div class="mrow"><span class="id">#'+r.i+'</span><span>'+badges(r)+'</span></div>'
   + '<div class="m2" style="color:var(--fg3)">'+esc(r.f)+'</div>'
   + '<div class="m2"><span class="deffdot" style="background:'+dcolS(r.de)+'"></span>'
   + 'N=<b>'+r.N+'</b> · d '+r.d+' · deff <b>'+r.de+'</b>/'+budget(r)+' · A <b>'+fmt(r.A,3)+'</b></div>'
   + '<div class="m2">morse '+r.mo+' · sep '+fmt(r.ms,3)+' · tw '+(Math.abs(r.tw)<1e-10?'0':fmt(r.tw,2))+'</div>';
}

function makeTile(r){
  var el=document.createElement('div');
  el.className='tile'; el.tabIndex=0;
  el.title=r.f+' #'+r.i+' — click for detail';
  var m=document.createElement('div'); m.className='cw main';
  var s=document.createElement('div'); s.className='cw strip';
  var b=document.createElement('div'); b.className='rb';
  var t={r:r, el:el, main:m, strip:s, rb:b, active:false, cvs:null, lastF:-1};
  el._t=t;
  el.appendChild(m); el.appendChild(s); el.appendChild(b);
  var meta=document.createElement('div'); meta.innerHTML=metaHTML(r);
  while(meta.firstChild) el.appendChild(meta.firstChild);
  el.addEventListener('click', function(){ openDetail(r); });
  el.addEventListener('keydown', function(e){ if(e.key==='Enter') openDetail(r); });
  sizeTile(t);
  r._t=t; TILES.push(t);
  prepObs.observe(el); visObs.observe(el);
  return t;
}
function sizeTile(t){
  var sh=PW+2+MH;
  t.main.style.width=TW+'px';  t.main.style.height=TW+'px';
  t.strip.style.width=TW+'px'; t.strip.style.height=sh+'px';
  t.rb.style.width=TW+'px';    t.rb.style.height=RH+'px';
}

function activate(t){
  if(t.active) return;
  var r=t.r, g;
  try { g=geom(r); } catch(e){ t.el.style.opacity=.4; t.active=true; return; }
  t.g=g; t.W=TW; t.pw=PW;
  var W=TW, sh=PW+2+MH;
  var mc=mkCanvas(W,W,'b'), dc=mkCanvas(W,W);
  t.main.appendChild(mc); t.main.appendChild(dc);
  t.mscale=paintMain(r,g,mc.getContext('2d'),W);
  t.dctx=dc.getContext('2d');
  var sc=mkCanvas(W,sh,'b'), sd=mkCanvas(W,sh);
  t.strip.appendChild(sc); t.strip.appendChild(sd);
  paintStrip(r,g,sc.getContext('2d'),W,PW);
  t.sdctx=sd.getContext('2d');
  var rc=mkCanvas(W,RH,'b'); t.rb.appendChild(rc);
  paintRibbon(r,g,rc.getContext('2d'),W,RH);
  t.cvs=[mc,dc,sc,sd,rc]; t.active=true; t.lastF=-1;
  ACT.push(t); trimActive(t);
  paintDots(t,0); paintSDots(t,0);
}
function deactivate(t){
  if(!t.active) return;
  if(t.cvs){ for(var i=0;i<t.cvs.length;i++) killCanvas(t.cvs[i]); }
  t.cvs=null; t.dctx=null; t.sdctx=null; t.active=false;
  if(t.r._g && t.r._g!==(OV.r&&OV.r._g)) t.r._g=null;   /* free ~100 KB per tile */
  t.g=null;
  var k=ACT.indexOf(t); if(k>=0) ACT.splice(k,1);
}
/* LRU trim.  `keep` is the tile currently being built: it is not in VIS yet
   (prepObs activates inside a 500px margin, i.e. below the fold), so without
   the guard a full ACT would evict it the instant it was pushed and the
   paintDots that follows would throw on a null context. */
function trimActive(keep){
  var i=0;
  while(ACT.length>MAXACT && i<ACT.length){
    var t=ACT[i];
    if(t===keep||VIS.has(t)){ i++; continue; }
    deactivate(t);
  }
}

var prepObs=new IntersectionObserver(function(es){
  for(var i=0;i<es.length;i++){
    var t=es[i].target._t; if(!t) continue;
    if(es[i].isIntersecting) activate(t);
    else if(!VIS.has(t)) deactivate(t);
  }
},{rootMargin:'500px 0px'});
var visObs=new IntersectionObserver(function(es){
  for(var i=0;i<es.length;i++){
    var t=es[i].target._t; if(!t) continue;
    if(es[i].isIntersecting){ VIS.add(t); if(!t.active) activate(t); }
    else VIS.delete(t);
  }
},{rootMargin:'0px'});


/* --------------------------------------------------- one global rAF loop */
var playing=!REDUCED, speed=1, phase=0, lastTs=0, rafId=0;
function loop(ts){
  rafId=requestAnimationFrame(loop);
  if(!lastTs) lastTs=ts;
  var dt=(ts-lastTs)/1000; lastTs=ts;
  if(dt>0.25) dt=0.25;
  if(playing){ phase+=dt*speed/PERIOD; phase-=Math.floor(phase); }
  VIS.forEach(function(t){
    if(!t.active) activate(t);
    if(!t.active||!t.dctx) return;
    var f=(phase*t.r.S)|0; if(f>=t.r.S) f=t.r.S-1;
    if(f===t.lastF) return;
    t.lastF=f; paintDots(t,f); paintSDots(t,f);
  });
}
function startLoop(){ if(!rafId){ lastTs=0; rafId=requestAnimationFrame(loop); } }
if(!REDUCED) startLoop();

/* -------------------------------------------------------- filter / sort */
var F={q:'',d:new Set(),N:new Set(),de:new Set(),rig:'all',rot:'all',
       sort:'A',asc:true,group:1};
function pass(r){
  if(F.d.size && !F.d.has(r.d)) return false;
  if(F.N.size && !F.N.has(r.N)) return false;
  if(F.de.size && !F.de.has(r.de)) return false;
  if(F.rig==='no' && isRigid(r)) return false;
  if(F.rig==='only' && !isRigid(r)) return false;
  if(F.rot==='rot' && !r.rotf) return false;
  if(F.rot==='in' && r.rotf) return false;
  if(F.q){
    var s=(r.f+' #'+r.i+' '+r.sym+' d'+r.d+' N'+r.N+' deff'+r.de).toLowerCase();
    if(s.indexOf(F.q)<0) return false;
  }
  return true;
}
function cmp(a,b){
  var k=F.sort, va,vb;
  if(k==='tw'){ va=Math.abs(a.tw); vb=Math.abs(b.tw); } else { va=a[k]; vb=b[k]; }
  if(va<vb) return F.asc?-1:1;
  if(va>vb) return F.asc?1:-1;
  return a.f<b.f?-1:a.f>b.f?1:a.i-b.i;
}
function render(){
  var list=RECS.filter(pass);
  list.sort(cmp);
  var host=document.getElementById('grid');
  host.textContent='';
  var frag=document.createDocumentFragment();
  if(F.group){
    var keys=[], map={};
    for(var i=0;i<list.length;i++){
      var r=list[i], k=r.d+'|'+r.N;
      if(!map[k]){ map[k]=[]; keys.push([r.d,r.N,k]); }
      map[k].push(r);
    }
    keys.sort(function(a,b){ return (b[0]-a[0]) || (a[1]-b[1]); });
    for(i=0;i<keys.length;i++){
      var kk=keys[i], arr=map[kk[2]];
      var h=document.createElement('div'); h.className='sect';
      h.innerHTML='d = '+kk[0]+' &nbsp; N = '+kk[1]+' <em>&nbsp;— '+arr.length+' record'+(arr.length===1?'':'s')+'</em>';
      frag.appendChild(h);
      frag.appendChild(gridOf(arr));
    }
  } else { frag.appendChild(gridOf(list)); }
  host.appendChild(frag);
  document.getElementById('count').textContent='showing '+list.length+' of '+RECS.length+' records';
}
function gridOf(arr){
  var g=document.createElement('div'); g.className='grid';
  for(var i=0;i<arr.length;i++){
    var r=arr[i], t=r._t||makeTile(r);
    g.appendChild(t.el);
  }
  return g;
}
"""


JS3 = r"""
/* ------------------------------------------------------- detail overlay */
var OV={r:null,g:null,ctx:null,W:420,yaw:YAW,pit:PIT,spin:true,ax:[0,1,2],
        raf:0,drag:false,px:0,py:0,ext:1,last:0,sd:null,spw:0,sww:0,lastF:-1,
        sx:null,sy:null,sz:null};

function axisOptions(r){
  var out=[], de=Math.max(3,r.de), k;
  for(k=0;k+2<Math.min(de,r.d);k++) out.push([k,k+1,k+2]);
  if(!out.length) out.push([0,1,2]);
  if(r.de>=6){ out.push([0,2,4]); out.push([1,3,5]); }
  if(r.de>=8){ out.push([0,3,6]); }
  return out;
}
function orbExtent(){
  var g=OV.g, a=OV.ax, m=1e-12, c,f;
  for(c=0;c<g.nc;c++){ var A=g.arrs[c];
    for(f=0;f<g.S;f++){ var o=f*g.d;
      var x=A[o+a[0]]||0, y=(a[1]<g.d?A[o+a[1]]:0)||0, z=(a[2]<g.d?A[o+a[2]]:0)||0;
      var n=Math.sqrt(x*x+y*y+z*z); if(n>m) m=n; } }
  OV.ext=m;
}
function paintOrbit(f){
  var r=OV.r, g=OV.g, ctx=OV.ctx, W=OV.W;
  ctx.setTransform(DPR,0,0,DPR,0,0);
  ctx.fillStyle=TH.canvas; ctx.fillRect(0,0,W,W);
  var cy=Math.cos(OV.yaw),sy=Math.sin(OV.yaw),cp=Math.cos(OV.pit),sp=Math.sin(OV.pit);
  var a=OV.ax, s=0.44*W/OV.ext, cx=W/2, col=dcol(r.de), S=g.S, c,i;
  var a0=r.rot?TH.ra0:TH.ta0, a1v=r.rot?TH.ra1:TH.ta1, NB=6;
  ctx.lineWidth=1.1; ctx.lineCap='round';
  if(!OV.sx||OV.sx.length<S){ OV.sx=new Float64Array(S); OV.sy=new Float64Array(S); OV.sz=new Float64Array(S); }
  var X=OV.sx, Y=OV.sy, Z=OV.sz;
  for(c=0;c<g.nc;c++){
    var A=g.arrs[c], dmin=1e30,dmax=-1e30;
    for(i=0;i<S;i++){
      var o=i*g.d;
      var p1=A[o+a[0]]||0, p2=(a[1]<g.d?A[o+a[1]]:0)||0, p3=(a[2]<g.d?A[o+a[2]]:0)||0;
      var xx=p1*cy-p2*sy, yt=p1*sy+p2*cy;
      X[i]=xx; Y[i]=yt*sp+p3*cp; Z[i]=yt*cp-p3*sp;
      if(Z[i]<dmin)dmin=Z[i]; if(Z[i]>dmax)dmax=Z[i];
    }
    var span=(dmax-dmin)||1, paths=[], t;
    for(t=0;t<NB;t++) paths.push(new Path2D());
    var prv=-1;
    for(i=0;i<S;i++){
      var i2=(i+1)%S, bk=(((Z[i]+Z[i2])*0.5-dmin)/span*NB)|0;
      if(bk<0)bk=0; else if(bk>=NB)bk=NB-1;
      if(bk!==prv){ paths[bk].moveTo(cx+X[i]*s, cx-Y[i]*s); prv=bk; }
      paths[bk].lineTo(cx+X[i2]*s, cx-Y[i2]*s);
    }
    for(t=0;t<NB;t++){ ctx.strokeStyle=rgba(col,a0+(a1v-a0)*(t/(NB-1))); ctx.stroke(paths[t]); }
  }
  for(var j=0;j<r.N;j++){
    var Ab=g.arrs[ci(g,j)], oo=fi(g,j,f)*g.d;
    var q1=Ab[oo+a[0]]||0, q2=(a[1]<g.d?Ab[oo+a[1]]:0)||0, q3=(a[2]<g.d?Ab[oo+a[2]]:0)||0;
    var bx=q1*cy-q2*sy, byt=q1*sy+q2*cy, by=byt*sp+q3*cp;
    var px=cx+bx*s, py=cx-by*s;
    ctx.fillStyle=rgba(TH.ink,0.55+0.45*(1-j/r.N));
    ctx.beginPath(); ctx.arc(px,py,4.2,0,TAU); ctx.fill();
    if(j===0){ ctx.lineWidth=1.4; ctx.strokeStyle=rgba(col,0.95);
      ctx.beginPath(); ctx.arc(px,py,7.5,0,TAU); ctx.stroke(); }
  }
}
function ovLoop(ts){
  OV.raf=requestAnimationFrame(ovLoop);
  if(!OV.last) OV.last=ts;
  var dt=(ts-OV.last)/1000; OV.last=ts; if(dt>0.25) dt=0.25;
  if(OV.spin && !OV.drag) OV.yaw+=dt*0.35;
  var f=(phase*OV.g.S)|0; if(f>=OV.g.S) f=OV.g.S-1;
  paintOrbit(f);
  if(OV.sd && f!==OV.lastF){ OV.lastF=f;
    paintSDots({g:OV.g,r:OV.r,sdctx:OV.sd,pw:OV.spw,W:OV.sww}, f); }
}
function closeDetail(){
  document.getElementById('ov').classList.remove('on');
  if(OV.raf) cancelAnimationFrame(OV.raf);
  OV.raf=0; OV.last=0;
  if(OV.r){ var keep=OV.r._t && OV.r._t.active; if(!keep) OV.r._g=null; }
  OV.r=null; OV.g=null; OV.ctx=null; OV.sd=null; OV.lastF=-1;
}
function famOf(r){ return RECS.filter(function(x){ return r.fam>0 && x.fam===r.fam && x!==r; }); }

function openDetail(r){
  var ov=document.getElementById('ov'), g;
  try{ g=geom(r); }catch(e){ return; }
  OV.r=r; OV.g=g; OV.spin=!REDUCED; OV.yaw=YAW; OV.pit=PIT;
  OV.ax=axisOptions(r)[0];
  var W=Math.min(420, Math.max(260, window.innerWidth-120)); OV.W=W;
  var i,k;
  var opts=axisOptions(r).map(function(a,ix){
     return '<option value="'+ix+'">principal axes '+(a[0]+1)+', '+(a[1]+1)+', '+(a[2]+1)+'</option>'; }).join('');
  var sumev=0; for(i=0;i<r.ev.length;i++) sumev+=Math.max(0,r.ev[i]);
  var cum=0, spec='';
  for(i=0;i<r.d;i++){
    var lv=Math.max(0,r.ev[i]); cum+=lv;
    var frac=lv/(r.ev[0]||1);
    spec+='<div class="specrow"><span style="width:22px">'+(i+1)+'</span>'
        + '<span class="bar" style="width:'+(Math.max(1,frac*160)).toFixed(1)+'px;background:'
        + (g.live[i]?dcolS(r.de):'var(--rule)')+'"></span>'
        + '<span class="mono">&sigma;='+Math.sqrt(lv).toPrecision(3)
        + ' &nbsp;&Sigma;='+(100*cum/(sumev||1)).toFixed(1)+'%</span></div>';
  }
  var fam=famOf(r);
  var path=META.paths[r.f]||('catalog/'+r.f);
  var kv=[
    ['file / id', esc(path)+' &nbsp; #'+r.i],
    ['N, d, deff', r.N+', '+r.d+', <b>'+r.de+'</b> of budget '+budget(r)+' (=2&lfloor;N/2&rfloor;)'],
    ['action', r.A],['energy', r.E],
    ['morse / nullity', r.mo+' / '+r.nu],
    ['min separation', r.ms],['rms / max r', fmt(r.rms,4)+' / '+fmt(r.mxr,4)],
    ['|L|', r.Ln],['twist / rel', fmt(r.tw,4)+' / '+fmt(r.twr,4)],
    ['rigidity', r.rig+(isRigid(r)?' &nbsp;<b style="color:var(--b-rig)">relative equilibrium</b>':'')],
    ['K / M / cover', r.K+' / '+r.M+' / '+r.cov],
    ['calib k / jet rel', r.ck+' / '+fmt(r.jr,3)],
    ['state / coef residual', fmt(r.re,3)+' / '+fmt(r.ce,3)],
    ['frame', r.rotf?('rotating, rates {'+r.rates.join(', ')+'}'
        +(r.closed?'':' — inertial motion NOT 2&pi;-periodic, so the co-rotating frame is drawn')):'inertial'],
    ['sym', r.sym?('<span class="mono">'+esc(r.sym)+'</span>'):'—'],
    ['seed / trial / hits', r.seed+' / '+r.trial+' / '+r.hits],
    ['loop pca[]', '<span class="mono" style="font-size:10px">'+r.pca.map(function(v){return (+v).toPrecision(3);}).join(' ')+'</span>'],
    ['Lsv[]', '<span class="mono" style="font-size:10px">'+r.Lsv.map(function(v){return (+v).toPrecision(3);}).join(' ')+'</span>']
  ].map(function(x){ return '<tr><td>'+x[0]+'</td><td>'+x[1]+'</td></tr>'; }).join('');

  ov.innerHTML='<div class="card"><button class="x" id="ovx">close ✕</button>'
    + '<h2>'+esc(r.f)+' &nbsp;#'+r.i+' '+badges(r)+'</h2>'
    + '<div style="color:var(--fg2);font-size:11px">N = '+r.N+' bodies in d = '+r.d
    + ', effective dimension <b style="color:'+dcolS(r.de)+'">'+r.de+'</b>, action '+fmt(r.A,6)+'</div>'
    + '<div class="dcols"><div>'
    +   '<div id="ovcv"></div>'
    +   '<div class="ctl" style="margin-top:6px">'
    +     '<select id="ovax">'+opts+'</select>'
    +     '<button id="ovspin"'+(REDUCED?'':' class="on"')+'>auto-spin</button>'
    +     '<span class="hint">drag to orbit</span></div>'
    +   '<div style="margin-top:10px"><div class="hint">principal-plane strip — '
    +     'deff = 2&times;(full panels) + (line panels)</div><div id="ovstrip"></div></div>'
    +   '<div style="margin-top:10px"><div class="hint">mutual distances |q<sub>0</sub>−q<sub>k</sub>|(t); '
    +     'flat lines = relative equilibrium. dashed = minsep '+fmt(r.ms,4)+'</div>'
    +     '<div id="ovrib"></div></div>'
    + '</div><div>'
    +   '<table class="kv">'+kv+'</table>'
    +   '<div class="spec"><div class="hint">pooled principal spectrum (over all '+r.N+' bodies)</div>'+spec+'</div>'
    +   (fam.length? '<div class="fam" style="margin-top:8px"><span class="hint">family siblings: </span>'
        + fam.map(function(x){ return '<a href="#" data-f="'+esc(x.f)+'" data-i="'+x.i+'">'+esc(x.f)+'#'+x.i+'</a>'; }).join('')+'</div>' : '')
    +   '<div class="repro">./hyperchoreography show '+esc(path)+' --id '+r.i+'\n'
    +     './hyperchoreography export '+esc(path)+' --id '+r.i+' --samples 720 --out curve.csv</div>'
    + (isRigid(r)? '<div class="warnbox">This is a <b>relative equilibrium</b> (rigidity '+r.rig
        +'): the configuration is frozen and merely rotates. However high its deff, it is dynamically trivial.</div>':'')
    + '</div></div></div>';
  ov.classList.add('on');

  var cv=mkCanvas(W,W); document.getElementById('ovcv').appendChild(cv);
  OV.ctx=cv.getContext('2d'); orbExtent();
  var PWo=Math.floor((W-2*(P-1))/P), sh=PWo+2+MH;
  var sc=mkCanvas(P*(PWo+2)-2,sh,'b'), sd=mkCanvas(P*(PWo+2)-2,sh);
  var sw=document.getElementById('ovstrip'); sw.className='cw';
  sw.style.width=(P*(PWo+2)-2)+'px'; sw.style.height=sh+'px'; sw.style.maxWidth='100%';
  sw.appendChild(sc); sw.appendChild(sd);
  paintStrip(r,g,sc.getContext('2d'),P*(PWo+2)-2,PWo);
  OV.sd=sd.getContext('2d'); OV.spw=PWo; OV.sww=P*(PWo+2)-2;
  var rw=document.getElementById('ovrib');
  var rc=mkCanvas(W,110,'b'); rw.appendChild(rc); rw.style.lineHeight=0;
  var rx=rc.getContext('2d'); paintRibbon(r,g,rx,W,110);
  rx.setTransform(DPR,0,0,DPR,0,0);
  rx.strokeStyle=rgba(TH.ink,0.35); rx.setLineDash([3,3]); rx.lineWidth=1;
  var yms=110-1.5-r.ms/g.ribMax*107;
  rx.beginPath(); rx.moveTo(0,yms); rx.lineTo(W,yms); rx.stroke(); rx.setLineDash([]);

  document.getElementById('ovx').onclick=closeDetail;
  document.getElementById('ovax').onchange=function(){ OV.ax=axisOptions(r)[+this.value]; orbExtent(); };
  document.getElementById('ovspin').onclick=function(){ OV.spin=!OV.spin; this.classList.toggle('on',OV.spin); };
  cv.style.cursor='grab';
  cv.addEventListener('pointerdown',function(e){ OV.drag=true; OV.px=e.clientX; OV.py=e.clientY;
     cv.setPointerCapture(e.pointerId); cv.style.cursor='grabbing'; });
  cv.addEventListener('pointermove',function(e){ if(!OV.drag) return;
     OV.yaw+=(e.clientX-OV.px)*0.01; OV.pit+=(e.clientY-OV.py)*0.01;
     if(OV.pit>1.5)OV.pit=1.5; if(OV.pit<-1.5)OV.pit=-1.5;
     OV.px=e.clientX; OV.py=e.clientY; });
  cv.addEventListener('pointerup',function(e){ OV.drag=false; cv.style.cursor='grab'; });
  var links=ov.querySelectorAll('.fam a');
  for(i=0;i<links.length;i++) links[i].onclick=function(e){
    e.preventDefault(); var ff=this.getAttribute('data-f'), ii=+this.getAttribute('data-i');
    for(var q=0;q<RECS.length;q++) if(RECS[q].f===ff && RECS[q].i===ii){ closeDetail(); openDetail(RECS[q]); return; }
  };
  OV.last=0; if(!OV.raf) OV.raf=requestAnimationFrame(ovLoop);
  paintOrbit(0);
}
"""


JS4 = r"""
/* ------------------------------------------------------------- controls */
function chipRow(host, vals, set, label){
  var el=document.getElementById(host); el.textContent='';
  var lab=document.createElement('span');
  lab.style.cssText='color:var(--fg2);font-size:11px;margin-right:2px'; lab.textContent=label;
  el.appendChild(lab);
  vals.forEach(function(v){
    var c=document.createElement('span'); c.className='chip'; c.textContent=v;
    c.onclick=function(){ if(set.has(v)) set.delete(v); else set.add(v);
      c.classList.toggle('on', set.has(v)); render(); };
    el.appendChild(c);
  });
}
function setTileWidth(w){
  TW=w; PW=Math.floor((TW-2*(P-1))/P); if(PW<12) PW=12;
  document.documentElement.style.setProperty('--tw', TW+'px');
  for(var i=0;i<TILES.length;i++){ deactivate(TILES[i]); sizeTile(TILES[i]); }
  requestAnimationFrame(function(){
    for(var i=0;i<TILES.length;i++) if(VIS.has(TILES[i])) activate(TILES[i]);
  });
}
function repaintAll(){
  var was=ACT.slice();
  for(var i=0;i<was.length;i++) deactivate(was[i]);
  for(i=0;i<was.length;i++) activate(was[i]);
}

/* devicePixelRatio changes when the page is zoomed or dragged to another
   display.  Every existing canvas keeps its old backing store, so rebuild them
   all rather than leaving 1x bitmaps upscaled on a 2x screen. */
function checkDPR(){
  var v=dprNow(); if(v===DPR) return;
  DPR=v; repaintAll();
  if(OV.r){ var r=OV.r; closeDetail(); openDetail(r); }
  armDPR();
}
function armDPR(){
  if(!window.matchMedia) return;
  var m=window.matchMedia('(resolution: '+(window.devicePixelRatio||1)+'dppx)');
  if(m.addEventListener) m.addEventListener('change', checkDPR, {once:true});
}
window.addEventListener('resize', checkDPR);

function init(){
  readTheme();
  armDPR();
  PW=Math.floor((TW-2*(P-1))/P); if(PW<12) PW=12;
  document.documentElement.style.setProperty('--tw', TW+'px');

  var ds={},ns={},des={};
  RECS.forEach(function(r){ ds[r.d]=1; ns[r.N]=1; des[r.de]=1; });
  var num=function(o){ return Object.keys(o).map(Number).sort(function(a,b){return a-b;}); };
  chipRow('fd', num(ds), F.d, 'd');
  chipRow('fn', num(ns), F.N, 'N');
  chipRow('fe', num(des), F.de, 'deff');

  document.getElementById('q').addEventListener('input', function(){
    F.q=this.value.trim().toLowerCase(); render(); });
  document.getElementById('sort').addEventListener('change', function(){ F.sort=this.value; render(); });
  document.getElementById('dir').addEventListener('click', function(){
    F.asc=!F.asc; this.textContent=F.asc?'▲ asc':'▼ desc'; render(); });
  document.getElementById('group').addEventListener('change', function(){ F.group=+this.value; render(); });
  document.getElementById('frig').addEventListener('change', function(){ F.rig=this.value; render(); });
  document.getElementById('frot').addEventListener('change', function(){ F.rot=this.value; render(); });
  document.getElementById('reset').addEventListener('click', function(){
    F.d.clear(); F.N.clear(); F.de.clear(); F.q=''; F.rig='all'; F.rot='all';
    document.getElementById('q').value='';
    document.getElementById('frig').value='all'; document.getElementById('frot').value='all';
    var cs=document.querySelectorAll('.chip'); for(var i=0;i<cs.length;i++) cs[i].classList.remove('on');
    render(); });
  var WS=[128,168,224];
  ['s','m','l'].forEach(function(k,ix){
    document.getElementById('sz'+k).classList.toggle('on', WS[ix]===TW); });
  ['s','m','l'].forEach(function(k,ix){
    document.getElementById('sz'+k).addEventListener('click', function(){
      var bs=document.querySelectorAll('.szb'); for(var i=0;i<bs.length;i++) bs[i].classList.remove('on');
      this.classList.add('on'); setTileWidth(WS[ix]); }); });
  var pb=document.getElementById('play');
  pb.textContent = playing?'⏸ pause':'▶ play';
  pb.addEventListener('click', function(){ playing=!playing;
    this.textContent=playing?'⏸ pause':'▶ play';
    if(playing) startLoop(); });
  document.getElementById('spd').addEventListener('input', function(){
    speed=+this.value; document.getElementById('spdv').textContent=speed.toFixed(1)+'×'; });

  document.getElementById('ov').addEventListener('click', function(e){ if(e.target===this) closeDetail(); });
  document.addEventListener('keydown', function(e){ if(e.key==='Escape') closeDetail(); });

  var mq=window.matchMedia('(prefers-color-scheme: dark)');
  var onTheme=function(){ readTheme(); repaintAll(); if(OV.r){ /* overlay repaints itself */ } };
  if(mq.addEventListener) mq.addEventListener('change', onTheme); else if(mq.addListener) mq.addListener(onTheme);
  // a host may switch theme by stamping data-theme on the root instead of changing the OS setting,
  // which fires no media query -- watch the attribute so the cached canvases repaint either way
  if(window.MutationObserver) new MutationObserver(onTheme).observe(document.documentElement,
    {attributes:true, attributeFilter:['data-theme']});

  render();
  if(REDUCED) document.getElementById('rmnote').style.display='inline';
}
init();
"""


PAGE = r"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>__CSS__</style>
</head><body>
<header>
  <h1>Hyperchoreography catalogue</h1><span class="sub">__SUB__</span>
  <div class="prov">__PROV__</div>
  <div class="ctl">
    <input type="search" id="q" placeholder="search file, id, sym…">
    <label>sort
      <select id="sort">
        <option value="A">action</option><option value="tw">|twist|</option>
        <option value="mo">morse index</option><option value="de">deff</option>
        <option value="N">N</option><option value="d">d</option>
        <option value="ms">min separation</option><option value="E">energy</option>
        <option value="hits">hits</option><option value="rig">rigidity</option>
        <option value="ce">coef residual</option><option value="i">id</option>
      </select></label>
    <button id="dir">▲ asc</button>
    <label>group
      <select id="group"><option value="1">by (d, N)</option><option value="0">none</option></select></label>
    <label>rigid
      <select id="frig"><option value="all">all</option><option value="no">hide</option>
        <option value="only">only</option></select></label>
    <label>frame
      <select id="frot"><option value="all">all</option><option value="rot">rotating</option>
        <option value="in">inertial</option></select></label>
    <button id="reset">reset</button>
    <span class="count" id="count"></span>
  </div>
  <div class="ctl">
    <span class="chips" id="fd"></span><span class="chips" id="fn"></span><span class="chips" id="fe"></span>
  </div>
  <div class="ctl">
    <span style="color:var(--fg2);font-size:11px">tile</span>
    <button class="szb" id="szs">S</button><button class="szb" id="szm">M</button><button class="szb" id="szl">L</button>
    <button id="play">⏸ pause</button>
    <label>speed <input type="range" id="spd" min="0.1" max="4" step="0.1" value="1" style="width:90px">
      <span id="spdv" class="mono" style="font-size:10.5px">1.0×</span></label>
    <span id="rmnote" style="display:none;color:var(--fg3);font-size:10.5px">reduced-motion: animation off by default</span>
  </div>
  <details class="legend"><summary>how to read a tile</summary><div class="lg">
    <div><b>Main view</b> — a static orthographic 3-D shadow of principal axes 1, 2, 3 of the
      <i>pooled</i> second moment (all N bodies, all samples). Depth is encoded as line opacity.
      The stored coordinate basis is meaningless (the action is O(d)-invariant), and a plain
      (1,2)-plane plot renders a large class of high-deff records as a featureless circle — hence 3-D.
      Dots are the bodies, brightest first, body 0 ringed.</div>
    <div><b>Strip</b> — one panel per principal plane (1,2), (3,4), (5,6)…, in three states:
      a <b>curve</b> = both axes live; a <b>horizontal line</b> = exactly one axis live;
      a <b>short grey rule</b> = neither (or the axis does not exist). So
      <b>deff = 2×(curve panels) + (line panels)</b>, countable by eye.
      Panels are normalised per axis, so shape here is qualitative only.</div>
    <div><b>Meter</b> — under each panel, two bars of height √(λ<sub>k</sub>/λ<sub>1</sub>):
      the magnitudes the strip's per-axis normalisation throws away. Grey = dead axis;
      a hairline = that axis does not exist at this d.</div>
    <div><b>Ribbon</b> — the mutual distances |q<sub>0</sub>−q<sub>k</sub>|(t), k = 1…⌊N/2⌋, over one
      period, floor = 0. Dead-flat lines mean a <b>relative equilibrium</b> (rigid body rotation);
      wavy lines mean a genuine choreography. The lowest minimum is <i>minsep</i>.</div>
    <div><b>deff colour</b>, fixed to 2…10 so it is stable across harvests:
      <span class="ramp" id="rampsw"></span></div>
    <div><b>Badges</b> — <span style="color:var(--b-rig)">rigid</span> relative equilibrium
      (rigidity &lt; 1e-4, a clean gap in the data — trivial however high its deff); <span style="color:var(--b-rot)">frame {…}</span> rotating
      frame with its rotation-rate multiset; <span style="color:var(--b-cov)">cover</span> multiple
      cover; <span style="color:var(--b-fam)">family</span> a continuous family (same N, same d, same action);
      <span style="color:var(--b-warn)">residual</span> a record whose certified state exceeds 1e-9 or
      whose stored coefficients exceed 1e-6.
      Click any tile for the full record.</div>
  </div></details>
</header>
<main><div id="grid"></div></main>
<div id="ov"></div>
<script type="application/json" id="DATA">__DATA__</script>
<script>__JS__
(function(){var s=document.getElementById('rampsw');if(!s)return;var h='';
for(var k=2;k<=11;k++)h+='<i style="background:var(--c'+k+')" title="deff '+k+'"></i>';
s.innerHTML=h+'<span style="margin-left:4px;color:var(--fg3)">2 → 11</span>';})();
</script>
</body></html>
"""


def to_fragment(html):
    """Strip the document wrapper, leaving <title>/<style>/body content/scripts.

    A published Artifact supplies its own <!doctype>/<html>/<head>/<body>, so the standalone
    document that is right for a local file is wrong there. Everything that carries meaning --
    the title, the stylesheet, the markup and the scripts -- survives; only the shell goes.
    """
    i = html.index("<title>")
    head = html[i:html.index("</head>")]
    body = html[html.index("</head><body>") + len("</head><body>"):html.rindex("</body></html>")]
    return head + body


def emit_page(records, meta, out_path, title, fragment=False):
    payload = json.dumps({"recs": records, "meta": meta},
                         separators=(",", ":"), sort_keys=False)
    payload = payload.replace("</", "<\\/")
    js = JS + JS2 + JS3 + JS4
    html = (PAGE.replace("__CSS__", CSS)
                .replace("__TITLE__", title)
                .replace("__SUB__", meta["sub"])
                .replace("__PROV__", meta["prov"])
                .replace("__DATA__", payload)
                .replace("__JS__", js))
    if fragment:
        html = to_fragment(html)
    parent = os.path.dirname(os.path.abspath(out_path))
    if parent and not os.path.isdir(parent):
        try:
            os.makedirs(parent)
        except OSError as e:
            raise CliError("cannot create output directory %s (%s)" % (parent, e))
    try:
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(html)
    except OSError as e:
        raise CliError("cannot write %s (%s)" % (out_path, e))
    return len(html.encode("utf-8"))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def find_binary(explicit):
    if explicit:
        return os.path.abspath(explicit)
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, os.pardir, "hyperchoreography"),
              os.path.join(os.getcwd(), "hyperchoreography")):
        c = os.path.abspath(c)
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return "hyperchoreography"


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build one self-contained HTML gallery from the choreography catalogue.")
    ap.add_argument("catalogs", nargs="*",
                    help="catalogue .bin files (default: catalog/*.bin)")
    ap.add_argument("--out", default="gallery.html")
    ap.add_argument("--samples", type=int, default=240,
                    help="target samples per period; rounded down to a multiple of N (default 240)")
    ap.add_argument("--max-records", type=int, default=0, help="0 = no limit")
    ap.add_argument("--min-deff", type=int, default=0)
    ap.add_argument("--no-rigid", action="store_true", help="drop relative equilibria")
    ap.add_argument("--tile", type=int, default=168, help="default tile width in px")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--bin", default=None, help="path to ./hyperchoreography")
    ap.add_argument("--title", default="Hyperchoreography catalogue")
    ap.add_argument("--fragment", action="store_true",
                    help="emit without the <!doctype>/<html>/<head>/<body> shell, for hosts that supply one")
    a = ap.parse_args(argv)

    t0 = time.time()
    binpath = find_binary(a.bin)
    paths = a.catalogs or sorted(glob.glob(os.path.join("catalog", "*.bin")))
    paths = [p for p in paths if p.endswith(".bin")]
    seen_real = set()
    uniq = []
    for p in paths:                      # the same file named twice is one file
        rp = os.path.realpath(p)
        if rp in seen_real:
            continue
        seen_real.add(rp)
        uniq.append(p)
    paths = uniq
    if not paths:
        print("gallery: no catalogue files found (looked for catalog/*.bin)", file=sys.stderr)
        return 2
    labels = path_labels(paths)

    warns = []
    def warn(msg):
        warns.append(msg)
        print("  warn: " + msg, file=sys.stderr)

    # ---- enumerate + build ------------------------------------------------
    # Ids are file positions: a search that saves a catalogue between the `list`
    # pass and the last `show`/`export` renumbers them.  `cat_show` rejects a
    # record whose id does not come back, and every file is fingerprinted around
    # its build so a rewrite that happens to keep ids valid is still caught and
    # that file is rebuilt from a fresh `list`.
    tmpdir = tempfile.mkdtemp(prefix="hcgal")
    built = {}                    # path -> [record, ...]     (replaced on retry)
    skips = {}                    # path -> ["<label>#id: why", ...]
    counts_by_path = {}           # path -> rows in `list`
    todo = list(paths)

    def work(job):
        p, rid, rigid = job
        try:
            return build_record(binpath, p, rid, rigid, a.samples, tmpdir, warn,
                                labels[p])
        except Exception as e:
            skips[p].append("%s#%s: %s" % (labels[p], rid, e))
            return None

    try:
        for attempt in range(1, 4):
            fps = {}
            jobs = []
            for p in todo:
                built[p] = []
                skips[p] = []
                fps[p] = fingerprint(p)
                try:
                    rows = cat_list(binpath, p)
                except Exception as e:
                    warn("%s: cannot list (%s)" % (p, e))
                    counts_by_path[p] = 0
                    continue
                counts_by_path[p] = len(rows)
                for row in rows:
                    jobs.append((p, row["id"], row["rigid"]))
            if a.max_records > 0:
                jobs = jobs[:a.max_records]
            if attempt == 1:
                print("gallery: %d records in %d catalogue file(s)"
                      % (len(jobs), len(paths)), file=sys.stderr)
            else:
                print("gallery: re-reading %d changed catalogue file(s), pass %d"
                      % (len(todo), attempt), file=sys.stderr)

            if a.jobs > 1:
                with ThreadPoolExecutor(max_workers=a.jobs) as ex:
                    out = list(ex.map(work, jobs))
            else:
                out = [work(j) for j in jobs]
            for job, r in zip(jobs, out):
                if r is not None:
                    built[job[0]].append(r)

            todo = [p for p in todo if fingerprint(p) != fps[p]]
            if not todo:
                break
            for p in todo:
                warn("%s CHANGED DURING THE BUILD (concurrent search?) -- "
                     "discarding its %d record(s)%s"
                     % (p, len(built[p]), " and re-reading it" if attempt < 3 else ""))
        else:
            for p in todo:
                warn("%s keeps changing; its records are LEFT OUT of the page" % p)
                built[p] = []
                counts_by_path[p] = 0
                skips[p] = ["%s: catalogue rewritten repeatedly during the build"
                            % labels[p]]
    finally:
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    recs = [r for p in paths for r in built.get(p, [])]
    skipped = [m for p in paths for m in skips.get(p, [])]
    counts = [(labels[p], counts_by_path.get(p, 0)) for p in paths]

    if a.min_deff:
        recs = [r for r in recs if r["de"] >= a.min_deff]
    if a.no_rigid:
        recs = [r for r in recs if not (r["rig"] < RIGID_TOL)]
    if not recs:
        print("gallery: no records survived", file=sys.stderr)
        return 1

    # deterministic order (threads finish out of order)
    recs.sort(key=lambda r: (r["f"], r["i"]))
    if a.max_records > 0:
        recs = recs[:a.max_records]

    # ---- continuous families: same N, same d, same action (README 13.1) ---
    # Grouped on the FULL-precision action, not the rounded display value, and
    # chained so a genuine one-parameter family stays one group.  The AMBIENT
    # dimension is part of the key: the C++ find_duplicate only dedupes inside
    # one catalogue file, so the same orbit re-found in a higher d (identical
    # action, minsep and deff, just embedded) shows up twice across files.  That
    # is a re-discovery, not a one-parameter family.
    order = sorted(range(len(recs)),
                   key=lambda k: (recs[k]["N"], recs[k]["d"], recs[k]["_a"]))
    fam = 0
    k = 0
    while k < len(order):
        j = k + 1
        while (j < len(order)
               and recs[order[j]]["N"] == recs[order[k]]["N"]
               and recs[order[j]]["d"] == recs[order[k]]["d"]
               and abs(recs[order[j]]["_a"] - recs[order[j - 1]]["_a"])
                   <= 1e-9 * max(1.0, abs(recs[order[k]]["_a"]))):
            j += 1
        if j - k > 1:
            fam += 1
            for t in range(k, j):
                recs[order[t]]["fam"] = fam
        k = j
    for r in recs:
        r.setdefault("fam", 0)
        del r["_a"]

    # ---- page metadata ---------------------------------------------------
    dmax = max(r["d"] for r in recs)
    dvals = sorted(set(r["d"] for r in recs))
    nvals = sorted(set(r["N"] for r in recs))
    evals = sorted(set(r["de"] for r in recs))
    nrot = sum(1 for r in recs if r["rot"])
    nrig = sum(1 for r in recs if r["rig"] < RIGID_TOL)
    nfam = fam
    sub = ("%d records · d %s · N %s · deff %s · %d rotating-frame · %d relative equilibria"
           % (len(recs), "/".join(map(str, dvals)), "/".join(map(str, nvals)),
              "%d–%d" % (evals[0], evals[-1]), nrot, nrig))
    prov = ("sources: " + ", ".join("%s (%d)" % (n, c) for n, c in counts if c)
            + " &nbsp;|&nbsp; %d samples/period, one curve + frame generator per record"
            % a.samples
            + (" &nbsp;|&nbsp; %d continuous famil%s" % (nfam, "y" if nfam == 1 else "ies") if nfam else "")
            + (" &nbsp;|&nbsp; %d skipped" % len(skipped) if skipped else "")
            + "<br>built by: python3 tools/gallery.py "
            + " ".join(sys.argv[1:] if argv is None else argv))
    meta = {"dmax": dmax, "rigtol": RIGID_TOL, "tile": a.tile, "sub": sub, "prov": prov,
            "paths": {labels[p]: p for p in paths}}

    try:
        nbytes = emit_page(recs, meta, a.out, a.title, fragment=a.fragment)
    except CliError as e:
        print("gallery: %s" % e, file=sys.stderr)
        return 1
    dt = time.time() - t0

    print("gallery: wrote %s -- %d records, %s bytes (%.2f MB), %.1f s"
          % (a.out, len(recs), format(nbytes, ","), nbytes / 1048576.0, dt), file=sys.stderr)
    if skipped:
        print("gallery: SKIPPED %d record(s):" % len(skipped), file=sys.stderr)
        for s in skipped:
            print("   " + s, file=sys.stderr)
    if warns:
        print("gallery: %d warning(s)" % len(warns), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
