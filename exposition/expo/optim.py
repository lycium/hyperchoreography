"""L-BFGS and the Levenberg-Marquardt Newton, ported from src/optim.hpp.

These are the two phases of every search trial, and they are reproduced here line
for line -- same ring buffers, same Armijo constant, same damping schedule -- so
the presentation animates the solver's own descent rather than an impression of it.
Each call records every iterate, which is what the scenes play back.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Trace:
    """Everything one optimisation phase did, one entry per iteration."""
    phase: str
    x: list = field(default_factory=list)          # iterates, flattened
    f: list = field(default_factory=list)          # objective
    gnorm: list = field(default_factory=list)      # |grad|
    action: list = field(default_factory=list)     # A, even when descending on |grad|^2
    extra: list = field(default_factory=list)      # per-phase: mu, or the step length

    def add(self, x, f, gnorm, action=None, extra=None):
        self.x.append(np.array(x, dtype=float, copy=True))
        self.f.append(float(f))
        self.gnorm.append(float(gnorm))
        self.action.append(float(action if action is not None else f))
        self.extra.append(extra)

    def __len__(self):
        return len(self.x)


def lbfgs(x0, fn, max_iter=200, gtol=1e-10, mem=10, trace=None, action_of=None):
    """Limited-memory BFGS with an Armijo backtracking line search.

    `fn(x) -> (f, g)`. Returns (x, converged). The two-loop recursion, the
    Barzilai-Borwein initial scaling and the 1e-4 sufficient-decrease constant are
    the ones in optim.hpp.
    """
    x = np.array(x0, dtype=float).reshape(-1).copy()
    n = x.size
    S = np.zeros((mem, n))
    Y = np.zeros((mem, n))
    rho = np.zeros(mem)
    m, head = 0, 0

    f, g = fn(x)
    if not np.isfinite(f):
        return x, False
    gn = float(np.linalg.norm(g))
    if trace is not None:
        trace.add(x, f, gn, action_of(x) if action_of else f, 1.0)

    for _ in range(max_iter):
        if gn <= gtol:
            break
        # two-loop recursion, oldest to newest
        idx = [(head - m + i) % mem for i in range(m)]
        q = g.copy()
        al = np.zeros(m)
        for i in reversed(range(m)):
            k = idx[i]
            al[i] = rho[k] * S[k] @ q
            q -= al[i] * Y[k]
        if m:
            k = idx[-1]
            q *= (S[k] @ Y[k]) / (Y[k] @ Y[k])
        else:
            q /= max(gn, 1e-30)
        for i in range(m):
            k = idx[i]
            q += (al[i] - rho[k] * (Y[k] @ q)) * S[k]
        d = -q

        gd = float(g @ d)
        if not gd < 0:                              # not a descent direction: reset
            d = -g / max(gn, 1e-30)
            gd = -gn
            m = 0

        step, ok = 1.0, False
        for _ls in range(40):
            xn = x + step * d
            fn_, gnw = fn(xn)
            if np.isfinite(fn_) and fn_ <= f + 1e-4 * step * gd:
                ok = True
                break
            step *= 0.5
        if not ok:
            break

        s, y = xn - x, gnw - g
        sy, ss, yy = float(s @ y), float(s @ s), float(y @ y)
        if sy > 1e-12 * np.sqrt(ss * yy):
            S[head], Y[head], rho[head] = s, y, 1.0 / sy
            head = (head + 1) % mem
            m = min(m + 1, mem)
        x, f, g = xn, fn_, gnw
        gn = float(np.linalg.norm(g))
        if trace is not None:
            trace.add(x, f, gn, action_of(x) if action_of else f, step)

    return x, gn <= gtol


def newton_lm(x0, fn, hess, max_iter=60, gtol=1e-10, trace=None, action_of=None):
    """Newton on grad A = 0, damped in the Hessian eigenbasis.

        delta = - sum_k  q_k * lambda_k / (lambda_k^2 + mu) * (q_k . g)

    The sign of lambda_k never appears, so this converges quadratically to a
    critical point of *any* Morse index -- which is the whole point, since most
    choreographies are saddles. Directions with lambda_k ~ 0, the gauge freedoms,
    receive almost no step and so are ignored without being identified.

    mu starts at 1e-4 lambda_max^2, shrinks by 5 on an accepted step and grows by
    8 on a rejected one; a step is accepted when it reduces |grad|, not A.
    """
    x = np.array(x0, dtype=float).reshape(-1).copy()
    f, g = fn(x)
    if not np.isfinite(f):
        return x, False
    gn = float(np.linalg.norm(g))
    mu = -1.0
    if trace is not None:
        trace.add(x, f, gn, action_of(x) if action_of else f, None)

    for _ in range(max_iter):
        if gn <= gtol:
            break
        H = hess(x)
        w, V = np.linalg.eigh(H)
        lmax = float(np.abs(w).max())
        if mu < 0:
            mu = 1e-4 * lmax * lmax
        gq = V.T @ g
        accepted = False
        for _try in range(25):
            delta = -(V @ (w / (w * w + mu) * gq))
            xn = x + delta
            fnw, gnw = fn(xn)
            gnn = float(np.linalg.norm(gnw)) if np.isfinite(fnw) else np.inf
            if gnn < gn:
                accepted = True
                x, f, g, gn = xn, fnw, gnw, gnn
                mu = max(mu / 5.0, 1e-18 * lmax * lmax)
                break
            mu *= 8.0
        if not accepted:
            break
        if trace is not None:
            trace.add(x, f, gn, action_of(x) if action_of else f,
                      {"mu": mu, "eig": w.copy()})

    return x, gn <= gtol


def morse(H, rel_tol=1e-8):
    """(negative, zero, positive) eigenvalue counts -- the Morse index and nullity."""
    w = np.linalg.eigvalsh(H)
    lmax = float(np.abs(w).max())
    neg = int(np.sum(w < -rel_tol * lmax))
    zero = int(np.sum(np.abs(w) <= rel_tol * lmax))
    return neg, zero, len(w) - neg - zero, w


# ---------------------------------------------------------------------------
# the gauge directions, and the Morse index that ignores them
# ---------------------------------------------------------------------------
def gauge_basis(P, x, omega=None):
    """Orthonormal directions along which the action cannot change.

    Two families. Shifting the time origin moves the loop along itself, and
    rotating it rigidly in R^d leaves every mutual distance alone, so both are
    exact zero modes of the Hessian at any critical point. In a rotating frame only
    the rotations that commute with Omega survive; that centraliser is the kernel of
    the map L -> [Omega, L], not the coordinate generators taken one at a time.
    """
    P_ = P
    x = P_.shaped(x)
    d, nm = P_.d, P_.nm
    m = P_.modes[:, None].astype(float)

    # the time shift: d/dtau of (c cos m(t+tau) + s sin m(t+tau)) at tau = 0
    shift = np.empty_like(x)
    shift[:, 0, :] = m * x[:, 1, :]
    shift[:, 1, :] = -m * x[:, 0, :]
    gens = [shift.reshape(-1)]

    pairs = [(a, b) for a in range(d) for b in range(a + 1, d)]
    if omega is None or not np.any(omega):
        basis = np.eye(len(pairs))
    else:                                   # kernel of ad_Omega in the E_ab basis
        D = len(pairs)
        M = np.zeros((D, D))
        for r, (a, b) in enumerate(pairs):
            for c, (e, f) in enumerate(pairs):
                v = 0.0
                if b == f:
                    v += omega[a, e]
                if b == e:
                    v -= omega[a, f]
                if a == e:
                    v -= omega[f, b]
                if a == f:
                    v += omega[e, b]
                M[r, c] = v
        w, V = np.linalg.eigh(M.T @ M)
        basis = V[:, np.abs(w) <= 1e-12 * max(np.abs(w).max(), 1e-300)].T

    for row in np.atleast_2d(basis):
        v = np.zeros_like(x)
        for c, (a, b) in enumerate(pairs):
            if row[c] == 0.0:
                continue
            v[:, :, a] -= row[c] * x[:, :, b]
            v[:, :, b] += row[c] * x[:, :, a]
        gens.append(v.reshape(-1))

    # Gram-Schmidt, dropping directions the loop does not actually have
    xn = float(np.linalg.norm(x))
    out = []
    for g in gens:
        for u in out:
            g = g - (u @ g) * u
        n = float(np.linalg.norm(g))
        if n < 1e-8 * xn:
            continue
        out.append(g / n)
    return np.array(out).T if out else np.zeros((P_.n, 0))


def inertia_gauge(P, x, H=None, omega=None):
    """(Morse index, nullity, eigenvalues) with the gauge directions lifted out.

    Adding sigma G G^T pushes the exact zero modes to the top of the spectrum
    instead of leaving them near zero to be classified by a magnitude threshold --
    which is fragile, because a genuinely soft direction and an exact gauge
    direction look the same to a tolerance. What remains is judged against |H G|,
    the size of the numerical error in those directions.
    """
    if H is None:
        H = P.hessian(x)
    G = gauge_basis(P, x, omega)
    ng = G.shape[1]
    sg = float(np.abs(np.diag(H)).max())
    tol = float(np.abs(H @ G).max()) if ng else 0.0
    tol = max(tol, 16 * 2.220446049250313e-16 * sg)
    w = np.linalg.eigvalsh(H + sg * (G @ G.T))
    neg = int(np.sum(w < -tol))
    zero = int(np.sum(np.abs(w) <= tol)) + ng
    return neg, zero, w
