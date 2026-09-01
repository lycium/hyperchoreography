#!/usr/bin/env python3
"""Checks that the film is telling the truth.

Three things could go wrong quietly. The numerics could drift from the solver's,
and every number on screen would be subtly wrong. The Morse index could be
computed the naive way, and the counts would disagree with the catalogue. Or the
time filter could stop being a tent, and the motion blur would be wrong in a way
nobody would notice by eye.

    .venv/bin/python selftest.py
"""

import subprocess
import sys
import tempfile
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from expo import catalog, nbody, optim, shoot, targets
from expo.pipeline import filter_chain, tent_weights

fails = []


def check(name, ok, detail=""):
    print("  %-46s %s  %s" % (name, "ok " if ok else "FAIL", detail))
    if not ok:
        fails.append(name)


RECORDS = ("eight", "n4_planar", "hiphop", "n4_deff4", "n6_deff6",
           "d7_inertial", "d7_twist", "d11")

print("action and gradient against the catalogue")
for name in RECORDS:
    rec = catalog.load(name)
    P = nbody.Action(rec.N, rec.d, int(rec.modes.max()), modes=rec.modes,
                     omega=rec.omega)
    A, g = P.value(rec.coef, grad=True)
    check("%s: action" % name, abs(A - rec.action) < 1e-8 * abs(rec.action),
          "%.9f vs %.9f" % (A, rec.action))
    check("%s: it is a critical point" % name, np.linalg.norm(g) < 1e-9,
          "|grad| = %.1e" % np.linalg.norm(g))

print("Morse index and nullity against the catalogue")
for name in RECORDS:
    rec = catalog.load(name)
    P = nbody.Action(rec.N, rec.d, int(rec.modes.max()), modes=rec.modes,
                     omega=rec.omega)
    neg, zero, _ = optim.inertia_gauge(P, rec.coef, omega=rec.omega)
    check("%s: morse %d nullity %d" % (name, neg, zero),
          neg == rec.morse and zero == rec.nullity,
          "stored %d / %d" % (rec.morse, rec.nullity))

print("the Hessian against finite differences of the gradient")
P = nbody.Action(3, 2, 12)
rng = np.random.default_rng(3)
x = rng.normal(size=P.n) * 0.3
x *= P.optimal_scale(x)
H = P.hessian(x)
Hn = np.empty_like(H)
h = 1e-6
for i in range(P.n):
    xp, xm = x.copy(), x.copy()
    xp[i] += h
    xm[i] -= h
    Hn[:, i] = (P.flat(P.gradient(xp)) - P.flat(P.gradient(xm))) / (2 * h)
rel = np.abs(H - Hn).max() / np.abs(H).max()
check("relative agreement", rel < 1e-8, "%.1e" % rel)

print("the featured descents land where the catalogue says")
for feat, want, morse in ((targets.pentagon(), 17.447668440, 2),
                          (targets.hiphop(), 26.761760441, 8),
                          (targets.saddle_n3(), 11.152080126, 1)):
    got = feat.run.phase2.action[-1]
    check("%s: action %.9f" % (feat.key, got), abs(got - want) < 1e-8, "want %.9f" % want)
    check("%s: Morse index %d" % (feat.key, feat.run.morse_index),
          feat.run.morse_index == morse, "want %d" % morse)

print("k-fold covers scale the action by k^(2/3)")
P = nbody.Action(3, 2, 24)
base = None
from expo.trial import run as trial_run
for k in (1, 2, 4):
    r = trial_run(N=3, d=2, K=24, phase1="action", n1=100, newton=60,
                  start=targets.circle_start(P, k))
    A = r.phase2.action[-1]
    if base is None:
        base = A
    check("%d-fold cover" % k, abs(A / base - k ** (2 / 3.)) < 1e-7,
          "ratio %.9f vs %.9f" % (A / base, k ** (2 / 3.)))

print("the shooting Newton reproduces the record's own residual")
rec = catalog.load("eight")
P = nbody.Action(rec.N, 2, int(rec.modes.max()), modes=rec.modes)
Z = shoot.state_from_loop(P, np.array(rec.coef[:, :, :2]))
res = float(np.abs(shoot.shift_residual(Z, 3, 2, 2 * np.pi, 3000)).max())
check("shift residual of the stored coefficients", 1e-11 < res < 1e-10,
      "%.1e (the record stores 3.5e-11)" % res)

print("the time filter is exactly a tent over overlapping windows")
check("tent weights", tent_weights(4) == [1, 2, 3, 4, 3, 2, 1], str(tent_weights(4)))
with tempfile.TemporaryDirectory() as d:
    from PIL import Image
    for i in range(8):                       # a ramp, so the answer is checkable
        Image.fromarray(np.full((16, 16, 3), i * 36, np.uint8)).save(
            os.path.join(d, "f%02d.png" % i))
    out = os.path.join(d, "out.mkv")
    got = os.path.join(d, "got.png")
    vf = filter_chain(4, 30, 16, 16, 16, 16)
    subprocess.run(["ffmpeg", "-hide_banner", "-loglevel", "error", "-y",
                    "-framerate", "120", "-i", os.path.join(d, "f%02d.png"),
                    "-vf", vf, "-r", "30", "-c:v", "ffv1", out], check=True)
    subprocess.run(["ffmpeg", "-hide_banner", "-loglevel", "error", "-y",
                    "-i", out, "-frames:v", "1", got], check=True)
    v = float(np.array(Image.open(got))[0, 0, 0])
    w = tent_weights(4)
    want = sum(a * b for a, b in zip(w, [i * 36 for i in range(7)])) / sum(w)
    check("first output frame is the tent mean", abs(v - want) < 0.5,
          "%g, want %g" % (v, want))

print()
if fails:
    print("%d check(s) FAILED: %s" % (len(fails), ", ".join(fails)))
    sys.exit(1)
print("all checks passed")
