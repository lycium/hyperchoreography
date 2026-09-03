"""Checks that the film is telling the truth."""

import subprocess
import sys
import tempfile
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from expo import catalog, nbody, optim, shoot, targets
from expo.pipeline import TimeFilter, tent_weights

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
for ss in (4, 16):
    tf = TimeFilter(ss, "tent")
    w = tent_weights(ss)
    n = 3 * ss + len(w)
    got = []
    for i in range(n):
        got += tf.push(np.full((4, 4, 4), (i * 7) % 256, np.uint8))
    want_n = (n - len(w)) // ss + 1
    check("ss %d: %d supersampled frames make %d" % (ss, n, want_n),
          len(got) == want_n, "%d" % len(got))
    worst = 0.0
    for k, data in enumerate(got):
        v = float(np.frombuffer(data, np.uint16)[0]) * 255.0 / 65535.0
        want = sum(a * ((i + k * ss) * 7 % 256)
                   for a, i in zip(w, range(len(w)))) / sum(w)
        worst = max(worst, abs(v - want))
    check("ss %d: every window is the tent mean of its frames" % ss, worst < 0.02,
          "worst %.4f code values" % worst)

print()
if fails:
    print("%d check(s) FAILED: %s" % (len(fails), ", ".join(fails)))
    sys.exit(1)
print("all checks passed")
