"""One search trial, end to end, with everything recorded.

This is src/search.hpp's run_trial minus the certification: a start, a Kepler
rescaling, phase 1 (L-BFGS), phase 2 (Newton-LM). It is deterministic in the seed,
so a scene can ask for the same descent every render.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nbody import Action
from .optim import Trace, lbfgs, newton_lm, inertia_gauge


@dataclass
class TrialRun:
    P: Action
    x0: np.ndarray                # the start, after Kepler rescaling
    x: np.ndarray                 # the critical point
    phase1: Trace
    phase2: Trace
    converged: bool
    morse_index: int = -1
    nullity: int = -1
    eig: np.ndarray | None = None


def random_start(P: Action, rng, K0=3, gamma=1.0) -> np.ndarray:
    """Low-mode random coefficients with amplitude m^-gamma, as random_guess does."""
    x = P.zeros()
    for mu, m in enumerate(P.modes):
        if m > K0:
            continue
        x[mu] = rng.normal(size=(2, P.d)) * m ** (-gamma)
    return x


def vertical_start(P: Action, rng, m1=1, m2=2) -> np.ndarray:
    """A circle in the first plane plus one transverse oscillation -- the
    unchained-polygon ansatz, the start that finds spatial orbits."""
    x = P.zeros()
    i1 = int(np.where(P.modes == m1)[0][0])
    x[i1, 0, 0] = 1.0
    x[i1, 1, 1] = 1.0
    if P.d > 2 and m2 in P.modes:
        i2 = int(np.where(P.modes == m2)[0][0])
        x[i2, 0, 2] = 0.35 * rng.uniform(0.6, 1.4)
    return x


def run(N=3, d=2, K=16, seed=0, phase1="action", n1=120, newton=40,
        start=None, gtol=1e-11, M=None) -> TrialRun:
    """A trial. `phase1` is "action" (descend on A) or "gradnorm" (descend on
    |grad A|^2/2, which is what reaches saddles); `start` overrides the random one."""
    P = Action(N, d, K, M=M)
    rng = np.random.default_rng(seed)

    x = random_start(P, rng) if start is None else np.array(start, dtype=float)
    for _ in range(20):                       # reject colliding starts, as the C++ does
        if np.isfinite(P.value(x)):
            break
        x = random_start(P, rng)
    x = x * P.optimal_scale(x)                # Kepler rescaling to the stationary size
    x0 = P.flat(x).copy()

    def fA(v):
        A, g = P.value(v, grad=True)
        return A, P.flat(g)

    def fG(v):
        """half |grad A|^2, whose gradient is H grad A -- reaches critical points of any index"""
        A, g = P.value(v, grad=True)
        if not np.isfinite(A):
            return np.inf, np.zeros(P.n)
        gf = P.flat(g)
        return 0.5 * float(gf @ gf), P.hessian(v) @ gf

    act = lambda v: P.value(v)
    t1 = Trace("L-BFGS on " + ("A" if phase1 == "action" else "|grad A|^2 / 2"))
    x1, _ = lbfgs(x0, fA if phase1 == "action" else fG,
                  max_iter=n1, gtol=gtol if phase1 == "action" else gtol ** 2,
                  trace=t1, action_of=act)

    t2 = Trace("Newton-LM on grad A = 0")
    x2, ok = newton_lm(x1, fA, lambda v: P.hessian(v), max_iter=newton, gtol=gtol,
                       trace=t2, action_of=act)

    neg, zero, w = inertia_gauge(P, x2)
    return TrialRun(P=P, x0=x0, x=x2, phase1=t1, phase2=t2, converged=ok,
                    morse_index=neg, nullity=zero, eig=w)


def landscape(P: Action, x, u, v, span=1.0, n=81):
    """A[x + a u + b v] on a grid -- the real surface, not a sketch."""
    a = np.linspace(-span, span, n)
    b = np.linspace(-span, span, n)
    Z = np.empty((n, n))
    xf = P.flat(x)
    for i, bi in enumerate(b):
        for j, aj in enumerate(a):
            Z[i, j] = P.value(xf + aj * u + bi * v)
    return a, b, Z
