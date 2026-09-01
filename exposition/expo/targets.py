"""The orbits the presentation follows, and the optimiser runs that reach them.

Each entry is a real descent: a start that means something, the two phases the
solver runs, and the critical point they land on. Nothing here is staged -- the
numbers printed on screen are the ones these functions return.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import catalog
from .nbody import Action
from .optim import inertia_gauge
from .trial import TrialRun, run


# ---------------------------------------------------------------------------
# putting a catalogued record into a working problem
# ---------------------------------------------------------------------------
def embed(orbit, P: Action) -> np.ndarray:
    """A record's coefficients in P's basis, dropping modes P does not carry."""
    x = P.zeros()
    for k, m in enumerate(orbit.modes):
        j = np.where(P.modes == m)[0]
        if len(j):
            x[j[0], :, :] = orbit.coef[k][:, :P.d]
    return x


def principal_frame(P: Action, x) -> np.ndarray:
    """Rows are the principal axes of the loop -- the frame every picture uses."""
    X = P.shaped(x).reshape(-1, P.d)
    w, V = np.linalg.eigh(X.T @ X)
    V = V[:, ::-1]
    for c in range(V.shape[1]):
        k = np.argmax(np.abs(V[:, c]))
        if V[k, c] < 0:
            V[:, c] = -V[:, c]
    return V.T


# ---------------------------------------------------------------------------
# starts that mean something
# ---------------------------------------------------------------------------
def circle_start(P: Action, m: int = 1, plane=(0, 1)) -> np.ndarray:
    """A body turning at constant rate in one plane: the rotating N-gon."""
    x = P.zeros()
    i = int(np.where(P.modes == m)[0][0])
    x[i, 0, plane[0]] = 1.0
    x[i, 1, plane[1]] = 1.0
    return x


def epicycle_start(P: Action, m1: int = 2, m2: int = 3, a: float = 1.0,
                   b: float = 0.35) -> np.ndarray:
    """Two circles turning opposite ways: exp(i m1 t) plus b exp(-i m2 t).

    Their difference of rates is m1 + m2, so the figure closes up after that many
    lobes -- which is where the five-fold pentagon comes from at (2, 3).
    """
    x = P.zeros()
    i1 = int(np.where(P.modes == m1)[0][0])
    i2 = int(np.where(P.modes == m2)[0][0])
    x[i1, 0, 0] = a
    x[i1, 1, 1] = a
    x[i2, 0, 0] = b
    x[i2, 1, 1] = -b
    return x


def vertical_start(P: Action, m1: int = 1, m2: int = 2, amp: float = 0.35,
                   rng=None) -> np.ndarray:
    """A rotating polygon plus one oscillation out of its plane.

    This is the unchained-polygon ansatz; it converges when m2/m1 is near one of the
    N-gon's transverse frequencies.
    """
    x = circle_start(P, m1)
    i2 = int(np.where(P.modes == m2)[0][0])
    x[i2, 0, 2] = amp
    return x


# ---------------------------------------------------------------------------
# the featured descents
# ---------------------------------------------------------------------------
@dataclass
class Featured:
    key: str
    title: str
    blurb: str
    run: TrialRun
    orbit: object | None = None

    @property
    def P(self):
        return self.run.P

    @property
    def x(self):
        return self.run.x

    def frame(self):
        return principal_frame(self.P, self.x)

    def bodies(self, x=None):
        P, xx = self.P, (self.x if x is None else x)
        return lambda ts: P.bodies(xx, ts)


def eight(seed: int = 1, K: int = 24) -> Featured:
    """N = 3, from a random low-mode start, descending on the action."""
    r = run(N=3, d=2, K=K, seed=seed, phase1="action", n1=200, newton=60)
    return Featured("eight", "the figure eight",
                    "N = 3, from a random start, descending on the action", r)


def saddle_n3(seed: int = 39, K: int = 24) -> Featured:
    """N = 3, index 1: a saddle, so descent on the action alone cannot find it."""
    r = run(N=3, d=2, K=K, seed=seed, phase1="gradnorm", n1=400, newton=120)
    return Featured("saddle_n3", "an N = 3 saddle",
                    "index 1: the action goes down in one direction", r)


def pentagon(K: int = 48, b: float = 0.35) -> Featured:
    """N = 4, five-fold, index 2, from two counter-rotating epicycles."""
    P = Action(4, 2, K)
    r = run(N=4, d=2, K=K, phase1="gradnorm", n1=200, newton=150,
            start=epicycle_start(P, 2, 3, 1.0, b))
    return Featured("pentagon", "the five-fold N = 4 orbit",
                    "index 2, from a mode 2 and a mode 3 turning opposite ways", r)


def hiphop(K: int = 32, seed: int = 0) -> Featured:
    """N = 4 out of the plane, from the unchained-polygon start at the 5:6 resonance."""
    P = Action(4, 3, K)
    r = run(N=4, d=3, K=K, phase1="gradnorm", n1=300, newton=150,
            start=vertical_start(P, 5, 6, 0.45))
    return Featured("hiphop", "the N = 4 hip-hop",
                    "modes 5 and 6, at a resonance of the rotating square", r)


# ---------------------------------------------------------------------------
# the two motions that cost nothing
# ---------------------------------------------------------------------------
def time_shift(P: Action, x, tau: float) -> np.ndarray:
    """The same loop, started at a different moment.

    q(t + tau) has the same modes with the coefficient pair turned by m*tau, so this
    is a rotation of each mode's (cos, sin) pair -- and the action does not move.
    """
    x = P.shaped(x)
    m = P.modes[:, None].astype(float)
    c, s = np.cos(m * tau), np.sin(m * tau)
    out = np.empty_like(x)
    out[:, 0, :] = c * x[:, 0, :] + s * x[:, 1, :]
    out[:, 1, :] = -s * x[:, 0, :] + c * x[:, 1, :]
    return out


def rotate(P: Action, x, angle: float, plane=(0, 1)) -> np.ndarray:
    """The same loop, turned in space. Every mutual distance is unchanged."""
    x = np.array(P.shaped(x), copy=True)
    a, b = plane
    c, s = np.cos(angle), np.sin(angle)
    xa, xb = x[:, :, a].copy(), x[:, :, b].copy()
    x[:, :, a] = c * xa - s * xb
    x[:, :, b] = s * xa + c * xb
    return x


def eight_newton(seed: int = 3, K: int = 24, n1: int = 40) -> Featured:
    """The eight again, but with phase one cut short so Newton has real work to do."""
    r = run(N=3, d=2, K=K, seed=seed, phase1="gradnorm", n1=n1, newton=60)
    return Featured("eight_newton", "the figure eight",
                    "phase two, starting from a half-finished phase one", r)
