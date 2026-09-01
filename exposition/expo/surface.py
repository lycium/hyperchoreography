"""A wireframe of the action over a two-dimensional slice of coefficient space.

The action lives on a space of a hundred-odd dimensions, so any picture of it is a
slice; these slices are real values of A, evaluated on a grid, not a sketch of what
a saddle looks like. The camera is a plain orthographic projection so it can be
turned smoothly while the rest of a scene runs.
"""

from __future__ import annotations

import numpy as np
from manim import VGroup, VMobject, ORIGIN, interpolate_color

from . import theme as T


def camera_matrix(elev: float, azim: float) -> np.ndarray:
    """Rows: right, up. A z-up world seen from `elev` above and `azim` around."""
    ce, se = np.cos(elev), np.sin(elev)
    ca, sa = np.cos(azim), np.sin(azim)
    right = np.array([-sa, ca, 0.0])
    up = np.array([-ca * se, -sa * se, ce])
    return np.stack([right, up])


class Wireframe(VGroup):
    """A grid of lines over Z(a, b), coloured by height."""

    def __init__(self, a, b, Z, width: float = 5.2, depth: float = 5.2,
                 height: float = 2.0, elev: float = 0.62, azim: float = -0.9,
                 center=ORIGIN, lines: int = 21, stroke: float = 1.7,
                 low=None, high=None, clip: float = None):
        super().__init__()
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        Z = np.asarray(Z, dtype=float)
        # Near a collision the action is infinite; a slice that clips a few corners
        # is honest and a slice that does not is unplottable.
        finite = Z[np.isfinite(Z)]
        if clip is None:
            lo = float(finite.min())
            clip = lo + 6.0 * (float(np.percentile(finite, 90)) - lo + 1e-12)
        self.Z = np.minimum(np.nan_to_num(Z, nan=clip, posinf=clip), clip)
        self.w, self.dp, self.h = width, depth, height
        self.center = np.asarray(center, dtype=float)
        self.low = low or T.COOL_DEEP
        self.high = high or T.CURVE

        z0, z1 = float(np.nanmin(self.Z)), float(np.nanmax(self.Z))
        self.z0, self.z1 = z0, z1 + (1e-12 if z1 == z0 else 0.0)

        self.rows = VGroup()
        self.cols = VGroup()
        ny, nx = self.Z.shape
        ri = np.unique(np.linspace(0, ny - 1, min(lines, ny)).astype(int))
        ci = np.unique(np.linspace(0, nx - 1, min(lines, nx)).astype(int))
        self._ri, self._ci = ri, ci
        for _ in ri:
            m = VMobject(stroke_width=stroke)
            m.set_fill(opacity=0)
            self.rows.add(m)
        for _ in ci:
            m = VMobject(stroke_width=stroke)
            m.set_fill(opacity=0)
            self.cols.add(m)
        self.add(self.rows, self.cols)
        # The world coordinates never change; only the camera does. Precomputing them
        # turns each frame of a turning shot into one matrix multiply per strand
        # instead of a Python loop over every vertex.
        A, B = self.a, self.b
        self._rows_w = [np.array([self.world(A[i], B[j], self.Z[j, i])
                                  for i in range(len(A))]) for j in ri]
        self._cols_w = [np.array([self.world(A[i], B[j], self.Z[j, i])
                                  for j in range(len(B))]) for i in ci]
        self._rows_c = [self._ramp(self.Z[j, :]) for j in ri]
        self._cols_c = [self._ramp(self.Z[:, i]) for i in ci]
        self.set_camera(elev, azim)

    # -- geometry ---------------------------------------------------------
    def world(self, ai: float, bi: float, z: float) -> np.ndarray:
        """Data coordinates to the wireframe's own 3-D box."""
        x = self.w * (ai - self.a[0]) / max(self.a[-1] - self.a[0], 1e-30) - self.w / 2
        y = self.dp * (bi - self.b[0]) / max(self.b[-1] - self.b[0], 1e-30) - self.dp / 2
        zz = self.h * (min(z, self.z1) - self.z0) / max(self.z1 - self.z0, 1e-30)
        return np.array([x, y, zz])

    def project(self, p3) -> np.ndarray:
        p = np.atleast_2d(np.asarray(p3, dtype=float))
        out = np.zeros((len(p), 3))
        out[:, :2] = p @ self.M.T
        return out + self.center

    def point(self, ai: float, bi: float, z: float) -> np.ndarray:
        return self.project(self.world(ai, bi, z))[0]

    def color_at(self, z: float):
        f = float(np.clip((z - self.z0) / max(self.z1 - self.z0, 1e-30), 0, 1))
        return interpolate_color(self.low, self.high, f ** 0.65)

    # -- drawing ----------------------------------------------------------
    def set_camera(self, elev: float, azim: float):
        self.elev, self.azim = elev, azim
        self.M = camera_matrix(elev, azim)
        for k, W in enumerate(self._rows_w):
            self.rows[k].set_points_smoothly(self.project(W))
            self.rows[k].set_stroke(self._rows_c[k], opacity=0.9)
        for k, W in enumerate(self._cols_w):
            self.cols[k].set_points_smoothly(self.project(W))
            self.cols[k].set_stroke(self._cols_c[k], opacity=0.55)
        return self

    def _ramp(self, zs, n: int = 10):
        """Colours sampled along a line, so each strand shades with its own height."""
        idx = np.linspace(0, len(zs) - 1, n).astype(int)
        return [self.color_at(float(zs[i])) for i in idx]


def turn(wire: Wireframe, d_azim: float = 0.9, d_elev: float = 0.0, extra=()):
    """Rotate the camera; `extra` are (mobject, (a, b, z)) markers to carry along."""
    from manim import UpdateFromAlphaFunc
    e0, z0 = wire.elev, wire.azim
    items = list(extra)

    def f(m, alpha):
        m.set_camera(e0 + d_elev * alpha, z0 + d_azim * alpha)
        for mob, (ai, bi, zi) in items:
            mob.move_to(m.point(ai, bi, zi))

    return UpdateFromAlphaFunc(wire, f)


def scaled_directions(H, x=None, kind=("neg", "pos"), frac: float = 0.22,
                      tol_rel: float = 1e-6):
    """Hessian eigenvectors scaled so a slice through them is readable.

    Two scalings, both necessary. Dividing by sqrt|lambda| makes the action change
    by about the same amount along each direction -- without it the stiffest
    direction of this Hessian is a thousand times stiffer than the softest, and a
    slice spanning both is a cliff beside a plain. Then the extent is capped at a
    fraction of the loop's own size, because far enough along any direction the loop
    runs into a collision and the action goes to infinity, which is true but not
    what the picture is about.

    Returns one (direction, eigenvalue, half-width) per entry of `kind`, where
    "neg" takes the steepest downhill direction and "pos" the gentlest uphill one.
    """
    w, V = np.linalg.eigh(H)
    lmax = float(np.abs(w).max())
    idx = [i for i in range(len(w)) if abs(w[i]) > tol_rel * lmax]
    neg = [i for i in idx if w[i] < 0]                    # most negative first
    pos = [i for i in idx if w[i] > 0]                    # gentlest uphill first
    xn = float(np.linalg.norm(x)) if x is not None else None
    used, out = set(), []
    for k in kind:
        pool = neg if k == "neg" else pos
        i = next((p for p in pool if p not in used), None)
        if i is None:
            i = next(p for p in idx if p not in used)
        used.add(i)
        u = V[:, i] / np.sqrt(abs(w[i]))
        span = 1.7 if xn is None else frac * xn / float(np.linalg.norm(u))
        out.append((u, float(w[i]), float(span)))
    return out
