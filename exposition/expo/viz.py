"""Drawing loops, bodies, spectra and eigenvalues."""

from __future__ import annotations

import numpy as np
from manim import (VGroup, VMobject, Dot, Circle, Line, Rectangle, DOWN, RIGHT,
                   ORIGIN, PI)

from . import theme as T
from .mathtext import C


class Projector:
    """An orthographic camera in R^d: two orthonormal rows, freely rotatable."""

    def __init__(self, d: int, frame: np.ndarray | None = None):
        self.d = d
        self.R = np.eye(d) if frame is None else np.array(frame, dtype=float)

    def matrix(self) -> np.ndarray:
        return self.R[:2]

    def rotated(self, i: int, j: int, angle: float) -> "Projector":
        """A copy of this camera turned by `angle` in the (i, j) plane."""
        G = np.eye(self.d)
        c, s = np.cos(angle), np.sin(angle)
        G[i, i] = G[j, j] = c
        G[i, j], G[j, i] = -s, s
        p = Projector(self.d)
        p.R = G @ self.R
        return p

    def __call__(self, X: np.ndarray) -> np.ndarray:
        """(..., d) -> (..., 2)."""
        return np.asarray(X) @ self.matrix().T


def fit(points2: np.ndarray, radius: float = 2.4, center=ORIGIN):
    """Scale a set of 2-D points so the widest extent spans 2*radius."""
    p = np.asarray(points2)
    span = float(np.abs(p).max()) if p.size else 1.0
    s = radius / max(span, 1e-12)
    return s, np.asarray(center, dtype=float)


def to_screen(points2: np.ndarray, scale: float, center) -> np.ndarray:
    p = np.asarray(points2) * scale
    out = np.zeros((len(p), 3))
    out[:, :2] = p
    return out + np.asarray(center, dtype=float)


def loop_curve(pts3: np.ndarray, color=None, width: float = 3.2,
               opacity: float = 1.0, gradient: bool = False) -> VMobject:
    """A closed curve through the given screen points."""
    m = VMobject(stroke_width=width)
    m.set_points_smoothly([*pts3, pts3[0]])
    if gradient:
        m.set_stroke(color=[T.CURVE, T.ROSE, T.COOL, T.CURVE], width=width, opacity=opacity)
    else:
        m.set_stroke(color=color or T.CURVE, width=width, opacity=opacity)
    m.set_fill(opacity=0)
    return m


HALO_LAYERS = 40
HALO_PEAK = 0.24
HALO_POWER = 1.7


def halo_alphas(layers: int = HALO_LAYERS, peak: float = HALO_PEAK,
                power: float = HALO_POWER):
    """Layer opacities that composite to the profile above; index 0 is the innermost"""
    t = [peak * max(0.0, 1.0 - (j - 0.5) / layers) ** power
         for j in range(1, layers + 1)] + [0.0]
    return [(t[j - 1] - t[j]) / (1.0 - t[j]) for j in range(1, layers + 1)]


def body_dot(color, r: float = 0.085, glow: bool = True,
             layers: int = HALO_LAYERS, reach: float = 3.4) -> VGroup:
    """A body: a solid dot inside a soft halo."""
    parts = []
    if glow:
        alphas = halo_alphas(layers)
        for j in range(layers, 0, -1):
            parts.append(Circle(radius=r * (1.0 + reach * j / layers), stroke_width=0,
                                fill_color=color, fill_opacity=alphas[j - 1]))
    parts.append(Circle(radius=r, stroke_width=0, fill_color=color, fill_opacity=1.0))
    return VGroup(*parts)


class OrbitView(VGroup):
    """A loop with its N bodies, ready to be animated in time."""

    def __init__(self, bodies_fn, N: int, radius: float = 2.4, center=ORIGIN,
                 projector: Projector = None, samples: int = 720,
                 curve_width: float = 3.0, gradient: bool = False,
                 show_bodies: bool = True, trail: float = 0.0,
                 curve_color=None, dot_radius: float = 0.085,
                 fixed_scale: bool = False):
        super().__init__()
        self.bodies_fn = bodies_fn
        self.N = N
        self.samples = samples
        self.radius = radius
        self.center = np.asarray(center, dtype=float)
        self.proj = projector
        self.gradient = gradient
        self.curve_width = curve_width
        self.curve_color = curve_color
        self.trail = trail
        self.dot_radius = dot_radius
        self.fixed_scale = fixed_scale

        ts = np.linspace(0, 2 * PI, samples, endpoint=False)
        self.raw = bodies_fn(ts)
        if self.proj is None:
            self.proj = Projector(self.raw.shape[2])
        self._rescale()

        self.curve = loop_curve(self._screen(self.raw[0]), color=curve_color,
                                width=curve_width, gradient=gradient)
        self.add(self.curve)
        self.dots = VGroup()
        if show_bodies:
            for k in range(N):
                self.dots.add(body_dot(T.body_color(k, N), r=dot_radius))
            self.add(self.dots)
        self.trails = VGroup()
        if trail > 0:
            for k in range(N):
                tr = VMobject(stroke_width=2.0)
                tr.set_stroke(T.body_color(k, N), opacity=0.55)
                self.trails.add(tr)
            self.add(self.trails)
        self.set_time(0.0)

    def _rescale(self):
        """Fit the drawing to `radius`."""
        X = self.raw.reshape(-1, self.raw.shape[2])
        flat = X if self.fixed_scale else self.proj(X)
        self.scale_factor, _ = fit(flat, self.radius)

    def _screen(self, X: np.ndarray) -> np.ndarray:
        return to_screen(self.proj(X), self.scale_factor, self.center)

    def redraw(self):
        """Rebuild every part from the current camera, scale, centre and time."""
        pts = self._screen(self.raw[0])
        self.curve.set_points_smoothly([*pts, pts[0]])
        self.set_time(self._t)
        return self

    def place(self, radius: float = None, center=None, rescale: bool = True):
        if radius is not None:
            self.radius = float(radius)
        if center is not None:
            self.center = np.asarray(center, dtype=float)
        if rescale:
            self._rescale()
        return self.redraw()

    def set_projector(self, p: Projector, rescale: bool = False):
        self.proj = p
        if rescale:
            self._rescale()
        return self.redraw()

    def set_coefficients(self, bodies_fn, rescale: bool = False):
        """Swap in a new loop (a new optimiser iterate) keeping the camera."""
        ts = np.linspace(0, 2 * PI, self.samples, endpoint=False)
        self.raw = bodies_fn(ts)
        if rescale:
            self._rescale()
        return self.redraw()

    def set_time(self, t: float):
        self._t = t
        j = int(round(t / (2 * PI) * self.samples)) % self.samples
        for k, dot in enumerate(self.dots):
            dot.move_to(self._screen(self.raw[k, j:j + 1])[0])
        if self.trail > 0:
            n = max(2, int(self.trail * self.samples))
            for k, tr in enumerate(self.trails):
                idx = (np.arange(j - n, j + 1)) % self.samples
                tr.set_points_smoothly(self._screen(self.raw[k][idx]))
        return self


def spectrum_bars(modes, power, width: float = 5.6, height: float = 1.5,
                  color=None, log: bool = True, floor: float = 1e-10,
                  label_every: int = 0) -> VGroup:
    """Per-mode power |c_m|^2 + |s_m|^2, the loop's Fourier fingerprint."""
    power = np.asarray(power, dtype=float)
    if log:
        v = np.log10(np.maximum(power, floor))
        v = (v - np.log10(floor)) / (-np.log10(floor) + np.log10(max(power.max(), floor)) + 1e-12)
        v = np.clip(v, 0, 1)
    else:
        v = power / max(power.max(), 1e-300)
    g = VGroup()
    n = len(modes)
    bw = width / max(n, 1)
    for i, (m, h) in enumerate(zip(modes, v)):
        bar = Rectangle(width=bw * 0.66, height=max(h * height, 0.006),
                        stroke_width=0, fill_color=color or T.COOL, fill_opacity=0.92)
        bar.move_to([-width / 2 + (i + 0.5) * bw, max(h * height, 0.006) / 2, 0])
        g.add(bar)
    axis = Line([-width / 2, 0, 0], [width / 2, 0, 0], stroke_width=1.4, color=T.RULE)
    bars = VGroup(*g)
    g.add(axis)
    if label_every:
        for i, m in enumerate(modes):
            if i % label_every:
                continue
            t = C(str(int(m)), size=T.SZ_TINY, color=T.RULE)
            t.move_to([-width / 2 + (i + 0.5) * bw, -0.22, 0])
            g.add(t)
    g.bars = bars
    g.axis = axis
    return g


def eigen_strip(w, width: float = 6.4, height: float = 0.5, rel_tol: float = 1e-8,
                label: bool = True, window: int = 0) -> VGroup:
    """The Hessian spectrum as a signed strip: red below zero, green above, violet"""
    w = np.asarray(w, dtype=float)
    lmax = float(np.abs(w).max())
    ws = np.sort(w)
    trimmed = 0
    if window and window < len(ws):
        trimmed = len(ws) - window
        ws = ws[:window]
    n = len(ws)
    g = VGroup()
    bw = width / n
    neg = zero = 0
    for i, lam in enumerate(ws):
        if abs(lam) <= rel_tol * lmax:
            col, zero = T.ZERO, zero + 1
        elif lam < 0:
            col, neg = T.BAD, neg + 1
        else:
            col = T.GOOD
        r = Rectangle(width=bw * 0.86, height=height, stroke_width=0,
                      fill_color=col, fill_opacity=0.92)
        r.move_to([-width / 2 + (i + 0.5) * bw, 0, 0])
        g.add(r)
    g.neg, g.zero, g.pos = neg, zero, n - neg - zero
    g.trimmed = trimmed
    if label:
        txt = "Morse index %d   nullity %d" % (neg, zero)
        if trimmed:
            txt += "    (%d larger positive values not shown)" % trimmed
        tag = C(txt, size=T.SZ_TINY, color=T.INK_DIM)
        tag.next_to(g, DOWN, buff=0.18)
        g.add(tag)
    return g


def value_readout(label: str, value: str, color=None, size=None) -> VGroup:
    """A small `name  value` pair for live numbers."""
    lab = C(label, size=size or T.SZ_SMALL, color=T.INK_DIM)
    val = C(value, size=size or T.SZ_SMALL, color=color or T.INK)
    g = VGroup(lab, val).arrange(RIGHT, buff=0.22)
    g.label, g.value = lab, val
    return g


def axes_box(width: float, height: float, color=None) -> VGroup:
    """A plain framed plotting area -- two rules, no chrome."""
    x = Line([0, 0, 0], [width, 0, 0], stroke_width=1.4, color=color or T.RULE)
    y = Line([0, 0, 0], [0, height, 0], stroke_width=1.4, color=color or T.RULE)
    g = VGroup(x, y)
    g.x_axis, g.y_axis = x, y
    return g


def log_plot(values, width: float = 5.2, height: float = 2.4, lo: float = None,
             hi: float = None, color=None, width_stroke: float = 2.6):
    """A |grad| history on a log scale, returned with the mapping so points can be"""
    v = np.maximum(np.asarray(values, dtype=float), 1e-300)
    lo = np.log10(lo if lo else max(v.min() * 0.5, 1e-16))
    hi = np.log10(hi if hi else v.max() * 2)

    def xy(i, val, n):
        x = width * (i / max(n - 1, 1))
        y = height * (np.log10(max(val, 1e-300)) - lo) / max(hi - lo, 1e-9)
        return np.array([x, np.clip(y, 0, height), 0.0])

    pts = [xy(i, val, len(v)) for i, val in enumerate(v)]
    m = VMobject(stroke_width=width_stroke)
    m.set_points_smoothly(pts) if len(pts) > 2 else m.set_points_as_corners(pts or [ORIGIN])
    m.set_stroke(color or T.COOL)
    m.set_fill(opacity=0)
    m.xy = xy
    m.lo, m.hi = lo, hi
    return m


def spin(view: OrbitView, turns: float = 1.0, t0: float = None):
    """An animation that runs the orbit clock forward."""
    from manim import UpdateFromAlphaFunc
    t0 = view._t if t0 is None else t0
    return UpdateFromAlphaFunc(view, lambda m, a: m.set_time(t0 + a * 2 * PI * turns))


def spin_and_place(view: OrbitView, radius: float = None, center=None,
                   turns: float = 1.0, t0: float = None):
    """Move and rescale the view while the orbit keeps running."""
    from manim import UpdateFromAlphaFunc
    t0 = view._t if t0 is None else t0
    r0 = view.radius
    c0 = np.array(view.center, dtype=float)
    r1 = r0 if radius is None else float(radius)
    c1 = c0 if center is None else np.asarray(center, dtype=float)

    def f(m, a):
        m.radius = r0 + (r1 - r0) * a
        m.center = c0 + (c1 - c0) * a
        m._rescale()
        m._t = t0 + a * 2 * PI * turns
        m.redraw()

    return UpdateFromAlphaFunc(view, f)


def default_tumble_planes(d: int):
    """Which planes to turn the camera in."""
    if d < 4:
        return [(0, 2)] if d >= 3 else []
    planes = []
    for a in range(0, d - 3, 4):
        planes += [(a, a + 3), (a + 1, a + 2)]
    planes = [(u, v) for u, v in planes if v < d]
    return planes[:4] or [(0, d - 1)]


def tumble(view: OrbitView, turns: float = 1.0, sweep: float = 2 * PI,
           planes=None, t0: float = None, base: Projector = None):
    """Turn the camera through the orbit's own principal planes while it runs."""
    from manim import UpdateFromAlphaFunc
    base = base or view.proj
    t0 = view._t if t0 is None else t0
    d = base.d
    if planes is None:
        planes = default_tumble_planes(d)

    def f(m, a):
        p = Projector(d, base.R)
        for i, (u, v) in enumerate(planes):
            p = p.rotated(u, v, sweep * a / (i + 1))
        m.proj = p
        m._t = t0 + a * 2 * PI * turns
        m.redraw()

    return UpdateFromAlphaFunc(view, f)


class PlaneGrid(VGroup):
    """One small panel per principal plane -- the gallery's own way of drawing a"""

    def __init__(self, orbit, cols: int = 3, panel: float = 1.55, gap: float = 0.42,
                 center=ORIGIN, samples: int = 480, show_bodies: bool = False,
                 label: bool = True):
        super().__init__()
        d = orbit.d
        ts = np.linspace(0, 2 * PI, samples, endpoint=False)
        F = orbit.principal_frame()
        X = np.einsum("ab,ktb->kta", F, orbit.bodies(ts))
        sv = orbit.principal_values()
        scale = panel * 0.44 / max(float(np.abs(X).max()), 1e-12)
        npl = (d + 1) // 2
        rows = (npl + cols - 1) // cols
        for i in range(npl):
            a, b = 2 * i, min(2 * i + 1, d - 1)
            r, c = divmod(i, cols)
            off = np.array([(c - (cols - 1) / 2) * (panel + gap),
                            ((rows - 1) / 2 - r) * (panel + gap + 0.28), 0.0])
            off = off + np.asarray(center, dtype=float)
            box = Rectangle(width=panel, height=panel, stroke_width=1.0,
                            color=T.RULE, fill_opacity=0)
            box.move_to(off)
            self.add(box)
            thin = (b == a) or sv[b] < 1e-4 * sv[0]
            pts = np.zeros((samples, 3))
            pts[:, 0] = X[0, :, a] * scale
            pts[:, 1] = (0.0 if thin else X[0, :, b] * scale)
            pts += off
            col = T.CURVE if not thin else T.GHOST
            m = VMobject(stroke_width=2.2)
            m.set_points_smoothly([*pts, pts[0]])
            m.set_stroke(col)
            m.set_fill(opacity=0)
            self.add(m)
            if show_bodies:
                for k in range(orbit.N):
                    dot = Dot(radius=0.035, color=T.body_color(k, orbit.N))
                    dot.move_to(off + np.array([X[k, 0, a] * scale,
                                                0.0 if thin else X[k, 0, b] * scale, 0]))
                    self.add(dot)
            if label:
                name = "%d" % (a + 1) if b == a else "%d,%d" % (a + 1, b + 1)
                t = C(name, size=T.SZ_TINY, color=T.RULE)
                t.next_to(box, DOWN, buff=0.10)
                self.add(t)
