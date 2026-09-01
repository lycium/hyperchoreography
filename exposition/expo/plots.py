"""A small plotting kit.

Manim's own Axes label through LaTeX, which is not in this render path, and they
carry more chrome than these pictures want. This is the minimum: a frame, ticks,
and a way to turn data coordinates into points so a curve can be drawn, grown, or
have a marker walked along it while the rest of the scene animates.
"""

from __future__ import annotations

import numpy as np
from manim import (VGroup, VMobject, Line, Dot, DashedLine, Rectangle, UP, DOWN,
                   LEFT, RIGHT, ORIGIN)

from . import theme as T
from .mathtext import C


class Plot(VGroup):
    """Axes over [x0, x1] x [y0, y1] with the origin at `center` minus half the size."""

    def __init__(self, x_range, y_range, width: float = 5.4, height: float = 3.0,
                 center=ORIGIN, log_y: bool = False, x_ticks=None, y_ticks=None,
                 x_label: str = None, y_label: str = None, frame: bool = False,
                 tick_fmt=None, color=None):
        super().__init__()
        self.w, self.h = width, height
        self.log_y = log_y
        self.x0, self.x1 = float(x_range[0]), float(x_range[1])
        y0, y1 = float(y_range[0]), float(y_range[1])
        self.y0, self.y1 = (np.log10(y0), np.log10(y1)) if log_y else (y0, y1)
        self.origin = np.asarray(center, dtype=float) + np.array([-width / 2, -height / 2, 0])
        col = color or T.RULE

        if frame:
            box = Rectangle(width=width, height=height, stroke_width=1.1, color=col,
                            fill_opacity=0).move_to(np.asarray(center, dtype=float))
            self.add(box)
        else:
            self.add(Line(self.origin, self.origin + [width, 0, 0],
                          stroke_width=1.3, color=col))
            self.add(Line(self.origin, self.origin + [0, height, 0],
                          stroke_width=1.3, color=col))

        fmt = tick_fmt or (lambda v: ("%g" % v))
        for v in (x_ticks or []):
            p = self.at(v, None)
            self.add(Line(p + DOWN * 0.06, p + UP * 0.06, stroke_width=1.2, color=col))
            t = C(fmt(v), size=T.SZ_TINY, color=T.INK_DIM)
            t.next_to(p, DOWN, buff=0.14)
            self.add(t)
        for v in (y_ticks or []):
            p = self.at(None, v)
            self.add(Line(p + LEFT * 0.06, p + RIGHT * 0.06, stroke_width=1.2, color=col))
            t = C(self.y_tick_text(v), size=T.SZ_TINY, color=T.INK_DIM)
            t.next_to(p, LEFT, buff=0.16)
            self.add(t)
        if x_label:
            t = C(x_label, size=T.SZ_TINY, color=T.INK_DIM)
            t.move_to(self.origin + [width / 2, -0.55, 0])
            self.add(t)
        if y_label:
            t = C(y_label, size=T.SZ_TINY, color=T.INK_DIM)
            t.move_to(self.origin + [0, height + 0.32, 0])
            t.align_to(self.origin + [0, 0, 0], LEFT)
            self.add(t)

    # -- coordinates ------------------------------------------------------
    def y_tick_text(self, v):
        if not self.log_y:
            return "%g" % v
        e = int(round(np.log10(v)))
        return "1e%d" % e

    def at(self, x=None, y=None) -> np.ndarray:
        p = np.array(self.origin, dtype=float)
        if x is not None:
            p = p + np.array([self.w * (x - self.x0) / max(self.x1 - self.x0, 1e-30), 0, 0])
        if y is not None:
            yv = np.log10(max(float(y), 1e-300)) if self.log_y else float(y)
            p = p + np.array([0, self.h * (yv - self.y0) / max(self.y1 - self.y0, 1e-30), 0])
        return p

    def clamp(self, p) -> np.ndarray:
        p = np.array(p, dtype=float)
        p[0] = np.clip(p[0], self.origin[0], self.origin[0] + self.w)
        p[1] = np.clip(p[1], self.origin[1], self.origin[1] + self.h)
        return p

    # -- content ----------------------------------------------------------
    def line(self, xs, ys, color=None, width: float = 2.6, smooth: bool = True,
             opacity: float = 1.0) -> VMobject:
        pts = [self.clamp(self.at(x, y)) for x, y in zip(xs, ys)]
        m = VMobject(stroke_width=width)
        if len(pts) > 2 and smooth:
            m.set_points_smoothly(pts)
        else:
            m.set_points_as_corners(pts if len(pts) > 1 else pts * 2)
        m.set_stroke(color or T.COOL, opacity=opacity)
        m.set_fill(opacity=0)
        return m

    def marker(self, x, y, color=None, r: float = 0.055) -> Dot:
        return Dot(self.clamp(self.at(x, y)), radius=r, color=color or T.CURVE)

    def hline(self, y, color=None, dashed: bool = True, opacity: float = 0.6):
        a, b = self.at(self.x0, y), self.at(self.x1, y)
        cls = DashedLine if dashed else Line
        return cls(a, b, stroke_width=1.2, color=color or T.RULE,
                   stroke_opacity=opacity)

    def vline(self, x, color=None, dashed: bool = True, opacity: float = 0.6):
        a, b = self.at(x, self.y0 if not self.log_y else 10 ** self.y0), \
               self.at(x, self.y1 if not self.log_y else 10 ** self.y1)
        cls = DashedLine if dashed else Line
        return cls(a, b, stroke_width=1.2, color=color or T.RULE,
                   stroke_opacity=opacity)


def growing_line(plot: Plot, xs, ys, color=None, width: float = 2.6):
    """A curve plus an updater that reveals it as an animation runs.

    Returns (mobject, setter); `setter(k)` shows the first k points, so the plot can
    be drawn in step with whatever else is happening.
    """
    xs = list(xs)
    ys = list(ys)
    m = VMobject(stroke_width=width)
    m.set_stroke(color or T.COOL)
    m.set_fill(opacity=0)

    def setter(k: int):
        k = int(max(2, min(len(xs), k)))
        pts = [plot.clamp(plot.at(x, y)) for x, y in zip(xs[:k], ys[:k])]
        m.set_points_as_corners(pts)
        return m

    setter(2)
    return m, setter
