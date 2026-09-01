"""Numbers that change every frame.

Building a text mobject costs a Pango layout and an SVG parse, which is fine once
and ruinous sixty times a second at four times the frame rate. So the glyphs are
built once into an atlas and each character slot is re-pointed at the glyph it
needs; the monospaced face makes the layout a fixed pitch, so nothing has to be
measured again either.
"""

from __future__ import annotations

import numpy as np
from manim import VGroup, VMobject, LEFT, RIGHT

from . import theme as T
from .mathtext import C

ALPHABET = "0123456789.+-e ×"


class Readout(VGroup):
    """A fixed-width field of `width` characters, updatable with `set`."""

    def __init__(self, width: int, size: float = None, color=None, align=LEFT):
        super().__init__()
        size = size or T.SZ_SMALL
        self.color_ = color or T.INK
        # Each glyph is measured next to a reference "0" so it keeps that glyph's
        # baseline and advance: centring a full stop in its own box would float it
        # to the middle of the line, which is where a decimal point must not be.
        self.atlas = {}
        ref = C("00", size=size, color=self.color_)
        self.pitch = float(ref[1].get_left()[0] - ref[0].get_left()[0])
        for ch in ALPHABET:
            if ch == " ":
                self.atlas[ch] = None
                continue
            t = C("0" + ch, size=size, color=self.color_)
            if len(t.submobjects) < 2:
                self.atlas[ch] = None
                continue
            anchor = np.array([t[0].get_left()[0], t[0].get_bottom()[1], 0.0])
            g = VGroup(*[m.copy() for m in t.submobjects[1:]])
            g.shift(-anchor)
            self.atlas[ch] = g
        self.slots = []
        for i in range(width):
            m = VMobject()
            self.slots.append(m)
            self.add(m)
        self.width_chars = width
        self.align = align
        self._origin = np.zeros(3)
        self.set("")

    def anchor(self, point, align=None):
        """`point` is the left end of the baseline."""
        self._origin = np.asarray(point, dtype=float)
        if align is not None:
            self.align = align
        return self.set(self._text)

    def set(self, text: str):
        self._text = text
        n = self.width_chars
        text = text[:n]
        pad = n - len(text)
        if self.align is RIGHT:
            text = " " * pad + text
        else:
            text = text + " " * pad
        base = self._origin
        if self.align is RIGHT:
            base = base - np.array([self.pitch * n, 0, 0])
        for i, ch in enumerate(text):
            g = self.atlas.get(ch)
            slot = self.slots[i]
            if g is None:
                slot.become(VMobject())
                continue
            c = g.copy()
            c.shift(base + np.array([self.pitch * i, 0, 0]))
            slot.become(c)
        return self


class Field(VGroup):
    """A label with a live value beside it."""

    def __init__(self, label: str, width: int = 12, size: float = None,
                 color=None, label_color=None):
        super().__init__()
        self.label = C(label, size=size or T.SZ_SMALL, color=label_color or T.INK_DIM)
        self.value = Readout(width, size=size, color=color)
        self.add(self.label, self.value)
        self.value.anchor(self.label.get_right() + np.array([0.22, 0, 0]))

    def set(self, text: str):
        self.value.set(text)
        return self

    def place(self, point, align=LEFT):
        self.label.move_to(np.asarray(point, dtype=float), aligned_edge=align)
        self.value.anchor(self.label.get_right() + np.array([0.22, 0, 0]))
        return self


def sci(v: float, digits: int = 2) -> str:
    """1.2e-07 rather than 1.20000e-07: the exponent is the news."""
    if v == 0 or not np.isfinite(v):
        return "0"
    e = int(np.floor(np.log10(abs(v))))
    m = v / 10.0 ** e
    return ("%." + str(digits) + "fe%+03d") % (m, e)
