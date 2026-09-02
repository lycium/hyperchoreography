"""Typesetting the formulas.

There is no LaTeX in the render path, so display maths is built from Pango text:
Unicode symbols in Palatino, with a small layout kit for the pieces Unicode cannot
do on its own -- real fractions, and big operators with limits above and below.

`M()` is the workhorse and takes Pango markup, so <i>q</i> italicises a variable
and <sup>2</sup> raises an exponent. Everything is a VGroup, so any part can be
coloured, faded or transformed on its own.
"""

from __future__ import annotations

from manim import (VGroup, Line, MarkupText, UP, DOWN, LEFT, RIGHT, ORIGIN)

from . import theme as T


# Two Pango facts to work around. Glyphs are laid out at whole pixels, so text asked
# for at 15 px comes back with visibly uneven letter spacing -- everything is
# therefore built at one large size and scaled to fit, which costs nothing since the
# glyphs are vector paths. But manim hands Pango a fixed layout width, so a long line
# built too large silently wraps; `_text` notices that and steps the build size down
# until the line is a single line again.
_BASE = 44.0
_MIN_BASE = 20.0


def _line_height(font: str, base: float) -> float:
    key = (font, base)
    if key not in _line_height.cache:
        _line_height.cache[key] = MarkupText("Hxg", font=font, font_size=base).height
    return _line_height.cache[key]


_line_height.cache = {}


def _text(markup: str, font: str, size: float, color, **kw) -> MarkupText:
    base = _BASE
    while True:
        m = MarkupText(markup, font=font, font_size=base, color=color, **kw)
        wrapped = ("\n" not in markup
                   and m.height > 1.6 * _line_height(font, base))
        if not wrapped or base <= _MIN_BASE:
            m.scale(size * 48.0 / base)
            return m
        base = max(_MIN_BASE, base * 0.7)


def M(markup: str, size: float = None, color=None, font: str = None, **kw):
    """One line of maths. Pango markup; the font defaults to the maths face."""
    return _text(markup, font or T.FONT_MATH, size or T.SZ_MATH, color or T.INK, **kw)


def B(markup: str, size: float = None, color=None, **kw):
    """One line of body text, in the narration face."""
    return _text(markup, T.FONT_BODY, size or T.SZ_BODY, color or T.INK, **kw)


def C(markup: str, size: float = None, color=None, **kw):
    """A command line or a record field, in the monospaced face."""
    return _text(markup, T.FONT_CODE, size or T.SZ_SMALL, color or T.INK_DIM, **kw)


def Crow(markup: str, size: float = None, color=None, **kw):
    """A row of a monospaced table. Rows are separate objects (so they can be coloured and
    revealed one by one), and manim aligns objects by their ink: two rows whose first glyphs
    have different side bearings land a pixel or two apart, and their columns with them. A
    hidden leading glyph gives every row the same bearing, so left-aligned rows share one
    character grid exactly."""
    m = C("|" + markup, size, color, **kw)
    m[0].set_color(T.BG).set_opacity(0)          # invisible whichever of the two a later call resets
    return m


def var(name: str) -> str:
    """A variable: italic, as it would be set in print."""
    return "<i>%s</i>" % name


def sb(s) -> str:
    """A subscript. Pango's own is set rather small for a formula, so it is
    enlarged inside the subscript baseline where it still sits correctly."""
    return '<sub><span size="152%%">%s</span></sub>' % s


def sp(s) -> str:
    """A superscript, enlarged the same way."""
    return '<sup><span size="152%%">%s</span></sup>' % s


def frac(num, den, width_pad: float = 0.12, rule_color=None) -> VGroup:
    """A real fraction: numerator, rule, denominator, centred on the rule."""
    num = M(num) if isinstance(num, str) else num
    den = M(den) if isinstance(den, str) else den
    w = max(num.width, den.width) + width_pad
    bar = Line(LEFT * w / 2, RIGHT * w / 2,
               stroke_width=2.2, color=rule_color or T.INK)
    g = VGroup(num, bar, den)
    num.next_to(bar, UP, buff=0.10)
    den.next_to(bar, DOWN, buff=0.12)
    g.rule = bar
    return g


def bigop(symbol: str, body, sub: str = None, sup: str = None,
          op_scale: float = 2.1, limit_scale: float = 0.58,
          limits_beside: bool = False, color=None) -> VGroup:
    """A big operator with its limits: sum, integral, product.

    `limits_beside` puts them to the right of the sign, the way an integral is
    usually set; the default stacks them, the way a sum is.
    """
    sign = M(symbol, size=T.SZ_MATH * op_scale, color=color)
    parts = [sign]
    lo = M(sub, size=T.SZ_MATH * limit_scale, color=color) if sub else None
    hi = M(sup, size=T.SZ_MATH * limit_scale, color=color) if sup else None
    if limits_beside:
        if hi:
            hi.next_to(sign, RIGHT, buff=0.04).align_to(sign, UP)
        if lo:
            lo.next_to(sign, RIGHT, buff=0.04).align_to(sign, DOWN)
    else:
        if hi:
            hi.next_to(sign, UP, buff=0.06)
        if lo:
            lo.next_to(sign, DOWN, buff=0.06)
    for p in (lo, hi):
        if p is not None:
            parts.append(p)
    body = M(body) if isinstance(body, str) else body
    body.next_to(sign, RIGHT, buff=0.14)
    parts.append(body)
    g = VGroup(*parts)
    g.sign, g.lo, g.hi, g.body = sign, lo, hi, body
    return g


def row(*items, buff: float = 0.16, align=None) -> VGroup:
    """Lay pieces out left to right on a common baseline-ish centre."""
    mobs = [M(i) if isinstance(i, str) else i for i in items]
    g = VGroup(*mobs)
    g.arrange(RIGHT, buff=buff, aligned_edge=align if align is not None else ORIGIN)
    return g


def stack(*items, buff: float = 0.30, align=LEFT) -> VGroup:
    mobs = [B(i) if isinstance(i, str) else i for i in items]
    g = VGroup(*mobs)
    g.arrange(DOWN, buff=buff, aligned_edge=align)
    return g


# ---------------------------------------------------------------------------
# the symbols the deck uses, in one place so they stay consistent
# ---------------------------------------------------------------------------
QDOT = "q̇"          # q with a dot above
ALPHA, BETA, GAMMA = "α", "β", "γ"
DELTA, LAMBDA, MU = "δ", "λ", "μ"
OMEGA, omega = "Ω", "ω"
CHI, PSI, PHI = "χ", "ψ", "φ"
SIGMA, THETA, EPS = "σ", "θ", "ε"
NABLA = "∇"
INT = "∫"
SUM = "∑"
OINT = "∮"
WEDGE = "∧"
APPROX, NEQ, LEQ, GEQ = "≈", "≠", "≤", "≥"
IN = "∈"
TO = "→"
TIMES = "×"
HALF = "½"
FLOOR_L, FLOOR_R = "⌊", "⌋"
LANGLE, RANGLE = "⟨", "⟩"
RR = "ℝ"
INFTY = "∞"
DOT = "·"
