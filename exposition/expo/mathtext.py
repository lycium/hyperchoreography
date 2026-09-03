"""Typesetting the formulas."""

from __future__ import annotations

from manim import (VGroup, Line, MarkupText, UP, DOWN, LEFT, RIGHT, ORIGIN)

from . import theme as T


_BASE = 44.0
_MIN_BASE = 20.0

SUP_SIZE, SUP_RISE = 62, 0.60
SUB_SIZE, SUB_DROP = 62, 0.30

# Pango's rise is an absolute displacement, 1/20480 of a manim unit whatever the
# font size, so it can only be filled in once the final size is known.
_RISE_PER_UNIT = 20480.0
_SUP_TOKEN, _SUB_TOKEN = "@@sup-rise@@", "@@sub-rise@@"


def _measure(markup: str, font: str, base: float) -> float:
    key = (markup, font, base)
    if key not in _measure.cache:
        _measure.cache[key] = MarkupText(markup, font=font, font_size=base).height
    return _measure.cache[key]


_measure.cache = {}


def _line_height(font: str, base: float) -> float:
    return _measure("Hxg", font, base)


def _resolve_rise(markup: str, font: str, base: float) -> str:
    """Fill in the sub/superscript displacements for the size this line is built at."""
    if _SUP_TOKEN not in markup and _SUB_TOKEN not in markup:
        return markup
    cap = _measure("H", font, base) * _RISE_PER_UNIT
    return (markup.replace(_SUP_TOKEN, "%d" % round(SUP_RISE * cap))
                  .replace(_SUB_TOKEN, "%d" % -round(SUB_DROP * cap)))


def _text(markup: str, font: str, size: float, color, **kw) -> MarkupText:
    base = _BASE
    while True:
        m = MarkupText(_resolve_rise(markup, font, base), font=font, font_size=base,
                       color=color, **kw)
        wrapped = ("\n" not in markup
                   and m.height > 1.75 * _line_height(font, base))
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
    """A row of a monospaced table. Rows are separate objects (so they can be coloured and"""
    m = C("|" + markup, size, color, **kw)
    m[0].set_color(T.BG).set_opacity(0)
    return m


_BAR_OVER = 1.15


def _cell(text: str, size, color, maker, bar: float):
    size = size or T.SZ_SMALL
    m = maker('<span size="%d%%">|</span>%s' % (round(100.0 * bar / size), text),
              size, color)
    m[0].set_color(T.BG).set_opacity(0)
    return m


def table(rows, align: str = None, sizes=None, colors=None, maker=None,
          col_colors=None, col_sizes=None,
          col_buff: float = 0.40, row_buff: float = 0.26,
          head_buff: float = None) -> VGroup:
    """A table whose columns line up whatever size each row is set at."""
    maker = maker or C
    rows = [tuple(r) for r in rows]
    ncol = max(len(r) for r in rows)
    align = (align or "l" * ncol).ljust(ncol, "l")

    def per(v, i):
        return v[i] if isinstance(v, (list, tuple)) else v

    def spec(by_col, by_row, i, j):
        v = per(by_col, j)
        return per(by_row, i) if v is None else v

    at = [[spec(col_sizes, sizes, i, j) or T.SZ_SMALL for j in range(len(r))]
          for i, r in enumerate(rows)]
    bar = _BAR_OVER * max(max(row) for row in at)

    grid = [[_cell(t, at[i][j], spec(col_colors, colors, i, j), per(maker, j), bar)
             if t != "" else None
             for j, t in enumerate(r)] + [None] * (ncol - len(r))
            for i, r in enumerate(rows)]

    widths = [max((c.width for c in col if c is not None), default=0.0)
              for col in zip(*grid)]
    x, left = [], 0.0
    for w in widths:
        x.append(left)
        left += w + col_buff
    pitch = max(c.height for line in grid for c in line if c is not None) + row_buff

    out = VGroup()
    y = 0.0
    for i, line in enumerate(grid):
        g = VGroup()
        for j, c in enumerate(line):
            if c is None:
                continue
            cx = (x[j] + c.width / 2 if align[j] == "l" else
                  x[j] + widths[j] - c.width / 2 if align[j] == "r" else
                  x[j] + widths[j] / 2)
            c.move_to([cx, y, 0])
            g.add(c)
        out.add(g)
        y -= pitch
        if i == 0 and head_buff:
            y -= head_buff
    out.rows = list(out.submobjects)
    out.move_to(ORIGIN)
    return out


def var(name: str) -> str:
    """A variable: italic, as it would be set in print."""
    return "<i>%s</i>" % name


def sb(s) -> str:
    """A subscript: an index, dropped clear of the baseline and set in the same face."""
    return '<span size="%d%%" rise="%s">%s</span>' % (SUB_SIZE, _SUB_TOKEN, s)


def sp(s) -> str:
    """A superscript: a power, raised clear of the x-height."""
    return '<span size="%d%%" rise="%s">%s</span>' % (SUP_SIZE, _SUP_TOKEN, s)


def idx(s) -> str:
    """An index in a subscript: italic, the way a variable is set everywhere else."""
    return sb(var(s))


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
    """A big operator with its limits: sum, integral, product."""
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


QDOT = "q̇"
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
