"""Palette, type and spacing for the whole presentation.

One place to change the look. Colours are chosen so that the loop, the bodies and
the diagnostics never fight each other: the curve is warm, everything measured is
cool, and red is reserved for the one thing red should mean here -- a negative
eigenvalue, a direction in which the action goes down.
"""

from manim import ManimColor

# -- ground and ink -------------------------------------------------------
BG        = ManimColor("#0A0C11")
BG_PANEL  = ManimColor("#11141B")
INK       = ManimColor("#EDEAE3")      # primary text
INK_DIM   = ManimColor("#98A0AB")      # secondary text, axis labels
RULE      = ManimColor("#2A3039")      # hairlines, grids

# -- the objects ----------------------------------------------------------
CURVE     = ManimColor("#FFB454")      # the loop q(t)
CURVE_DIM = ManimColor("#7A5A2E")      # the loop, de-emphasised
GHOST     = ManimColor("#3C4450")      # previous iterate, construction lines

# -- measurement ----------------------------------------------------------
COOL      = ManimColor("#6FD3E0")      # actions, residuals, anything measured
COOL_DEEP = ManimColor("#2E7C8C")
GOOD      = ManimColor("#7BE495")      # positive eigenvalue, a passed gate
BAD       = ManimColor("#FF6B6B")      # negative eigenvalue, a failed gate
ZERO      = ManimColor("#B08CE0")      # the gauge null space
GOLD      = ManimColor("#E8C170")      # emphasis in text
ROSE      = ManimColor("#E08CA8")

# the N bodies: one hue ramp, dark to light, so body 0 is always identifiable
BODY_COLORS = [
    ManimColor(c) for c in
    ("#FFD79A", "#FFB454", "#F58E5C", "#E0687A", "#C05A9E", "#8C63C0",
     "#5C7AD0", "#3E9FC8", "#3FBFAE", "#6FD08A", "#A8D76A", "#D8C95E")
]


def body_color(k: int, n: int):
    """Evenly spaced colours around the ramp, so N = 3 and N = 12 both read."""
    if n <= 1:
        return BODY_COLORS[1]
    i = int(round(k * (len(BODY_COLORS) - 1) / (n - 1)))
    return BODY_COLORS[i]


# -- type -----------------------------------------------------------------
FONT_BODY = "Avenir Next"     # narration and labels
FONT_MATH = "Palatino"        # formulas
FONT_CODE = "Menlo"           # command lines and record fields

# sizes are in manim units at the default 8-unit-high frame
SZ_TITLE   = 0.62
SZ_HEAD    = 0.42
SZ_BODY    = 0.34
SZ_CAPTION = 0.32
SZ_MATH    = 0.44
SZ_SMALL   = 0.26
SZ_TINY    = 0.21

# -- spacing --------------------------------------------------------------
CAPTION_Y  = -3.30            # the bottom of the narration band
CAPTION_MAX_HEIGHT = 1.30     # three lines; longer text is eased down to fit
MARGIN     = 0.55

# -- pacing ---------------------------------------------------------------
# Captions are held long enough to be read aloud: this is the reading rate the
# dwell time is computed from, plus a breath at the end of every sentence.
WORDS_PER_MINUTE = 145.0
BREATH = 0.75                 # seconds added to every caption
MIN_DWELL = 1.8
