"""Static plates for the Allura Studio comp: the title card and the end card.

Each scene is one full-frame transparent PNG, rendered by `make_plates.sh` with
manim's `-s -t` (last frame, alpha kept). Everything is set in the film's own
faces and palette, and the orbit art is drawn by the same OrbitView the film
uses, so the cards read as part of the picture rather than as a wrapper around
it. Positioning happens HERE, not in the comp: Allura stacks the plates at
identity and owns only their timing, fades and processing.
"""

import numpy as np
from manim import Scene, VGroup, ORIGIN, UP, DOWN, PI

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from expo import catalog, theme as T
from expo.mathtext import B, C
from expo.viz import OrbitView, Projector


def glow_orbit(name: str, radius: float, center, alpha: float, turn=(),
               show_bodies=False, samples=2000):
    """The loop drawn through the film's own view, with a baked halo: the same
    curve restated at falling widths and opacities underneath a bright core.
    `turn` mixes principal axes into the camera — the film's own gesture for
    showing that a high-deff orbit is not a circle seen edge-on."""
    orbit = catalog.load(name)
    proj = Projector(orbit.d, orbit.principal_frame())
    for i, j, angle in turn:
        proj = proj.rotated(i, j, angle)
    g = VGroup()
    for width, a in ((16.0, 0.05), (10.0, 0.09), (6.0, 0.16), (3.2, 0.34), (1.7, 1.0)):
        view = OrbitView(orbit.bodies, orbit.N, radius=radius, center=center,
                         projector=proj, samples=samples, curve_width=width,
                         show_bodies=show_bodies and width < 2.0,
                         curve_color=T.CURVE, dot_radius=0.075)
        view.curve.set_stroke(opacity=a * alpha)
        g.add(view)
    return g


class TitleOrbit(Scene):
    """The figure eight, softly glowing, with its three bodies frozen in place."""
    def construct(self):
        self.add(glow_orbit("eight", radius=2.9, center=UP * 0.55, alpha=0.55,
                            show_bodies=True))


class TitleName(Scene):
    def construct(self):
        self.add(B("hyperchoreography", size=T.SZ_TITLE, color=T.INK)
                 .move_to(DOWN * 1.55))


class EndOrbit(Scene):
    """A spatial three-body choreography, faint behind the closing text: the
    film opens on three bodies tracing the planar eight and closes on the same
    three somewhere richer. (The d = 11 champion was auditioned and declined —
    every two-axis shadow of it is an ellipse, radial spread 0.1% on the leading
    pair, which is exactly the edge-on-circle trap s11 warns about.)"""
    def construct(self):
        self.add(glow_orbit("n3_spatial", radius=3.4, center=ORIGIN, alpha=0.32,
                            samples=3000))


class EndSubscribe(Scene):
    def construct(self):
        self.add(B("subscribe — the search continues at d = 13",
                   size=T.SZ_HEAD, color=T.GOLD).move_to(UP * 1.55))


class EndCatalogue(Scene):
    def construct(self):
        g = VGroup(B("the catalogue", size=T.SZ_SMALL, color=T.INK_DIM),
                   C("lycium.github.io/hyperchoreography", size=T.SZ_BODY, color=T.INK))
        g.arrange(DOWN, buff=0.16).move_to(DOWN * 0.15)
        self.add(g)


class EndCode(Scene):
    def construct(self):
        g = VGroup(B("code and data", size=T.SZ_SMALL, color=T.INK_DIM),
                   C("github.com/lycium/hyperchoreography", size=T.SZ_BODY, color=T.INK))
        g.arrange(DOWN, buff=0.16).move_to(DOWN * 1.45)
        self.add(g)


class EndAllura(Scene):
    def construct(self):
        self.add(B("composed with allura studio", size=T.SZ_SMALL, color=T.INK_DIM)
                 .move_to(DOWN * 3.0))
