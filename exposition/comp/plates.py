"""Plates for the Allura Studio comp: the end card, and the title art it kept."""

import numpy as np
from manim import Scene, VGroup, ORIGIN, UP, DOWN, PI, linear

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from expo import catalog, theme as T
from expo.mathtext import B, C
from expo.viz import OrbitView, Projector


def glow_orbit(name: str, radius: float, center, alpha: float, turn=(),
               show_bodies=False, samples=2000):
    """The loop drawn through the film's own view, with a baked halo: the same"""
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
    """A spatial three-body choreography, faint behind the closing text: the"""
    def construct(self):
        self.add(glow_orbit("n3_spatial", radius=3.4, center=ORIGIN, alpha=0.32,
                            samples=3000))


class EndOrbitTurning(Scene):
    """EndOrbit, turning: the same curve, the same baked halo, the same faintness,"""

    SECONDS = 14.6
    SWEEP = PI

    def construct(self):
        from manim import UpdateFromAlphaFunc
        from expo.viz import default_tumble_planes

        self.camera.background_color = T.BG
        g = glow_orbit("n3_spatial", radius=3.4, center=ORIGIN, alpha=0.32,
                       samples=3000)
        self.add(g)
        orbit = catalog.load("n3_spatial")
        base = Projector(orbit.d, orbit.principal_frame())
        planes = default_tumble_planes(orbit.d)

        def turn(m, a):
            p = base
            for i, (u, v) in enumerate(planes):
                p = p.rotated(u, v, self.SWEEP * a / (i + 1))
            for view in m:
                view.set_projector(p)

        self.play(UpdateFromAlphaFunc(g, turn), run_time=self.SECONDS,
                  rate_func=linear)


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


class EndCredits(Scene):
    """Who made it and what with, quiet under the links. Replaces the older single"""
    def construct(self):
        g = VGroup(B("direction: thomas ludwig / lycium", size=T.SZ_SMALL, color=T.INK),
                   B("heavy lifting: claude code", size=T.SZ_SMALL, color=T.INK),
                   B("rendered with manim and allura studio", size=T.SZ_SMALL,
                     color=T.INK_DIM))
        g.arrange(DOWN, buff=0.20).move_to(DOWN * 2.95)
        self.add(g)
