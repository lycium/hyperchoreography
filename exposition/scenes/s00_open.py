"""Opening: three bodies on one curve."""

from manim import (VGroup, FadeIn, FadeOut, UP, DOWN, ORIGIN, PI,
                   linear, rate_functions)

from expo import catalog, theme as T
from expo.base import ExpoScene
from expo.mathtext import B
from expo.viz import OrbitView, Projector, spin, spin_and_place


class Opening(ExpoScene):
    section_title = ""

    def story(self):
        eight = catalog.load("eight")
        view = OrbitView(eight.bodies, eight.N, radius=2.15, center=ORIGIN,
                         projector=Projector(eight.d, eight.principal_frame()),
                         samples=900, curve_width=3.4)
        view.curve.set_stroke(opacity=0.0)
        self.add(view)

        self.play(*[FadeIn(d, scale=0.4) for d in view.dots], spin(view, 1.6 / 7.0),
                  run_time=1.6)
        self.play(spin(view), run_time=7.0, rate_func=linear)

        self.play(view.curve.animate.set_stroke(opacity=0.9), spin(view, 1.8 / 7.0),
                  run_time=1.8)
        self.play(spin(view), run_time=6.0, rate_func=linear)

        self.say_with("Three equal masses, pulling on each other by Newton's law of "
                      "gravity, and all three running along the same closed curve.",
                      spin(view), run_time=7.0)

        self.say_with("They are evenly spaced in time: each body is a third of a period "
                      "behind the one in front. An orbit like this is called a "
                      "choreography.",
                      spin(view), run_time=8.0)

        self.play(spin_and_place(view, radius=1.30, center=UP * 1.35, turns=0.5),
                  run_time=2.4, rate_func=rate_functions.ease_in_out_sine)

        title = B("hyperchoreography", size=T.SZ_TITLE, color=T.INK)
        sub = B("a search engine for N-body choreographies in any dimension",
                size=T.SZ_BODY, color=T.INK_DIM)
        g = VGroup(title, sub).arrange(DOWN, buff=0.34).move_to(DOWN * 1.15)
        self.play(FadeIn(title, shift=UP * 0.18), spin(view, 1.2 / 7.0), run_time=1.2)
        self.play(FadeIn(sub, shift=UP * 0.10), spin(view, 1.0 / 7.0), run_time=1.0)

        self.say_with("This one is the figure eight, found by Cristopher Moore in 1993 "
                      "and proved to exist by Chenciner and Montgomery in 2000.",
                      spin(view), run_time=8.0)

        self.say_with("This film is about a program that goes looking for the others.",
                      spin(view), run_time=7.0, extra=1.0)

        self.play(FadeOut(VGroup(view, g)), *self.caption_anims(None), run_time=1.6)
