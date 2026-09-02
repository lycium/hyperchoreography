"""Closing: what is in the catalogue, and what is not."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, UP, DOWN, LEFT, RIGHT, ORIGIN, PI,
                   linear, rate_functions)

from expo import catalog, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C
from expo.viz import OrbitView, Projector, spin, tumble


class Closing(ExpoScene):
    section_number = 15
    section_title = "what is there"

    def story(self):
        # — a last look ----------------------------------------------------
        gallery = [("eight", "N = 3, in a plane"),
                   ("hiphop", "N = 4, out of the plane"),
                   ("n4_deff4", "N = 4, four dimensions"),
                   ("n6_deff6", "N = 6, six dimensions"),
                   ("d7_inertial", "N = 10, seven dimensions, no frame"),
                   ("d7_twist", "N = 10, seven dimensions, in a frame"),
                   ("d11", "N = 12, eleven dimensions")]

        lines = [
            "Everything in this film was one of these: a start, two phases, a "
            "certification, and a check that it was not already in the file.",
            None,
            "The two-dimensional ones were known. The four-dimensional ones were "
            "mostly known.",
            None,
            "Above six, nothing was, because the search that finds them did not "
            "exist.",
            None,
            "This one fills eleven dimensions, needed twelve bodies and a frame built "
            "from an exceptional group to be found at all — and it is a theorem.",
        ]

        prev = None
        for i, (name, text) in enumerate(gallery):
            o = catalog.load(name)
            v = OrbitView(o.bodies, o.N, radius=2.0, center=UP * 0.5,
                          projector=Projector(o.d, o.principal_frame()),
                          samples=720, curve_width=2.6, dot_radius=0.065,
                          fixed_scale=True)
            lab = C(text, size=T.SZ_SMALL, color=T.INK_DIM).move_to(DOWN * 1.95)
            if prev is None:
                self.play(FadeIn(v), FadeIn(lab), run_time=1.0)
            else:
                self.play(FadeOut(prev[0]), FadeOut(prev[1]), FadeIn(v), FadeIn(lab),
                          run_time=0.9)
            anim = tumble(v, turns=1.0, sweep=2 * PI) if o.d > 2 else spin(v, 1.0)
            line = lines[i] if i < len(lines) else None
            if line:
                self.say_with(line, anim, run_time=7.0, rate_func=linear)
            else:
                self.play(anim, run_time=5.5, rate_func=linear)
            prev = (v, lab)

        self.say_with("Between them they are about eighteen hundred certified orbits, "
                      "in dimensions two through eleven.",
                      tumble(prev[0], turns=1.0, sweep=2 * PI), run_time=6.5,
                      rate_func=linear)

        self.play(FadeOut(prev[0]), FadeOut(prev[1]), run_time=0.9)

        # — the honest part -------------------------------------------------
        open_ = VGroup(
            VGroup(B("linear stability", size=T.SZ_BODY, color=T.INK),
                   B("the catalogue reports a Morse index, which is a statement about "
                     "the action, not about whether the orbit survives a nudge",
                     size=T.SZ_SMALL, color=T.INK_DIM))
            .arrange(DOWN, buff=0.20, aligned_edge=LEFT),
            VGroup(B("existence, record by record", size=T.SZ_BODY, color=T.INK),
                   B("every orbit the interval method has been run on is a theorem; "
                     "the rest of the file is still only certified, and the sweep is "
                     "machine time", size=T.SZ_SMALL, color=T.INK_DIM))
            .arrange(DOWN, buff=0.20, aligned_edge=LEFT),
            VGroup(B("the search above six dimensions needs a turning frame",
                     size=T.SZ_BODY, color=T.INK),
                   B("and not one high-dimensional orbit found that way has been "
                     "deformed back to one that does not", size=T.SZ_SMALL,
                     color=T.INK_DIM))
            .arrange(DOWN, buff=0.20, aligned_edge=LEFT),
        ).arrange(DOWN, buff=0.62, aligned_edge=LEFT).move_to(UP * 0.2)

        self.say_with("It is worth being clear about what has not been done.",
                      FadeIn(open_[0]), run_time=1.2)
        self.say("The Morse index says how the action behaves nearby. It does not say "
                 "whether the orbit is stable, and most of these are violently unstable.")
        self.play(FadeIn(open_[1]), run_time=0.9)
        self.say("The existence proofs are the newest part. Every record the interval "
                 "method has been given has passed, and the rest of the file is still "
                 "only certified until the sweep reaches it.")
        self.play(FadeIn(open_[2]), run_time=0.9)
        self.say("The largest open question is the last one: whether the high "
                 "dimensional orbits found in a turning frame have relatives that need "
                 "no frame at all.")
        self.say("So far, following them has always led somewhere else.")

        self.play(FadeOut(open_), run_time=1.0)

        # — sign off --------------------------------------------------------
        o = catalog.load("eight")
        v = OrbitView(o.bodies, o.N, radius=1.45, center=UP * 1.15,
                      projector=Projector(o.d, o.principal_frame()),
                      samples=900, curve_width=3.0)
        end = VGroup(
            B("hyperchoreography", size=T.SZ_HEAD, color=T.INK),
            C("github.com/lycium/hyperchoreography", size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.34).move_to(DOWN * 1.4)
        self.play(FadeIn(v), run_time=1.0)
        self.play(FadeIn(end), run_time=1.0)
        self.say_with("The whole thing is one binary and about five thousand lines.",
                      spin(v, 1.0), run_time=7.0, extra=1.2)
        self.play(FadeOut(VGroup(v, end)), *self.caption_anims(None), run_time=1.6)
