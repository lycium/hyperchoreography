"""What a choreography is: N bodies, N curves, and then only one."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, Transform, UP, DOWN, LEFT,
                   RIGHT, ORIGIN, PI, linear, rate_functions)

from expo import catalog, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, frac
from expo.viz import OrbitView, Projector, spin, spin_and_place, tumble


class WhatIsAChoreography(ExpoScene):
    section_number = 1
    section_title = "the problem"

    def story(self):
        # — Newton ------------------------------------------------------
        law = M("<i>q̈</i>%s  =  " % sb("<i>k</i>"), size=0.52)
        rhs = M("∑%s  <i>m</i>%s (<i>q</i>%s − <i>q</i>%s) / "
                "|<i>q</i>%s − <i>q</i>%s|%s"
                % (sb("<i>j</i> ≠ <i>k</i>"), sb("<i>j</i>"), sb("<i>j</i>"),
                   sb("<i>k</i>"), sb("<i>j</i>"), sb("<i>k</i>"), sp("3")),
                size=0.52)
        eq = VGroup(law, rhs).arrange(RIGHT, buff=0.12).move_to(UP * 0.9)

        self.say_with("We start where everyone starts: N point masses, each one "
                      "accelerated by the pull of all the others.",
                      FadeIn(eq, shift=UP * 0.15), run_time=1.6)

        note = B("N unknown curves, coupled to each other", size=T.SZ_SMALL,
                 color=T.INK_DIM).next_to(eq, DOWN, buff=0.7)
        self.say_with("That is N unknown curves in space, all tangled together. "
                      "Solving it in general is hopeless.",
                      FadeIn(note), run_time=1.4)

        # — the choreography constraint ---------------------------------
        self.play(FadeOut(note), eq.animate.move_to(UP * 2.55).scale(0.72),
                  run_time=1.2)

        cons = M("<i>q</i>%s(<i>t</i>)  =  <i>q</i>(<i>t</i> + <i>kT</i>/<i>N</i>)"
                 % sb("<i>k</i>"), size=0.62, color=T.GOLD).move_to(UP * 1.35)
        self.say_with("A choreography is the boldest possible guess about the answer: "
                      "that all N bodies trace out one and the same curve.",
                      FadeIn(cons, scale=1.06), run_time=1.6)

        eight = catalog.load("eight")
        view = OrbitView(eight.bodies, eight.N, radius=1.55, center=DOWN * 1.35,
                         projector=Projector(eight.d, eight.principal_frame()),
                         samples=900, curve_width=3.2)
        self.play(FadeIn(view), run_time=1.2)
        self.say_with("One curve q, and body number k is simply body zero seen a "
                      "fraction k over N of a period later.",
                      spin(view, 1.0), run_time=7.0)

        self.say_with("So instead of N unknown functions there is one, and instead of "
                      "a hopeless problem there is a hard but finite one.",
                      spin(view, 1.0), run_time=7.5)

        # — it is a real restriction, and it still leaves plenty ---------
        self.play(FadeOut(VGroup(eq, cons)), spin_and_place(view, radius=2.0,
                  center=UP * 0.35, turns=0.4), run_time=1.6)

        self.say_with("It is a genuine restriction. Most solutions of the N-body "
                      "problem are not choreographies at all.",
                      spin(view, 0.9), run_time=6.5)

        gallery = [("hiphop", "N = 4, out of the plane"),
                   ("n4_deff4", "N = 4, filling four dimensions"),
                   ("d7_twist", "N = 10, filling seven"),
                   ("d11", "N = 12, filling eleven")]
        self.play(FadeOut(view), run_time=0.8)

        prev = None
        for name, text in gallery:
            o = catalog.load(name)
            v = OrbitView(o.bodies, o.N, radius=2.05, center=UP * 0.45,
                          projector=Projector(o.d, o.principal_frame()),
                          samples=720, curve_width=2.6, dot_radius=0.070,
                          fixed_scale=True)
            lab = C(text, size=T.SZ_SMALL, color=T.INK_DIM)
            lab.move_to(DOWN * 1.95)
            if prev is None:
                self.play(FadeIn(v), FadeIn(lab), run_time=0.9)
            else:
                self.play(FadeOut(prev[0]), FadeOut(prev[1]),
                          FadeIn(v), FadeIn(lab), run_time=0.9)
            self.play(tumble(v, turns=1.0, sweep=2 * PI), run_time=6.0,
                      rate_func=linear)
            prev = (v, lab)

        self.say("Turning the picture is the only honest way to look at these: in any "
                 "two directions you happen to pick, a loop that fills eleven "
                 "dimensions looks exactly like a circle.")
        self.play(tumble(prev[0], turns=1.0, sweep=2 * PI), run_time=6.0,
                  rate_func=linear)

        self.say_with("There are enough of them to fill a catalogue, in dimensions "
                      "nobody had looked in before.",
                      tumble(prev[0], turns=1.0, sweep=2 * PI),
                      run_time=6.5, extra=0.4)

        self.play(FadeOut(prev[0]), FadeOut(prev[1]),
                  *self.caption_anims(None), run_time=1.2)
