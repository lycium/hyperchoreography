"""The rotating frame, and what it buys."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, UP, DOWN, LEFT, RIGHT,
                   ORIGIN, PI, linear, rate_functions, UpdateFromAlphaFunc)

from expo import catalog, nbody, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, Crow
from expo.viz import OrbitView, Projector, PlaneGrid, spin, tumble


class TheRotatingFrame(ExpoScene):
    section_number = 12
    section_title = "the rotating frame"

    def story(self):
        self.say("Up to six dimensions the search works as described. Above six it "
                 "stops working, and the fix is a change of viewpoint.")

        wall = VGroup(
            Crow("eight dimensions, no rotating frame", size=T.SZ_SMALL, color=T.INK_DIM),
            Crow("60 331 trials     0 orbits", size=T.SZ_BODY, color=T.BAD),
            Crow("eight dimensions, a calibrated frame", size=T.SZ_SMALL, color=T.INK_DIM),
            Crow("8 799 trials      6 orbits", size=T.SZ_BODY, color=T.GOOD),
        ).arrange(DOWN, buff=0.28, aligned_edge=LEFT).move_to(UP * 0.6)
        self.play(FadeIn(wall[0]), FadeIn(wall[1]), run_time=1.0)
        self.say("Sixty thousand trials in eight dimensions, with the relative "
                 "equilibria filtered out, found nothing at all.")
        self.play(FadeIn(wall[2]), FadeIn(wall[3]), run_time=1.0)
        self.say("Nine thousand trials of the same search, in a frame that is itself "
                 "turning, found six.")

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        ans = M("<i>q</i>%s(<i>t</i>)  =  exp(Ω<i>t</i>) · <i>q</i>(<i>t</i> + "
                "2π<i>j</i>/<i>N</i>)" % sb("<i>j</i>"), size=0.52, color=T.GOLD)
        ans.move_to(UP * 2.35)
        self.play(FadeIn(ans), run_time=1.0)
        self.say("The change is small to write down. Ask for a loop that closes up only "
                 "after a rotation, rather than exactly.")

        rec = catalog.load("d7_twist")
        vf = OrbitView(rec.frame_bodies, rec.N, radius=1.7,
                       center=LEFT * 3.4 + DOWN * 0.3,
                       projector=Projector(rec.d, rec.principal_frame()),
                       samples=720, curve_width=2.6, dot_radius=0.055,
                       fixed_scale=True)
        vi = OrbitView(rec.bodies, rec.N, radius=1.7,
                       center=RIGHT * 3.4 + DOWN * 0.3,
                       projector=Projector(rec.d, rec.principal_frame()),
                       samples=720, curve_width=2.6, dot_radius=0.055,
                       fixed_scale=True)
        lf = C("in the turning frame", size=T.SZ_TINY, color=T.INK_DIM)
        lf.move_to(vf.center + DOWN * 1.95)
        li = C("what actually happens", size=T.SZ_TINY, color=T.INK_DIM)
        li.move_to(vi.center + DOWN * 1.95)
        self.play(FadeIn(vf), FadeIn(vi), FadeIn(lf), FadeIn(li), run_time=1.2)
        self.say_with("A ten-body orbit in seven dimensions, drawn both ways. On the "
                      "left, seen from the turning frame; on the right, what an "
                      "observer sitting still would see.",
                      spin(vf, 1.0), spin(vi, 1.0), run_time=9.0)

        self.say_with("It is the same solution. The frame is a way of asking the "
                      "question, not a property of the answer.",
                      tumble(vf, turns=1.0, sweep=PI), tumble(vi, turns=1.0, sweep=PI),
                      run_time=8.0, rate_func=linear)

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        f1 = VGroup(
            B("the potential does not change at all", size=T.SZ_BODY, color=T.GOOD),
            B("a rotation preserves every distance, so the term that costs the most to "
              "evaluate is untouched", size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        f2 = VGroup(
            M("½ ∫ |<i>q̇</i>|%s   →   ½ ∫ |<i>q̇</i> + Ω<i>q</i>|%s"
              % (sp("2"), sp("2")), size=0.48, color=T.GOLD),
            B("the kinetic term picks up one extra matrix, and stays quadratic",
              size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.24)
        g = VGroup(f1, f2).arrange(DOWN, buff=0.85).move_to(UP * 0.3)
        self.say_with("Two facts make this cheap.", FadeIn(f1), run_time=1.2)
        self.say("The potential term is untouched, because rotating everything leaves "
                 "every distance between bodies exactly where it was.")
        self.play(FadeIn(f2), run_time=0.9)
        self.say("And the kinetic term picks up a single extra matrix while staying "
                 "quadratic, so the gradient and the Hessian stay closed form.")

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        om = VGroup(
            B("Omega is a rotation rate in each of several orthogonal planes at once",
              size=T.SZ_BODY, color=T.INK),
            Crow("--omega \"1,2\"        rates in the coordinate planes",
              size=T.SZ_SMALL, color=T.INK_DIM),
            Crow("--omega su:1,2       rates that sum to zero",
              size=T.SZ_SMALL, color=T.INK_DIM),
            Crow("--omega g2:1,6       the maximal torus of the group G2",
              size=T.SZ_SMALL, color=T.GOLD),
        ).arrange(DOWN, buff=0.30, aligned_edge=LEFT).move_to(UP * 0.9)
        self.say_with("What is left is choosing the rotation, and that turns out to "
                      "matter enormously.", FadeIn(om), run_time=1.4)
        self.say("The rates that work are not arbitrary. The ones that opened seven "
                 "dimensions, and then nine and eleven, come from the torus of an "
                 "exceptional group.")

        res = VGroup(
            C("d = 7,  N = 10     142 records, 76 of them filling all seven dimensions",
              size=T.SZ_SMALL, color=T.INK),
            C("d = 9,  N = 10      31 records, 13 filling all nine",
              size=T.SZ_SMALL, color=T.INK),
            C("d = 11, N = 12     122 records, 38 filling all eleven",
              size=T.SZ_SMALL, color=T.INK),
        ).arrange(DOWN, buff=0.26, aligned_edge=LEFT).move_to(DOWN * 1.6)
        self.play(FadeIn(res), run_time=1.1)
        self.say("That is what is in the catalogue above six dimensions, and none of it "
                 "exists without the frame.", extra=0.5)

        self.play(FadeOut(VGroup(om, res)), *self.caption_anims(None), run_time=1.1)
