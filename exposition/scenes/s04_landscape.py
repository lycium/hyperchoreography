"""The landscape: critical points, and why most of them are not minima."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, UP, DOWN, LEFT, RIGHT,
                   ORIGIN, PI, Dot, linear, rate_functions)

from expo import catalog, nbody, optim, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.surface import Wireframe, turn, scaled_directions
from expo.viz import OrbitView, Projector, eigen_strip, spin


def slice_of_action(P, x, u, v, span_a=0.5, span_b=None, n=41):
    """A[x + a u + b v] on a grid: real values, not a drawing of a saddle."""
    span_b = span_a if span_b is None else span_b
    a = np.linspace(-span_a, span_a, n)
    b = np.linspace(-span_b, span_b, n)
    Z = np.empty((n, n))
    xf = P.flat(x)
    for j, bj in enumerate(b):
        for i, ai in enumerate(a):
            Z[j, i] = P.value(xf + ai * u + bj * v)
    return a, b, Z


class TheLandscape(ExpoScene):
    section_number = 4
    section_title = "the landscape"

    def story(self):
        self.say("So: a function of a hundred or so numbers, and we want the places "
                 "where its gradient vanishes.")

        grad = M("∇<i>A</i>(<i>x</i>)  =  0", size=0.60, color=T.GOLD)
        grad.move_to(UP * 2.2)
        self.play(FadeIn(grad), run_time=0.9)

        self.say("The obvious move is to roll downhill. That finds the bottoms of "
                 "valleys, and it will find some real orbits. But it is not enough, "
                 "and it is worth seeing exactly why.")

        # -- the eight: a bowl ---------------------------------------------
        rec = catalog.load("eight")
        P = nbody.Action(3, 2, 16)
        x = P.flat(targets.embed(rec, P))
        H = P.hessian(x)
        w, V = np.linalg.eigh(H)

        # the spectrum shown later is the record's own, at its full 34 modes: at 16
        # the direction that makes the eight a saddle is below the truncation error
        Pf = nbody.Action(rec.N, 2, int(rec.modes.max()), modes=rec.modes)
        xf = np.array(rec.coef[:, :, :2])
        Hf = Pf.hessian(xf)
        neg, zero, _lifted = optim.inertia_gauge(Pf, xf, Hf)
        wraw = np.linalg.eigvalsh(Hf)

        # two uphill directions, each scaled by its own curvature -> a clean bowl
        (u, lu, sa), (v, lv, sbv) = scaled_directions(H, x, kind=("pos", "pos"))
        a, b, Z = slice_of_action(P, x, u, v, sa, sbv, n=41)
        wire = Wireframe(a, b, Z, width=4.6, depth=4.6, height=2.1,
                         center=LEFT * 3.2 + DOWN * 0.55, lines=19)
        A0 = P.value(x)
        mark = Dot(wire.point(0, 0, A0), radius=0.055, color=T.CURVE)

        self.play(FadeOut(grad), run_time=0.6)
        self.say_with("Here is the action itself, on a two-dimensional slice through "
                      "the figure eight. The orbit sits at the marked point.",
                      FadeIn(wire), FadeIn(mark), run_time=1.6)
        self.play(turn(wire, d_azim=1.5, extra=[(mark, (0, 0, A0))]),
                  run_time=7.0, rate_func=rate_functions.ease_in_out_sine)

        lab1 = C("the figure eight: a minimum in these directions",
                 size=T.SZ_TINY, color=T.INK_DIM)
        lab1.next_to(wire, DOWN, buff=0.25)
        self.play(FadeIn(lab1), run_time=0.6)

        # -- the saddle -----------------------------------------------------
        s = targets.saddle_n3()
        Ps, xs = s.P, s.x
        Hs = Ps.hessian(xs)
        ws, Vs = np.linalg.eigh(Hs)
        negs, zeros, _ = optim.inertia_gauge(Ps, xs, Hs)
        (us, lus, sa2), (vs, lvs, sb2) = scaled_directions(Hs, xs, kind=("neg", "pos"))
        a2, b2, Z2 = slice_of_action(Ps, xs, us, vs, sa2, sb2, n=41)
        wire2 = Wireframe(a2, b2, Z2, width=4.6, depth=4.6, height=2.1,
                          center=RIGHT * 3.2 + DOWN * 0.55, lines=19)
        As = Ps.value(xs)
        mark2 = Dot(wire2.point(0, 0, As), radius=0.055, color=T.CURVE)
        lab2 = C("another N = 3 orbit: a pass, not a basin",
                 size=T.SZ_TINY, color=T.INK_DIM)
        lab2.next_to(wire2, DOWN, buff=0.25)

        self.say_with("And here is another genuine three-body orbit, drawn the same "
                      "way. This one is not at the bottom of anything.",
                      FadeIn(wire2), FadeIn(mark2), FadeIn(lab2), run_time=1.8)
        self.play(turn(wire, d_azim=1.2, extra=[(mark, (0, 0, A0))]),
                  turn(wire2, d_azim=1.2, extra=[(mark2, (0, 0, As))]),
                  run_time=8.0, rate_func=rate_functions.ease_in_out_sine)

        self.say_with("It is a mountain pass. Walk one way and the action rises; walk "
                      "the other way and it falls. Rolling downhill will never stop "
                      "here.",
                      turn(wire, d_azim=1.2, extra=[(mark, (0, 0, A0))]),
                      turn(wire2, d_azim=1.2, extra=[(mark2, (0, 0, As))]),
                      run_time=8.0)

        self.wipe(run_time=1.0)

        # -- Morse index ----------------------------------------------------
        self.say("The number of downhill directions has a name: the Morse index. "
                 "Zero means a genuine minimum.")

        strip = eigen_strip(wraw, width=6.6, height=0.44, label=False, window=36)
        strip.move_to(UP * 1.0)
        cap = C("the figure eight:  Morse index %d,  nullity %d" % (neg, zero),
                size=T.SZ_SMALL, color=T.INK_DIM)
        cap.next_to(strip, DOWN, buff=0.28)
        self.say_with("These are the smallest second derivatives at the figure eight -- "
                      "green where the action curves up, red where it curves down, and "
                      "violet for the directions that change nothing at all.",
                      FadeIn(strip), FadeIn(cap), run_time=1.2)

        stat = VGroup(
            B("Across the catalogue:", size=T.SZ_SMALL, color=T.INK_DIM),
            B("76 of about 1800 records are minima", size=T.SZ_BODY, color=T.INK),
            B("the median Morse index is 20", size=T.SZ_BODY, color=T.INK),
        ).arrange(DOWN, buff=0.30).move_to(DOWN * 1.4)
        self.play(FadeIn(stat), run_time=1.2)
        self.say("And the catalogue is blunt about it. Almost nothing in there is a "
                 "minimum; the typical orbit has twenty directions in which the action "
                 "goes down.")
        self.say("A search that can only find minima would miss nearly all of them.",
                 extra=0.5)

        self.play(FadeOut(VGroup(strip, cap, stat)),
                  *self.caption_anims(None), run_time=1.0)
