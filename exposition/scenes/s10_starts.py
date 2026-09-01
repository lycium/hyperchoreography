"""Where the starting guesses come from."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, UP, DOWN, LEFT, RIGHT,
                   ORIGIN, PI, linear, rate_functions, UpdateFromAlphaFunc)

from expo import catalog, nbody, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.plots import Plot
from expo.viz import OrbitView, Projector, spin, tumble


def ngon_frequencies(N: int):
    """Transverse frequencies of the rotating N-gon, in units of the rotation rate.

    Linearise the equations about the polygon and the motion across its plane
    decouples into patterns numbered by k; each one has its own frequency, and the
    formula is a finite sum over the other bodies.
    """
    R = sum(1.0 / (4.0 * np.sin(PI * k / N)) for k in range(1, N)) ** (1.0 / 3.0)
    dl = [2 * R * np.sin(PI * l / N) for l in range(N)]
    w = {}
    for k in range(1, N):
        w[k] = np.sqrt(sum((1 - np.cos(2 * PI * k * l / N)) / dl[l] ** 3
                           for l in range(1, N)))
    return R, {k: w[k] / w[N - 1] for k in w}


class WhereStartsComeFrom(ExpoScene):
    section_number = 10
    section_title = "starts"

    def story(self):
        self.say("Everything so far has been about what happens after the start. The "
                 "start turns out to matter more than any of it.")

        stat = VGroup(
            C("120 random low-mode starts, three bodies in a plane",
              size=T.SZ_TINY, color=T.INK_DIM),
            C("circle   82", size=T.SZ_BODY, color=T.INK),
            C("eight    36", size=T.SZ_BODY, color=T.INK),
            C("covers    2", size=T.SZ_BODY, color=T.INK),
        ).arrange(DOWN, buff=0.26, aligned_edge=LEFT).move_to(UP * 1.0)
        self.play(FadeIn(stat), run_time=1.1)
        self.say("A hundred and twenty random starts in the plane give the circle, the "
                 "figure eight, and covers of the circle. Nothing else.")
        self.say("In three dimensions it is worse: fifty-eight thousand random starts "
                 "produced the circle and the eight, and not one genuinely "
                 "three-dimensional orbit.")
        self.say("Random guessing does not reach into the interesting part of the "
                 "space. Something has to point at it.")

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        P = nbody.Action(5, 3, 24)
        x = targets.circle_start(P, 1)
        x = x * P.optimal_scale(x)
        view = OrbitView(lambda ts: P.bodies(x, ts), 5, radius=1.55,
                         center=LEFT * 3.4 + UP * 0.2, projector=Projector(3),
                         samples=600, curve_width=2.6, fixed_scale=True)
        self.say_with("The something is a linear calculation about the simplest orbit "
                      "there is.", FadeIn(view), run_time=1.2)
        self.say_with("Take N bodies at the corners of a rotating polygon, and push one "
                      "of them a little way out of the plane.",
                      spin(view, 1.0), run_time=7.0)

        f = M("ω%s%s  =  ∑%s  ( 1 − cos 2π<i>kl</i>/<i>N</i> ) / <i>d</i>%s"
              % (sb("k"), sp("2"), sb("<i>l</i> ≠ 0"), sb("<i>l</i>") + sp("3")),
              size=0.46, color=T.GOLD).move_to(RIGHT * 3.1 + UP * 1.6)
        self.play(FadeIn(f), run_time=1.0)
        self.say_with("Linearise, and the wobble splits into independent patterns, one "
                      "for each way of distributing the push around the ring. Each "
                      "pattern has its own frequency, in closed form.",
                      spin(view, 1.0), run_time=10.0)

        tab = VGroup(C("N   pattern   frequency", size=T.SZ_TINY, color=T.INK_DIM))
        for N in (4, 5, 6, 7):
            _R, w = ngon_frequencies(N)
            best = max(w, key=lambda k: w[k])
            tab.add(C("%-3d    k = %d      %.4f" % (N, best, w[best]),
                      size=T.SZ_SMALL, color=T.INK))
        tab.arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        tab.move_to(RIGHT * 3.1 + DOWN * 0.9)
        self.play(FadeIn(tab), run_time=1.1)

        self.say_with("Now the useful part. If that frequency is a simple ratio of "
                      "whole numbers, the wobble closes up after a whole number of "
                      "turns.",
                      spin(view, 1.0), run_time=7.0)
        self.say_with("And then the whole motion is periodic again -- which is the "
                      "one thing a choreography has to be.",
                      spin(view, 1.0), run_time=6.5)
        self.say_with("So the program builds starts at simple ratios near those "
                      "frequencies: a circle at one mode, plus a transverse "
                      "oscillation at another.",
                      spin(view, 1.0), run_time=9.0)
        self.say("That is where the hip-hop's modes five and six came from -- and it is "
                 "how the search gets out of the plane at all.")

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        fam = [
            ("random", "low modes, random amplitudes"),
            ("torus", "a rotation in each of several orthogonal planes at once"),
            ("vertical", "a rotating polygon plus one transverse resonance"),
            ("hyper", "several transverse resonances at once, in polarised pairs"),
            ("inplane", "the other half of the same linearisation"),
            ("kick", "a catalogued orbit, embedded and pushed along its softest "
                     "directions"),
        ]
        rows = VGroup()
        for name, desc in fam:
            rows.add(VGroup(C(name, size=T.SZ_BODY, color=T.GOLD),
                            B(desc, size=T.SZ_SMALL, color=T.INK_DIM))
                     .arrange(RIGHT, buff=0.45, aligned_edge=DOWN))
        rows.arrange(DOWN, buff=0.34, aligned_edge=LEFT).move_to(UP * 0.2)
        self.say_with("There are six families of start in the program, and this "
                      "reasoning is behind most of them.",
                      FadeIn(rows[0], shift=RIGHT * 0.1), run_time=1.2)
        for r in rows[1:]:
            self.play(FadeIn(r, shift=RIGHT * 0.1), run_time=0.40)
        self.say("The one called hyper is the one that matters most: it stacks several "
                 "resonant patterns at once, in pairs that turn rather than oscillate.")
        self.say("Every record in the catalogue with an effective dimension of five or "
                 "more came out of that one family.", extra=0.5)

        self.play(FadeOut(rows), *self.caption_anims(None), run_time=1.0)
