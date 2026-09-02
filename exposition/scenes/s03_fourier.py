"""The loop as a Fourier series, and the finite problem that comes out of it."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, Transform, UP, DOWN, LEFT,
                   RIGHT, ORIGIN, PI, linear, UpdateFromAlphaFunc, Line)

from expo import catalog, nbody, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.viz import OrbitView, Projector, spin, spectrum_bars
from expo.plots import Plot


class TheFourierLoop(ExpoScene):
    section_number = 3
    section_title = "the unknown"

    def story(self):
        rec = catalog.load("eight")
        P = nbody.Action(3, 2, 24)
        x = targets.embed(rec, P)
        frame = targets.principal_frame(P, x)

        # ------------------------------------------------------------------
        series = M("<i>q</i>(<i>t</i>)  =  ∑%s  <i>c</i>%s cos <i>mt</i>  +  "
                   "<i>s</i>%s sin <i>mt</i>" % (sb("m"), sb("m"), sb("m")),
                   size=0.50).move_to(UP * 2.35)
        self.say_with("The unknown is a closed curve, so write it as a Fourier series. "
                      "Each mode is a vector coefficient, one for the cosine and one "
                      "for the sine.",
                      FadeIn(series, shift=UP * 0.14), run_time=1.4)

        # build the eight mode by mode
        view = OrbitView(lambda ts: P.bodies(x, ts), P.N, radius=1.65,
                         center=LEFT * 3.5 + DOWN * 0.35,
                         projector=Projector(2, frame), samples=720,
                         curve_width=3.0, show_bodies=False, fixed_scale=True)
        power = np.einsum("mca,mca->m", P.shaped(x), P.shaped(x))
        bars = spectrum_bars(P.modes, power, width=5.4, height=1.7)
        bars.move_to(RIGHT * 3.3 + DOWN * 1.0)
        blab = C("power in each mode", size=T.SZ_TINY, color=T.INK_DIM)
        blab.next_to(bars, UP, buff=0.35).align_to(bars, LEFT)

        partial = np.zeros_like(P.shaped(x))
        view.set_coefficients(lambda ts: P.bodies(partial, ts), rescale=False)
        self.play(FadeIn(view), FadeIn(bars), FadeIn(blab), run_time=1.0)

        order = np.argsort(-power)
        keep = [int(i) for i in order[:6]]
        self.say("Adding them one at a time, in order of how much they matter:")
        for step, i in enumerate(keep):
            partial[i] = P.shaped(x)[i]
            hl = bars.bars[i].copy().set_fill(T.CURVE, opacity=1.0)
            self.play(view.animate.set_coefficients(
                          lambda ts, pp=partial.copy(): P.bodies(pp, ts)),
                      Transform(bars.bars[i], hl),
                      run_time=0.85)
            self.wait(0.35)
        self.play(view.animate.set_coefficients(lambda ts: P.bodies(x, ts)),
                  run_time=1.0)
        self.say("Six modes already draw the figure eight. The catalogue's copy runs "
                 "out to mode fifty, because the last few digits of the answer live in "
                 "the tail.")

        self.wipe(run_time=0.9)

        # ------------------------------------------------------------------
        # modes that are multiples of N do nothing but move the centre of mass
        self.say("One family of modes has to go.")

        drop = M("<i>m</i>  ≡  0   (mod <i>N</i>)", size=0.54, color=T.BAD)
        drop.move_to(UP * 2.4)
        self.play(FadeIn(drop), run_time=0.9)

        def with_cm(a):
            def f(ts):
                base = P.bodies(x, ts)
                extra = np.stack([a * np.cos(P.N * ts), a * np.sin(P.N * ts)], axis=1)
                return base + extra[None, :, :]
            return f

        view2 = OrbitView(with_cm(0.0), P.N, radius=1.5, center=DOWN * 0.55,
                          projector=Projector(2, frame), samples=720,
                          curve_width=2.8, fixed_scale=True)
        view2.radius = 1.5
        self.play(FadeIn(view2), run_time=0.9)

        def push(m, a):
            amp = 0.55 * np.sin(2 * PI * a)
            m.set_coefficients(with_cm(amp), rescale=False)
            m.set_time(a * 2 * PI)

        self.say_with("A mode whose number is a multiple of N shifts every body by the "
                      "same vector at the same instant.",
                      UpdateFromAlphaFunc(view2, push), run_time=7.0)

        self.say_with("The shape never changes, only where it sits. That is the centre "
                      "of mass drifting, and it is not a new orbit — so those modes "
                      "are simply left out of the basis.",
                      UpdateFromAlphaFunc(view2, push), run_time=8.0)

        self.wipe(run_time=0.9)

        # ------------------------------------------------------------------
        self.say("One more freedom is worth removing before we start.")
        kep = M("<i>q</i> → λ<i>q</i>,    <i>t</i> → λ%s <i>t</i>"
                % sp("(α+2)/2"), size=0.52, color=T.GOLD).move_to(UP * 1.7)
        self.play(FadeIn(kep), run_time=1.0)
        self.say("Kepler's scaling turns any solution into a whole family of larger, "
                 "slower ones. So the period is nailed down at two pi, and the size is "
                 "fixed by rescaling each start to where the action is stationary.")

        # ------------------------------------------------------------------
        fin = VGroup(
            M("<i>x</i>  =  ( <i>c</i>%s, <i>s</i>%s, … )   ∈   ℝ%s"
              % (sb("1"), sb("1"), sp("n")), size=0.50),
            M("<i>n</i>  =  2 · (modes) · <i>d</i>", size=0.42, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.45).move_to(DOWN * 0.5)
        self.play(FadeIn(fin), run_time=1.2)

        self.say("What is left is a list of numbers. The action is now an ordinary "
                 "function of that list, and its gradient and its second derivatives "
                 "are all available in closed form.")
        self.say("Finding a choreography has become the problem of finding a point "
                 "where that gradient is zero.", extra=0.6)

        self.play(FadeOut(VGroup(kep, fin)), *self.caption_anims(None), run_time=1.1)
