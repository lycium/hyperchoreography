"""Phase two: Newton on the gradient, which does not care about signs."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, UP, DOWN, LEFT, RIGHT,
                   ORIGIN, PI, Dot, linear, rate_functions, UpdateFromAlphaFunc)

from expo import nbody, optim, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, table
from expo.plots import Plot, growing_line
from expo.readout import Field, sci
from expo.viz import OrbitView, Projector, eigen_strip


class PhaseTwo(ExpoScene):
    section_number = 6
    section_title = "phase two"

    def story(self):
        feat = targets.eight_newton()
        P, tr = feat.P, feat.run.phase2

        self.say("Phase two drops the idea of going downhill altogether, and solves "
                 "the equation we actually want.")

        eq = M("∇<i>A</i>(<i>x</i>)  =  0", size=0.62, color=T.GOLD).move_to(UP * 2.3)
        sub = B("n equations, n unknowns", size=T.SZ_SMALL, color=T.INK_DIM)
        sub.next_to(eq, DOWN, buff=0.32)
        self.play(FadeIn(eq), FadeIn(sub), run_time=1.0)
        self.say("Nothing in that equation says minimum. It is satisfied at a peak, at "
                 "a pass and at a basin alike, and Newton's method solves it without "
                 "ever being told which is which.")

        newton = M("<i>H</i> δ  =  − ∇<i>A</i>", size=0.54).move_to(UP * 0.75)
        self.play(FadeIn(newton), run_time=0.9)
        self.say("Linearise, solve for the step, repeat. H is the matrix of second "
                 "derivatives, which this program has in closed form.")

        self.play(FadeOut(VGroup(eq, sub)), newton.animate.move_to(UP * 2.5).scale(0.8),
                  run_time=1.0)

        step = M("δ  =  − ∑%s  <i>q</i>%s  "
                 "λ%s / ( λ%s + μ )  ( <i>q</i>%s · ∇<i>A</i> )"
                 % (sb("k"), sb("k"), sb("k"), sb("k") + sp("2"), sb("k")),
                 size=0.46, color=T.GOLD).move_to(UP * 1.75)
        self.play(FadeIn(step), run_time=1.2)
        self.say("Written in the eigenvectors of H, the step falls apart into "
                 "independent pieces: one number per direction, each scaled by that "
                 "direction's own curvature.")

        lam = np.linspace(-3.0, 3.0, 601)
        mu = 0.25
        plot = Plot((-3, 3), (-2.4, 2.4), width=5.2, height=2.0,
                    center=DOWN * 0.55, x_ticks=[-3, 0, 3], y_ticks=[],
                    x_label="curvature of the direction")
        pure = plot.line(lam[np.abs(lam) > 0.34],
                         np.clip(1.0 / lam[np.abs(lam) > 0.34], -2.4, 2.4),
                         color=T.GHOST, width=2.0, smooth=False)
        damped = plot.line(lam, lam / (lam ** 2 + mu), color=T.COOL, width=3.0)
        zeroline = plot.hline(0, dashed=False, opacity=0.35)
        self.play(FadeIn(plot), FadeIn(zeroline), Create(pure), run_time=1.4)
        self.say("Undamped, the size of the step in a direction is one over its "
                 "curvature — which explodes wherever the curvature is near zero.")
        self.play(Create(damped), run_time=1.6)
        self.say("The damping replaces that by this. Two things about the shape of "
                 "the curve matter, and they are exactly the two problems from before.")

        p1 = B("it is odd: the sign of the curvature only turns the step around",
               size=T.SZ_SMALL, color=T.GOOD)
        p2 = B("it vanishes at zero: flat directions get no step at all",
               size=T.SZ_SMALL, color=T.ZERO)
        pts = VGroup(p1, p2).arrange(DOWN, buff=0.24).move_to(UP * 0.95)
        self.play(FadeIn(p1), run_time=0.8)
        self.say("The curve is odd. Flip the sign of the curvature and the step turns "
                 "around, but nothing else changes. A pass is no harder than a basin, "
                 "which is the whole reason this method reaches most of the catalogue.")
        self.play(FadeIn(p2), run_time=0.8)
        self.say("And it goes to zero at zero curvature, which handles a nuisance the "
                 "problem has built into it.")

        self.wipe(run_time=1.0)

        x = feat.x
        frame = targets.principal_frame(P, x)
        view = OrbitView(lambda ts: P.bodies(x, ts), P.N, radius=1.85,
                         center=LEFT * 3.2 + UP * 0.3,
                         projector=Projector(2, frame), samples=600,
                         curve_width=3.0, fixed_scale=True)
        fA = Field("action", 12, color=T.CURVE).place(LEFT * 5.6 + DOWN * 1.75)
        fA.set("%.9f" % P.value(x))
        what = B("", size=T.SZ_BODY, color=T.INK).move_to(RIGHT * 3.2 + UP * 0.9)
        self.play(FadeIn(view), FadeIn(fA), run_time=1.0)

        w1 = B("start the same loop at a different moment",
               size=T.SZ_BODY, color=T.INK).move_to(RIGHT * 3.1 + UP * 0.8)
        self.play(FadeIn(w1), run_time=0.7)

        def slide(m, a):
            xs = targets.time_shift(P, x, a * 2 * PI)
            m.set_coefficients(lambda ts, xx=xs: P.bodies(xx, ts), rescale=False)
            m.set_time(a * 2 * PI)
            fA.set("%.9f" % P.value(xs))

        self.say_with("Slide the whole thing forward in time and it is the same orbit, "
                      "so the action does not move at all.",
                      UpdateFromAlphaFunc(view, slide), run_time=7.0)

        w2 = B("turn it in space", size=T.SZ_BODY, color=T.INK)
        w2.move_to(RIGHT * 3.1 + UP * 0.1)
        self.play(FadeIn(w2), run_time=0.7)

        def turn_it(m, a):
            xs = targets.rotate(P, x, a * 2 * PI)
            m.set_coefficients(lambda ts, xx=xs: P.bodies(xx, ts), rescale=False)
            m.set_time(a * 2 * PI)
            fA.set("%.9f" % P.value(xs))

        self.say_with("Turn it bodily in space and every distance between bodies is "
                      "unchanged, so again the action does not move.",
                      UpdateFromAlphaFunc(view, turn_it), run_time=7.0)

        self.say("Those are perfectly flat directions in the landscape, and they are "
                 "not telling us anything about the orbit. They are telling us that "
                 "one orbit has many descriptions.")
        self.say("The damping deals with them without having to be told they are "
                 "there: no curvature, no step.")

        self.wipe(run_time=1.0)

        g = tr.gnorm
        plot2 = Plot((0, len(g) - 1), (1e-13, 10.0), width=5.6, height=3.0,
                     center=RIGHT * 2.9 + DOWN * 0.2, log_y=True,
                     x_ticks=list(range(0, len(g))),
                     y_ticks=[1e-12, 1e-8, 1e-4, 1e0], x_label="Newton step")
        line = plot2.line(range(len(g)), g, color=T.COOL, width=3.0, smooth=False)
        dots = VGroup(*[plot2.marker(i, v, color=T.CURVE, r=0.06)
                        for i, v in enumerate(g)])
        tab = table([("step", "size of the gradient")]
                    + [("%d" % i, sci(v, 2)) for i, v in enumerate(g)],
                    align="rr", sizes=[T.SZ_TINY] + [T.SZ_SMALL] * len(g),
                    colors=[T.INK_DIM] + [T.INK_DIM if i == 0 else T.INK
                                          for i in range(len(g))],
                    col_buff=0.55, row_buff=0.20, head_buff=0.14)
        tab.move_to(LEFT * 4.3 + DOWN * 0.1)
        head, rows = tab.rows[0], tab.rows[1:]

        self.play(FadeIn(plot2), FadeIn(head), run_time=1.0)
        self.say_with("This is the real thing, finishing the figure eight off.",
                      Create(line), run_time=2.0)
        for i in range(len(g)):
            self.play(FadeIn(dots[i], scale=0.5), FadeIn(rows[i]), run_time=0.28)
        self.say("Three, two, one, then a tenth, then two thousandths — and once it is "
                 "close, the number of correct digits doubles at every step.")
        self.say("Eight steps take it from a rough shape to twelve digits. That "
                 "doubling is what makes a catalogue of a couple of thousand certified "
                 "orbits possible at all.", extra=0.6)

        self.play(FadeOut(VGroup(plot2, line, dots, rows, head)),
                  *self.caption_anims(None), run_time=1.1)
