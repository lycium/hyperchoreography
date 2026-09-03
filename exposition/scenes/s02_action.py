"""The action principle, and the reduction to a single curve."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, UP, DOWN, LEFT, RIGHT,
                   ORIGIN, PI, linear, rate_functions, UpdateFromAlphaFunc)

from expo import catalog, nbody, optim, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, var
from expo.plots import Plot
from expo.surface import scaled_directions
from expo.viz import OrbitView, Projector, spin


class TheActionPrinciple(ExpoScene):
    section_number = 2
    section_title = "the action"

    def story(self):
        self.say("Newton's law has a second life as a statement about whole paths "
                 "rather than about instants.")

        act = VGroup(
            M("<i>A</i>[<i>q</i>]  =  ∫ <i>L</i> d<i>t</i>,     "
              "<i>L</i>  =  <i>T</i>  −  <i>V</i>", size=0.52),
            B("the Lagrangian: the kinetic energy less the potential energy",
              size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.36).move_to(UP * 1.9)
        self.play(FadeIn(act, shift=UP * 0.15), run_time=1.2)

        self.say("Give a path a single number, its action: the integral over the whole "
                 "motion of one quantity, the Lagrangian — the kinetic energy less the "
                 "potential energy.")

        self.say("A path is a solution of Newton's equations exactly when that number "
                 "is stationary: nudge the path a little, and the action does not "
                 "change to first order.")

        rec = catalog.load("eight")
        P = nbody.Action(rec.N, 2, 24)
        x = P.flat(targets.embed(rec, P))
        (eta, _lam, _span) = scaled_directions(P.hessian(x), x, kind=("pos",))[0]
        eta *= 0.26 * np.linalg.norm(x) / np.linalg.norm(eta)

        eps = np.linspace(-1.0, 1.0, 81)
        A = np.array([P.value(x + e * eta) for e in eps])
        A0 = P.value(x)
        span = float(A.max() - A0)

        plot = Plot((-1.0, 1.0), (A0 - 0.08 * span, A0 + 1.06 * span),
                    width=5.0, height=2.9,
                    center=DOWN * 0.85 + RIGHT * 3.05,
                    x_ticks=[-1, 0, 1], y_ticks=[],
                    x_label="size of the nudge")
        lab = C("action", size=T.SZ_TINY, color=T.INK_DIM)
        lab.next_to(plot, UP, buff=0.18).align_to(plot, LEFT)
        curve = plot.line(eps, A, color=T.COOL)
        mark = plot.marker(0, A0, color=T.CURVE)

        view = OrbitView(lambda ts: P.bodies(x, ts), P.N, radius=1.5,
                         center=DOWN * 0.7 + LEFT * 3.4,
                         projector=Projector(2, targets.principal_frame(P, x)),
                         samples=600, curve_width=3.0, show_bodies=False,
                         fixed_scale=True)
        ghost = view.curve.copy().set_stroke(T.GHOST, opacity=0.55, width=2.0)

        self.play(FadeOut(act), run_time=0.7)
        self.play(FadeIn(view), FadeIn(ghost), FadeIn(plot), FadeIn(lab), run_time=1.2)

        def wobble(m, a):
            e = 0.92 * np.sin(2 * PI * a * 2.0)
            xa = x + e * eta
            m.set_coefficients(lambda ts: P.bodies(xa, ts), rescale=False)

        def walk(m, a):
            e = 0.92 * np.sin(2 * PI * a * 2.0)
            m.move_to(plot.clamp(plot.at(e, P.value(x + e * eta))))

        self.play(Create(curve), run_time=1.6)
        self.add(mark)
        self.say_with("Here is the figure eight, and here is its action as the curve is "
                      "pushed away from it in one particular direction.",
                      UpdateFromAlphaFunc(view, wobble),
                      UpdateFromAlphaFunc(mark, walk), run_time=8.0, extra=0.4)

        self.say_with("The bottom of that curve is flat. That flatness, in every "
                      "direction at once, is the whole content of Newton's law.",
                      UpdateFromAlphaFunc(view, wobble),
                      UpdateFromAlphaFunc(mark, walk), run_time=8.0)

        self.wipe(run_time=0.9)

        tv = VGroup(
            M("<i>T</i>  =  ½ ∑%s |<i>q̇</i>%s|%s"
              % (sb(var("j")), sb(var("j")), sp("2")), size=0.50),
            M("<i>V</i>  =  − ∑%s |<i>q</i>%s − <i>q</i>%s|%s"
              % (sb(var("i") + " &lt; " + var("j")), sb(var("i")), sb(var("j")),
                 sp("−α")), size=0.50),
        ).arrange(DOWN, buff=0.46, aligned_edge=LEFT).move_to(UP * 1.7)
        sign = B("the potential is negative — the bodies pull on each other — so "
                 "subtracting it puts a plus in the action",
                 size=T.SZ_SMALL, color=T.INK_DIM).move_to(UP * 0.35)
        self.play(FadeIn(tv), run_time=1.2)
        self.say("Written out for N bodies: the kinetic energy is a sum over the "
                 "bodies, and the potential energy a sum over the pairs of them.")
        self.play(FadeIn(sign), run_time=0.8)
        self.say("The potential is negative, because the bodies pull on each other. So "
                 "subtracting it puts a plus sign in the action.")

        self.play(FadeOut(tv), FadeOut(sign), run_time=0.7)
        self.say("Now put the choreography constraint in.")

        full = M("<i>A</i>  =  ∫ [ ½ ∑%s |<i>q̇</i>%s|%s  +  ∑%s "
                 "|<i>q</i>%s − <i>q</i>%s|%s ] d<i>t</i>"
                 % (sb(var("j")), sb(var("j")), sp("2"),
                    sb(var("i") + " &lt; " + var("j")), sb(var("i")), sb(var("j")),
                    sp("−α")),
                 size=0.46).move_to(UP * 1.9)
        self.play(FadeIn(full), run_time=1.2)
        self.say("Every body is the same curve at a different time, so every term in "
                 "these sums is the same integral written N times over.")

        red = M("<i>A</i>[<i>q</i>]  =  ½ ∫ |<i>q̇</i>|%s d<i>t</i>  +  ½ ∑%s "
                "∫ |<i>q</i>(<i>t</i>) − <i>q</i>(<i>t</i> + 2π<i>k</i>/<i>N</i>)|%s "
                "d<i>t</i>"
                % (sp("2"), sb(var("k") + " = 1 … <i>N</i>−1"), sp("−α")),
                size=0.46, color=T.GOLD).move_to(UP * 0.6)
        self.play(FadeIn(red, shift=UP * 0.12), run_time=1.4)

        self.say("What is left is a functional of one closed curve. The first term is "
                 "its kinetic energy; the second is how close it comes to its own "
                 "time-shifted copies.")

        gloss = VGroup(
            VGroup(M("here <i>k</i> counts the time shifts, not the bodies;",
                     size=T.SZ_SMALL, color=T.INK_DIM),
                   M("the <i>N</i> bodies contribute the same integral, so this is "
                     "the action per body", size=T.SZ_SMALL, color=T.INK_DIM),
                   ).arrange(DOWN, buff=0.22),
            M("α = 1 is Newtonian gravity; the program takes any α",
              size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.40).move_to(DOWN * 0.60)
        self.play(FadeIn(gloss[0]), run_time=0.8)
        self.say("The index k there counts the shifts in time, not the bodies — and "
                 "since the N bodies all contribute the same integral, this is the "
                 "action of a single one of them.")
        self.play(FadeIn(gloss[1]), run_time=0.8)
        self.say("The exponent is left free. Newtonian gravity is alpha equal to one, "
                 "and other values are useful as a way of walking from an easier "
                 "problem to this one.")

        self.play(FadeOut(gloss), *self.caption_anims(None), run_time=0.6)

        palais = B("Palais' principle of symmetric criticality", size=T.SZ_BODY,
                   color=T.INK).move_to(DOWN * 1.35)
        sub = B("a critical point of the reduced functional is a genuine solution, "
                "not merely the best symmetric one", size=T.SZ_SMALL,
                color=T.INK_DIM).next_to(palais, DOWN, buff=0.28)
        self.play(FadeIn(palais), FadeIn(sub), run_time=1.2)
        self.say("There is a subtlety worth naming. We have thrown away every "
                 "non-symmetric direction, so why should a critical point of what is "
                 "left be a critical point of the whole thing?")
        self.say("Because the symmetry acts by isometries, and a theorem of Palais "
                 "says that is enough.")
        self.say("Anything critical here is a real solution of the N-body problem.",
                 extra=0.6)

        self.play(FadeOut(VGroup(full, red, palais, sub)),
                  *self.caption_anims(None), run_time=1.2)
