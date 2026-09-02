"""Three orbits that need more than a random start."""

import numpy as np
from manim import (VGroup, VMobject, FadeIn, FadeOut, Create, Line, Dot, UP, DOWN,
                   LEFT, RIGHT, ORIGIN, PI, linear, rate_functions,
                   UpdateFromAlphaFunc)

from expo import catalog, nbody, optim, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, Crow
from expo.plots import Plot
from expo.readout import Field, sci
from expo.viz import OrbitView, Projector, spin, tumble, spectrum_bars


class Epicycles(VGroup):
    """Two arms, one turning forward and one back, with the tip's trace.

    This is the whole content of a two-mode Fourier series drawn as a machine: the
    first arm turns at rate m1, the second is carried by it and turns at rate -m2,
    and the tip draws whatever the sum of the two is.
    """

    def __init__(self, m1=2, m2=3, a=1.0, b=0.35, scale=1.15, center=ORIGIN,
                 samples=720):
        super().__init__()
        self.m1, self.m2, self.a, self.b = m1, m2, a, b
        self.scale_, self.center_ = scale, np.asarray(center, dtype=float)
        self.samples = samples

        self.arm1 = Line(ORIGIN, RIGHT, stroke_width=2.6, color=T.COOL)
        self.arm2 = Line(ORIGIN, RIGHT, stroke_width=2.6, color=T.ROSE)
        self.circle1 = VMobject(stroke_width=1.2)
        self.circle1.set_stroke(T.COOL, opacity=0.35)
        self.circle1.set_fill(opacity=0)
        self.trace = VMobject(stroke_width=3.0)
        self.trace.set_stroke(T.CURVE)
        self.trace.set_fill(opacity=0)
        self.tip = Dot(ORIGIN, radius=0.062, color=T.CURVE)
        self.joint = Dot(ORIGIN, radius=0.045, color=T.COOL)
        self.add(self.circle1, self.trace, self.arm1, self.arm2, self.joint, self.tip)
        self.set_state(0.0, b)

    def p1(self, t):
        return self.center_ + self.scale_ * self.a * np.array(
            [np.cos(self.m1 * t), np.sin(self.m1 * t), 0.0])

    def p2(self, t, b):
        return self.p1(t) + self.scale_ * b * np.array(
            [np.cos(self.m2 * t), -np.sin(self.m2 * t), 0.0])

    def set_state(self, t: float, b: float = None, drawn: float = 1.0):
        if b is not None:
            self.b = b
        ts = np.linspace(0, 2 * PI, self.samples, endpoint=False)
        pts = np.array([self.p2(s, self.b) for s in ts])
        n = max(3, int(drawn * len(pts)))
        self.trace.set_points_smoothly(pts[:n] if drawn < 1.0
                                       else np.vstack([pts, pts[:1]]))
        c1 = np.array([self.p1(s) for s in ts])
        self.circle1.set_points_smoothly(np.vstack([c1, c1[:1]]))
        self.arm1.put_start_and_end_on(self.center_, self.p1(t))
        if self.b > 1e-6:                       # a zero-length arm is not a line
            self.arm2.put_start_and_end_on(self.p1(t), self.p2(t, self.b))
            self.arm2.set_stroke(opacity=1.0)
        else:
            self.arm2.set_stroke(opacity=0.0)
        self.joint.move_to(self.p1(t))
        self.tip.move_to(self.p2(t, self.b))
        return self


class HardCases(ExpoScene):
    section_number = 7
    section_title = "harder cases"

    def story(self):
        self.say("A random start finds the easy orbits and nothing else. Here are "
                 "three that need a start with an idea behind it.")

        # ==================================================================
        # 1. the five-fold N = 4 orbit
        # ==================================================================
        head = B("one: two circles turning opposite ways", size=T.SZ_HEAD,
                 color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(head), run_time=0.8)

        ep = Epicycles(2, 3, 1.0, 0.0, scale=1.25, center=LEFT * 3.3 + DOWN * 0.3)
        self.play(FadeIn(ep), run_time=1.0)

        self.say_with("Start with one arm turning twice per period. Its tip goes round "
                      "a circle — which is a perfectly good, if dull, choreography.",
                      UpdateFromAlphaFunc(ep, lambda m, a: m.set_state(a * 2 * PI, 0.0)),
                      run_time=7.0, rate_func=linear)

        f2 = M("<i>q</i>(<i>t</i>)  =  <i>a</i> e%s  +  <i>b</i> e%s"
               % (sp("2<i>it</i>"), sp("−3<i>it</i>")), size=0.50)
        f2.move_to(RIGHT * 3.2 + UP * 1.4)
        self.play(FadeIn(f2), run_time=0.9)

        self.say_with("Now hang a second arm off its tip, turning three times per "
                      "period the other way.",
                      UpdateFromAlphaFunc(ep, lambda m, a: m.set_state(a * 2 * PI,
                                                                       0.35 * a)),
                      run_time=7.0, rate_func=linear)

        self.say_with("The tip stops drawing a circle and starts drawing petals.",
                      UpdateFromAlphaFunc(ep, lambda m, a: m.set_state(a * 2 * PI, 0.35)),
                      run_time=7.0, rate_func=linear)

        why = VGroup(
            M("2 − (−3)  =  5", size=0.52, color=T.GOLD),
            B("the two arms line up five times per period", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.30).move_to(RIGHT * 3.2 + DOWN * 0.9)
        self.play(FadeIn(why), run_time=1.0)
        self.say_with("Five petals, and it is the difference of the two rates that says "
                      "so: the arms come back into line five times in a period, and the "
                      "figure closes up.",
                      UpdateFromAlphaFunc(ep, lambda m, a: m.set_state(a * 2 * PI, 0.35)),
                      run_time=8.0, rate_func=linear)

        self.say_with("That is a guess, not an orbit. Nothing about it satisfies "
                      "Newton's equations yet.",
                      UpdateFromAlphaFunc(ep, lambda m, a: m.set_state(a * 2 * PI, 0.35)),
                      run_time=7.0, rate_func=linear)

        # — hand it to the optimiser ---------------------------------------
        self.play(FadeOut(VGroup(f2, why)), run_time=0.7)
        feat = targets.pentagon()
        P, tr1, tr2 = feat.P, feat.run.phase1, feat.run.phase2
        frame = targets.principal_frame(P, feat.x)
        xs = tr1.x + tr2.x
        acts = tr1.action + tr2.action
        gs = tr1.gnorm + tr2.gnorm

        view = OrbitView(lambda ts: P.bodies(xs[0], ts), P.N, radius=1.85,
                         center=LEFT * 3.3 + DOWN * 0.3,
                         projector=Projector(2, frame), samples=720,
                         curve_width=3.0, fixed_scale=True)
        fA = Field("action", 12, color=T.CURVE).place(RIGHT * 1.3 + UP * 0.9)
        fG = Field("gradient", 12, color=T.COOL).place(RIGHT * 1.3 + UP * 0.35)
        fM = Field("Morse index", 6, color=T.BAD).place(RIGHT * 1.3 + DOWN * 0.2)
        fA.set("%.9f" % acts[0])
        fG.set(sci(gs[0]))
        fM.set("--")
        self.play(FadeOut(ep), FadeIn(view), FadeIn(fA), FadeIn(fG), FadeIn(fM),
                  run_time=1.1)

        n = len(xs)

        def drive(m, alpha):
            k = int(round((n - 1) * alpha ** 2.0))
            m.set_coefficients(lambda ts, xx=xs[k]: P.bodies(xx, ts), rescale=False)
            m.set_time(alpha * 2 * PI * 1.5)
            fA.set("%.9f" % acts[k])
            fG.set(sci(gs[k]))

        self.say_with("So hand it to the two phases and let them work.",
                      UpdateFromAlphaFunc(view, drive), run_time=11.0, rate_func=linear)
        self.say_with("The petals stay — the start put the orbit in the right family "
                      "— and everything else moves until the gradient is gone.",
                      UpdateFromAlphaFunc(view, drive), run_time=9.0, rate_func=linear)

        fM.set("%d" % feat.run.morse_index)
        self.play(FadeIn(fM), run_time=0.4)
        self.say("Action seventeen point four four seven, matching the catalogue to "
                 "every digit it stores.")
        self.say("And a Morse index of two: the action goes down in two directions "
                 "here, so no amount of rolling downhill would have found it.")

        sym = VGroup(
            M("<i>q</i>(<i>t</i> + 2π/5)  =  <i>R</i> <i>q</i>(<i>t</i>)", size=0.48,
              color=T.GOLD),
            B("every mode in the answer sits in one class modulo 5", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.28).move_to(RIGHT * 2.9 + DOWN * 1.6)
        self.play(FadeIn(sym), run_time=1.0)
        self.say_with("Afterwards the program detects the symmetry rather than "
                      "assuming it, and finds exactly the five-fold one the two arms "
                      "predicted.",
                      spin(view, 1.0), run_time=7.5)

        self.wipe(run_time=1.0)

        # ==================================================================
        # 2. the saddle only the second objective finds
        # ==================================================================
        head2 = B("two: an orbit that downhill cannot reach", size=T.SZ_HEAD,
                  color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(head2), run_time=0.8)

        sd = targets.saddle_n3()
        v2 = OrbitView(sd.bodies(), sd.P.N, radius=1.7, center=LEFT * 3.4 + UP * 0.1,
                       projector=Projector(2, sd.frame()), samples=600,
                       curve_width=3.0, fixed_scale=True)
        self.play(FadeIn(v2), run_time=0.9)

        obj = M("½ |∇<i>A</i>|%s" % sp("2"), size=0.56, color=T.COOL)
        obj.move_to(RIGHT * 3.2 + UP * 1.7)
        note = B("descend on this instead, and every critical point is a minimum of it",
                 size=T.SZ_SMALL, color=T.INK_DIM)
        note.next_to(obj, DOWN, buff=0.36)
        self.play(FadeIn(obj), FadeIn(note), run_time=1.0)
        self.say_with("Phase one has a second setting. Instead of descending on the "
                      "action, descend on the size of its gradient.",
                      spin(v2, 1.0), run_time=7.0)
        self.say_with("Every critical point of the action, whatever its index, is a "
                      "minimum of that — so downhill now works, and lands wherever "
                      "the gradient happens to vanish.",
                      spin(v2, 1.0), run_time=8.5)

        tbl = VGroup(
            Crow("120 random starts, N = 3 in the plane", size=T.SZ_TINY, color=T.INK_DIM),
            Crow("descending on the action     circle 82   eight 36   double cover 2",
              size=T.SZ_SMALL, color=T.INK),
            Crow("descending on |grad A|^2     circle 57   eight 35   two saddles",
              size=T.SZ_SMALL, color=T.INK),
        ).arrange(DOWN, buff=0.26, aligned_edge=LEFT)
        tbl.move_to(RIGHT * 2.6 + DOWN * 1.4)
        self.play(FadeIn(tbl), run_time=1.1)
        self.say_with("A hundred and twenty random starts each way. The action found "
                      "nothing but minima.",
                      spin(v2, 1.0), run_time=6.0)
        self.say_with("The gradient found this one, at action eleven point one five, "
                      "with a single downhill direction.",
                      spin(v2, 1.0), run_time=6.5, extra=0.4)

        self.wipe(run_time=1.0)

        # ==================================================================
        # 3. out of the plane
        # ==================================================================
        head3 = B("three: leaving the plane", size=T.SZ_HEAD,
                  color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(head3), run_time=0.8)

        hh = targets.hiphop()
        v3 = OrbitView(hh.bodies(), hh.P.N, radius=1.9, center=LEFT * 3.2 + DOWN * 0.1,
                       projector=Projector(3, hh.frame()), samples=720,
                       curve_width=2.8, fixed_scale=True)
        self.play(FadeIn(v3), run_time=0.9)

        st = VGroup(
            B("a rotating square, plus one oscillation across it", size=T.SZ_BODY,
              color=T.INK),
            M("modes 5 and 6", size=0.44, color=T.GOLD),
            B("chosen because 6/5 is close to a resonance of the square", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.32).move_to(RIGHT * 3.1 + UP * 0.4)
        self.play(FadeIn(st), run_time=1.1)

        self.say_with("The last one is the reason the program has structured starts at "
                      "all. Four bodies at the corners of a rotating square, pushed a "
                      "little out of their plane.",
                      tumble(v3, turns=1.0, sweep=2 * PI), run_time=9.0,
                      rate_func=linear)
        self.say_with("Push at the wrong rate and the wobble decays back into the "
                      "plane. Push near a rate the square resonates at, and it does "
                      "not.",
                      tumble(v3, turns=1.0, sweep=2 * PI), run_time=8.5,
                      rate_func=linear)
        self.say_with("Six to five is such a rate, and the orbit it converges on is "
                      "genuinely three-dimensional — the hip-hop, at action twenty-six "
                      "point seven six.",
                      tumble(v3, turns=1.0, sweep=2 * PI), run_time=9.0,
                      rate_func=linear, extra=0.4)

        self.play(FadeOut(VGroup(v3, st, head3)), *self.caption_anims(None),
                  run_time=1.2)
