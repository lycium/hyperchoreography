"""Recognising the same orbit in a different disguise."""

import numpy as np
from manim import (VGroup, VMobject, FadeIn, FadeOut, Create, Transform, UP, DOWN,
                   LEFT, RIGHT, ORIGIN, PI, linear, rate_functions,
                   UpdateFromAlphaFunc)

from expo import catalog, nbody, targets, trial, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.plots import Plot
from expo.readout import Field, sci
from expo.viz import OrbitView, Projector, spin, tumble


class TheSameOrbitTwice(ExpoScene):
    section_number = 9
    section_title = "the same orbit twice"

    def story(self):
        self.say("A search runs the same trial tens of thousands of times, so the same "
                 "orbit turns up again and again wearing different clothes. Four "
                 "disguises have to be seen through.")

        # ==================================================================
        # 1. covers
        # ==================================================================
        h = B("one: going round twice", size=T.SZ_HEAD, color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(h), run_time=0.8)

        P = nbody.Action(3, 2, 24)
        runs = []
        for m in (1, 2, 4):          # k must be coprime to N, so 3 is not a cover here
            r = trial.run(N=3, d=2, K=24, phase1="action", n1=100, newton=60,
                          start=targets.circle_start(P, m))
            runs.append(r)

        views, labs = VGroup(), VGroup()
        for i, r in enumerate(runs):
            v = OrbitView(lambda ts, rr=r: rr.P.bodies(rr.x, ts), 3, radius=1.15,
                          center=LEFT * 4.0 + RIGHT * 4.0 * i + UP * 0.35,
                          projector=Projector(2, targets.principal_frame(r.P, r.x)),
                          samples=600, curve_width=2.8, dot_radius=0.07,
                          fixed_scale=True)
            lab = C("action %.9f" % r.phase2.action[-1], size=T.SZ_TINY,
                    color=T.INK_DIM)
            lab.move_to(v.center + DOWN * 1.75)
            views.add(v)
            labs.add(lab)

        self.play(FadeIn(views[0]), FadeIn(labs[0]), run_time=0.9)
        self.say_with("Here is the simplest choreography there is: three bodies at the "
                      "corners of a rotating triangle.",
                      spin(views[0], 1.0), run_time=6.5)
        self.play(FadeIn(views[1]), FadeIn(labs[1]), FadeIn(views[2]), FadeIn(labs[2]),
                  run_time=1.0)
        self.say_with("And here it is going round twice in a period, and four times. "
                      "Every one of those is a genuine critical point, and every one is "
                      "the same circle.",
                      *[spin(v, 1.0) for v in views], run_time=8.0)

        law = M("<i>A</i>(<i>k</i> times round)  =  <i>k</i>%s <i>A</i>"
                % sp("2/3"), size=0.50, color=T.GOLD).move_to(DOWN * 1.95)
        self.play(FadeIn(law), run_time=0.9)
        self.say_with("The actions are in exactly the ratio k to the two thirds, which "
                      "is Kepler's third law wearing a different hat.",
                      *[spin(v, 1.0) for v in views], run_time=7.5)
        self.say("They are spotted from the highest common factor of the modes that "
                 "carry any power, and unwound before anything else happens.")

        self.wipe(run_time=1.0)

        # ==================================================================
        # 2. the frame
        # ==================================================================
        h2 = B("two: the same loop seen from a different angle", size=T.SZ_HEAD,
               color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(h2), run_time=0.8)

        rec = catalog.load("eight")
        P3 = nbody.Action(3, 3, 24)
        x3 = targets.embed(rec, P3)
        R = np.eye(3)
        th = 0.9
        R[[0, 2]] = [[np.cos(th), 0, -np.sin(th)], [np.sin(th), 0, np.cos(th)]]
        xt = np.einsum("ab,mcb->mca", R, P3.shaped(x3))

        v = OrbitView(lambda ts: P3.bodies(xt, ts), 3, radius=1.8,
                      center=LEFT * 3.2 + DOWN * 0.1, projector=Projector(3),
                      samples=600, curve_width=3.0, fixed_scale=True)
        self.play(FadeIn(v), run_time=0.9)
        self.say_with("A search in three dimensions can find the figure eight lying at "
                      "any angle at all, and every one of them is a different list of "
                      "numbers.",
                      tumble(v, turns=1.0, sweep=PI), run_time=8.0, rate_func=linear)

        pf = M("principal axes of  <i>X</i>%s<i>X</i>" % sp("T"), size=0.48,
               color=T.GOLD).move_to(RIGHT * 3.2 + UP * 0.9)
        note = B("rotate the coefficients into the frame the loop itself picks out",
                 size=T.SZ_SMALL, color=T.INK_DIM)
        note.next_to(pf, DOWN, buff=0.34)
        self.play(FadeIn(pf), FadeIn(note), run_time=1.0)

        base = Projector(3, targets.principal_frame(P3, xt))

        def straighten(m, a):
            p = Projector(3, base.R)
            m.proj = Projector(3, (1 - a) * v0.R + a * base.R)
            m.redraw()

        v0 = Projector(3, v.proj.R)
        self.say_with("So every record is turned into the frame its own motion picks "
                      "out, before anything is compared to anything.",
                      UpdateFromAlphaFunc(v, straighten), run_time=4.0,
                      rate_func=rate_functions.ease_in_out_sine)

        deff = VGroup(
            M("<i>d</i>%s  =  how many of those axes it actually uses" % sb("eff"),
              size=0.44),
            B("the figure eight found in three dimensions is flat, and says so",
              size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.30).move_to(RIGHT * 3.1 + DOWN * 1.5)
        self.play(FadeIn(deff), run_time=1.0)
        self.say_with("Counting how many of those axes carry anything gives the number "
                      "the catalogue calls the effective dimension -- which is a "
                      "property of the orbit, not of the search that found it.",
                      spin(v, 1.0), run_time=9.0)

        self.wipe(run_time=1.0)

        # ==================================================================
        # 3. Procrustes
        # ==================================================================
        h3 = B("three: shifted, reflected, relabelled", size=T.SZ_HEAD,
               color=T.GOLD).move_to(UP * 3.0)
        rows = VGroup(
            B("start it at a different moment", size=T.SZ_BODY, color=T.INK),
            B("run it backwards", size=T.SZ_BODY, color=T.INK),
            B("rotate or reflect it", size=T.SZ_BODY, color=T.INK),
            B("renumber the bodies", size=T.SZ_BODY, color=T.INK),
        ).arrange(DOWN, buff=0.34).move_to(UP * 0.5)
        self.play(FadeIn(h3), run_time=0.7)
        for r in rows:
            self.play(FadeIn(r, shift=RIGHT * 0.12), run_time=0.45)
        self.say("Four more ways to write the same orbit down. Two records are the same "
                 "choreography when some combination of these turns one into the other.")

        proc = VGroup(
            B("cheap invariants first, then a Procrustes fit", size=T.SZ_BODY,
              color=T.GOLD),
            B("the best rotation between two sampled loops is a singular value "
              "decomposition away", size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.28).move_to(DOWN * 1.72)
        self.play(FadeIn(proc), run_time=1.0)
        self.say("Candidates are filtered on invariants that are quick to compute, and "
                 "the survivors are compared properly: for a fixed time shift the best "
                 "rotation between two loops has a closed form.")

        self.wipe(run_time=1.0)

        # ==================================================================
        # 4. relative equilibria
        # ==================================================================
        h4 = B("four: the ones that are not really moving", size=T.SZ_HEAD,
               color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(h4), run_time=0.8)

        circ = trial.run(N=3, d=2, K=24, phase1="action", n1=100, newton=60,
                         start=targets.circle_start(P, 1))
        e8 = targets.eight(seed=1)

        pair = []
        for i, (r, name) in enumerate(((circ, "a rotating triangle"),
                                       (e8, "the figure eight"))):
            v = OrbitView(lambda ts, rr=r: rr.P.bodies(rr.x, ts), 3, radius=1.15,
                          center=LEFT * 3.6 + RIGHT * 2.6 * i + UP * 0.55,
                          projector=Projector(2, targets.principal_frame(r.P, r.x)),
                          samples=600, curve_width=2.6, dot_radius=0.065,
                          fixed_scale=True)
            lab = C(name, size=T.SZ_TINY, color=T.INK_DIM)
            lab.move_to(v.center + DOWN * 1.6)
            pair.append((v, lab, r))
            self.play(FadeIn(v), FadeIn(lab), run_time=0.6)

        # the mutual-distance ribbon
        ts = np.linspace(0, 2 * PI, 400)
        plot = Plot((0, 2 * PI), (0.0, 2.2), width=4.6, height=2.4,
                    center=RIGHT * 3.4 + UP * 0.25, x_ticks=[0], y_ticks=[0, 1, 2],
                    x_label="one period")
        lines = VGroup()
        for i, (v, lab, r) in enumerate(pair):
            q0 = r.P.curve(r.x, ts)
            q1 = r.P.curve(r.x, ts + 2 * PI / 3)
            dist = np.linalg.norm(q0 - q1, axis=1)
            lines.add(plot.line(ts, dist,
                                color=T.GOOD if i == 0 else T.CURVE, width=2.8))
        dlab = C("distance between two of the bodies", size=T.SZ_TINY, color=T.INK_DIM)
        dlab.next_to(plot, UP, buff=0.30).align_to(plot, LEFT)
        self.play(FadeIn(plot), FadeIn(dlab), Create(lines[0]), Create(lines[1]),
                  run_time=1.8)

        self.say_with("A rotating polygon is a choreography, technically. Every body "
                      "does follow the same circle. But nothing about the shape ever "
                      "changes.",
                      *[spin(v, 1.0) for v, _, _ in pair], run_time=9.0)
        self.say_with("Measure the distance between two bodies and it is a flat line. "
                      "For the figure eight it is not.",
                      *[spin(v, 1.0) for v, _, _ in pair], run_time=8.0)

        gate = VGroup(
            B("the rigidity defect: how far from flat that line is", size=T.SZ_BODY,
              color=T.INK),
            C("clustered below 2e-11, then a 4300-fold gap, then real orbits at 1.7e-2",
              size=T.SZ_TINY, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.26).move_to(DOWN * 1.95)
        self.play(FadeIn(gate), run_time=1.0)
        self.say("The catalogue stores that number for every record, and it comes out "
                 "in two clumps with an enormous empty gap between them.")
        self.say("So the threshold that throws the rigid ones away is not a judgement "
                 "call.")
        self.say("Which matters, because in high dimensions those rigid solutions "
                 "outnumber everything else.", extra=0.5)

        self.play(FadeOut(VGroup(*[v for v, _, _ in pair], *[l for _, l, _ in pair],
                                 plot, lines, dlab, gate, h4)),
                  *self.caption_anims(None), run_time=1.2)
