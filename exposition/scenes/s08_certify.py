"""From a critical point to a certified orbit."""

import numpy as np
from manim import (VGroup, VMobject, FadeIn, FadeOut, Create, Circle, Dot, Line,
                   UP, DOWN, LEFT, RIGHT, ORIGIN, PI, linear, rate_functions,
                   UpdateFromAlphaFunc)

from expo import catalog, nbody, shoot, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.plots import Plot
from expo.readout import Field, sci
from expo.viz import OrbitView, Projector, body_dot, spin


def prepare():
    """A deliberately under-resolved eight, and the shooting Newton that fixes it."""
    rec = catalog.load("eight")
    P = nbody.Action(rec.N, 2, int(rec.modes.max()), modes=rec.modes)
    x = np.array(rec.coef[:, :, :2])

    def build():
        pw = np.einsum("mca,mca->m", x, x)
        xt = np.zeros_like(x)
        for i in np.argsort(-pw)[:5]:
            xt[i] = x[i]
        Z0 = shoot.state_from_loop(P, xt)
        Zc, hist = shoot.shoot_newton(Z0, P.N, P.d, steps=4000, iters=6)
        return {"Z0": hist[0][0], "Zc": Zc,
                "res": np.array([h[1] for h in hist]),
                "Zs": np.array([h[0] for h in hist]),
                "xt": xt.reshape(-1)}

    return P, x, shoot.cached_run("certify_eight", build)


def arc(Z, N, d, T, steps=900):
    """Where the bodies actually go over one period divided by N."""
    nd = N * d
    p = Z[:nd].reshape(N, d).copy()
    v = Z[nd:].reshape(N, d).copy()
    out = [p.copy()]
    sub = max(1, steps // 180)
    for _ in range(180):
        p, v = shoot.rk4(p, v, T / N / 180.0, sub)
        out.append(p.copy())
    return np.array(out)                       # (frames, N, d)


class Certification(ExpoScene):
    section_number = 8
    section_title = "certification"

    def story(self):
        P, xfull, cache = prepare()

        self.say("At this point we have a curve whose gradient is zero. That is not "
                 "the same thing as a solution of the N-body problem, and the "
                 "difference matters.")

        pt = VGroup(
            B("the series was cut off after K modes", size=T.SZ_BODY, color=T.INK),
            B("the potential was integrated by a finite sum", size=T.SZ_BODY,
              color=T.INK),
        ).arrange(DOWN, buff=0.34).move_to(UP * 1.1)
        self.play(FadeIn(pt), run_time=1.0)
        self.say("Two approximations went in: the Fourier series was truncated, and "
                 "the integral in the action was replaced by a sum over sample points.")
        self.say("So the program stops trusting its own formulation and goes and "
                 "checks against the equations of motion themselves.")
        self.wipe(run_time=0.9)

        # ------------------------------------------------------------------
        self.say("Read the positions and the velocities of all N bodies at time zero "
                 "straight off the series. That is the only thing the ordinary "
                 "differential equation needs.")

        cond = M("Φ%s(<i>Z</i>)  =  <i>S</i> <i>Z</i>" % sb("<i>T</i>/<i>N</i>"),
                 size=0.58, color=T.GOLD).move_to(UP * 2.45)
        cw = B("integrate for one period over N, and every body should land exactly "
               "where the next one started", size=T.SZ_SMALL, color=T.INK_DIM)
        cw.next_to(cond, DOWN, buff=0.30)
        self.play(FadeIn(cond), FadeIn(cw), run_time=1.1)
        self.say("Then integrate forward by one period divided by N, and check the one "
                 "condition that makes a choreography a choreography: every body must "
                 "land exactly where the body in front of it started.")

        # -- the under-resolved attempt --------------------------------------
        N, d = P.N, P.d
        Z0, Zc = cache["Z0"], cache["Zc"]
        res = cache["res"]
        A0 = arc(Z0, N, d, 2 * PI)
        A1 = arc(Zc, N, d, 2 * PI)

        scale = 1.75 / max(np.abs(A1).max(), 1e-9)
        centre = LEFT * 3.3 + DOWN * 0.55

        def screen(p):
            return centre + np.array([p[0] * scale, p[1] * scale, 0.0])

        def panel(traj, Z):
            g = VGroup()
            trails = []
            for k in range(N):
                tr = VMobject(stroke_width=2.2)
                tr.set_stroke(T.body_color(k, N), opacity=0.7)
                tr.set_fill(opacity=0)
                trails.append(tr)
                g.add(tr)
            dots = VGroup(*[body_dot(T.body_color(k, N), r=0.075) for k in range(N)])
            targ = VGroup(*[Circle(radius=0.115, stroke_width=2.0,
                                   color=T.body_color((k + 1) % N, N))
                            for k in range(N)])
            p0 = Z[:N * d].reshape(N, d)
            for k in range(N):
                targ[k].move_to(screen(p0[(k + 1) % N]))
            g.add(targ, dots)
            g.trails, g.dots, g.targets, g.traj = trails, dots, targ, traj
            return g

        pan = panel(A0, Z0)
        fR = Field("shift residual", 10, color=T.BAD).place(RIGHT * 0.9 + UP * 0.6)
        fR.set(sci(res[0]))
        lab = C("hollow rings: where each body has to arrive",
                size=T.SZ_TINY, color=T.INK_DIM)
        lab.next_to(pan, DOWN, buff=1.15)
        self.play(FadeOut(VGroup(cond, cw)), FadeIn(pan), FadeIn(fR), FadeIn(lab),
                  run_time=1.2)

        def march(traj):
            def f(m, a):
                j = int(round(a * (len(traj) - 1)))
                for k in range(N):
                    pts = np.array([screen(p) for p in traj[:max(j + 1, 2), k]])
                    m.trails[k].set_points_smoothly(pts)
                    m.dots[k].move_to(screen(traj[j, k]))
            return f

        self.say_with("Here is a version of the figure eight cut down to five modes. "
                      "Integrate it and the bodies miss.",
                      UpdateFromAlphaFunc(pan, march(A0)), run_time=8.0,
                      rate_func=linear)
        self.say("Five parts in a hundred out. Small, but the equations of motion do "
                 "not grade on a curve.")

        # -- the shooting Newton ---------------------------------------------
        plot = Plot((0, len(res) - 1), (1e-15, 1.0), width=4.5, height=2.5,
                    center=RIGHT * 3.4 + DOWN * 1.0, log_y=True,
                    x_ticks=list(range(len(res))),
                    y_ticks=[1e-15, 1e-10, 1e-5, 1e0], x_label="shooting step")
        line = plot.line(range(len(res)), res, color=T.COOL, width=3.0, smooth=False)
        dots = VGroup(*[plot.marker(i, v, color=T.CURVE, r=0.055)
                        for i, v in enumerate(res)])
        self.play(FadeIn(plot), run_time=0.8)
        self.say("So correct the starting state until they do not. That is a Newton "
                 "method again, this time on the missing amount as a function of the "
                 "starting state.")
        self.play(Create(line), run_time=1.6)
        for i in range(len(res)):
            fR.set(sci(res[i]))
            self.play(FadeIn(dots[i], scale=0.5), run_time=0.30)
        self.say("Four steps, and the bodies arrive to fourteen digits.")

        self.play(UpdateFromAlphaFunc(pan, lambda m, a: None), run_time=0.1)
        pan2 = panel(A1, Zc)
        self.play(FadeOut(pan), FadeIn(pan2), run_time=0.8)
        self.say_with("Now the same picture with the corrected state.",
                      UpdateFromAlphaFunc(pan2, march(A1)), run_time=8.0,
                      rate_func=linear)

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        self.say("Two details of how this is done are worth pulling out.")

        rows = VGroup(
            VGroup(B("the series is re-read off the certified orbit",
                     size=T.SZ_BODY, color=T.INK),
                   B("over one T/N segment: exact, N times cheaper, and it never "
                     "accumulates the error a full period would", size=T.SZ_SMALL,
                     color=T.INK_DIM)).arrange(DOWN, buff=0.20, aligned_edge=LEFT),
            VGroup(B("where double precision stalls, the solve is redone in MPFR",
                     size=T.SZ_BODY, color=T.INK),
                   B("the figure eight goes to 1e-71 that way", size=T.SZ_SMALL,
                     color=T.INK_DIM)).arrange(DOWN, buff=0.20, aligned_edge=LEFT),
        ).arrange(DOWN, buff=0.70, aligned_edge=LEFT).move_to(UP * 0.6)
        self.play(FadeIn(rows[0]), run_time=0.9)
        self.say("First, the Fourier series that gets stored is not the one we started "
                 "with. It is read back off the certified orbit, over a single segment "
                 "of length T over N.")
        self.play(FadeIn(rows[1]), run_time=0.9)
        self.say("Second, when double precision cannot get there, the whole solve is "
                 "repeated in arbitrary precision and rounded back.")

        stat = VGroup(
            C("over about 1800 catalogued records", size=T.SZ_TINY, color=T.INK_DIM),
            C("median shift residual   2.2e-15", size=T.SZ_SMALL, color=T.GOOD),
            C("worst                   6.0e-11", size=T.SZ_SMALL, color=T.INK),
        ).arrange(DOWN, buff=0.24, aligned_edge=LEFT).move_to(DOWN * 1.72)
        self.play(FadeIn(stat), run_time=1.0)
        self.say("Across the whole catalogue the typical orbit satisfies the shift "
                 "condition to fifteen digits.")
        self.say("The worst record in the file manages eleven.", extra=0.5)

        self.play(FadeOut(VGroup(rows, stat)), *self.caption_anims(None), run_time=1.1)
