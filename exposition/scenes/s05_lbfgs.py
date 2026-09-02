"""Phase one: rolling downhill, and watching a curve become an orbit."""

import numpy as np
from manim import (VGroup, VMobject, FadeIn, FadeOut, Create, UP, DOWN, LEFT,
                   RIGHT, ORIGIN, PI, Dot, linear, rate_functions,
                   UpdateFromAlphaFunc)

from expo import nbody, optim, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.plots import Plot, growing_line
from expo.surface import Wireframe, turn
from expo.readout import Field, sci
from expo.viz import OrbitView, Projector, spectrum_bars


class PhaseOne(ExpoScene):
    section_number = 5
    section_title = "phase one"

    def story(self):
        feat = targets.eight(seed=1, K=24)
        P, tr = feat.P, feat.run.phase1
        xs = tr.x
        frame = targets.principal_frame(P, feat.x)

        # — the stage ------------------------------------------------------
        view = OrbitView(lambda ts: P.bodies(xs[0], ts), P.N, radius=1.85,
                         center=LEFT * 3.4 + UP * 0.55,
                         projector=Projector(2, frame), samples=600,
                         curve_width=3.2, fixed_scale=True)
        view.place(radius=1.85, center=LEFT * 3.4 + UP * 0.55)

        gmax = max(tr.gnorm[0], 1e-12)
        plot = Plot((0, len(tr) - 1), (1e-8, 10 ** np.ceil(np.log10(gmax))),
                    width=4.9, height=2.15, center=RIGHT * 3.3 + DOWN * 0.55,
                    log_y=True, x_ticks=[0, len(tr) - 1],
                    y_ticks=[1e-8, 1e-4, 1e0], x_label="L-BFGS iteration")
        glab = C("size of the gradient", size=T.SZ_TINY, color=T.INK_DIM)
        glab.next_to(plot, UP, buff=0.30).align_to(plot, LEFT)
        gline, gset = growing_line(plot, range(len(tr)), tr.gnorm, color=T.COOL)

        power0 = np.einsum("mca,mca->m", P.shaped(xs[0]), P.shaped(xs[0]))
        bars = spectrum_bars(P.modes, power0, width=4.9, height=1.15,
                             label_every=4)
        bars.move_to(RIGHT * 3.3 + UP * 2.15)
        blab = C("power in each mode", size=T.SZ_TINY, color=T.INK_DIM)
        blab.next_to(bars, UP, buff=0.42).align_to(bars, LEFT)

        fA = Field("action", 11, color=T.CURVE).place(LEFT * 5.9 + DOWN * 1.55)
        fG = Field("gradient", 11, color=T.COOL).place(LEFT * 5.9 + DOWN * 2.05)
        fA.set("%.6f" % tr.action[0])
        fG.set(sci(tr.gnorm[0]))

        self.say_with("Every trial starts the same way: a random curve, made of a few "
                      "low modes, scaled to the size where the action is stationary.",
                      FadeIn(view), FadeIn(plot), FadeIn(glab), FadeIn(bars),
                      FadeIn(blab), FadeIn(fA), FadeIn(fG), run_time=1.6)
        self.add(gline)

        # — the descent ----------------------------------------------------
        n = len(tr)
        pmax = float(power0.max())

        def step(alpha):
            """Early iterations change everything, so give them most of the time."""
            return int(round((n - 1) * alpha ** 2.2))

        def drive(m, alpha):
            k = step(alpha)
            xk = xs[k]
            view.set_coefficients(lambda ts, xx=xk: P.bodies(xx, ts), rescale=False)
            view.set_time(alpha * 2 * PI * 2.0)
            gset(k + 1)
            fA.set("%.6f" % tr.action[k])
            fG.set(sci(tr.gnorm[k]))
            pw = np.einsum("mca,mca->m", P.shaped(xk), P.shaped(xk))
            for i, bar in enumerate(bars.bars):
                h = max(np.log10(max(pw[i], 1e-10)) - np.log10(1e-10), 0.0)
                h = max(h / 10.0 * 1.15, 0.006)
                bar.stretch_to_fit_height(h)
                bar.move_to([bar.get_center()[0], bars.axis.get_center()[1] + h / 2, 0])

        self.say_with("L-BFGS then rolls downhill. It never builds the matrix of second "
                      "derivatives; it infers the curvature it needs from the last ten "
                      "steps it took.",
                      UpdateFromAlphaFunc(view, drive), run_time=12.0,
                      rate_func=linear)

        self.say_with("The action falls, the gradient falls, and the spectrum fills in "
                      "from the low modes upward.",
                      UpdateFromAlphaFunc(view, drive), run_time=10.0,
                      rate_func=linear)

        self.wait(1.0)
        self.say("After a couple of hundred steps the curve is the figure eight, and "
                 "the gradient is down to about a millionth.")

        self.wipe(run_time=1.0)

        # — the same descent, drawn on the action itself --------------------

        # the two directions the search actually moved in, over its last stretch
        tail = np.array([P.flat(v) for v in xs[-70:]])
        D = tail - P.flat(feat.x)
        _u, _s, Vt = np.linalg.svd(D, full_matrices=False)
        u, v = Vt[0], Vt[1]
        proj = np.array([[d @ u, d @ v] for d in D])
        rad = float(np.abs(proj).max()) * 1.25

        a = np.linspace(-rad, rad, 45)
        xf = P.flat(feat.x)
        Z = np.array([[P.value(xf + ai * u + bj * v) for ai in a] for bj in a])

        wire = Wireframe(a, a, Z, width=6.4, depth=6.4, height=2.4,
                         center=DOWN * 0.30, lines=23, elev=0.58, azim=-0.85)
        lab = C("the action over the two directions the search actually moved in",
                size=T.SZ_TINY, color=T.INK_DIM)
        lab.to_corner(UP + LEFT, buff=0.45)
        self.say_with("It is worth seeing that descent the other way round: here is the "
                      "action itself, over the two directions the search moved in most.",
                      FadeIn(wire), FadeIn(lab), run_time=1.8)

        path_pts = [(float(pa), float(pb),
                     float(P.value(xf + pa * u + pb * v))) for pa, pb in proj]
        trail = VMobject(stroke_width=3.0)
        trail.set_stroke(T.CURVE)
        trail.set_fill(opacity=0)
        walker = Dot(wire.point(*path_pts[0]), radius=0.065, color=T.INK)

        def walk_path(m, alpha):
            k = max(2, int(round(alpha * (len(path_pts) - 1))) + 1)
            trail.set_points_as_corners([wire.point(*q) for q in path_pts[:k]])
            walker.move_to(wire.point(*path_pts[k - 1]))

        self.add(trail, walker)
        self.say_with("And here is where L-BFGS went: down the wall, and then a long "
                      "slow crawl along the floor of the valley.",
                      UpdateFromAlphaFunc(wire, walk_path), run_time=9.0,
                      rate_func=linear)

        def turn_all(m, alpha):
            e0, z0 = 0.58, -0.85
            m.set_camera(e0, z0 + 1.6 * alpha)
            trail.set_points_as_corners([m.point(*q) for q in path_pts])
            walker.move_to(m.point(*path_pts[-1]))

        self.say_with("That crawl is the whole problem with the first phase. The valley "
                      "floor is nearly flat, and a method that only knows which way is "
                      "downhill has almost nothing left to go on.",
                      UpdateFromAlphaFunc(wire, turn_all), run_time=11.0,
                      rate_func=rate_functions.ease_in_out_sine)

        self.play(FadeOut(VGroup(wire, trail, walker, lab)), run_time=1.0)

        self.say("So that is what a good case looks like. Two things are still wrong "
                 "with it.")

        # — what is wrong with it ------------------------------------------
        pts = [
            ("Downhill only finds minima.",
             "and almost no choreography is a minimum"),
            ("Convergence slows to a crawl.",
             "the last few digits would take thousands more steps"),
        ]
        rows = VGroup()
        for a, b in pts:
            rows.add(VGroup(B(a, size=T.SZ_BODY, color=T.INK),
                            B(b, size=T.SZ_SMALL, color=T.INK_DIM))
                     .arrange(DOWN, buff=0.20, aligned_edge=LEFT))
        rows.arrange(DOWN, buff=0.75, aligned_edge=LEFT).move_to(UP * 0.35)
        self.play(FadeIn(rows[0]), run_time=0.9)
        self.say("The first is fatal on its own: rolling downhill can only stop at the "
                 "bottom of a valley.")
        self.play(FadeIn(rows[1]), run_time=0.9)
        self.say("The second is merely expensive. Gradient methods converge linearly, "
                 "so each extra digit costs as much as the last one did.")
        self.say("Both are fixed by the same second phase.", extra=0.5)

        self.play(FadeOut(rows), *self.caption_anims(None), run_time=1.0)
