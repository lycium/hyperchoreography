"""From a certified orbit to a theorem: interval arithmetic, the validated flow, Krawczyk's test."""

import json
import os
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING

import numpy as np
from manim import (VGroup, VMobject, FadeIn, FadeOut, Create, Square, Polygon, Line, Dot,
                   Arrow, DashedLine, UpdateFromAlphaFunc, Transform,
                   UP, DOWN, LEFT, RIGHT, ORIGIN, PI, linear, rate_functions)

from expo import catalog, nbody, shoot, theme as T, viz
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, Crow
from expo.plots import Plot, growing_line
from expo.readout import Field, sci
from expo.viz import OrbitView, Projector, tumble

T3 = 2 * PI / 3


def load_proof(name):
    with open(os.path.join(catalog.DATA, name + ".json")) as f:
        return json.load(f)


def prepare_ellipse():
    """Body 0's position over one T/N, and how it moves with its own starting position:"""
    rec = catalog.load("eight")
    P = nbody.Action(rec.N, 2, int(rec.modes.max()), modes=rec.modes)
    x = np.array(rec.coef[:, :, :2])

    def build():
        Z0 = shoot.state_from_loop(P, x)
        N, d = P.N, 2
        nd = N * d
        Mn = 120

        def run(Z):
            p = Z[:nd].reshape(N, d).copy()
            v = Z[nd:].reshape(N, d).copy()
            out = [p[0].copy()]
            for _ in range(Mn):
                p, v = shoot.rk4(p, v, T3 / Mn, 8)
                out.append(p[0].copy())
            return np.array(out)

        base = run(Z0)
        eps = 1e-5
        J = np.zeros((Mn + 1, 2, 2))
        for c in range(2):
            Zp, Zm = Z0.copy(), Z0.copy()
            Zp[c] += eps
            Zm[c] -= eps
            J[:, :, c] = (run(Zp) - run(Zm)) / (2 * eps)
        return {"ts": np.linspace(0, T3, Mn + 1), "pos": base, "J": J}

    return rec, shoot.cached_run("prove_ellipse", build)


def interval_markup(lo: str, hi: str, places: int, tail_color: str):
    """[lo, hi] printed to `places` decimals, rounded outward so the printed interval still"""
    q = Decimal(1).scaleb(-places)
    a = format(Decimal(lo).quantize(q, rounding=ROUND_FLOOR), "f")
    b = format(Decimal(hi).quantize(q, rounding=ROUND_CEILING), "f")
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    span = '<span foreground="%s">%s</span>'
    return "[ %s%s ,  %s%s ]" % (a[:n], span % (tail_color, a[n:]), b[:n], span % (tail_color, b[n:]))


class GlowSquare(VGroup):
    """A box inside a soft halo, resizable in place."""

    REACH = 1.2

    def __init__(self, color, side: float, layers: int = viz.HALO_LAYERS):
        super().__init__()
        self.layers = layers
        self.side = side
        alphas = viz.halo_alphas(layers)
        for k in range(layers, 0, -1):
            s = Square(side_length=side * self._span(k / layers), stroke_width=0,
                       fill_color=color, fill_opacity=alphas[k - 1])
            self.add(s)
        self.core = Square(side_length=side, stroke_width=1.8, stroke_color=color,
                           fill_color=color, fill_opacity=0.30)
        self.add(self.core)

    def _span(self, f: float) -> float:
        return 1.0 + self.REACH * f

    def resize(self, side: float, at):
        """Scale the whole stack at once rather than stretching each layer."""
        if side != self.side:
            self.scale(side / self.side)
            self.side = side
        self.move_to(at)
        return self


class TheProof(ExpoScene):
    section_number = 9
    section_title = "existence"

    def story(self):
        pf = load_proof("prove_eight")
        pg = load_proof("prove_d11")
        rec, ell = prepare_ellipse()
        cool, gold, zero = T.COOL.to_hex(), T.GOLD.to_hex(), T.ZERO.to_hex()

        row1 = VGroup(M("|Φ%s(<i>Z</i>) − <i>S</i> <i>Z</i>|  ≤  1.2·10%s" % (sb("<i>T</i>/<i>N</i>"), sp("−15")),
                        size=0.50, color=T.COOL),
                      B("a measurement, in floating point", size=T.SZ_SMALL, color=T.INK_DIM)
                      ).arrange(DOWN, buff=0.22)
        row2 = VGroup(M("Φ%s(<i>Z</i>*) − <i>S</i> <i>Z</i>*  =  0" % sb("<i>T</i>/<i>N</i>"),
                        size=0.50, color=T.GOLD),
                      B("a theorem: such a <i>Z</i>* exists", size=T.SZ_SMALL, color=T.INK_DIM)
                      ).arrange(DOWN, buff=0.22)
        rows = VGroup(row1, row2).arrange(DOWN, buff=0.85).move_to(UP * 0.55)
        arrow = Arrow(row1.get_bottom() + DOWN * 0.08, row2.get_top() + UP * 0.08, buff=0,
                      stroke_width=2.0, color=T.RULE, max_tip_length_to_length_ratio=0.35)
        self.say_with("A residual of ten to the minus fifteen is a statement about a "
                      "computation. It is not a statement that an orbit exists. Between the "
                      "two sits every rounding error the computation made.",
                      FadeIn(row1), run_time=1.2)
        self.say_with("To get from the first line to the second, the arithmetic itself has "
                      "to change.", Create(arrow), FadeIn(row2), run_time=1.1)
        self.wipe(run_time=0.9)

        axis_y = 1.35
        axis = Line(LEFT * 5.2 + UP * axis_y, RIGHT * 5.2 + UP * axis_y, stroke_width=1.2,
                    color=T.RULE)

        def bar(x0, x1, y, color, h=0.34, op=0.22):
            g = VGroup()
            body = Polygon([x0, y - h / 2, 0], [x1, y - h / 2, 0], [x1, y + h / 2, 0],
                           [x0, y + h / 2, 0], stroke_width=0, fill_color=color,
                           fill_opacity=op)
            g.add(body)
            for x in (x0, x1):
                g.add(Line([x, y - h / 2, 0], [x, y + h / 2, 0], stroke_width=2.2, color=color))
            return g

        truth = Dot([-1.35, axis_y, 0], radius=0.06, color=T.GOLD)
        ia = bar(-2.6, -0.4, axis_y, T.COOL)
        la = C("[a]", size=T.SZ_SMALL, color=T.COOL).next_to(ia, UP, buff=0.18)
        lt = C("the true value", size=T.SZ_TINY, color=T.GOLD).next_to(truth, DOWN, buff=0.16)
        self.play(FadeIn(axis), FadeIn(ia), FadeIn(la), run_time=0.9)
        self.play(FadeIn(truth, scale=0.4), FadeIn(lt), run_time=0.7)
        self.say("Every number becomes a pair of bounds that provably contain the truth, "
                 "and every operation rounds its result outward. Whatever happens, the "
                 "truth never leaves the box.")

        ib = bar(1.1, 2.4, axis_y, T.ROSE)
        lb = C("[b]", size=T.SZ_SMALL, color=T.ROSE).next_to(ib, UP, buff=0.18)
        self.play(FadeIn(ib), FadeIn(lb), run_time=0.7)
        y2 = -0.55
        axis2 = Line(LEFT * 5.2 + UP * y2, RIGHT * 5.2 + UP * y2, stroke_width=1.2, color=T.RULE)
        eps = 0.09
        isum = bar(-2.6 + 1.1 - eps, -0.4 + 2.4 + eps, y2, T.CURVE)
        lsum = C("[a] + [b]", size=T.SZ_SMALL, color=T.CURVE).next_to(isum, UP, buff=0.18)
        rnd = VGroup(*[Line([x, y2 - 0.30, 0], [x, y2 + 0.30, 0], stroke_width=1.0,
                            color=T.INK_DIM).set_opacity(0.7)
                       for x in (-2.6 + 1.1, -0.4 + 2.4)])
        lr = C("rounded outward, one unit in the last place", size=T.SZ_TINY,
               color=T.INK_DIM).next_to(isum, DOWN, buff=0.20)
        self.play(FadeIn(axis2), run_time=0.4)
        self.play(Transform(VGroup(ia.copy(), ib.copy()), isum), FadeIn(lsum), run_time=1.3)
        self.play(FadeIn(rnd), FadeIn(lr), run_time=0.7)
        self.say("That is interval arithmetic. It is slow, it is pessimistic, and it is the "
                 "only arithmetic a theorem can be built on.")
        self.wipe(run_time=0.9)

        view = OrbitView(rec.bodies, rec.N, radius=2.05, center=LEFT * 3.4 + DOWN * 0.3,
                         projector=Projector(rec.d, rec.principal_frame()), samples=720,
                         curve_width=2.4, dot_radius=0.07)

        def screen(xy):
            return view._screen(np.array([[xy[0], xy[1], 0.0]]))[0]

        rects = VGroup()
        wb = np.array(pf["wbox"])
        for t0, h, xlo, xhi, ylo, yhi in wb:
            poly = Polygon(screen((xlo, ylo)), screen((xhi, ylo)), screen((xhi, yhi)),
                           screen((xlo, yhi)), stroke_width=1.2, stroke_color=T.COOL,
                           fill_color=T.COOL, fill_opacity=0.0)
            poly.set_stroke(opacity=0.0)
            rects.add(poly)
        trace = np.array(pf["box"])
        tt, ww = trace[:, 0], trace[:, 2]
        lw0, lw1 = np.log10(ww.min()), np.log10(ww.max())

        def side_of(w):
            return 0.09 + 0.55 * (np.log10(w) - lw0) / max(lw1 - lw0, 1e-9)

        glow = GlowSquare(T.GOLD, side_of(ww[0]))
        p0 = screen(rec.bodies([0.0])[0, 0, :2])
        glow.resize(side_of(ww[0]), p0)

        plot = Plot((0, T3), (10 ** (lw0 - 0.5), 10 ** (lw1 + 0.8)), width=4.4, height=2.3,
                    center=RIGHT * 3.45 + DOWN * 0.55, log_y=True,
                    x_ticks=[0, 1, 2], y_ticks=[1e-20, 1e-18, 1e-16],
                    x_label="time, over one period divided by three")
        gline, setter = growing_line(plot, tt, ww, color=T.COOL, width=2.6)
        fw = Field("width of the box", 10, color=T.COOL).place(RIGHT * 1.35 + UP * 1.75)
        fw.set(sci(ww[0]))
        stage = VGroup(rects, view, glow)
        self.play(FadeIn(view), FadeIn(glow), FadeIn(plot), FadeIn(fw), run_time=1.1)
        self.add(gline)

        wb_t0 = wb[:, 0]

        def march(m, a):
            t = a * T3
            view.set_time(t)
            w = np.exp(np.interp(t, tt, np.log(ww)))
            q = screen(rec.bodies([t])[0, 0, :2])
            glow.resize(side_of(w), q)
            k = int(np.searchsorted(wb_t0, t, side="right"))
            for i in range(k):
                if i < len(rects):
                    rects[i].set_fill(opacity=0.26)
                    rects[i].set_stroke(opacity=0.85)
            setter(int(np.searchsorted(tt, t, side="right")))
            fw.set(sci(w))

        self.say_with("The integrator that certified every orbit is a Taylor series "
                      "recurrence. Run unchanged on intervals it becomes a validated flow: "
                      "a box around the starting state, and a box that provably contains "
                      "wherever that state goes.",
                      UpdateFromAlphaFunc(stage, march), run_time=10.0, rate_func=linear)
        self.say("Each step first lays down a rough enclosure, the pale rectangles, inside "
                 "which the true solution is trapped for the whole step; the tail of the "
                 "series is bounded there. Nothing in the flow is estimated.")
        self.say("The tight box grows, here by four orders of magnitude over a third of a "
                 "period. The cure is not cleverness but precision: three hundred and "
                 "eighty-eight bits, so that a box which starts at ten to the minus twenty "
                 "can afford to.")

        self.play(FadeOut(rects), FadeOut(glow), FadeOut(plot), FadeOut(gline), FadeOut(fw),
                  run_time=0.9)
        view.set_time(0.0)
        r0 = 0.10
        th = np.linspace(0, 2 * PI, 72, endpoint=False)
        unit = np.stack([np.cos(th), np.sin(th)]) * r0
        ets, epos, eJ = ell["ts"], ell["pos"], ell["J"]

        def ellipse_pts(t):
            i = int(np.clip(np.searchsorted(ets, t), 1, len(ets) - 1))
            a = (t - ets[i - 1]) / (ets[i] - ets[i - 1])
            J = (1 - a) * eJ[i - 1] + a * eJ[i]
            c = (1 - a) * epos[i - 1] + a * epos[i]
            return (c[:, None] + J @ unit).T

        disc = VMobject(stroke_width=1.6)
        disc.set_stroke(T.ZERO, opacity=0.95)
        disc.set_fill(T.ZERO, opacity=0.32)
        box = VMobject(stroke_width=1.6)
        box.set_stroke(T.COOL, opacity=0.9)
        box.set_fill(T.COOL, opacity=0.06)

        def shape(t):
            pts = ellipse_pts(t)
            sp_ = [screen(p) for p in pts]
            disc.set_points_as_corners(sp_ + [sp_[0]])
            lo, hi = pts.min(axis=0), pts.max(axis=0)
            corners = [screen((lo[0], lo[1])), screen((hi[0], lo[1])), screen((hi[0], hi[1])),
                       screen((lo[0], hi[1]))]
            box.set_points_as_corners(corners + [corners[0]])

        shape(0.0)
        ghosts = VGroup()
        for _ in range(14):
            g_ = VMobject(stroke_width=1.0)
            g_.set_stroke(T.ZERO, opacity=0.0)
            g_.set_fill(T.ZERO, opacity=0.0)
            g_.set_points_as_corners([ORIGIN, ORIGIN])
            ghosts.add(g_)
        placed = [False] * len(ghosts)

        def leave_ghost(a):
            i = int(a * len(ghosts))
            if 0 < i <= len(ghosts) and not placed[i - 1]:
                placed[i - 1] = True
                pts = [screen(p) for p in ellipse_pts((i - 0.5) / len(ghosts) * T3)]
                g_ = ghosts[i - 1]
                g_.set_points_as_corners(pts + [pts[0]])
                g_.set_stroke(T.ZERO, opacity=0.35)
                g_.set_fill(T.ZERO, opacity=0.10)

        fs = Field("stretch of the cloud, long over short", 8, color=T.ZERO).place(RIGHT * 1.35 + UP * 1.75)
        fb = Field("area of its box, relative to the disc", 8, color=T.COOL).place(RIGHT * 1.35 + UP * 1.15)
        leg = VGroup(C("the cloud: every start within a small disc of the first body",
                       size=T.SZ_TINY, color=T.ZERO),
                     C("the box: the smallest axis-aligned one that holds it",
                       size=T.SZ_TINY, color=T.COOL)
                     ).arrange(DOWN, buff=0.16, aligned_edge=LEFT).move_to(RIGHT * 3.4 + DOWN * 0.9)
        area0 = None

        def sv(J):
            return np.linalg.svd(J, compute_uv=False)

        self.play(FadeIn(box), FadeIn(disc), FadeIn(fs), FadeIn(fb), FadeIn(leg), run_time=0.9)
        stage2 = VGroup(ghosts, view, box, disc)

        def march2(m, a):
            t = a * T3
            view.set_time(t)
            shape(t)
            leave_ghost(a)
            pts = ellipse_pts(t)
            i = int(np.clip(np.searchsorted(ets, t), 1, len(ets) - 1))
            s = sv(eJ[i])
            fs.set("%.1f" % (s[0] / max(s[1], 1e-12)))
            lo, hi = pts.min(axis=0), pts.max(axis=0)
            fb.set("%.1f" % ((hi[0] - lo[0]) * (hi[1] - lo[1]) / (PI * r0 * r0)))

        self.say_with("Why it grows. Perturb where the first body starts, in a small disc, "
                      "and watch where those starts end up: a cloud that stretches and "
                      "turns. A box is axis-aligned, so it must hold the whole tilted cloud, "
                      "and the next step must hold the box's image, not the cloud's.",
                      UpdateFromAlphaFunc(stage2, march2), run_time=11.0, rate_func=linear)
        self.say("That wrapping is the price of the simplest possible enclosure, and the "
                 "whole method is a bet that paying it in bits is cheaper than avoiding it "
                 "in code.")
        self.wipe(run_time=0.9)

        formula = M("<i>K</i>(<i>B</i>)  =  <i>Z</i>%s − <i>Y</i> <i>F</i>(<i>Z</i>%s) + (<i>I</i> − <i>Y</i> <i>DF</i>(<i>B</i>)) (<i>B</i> − <i>Z</i>%s)"
                    % (sb("0"), sb("0"), sb("0")), size=0.46, color=T.INK).move_to(UP * 2.45)
        self.say_with("The same recurrence, differentiated term by term, encloses the "
                      "derivative of the flow over the whole box: how the end state moves "
                      "when the start does.", FadeIn(formula), run_time=1.0)

        bc = LEFT * 3.3 + DOWN * 0.45
        Bsq = Square(side_length=2.9, stroke_width=2.2, stroke_color=T.GOLD, fill_color=T.GOLD,
                     fill_opacity=0.04).move_to(bc)
        lB = M("<i>B</i>", size=0.44, color=T.GOLD).next_to(Bsq, UP + LEFT, buff=-0.55)
        lB.shift(RIGHT * 0.25 + DOWN * 0.15)
        z0 = Dot(bc, radius=0.05, color=T.INK)
        lz = M("<i>Z</i>%s" % sb("0"), size=0.34, color=T.INK).next_to(z0, DOWN + RIGHT, buff=0.08)
        kside = 2.9 * 0.11
        Ksq = GlowSquare(T.COOL, kside)
        Ksq.resize(kside, bc + RIGHT * 0.07 + UP * 0.05)
        lK = M("<i>K</i>(<i>B</i>)", size=0.40, color=T.COOL).next_to(Ksq, RIGHT, buff=0.30)
        facts = VGroup()
        for i, (lab, val, col) in enumerate((("radius of B", pf["radius"], T.GOLD),
                                             ("|Y F(Z0)|", pf["newton"], T.INK),
                                             ("|I − Y DF(B)|", pf["kappa"], T.COOL))):
            y = 1.25 - 0.62 * i
            row = VGroup(Crow(lab, size=T.SZ_SMALL, color=T.INK_DIM),
                         Crow(sci(val, 1), size=T.SZ_SMALL, color=col))
            row[0].move_to(RIGHT * 0.9 + UP * y, aligned_edge=LEFT)
            row[1].move_to(RIGHT * 3.9 + UP * y, aligned_edge=LEFT)
            facts.add(row)
        verdict = M("<i>K</i>(<i>B</i>)  ⊂  int <i>B</i>   ⇒   exactly one zero in <i>B</i>",
                    size=0.40, color=T.GOOD).move_to(RIGHT * 3.3 + DOWN * 1.05)
        self.play(FadeIn(Bsq), FadeIn(lB), FadeIn(z0), FadeIn(lz), run_time=0.9)
        self.say_with("Then Krawczyk's test: apply one Newton step to the entire box at once. "
                      "If the image lands strictly inside the box, there is exactly one "
                      "solution in it. Every symbol in that formula is now an interval.",
                      FadeIn(Ksq), FadeIn(lK), run_time=1.4)
        for f in facts:
            self.play(FadeIn(f, shift=RIGHT * 0.08), run_time=0.45)
        self.play(FadeIn(verdict), run_time=0.8)
        self.say("For the eight the box has radius ten to the minus twenty, the Newton step "
                 "is ten to the minus forty-three, and the image is smaller than the box by "
                 "thirteen orders of magnitude. The small square is drawn a trillion times "
                 "too large.")
        self.wipe(run_time=0.9)

        left = ["time shift", "translation, twice", "rotation"]
        right = ["energy", "momentum, twice", "angular momentum"]
        hl = B("flat directions", size=T.SZ_SMALL, color=T.INK_DIM).move_to(LEFT * 3.2 + UP * 2.15)
        hr = B("conservation laws", size=T.SZ_SMALL, color=T.INK_DIM).move_to(RIGHT * 3.2 + UP * 2.15)
        pairs = VGroup()
        for i, (a, b) in enumerate(zip(left, right)):
            y = 1.25 - 0.85 * i
            ta = B(a, size=T.SZ_BODY, color=T.ZERO).move_to(LEFT * 3.2 + UP * y)
            tb = B(b, size=T.SZ_BODY, color=T.GOOD).move_to(RIGHT * 3.2 + UP * y)
            ln = DashedLine(ta.get_right() + RIGHT * 0.3, tb.get_left() + LEFT * 0.3,
                            stroke_width=1.4, color=T.RULE)
            pairs.add(VGroup(ta, tb, ln))
        foot = C("four singular directions of DF, four equations dropped, four laws that put them back",
                 size=T.SZ_TINY, color=T.INK_DIM).move_to(DOWN * 1.65)
        self.say_with("One obstacle. The shift condition is flat in four directions: slide "
                      "the orbit in time, translate it, rotate it. So its derivative is "
                      "singular, and there is nothing to invert.",
                      FadeIn(hl), *[FadeIn(p[0]) for p in pairs], run_time=1.2)
        self.play(FadeIn(hr), *[FadeIn(p[1]) for p in pairs], *[Create(p[2]) for p in pairs],
                  FadeIn(foot), run_time=1.2)
        self.say("The test runs on a slice across those directions with four equations "
                 "dropped, and the dropped four come back from conservation: the flow keeps "
                 "energy, momentum and angular momentum, and so does the shift, so agreeing "
                 "everywhere else forces agreement there.")
        self.say("That closure is a second, tiny interval Newton argument, and it is checked "
                 "too.")
        self.wipe(run_time=0.9)

        head = Crow("PROVEN    N = %d    d = %d    period 2π" % (pf["N"], pf["d"]),
                 size=T.SZ_BODY, color=T.GOOD)
        e_lo, e_hi = pf["energy"]
        a_lo, a_hi = pf["action"]
        cert = VGroup(
            head,
            Crow("energy  ∈  " + interval_markup(e_lo, e_hi, 20, cool), size=T.SZ_SMALL, color=T.INK),
            Crow("action  ∈  " + interval_markup(a_lo, a_hi, 20, cool), size=T.SZ_SMALL, color=T.INK),
            Crow("radius %s     contraction %s     closure %s     %d validated steps     %.1f s"
              % (sci(pf["radius"], 1), sci(pf["kappa"], 1), sci(pf["closure"], 1), pf["steps"],
                 pf["seconds"]), size=T.SZ_TINY, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.36, aligned_edge=LEFT).move_to(UP * 0.95)
        self.play(FadeIn(cert[0]), run_time=0.8)
        self.say_with("What comes out is a theorem with numbers in it: a choreography of "
                      "three bodies exists within ten to the minus twenty of the refined "
                      "state, unique there up to the symmetries, with its energy and its "
                      "action pinned to twenty digits.",
                      FadeIn(cert[1], shift=RIGHT * 0.08), FadeIn(cert[2], shift=RIGHT * 0.08),
                      FadeIn(cert[3]), run_time=1.8, lag=0.3)

        fr = pg["frame"] or {}
        rates = fr.get("rates", [])
        try:
            spec = rates[:2] + rates[3:] if len(rates) >= 3 and float(rates[2]) == float(rates[0]) + float(rates[1]) else rates
            spec = ",".join("%g" % float(r) for r in spec)
        except ValueError:
            spec = ",".join(rates)
        secs = pg["seconds"]
        took = "%d min %d s" % (int(secs // 60), int(round(secs % 60))) if secs >= 90 else "%.0f s" % secs
        words = "one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen nineteen twenty".split()
        num = lambda n: words[n - 1] if 1 <= n <= 20 else str(n)
        spoken = ("%s minutes" % num(int(round(secs / 60)))) if secs >= 90 else ("%s seconds" % num(int(round(secs))))
        o11 = catalog.load("d11")
        view11 = OrbitView(o11.bodies, o11.N, radius=1.6, center=LEFT * 4.5 + DOWN * 0.1,
                           projector=Projector(o11.d, o11.principal_frame()), samples=720,
                           curve_width=2.4, dot_radius=0.06, fixed_scale=True)
        cert11 = VGroup(
            Crow("PROVEN    N = %d    d = %d    frame g2:%s" % (pg["N"], pg["d"], spec),
                 size=T.SZ_SMALL, color=T.GOOD),
            Crow("energy  ∈  " + interval_markup(pg["energy"][0], pg["energy"][1], 18, cool),
                 size=T.SZ_TINY, color=T.INK),
            Crow("action  ∈  " + interval_markup(pg["action"][0], pg["action"][1], 18, cool),
                 size=T.SZ_TINY, color=T.INK),
            Crow("radius %s   contraction %s   closure %s   %d validated steps   %s"
                 % (sci(pg["radius"], 1), sci(pg["kappa"], 1), sci(pg["closure"], 1), pg["steps"], took),
                 size=T.SZ_TINY, color=T.INK_DIM),
            Crow("slice dimension %d, %d gauge generators; the same code, in eleven dimensions"
                 % (pg["slice"], pg["gauge"]), size=T.SZ_TINY, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.30, aligned_edge=LEFT).move_to(RIGHT * 2.45 + DOWN * 0.25)
        self.play(FadeOut(cert), run_time=0.8)
        self.say_with("In a turning frame the same argument goes through: the frame's rotation "
                      "becomes an interval matrix, and its commuting rotations are the gauge. "
                      "The eleven-dimensional orbit that closes this film is a theorem: "
                      "%s." % spoken,
                      FadeIn(view11), FadeIn(cert11, shift=RIGHT * 0.08),
                      tumble(view11, turns=1.0, sweep=2 * PI), run_time=10.0, rate_func=linear)
        self.say_with("Every record this has been run on has passed. Proving the rest of the "
                      "catalogue is machine time, not new ideas.",
                      tumble(view11, turns=1.0, sweep=2 * PI), run_time=7.0, rate_func=linear,
                      extra=0.6)

        self.play(FadeOut(VGroup(view11, cert11)), *self.caption_anims(None), run_time=1.2)
