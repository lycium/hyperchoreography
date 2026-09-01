"""A look at the structure underneath: redundancy, saturation, and the twist."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, Create, UP, DOWN, LEFT, RIGHT,
                   ORIGIN, PI, linear, rate_functions, UpdateFromAlphaFunc)

from expo import catalog, nbody, optim, targets, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp
from expo.plots import Plot
from expo.viz import OrbitView, Projector, PlaneGrid, eigen_strip, spin, tumble


class TheHigherStructure(ExpoScene):
    section_number = 13
    section_title = "the structure underneath"

    def story(self):
        self.say("The last part of this is the part the program is really about, and it "
                 "is worth at least a glimpse.")
        self.say("Almost every hard problem in it is the same problem: one orbit has "
                 "many descriptions, and the redundancy has to be taken out by hand, "
                 "carefully, in the right place.")

        # ==================================================================
        # 1. counting the free directions
        # ==================================================================
        h = B("counting the directions that mean nothing", size=T.SZ_HEAD,
              color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(h), run_time=0.8)

        rec = catalog.load("hiphop")
        P = nbody.Action(rec.N, rec.d, int(rec.modes.max()), modes=rec.modes)
        x = np.array(rec.coef)
        H = P.hessian(x)
        neg, zero, _lifted = optim.inertia_gauge(P, x, H)
        # the strip shows the *raw* spectrum, where the flat directions are still
        # sitting at zero; the two counts come from the careful computation, which
        # lifts them out rather than trying to recognise them by size
        raw = np.linalg.eigvalsh(H)
        strip = eigen_strip(raw, width=7.2, height=0.5, label=False,
                            window=40).move_to(UP * 1.1)
        cap = C("the 40 smallest second derivatives at one N = 4 orbit"
                "      Morse index %d, nullity %d" % (neg, zero),
                size=T.SZ_TINY, color=T.INK_DIM)
        cap.next_to(strip, DOWN, buff=0.32)
        self.say_with("We met the flat directions already: sliding the loop in time, "
                      "and turning it in space. There is one of the first and, in d "
                      "dimensions, d times d minus one over two of the second.",
                      FadeIn(strip), FadeIn(cap), run_time=1.2)

        count = VGroup(
            M("1  +  <i>d</i>(<i>d</i>−1)/2", size=0.50, color=T.ZERO),
            B("for this orbit, in three dimensions: 1 + 3 = 4", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.28).move_to(DOWN * 0.6)
        self.play(FadeIn(count), run_time=1.0)
        self.say("So a correct answer has a nullity of exactly four here, and a stored "
                 "nullity of zero is impossible: the time shift alone forces one.")

        prob = VGroup(
            B("but a genuinely soft direction and an exactly flat one look the same "
              "to a tolerance", size=T.SZ_BODY, color=T.INK),
            B("so the flat ones are constructed and lifted out of the spectrum, rather "
              "than recognised in it", size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.26).move_to(DOWN * 1.72)
        self.play(FadeIn(prob), run_time=1.1)
        self.say("Which is the whole difficulty. A direction that is nearly flat "
                 "because the orbit is nearly degenerate looks exactly like one that is "
                 "flat by construction.")
        self.say("So the program builds the flat directions explicitly, pushes them "
                 "out of the spectrum, and counts what is left.")
        self.say("In a turning frame only the rotations that commute with the frame "
                 "survive -- a smaller set, and not one you can find a generator at a "
                 "time.")

        self.wipe(run_time=1.0)

        # ==================================================================
        # 2. why the dimension saturates
        # ==================================================================
        h2 = B("why the dimension runs out", size=T.SZ_HEAD,
               color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(h2), run_time=0.8)

        rec11 = catalog.load("d11")
        grid = PlaneGrid(rec11, cols=3, panel=1.30, gap=0.40,
                         center=DOWN * 0.15, show_bodies=True)
        glab = C("twelve bodies, drawn in each of the orbit's own principal planes",
                 size=T.SZ_TINY, color=T.INK_DIM)
        glab.move_to(UP * 2.35)
        self.say_with("Here is a question with a sharp answer. Given N bodies, how many "
                      "dimensions can a choreography actually fill?",
                      FadeIn(grid), FadeIn(glab), run_time=1.4)
        self.say("This is the eleven-dimensional one, drawn plane by plane. Five of "
                 "them carry a curve and the last carries a line, and eleven is what "
                 "that adds up to.")
        self.play(FadeOut(VGroup(grid, glab)), run_time=0.9)

        bud = VGroup(
            M("<i>d</i>%s  ≤  2 ⌊<i>N</i>/2⌋" % sb("eff"), size=0.60, color=T.GOLD),
            B("two directions oscillating at the same rate only span a plane, so extra "
              "dimensions have to be bought with extra modes", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.34).move_to(UP * 1.1)
        self.play(FadeIn(bud), run_time=1.2)
        self.say("Twice N over two, rounded down. The reason is a counting argument: "
                 "two directions oscillating at the same rate together span only a "
                 "plane, so it is extra modes, not extra directions, that buy "
                 "dimension.")

        tbl = VGroup(C("N       4    5    6    7    8    9   10   12",
                       size=T.SZ_SMALL, color=T.INK_DIM),
                     C("limit   4    4    6    6    8    8   10   12",
                       size=T.SZ_SMALL, color=T.INK))
        tbl.arrange(DOWN, buff=0.24).move_to(DOWN * 0.5)
        self.play(FadeIn(tbl), run_time=1.0)
        self.say("The table it predicts was checked against the catalogue, and then "
                 "used the other way round: eleven dimensions needs at least twelve "
                 "bodies, so that is where the search for eleven was pointed.")

        strat = VGroup(
            B("and a second fact, which is why starts matter so much", size=T.SZ_BODY,
              color=T.INK),
            B("the set of loops lying in a smaller subspace is invariant under both "
              "phases: a flat orbit can be entered, never left", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.26).move_to(DOWN * 1.72)
        self.play(FadeIn(strat), run_time=1.1)
        self.say("And there is a companion fact with teeth. If a loop ever lies "
                 "exactly in a smaller subspace, neither phase can take it out again.")
        self.say("The action is symmetric under reflection in that subspace, so both "
                 "methods preserve it exactly.")
        self.say("Measured: kicking a converged orbit along its softest directions "
                 "raised the effective dimension in zero of eight hundred and sixteen "
                 "tries. Starts are everything.")

        self.wipe(run_time=1.0)

        # ==================================================================
        # 3. the twist
        # ==================================================================
        h3 = B("a number that sees the group", size=T.SZ_HEAD,
               color=T.GOLD).move_to(UP * 3.0)
        self.play(FadeIn(h3), run_time=0.8)

        jet = VGroup(
            M("<i>A</i>%s  =  (1/2π) ∮ <i>q</i> ∧ <i>q̇</i> ∧ … ∧ "
              "<i>q</i>%s d<i>t</i>" % (sb("k"), sp("(<i>k</i>−1)")), size=0.48),
            B("the jet moment: a k-form built from the loop and its derivatives",
              size=T.SZ_SMALL, color=T.INK_DIM),
        ).arrange(DOWN, buff=0.28).move_to(UP * 1.3)
        tw = VGroup(
            M("χ*  =  max over rotations of  ⟨ <i>A</i>%s , <i>R</i>·ψ ⟩" % sb("k"),
              size=0.48, color=T.GOLD),
            B("paired against the form that the group preserves", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.28).move_to(DOWN * 0.35)
        self.say_with("Last one. Once a search runs in a turning frame built out of an "
                      "exceptional group, you would like to know whether the orbits it "
                      "finds know anything about that group.",
                      FadeIn(jet), run_time=1.4)
        self.say("Build a k-form out of the loop and its first few derivatives, "
                 "integrated over a period. It rotates with the loop, so it is not an "
                 "invariant on its own.")
        self.play(FadeIn(tw), run_time=1.1)
        self.say("Then pair it against the differential form that the group leaves "
                 "alone, and maximise over rotations. Now it is invariant under "
                 "everything a record is allowed to be rewritten by, and it still sees "
                 "the group.")

        van = VGroup(
            B("in the complex Fourier basis this is a sum over mode triples that add "
              "up to zero", size=T.SZ_BODY, color=T.INK),
            B("so only resonant orbits can have it, spread-out modes dominate, and a "
              "non-zero value certifies at least k dimensions", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.26).move_to(DOWN * 1.72)
        self.play(FadeIn(van), run_time=1.2)
        self.say("Written out in the Fourier basis it collapses to a sum over sets of "
                 "modes that add up to zero, weighted by a Vandermonde determinant.")
        self.say("So it is non-zero only for orbits whose modes are in arithmetic "
                 "resonance -- and being non-zero certifies, on its own, that the orbit "
                 "fills at least k dimensions.")

        self.wipe(run_time=1.0)

        # -- the ladder ------------------------------------------------------
        rows = VGroup(C("d     k     group", size=T.SZ_TINY, color=T.INK_DIM))
        for d, k, g in (("3", "3", "SO(3)"), ("4", "2", "SU(2)·U(1)"),
                        ("6", "3", "SU(3)"), ("7", "3", "G2"),
                        ("8", "4", "Spin(7)"), ("10", "5", "SU(5)")):
            rows.add(C("%-5s %-5s %s" % (d, k, g), size=T.SZ_SMALL,
                       color=T.GOLD if g in ("G2", "Spin(7)") else T.INK))
        rows.arrange(DOWN, buff=0.24, aligned_edge=LEFT).move_to(LEFT * 3.0 + UP * 0.2)
        note = VGroup(
            B("seven is the only odd dimension with a proper reduction",
              size=T.SZ_BODY, color=T.INK),
            B("which is why the seven-dimensional case is the one that opened, and why "
              "nine and eleven are built on top of it", size=T.SZ_SMALL,
              color=T.INK_DIM),
        ).arrange(DOWN, buff=0.28).move_to(RIGHT * 2.6 + UP * 0.2)
        self.play(FadeIn(rows), run_time=1.0)
        self.say("There is one of these in each of a handful of dimensions, and the "
                 "list is short because the groups are.")
        self.play(FadeIn(note), run_time=1.0)
        self.say("Seven is the only odd dimension on it, and that is not a "
                 "coincidence of this program.")
        self.say("It is why seven was the dimension that opened, and why nine and "
                 "eleven were built on top of the same structure.", extra=0.5)

        self.play(FadeOut(VGroup(rows, note)), *self.caption_anims(None), run_time=1.1)
