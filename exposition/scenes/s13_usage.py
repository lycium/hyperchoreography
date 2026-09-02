"""What you type, and what comes back."""

import numpy as np
from manim import (VGroup, FadeIn, FadeOut, UP, DOWN, LEFT, RIGHT, ORIGIN, PI,
                   linear)

from expo import catalog, theme as T
from expo.base import ExpoScene
from expo.mathtext import B, M, C, sb, sp, Crow
from expo.viz import OrbitView, Projector, PlaneGrid, spin, tumble


def annotated(pairs, size=None, buff=0.26):
    rows = VGroup()
    for flag, meaning in pairs:
        rows.add(VGroup(Crow(flag, size=size or T.SZ_SMALL, color=T.GOLD),
                        B(meaning, size=T.SZ_SMALL, color=T.INK_DIM))
                 .arrange(RIGHT, buff=0.5, aligned_edge=DOWN))
    rows.arrange(DOWN, buff=buff, aligned_edge=LEFT)
    for r in rows:
        r[0].align_to(rows[0][0], LEFT)
        r[1].align_to(rows, LEFT).shift(RIGHT * 3.6)
    return rows


class ParametersAndOutputs(ExpoScene):
    section_number = 13
    section_title = "using it"

    def story(self):
        self.say("So much for how it works. Here is what you actually type, and what "
                 "comes back.")

        cmd = C("hyperchoreography search --d 7 --N 10 --K 24 --omega g2:1,6",
                size=T.SZ_BODY, color=T.INK).move_to(UP * 2.6)
        self.play(FadeIn(cmd), run_time=1.0)

        rows = annotated([
            ("--d", "the dimension of the space to search in"),
            ("--N", "how many bodies"),
            ("--K", "how many Fourier modes the search carries"),
            ("--omega", "the rotating frame, if any"),
            ("--starts", "which families of starting guess to draw from"),
            ("--sym", "restrict the search to loops with a given symmetry"),
            ("--threads / --minutes", "how hard, and for how long"),
        ])
        rows.move_to(DOWN * 0.35)
        for r in rows:
            self.play(FadeIn(r, shift=RIGHT * 0.1), run_time=0.38)

        self.say("Dimension, number of bodies, and how many Fourier modes to carry. "
                 "Then the frame, the families of start to draw from, and how long to "
                 "run.")
        self.say("It writes a catalogue file as it goes and checkpoints every thirty "
                 "seconds. Stop it with control-C and run the same line again, and it "
                 "picks up at the trial it was on.")
        self.say("Trials are a deterministic function of the seed and the trial number, "
                 "so a resumed run is the same run — whatever order the threads "
                 "happen to finish in.")

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        three = annotated([
            ("--K", "too few modes and a close approach cannot be represented; "
                    "too many and every trial costs more"),
            ("--min-rigid", "the rigidity gate: how far from a rotating rigid "
                            "configuration a record has to be"),
            ("--min-deff", "throw away anything that does not fill enough dimensions, "
                           "before the expensive checks"),
        ], buff=0.55)
        three.move_to(ORIGIN)
        self.say_with("Three parameters are worth a sentence each, because they are "
                      "the ones that decide what you get.",
                      FadeIn(three[0]), run_time=1.2)
        for r in three[1:]:
            self.play(FadeIn(r), run_time=0.5)
        self.say("The number of modes is a trade: too few and an orbit with a close "
                 "approach cannot be written down at all, too many and every trial "
                 "costs more than it needs to.")
        self.say("The other two are filters. The rigidity gate throws away the "
                 "rotating rigid configurations, and the dimension filter throws away "
                 "flat orbits before the expensive checking starts.")

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        rec = catalog.load("d7_twist")
        view = OrbitView(rec.bodies, rec.N, radius=1.55,
                         center=LEFT * 4.0 + UP * 0.4,
                         projector=Projector(rec.d, rec.principal_frame()),
                         samples=720, curve_width=2.4, dot_radius=0.05,
                         fixed_scale=True)
        self.say_with("And here is one record, field by field.",
                      FadeIn(view), run_time=1.2)

        fields = annotated([
            ("action", "the value of A at the orbit, per body — its name, in practice"),
            ("energy", "fixed by the action through the virial theorem"),
            ("deff / d", "dimensions the motion uses, out of the ones searched"),
            ("morse", "how many directions the action goes down in"),
            ("nullity", "how many it is flat in"),
            ("ret_err", "how well the certified state satisfies the shift condition"),
            ("rigid", "distance from being a rotating rigid configuration"),
            ("twist", "a number that sees the symmetry group of the frame"),
            ("hits", "how many separate trials found this same orbit"),
        ], buff=0.24)
        fields.scale(0.86).move_to(RIGHT * 1.9 + UP * 0.1)
        for r in fields:
            self.play(FadeIn(r, shift=RIGHT * 0.08), run_time=0.30)
        self.say_with("Action and energy name the orbit. The two dimensions say how "
                      "much of the space it uses, against how much it was allowed.",
                      spin(view, 1.0), run_time=8.0)
        self.say_with("Morse index and nullity describe the landscape around it. The "
                      "return error is the certification. The last three are for "
                      "sorting and for sifting.",
                      spin(view, 1.0), run_time=9.0)

        self.wipe(run_time=1.0)

        # ------------------------------------------------------------------
        cmds = annotated([
            ("list", "everything in a file, sorted how you like"),
            ("show", "one record as JSON, coefficients and all"),
            ("export", "body positions over one period, as a CSV"),
            ("verify", "re-integrate and re-check, one record or the whole file"),
            ("refine", "redo the solve in arbitrary precision, to any number of digits"),
            ("prove", "turn a certified record into a theorem, by interval arithmetic"),
            ("merge", "union two catalogues, de-duplicating and re-applying the gates"),
            ("symmetry", "detect each loop's symmetry group rather than assume it"),
            ("continue", "follow a solution as a parameter is changed"),
        ], buff=0.30)
        cmds.move_to(ORIGIN)
        for r in cmds:
            self.play(FadeIn(r, shift=RIGHT * 0.08), run_time=0.32)
        self.say("The rest of the program is about reading catalogues back: listing "
                 "them, dumping a record, exporting the trajectories, re-verifying, "
                 "refining to as many digits as you like — and proving.")
        self.say("Merging is how a long search is spread over several machines: run "
                 "them with different seeds and union the results, and duplicates fold "
                 "together with their hit counts added.", extra=0.5)

        self.play(FadeOut(cmds), *self.caption_anims(None), run_time=1.0)
