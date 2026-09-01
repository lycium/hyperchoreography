"""The stage: a scene base class with a narration band, titles and section marks.

Every line of explanation goes through `say()`, which holds the text on screen long
enough to be read aloud. The length of that hold is measured from the synthesised
narration when there is any, and estimated from theme.WORDS_PER_MINUTE when there
is not -- so the deck is paced by the reading rather than by guesses at animation
lengths.
"""

from __future__ import annotations

import re

from manim import (Scene, VGroup, Line, FadeIn, FadeOut, Create, Succession,
                   UP, DOWN, LEFT, RIGHT, ORIGIN)

from . import narrate
from . import theme as T
from .mathtext import B, C

_TAG = re.compile(r"<[^>]+>")


def visible_len(s: str) -> int:
    return len(_TAG.sub("", s))


def wrap_markup(text: str, width: int = 78):
    """Wrap on spaces, counting only visible characters, so Pango tags are free.

    Tags are assumed not to span a line break, which is true of everything the
    deck writes: emphasis is applied to single words or short phrases.
    """
    words, lines, cur, n = text.split(" "), [], [], 0
    for w in words:
        lw = visible_len(w)
        if cur and n + 1 + lw > width:
            lines.append(" ".join(cur))
            cur, n = [w], lw
        else:
            cur.append(w)
            n += (1 if n else 0) + lw
    if cur:
        lines.append(" ".join(cur))
    return lines


def spoken(text: str) -> str:
    """The line as it is read aloud: markup gone, and nothing else changed."""
    return _TAG.sub("", text)


def reading_time(text: str) -> float:
    """How long the line takes to read aloud, plus a breath.

    When a narration track is being made the answer is measured from the synthesised
    audio, so the film is paced by the reading rather than by an estimate of it.
    """
    plain = spoken(text)
    d = narrate.duration(plain)
    if d is None:
        words = max(1, len(plain.split()))
        d = words / T.WORDS_PER_MINUTE * 60.0
    return max(T.MIN_DWELL, d + T.BREATH)


class ExpoScene(Scene):
    """Shared furniture. Subclasses override `story()`, not `construct()`."""

    section_number = 0
    section_title = ""
    caption_width = 64
    caption_max_width = 84

    def construct(self):
        self.camera.background_color = T.BG
        self._caption = None
        self._mark = None
        self.cues = []
        if self.section_title:
            self.place_section_mark()
        self.story()
        self.wait(1.0)

    # -- narration --------------------------------------------------------
    def _build_caption(self, markup: str) -> VGroup:
        """The narration band, wrapped to at most three lines.

        The band grows upward from a fixed bottom edge, so an unusually long line
        would climb into the picture. Widening the wrap and, past a point, easing
        the size down keeps the top of the band where the scenes expect it.
        """
        n = visible_len(markup)
        width = min(self.caption_max_width,
                    max(self.caption_width, n // 3 + 2))
        lines = wrap_markup(markup, width)
        if len(lines) > 3:
            lines = wrap_markup(markup, self.caption_max_width)
        mobs = [B(l, size=T.SZ_CAPTION) for l in lines]
        g = VGroup(*mobs).arrange(DOWN, buff=0.16)
        if g.height > T.CAPTION_MAX_HEIGHT:          # never climb into the picture
            g.scale(T.CAPTION_MAX_HEIGHT / g.height)
        g.move_to([0, T.CAPTION_Y + g.height / 2 - 0.2, 0])
        return g

    def now(self) -> float:
        """Where we are on the scene clock, for lining the narration up."""
        return float(getattr(self.renderer, "time", 0.0))

    def caption_anims(self, markup: str | None):
        """Animations that swap the narration band, without waiting.

        The old line is faded out *before* the new one arrives rather than at the
        same time -- two lines dissolving through each other in the same place is
        unreadable. Returned as one animation so that a caller playing it alongside
        a long animation cannot stretch it: manim's `play(run_time=...)` overrides
        every animation it is given, so the swap has to be immune to that.
        """
        if markup:
            self.cues.append({"t": self.now() + 0.5, "text": spoken(markup)})
        out = []
        if self._caption is not None:
            out.append(FadeOut(self._caption, shift=DOWN * 0.10, run_time=0.40))
        self._caption = self._build_caption(markup) if markup else None
        if self._caption is not None:
            out.append(FadeIn(self._caption, shift=UP * 0.10, run_time=0.55))
        if not out:
            return []
        if len(out) == 1:
            return out
        return [Succession(*out)]

    def _play_timed(self, caption_anims, anims, run_time, rate_func=None,
                    lag: float = 0.0):
        """Play a caption swap alongside animations, each keeping its own length."""
        for a in anims:
            a.run_time = run_time
            if rate_func is not None:
                a.rate_func = rate_func
        group = list(caption_anims) + list(anims)
        if not group:
            return
        if lag:
            self.play(*group, lag_ratio=lag)
        else:
            self.play(*group)

    def say(self, markup: str, extra: float = 0.0):
        """Show a line and hold it for a spoken reading."""
        anims = self.caption_anims(markup)
        rt = 0.95
        if anims:
            self.play(*anims)
        self.hold(reading_time(markup) - rt + extra)

    def say_with(self, markup: str, *anims, run_time: float = 2.0, extra: float = 0.0,
                 lag: float = 0.0, rate_func=None):
        """Show a line while something happens, then hold for the rest of the reading."""
        caps = self.caption_anims(markup)
        self._play_timed(caps, anims, run_time, rate_func, lag)
        self.hold(reading_time(markup) - max(run_time, 0.95) + extra)

    def clear_caption(self):
        anims = self.caption_anims(None)
        if anims:
            self.play(*anims)

    def hold(self, t: float):
        """Wait, unless the animation already outlasted the line."""
        if t > 0.02:
            self.wait(t)

    def beat(self, t: float = 0.9):
        self.hold(t)

    # -- section furniture ------------------------------------------------
    def place_section_mark(self):
        num = C("%02d" % self.section_number, size=T.SZ_TINY, color=T.RULE)
        name = C(self.section_title.upper(), size=T.SZ_TINY, color=T.RULE)
        g = VGroup(num, name).arrange(RIGHT, buff=0.22)
        g.to_corner(UP + RIGHT, buff=0.42)
        self.add(g)
        self._mark = g

    def title_card(self, title: str, subtitle: str = None, hold: float = 2.6):
        t = B(title, size=T.SZ_TITLE, color=T.INK)
        parts = [t]
        if subtitle:
            s = B(subtitle, size=T.SZ_BODY, color=T.INK_DIM)
            parts.append(s)
        g = VGroup(*parts).arrange(DOWN, buff=0.42).move_to(ORIGIN)
        rule = Line(LEFT * 1.6, RIGHT * 1.6, stroke_width=1.6, color=T.RULE)
        rule.next_to(g, DOWN, buff=0.55)
        self.play(FadeIn(t, shift=UP * 0.2), run_time=1.0)
        if subtitle:
            self.play(FadeIn(parts[1], shift=UP * 0.1), Create(rule), run_time=0.9)
        self.wait(hold)
        self.play(FadeOut(VGroup(g, rule)), run_time=0.9)

    def heading(self, text: str, at=UP * 3.05):
        h = B(text, size=T.SZ_HEAD, color=T.GOLD)
        h.move_to(at)
        return h

    # -- housekeeping -----------------------------------------------------
    def wipe(self, *keep, run_time: float = 0.8):
        """Fade everything except the narration band and whatever is kept."""
        keepers = set(id(m) for m in keep)
        if self._caption is not None:
            keepers.add(id(self._caption))
        if self._mark is not None:
            keepers.add(id(self._mark))
        gone = [m for m in self.mobjects if id(m) not in keepers]
        if gone:
            self.play(*[FadeOut(m) for m in gone], run_time=run_time)

    def story(self):
        raise NotImplementedError
