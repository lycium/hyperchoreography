#!/usr/bin/env python3
"""Render an animated plate to a movie through the film's own pipeline.

    ./make_plate_movie.py EndOrbitTurning plates/EndOrbit.mkv [--res 4k] [--ss 16]

Opaque, on the background colour, so the comp can fade it without an alpha channel.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from expo.pipeline import RESOLUTIONS, open_encoder, attach, close  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("scene", help="a class in plates.py")
    ap.add_argument("out", help="where to write it, relative to comp/")
    ap.add_argument("--res", default="4k", choices=sorted(RESOLUTIONS))
    ap.add_argument("--fps", type=float, default=30.0)
    ap.add_argument("--ss", type=int, default=16)
    a = ap.parse_args()

    W, H = RESOLUTIONS[a.res]
    out = os.path.join(HERE, a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)

    from manim import tempconfig

    proc = open_encoder(out, W, H, a.fps, W, H, codec="ffv1", depth=8)
    state = attach(proc, a.ss, "tent")
    settings = dict(
        pixel_width=W, pixel_height=H, frame_rate=a.fps * a.ss,
        format="mp4", write_to_movie=True, disable_caching=True,
        save_last_frame=False, preview=False, verbosity="ERROR",
        media_dir=os.path.join(ROOT, "media", "plates"), progress_bar="none",
    )
    t0 = time.time()
    with tempconfig(settings):
        import importlib
        plates = importlib.import_module("plates")
        getattr(plates, a.scene)().render()
    close(proc)

    dt = time.time() - t0
    print("  %s: %d supersampled frames -> %d frames, %.2f s in %.0f s"
          % (a.scene, state["frames"], state["out"], state["out"] / a.fps, dt))
    return 0


if __name__ == "__main__":
    sys.exit(main())