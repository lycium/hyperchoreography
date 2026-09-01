"""Render one scene straight into ffmpeg. Invoked by render.py as a subprocess.

A subprocess per scene, because manim's configuration is global and a long deck
should not accumulate one scene's state into the next.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("scene")
    ap.add_argument("--out", required=True)
    ap.add_argument("--width", type=int, required=True)
    ap.add_argument("--height", type=int, required=True)
    ap.add_argument("--render-width", type=int, required=True)
    ap.add_argument("--render-height", type=int, required=True)
    ap.add_argument("--fps", type=float, required=True)
    ap.add_argument("--ss", type=int, default=4)
    ap.add_argument("--codec", default="ffv1")
    ap.add_argument("--depth", type=int, default=8)
    ap.add_argument("--shape", default="tent")
    ap.add_argument("--linear-light", action="store_true")
    ap.add_argument("--no-narrate", action="store_true")
    a = ap.parse_args()

    if a.no_narrate:
        os.environ["EXPO_NARRATE"] = "0"

    from manim import config, tempconfig
    from expo import pipeline

    from expo import narrate

    fps_hi = a.fps * a.ss
    silent = a.out + ".silent" + os.path.splitext(a.out)[1]
    proc = pipeline.open_encoder(
        silent, a.render_width, a.render_height, fps_hi, a.ss, a.fps,
        a.width, a.height, codec=a.codec, depth=a.depth,
        linear_light=a.linear_light, shape=a.shape)
    state = pipeline.attach(proc)

    settings = dict(
        pixel_width=a.render_width,
        pixel_height=a.render_height,
        frame_rate=fps_hi,
        format="mp4",
        write_to_movie=True,
        disable_caching=True,          # every animation must actually produce frames
        save_last_frame=False,
        preview=False,
        verbosity="ERROR",
        media_dir=os.path.join(ROOT, "media"),
        progress_bar="none",
    )

    t0 = time.time()
    with tempconfig(settings):
        mod = importlib.import_module("scenes." + a.module)
        cls = getattr(mod, a.scene)
        scene = cls()
        scene.render()
    pipeline.close(proc)

    dt = time.time() - t0
    n = state["frames"]
    total = n / fps_hi

    cues = getattr(scene, "cues", [])
    with open(a.out + ".cues.json", "w") as f:               # for the spoken script
        json.dump({"scene": a.scene, "module": a.module,
                   "duration": total, "cues": cues}, f, indent=1)
    voiced = False
    if cues and not a.no_narrate and narrate.enabled():
        # kept, not deleted: a comp wants the speech as its own layer rather than
        # cross-faded along with the picture
        wav = os.path.splitext(a.out)[0] + ".narration.wav"
        if narrate.build_track(cues, total, wav):
            lossless = a.codec in ("ffv1", "x264rgb", "utvideo", "prores")
            voiced = narrate.mux(silent, wav, a.out, lossless=lossless)
    if voiced:
        os.remove(silent)
    else:
        os.replace(silent, a.out)

    print("    %s: %d supersampled frames -> %.1f s%s in %.0f s (%.0f fps)"
          % (a.scene, n, total, " with %d narration cues" % len(cues) if voiced else "",
             dt, n / max(dt, 1e-9)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
