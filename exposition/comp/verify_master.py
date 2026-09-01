#!/usr/bin/env python3
"""Verify the master against what the timeline says it must be.

Measured, not assumed: the frame count and duration; every scene's video
frames bit-identical to its source file at its exact placement (framemd5 of
the master, split by the timeline's ranges, against framemd5 of each scene);
every scene's audio equal to its narration WAV at its exact sample offset
under the spec's stated uniform gain -- bit-identical at gain 1, within one
code value of the quantised multiply otherwise (the overlapping tails are
all-zero, so a slice checks cleanly even where the next scene's head is
summed in); and silence in the gaps the timeline claims are silent (the
title and the end card).
"""

import json
import os
import subprocess
import sys
import wave

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MASTER = "out/master.mkv"
T = json.load(open("comp/timeline.json"))
FPS = T["fps"]
TITLE_F = int(round(T["title_seconds"] * FPS))
failures = 0


def check(name, ok, detail=""):
    global failures
    print(("  ok   " if ok else "  FAIL ") + name + (f"  ({detail})" if detail else ""))
    failures += 0 if ok else 1


def framemd5(path, select="0:v:0"):
    r = subprocess.run(["ffmpeg", "-v", "error", "-i", path, "-map", select,
                        "-f", "framemd5", "-"], capture_output=True, text=True, check=True)
    return [l.split(",")[-1].strip() for l in r.stdout.splitlines() if not l.startswith("#")]


# ---- geometry ------------------------------------------------------------------------------
master = framemd5(MASTER)
check("master frame count", len(master) == T["master_frames"],
      f"{len(master)} vs {T['master_frames']}")

# ---- every scene's picture, in place -------------------------------------------------------
cursor = TITLE_F
all_scenes = True
for s in T["scenes"]:
    src = framemd5("out/" + os.path.basename(s["video"]))
    seg = master[cursor:cursor + s["frames"]]
    if seg != src:
        all_scenes = False
        bad = sum(1 for a, b in zip(seg, src) if a != b)
        check(f"scene {s['module']} placed bit-exactly", False, f"{bad} frames differ")
    cursor += s["frames"]
check("all 15 scenes bit-exact at their placements", all_scenes)
check("title + end card fill the remainder",
      len(master) - (cursor - TITLE_F) - TITLE_F == int(round(T["end_seconds"] * FPS)))

# ---- every scene's sound, at its sample ----------------------------------------------------
r = subprocess.run(["ffmpeg", "-v", "error", "-i", MASTER, "-map", "0:a:0",
                    "-f", "s16le", "-"], capture_output=True, check=True)
pcm = r.stdout
expect_samples = T["master_frames"] * 1600
check("audio length matches the picture", len(pcm) == expect_samples * 2,
      f"{len(pcm) // 2} vs {expect_samples} samples")

# The gain every voice layer carries, from the spec itself; uniform by assertion.
spec = json.load(open("comp/spec.json"))
masterComp = next(c for c in spec["comps"] if c["name"] == "hyperchoreography")
gains = {l["gain"] for l in masterComp["layers"]
         if l.get("name", "").endswith(".voice")}
assert len(gains) == 1, gains
GAIN = gains.pop()

import array
cursor = TITLE_F * 1600
sound_ok = True
worst = 0
for s in T["scenes"]:
    with wave.open("out/" + os.path.basename(s["narration"])) as w:
        wav = array.array("h")
        wav.frombytes(w.readframes(w.getnframes()))
    got = array.array("h")
    got.frombytes(pcm[cursor * 2:(cursor + len(wav)) * 2])
    if GAIN == 1.0:
        ok = got == wav
    else:
        bad = 0
        for a, b in zip(wav, got):
            e = abs(b - a * GAIN)
            worst = max(worst, e)
            bad += e > 1.0 + 1e-9
        ok = bad == 0
    if not ok:
        sound_ok = False
        check(f"narration {s['module']} at its sample", False)
    cursor += s["frames"] * 1600
check(f"all 15 narrations exact under gain x{GAIN} (worst residual {worst:.2f} LSB)",
      sound_ok)
tail = pcm[cursor * 2:]
check("the end card is silent", all(b == 0 for b in tail),
      f"{len(tail)//2} samples")
check("the title is silent", all(b == 0 for b in pcm[:TITLE_F * 1600 * 2]))

# ---- the container -------------------------------------------------------------------------
r = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v", "-show_entries",
                    "stream=codec_name,pix_fmt,color_range", "-of", "csv=p=0", MASTER],
                   capture_output=True, text=True)
check("container is ffv1 bgr0 full-range", r.stdout.strip() == "ffv1,bgr0,pc", r.stdout.strip())

print()
if failures:
    print(f"FAIL: {failures} check(s) failed")
    sys.exit(1)
print(f"PASS: the master is {len(master)} frames = {len(master)/FPS:.2f} s, every scene's "
      "picture and sound bit-exact in place, silence where silence belongs")
