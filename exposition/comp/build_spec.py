#!/usr/bin/env python3
"""Build the Allura Studio conform spec for the exposition film.

Reads out/manifest.json for the running order and the cue map, but takes every
number that places a frame or a sample from the FILES: the manifest's frame
counts are the planned manim timeline, and the encoded scenes run 1-2 frames
shorter (the tent filter's terminal windows). Writes:

  comp/spec.json      three comps -- title, endcard, hyperchoreography
  comp/chapters.txt   YouTube chapter list for the master
  comp/timeline.json  what was measured (frames, samples, loudness), for verify

The double-precision discipline the trims depend on: a layer's "until" and the
next layer's "at" are printed from the SAME repr of the SAME python double, so
the loader parses the SAME IEEE value on both sides and every boundary frame
lands in exactly one trim (DESIGN-media-and-audio, Conform).
"""

import json
import os
import subprocess
import sys
import wave

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "out")
FPS = 30

TITLE_SECONDS = 6.0
END_SECONDS = 14.0

# One human title per scene, for the chapter list.
CHAPTERS = {
    "s00_open": "Opening: three bodies on one curve",
    "s01_choreo": "The choreography constraint",
    "s02_action": "Least action",
    "s03_fourier": "The Fourier basis",
    "s04_landscape": "The landscape and the Morse index",
    "s05_lbfgs": "Phase one: L-BFGS",
    "s06_newton": "Phase two: the damped Newton step",
    "s07_hard": "Three hard orbits",
    "s08_certify": "Certification: the shooting Newton",
    "s09_prove": "Existence: the interval proof",
    "s10_sameorbit": "Covers, frames and the rigidity gate",
    "s11_starts": "Where to start: the N-gon spectrum",
    "s12_frame": "The rotating frame",
    "s13_usage": "The command line and a record",
    "s14_structure": "Redundancy, budgets and the twist",
    "s15_close": "What is in the catalogue",
}


def video_frames(path):
    r = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v", "-count_packets",
         "-show_entries", "stream=nb_read_packets", "-of", "csv=p=0", path],
        capture_output=True, text=True, check=True)
    return int(r.stdout.strip())


def wav_samples(path):
    with wave.open(path) as w:
        assert w.getframerate() == 48000 and w.getnchannels() == 1
        return w.getnframes()


def wav_peak(path):
    import array
    with wave.open(path) as w:
        a = array.array("h")
        a.frombytes(w.readframes(w.getnframes()))
    return max(abs(min(a)), max(a)) / 32768.0 if a else 0.0


def loudness(path):
    r = subprocess.run(
        ["ffmpeg", "-v", "info", "-i", path, "-filter_complex", "ebur128",
         "-f", "null", "-"], capture_output=True, text=True)
    for line in r.stderr.splitlines():
        if line.strip().startswith("I:") and "LUFS" in line:
            return float(line.split()[1])
    return None


def main():
    with open(os.path.join(OUT, "manifest.json")) as f:
        manifest = json.load(f)

    scenes = []
    for s in manifest["scenes"]:
        video = s["video"]
        frames = video_frames(os.path.join(OUT, video))
        samples = wav_samples(os.path.join(OUT, s["narration"]))
        scenes.append({
            "module": s["module"],
            "video": os.path.join("..", "out", video),
            "narration": os.path.join("..", "out", s["narration"]),
            "frames": frames,
            "planned_frames": s["frames"],
            "samples": samples,
            "lufs": loudness(os.path.join(OUT, s["narration"])),
            "peak": wav_peak(os.path.join(OUT, s["narration"])),
        })

    # One uniform narration gain, measured rather than chosen: lift toward -16 LUFS
    # integrated, clamped so the loudest track's true peak stays under -0.13 dBFS.
    # Uniform because the tracks sit within a fraction of an LU of one another and a
    # per-scene trim would repair a consistency they already have.
    mean_lufs = sum(x["lufs"] for x in scenes) / len(scenes)
    max_peak = max(x["peak"] for x in scenes)
    gain = round(min(10 ** ((-16.0 - mean_lufs) / 20.0), 0.985 / max_peak), 3)

    # ---- the title and the end card ---------------------------------------------------------
    title_comp = {
        "name": "title",
        "size": [1920, 1080], "fps": [FPS, 1], "duration": TITLE_SECONDS,
        "layers": [
            {"footage": "plates/Backdrop.png", "name": "backdrop", "at": 0.0,
             "until": TITLE_SECONDS, "hold": True, "transfer": "linear"},
            {"footage": "plates/TitleOrbit.png", "name": "orbit", "at": 0.0,
             "until": TITLE_SECONDS, "hold": True, "transfer": "linear",
             "fade_in": 1.8, "fade_out": 1.3},
            {"footage": "plates/TitleName.png", "name": "name", "at": 0.9,
             "until": TITLE_SECONDS - 0.4, "hold": True, "transfer": "linear",
             "fade_in": 1.2, "fade_out": 1.0},
        ],
    }
    end_layers = [
        {"footage": "plates/Backdrop.png", "name": "backdrop", "at": 0.0,
         "until": END_SECONDS, "hold": True, "transfer": "linear"},
        {"footage": "plates/EndOrbit.png", "name": "orbit", "at": 0.3,
         "until": END_SECONDS - 0.4, "hold": True, "transfer": "linear",
         "fade_in": 1.8, "fade_out": 1.6},
    ]
    for i, plate in enumerate(["EndSubscribe", "EndCatalogue", "EndCode", "EndAllura"]):
        end_layers.append({"footage": f"plates/{plate}.png", "name": plate.lower(),
                           "at": 1.0 + 0.6 * i, "until": END_SECONDS - 0.8,
                           "hold": True, "transfer": "linear",
                           "fade_in": 1.1, "fade_out": 1.1})
    end_comp = {
        "name": "endcard",
        "size": [1920, 1080], "fps": [FPS, 1], "duration": END_SECONDS,
        "layers": end_layers,
    }

    # ---- the master -------------------------------------------------------------------------
    # Scenes at their true lengths after the title; audio at each scene's own video start (the
    # narration tails run ~50 ms into the next scene's black head and are silent there); the
    # end card after the last scene. Gains stay 1.0: the fifteen tracks measure within 0.6 LU
    # of one another (timeline.json carries the numbers) and their true peaks reach -0.06 dBFS,
    # so any lift would clip -- and the whole track is a stand-in for a real recording anyway.
    title_frames = int(round(TITLE_SECONDS * FPS))
    master_layers = [{"comp": "title", "name": "title",
                      "at": 0.0, "until": title_frames / FPS}]
    cursor_f = title_frames
    marks = [{"t": 0.0, "text": CHAPTERS["s00_open"]}]
    for i, s in enumerate(scenes):
        at, until = cursor_f / FPS, (cursor_f + s["frames"]) / FPS
        master_layers.append({"footage": s["video"], "name": s["module"],
                              "at": at, "until": until, "transfer": "linear",
                              "gain": 0.0})  # narration comes from the WAV layers alone
        if i > 0:
            marks.append({"t": at, "text": CHAPTERS[s["module"]]})
        cursor_f += s["frames"]
    film_end = cursor_f / FPS
    cursor = title_frames
    for s in scenes:
        at = cursor / FPS
        master_layers.append({"footage": s["narration"], "name": s["module"] + ".voice",
                              "at": at, "until": at + s["samples"] / 48000,
                              "gain": gain})
        cursor += s["frames"]
    end_frames = int(round(END_SECONDS * FPS))
    master_layers.append({"comp": "endcard", "name": "endcard",
                          "at": film_end, "until": (cursor_f + end_frames) / FPS})
    marks.append({"t": film_end, "text": "The catalogue, and what is open"})
    master_comp = {
        "name": "hyperchoreography",
        "size": [1920, 1080], "fps": [FPS, 1],
        "duration": (cursor_f + end_frames) / FPS,
        "markers": marks,
        "layers": master_layers,
    }

    spec = {"comps": [title_comp, end_comp, master_comp]}
    with open(os.path.join(HERE, "spec.json"), "w") as f:
        json.dump(spec, f, indent=1)

    # ---- chapters (YouTube: first at 0:00, each >= 10 s) ------------------------------------
    lines = []
    for m in marks:
        t = int(round(m["t"]))
        lines.append("%d:%02d:%02d %s" % (t // 3600, t % 3600 // 60, t % 60, m["text"])
                     if t >= 3600 else "%d:%02d %s" % (t // 60, t % 60, m["text"]))
    with open(os.path.join(HERE, "chapters.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    report = {
        "fps": FPS,
        "title_seconds": TITLE_SECONDS, "end_seconds": END_SECONDS,
        "film_frames": cursor_f - title_frames,
        "master_frames": cursor_f + end_frames,
        "master_seconds": (cursor_f + end_frames) / FPS,
        "scenes": scenes,
    }
    with open(os.path.join(HERE, "timeline.json"), "w") as f:
        json.dump(report, f, indent=1)

    total = sum(s["frames"] for s in scenes)
    drift = sum(s["planned_frames"] for s in scenes) - total
    print(f"{len(scenes)} scenes, {total} film frames ({drift} short of the manifest's plan)")
    print(f"master: {report['master_frames']} frames = {report['master_seconds']:.2f} s")
    lufs = [s["lufs"] for s in scenes]
    print(f"narration: {min(lufs):.1f}..{max(lufs):.1f} LUFS, spread {max(lufs)-min(lufs):.1f} LU; "
          f"gain x{gain} -> ~{mean_lufs + 20*__import__('math').log10(gain):.1f} LUFS, "
          f"peak {max_peak * gain:.3f}")


if __name__ == "__main__":
    main()
