"""Build the Allura Studio conform spec for the exposition film."""

import json
import os
import subprocess
import sys
import wave

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "out")
FPS = 30

TITLE_SECONDS = 0.0        # no title card; above zero puts one back
END_SECONDS = 14.0

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


def video_size(path):
    r = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v", "-show_entries",
         "stream=width,height", "-of", "csv=p=0", path],
        capture_output=True, text=True, check=True)
    return [int(x) for x in r.stdout.strip().split(",")[:2]]


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
    SIZE = video_size(os.path.join(OUT, manifest["scenes"][0]["video"]))

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

    mean_lufs = sum(x["lufs"] for x in scenes) / len(scenes)
    max_peak = max(x["peak"] for x in scenes)
    gain = round(min(10 ** ((-16.0 - mean_lufs) / 20.0), 0.985 / max_peak), 3)

    title_comp = None if TITLE_SECONDS <= 0 else {
        "name": "title",
        "size": SIZE, "fps": [FPS, 1], "duration": TITLE_SECONDS,
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
        {"footage": "plates/EndOrbit.mkv", "name": "orbit", "at": 0.3,
         "until": END_SECONDS - 0.4, "transfer": "linear",
         "fade_in": 1.8, "fade_out": 1.6},
    ]
    for i, plate in enumerate(["EndSubscribe", "EndCatalogue", "EndCode", "EndCredits"]):
        end_layers.append({"footage": f"plates/{plate}.png", "name": plate.lower(),
                           "at": 1.0 + 0.6 * i, "until": END_SECONDS - 0.8,
                           "hold": True, "transfer": "linear",
                           "fade_in": 1.1, "fade_out": 1.1})
    end_comp = {
        "name": "endcard",
        "size": SIZE, "fps": [FPS, 1], "duration": END_SECONDS,
        "layers": end_layers,
    }

    title_frames = int(round(TITLE_SECONDS * FPS))
    master_layers = [] if title_comp is None else [
        {"comp": "title", "name": "title", "at": 0.0, "until": title_frames / FPS}]
    cursor_f = title_frames
    marks = [{"t": 0.0, "text": CHAPTERS["s00_open"]}]
    for i, s in enumerate(scenes):
        at, until = cursor_f / FPS, (cursor_f + s["frames"]) / FPS
        master_layers.append({"footage": s["video"], "name": s["module"],
                              "at": at, "until": until, "transfer": "linear",
                              "gain": 0.0})
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
        "size": SIZE, "fps": [FPS, 1],
        "duration": (cursor_f + end_frames) / FPS,
        "markers": marks,
        "layers": master_layers,
    }

    spec = {"comps": [c for c in (title_comp, end_comp, master_comp) if c]}
    with open(os.path.join(HERE, "spec.json"), "w") as f:
        json.dump(spec, f, indent=1)

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
