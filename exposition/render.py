"""Render the presentation."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from expo.manifest import SCENES
from expo.pipeline import CODECS, RESOLUTIONS, ffmpeg_binary


def write_manifest(parts, out_dir, name, width, height, fps, codec):
    """One file describing the assets, for whatever assembles them next."""
    import json
    from expo import narrate
    scenes = []
    for p in parts:
        cue_path = p + ".cues.json"
        if not os.path.exists(cue_path):
            continue
        with open(cue_path) as f:
            d = json.load(f)
        wav = os.path.splitext(p)[0] + ".narration.wav"
        if not os.path.exists(wav) and d["cues"] and narrate.enabled():
            narrate.build_track(d["cues"], d["duration"], wav)
        scenes.append({
            "module": d["module"], "scene": d["scene"],
            "video": os.path.basename(p),
            "narration": os.path.basename(wav) if os.path.exists(wav) else None,
            "duration": d["duration"],
            "frames": int(round(d["duration"] * fps)),
            "cues": d["cues"],
        })
    man = {
        "name": name, "width": width, "height": height, "fps": fps,
        "codec": codec, "colour": "lossless 8-bit RGB, full range, no subsampling",
        "opens_and_closes_on_black": True,
        "total_duration": sum(s["duration"] for s in scenes),
        "scenes": scenes,
    }
    path = os.path.join(out_dir, "manifest.json")
    with open(path, "w") as f:
        json.dump(man, f, indent=1)
    return path


def write_script(parts, out_dir, name):
    """Collect the cue files into one readable narration script."""
    import json
    lines = ["# " + name, "",
             "The spoken text, in order, with the time within each scene at which "
             "it appears.", ""]
    total = 0.0
    for p in parts:
        cue_path = p + ".cues.json"
        if not os.path.exists(cue_path):
            continue
        with open(cue_path) as f:
            d = json.load(f)
        lines.append("## %s  (%s, %.0f s)" % (d["scene"], d["module"], d["duration"]))
        lines.append("")
        for c in d["cues"]:
            lines.append("- `%5.1f`  %s" % (c["t"], c["text"]))
        lines.append("")
        total += d["duration"]
    lines.append("---")
    lines.append("")
    lines.append("Total running time: %d minutes %d seconds."
                 % (int(total) // 60, int(total) % 60))
    path = os.path.join(out_dir, "narration.md")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def scene_seconds(out_dir, mod, default=120.0):
    """How long this scene ran last time, so the long ones can be started first."""
    import json
    try:
        with open(os.path.join(out_dir, mod + ".mkv.cues.json")) as f:
            return float(json.load(f)["duration"])
    except (OSError, ValueError, KeyError):
        return default


def run_jobs(todo, jobs: int, out_dir: str) -> int:
    """Render the scenes, `jobs` at a time, longest first."""
    order = sorted(todo, key=lambda t: -scene_seconds(out_dir, t[0]))
    running, failed = {}, []
    while order or running:
        while order and len(running) < jobs:
            mod, dst, cmd = order.pop(0)
            print("  %-16s -> %s" % (mod, os.path.basename(dst)), flush=True)
            running[subprocess.Popen(cmd, cwd=HERE)] = mod
        done = [p for p in running if p.poll() is not None]
        if not done:
            time.sleep(0.5)
            continue
        for p in done:
            mod = running.pop(p)
            if p.returncode != 0:
                print("  ** %s failed (%d)" % (mod, p.returncode))
                failed.append(mod)
    return 1 if failed else 0


def concat(parts, out_path):
    """Join the finished scenes without re-encoding a single pixel."""
    lst = out_path + ".txt"
    with open(lst, "w") as f:
        for p in parts:
            f.write("file '%s'\n" % os.path.abspath(p).replace("'", r"'\''"))
    subprocess.run([ffmpeg_binary(), "-y", "-hide_banner", "-loglevel", "error",
                    "-f", "concat", "-safe", "0", "-i", lst,
                    "-c", "copy", out_path], check=True)
    os.remove(lst)


def transcode(src, dst, crf=16):
    subprocess.run([ffmpeg_binary(), "-y", "-hide_banner", "-loglevel", "error",
                    "-i", src, "-c:v", "libx264", "-crf", str(crf),
                    "-preset", "slow", "-pix_fmt", "yuv420p",
                    "-c:a", "aac", "-b:a", "192k",
                    "-movflags", "+faststart", dst], check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scenes", nargs="*", help="scene modules to render (default: all)")
    ap.add_argument("--res", default="1080p", choices=sorted(RESOLUTIONS),
                    help="output resolution (default 1080p)")
    ap.add_argument("--fps", type=float, default=30.0, help="output frame rate")
    ap.add_argument("--ss", type=int, default=4,
                    help="temporal supersampling: manim renders at fps*ss (default 4)")
    ap.add_argument("--jobs", "-j", type=int, default=1,
                    help="scenes to render at once (default 1)")
    ap.add_argument("--ssaa", type=int, default=1,
                    help="spatial supersampling: render at ssaa times the resolution")
    ap.add_argument("--shape", default="tent", choices=("tent", "box"),
                    help="time filter over the supersampled frames (default tent)")
    ap.add_argument("--codec", default="ffv1", choices=sorted(CODECS))
    ap.add_argument("--depth", type=int, default=8, choices=(8, 16),
                    help="bit depth of the finished video")
    ap.add_argument("--linear-light", action="store_true",
                    help="average in linear light rather than gamma-encoded values")
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--name", default="hyperchoreography")
    ap.add_argument("--preview", action="store_true",
                    help="480p, 24 fps, ss 2 -- for checking a scene quickly")
    ap.add_argument("--mp4", action="store_true",
                    help="also write an H.264 copy of the finished film")
    ap.add_argument("--no-narrate", action="store_true",
                    help="skip the synthesised narration track")
    ap.add_argument("--no-concat", action="store_true")
    ap.add_argument("--concat-only", action="store_true",
                    help="skip rendering; just join what is already in --out")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list:
        for mod, cls in SCENES:
            print("  %-16s %s" % (mod, cls))
        return 0

    if a.preview:
        a.res, a.fps, a.ss, a.ssaa, a.codec = "480p", 24.0, 2, 1, "h264"

    wanted = [s for s in SCENES if not a.scenes
              or s[0] in a.scenes or s[1] in a.scenes]
    if not wanted:
        print("no scene matches %s; --list shows them all" % ", ".join(a.scenes))
        return 1

    W, H = RESOLUTIONS[a.res]
    RW, RH = W * a.ssaa, H * a.ssaa
    suffix = CODECS[a.codec][1]
    os.makedirs(a.out, exist_ok=True)

    print("%d scene(s) at %dx%d %.0f fps; rendering %dx%d at %.0f fps, %s filter over %d"
          "%s"
          % (len(wanted), W, H, a.fps, RW, RH, a.fps * a.ss, a.shape,
             2 * a.ss - 1 if a.shape == "tent" else a.ss,
             ", %d at a time" % a.jobs if a.jobs > 1 else ""))

    parts, t0 = [], time.time()
    todo = []
    for mod, cls in wanted:
        dst = os.path.join(a.out, "%s%s" % (mod, suffix))
        parts.append(dst)
        if a.concat_only:
            if not os.path.exists(dst):
                print("  %-16s missing" % mod)
                parts.pop()
            continue
        cmd = [sys.executable, "-m", "expo.renderone", mod, cls,
               "--out", dst, "--width", str(W), "--height", str(H),
               "--render-width", str(RW), "--render-height", str(RH),
               "--fps", str(a.fps), "--ss", str(a.ss), "--codec", a.codec,
               "--depth", str(a.depth), "--shape", a.shape]
        if a.linear_light:
            cmd.append("--linear-light")
        if a.no_narrate:
            cmd.append("--no-narrate")
        todo.append((mod, dst, cmd))

    if todo:
        rc = run_jobs(todo, max(1, a.jobs), a.out)
        if rc:
            return rc

    if parts:
        print("  script -> %s" % os.path.basename(write_script(parts, a.out, a.name)))
        print("  manifest -> %s" % os.path.basename(
            write_manifest(parts, a.out, a.name, W, H, a.fps, a.codec)))

    if len(parts) > 1 and not a.no_concat:
        full = os.path.join(a.out, a.name + suffix)
        print("  joining %d scenes -> %s" % (len(parts), os.path.basename(full)))
        concat(parts, full)
        if a.mp4:
            mp4 = os.path.join(a.out, a.name + ".mp4")
            print("  transcoding -> %s" % os.path.basename(mp4))
            transcode(full, mp4)

    print("done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
