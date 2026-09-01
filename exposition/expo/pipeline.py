"""Rendering: manim's frames go straight into ffmpeg, and nothing lossy happens.

Manim's own writer encodes to H.264 at CRF 23, which is not what you want under a
temporal filter -- the compression artefacts would be averaged in along with the
picture. So the writer is patched to hand every frame to an ffmpeg process we own,
raw, and ffmpeg does all three jobs in one pass:

    supersampled frames  ->  tent-weighted tmix  ->  Lanczos downscale  ->  lossless

Manim renders at SS times the output frame rate; every output frame is then a
tent-weighted mean of 2*SS - 1 of those, and consecutive windows overlap by SS - 1
frames. That is a real reconstruction filter rather than a shutter: fast motion
gets smooth, continuous blur instead of the stepping a box filter over disjoint
blocks leaves behind. Space is supersampled the same way, by rendering large and
scaling down with Lanczos.
"""

from __future__ import annotations

import shutil
import subprocess
import sys

import numpy as np

RESOLUTIONS = {
    "480p":  (854, 480),
    "720p":  (1280, 720),
    "1080p": (1920, 1080),
    "1440p": (2560, 1440),
    "2160p": (3840, 2160),
    "4k":    (3840, 2160),
}

CODECS = {
    # name: (ffmpeg args, container suffix, lossless?)
    "ffv1":     (["-c:v", "ffv1", "-level", "3", "-coder", "1", "-context", "1",
                  "-g", "1", "-slices", "16", "-slicecrc", "1"], ".mkv", True),
    "x264rgb":  (["-c:v", "libx264rgb", "-qp", "0", "-preset", "slow"], ".mkv", True),
    "utvideo":  (["-c:v", "utvideo"], ".mkv", True),
    "prores":   (["-c:v", "prores_ks", "-profile:v", "4444", "-qscale:v", "4"], ".mov", False),
    "h264":     (["-c:v", "libx264", "-crf", "16", "-preset", "slow",
                  "-pix_fmt", "yuv420p"], ".mp4", False),
}


def ffmpeg_binary() -> str:
    exe = shutil.which("ffmpeg")
    if not exe:
        sys.exit("ffmpeg not found on PATH -- brew install ffmpeg")
    return exe


def tent_weights(ss: int):
    """Triangular weights over 2*ss - 1 frames: 1, 2, ... ss, ... 2, 1.

    A box filter over disjoint blocks is what a shutter does, but it leaves the
    seams between blocks visible on slow pans -- each output frame is built from
    a completely different set of input frames from its neighbour. A tent over
    *overlapping* windows shares half its frames with the frame either side, so
    motion reads as continuous rather than stepped, and the weighting means the
    frames nearest an output instant count most.
    """
    up = list(range(1, ss + 1))
    return up + up[-2::-1]


def filter_chain(ss: int, fps_out: float, out_w: int, out_h: int,
                 render_w: int, render_h: int, depth: int = 8,
                 linear_light: bool = False, shape: str = "tent") -> str:
    """Time filter, then a Lanczos downscale.

    tmix averages a *sliding* window of the last L frames, so the filter is applied
    at every input frame and the ones we keep are every ss-th. With L = 2 ss - 1 and
    tent weights, output frame k is the weighted mean of input frames
    [k ss, k ss + 2 ss - 2] -- windows that step by ss and overlap by ss - 1, which
    is the overlapping tent asked for. With a box shape L = ss and the windows abut.
    """
    work = "gbrp16le" if (depth > 8 or ss > 1 or linear_light) else "gbrp"
    parts = ["format=%s" % work]

    if linear_light:                       # average light, not gamma-encoded values
        parts.append("lutrgb=r='pow(val/65535,2.2)*65535'"
                     ":g='pow(val/65535,2.2)*65535'"
                     ":b='pow(val/65535,2.2)*65535'")

    if ss > 1:
        w = tent_weights(ss) if shape == "tent" else [1] * ss
        L = len(w)
        parts.append("tmix=frames=%d:weights='%s'" % (L, " ".join(str(v) for v in w)))
        # Keep one window in ss, starting at the first one that is completely
        # filled. The gte is not decoration: ffmpeg's mod follows fmod, so
        # mod(n - (L-1), ss) is also zero for small n, and without the guard the
        # first frame kept would be a half-filled window.
        parts.append("select='gte(n\\,%d)*not(mod(n-%d\\,%d))'" % (L - 1, L - 1, ss))
        parts.append("setpts=N/%.10g/TB" % fps_out)

    if linear_light:
        parts.append("lutrgb=r='pow(val/65535,1/2.2)*65535'"
                     ":g='pow(val/65535,1/2.2)*65535'"
                     ":b='pow(val/65535,1/2.2)*65535'")

    if (render_w, render_h) != (out_w, out_h):
        parts.append("scale=%d:%d:flags=lanczos" % (out_w, out_h))

    parts.append("format=%s" % ("gbrp16le" if depth > 8 else "gbrp"))
    return ",".join(parts)


def open_encoder(out_path, render_w: int, render_h: int, fps_hi: float,
                 ss: int, fps_out: float, out_w: int, out_h: int,
                 codec: str = "ffv1", depth: int = 8, linear_light: bool = False,
                 shape: str = "tent", quiet: bool = True) -> subprocess.Popen:
    """An ffmpeg process reading raw RGBA on stdin and writing the finished scene."""
    args, _suffix, _ll = CODECS[codec]
    cmd = [ffmpeg_binary(), "-y",
           "-f", "rawvideo", "-pix_fmt", "rgba",
           "-s", "%dx%d" % (render_w, render_h),
           "-framerate", "%.10g" % fps_hi,
           "-i", "-",
           "-an",
           "-vf", filter_chain(ss, fps_out, out_w, out_h, render_w, render_h,
                               depth, linear_light, shape),
           "-r", "%.10g" % fps_out]
    cmd += args
    if codec in ("ffv1", "x264rgb", "utvideo") and depth > 8:
        cmd += ["-pix_fmt", "gbrp16le"]
    elif codec in ("ffv1", "utvideo"):
        cmd += ["-pix_fmt", "gbrp"]
    cmd += [str(out_path)]
    if quiet:
        cmd[1:1] = ["-hide_banner", "-loglevel", "error"]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE)


# ---------------------------------------------------------------------------
# the patch
# ---------------------------------------------------------------------------
def attach(proc: subprocess.Popen):
    """Redirect manim's frame writer into `proc`, and disable its own encoding.

    Everything patched here is a writer method that only exists to make files;
    the renderer itself is untouched, so the frames are exactly the ones manim
    would have encoded.
    """
    from manim.scene.scene_file_writer import SceneFileWriter

    state = {"frames": 0}

    def write_frame(self, frame_or_renderer, num_frames: int = 1):
        frame = frame_or_renderer
        if not isinstance(frame, np.ndarray):
            frame = frame.get_frame()
        buf = np.ascontiguousarray(frame, dtype=np.uint8).tobytes()
        for _ in range(num_frames):
            proc.stdin.write(buf)
        state["frames"] += num_frames

    def noop(self, *a, **kw):
        return None

    SceneFileWriter.write_frame = write_frame
    SceneFileWriter.open_partial_movie_stream = noop
    SceneFileWriter.close_partial_movie_stream = noop
    SceneFileWriter.join_all_encode_jobs = noop
    SceneFileWriter.combine_to_movie = noop
    SceneFileWriter.combine_to_section_videos = noop
    SceneFileWriter.clean_cache = noop
    SceneFileWriter.flush_cache_directory = noop
    SceneFileWriter.print_file_ready_message = noop
    return state


def close(proc: subprocess.Popen):
    try:
        proc.stdin.close()
    except (BrokenPipeError, ValueError):
        pass
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError("ffmpeg exited with %d" % rc)
