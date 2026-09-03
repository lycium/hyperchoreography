"""Rendering: manim's frames go straight into ffmpeg, and nothing lossy happens."""

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
    """Triangular weights over 2*ss - 1 frames: 1, 2, ... ss, ... 2, 1."""
    up = list(range(1, ss + 1))
    return up + up[-2::-1]


def filter_weights(ss: int, shape: str = "tent"):
    return tent_weights(ss) if shape == "tent" else [1] * ss


class TimeFilter:
    """Tent-weighted means over overlapping windows of supersampled frames."""

    def __init__(self, ss: int, shape: str = "tent", linear_light: bool = False):
        self.ss = max(1, int(ss))
        self.w = filter_weights(self.ss, shape)
        self.L = len(self.w)
        self.total = float(sum(self.w))
        self.linear = linear_light
        # 255 * sum(w) is 65 280 at ss = 16; past that the accumulator must grow.
        self.dtype = np.uint16 if 255 * sum(self.w) <= 65535 else np.uint32
        self.n = 0
        self._open = {}
        self._free = []
        self._wide = self._tmp = None
        if linear_light:
            v = np.arange(256, dtype=np.float64) / 255.0
            self._to_linear = np.round((v ** 2.2) * 255.0).astype(np.uint8)

    def _sized(self, n: int):
        if self._wide is None or len(self._wide) != n:
            self._wide = np.empty(n, dtype=self.dtype)
            self._tmp = np.empty(n, dtype=self.dtype)
            self._free = []
            self._open = {}
        return n

    def _buffer(self, n: int):
        return self._free.pop() if self._free else np.empty(n, dtype=self.dtype)

    def push(self, frame: np.ndarray):
        """Add one supersampled frame; return the finished output frames, if any."""
        buf = frame.reshape(-1)
        if self.linear:
            buf = self._to_linear[buf]
        n = self._sized(len(buf))
        np.copyto(self._wide, buf, casting="unsafe")
        out = []
        for j in range(max(0, -(-(self.n - self.L + 1) // self.ss)),
                       self.n // self.ss + 1):
            i = self.n - j * self.ss
            if not 0 <= i < self.L:
                continue
            self._accumulate(j, self.w[i], n)
            if i == self.L - 1:
                out.append(self._finish(j))
        self.n += 1
        return out

    def _accumulate(self, j: int, weight: int, n: int):
        acc = self._open.get(j)
        if acc is None:
            acc = self._open[j] = self._buffer(n)
            if weight == 1:
                np.copyto(acc, self._wide)
                return
            np.multiply(self._wide, self.dtype(weight), out=acc)
            return
        if weight == 1:
            src = self._wide
        else:
            np.multiply(self._wide, self.dtype(weight), out=self._tmp)
            src = self._tmp
        nxt = self._buffer(n)
        np.add(acc, src, out=nxt)
        self._free.append(acc)
        self._open[j] = nxt

    def _finish(self, j: int) -> bytes:
        acc = self._open.pop(j)
        v = acc.astype(np.float32)
        v *= 65535.0 / (255.0 * self.total)
        if self.linear:
            np.power(v * (1.0 / 65535.0), 1.0 / 2.2, out=v)
            v *= 65535.0
        np.add(v, 0.5, out=v)
        data = v.astype(np.uint16).tobytes()
        self._free.append(acc)
        return data


def filter_chain(out_w: int, out_h: int, render_w: int, render_h: int,
                 depth: int = 8) -> str:
    """What is left for ffmpeg once the time filter has run: space, then depth."""
    parts = ["format=gbrp16le"]
    if (render_w, render_h) != (out_w, out_h):
        parts.append("scale=%d:%d:flags=lanczos" % (out_w, out_h))
    parts.append("format=%s" % ("gbrp16le" if depth > 8 else "gbrp"))
    return ",".join(parts)


def open_encoder(out_path, render_w: int, render_h: int, fps_out: float,
                 out_w: int, out_h: int, codec: str = "ffv1", depth: int = 8,
                 quiet: bool = True) -> subprocess.Popen:
    """An ffmpeg process reading the filtered 16-bit frames and writing the scene."""
    args, _suffix, _ll = CODECS[codec]
    cmd = [ffmpeg_binary(), "-y",
           "-f", "rawvideo", "-pix_fmt", "rgba64le",
           "-s", "%dx%d" % (render_w, render_h),
           "-framerate", "%.10g" % fps_out,
           "-i", "-",
           "-an",
           "-vf", filter_chain(out_w, out_h, render_w, render_h, depth),
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


def attach(proc: subprocess.Popen, ss: int = 1, shape: str = "tent",
           linear_light: bool = False):
    """Redirect manim's frame writer through the time filter into `proc`."""
    from manim.scene.scene_file_writer import SceneFileWriter

    state = {"frames": 0, "out": 0}
    tf = TimeFilter(ss, shape, linear_light)

    def write_frame(self, frame_or_renderer, num_frames: int = 1):
        frame = frame_or_renderer
        if not isinstance(frame, np.ndarray):
            frame = frame.get_frame()
        frame = np.ascontiguousarray(frame, dtype=np.uint8)
        for _ in range(num_frames):
            for done in tf.push(frame):
                proc.stdin.write(done)
                state["out"] += 1
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
