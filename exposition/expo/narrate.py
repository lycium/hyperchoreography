"""The narration track."""

from __future__ import annotations

import hashlib
import json
import re
import os
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CACHE = os.path.join(ROOT, "narration", "cache")
INDEX = os.path.join(ROOT, "narration", "durations.json")

VOICE = os.environ.get("EXPO_VOICE", "Daniel")
RATE = int(os.environ.get("EXPO_RATE", "168"))

_durations = None

TTS_PY = os.path.join(ROOT, "narration", ".venv-tts", "bin", "python")
_KOKORO_SNIPPET = """\
import sys
import numpy as np
import soundfile as sf
from kokoro import KPipeline
voice, speed, text, path = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
pipe = KPipeline(lang_code="b" if voice.startswith("b") else "a", repo_id="hexgrad/Kokoro-82M")
chunks = [audio for _gs, _ps, audio in pipe(text, voice=voice, speed=speed)]
sf.write(path, np.concatenate(chunks), 24000, subtype="PCM_16")
"""


def _kokoro_spec() -> tuple[str, float] | None:
    """(voice, speed) when VOICE names the neural backend, else None."""
    if not VOICE.startswith("kokoro:"):
        return None
    spec = VOICE[len("kokoro:"):] or "bm_george"
    voice, _, speed = spec.partition("@")
    return voice, float(speed) if speed else 1.0


def _kokoro_synth(text: str, path: str, speed: float | None = None) -> bool:
    voice, voice_speed = _kokoro_spec()
    speed = voice_speed if speed is None else voice_speed * speed
    env = dict(os.environ)
    ca = os.path.join(os.path.dirname(TTS_PY), "..", "ca-bundle.pem")
    if os.path.exists(ca):
        env.setdefault("REQUESTS_CA_BUNDLE", os.path.abspath(ca))
        env.setdefault("SSL_CERT_FILE", os.path.abspath(ca))
    env.setdefault("HF_HUB_DISABLE_XET", "1")
    r = subprocess.run([TTS_PY, "-c", _KOKORO_SNIPPET, voice, "%g" % speed, text,
                        path], capture_output=True, env=env)
    if r.returncode != 0:
        import sys
        sys.stderr.write(r.stderr.decode(errors="replace")[-500:])
    return r.returncode == 0 and os.path.exists(path)


def available() -> bool:
    if _kokoro_spec() is not None:
        return os.path.exists(TTS_PY) and shutil.which("ffmpeg") is not None
    return shutil.which("say") is not None and shutil.which("ffmpeg") is not None


def enabled() -> bool:
    """Narration is on unless the environment turns it off."""
    return os.environ.get("EXPO_NARRATE", "1") not in ("0", "", "no") and available()


def _load_index() -> dict:
    global _durations
    if _durations is None:
        try:
            with open(INDEX) as f:
                _durations = json.load(f)
        except (OSError, ValueError):
            _durations = {}
    return _durations


def _save_index():
    """Merge into the file on disk under a lock and replace it atomically, so several"""
    import fcntl
    os.makedirs(os.path.dirname(INDEX), exist_ok=True)
    with open(INDEX + ".lock", "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            with open(INDEX) as f:
                disk = json.load(f)
        except (OSError, ValueError):
            disk = {}
        mine = _load_index()
        disk.update(mine)
        mine.update(disk)
        tmp = INDEX + ".tmp"
        with open(tmp, "w") as f:
            json.dump(disk, f, indent=0, sort_keys=True)
        os.replace(tmp, INDEX)


SAID_AS = [
    (re.compile(r"\bKrawczyk\b"), "[Krawczyk](/k\u0279\u02c8avt\u0283\u026ak/)"),
    (re.compile(r"\bFourier\b"), "[Fourier](/f\u02c8\u028a\u0279ie\u026a/)"),
    (re.compile(r"\bChenciner\b"), "[Chenciner](/\u0283\u0252nsi\u02c8ne\u026a/)"),
    (re.compile(r"\bLagrangian\b"), "[Lagrangian](/l\u0250\u02c8\u0261\u0279\u0251\u02d0n\u0292i\u0259n/)"),
    (re.compile(r"\bindex\b"), "[index](/\u02c8\u026and\u0259ks/)"),
    (re.compile(r"\barithmetic\b(?!\s+resonance)"),
     "[arithmetic](/\u0250\u0279\u02c8\u026a\u03b8m\u0259t\u026ak/)"),
    # Two words: misaki joins the hyphen into one, and the n then assimilates into
    # the b with no closure between them -- "embody".
    (re.compile(r"\bN-body\b"), "[N-body](/\u02c8\u025b\u02d0n b\u02c8\u0252di/)"),
]


PACE = [
    (re.compile(r"^This film is about a program\b"), 0.94),
]


def pace_of(text: str) -> float:
    for pat, speed in PACE:
        if pat.search(text):
            return speed
    return 1.0


def said_as(text: str) -> str:
    """The line as the synthesiser should hear it; the caption is left alone."""
    for pat, rep in SAID_AS:
        text = pat.sub(rep, text)
    return text


def key(text: str) -> str:
    return hashlib.sha1(("%s|%d|%s" % (VOICE, RATE, text)).encode()).hexdigest()[:16]


def _key(text: str, speed: float) -> str:
    """Keys that carry the pace only when it is not the default, so adding PACE to one"""
    return key(text if speed == 1.0 else "%s|%g" % (text, speed))


def clip(text: str) -> str | None:
    """The audio file for one line, synthesised if it is not already cached."""
    if not available():
        return None
    speed = pace_of(text)
    text = said_as(text)
    k = _key(text, speed)
    kokoro = _kokoro_spec() is not None
    path = os.path.join(CACHE, k + (".wav" if kokoro else ".aiff"))
    if not os.path.exists(path):
        os.makedirs(CACHE, exist_ok=True)
        if kokoro:
            if not _kokoro_synth(text, path, speed):
                return None
        else:
            r = subprocess.run(["say", "-v", VOICE, "-r", str(RATE), "-o", path,
                                text], capture_output=True)
            if r.returncode != 0 or not os.path.exists(path):
                return None
    return path


def measured(text: str) -> bool:
    """Whether the line's length is already known, without synthesising it."""
    return enabled() and _key(text, pace_of(text)) in _load_index()


def duration(text: str) -> float | None:
    """How long the line takes to say, in seconds."""
    if not enabled():
        return None
    idx = _load_index()
    k = _key(text, pace_of(text))
    if k in idx:
        return float(idx[k])
    path = clip(text)
    if not path:
        return None
    out = subprocess.run(["ffprobe", "-v", "error", "-show_entries",
                          "format=duration", "-of", "default=nw=1:nk=1", path],
                         capture_output=True, text=True)
    try:
        d = float(out.stdout.strip())
    except ValueError:
        return None
    idx[k] = d
    _save_index()
    return d


def build_track(cues, total: float, out_path: str) -> bool:
    """One WAV holding every line at its own start time."""
    if not cues or not available():
        return False
    items = [(c["t"], clip(c["text"])) for c in cues]
    items = [(t, p) for t, p in items if p]
    if not items:
        return False

    cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "error",
           "-f", "lavfi", "-t", "%.4f" % max(total, 0.1),
           "-i", "anullsrc=channel_layout=mono:sample_rate=48000"]
    for _t, p in items:
        cmd += ["-i", p]

    parts = []
    labels = ["[0:a]"]
    for i, (t, _p) in enumerate(items):
        lab = "a%d" % i
        parts.append("[%d:a]aresample=48000,aformat=channel_layouts=mono,"
                     "adelay=%d[%s]" % (i + 1, int(round(t * 1000)), lab))
        labels.append("[%s]" % lab)
    parts.append("%samix=inputs=%d:normalize=0:duration=first[out]"
                 % ("".join(labels), len(labels)))

    cmd += ["-filter_complex", ";".join(parts), "-map", "[out]",
            "-t", "%.4f" % max(total, 0.1), "-c:a", "pcm_s16le", out_path]
    return subprocess.run(cmd).returncode == 0


def mux(video: str, audio: str, out: str, lossless: bool = True) -> bool:
    """Put the track alongside the picture without touching a single video frame."""
    acodec = ["-c:a", "flac"] if lossless else ["-c:a", "aac", "-b:a", "192k"]
    cmd = ["ffmpeg", "-y", "-hide_banner", "-loglevel", "error",
           "-i", video, "-i", audio, "-map", "0:v:0", "-map", "1:a:0",
           "-c:v", "copy"] + acodec + [out]
    return subprocess.run(cmd).returncode == 0
