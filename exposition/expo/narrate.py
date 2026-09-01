"""The narration track.

The captions are written to be read aloud, so the deck is paced by how long a
reading actually takes rather than by a guess. Two synthesisers are wired in:

* macOS `say` (the default) — instant, dependency-free, robotic; the original
  stand-in.
* Kokoro, a small neural TTS that runs locally — set EXPO_VOICE to
  `kokoro:<voice>` (e.g. `kokoro:bm_george`), with an optional speed suffix
  (`kokoro:bm_george@1.05`). First use downloads the model; en-GB voices are
  bm_daniel, bm_fable, bm_george, bm_lewis.

Either way the clips are ordinary audio files, so a real recording can still be
dropped in later by replacing them. Each line is synthesised once and cached by
the hash of (voice, rate, text), so a re-render costs nothing and the timing
never moves under you — and CHANGING the voice re-times the whole film, which is
the point: the film is paced by the reading it actually carries.
"""

from __future__ import annotations

import hashlib
import json
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

# Kokoro lives in its OWN venv on an older interpreter (torch and spacy have no
# wheels for the film venv's Python yet), so synthesis is a subprocess exactly as
# `say` is. The model loads per invocation (~3 s); every line is cached forever,
# so the cost is paid once per line ever.
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


def _kokoro_synth(text: str, path: str) -> bool:
    voice, speed = _kokoro_spec()
    env = dict(os.environ)
    # This network re-signs TLS with a CA that macOS trusts and Python's bundle does
    # not, and hf-xet (the hub's Rust transfer backend) reads neither REQUESTS_CA_BUNDLE
    # nor the keychain — so the model download runs on the plain backend with a bundle
    # exported from the system keychain (built beside the venv). Both are inert on a
    # network with honest certificates.
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
    os.makedirs(os.path.dirname(INDEX), exist_ok=True)
    with open(INDEX, "w") as f:
        json.dump(_load_index(), f, indent=0, sort_keys=True)


def key(text: str) -> str:
    return hashlib.sha1(("%s|%d|%s" % (VOICE, RATE, text)).encode()).hexdigest()[:16]


def clip(text: str) -> str | None:
    """The audio file for one line, synthesised if it is not already cached.

    The two backends cache under different extensions, and the sample rate need
    not match the mix: build_track resamples every clip to the 48 kHz bed.
    """
    if not available():
        return None
    k = key(text)
    kokoro = _kokoro_spec() is not None
    path = os.path.join(CACHE, k + (".wav" if kokoro else ".aiff"))
    if not os.path.exists(path):
        os.makedirs(CACHE, exist_ok=True)
        if kokoro:
            if not _kokoro_synth(text, path):
                return None
        else:
            r = subprocess.run(["say", "-v", VOICE, "-r", str(RATE), "-o", path,
                                text], capture_output=True)
            if r.returncode != 0 or not os.path.exists(path):
                return None
    return path


def duration(text: str) -> float | None:
    """How long the line takes to say, in seconds."""
    if not enabled():
        return None
    idx = _load_index()
    k = key(text)
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
    """One WAV holding every line at its own start time.

    `cues` are {"t": seconds, "text": line}. Each clip is delayed to its cue and the
    lot are mixed over a silent bed of the right length, so the track lines up with
    the video without any further trimming.
    """
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
