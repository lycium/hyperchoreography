# The Allura Studio production

The film in `../out/` is assembled, mixed and mastered by
[Allura Studio](../../..//allura_studio) — this directory is the comp, its
inputs, and the proofs. It is the studio's first real production, so the
pipeline itself is the deliverable: assets in, one committed `.allura`
project, lossless out, and a verification that the compositor changed nothing
it was not asked to change.

## Rebuild the master

    make master        # ../out/master.mkv  — FFV1 + FLAC, lossless
    make mp4           # ../out/hyperchoreography-upload.mp4 — the lossless YouTube copy
    make verify        # every scene and narration proven bit-exact in place

`ALLURA` points at the studio binary (default `../../../allura_studio/build/allura`).
The master renders under `--exact` — every sample on its pixel centre, 1 spp,
no resample at 1:1, linear in and linear out — so the sixteen scenes pass
through the comp **bit-exactly**: `verify_master.py` checks every scene's
frames and every narration's samples at their placements. (The pipeline's own
transparency — that a lossless source leaves an `--exact` export as its own
bytes — is held permanently by `unit_render_passthrough` in the studio's test
suite, so the film needs no separate null render to restate it.)

## What is in the comp

`hyperchoreography.allura` holds four comps:

* **title** (6 s) — the figure eight with its three bodies over the film's own
  background, the film's name in its own face; plates held still, fades and
  timing by the comp.
* **endcard** (14 s) — a spatial three-body choreography behind the catalogue
  and repository links and the subscribe line: the film opens on three bodies
  tracing the planar eight and closes on the same three somewhere richer.
* **hyperchoreography** (2052.4 s) — the master: title, the sixteen scenes at
  their true encoded lengths, each narration at its scene's first frame
  (sample-exact), the end card. Its markers are the YouTube chapters
  (`chapters.txt`).

**Transitions are butt joins, decided rather than defaulted.** Every scene
opens and closes on the film's background, so a cut at the boundary *is* a
dip to black — the fade is authored content, not a comp effect, and a
cross-dissolve on top of it would only double-fade. The alternative (true
overlaps) needs the terminal `FadeOut`s removed and those scenes re-rendered,
shifts every timing after the cut, and buys a grammar the film was not paced
for; the dips also give the chapter list clean edges. The title and end card
follow the same grammar — they fade to the scenes' own background colour,
carried as an exact-valued plate because ffmpeg's `color` source rounds it
one code value off through yuv (`make_plates.sh` says how).

## Regenerating the comp

    make plates        # render the text/orbit plates (needs ../.venv, ~1 min)
    make spec          # probe the scene files, write spec.json + chapters.txt
    make conform       # spec.json -> hyperchoreography.allura, canonically

**The narration is Kokoro TTS** (`bm_daniel`, a small neural model running
locally from `narration/.venv-tts` — a separate Python 3.12 venv, because torch
has no wheels for the film venv's 3.14 yet). `make film` re-renders the whole
film against it; `EXPO_VOICE=kokoro:<voice>[@speed]` picks another voice, and
plain `EXPO_VOICE=Daniel` falls back to the original macOS `say` stand-in. On a
network that intercepts TLS, the first model download needs the CA bundle beside
the venv (`expo/narrate.py` wires it); after that everything is offline. A real
human recording still drops in the same way it always did: replace the
`*.narration.wav` files and re-run the pipeline.

`build_spec.py` takes every placement from the FILES, not the manifest: the
encoded scenes run 1–2 frames shorter than the manifest's planned timeline
(the tent filter's terminal windows), and each scene's audio overruns its
picture by a frame or two of silence, which the master lets run under the
next scene's black head. The sixteen narration tracks measure within 0.1 LU
of one another; `build_spec.py` lifts them uniformly toward −16 LUFS as far
as the loudest true peak allows (×1.476 → ≈ −21 LUFS at peak −0.13 dBFS for
the current voice — the rest of the distance would need a limiter, which the
gain-only bounce deliberately is not), and `verify_master.py` re-derives the
gain from the spec and checks every sample against the quantised multiply.
Real recordings drop into `../out/` under the same names and everything
above reruns unchanged.
