# The Allura Studio production

The film in `../out/` is assembled, mixed and mastered by
[Allura Studio](../../..//allura_studio) — this directory is the comp, its
inputs, and the proofs. It is the studio's first real production, so the
pipeline itself is the deliverable: assets in, one committed `.allura`
project, lossless out, and a verification that the compositor changed nothing
it was not asked to change.

## Rebuild the film

    make upload        # ../out/hyperchoreography-upload.mp4 — the deliverable
    make verify-pipe   # prove the one-pass route carries the comp's own bytes
    make master        # ../out/master.mkv — FFV1 + FLAC, for the archive
    make verify        # every scene and narration proven bit-exact in that master

**The upload is the master.** `make upload` renders the comp straight into the
encoder: Allura writes FFV1 into a fifo and x265 reads it out the other side, so
the twenty-two gigabytes of intermediate never exist. Lossless HEVC in 10-bit
4:4:4 is a tenth the size of the FFV1 it came from and plays where the FFV1 does
not, and there is nothing in it that the FFV1 held and it does not — `make
verify-pipe` renders ninety frames both ways and compares the encodes frame by
frame, then decodes the result back to RGB and asks for the identity, which is
what the trip through YCbCr gives at 10 bits: every 8-bit RGB triple comes back
exactly.

**Studio range, not full**, and this is the one colour decision worth reading
twice. The upload used to be tagged full range, which cost a code value on 0.4 %
of samples and, much worse, could not be relied on to say so: two encodes from
the same command have differed in the `video_full_range_flag` and nothing else.
Full-range samples under a studio-range flag are expanded by 255/219 by every
conforming decoder, which clips everything below code 16 to black — and the
film's background is 10, 12, 17, so that is the floor of every frame in it. The
symptom is a crushed, banded dark end rather than a wrong-looking picture, which
is why it is worth naming. Studio range says what it is, is bit-exact, and
cannot be read wrong, so it is what the film ships as. (True RGB is the other
way to be bit-exact — x265 will carry the comp's GBR planes untouched and
smaller — but it has to signal `matrix_coefficients = identity`, and a player
that ignores that tag reads G-B-R as Y-Cb-Cr and paints the film green. Studio
range gives up nothing to buy that risk.)

That format is not, however, one a graphics card will decode: 4:4:4 is outside
every consumer hardware decoder's HEVC profile, so a player falls back to
software and a slower machine stutters on it. That is the file's job — it exists
to hand YouTube the picture intact, and YouTube re-encodes it — but it is worth
knowing before concluding the film is at fault. The picture itself is exactly
30 fps throughout, every frame present, which `make upload` checks against the
timeline before it finishes.

`ALLURA` points at the studio binary (default `../../../allura_studio/build/allura`).
Either route renders under `--exact` — every sample on its pixel centre, 1 spp,
no resample at 1:1, linear in and linear out — so the sixteen scenes pass
through the comp **bit-exactly**. `make master` writes that FFV1 out as a file
instead, and `verify_master.py` then checks every scene's frames and every
narration's samples at their placements, which needs a file it can seek in.
(The pipeline's own transparency — that a lossless source leaves an `--exact`
export as its own bytes — is held permanently by `unit_render_passthrough` in
the studio's test suite, so the film needs no separate null render to restate
it.)

## What is in the comp

`hyperchoreography.allura` holds two comps:

* **endcard** (14 s) — a spatial three-body choreography **turning** behind the
  catalogue and repository links and the subscribe line: the film opens on three
  bodies tracing the planar eight and closes on the same three somewhere richer.
  The orbit is footage (`plates/EndOrbit.mkv`), not a held still; everything in
  front of it is a plate.
* **hyperchoreography** — the film: the sixteen scenes at their true encoded
  lengths, each narration at its scene's first frame (sample-exact), the end
  card. Its markers are the YouTube chapters (`chapters.txt`).

**There is no title card.** There was one, for six seconds: the figure eight and
the film's name, both held still. Scene 00 shows the figure eight and the film's
name at 0:27, so the film opened by saying the whole thing twice with a cut
between the two — and the first six seconds of a film about orbits were a
photograph of one. The scene's own opening replaces it, and is better: three
bodies arrive out of nothing already turning, and the curve they share is drawn
through them before anything is claimed. `build_spec.py` puts one back at
`TITLE_SECONDS`, and the plates are still in `plates.py`.

**Transitions are butt joins, decided rather than defaulted.** Every scene
opens and closes on the film's background, so a cut at the boundary *is* a
dip to black — the fade is authored content, not a comp effect, and a
cross-dissolve on top of it would only double-fade. The alternative (true
overlaps) needs the terminal `FadeOut`s removed and those scenes re-rendered,
shifts every timing after the cut, and buys a grammar the film was not paced
for; the dips also give the chapter list clean edges. The end card follows the
same grammar — it fades to the scenes' own background colour, carried as an
exact-valued plate because ffmpeg's `color` source rounds it one code value off
through yuv (`make_plates.sh` says how), which is also why its turning orbit is
rendered opaque on that colour rather than with an alpha channel.

## Regenerating the comp

    make plates        # the text plates, and the end card's turning orbit
    make spec          # probe the scene files, write spec.json + chapters.txt
    make conform       # spec.json -> hyperchoreography.allura, canonically

`prewarm.py` reads the spoken lines out of the scene *sources* and synthesises
them before a render, so a rewritten line does not stall the scene that speaks it
for the length of a model load.

**The narration is Kokoro TTS** (`bm_daniel`, a small neural model running
locally from `narration/.venv-tts` — a separate Python 3.12 venv, because torch
has no wheels for the film venv's 3.14 yet). `make film` re-renders the whole
film against it; `EXPO_VOICE=kokoro:<voice>[@speed]` picks another voice, and
plain `EXPO_VOICE=Daniel` falls back to the original macOS `say` stand-in. On a
network that intercepts TLS, the first model download needs the CA bundle beside
the venv (`expo/narrate.py` wires it); after that everything is offline. A real
human recording still drops in the same way it always did: replace the
`*.narration.wav` files and re-run the pipeline.

`build_spec.py` takes every placement from the FILES, not the manifest, and
prints the drift between them (zero at the moment: the renderer reports the
frames the time filter actually wrote rather than the ones the animation
planned). Each scene's audio still overruns its picture by a frame or two of
silence, which the master lets run under the next scene's black head. The sixteen narration tracks measure within 0.1 LU
of one another; `build_spec.py` lifts them uniformly toward −16 LUFS as far
as the loudest true peak allows (×1.476 → ≈ −21 LUFS at peak −0.13 dBFS for
the current voice — the rest of the distance would need a limiter, which the
gain-only bounce deliberately is not), and `verify_master.py` re-derives the
gain from the spec and checks every sample against the quantised multiply.
Real recordings drop into `../out/` under the same names and everything
above reruns unchanged.
