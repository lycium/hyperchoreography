# exposition

A manim presentation of how the `hyperchoreography` solver works, what its
parameters mean, and what its outputs are. Sixteen scenes, about half an hour,
rendered losslessly with tent-filtered motion blur and a narration track.

```
./setup.sh                                  # .venv, manim, and the orbit data
.venv/bin/python render.py --list           # the running order
.venv/bin/python render.py --preview s05_lbfgs      # one scene, fast
.venv/bin/python render.py --res 1080p --ss 4       # the whole film
.venv/bin/python render.py --res 4k --ss 8 --ssaa 2 --mp4
```

Needs `ffmpeg` on PATH. No LaTeX: the formulas are set in Palatino through Pango,
so nothing has to be installed to render them.

---

## What is on screen is computed, not drawn

Everything numerical in the film comes from a NumPy re-implementation of the
solver's own kernels, in `expo/`:

* `nbody.py` — the action, its gradient and its exact Hessian, mirroring
  `src/action.hpp`, rotating frame included. On every record the film uses, from
  the figure eight to the eleven-dimensional one, it reproduces the stored action
  to nine decimals with a gradient around 1e-12.
* `optim.py` — L-BFGS and the Levenberg–Marquardt Newton, ported line for line
  from `src/optim.hpp`, plus the gauge-aware Morse index from `src/search.hpp`.
  It reproduces the catalogue's stored `morse` and `nullity` exactly on every
  record tested.
* `shoot.py` — a Runge–Kutta integrator and the shooting Newton. On the figure
  eight's stored coefficients it returns a shift residual of 3.7e-11 against the
  record's own stored 3.5e-11.
* `catalog.py` — records read through `hyperchoreography show`, cached in `data/`.
* the existence proofs the film shows (`s09_prove`) are the solver's own `prove --verbose` runs,
  parsed by `make_data.py` into `data/prove_*.json`: the certificate, every validated step's width
  and the rough enclosure it used. The uncertainty ellipse is `shoot.py`'s flow, differentiated.

So when a scene shows a gradient falling from 3.2 to 7.5e-12 in eight steps, that
is a real run of the real algorithm, and the orbit it lands on matches the
catalogue to every digit the catalogue stores.

`selftest.py` pins all of that:

```
.venv/bin/python selftest.py
```

It checks the action and gradient of eight catalogued records against their stored
values, the Morse index and nullity against the stored ones (including the
eleven-dimensional record, at 11 and 6), the Hessian against finite differences,
the three featured descents against the orbits they are supposed to reach, the
`k^(2/3)` law for k-fold covers, the shooting residual against the record's own,
and the time filter against the tent mean it claims to be — that last one on a
ramp of flat frames, where the answer is an integer you can check by hand.

`make_data.py` refreshes `data/` from the catalogue; re-run it after a harvest
renumbers record ids.

## Rendering

Manim's frames never reach a lossy encoder. `expo/pipeline.py` patches the frame
writer to hand raw RGBA straight to an ffmpeg process, which does three jobs in
one pass:

**Time.** Manim renders at `fps * ss` frames per second. Each output frame is a
tent-weighted mean over `2*ss - 1` of them, and consecutive windows overlap by
`ss - 1` frames — so with `--ss 4` the deck renders at 120 fps and every 30 fps
frame is a weighted average of seven of those, sharing three with the frame either
side. That is a reconstruction filter rather than a shutter: motion reads as
continuous instead of stepped. `--shape box` gives the disjoint-block average
instead, for comparison.

**Space.** `--ssaa 2` renders at twice the output resolution and scales back down
with Lanczos.

**Encoding.** FFV1 in Matroska by default, which is lossless; `--codec x264rgb`,
`utvideo`, `prores` and `h264` are also there, and `--depth 16` keeps sixteen bits
per channel through the whole chain. `--mp4` additionally writes a shareable
H.264 copy. Scenes are joined with a stream copy, so the finished film is bit for
bit the scenes that went into it.

`--linear-light` averages in linear light rather than in gamma-encoded values.

## What `out/` contains

One lossless MKV per scene (FFV1 video, FLAC audio) plus, beside each, the same
narration again as a standalone WAV — a compositor wants speech on its own layer
rather than cross-faded along with the picture.

* `manifest.json` — every scene's file, duration, frame count and narration cue
  times, with the frame size, rate and colour policy. Enough to build a timeline
  from without probing anything.
* `narration.md` — the spoken text in order, with the time within each scene at
  which each line appears.
* `<name>.mkv` — the scenes joined by a stream copy, so the film is bit for bit the
  scenes that went into it; `--mp4` adds a shareable H.264 copy.

Every scene opens and closes on black. A cross-dissolve between two of them will
therefore read as a dip to black rather than a true overlap — worth deciding
deliberately rather than discovering.

## Narration

The captions are written to be read aloud, and the film is paced by how long a
reading actually takes rather than by a guess: `expo/narrate.py` synthesises each
line with the system speech synthesiser, measures it, and `reading_time()` returns
that. The clips are then delayed to their cue times, mixed into one track and
muxed alongside the picture without re-encoding a frame.

The audio is a stand-in. To replace it with a real recording, drop the new files
into `narration/cache/` under the same names (the hash of voice, rate and text) or
mute it with `--no-narrate` and lay a recording over the finished video. Set
`EXPO_VOICE` and `EXPO_RATE` to change the synthesised voice.

Deleting `narration/durations.json` re-measures everything, which changes every
timing in the film — so it is worth keeping if you want a re-render to match.

## The running order

| | scene | what it covers |
|---|---|---|
| 00 | `s00_open` | three bodies on one curve |
| 01 | `s01_choreo` | the choreography constraint: N curves become one |
| 02 | `s02_action` | least action, and the reduction to a single loop |
| 03 | `s03_fourier` | the Fourier basis, the modes that are dropped, Kepler scaling |
| 04 | `s04_landscape` | critical points, the Morse index, why minima are not enough |
| 05 | `s05_lbfgs` | phase one: L-BFGS, on a real descent to the figure eight |
| 06 | `s06_newton` | phase two: the damped Newton step, and quadratic convergence |
| 07 | `s07_hard` | the five-fold N=4 orbit from two epicycles; a saddle; the hip-hop |
| 08 | `s08_certify` | leaving the Fourier world: the shooting Newton |
| 09 | `s09_prove` | existence: interval arithmetic, the validated flow, Krawczyk's test |
| 10 | `s10_sameorbit` | covers, the canonical frame, Procrustes, the rigidity gate |
| 11 | `s11_starts` | why random starts fail, and the N-gon's transverse spectrum |
| 12 | `s12_frame` | the rotating frame, and what it buys above six dimensions |
| 13 | `s13_usage` | the command line, and a record field by field |
| 14 | `s14_structure` | the redundancies, the dimension budget, the calibration twist |
| 15 | `s15_close` | what is in the catalogue and what is still open |

## Layout

```
expo/       theme, typesetting, plotting, the solver kernels, the render pipeline
scenes/     one module per scene, each an ExpoScene with a story() method
data/       orbits exported from the catalogue, plus cached computations
narration/  synthesised speech, cached by the hash of the line
```

`expo/base.py` holds the stage: `say()` puts a line in the narration band and
holds it for a spoken reading; `say_with()` does the same while something happens.
`expo/theme.py` is the one place to change the look.
