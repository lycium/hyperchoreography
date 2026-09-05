# hyperchoreography

A search engine for equal-mass N-body choreographies and **relative choreographies** in higher dimension,
under the attractive homogeneous `1/r^α` potential. In the inertial ansatz the bodies share one curve,
`q_k(t) = q(t + kT/N)`; rotating-frame records require the distinction in §6 below.
Numerical search/refinement and interval existence proofs are separate stages.

**Correctness and research update:** [AUDIT.md](AUDIT.md) records the mathematical derivations, corrected proof bounds,
regressions, performance measurements, remaining limitations, and discovery priorities. Historical proof
flags now require recomputation; the full catalogue has not been re-proved. The integrator is high-order
Taylor, not Bulirsch–Stoer. See [CREDITS.md](CREDITS.md) for prior work, including Li & Liao's 2025 spatial
three-body dataset, and [data/README.md](data/README.md) for the attributed reference/seed workflow.

```
make                  # native build (clang/gcc; -mcpu=native on arm64, -march=native elsewhere)
make test             # self-checks of every kernel (derivatives, symmetry, integrator, dedup, I/O)
make mpinfo           # which GMP this build links, and which mpn assembly it is actually running

./hyperchoreography search   --d 3 --N 4 --K 24 --threads 16 --minutes 600 --out d3n4.bin   # Ctrl-C any time; rerun to resume
./hyperchoreography continue --root circle --d 3 --N 4 --K 32 --covers 9 --depth 2 --out d3n4.bin
./hyperchoreography list     d3n4.bin                 ./hyperchoreography show   d3n4.bin --id 7
./hyperchoreography verify   d3n4.bin [--id 7]        ./hyperchoreography export d3n4.bin --id 7 --out orbit.csv
./hyperchoreography refine   d3n4.bin --id 7 --digits 100 --out d3n4_7.txt
./hyperchoreography prove    d3n4.bin --id 7 --digits 40 --write
./hyperchoreography merge    all.bin d3n4.bin other.bin [--min-rigid r --min-deff k]
./hyperchoreography extras   d3n4.bin                 ./hyperchoreography symmetry d3n4.bin
make gallery                                          # docs/index.html, served at lycium.github.io/hyperchoreography
make gallery-check                                    # packed geometry, metadata, provenance and script syntax
```

Dependencies: a C++20 compiler; MPFR + GMP for `refine` and `prove` (`make NOMPFR=1` drops them); on macOS,
Accelerate for LAPACK `dsyevd` above n = 64 (`make NOACCEL=1` drops it). Everything else — symmetric
eigensolver, pivoted LU, optimisers, Taylor integrator, MPFR wrapper — is in `src/`.
The gallery tools use Python 3's standard library. `gallery-check` also checks JavaScript syntax when Node
is installed; it reports explicitly when that optional check is unavailable. No browser code is executed.

The GMP you already have is probably the wrong one. GMP binds its `mpn` assembly to a CPU, and the table it
uses in 6.3.0 stops at Zen 3, so anything newer falls through to the baseline k8 path; MSYS2 ships that path
outright. Nothing reports it and the flags do not give it away — a package can record `-march=native` and
dispatch to k8 anyway. `make mpinfo` asks the library which path it took; `tools/deps.sh` builds a GMP and an
MPFR that know the machine into `$HOME/.local/opt/gmp-zen`, which the Makefile picks up (`MPPREFIX=` for
elsewhere). On a Zen 5: `prove --digits 40` in 6.4 s rather than 7.7 s, `refine --digits 100` in 21.0 s
rather than 26.6 s, to identical digits.

Naming: **dimension first, then bodies** (`d3n4.bin`); `list` sorts by (deff, N, action).

The paper is `paper/hyperchoreography.tex` (`tectonic hyperchoreography.tex`), with the
existence certificates it cites in `paper/certificates/`. The 35-minute exposition film is
`exposition/` — sixteen manim scenes, narrated and mastered into one lossless upload;
`exposition/comp/README.md` is its build.

---

## 1. What is in the catalogue

The counts and harvest comparisons below describe the original 1,800-record snapshot. The audit adds a
separate attributed Li–Liao reference subset; it does not retroactively reclassify the original records.

`deff` is the dimension the motion actually occupies; `d` is the ambient dimension it was found in. The
budget of §5, `2⌊N/2⌋`, bounds every record with `N ≥ 4` and is reached exactly at `N = 6`, `7` and `10`; it is
a count for one start family, not a theorem, and `N = 3` has three non-planar records (`d2-3_n3:2–4`,
with historical proof flags requiring revalidation).

Modern records store the **numerically refined initial state** as well as the Fourier series rendering it: `h.ret_err`
is the state's residual (median 2.2e-15, worst 6.0e-11 over 1800 records), `extra[6]` the series', which the
`m²` weighting tends to make looser. The stored finite-Fourier Morse counts identify 76 action-minimum
candidates and a median index of 20; these are not Floquet stability classifications. Unstable return maps
can cost precision without ruling out existence (§4).

| `d` | `N` | max `deff` | file | notes |
|---|---|---|---|---|
| 2–4 | 3–6 | 4 | `d2-3_n3`, `d2-4_n4`, `d2-4_n5`, `d3-4_n6` | the classical planar and spatial families |
| 2–4 | 7–11 | 4 | `d2_n7`, `d3-4_n7`, `d3-4_n8`, `d3-4_n9`, `d3_n10`, `d3_n11` | 256 records, all inertial; `d = 3` is spatial up to `N = 11` |
| 5–6 | 6–8 | 6 | `d5-6_n6`, `d5-6_n7`, `d5-6_n8` | 551 records, 128 at `deff = 6`; frames `su:1,2` (`N = 6, 7`) and `su:1,3` (`N = 8`) at `d = 6`, `"1,2"` at `d = 5` |
| 7 | 8, 10 | 7 | `d7_n8`, `d7_n10` | first `deff = 7`, inertial |
| 7 | 10 | 7 | `d7_n10_g2`, `d7_n10_g2_16` | 154 records, 76 at `deff = 7`, frame `g2:1,6` |
| 7 | 11, 12 | 7 | `d7_n11_g2`, `d7_n12_g2` | 336 and 219 records, frame `g2:2,3` |
| 9 | 10 | 9 | `d9_n10_g2` | 31 records, 13 at `deff = 9`; numerical isolation candidates |
| 11 | 12 | **11** | `d11_n12_g2` | 122 records, 38 at `deff = 11`; frame `g2:1,2,4,5` (§6.3) |

Many of the highest-dimensional records use a **rotating frame** (§6), although the table also includes
inertial seven-dimensional examples. Historically, with relative equilibria filtered out, the
inertial search finds nothing at `d = 8` in 60 331 trials while a calibrated frame finds six in 8 799.
`d = 6` is where that stops — **eight of the `deff = 6` records are inertial** (2, 1 and 5 at `N = 6, 7, 8`),
as are 13 of the `deff = 5` ones. These are empirical search results, not an upper dimension bound. No
high-`deff` relative choreography found so far deforms back to an inertial one (§11, item 1).

Arbitrary-precision refinement works in **any** frame: the MPFR Newton carries `G = exp(2πΩ/N)` and threads
its Jacobian one column per task. The figure eight reaches **1e-71**, `d11_n12_g2:106` 4.8e-4 → **1.4e-48**.

## 2. The variational formulation

The loop is a real Fourier series with period 2π (fixed by the Kepler scaling `q → λq, t → λ^{(α+2)/2} t`):

    q(t) = Σ_m  c_m cos(mt) + s_m sin(mt),      c_m, s_m ∈ R^d,   1 ≤ m ≤ K,   m ≢ 0 (mod N).

Modes `m ≡ 0 (mod N)` only move the centre of mass, so they are excluded. Using the `Z_N` symmetry the action
of the whole configuration reduces to a functional of the single curve (overall factor N dropped):

    A[q] = ½ ∫ |q̇|² dt  +  ½ Σ_{k=1}^{N-1} ∫ |q(t) − q(t + 2πk/N)|^{−α} dt .

By Palais' symmetric criticality, critical points are choreographic solutions. The kinetic term is analytic
in the coefficients; the potential is a trapezoid quadrature with `M` nodes, spectrally accurate for
collision-free loops, needing only `⌊(N−1)/2⌋` shifts. Sampled data are structure-of-arrays with a doubled,
wrap-free layout, so every inner loop is contiguous and vectorisable. Value, gradient, Hessian–vector
products, the full Hessian (two DFT bins per mode pair) and `∂∇A/∂α` are closed-form.

In a frame rotating with `Ω ∈ so(d)` the ansatz is `q_j(t) = exp(Ωt) q(t + 2πj/N)`. The potential is
**unchanged** (rotations are orthogonal); only the kinetic term becomes `½∫|q̇ + Ωq|²`, a per-mode `2d × 2d`
block `π[[m²I − Ω², −2mΩ], [2mΩ, m²I − Ω²]]` in place of the scalar `π m²`. At `(d, N, K) = (7, 7, 32)` the
Hessian costs 355 µs against 5.6 ms for its eigendecomposition, the dominant cost of a high-`d` trial.

## 3. Finding critical points of any Morse index

Most choreographies are saddles of the action, which plain minimisation cannot reach. Each trial runs

1. **a start** (§5), Kepler-rescaled to the optimal size `λ^{α+2} = αU/(2K)`;
2. **phase 1**: randomised-length L-BFGS on the action (→ minimisers) or on `½|∇A|²` (→ any critical point);
3. **phase 2**: Levenberg–Marquardt Newton on `∇A = 0` with the exact Hessian in its eigenbasis,
   `δ = −Σ_k q_k λ_k/(λ_k²+μ)(q_kᵀ∇A)` — quadratic convergence to critical points of any index, with the
   gauge null space ignored automatically because `λ_k ≈ 0` directions receive no step;
4. optional **α-continuation** from a strong-force exponent back to α = 1.

Symmetry groups (§4) restrict the search to a fixed subspace `x = By`.

## 4. Certification, canonical form, de-duplication

A Fourier critical point is never accepted on its own merits. `verify <file>` with no `--id` checks the
stored numerical states and their Fourier renderings; `refine` handles one state in MPFR. These are
numerical checks, distinct from the interval `prove` command below. The Newton runs in double, but where an orbit's
monodromy stalls it above `--mpfr-gate` (1e-12) the solve is repeated in MPFR and rounded back
(`d7_n11_g2:127`, 4.7e-12 → 1.3e-15); where rounding is amplified back up (`:142`, 9e-37 in MPFR, 3.7e-12
rounded) the double answer is kept.

* **ODE validation** — initial conditions are read off the series and integrated with a 22nd-order Taylor
  method over `T/N`; the shift residual must be ≤ `--ret-reject`. Its default 1e-1 is deliberately loose: the
  shooting Newton is the real gate, and tightening it to 1e-5 costs 2.5–4.5× in unique orbits per second.
* **Shooting Newton** — `Φ_{T/N}(Z) = G S Z` with `G = exp(2πΩ/N)`, solved to 1e-12 with Levenberg damping
  adapted on the residual (fixed damping abandoned 27 % of admissible candidates).
* **Fourier re-extraction** — the loop is re-sampled from the refined orbit over one `T/N` segment (using symmetry,
  N× cheaper, immune to the error an unstable orbit accumulates over a full period).
* **k-fold covers** — inertial covers are detected from the gcd of significant modes and unwound.
  Rotating-frame unwinding is disabled: it also changes `Ω`, and the Fourier gcd alone is insufficient.
* **Canonical frame** — coefficients rotated into the principal axes of `XᵀX`; `deff` is the number of
  non-negligible principal values, so a planar loop found in a 4-D search compares equal to its 2-D twin.
* **Rigidity defect** — a relative equilibrium moves rigidly, so every mutual distance is constant (`⌊N/2⌋`
  separations exhaust the test, and a common rotation cancels, so the answer is independent of `Ω`). These are
  trivial however large their `deff`, and they **dominate** the high-`d` search. The stored defect is
  trimodal: a cluster at `≤ 2e-11`, a second at `1.3e-6 … 4.0e-6`, then a **4300× empty gap** to the first
  genuine orbit at `1.7e-2`, so `--min-rigid` defaults to **1e-4**. `merge --min-rigid r --min-deff k` re-gates
  a file on disk, so a threshold change costs seconds rather than a re-harvest.
* **`deff` in a rotating frame** is not the loop's own rank: `exp(Ωt)` raises it when a fixed point sweeps a
  circle and lowers it when a circularly polarised mode at rate `−w` becomes a linear oscillation. At `Ω = 0`
  the two agree exactly, so no inertial record changes.
* **Equivalence** — candidates are filtered on invariants, then matched up to phase, reversal and `O(d)`.
  The proposed alignment must fit all Fourier coefficients and the frame. Equal action or mode powers do
  not establish equivalence; members of a possible continuous family are retained. This is a numerical
  tolerance test, not a rigorous distinctness theorem. Prior-work citations survive duplicate replacement.
* **Morse index / nullity** of the full Fourier Hessian, with the exact gauge directions (the time shift and
  the rotations commuting with `Ω`) lifted out by `H + σGGᵀ` rather than classified by magnitude. That
  centraliser is the kernel of `ad_Ω` on `so(d)`, **not** the coordinate generators `E_ab` tested one at a
  time: on the `g2` torus no single `E_ab` commutes, so the per-generator test leaves the torus directions in
  the spectrum as numerical zeros. A stored `nullity = 0` is impossible — the time shift alone forces 1.
* **The calibration twist `χ*`** (§7) and **the symmetry group of the loop**, detected rather than assumed
  (`symmetry`): for fixed `(ε, θ)` the best `R` is the Procrustes fit, so the group is the zero set of the
  loop's *self*-distance. The relative squared residual bottoms out near 1e-15, which sets the default
  tolerance at 1e-6 on the r.m.s.; a loop for which every shift works is reported as `S¹`.

A generator of `--sym` is `(ε, θ = 2πp/q, R ∈ O(d))` imposing `q(εt + θ) = R q(t)`, written `t+p/q s[±i,±j,…]`
(signed permutation) or `r(i,j,p/q)` (plane rotation). The figure-eight class is `"t+1/2 s[-1,2]; t-0
s[1,-2]"` — its half-period map is `diag(−1,+1)`, not `−I`; `s[-1,-2]` is the *circle's*. `cyc:p` cycles the
`⌊d/2⌋` coordinate planes, `fano:p` is the `d = 7` Fano 7-cycle; both are cyclic with trivial core, so they
do not cap `deff` (§5) — their value is cost, `cyc:1` at `d = 6, N = 5` fixing `r = 52` of `n = 156`
parameters.

### Existence proofs (`prove`)

Everything above is numerical evidence. `prove` turns a record into a theorem: *in a box of radius `r` of a
slice transversal to the gauge group there is exactly one initial state whose flow satisfies
`Φ_{T/N}(Z) = G S Z` exactly*, hence a relative choreography with rotating-frame period 2π. Inertial
periodicity additionally follows when `G^N = I`; a shared inertial curve is a further condition. On the same run it checks that the bodies
span `R^d` (the interval Gram matrix of the positions over the segment, relative to their centre of mass, is positive definite) and that the orbit
is not a relative equilibrium (a pair distance has disjoint enclosures at two times) — every quantity an
interval computed in MPFR with outward rounding (`interval.hpp`), nothing that matters a floating-point
approximation. The formal statements and proofs are in `paper/` (`tectonic hyperchoreography.tex`).

* **The test** is Krawczyk's: with `Y ≈ DF(Z₀)⁻¹` (any real matrix), `K = Z₀ − Y F(Z₀) + (I − Y [DF](B))(B − Z₀)`
  lands strictly inside the box `B` and `‖I − Y [DF](B)‖ < 1`. `F(Z₀)` needs a validated flow from a point,
  `[DF](B)` one of the variational equation over the whole box.
* **The validated flow** is the Taylor recurrence of `taylor.hpp` run unchanged on intervals — the same code,
  a different scalar — plus its linearisation (`Tangent`, differentiated term by term, checked against central
  differences). Each step is the degree-`p−1` polynomial on the tight box plus a Lagrange remainder
  `h^p c_p(W)` on a rough enclosure `W` from the high-order Picard inclusion
  `Σ_{k<q} [0,h]^k c_k + [0,h]^q c_q(W) ⊂ int W` (contraction `(h/ρ)^q`, so the step is tolerance-driven, not
  `h‖Df‖ < 1`). The step is set by that remainder *as the intervals see it on `W`*, which dependency inflates
  by tens of orders over the tight-box estimate.
* **Wrapping is paid for in precision, not in code.** Boxes are propagated naively, so widths grow like
  `exp(∫‖|Df|‖)` — 5e4 over `T/3` for the eight, 1e11 for a spatial `N = 3` orbit, 1e6 for `N = 10`. With
  `--digits D` the precision is `3.32 D + 256` bits, the box radius `10^{−D/2}` and the tolerance `10^{−D−4}`,
  so a radius that fails the contraction (the width of `[DF]` is linear in it) is retried smaller. No Lohner
  QR, no Taylor models — a shorter, checkable argument at the cost of a few hundred bits.
* **The gauge group.** `DF` is singular along the time shift, the translations `G` fixes and the rotations
  commuting with `G`, so the test runs on the slice orthogonal to those generators with as many equations
  dropped; the dropped components follow from conservation of energy, momentum and angular momentum, which
  the flow and `G S` both respect — `Q(Y) = Q(Y′)` forces two states agreeing outside the dropped components
  equal when `∂Q/∂Y_D` is verifiably nonsingular on their hull (the reported *closure*). What to drop is
  chosen greedily on `∂(P,E,L)/∂Y`; a rotation fixing the orbit (deff < d) is dropped with its law.
* **Rotating frames** are exact: `Ω` is decomposed into planes, rates snapped so that coinciding angles
  coincide exactly, and the theorem stated for that frame in its own axes with `G` a block matrix of interval
  cosines. The commutant of `G` is built structurally from the angle classes.
* **Output**: the box's half-width and slice radius, `|Y F|`, the contraction and closure norms, the two
  qualifiers, and interval enclosures of the energy and the action (carried in the flow with its own
  remainder). With `--write` the radius goes into `extra[7]` and layout 2 marks a recomputed proof run.
  Older flags are labelled `legacy-proof (rerun)`. The radius refers to the refined slice, not the stored
  binary64 state; a portable replay witness is not yet serialized. Historical timings below predate the
  corrections; current checks are in [AUDIT.md](AUDIT.md).

Measured at `--digits 40`: the eight in 1 s, every `N = 3` record in under a minute, `N = 10` in `d = 3` in
132 s, `N = 8` in `d = 6` in 216 s, `d = 5, N = 6` in a rotating frame in 41 s, the `deff = 7` champion
`d7_n10_g2:0` (frame `g2:1,2`, slice dimension 135) in 74 s on nine Taylor steps, and the eleven-dimensional
headline `d11_n12_g2:13` (`N = 12`, frame `g2:1,2,4,5`, slice dimension 257, seven gauge generators) in under
ten minutes: seven steps, contraction 8.5e-12, energy and action to twenty digits. Its frame turns one plane
by exactly a quarter turn, which is where an interval cosine must know a zero crossing from an extremum. Cost
per step is `p² · N² d · 2Nd` interval operations; close approaches cost steps. Only `α = 1` is wired; the
recurrences are generic, the frame code assumes nothing about `d`.

## 5. Where the high-dimensional choreographies come from

Random low-mode starts collapse onto the Lagrange circle in `d ≥ 3` (58 000 trials in 3-D gave circle and
eight, nothing else); the structured starts (`--starts`) and `continue` are what produce spatial solutions.
**torus** rotates in ⌊d/2⌋ orthogonal planes, **vertical** is a rotating circle plus one transverse
oscillation (the unchained-polygon ansatz), **hyper** puts one near-resonant transverse mode per pattern
`k = 2…⌊N/2⌋` circularly polarised in pairs — which produced every `deff ≥ 5` record — and **inplane**
(below), **fano** (`d = 7`) and **kick** (`--seed-from`, a catalogued solution perturbed along Hessian
eigenvectors) complete the set. The present kick selects eigenvalues in ascending algebraic order, not
the smallest magnitudes; gauge-aware transverse soft-mode selection remains an improvement to make.

Linearising the rotating N-gon gives the transverse frequencies in closed form,
`ω_k² = Σ_{l≠0} (1 − cos 2πkl/N)/d_l³`, `d_l = 2R sin(πl/N)`, so inertial choreographies of vertical type sit
at rational resonances `m₂/m₁ ≈ ω_k/ω_N` with pattern `k ≡ m₂·m₁⁻¹ (mod N)`. N=4, k=2 gives 1.2156 → the
(5,6) solution; N=5, k=2 gives 1.3277 → the (7,9) solution. For N=3 every transverse mode is a tilt, so this
family produces nothing spatial there; the non-planar `N = 3` records are not of this type.

The **in-plane** block is the other half of the same linearisation (`ngon_inplane_freq`): a quartic in the
same sums `C_k = Σ_{p≠0} g_p cos(2πkp/N)`, of which `C_0 − C_1 = 1` *is* the radius equation. Its real roots,
which exist for every `k ≥ 2` at `N ≥ 7`, give a second resonance family `m₂/m₁ = 1 ± ν̂_k` with
`k ≡ m₂·m₁⁻¹ − 1 (mod N)`, seeded by `--starts inplane`. At equal budget it gives fewer records than the
`vertical,hyper,torus` mixture (21 against 25 at `d = 3, N = 9`) but roughly half of them ones the mixture
never found — so it belongs **in** the mix, not in place of it.

**Why `deff` often saturates.** Invariant low-dimensional subspaces and solver basins bias the search.
For an inertial problem, an exactly lower-dimensional start remains in that subspace under the symmetric
gradient/Hessian operations. This does not prove transverse attraction, nor does it forbid a higher-dimensional
branch. Cubic transverse contraction requires additional nondegeneracy and undamped-Newton assumptions;
LM damping and a rotating frame change that analysis. The historical unsuccessful kicks are evidence about
these starts and solvers, not a nonexistence result.

*The budget.* A vertical-family loop has `deff = 2 + Σ_j rank[c_j, s_j]`: two directions oscillating at the
same mode span only a plane, so extra *modes*, not extra directions, buy dimension. With `⌊N/2⌋ − 1` usable
transverse patterns,

    deff_max(N) = 2⌊N/2⌋ :   N=3→2, 4→4, 5→4, 6→6, 7→6, 8→8, 9→8, 10→10, 12→12 .

This counts one seed family's resonances, not all solutions: it does not establish a minimum N for a
desired dimension. The general center-of-mass-frame bound is `deff ≤ min(d, 2(N−1))`, from the invariant
span of initial positions and velocities. Spatial `N = 3` choreographies already exceed the vertical-family
budget, including the prior work credited below.

## 6. The rotating frame (`--omega`)

`--omega "w1,w2,…"` puts rates in the successive coordinate planes; `su:w1,…` closes the list with
`w_n = −Σw`; `g2:p,q[,r,…]` is the maximal torus of `g₂` — rates `(p, q, p+q)` on the Fourier planes of the
Fano 7-cycle, plus one free rate per leftover coordinate plane. Validated exactly: at `Ω = −1` the Lagrange
circle reappears as the mode-2 curve with an identical action and `|∇A| = 1.3e-15`.

Calibration preservation is an optional geometric restriction, not a requirement of the isotropic N-body
equations. For a chosen form (§7), `su:` puts its rates on coordinate
planes rather than Fano planes, and its calibration defect reads 3.0 at `d = 7`, 6.0 at `d = 9` and 4.0 at
`d = 11`, against 1e-15 for the `g₂` frame.

The associative 3-form is **not** a `d = 7` object. It lives on `span(e₀…e₆)` in every `d ≥ 7`, and any
rotation of the remaining coordinates fixes it pointwise, so `stab(φ) = g₂ ⊕ so(d−7)` — dimension **14, 15,
20** at `d = 7, 9, 11`, computed and pinned by `make test`. This supplies useful frames, not a theorem
that calibration is needed to reach those dimensions.

With `G = exp(2πΩ/N)`, shooting implies `Φ_(2π)(Z) = G^N Z` and
`Q_j(t) = G^(-j) Q_0(t + 2πj/N)`. Integer rates imply `G^N = I`, but do not remove the factor `G^(-j)`:
closed body paths need not be the same inertial curve. `G = I` is sufficient for the usual choreography.
The prover prints the exact dyadic rates of its snapped frame, which can differ slightly from the input.

`continue --param omega` walks `Ω = sΩ₀` under the same pseudo-arclength corrector as `α`, allowing folds
to be traversed; `s = 0` is numerically checked against the genuine inertial problem. Because the kinetic operator is
exactly quadratic in `Ω`, the frame derivative is one kernel rather than a new solver.

### 6.1 Choosing a frame

The search problem depends on `Ω` only through the **multiset of `|rates|`**, which collapses the parameter
box. Raw record counts can overvalue repeated sampling of a possible family: `g2:3,4` at `d = 7, N = 10`
returns 21 same-action candidates. A continuous family needs continuation or a separate theorem, not
just coincident scalar invariants. Degenerate cells
collapse, the cell with the two largest rates returns little, and the smallest distinct rates win — but the
ranking **does not transfer between `N`** (`{2,3,5}` is seventh of eight at `N = 10` and first at `N = 11`
and `12`), so every new `N` needs its own sweep.

Harvest frames: `"1,2"` at `d = 5`; `su:1,2` at `d = 6, N = 6, 7` and `su:1,3` at `N = 8`; `g2:1,6` at
`d = 7, N = 10`; `g2:2,3` at `N = 11, 12`; `g2:1,3,2` at `d = 9`; `g2:1,2,4,5` at `d = 11`.

### 6.2 The leftover rates buy the last dimension

Setting the leftover rate to zero leaves the `d = 7` frame inside `R^d` with a fixed subspace. At `d = 9` that
reaches `deff = 8` but **never 9**; at `d = 8`, where `so(1) = 0` and a `g₂` frame has no leftover rate at
all, `deff = 8` is still reached — so a direction the frame does not rotate **can** be occupied. What survives
at every budget and rate set is narrower: **no `r = 0` run has ever reached the top dimension.** At `d = 11`:

| frame | `deff` 11 / 10 / 9 @ 1× | @ 4× |
|---|---|---|
| `g2:1,2,4,0`, one rate zeroed | 0 / 1 / 6 | 0 / 4 / 12 |
| `g2:1,2,0,0`, both zeroed | 0 / 0 / 0 | no records at all |

With both rates zeroed the frame is the bare `d = 7` torus in `R¹¹` with a fixed 4-space; it returns nothing
while running at the **highest** trial rate of any cell, dying cheaply on the `deff` pre-filter.

### 6.3 `d = 11`

`--omega g2:p,q,r₁,r₂` with `N = 12`, the cheapest `N` the budget allows for `deff = 11`; `K = 24`,
`--starts hyper`, ten cells scored by distinct `deff = 11` records:

| rates | `deff` 11 | 10 | 9 | best `χ*` | trials/s |
|---|---|---|---|---|---|
| {1,2,3,4,5} (`g2:1,2,4,5`) | **8** | 1 | 4 | 81.3 | 16.8 |
| {1,2,3,5,6} (`g2:1,2,5,6`) | 2 | 2 | 3 | 133.0 | 16.5 |
| {1,2,3,4,7} (`g2:3,4,1,2`) | 2 | 0 | 0 | 18.0 | 14.9 |
| {1,2,3,6,7} (`g2:1,6,2,3`) | 1 | 0 | 3 | 45.3 | 13.9 |
| four more, incl. the degenerate `g2:1,2,3,4` | 0 | 0–3 | 1–6 | — | 12.9–14.6 |

The ranking survives a second seed, so `g2:1,2,4,5` is the frame to harvest in. One hour (57 327 trials,
seeds 1, 2 and 11 merged) gives **122 records, 38 of them `deff = 11`**, with stored `nullity = 6`, matching
the expected gauge dimension of the 5-torus plus the time shift. These finite-resolution counts suggest
isolation but do not prove it. The historical checks found no missing Morse data, no near-collision,
every return error below 1e-9, and nothing inside the rigidity gate. The `χ*` champion is id
13, `A = 22.745582770`, `χ* = 81.3` (`tw_rel = 0.32`), numerically checked at shift residual 1.8e-12,
full-period return 1.6e-12 and energy drift 1.8e-15. The audit additionally re-proved existence, span 11 and
non-rigidity for this record; see [AUDIT.md](AUDIT.md) for the exact-frame qualification and bounds.

## 7. The calibration ladder and `χ*`

Let `G ⊂ SO(d)` be a proper subgroup stabilising a `k`-form `ψ`. Then

    A_k[q] = (1/2π) ∮ q ∧ q̇ ∧ … ∧ q^(k−1) dt ∈ Λ^k(R^d)        the jet moment
    χ*_ψ[q] = max over R ∈ O(d) of ⟨A_k, R·ψ⟩                    the twist

`A_k` is `O(d)`-equivariant, so `χ*` is invariant under the whole equivalence group of a record while still
*seeing* `G`, because the `G`-orbit of `ψ` is a proper subvariety of its `O(d)` orbit. The rungs, with
stabiliser dimensions computed rather than looked up and pinned by `make test`:

| `d` | `k` | `ψ` | `G` | `dim G` |
|-----|-----|-----|-----|---------|
| 3   | 3   | volume form            | SO(3)    | 3  |
| 4   | 2   | `Re(dz₁∧dz₂)`          | SU(2)·U(1)| 4 |
| 6   | 3   | `Re(dz₁∧dz₂∧dz₃)`      | SU(3)    | 8  |
| 7   | 3   | associative `φ`        | **G₂**   | 14 |
| 8   | 4   | Cayley `dx₀∧φ + ⋆φ`    | **Spin(7)**| 21 |
| 10  | 5   | `Re(dz₁∧…∧dz₅)`        | SU(5)    | 24 |

The G₂ rung is special in dimension seven, but it is not the only possible proper reduction in odd
dimension. In the complex Fourier basis the jet moment is exactly

    A_k = i^{k(k−1)/2} · Σ over {n₁<…<n_k} ⊂ Z∖{0} with Σn = 0  of  V(n) · z_{n₁} ∧ … ∧ z_{n_k},
    V(n) = Π_{r<s} (n_s − n_r)   (Vandermonde)

— verified against the quadrature to 1e-15 in `make test`. So only **resonant** mode tuples contribute; the
modes must be **distinct**, with spread-out sets dominating. An exactly nonzero k-form moment forces the
span of the represented loop to have dimension at least k. Floating-point nonzero values are not certified
rank bounds, and in a rotating frame the span of q is not necessarily the effective inertial dimension.

`χ*` is silenced by exactly two things. **No resonance** — a relative equilibrium sits at the rotation rates
of `Λ`, generically non-resonant, so every tuple is empty. **Time reversal at the wrong parity of `k`** — if
`q(−t + θ) = R q(t)` then `R·A_k = (−1)^{k(k−1)/2} A_k`, which at `k = 3` is `−1`, so a reversal whose `R`
fixes `A₃` annihilates the twist however resonant the modes are, while at `k = 4` it costs nothing. The data
agrees sharply: 13 of 13 records at `d = 8` carry twist, against 24 of 98 at `k = 3`.

Two cautions, both learned the hard way. **`χ*` is not invariant under the rotating-frame gauge**:
`exp(Ωt) q(t)` is unchanged when a rate moves by `N` and the modes compensate, so `χ*`, read off `q` alone,
moves with the gauge — by a factor of **21** on one orbit whose two records agreed on every other stored
scalar. Twist is comparable *within* a frame, never across them. And **`χ*` ranks orbits only downstream of
the rigidity gate**: a *resonant* relative equilibrium maximises exactly the pairing `χ*` measures (20.0
against 4e-14 for a non-resonant one), so a rigid record once looked like the project's twist champion at
446.5. That is one sharp trap, not a general twist/rigidity correlation.

Implementation qualification: `calib_max` is a multistart local optimization, not a certified global
maximum. The product-of-L² jet normalization is a heuristic scale, not a general integral Hadamard bound
for `k > 2`. Neither a large score nor an embedded G₂ form certifies full ambient dimension.

## 8. Files, resumability, tooling

* `catalog.bin` — native-ABI binary (little-endian on the tested machines), bit-exact doubles:
  `"HYPCHOR1"` then records (`RecHdr`, 168 bytes, then
  `int32 modes[nm]`, `double coef[nm·2·d]`, `double Lsv[d]`, `double pca[d]`, `char sym[]`). Records carry a
  trailing `double extra[]`: `[0]` `χ*`, `[1]` `χ*` over the jet scale, `[2]` the rung's jet order, `[3]`
  `‖A_k‖` over the jet scale, `[4]` the rigidity defect, `[5]` the layout, `[6]` the coefficients' residual,
  `[7]` the proof-run slice radius (`prove --write`, 0 = none), then `d²` entries of `Ω` and, in layout ≥ 1,
  the refined state. Layout 2 marks the corrected proof revision. An optional `SRC1` suffix retains arXiv
  orbit citations through merge/replacement; `show` exposes them as `sources`. Legacy files load without
  trusting their old proof flags. Accessors must read `extra.size()`, not the header's `nextra`.
* `catalog.bin.state` — seed and next trial index. Trials are deterministic functions of `(seed, trial)`, so a
  completed trials are reproducible. The current checkpoint stores the next issued trial, so a crash can
  skip in-flight trials; the time limit stops new work, not an already-running trial. Checkpoints every 30 s. Use a
  different `--seed` per machine and `merge` the results; the better-resolved record wins and hit counts add.
* `src/` — `action.hpp` (basis, kernels, symmetry), `calib.hpp` (the ladder), `g2.hpp` (Fano frame and torus),
  `optim.hpp`, `taylor.hpp`, `invariants.hpp`, `catalog.hpp`, `search.hpp`, `continue.hpp`, `mpreal.hpp`,
  `interval.hpp` (MPFR intervals), `prove.hpp` (validated flow, Krawczyk), `linalg.hpp`, `main.cpp`,
  `tests.cpp`.
* `tools/gallery.py` — the whole catalogue as one self-contained page, standard library only, organised by
  **effective dimension** (the ambient space is a search setting, not a property of the orbit); `--heroes N`,
  `--hero`/`--no-hero FILE#ID` pick the selection at the top, `--split` writes one page per dimension. A
  `deff = 9` orbit and a circle look identical in any two coordinates you pick, so each tile is drawn in the
  record's **own principal frame**: a 3-D orthographic shadow, one panel per principal plane in three states
  (curve, line, grey rule — so `deff = 2×curves + lines` is countable), a bar meter of `√(λₖ/λ₁)`, and a
  mutual-distance ribbon flat exactly for a relative equilibrium. Measured against plotting the first two
  coordinates: **2 of 316** tiles read as a near-circle, against **80 of 333**. Only body 0 is stored; the
  rest follow from `body_j(f) = Gʲ·body_0((f + jS/N) mod S)`.

## 9. Reference values

Period 2π, unit masses, G = 1, action per body; `E_total = −N·A/(6π)` by the virial theorem.

| deff | N | solution | action | energy | Morse | modes |
|---|---|----------|--------|--------|-------|-------|
| 2 | 3 | Lagrange circle | 6.534776057 | −1.0400419 | 0 | 1 |
| 2 | 3 | figure eight | 8.123975492 | −1.29297085712209404 | 0 | 59 |
| 2 | 4 | square RE | 9.153307580 | −1.9423922 | 0 | 1 |
| 4 | 4 | Clifford-torus RE, modes (2,3) | 16.579572165 | −3.5182940 | 6 | 3 |
| 3 | 4 | hip-hop resonance (5,6) | 26.761760441 | −5.6790220 | 8 | 66 |
| 3 | 5 | hip-hop resonance (7,9) | 42.644345524 | −11.311764 | 12 | 99 |
| 4 | 5 | Clifford-torus RE, modes (6,7) | 37.945568979 | −10.065375 | 18 | 7 |
| 7 | 10 | genuinely 7-dimensional, inertial | 61.724155215 | −32.7456814 | 48 | 248 |
| 10 | 10 | `d = 10`, frame `Σw = 0`, time-reversible | 40.995282456 | | | |
| 11 | 12 | `d = 11`, frame `g2:1,2,4,5`, `χ* = 81.3` | 22.745582770 | −14.4802877 | 11 | 23 |

## 10. Known issues

* **Some `continue` branches stall before α = 1**, mixing two causes: branches that genuinely end in a
  collision (the code should recognise the shrinking minimum separation and say so), and switching problems
  where a degenerate *pair* of eigenvalues crosses — rotational symmetry makes the bifurcating set a circle
  of solutions, and the bordered system needs an explicit phase condition.
* **Two records have a floor MPFR did not lift in the reported runs** — this is not a nonexistence proof:
  `d7_n12_g2:61` at 6.0e-11 and `d7_n10_g2_16:111` at 2.8e-12, the
  only records above 3.7e-12. `:61`'s canonical frame has principal values 0.0984/0.0981, a relative gap of
  9e-7, leaving the stored `Ω = R Ω Rᵀ` wrong at ~1e-10 — the frame the record names is not the one its orbit
  solves. Solving for `Ω` alongside `Z` would settle it.
* **Calibrated starts do not work below `N ≈ 12`, for a structural reason.** A twisted loop needs modes in
  *arithmetic* resonance (`Σn = 0`); a start that converges needs the *dynamical* resonance `m₂/m₁ ≈ ω_k`.
  The two are incompatible at small `N` — the closest `1 + ω_j = ω_k` is off by 0.64 at `N = 8` and 0.37 at
  `N = 10`, but only 0.087 at `N = 12`; forcing the resonance at `d = 7, N = 8` gives 4 records against
  `hyper`'s 12. Twist comes from the nonlinear mode cascade, not the start — a falsifiable prediction for
  `N ≳ 12`.
* Double precision was sufficient for many tested search candidates, not for every possible orbit or proof.
  The reported action agreement was 1–2 ulp against MPFR; unsuccessful rescue runs do not establish that
  higher precision cannot recover other candidates. Instability amplifies roundoff and integration error.
  Newton–Krylov was measured and dropped (0.88×–1.57×, no predictor). The `hits` counter is saved only at
  checkpoints, so a kill between them loses ≤ 30 s of counts.

## 11. Roadmap

1. **Branch-switch off the `Ω` isola.** Of 18 branches followed from the `deff = 7` records of
   `d7_n10_g2.bin`, 5 reach `s = 0` and **every one arrives on a relative equilibrium**; 1 is a **closed
   isola** (the `χ* = 16.2` champion) never leaving `s ∈ [0.98, 1.06]`. The Morse index does cross (6↔7) at
   `s ≈ 0.9829` and `s ≈ 0.9986` — what `--depth` is for, and the one untried escape.
2. **Linear stability (Floquet)** — the largest remaining gap, since a Morse index is not linear stability.
   For a relative orbit the return matrix is `(GS)⁻¹ DΦ_{T/N}`, not merely `S⁻¹ DΦ_{T/N}`.
   In the inertial case the latter's measured symplectic defect is 1.3e-07 against
   9.4e+02 for `DΦ_T` taken directly. The shooting Jacobian is already computed; the missing pieces are a
   reliable general eigensolver, symmetry reduction, and conditioning/error checks.
3. **Possible continuous families.** The `g2:3,4` cluster at `A = 19.827310867` has stored nullity 6
   (4 gauge + 2), suggesting two transverse degeneracies. Establish whether it is a family by continuation
   and validated analysis. The lack of analogous signals at `d = 9` and `d = 11` is not a nonexistence theorem.
4. **Mine the symmetry classes instead of guessing them.** `symmetry` emits exactly the text `--sym` consumes,
   so the loop closes: detect the groups that occur, keep those with a small fixed subspace, re-search inside.
5. **More cells.** The current hyper-start recipe suggests `d = 13`, `--omega g2:p,q,r₁,r₂,r₃`, `N ≥ 14`;
   this is not a necessity theorem. `d = 12, N = 12` gave 0 records in
   8 537 trials at a quarter of the trial rate, which is the warning.
6. Better branch switching for degenerate bifurcations; multi-threaded `continue`; FFT synthesis for
   `K ≳ 200`; strong-force homotopy classes; mountain-pass between catalogued minimisers; relative
   choreographies with several curves. For `prove`: a Lohner-type wrapping control would cut the precision the
   large-`N` orbits need, and `α ≠ 1` only needs the potential's interval power wired through.

## 12. References

* X. Li, S. Liao, *Discovery of 10,059 new three-dimensional periodic orbits of general three-body problem*
  (2025), [arXiv:2508.08568v1](https://arxiv.org/html/2508.08568v1). Their 21 equal-mass spatial
  choreographies are a reference set, not discoveries attributable to this engine. See [CREDITS.md](CREDITS.md).
* C. Simó, *New families of solutions in N-body problems*, 3rd European Congress of Mathematics (2000).
* G. Minton, *Choreographies* (gminton.org/choreo.html).
* A. Chenciner, R. Montgomery, *A remarkable periodic solution of the three-body problem* (2000).
* A. Chenciner, J. Féjoz, *Unchained polygons and the N-body problem* (2009).
* D. Ferrario, S. Terracini, *On the existence of collisionless equivariant minimizers* (2004).
