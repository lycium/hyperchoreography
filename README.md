# hyperchoreography

A search engine for **N-body choreographies in arbitrary dimension**: N equal masses chasing each other along
one closed curve, `q_k(t) = q(t + kT/N)`, under the Newtonian (or any homogeneous `1/r^α`) potential.
Built to run for weeks on many cores, to be resumable and reproducible, to never double-count, and to certify
every catalogued orbit against the true equations of motion — in double, and in MPFR to any number of digits.

```
make                  # native-optimised build (clang/gcc, -march=native; 512-bit vectors preferred on x86-64)
make test             # numerical self-checks of every kernel (derivatives, symmetry, integrator, dedup, I/O)

./hyperchoreography search   --d 3 --N 4 --K 24 --threads 16 --minutes 600 --out d3n4.bin     # Ctrl-C any time; rerun to resume
./hyperchoreography continue --root circle --d 3 --N 4 --K 32 --covers 9 --depth 2 --out d3n4.bin   # bifurcation-tree enumeration
./hyperchoreography list     d3n4.bin                      ./hyperchoreography show d3n4.bin --id 7          # table / JSON dump
./hyperchoreography verify   d3n4.bin --id 7               ./hyperchoreography export d3n4.bin --id 7 --out orbit.csv
./hyperchoreography refine   d3n4.bin --id 7 --digits 100 --out d3n4_7.txt                     # MPFR shooting Newton
./hyperchoreography merge    all.bin d3n4.bin d3n4_machine2.bin                                # union of catalogues
./hyperchoreography bench    --d 4 --N 5 --K 24
```
Dependencies: a C++20 compiler; MPFR + GMP for `refine` (`make NOMPFR=1` drops them); on macOS, Accelerate
for LAPACK `dsyevd` above n = 64 (`make NOACCEL=1` drops it). No other libraries — the linear algebra
(symmetric eigensolver, pivoted LU), the optimisers, the Taylor integrator and the MPFR wrapper are all in
`src/` (~2 500 lines).

Naming convention: **dimension first, then bodies** (`d3n4.bin`); `list` sorts by (deff, N, action).

---------------------------------------------------------------------------------------------------------

## 1. Status — what has been achieved

* A complete, tested pipeline: random or structured starts → L-BFGS → exact-Hessian Newton–LM (critical
  points of *any* Morse index) → ODE certification by shooting Newton → Fourier re-extraction → canonical
  form → symmetry-invariant de-duplication → binary catalogue. ~200–650 trials/s/8 threads for N=3–5,
  d=2–6 on an Apple M5 Pro. Resumable per `(seed, trial)`; multi-machine via `merge`.
* Arbitrary-precision refinement: the eight goes from 1e-11 to **1e-71** in three Newton iterations (1.5 s
  each at 295 bits); 400-bit integration costs 7 ms per `T/N` (allocation-free MPFR kernels).
* Genuinely high-dimensional choreographies, found within seconds and all ODE-certified to ~1e-13:
  the N=4 Clifford-torus relative equilibrium in **d=4** (modes (2,3)), the N=5 one (modes (6,7)),
  the N=4 "hip-hop resonance" spatial choreography in **d=3** (modes (5,6)), a second spatial N=4 solution
  with reflection–reversal symmetry, and the N=5 (7,9) spatial choreography. The linear analysis of the
  rotating N-gon predicts these resonances in closed form (§5) — and the α-continuation tree finds them
  systematically from the k-fold covers of the circle (§6).
* **Rotating-frame ("twisted") formulation**: `q_j(t) = exp(Ωt) q(t + 2πj/N)` for any `Ω ∈ so(d)` (§6).
  Only the kinetic term changes; the potential machinery is untouched. Validated exactly — at `Ω = −1`
  the Lagrange circle reappears as the mode-2 curve with an identical action to 10 digits.
* **The first `d = 6` and `d = 8` choreographies.** Both dimensions were empty simply because nobody had run
  them: 30 s per `N` at `d = 6` gives 43 records (up to `deff = 5`), and 20 minutes at `d = 8` gives 13, all
  `deff = 6` and all carrying non-zero Spin(7) Cayley twist. These runs predate the rigidity defect of §4
  and were not kept, so what fraction of them was rigid is unknown.
* **Genuinely seven-dimensional choreographies.** Two `deff = 7` orbits at `d = 7, N = 10`, certified to
  2e-15 and 6e-15, against a previous catalogue maximum of `deff = 4`. The `deff` budget of §5.1 predicted
  exactly which `N` would work, and their principal values come out as `1 + 2 + 2 + 2`, the polarisation
  structure the budget calls for. Eight of the twenty records in `catalog/d7_n{8,9,10}.bin` are genuine
  orbits; the other twelve are *relative equilibria* — rigid configurations, trivial however high their
  `deff`, which the rigidity defect of §4 now separates and rejects.
* **A ten-dimensional relative choreography.** `d = 10, N = 10` in a rotating frame on the calibrated line
  `Σw = 0` (§14): `deff = 10/10`, `A = 40.995282456`, shift residual 7e-13, full-period return 1e-12, energy
  drift 2e-14, all ten inertial principal values ≥ 0.64, rigidity defect 0.44 — and time-reversible,
  `q(−t) = R q(t)` with `R` negating five axes. This is the `deff` budget's ceiling `2⌊N/2⌋`, reached
  exactly. `d = 8, N = 10` gives six more at `deff = 8`. The frame is what opens them: with relative
  equilibria filtered out, the inertial search finds **none** at `d = 8` in 60 331 trials while `--omega
  su:1,2,3` finds six in 8 799, and `Σw = 0` beats an arbitrary rate vector (`1,2,3,4` gives one).
  `d = 12, N = 12` is not yet reachable: 0 records in 8 537 trials, at a quarter of the trial rate.
* **The calibration ladder** (§14): for the proper subgroup `G ⊂ SO(d)` stabilising a `k`-form `ψ`, the
  pairing of `ψ` with the loop's jet moment is an invariant no `O(d)` invariant of the loop can see. It is
  implemented for every dimension — `G₂` at `d = 7`, `Spin(7)` at `d = 8`, `SU(n)` at `d = 2n` — and the
  stored *twist* `χ*` now discriminates records in every dimension the catalogue reaches, not only `d = 7`.
* **Automatic symmetry detection** (§4, §7): every stored loop's group `q(εt + θ) = R q(t)` is recovered from
  the series, including continuous ones. It reproduces the figure eight's `Z₂ × Z₂` exactly, and it found a
  documentation error in this file that had stood since the first commit.
* Symmetry groups in any dimension via a small DSL.

## 2. The variational formulation

The loop is a real Fourier series with period 2π (fixed by the Kepler scaling `q → λq, t → λ^{(α+2)/2} t`):

    q(t) = Σ_m  c_m cos(mt) + s_m sin(mt),      c_m, s_m ∈ R^d,   1 ≤ m ≤ K,   m ≢ 0 (mod N).

Modes `m ≡ 0 (mod N)` only move the centre of mass (the potential is blind to them and the kinetic term forces
them to zero at a critical point), so they are excluded. Using the Z_N symmetry the action of the whole
configuration reduces to a functional of the single curve (overall factor N dropped):

    A[q] = ½ ∫ |q̇|² dt  +  ½ Σ_{k=1}^{N-1} ∫ |q(t) − q(t + 2πk/N)|^{−α} dt .

In a frame rotating with `Ω ∈ so(d)` (§13) the ansatz becomes `q_j(t) = exp(Ωt) q(t + 2πj/N)`; the potential
term is *unchanged* (rotations are orthogonal) and only the kinetic term becomes `½∫|q̇ + Ωq|²`, a per-mode
`2d × 2d` block `π[[m²I − Ω², −2mΩ], [2mΩ, m²I − Ω²]]` in place of the scalar `π m²`.
By Palais' symmetric criticality, critical points of `A` are choreographic solutions; the Euler–Lagrange
equation is literally `q̈ = Σ_k (q_k − q)/|q_k − q|³`. The kinetic term is analytic in the coefficients
(`½ Σ π m² (|c_m|²+|s_m|²)`); the potential is a uniform-sample trapezoid quadrature with `M` nodes (multiple
of N and of 8), spectrally accurate for analytic (collision-free) loops. Only `⌊(N−1)/2⌋` shifts (plus `N/2`
with weight ½) are needed. Sampled data are structure-of-arrays with a doubled, wrap-free layout so every inner
loop is a contiguous, gather-free, vectorisable loop over the sample index (`#pragma omp simd`).
Closed-form: value, gradient, Hessian–vector products (`O(M·N·d)`), full Hessian (§13 — assembled from two
DFT bins per mode pair, 3.3× faster at `d = 2` up to 8.7× at `d = 7, N = 7, K = 32`) and
`∂∇A/∂α` (for continuation). Kernel timings (N=3, d=2, K=16): action+gradient 1.0 µs, Hessian 6.2 µs,
Hessian eigendecomposition 75 µs, ODE shift residual 112 µs; at (d, N, K) = (7, 7, 32) the Hessian is 355 µs
against an eigendecomposition of 5.6 ms, which is now the dominant cost of a high-dimensional trial.

## 3. Finding critical points of any Morse index

Most choreographies are saddles of the action, which plain minimisation cannot reach (Simó's survey notes that
his gradient method "detects passages close to saddle points … this can be used in the future to try to locate
these solutions"; Minton uses L-SR1 to tolerate indefinite Hessians). Each trial here runs

1. **a start** (§5), Kepler-rescaled to the optimal size `λ^{α+2} = αU/(2K)`;
2. **phase 1**: a randomised-length L-BFGS descent on the action (→ minimisers) or on `½|∇A|²`
   (→ any critical point; gradient `H∇A` is one Hessian–vector product) — `--phase1`;
3. **phase 2**: Levenberg–Marquardt Newton on `∇A = 0` with the exact Hessian in its eigenbasis,
   `δ = −Σ_k q_k λ_k/(λ_k²+μ)(q_kᵀ∇A)`: quadratic convergence to critical points of any index; the gauge null
   space (time shift, rotations) is ignored automatically because `λ_k ≈ 0` directions receive no step;
4. optional **α-continuation** from a strong-force exponent (`--alpha-start 2`) back to α = 1.

Symmetry groups (§4) restrict the search to a fixed subspace `x = By` (reduced gradient/Hessian `Bᵀ∇A`, `BᵀHB`).

## 4. Certification, canonical form, de-duplication

A Fourier critical point is never accepted on its own merits:

* **ODE validation.** Initial conditions of all N bodies are read off the series and integrated with a
  22nd-order Taylor method (automatic differentiation of the N-body recurrences, adaptive step) over `T/N`; the
  shift residual `|Φ_{T/N}(Z) − SZ|` must be ≤ `--ret-reject` (one mode doubling is allowed if it is > 1e-4).
  The default 1e-1 is deliberately loose: the shooting Newton below is the real gate, and tightening it to
  1e-5 discards candidates that shooting would have corrected, costing 2.5–4.5x in unique orbits per second.
* **Shooting Newton.** The shooting map — `Φ_{T/N}(Z) = G S Z` with `G = exp(2πΩ/N)`, the identity in the
  inertial case — is solved to `1e-12` (finite-difference Jacobian, centre of mass removed, Levenberg
  damping adapted on the residual and retried against the same Jacobian; fixed damping abandoned 27 % of
  admissible candidates). Every catalogued orbit is a solution of the true
  equations to machine precision regardless of the Fourier truncation.
* **Fourier re-extraction.** The loop is re-sampled from the certified orbit — one `T/N` segment of all N
  bodies assembled by the choreography symmetry (exact, N× cheaper, and immune to the error growth an unstable
  orbit accumulates over a full period) — and its series taken by DFT down to an amplitude floor of `1e-13`.
  Records store exactly the modes the orbit needs (3 for the torus, 66 for the hip-hop, ~60 for the eight).
* **k-fold covers.** A loop traversed k times (gcd(k,N)=1) is the same choreography with period T/k; it is
  detected from the gcd of the significant modes and unwound (mode km → m, amplitude × k^{2/(α+2)}).
* **Canonical frame.** Coefficients are rotated into the principal axes of `XᵀX` with fixed column signs;
  `deff` = number of non-negligible principal values. A planar loop found in a 4-D search is catalogued with
  `deff = 2` and compares equal to its 2-D twin.
* **Rigidity defect.** A relative equilibrium moves rigidly, so every mutual distance is constant. The
  choreography shift sends the pair `(j,k)` to `(j−k,0)`, so `⌊N/2⌋` separations exhaust the test, and a
  common rotation cancels, so the answer is independent of `Ω`. Zero means a rigid rotating configuration —
  the high-dimensional analogue of the N-gon, and a *trivial* choreography however large its `deff`.
  These dominate the high-`d` search: in a rotating frame at `d = 8, N = 10` they were 9 of the first 13
  records, and 31 of the 170 records catalogued before this test existed are rigid (4 of 4 in
  `d7_n9.bin`, 5 of 8 in `d7_n10.bin`). `search` rejects them by default (`--min-rigid`, default 1e-6);
  the separation is clean, `1e-15…1e-13` against `0.1…0.7`.
* **Effective dimension in a rotating frame.** `deff` is the dimension the motion actually occupies, which
  is *not* the loop's own rank once `Ω ≠ 0`: `exp(Ωt)` raises it when a fixed point sweeps a circle, and
  lowers it when a circularly polarised mode at rate `−w` becomes a linear oscillation. Body `k` traces
  `exp(Ωt)q(t + 2πj/N) = exp(−2πΩk/N)·(`body 0's path`)`, so the occupied subspace is the sum of `N` rotated
  copies, and a plane collapses only when `w_p ≡ 0 (mod N)`. At `Ω = 0` this is the loop's rank exactly, so
  no inertial record changes. The cheap pre-filter still gates on the loop's own rank: measured over 49
  records it discards a genuinely high-`deff` candidate 2 % of the time and buys 18 % more trials.
* **Equivalence.** Two loops are the same choreography iff they agree up to time shift, time reversal,
  O(d) (rotations *and* reflections) and body relabelling. Candidates are filtered on invariants (action,
  energy, r.m.s. size; relative 1e-4) and then the Procrustes distance `min_{τ,ε,R} ‖A − R B(ε·+τ)‖` is
  computed: for each discrete shift the optimal rotation is given by the nuclear norm of the d×d cross-
  covariance, and the best shift is refined continuously. Below 1e-3 (relative) the loops are identified.
* **Morse index / nullity** of the full Fourier Hessian, computed only for new records. The series is
  re-polished to a genuine critical point of the truncation first (the index is undefined anywhere else) and
  the exact gauge directions — the time shift and the rotations that commute with `Ω` — are lifted out of
  the spectrum by `H + σGGᵀ` instead of being classified by magnitude. `(morse, nullity)` is then stable
  across `K_index = 24…80`; before this it was not, and `nullity = 0` always meant this bug.
  The nullity is 1 + the number of rotation generators acting non-trivially; a larger value flags a
  degenerate (bifurcating) solution. A virial check `N·A + 6π·E = 0` guards every certified orbit.
  `hyperchoreography extras` recomputes the index, nullity and twist over an existing catalogue; the values
  are unchanged over a 4× larger index basis (`--K-index 48` against `200`).
* **The calibration twist `χ*`** (§14), stored per record together with its scale-free form and the jet
  order `k` of the rung used.
* **The symmetry group of the loop**, detected rather than assumed (`hyperchoreography symmetry`). For a
  fixed `(ε, θ)` the best `R` is the Procrustes fit, so the group is exactly the set of zeros of the loop's
  *self*-distance over `(ε, θ)` — the same function `loop_distance` already minimises. The fit is done inside
  the `deff`-dimensional span, where `R` is determined and commutes with the inertia tensor, so it comes out
  as the signed permutation the `--sym` DSL wants whenever the principal values are simple. The relative
  squared residual bottoms out near `1e-15` (the `‖A‖²+‖B‖²−2‖AᵀB‖_*` cancellation), which sets the default
  tolerance at `1e-6` on the r.m.s.; a loop for which *every* shift works is reported as `S¹`.

## 5. Where the high-dimensional choreographies come from

Random low-mode starts collapse onto the Lagrange circle in d ≥ 3 (58 000 trials in 3-D gave circle + eight
and nothing else). The structured start families (`--starts`) and `continue` are what produce spatial
solutions:

* **fano** (`d = 7`): a resonant triple `m₁ + m₂ = m₃` placed in three distinct σ-eigenplanes of the Fano
  7-cycle, plus a mode `≡ 0 (mod 7)` on the fixed axis (§14).
* **torus** (even d): rotations in ⌊d/2⌋ orthogonal planes with distinct modes coprime to N. Their relative
  equilibria — N bodies on a Clifford torus — exist and are found within seconds.
* **vertical** (d ≥ 3): a rotating circle (mode m₁) plus one transverse oscillation (mode m₂) — the
  "unchained polygon" ansatz (Chenciner–Féjoz). Linearising the rotating N-gon gives its transverse mode
  frequencies in closed form, `ω_k² = Σ_{l≠0} (1 − cos 2πkl/N)/d_l³`, `d_l = 2R sin(πl/N)` (k = ±1 is the
  trivial tilt, ω = ω_N), so inertial choreographies of this type sit at rational resonances
  `m₂/m₁ ≈ ω_k/ω_N` with phase pattern `k ≡ m₂·m₁⁻¹ (mod N)`; the draw is biased to those. N=4, k=2 has
  ω₂/ω_N = 1.2156 → the (5,6) solution; N=5, k=2 has 1.3277 → the (7,9) solution. For N=3 all transverse modes
  are tilts (ω = ω_N), so no such family exists — consistent with nothing spatial being found for N=3.
* **hyper** (`d ≥ 3`, the family that produced every `deff ≥ 5` record): an `m₁`-fold circle plus one
  near-resonant transverse mode per pattern `k = 2…⌊N/2⌋`, drawn without replacement and circularly
  polarised in pairs while transverse directions remain. See §5.1 for why this, and not kicks, is the lever.
* **kick** (`--seed-from`): a catalogued solution embedded in the current dimension and kicked along its
  softest Hessian directions (or randomly), then Newton. The eight has Morse index 1 in 3-D (0 in 2-D):
  such transverse negative directions are what a kick exploits.
* **symmetry groups** (`--sym`): useful for *organising, re-finding and cheapening* a search — the reduced
  dimension `r ≪ n` makes the Newton eigendecomposition ~100× cheaper at `d = 7` — but they cannot force
  `deff = d`. A generator is `(ε, θ, R)`; forcing `deff = d` needs a trivial core `G ∩ ({+1}×{0}×O(d))`,
  which makes `θ` injective on `{ε = +1}` and hence `G` cyclic or dihedral; real irreps of those have
  dimension ≤ 2, so `Fix(G)` always contains loops in a plane — generically the N-gon itself. Measured:
  500 `--sym random` draws at `d = 5, N = 4` and `d = 7, N = 6` produced **not one** loop with `deff > 2`.
  Conversely an *irreducible* group has an **empty** fixed subspace (§14). The (5,6) solution was
  re-found with `t+2/3 s[2,3,1]; t-0/2 s[1,3,2]`, which is what symmetry is good for.

### 5.1 Why `deff` saturates, and how far it can go

The effective dimension is not limited by the method but by two exact facts.

**The reflection stratum is attracting.** Let `σ_r ∈ O(d)` fix `e_1…e_r`. Since `A∘σ_r = A`, we have
`∇A(σ_r x) = σ_r ∇A(x)` and `H(σ_r x) = σ_r H(x) σ_r`, so `S_r = {loops in R^r}` is invariant under L-BFGS
*and* under the Newton–LM step: a `deff = r` critical point can never be left, only entered. At `x ∈ S_r`
the Hessian is exactly block diagonal, `H = H_∥ ⊕ (I_{d−r} ⊗ H_⊥)` (measured bit-exact), and
`∇_⊥A = H_⊥ξ_⊥ + O(‖ξ_⊥‖³)` is odd in `ξ_⊥`, so the Newton step contracts transversally like `‖ξ_⊥‖³`
*whatever the sign of `H_⊥`*. Instrumented: phase-1 L-BFGS keeps ~50 % of `d`-dimensional starts
`d`-dimensional and the Newton step destroys 99 % of those. Kicking along the softest, the softest
*negative*, or random Hessian directions all raise `deff` in **0 of 816** trials. **Starts are everything.**

**The `deff` budget.** A vertical-family loop has `deff = 2 + Σ_j rank[c_j, s_j]`: two directions
oscillating at the *same* mode span only a plane, so extra directions buy nothing and only extra *modes*
do. With `⌊N/2⌋ − 1` usable transverse patterns `k = 2…⌊N/2⌋`,

    deff_max(N) = 2⌊N/2⌋ :   N=3→2, 4→4, 5→4, 6→6, 7→6, 8→8, 9→8, 10→10 .

This reproduces the whole previous catalogue (max `deff` 4 at N = 4 and 5) and predicted the rest before the
runs: `d = 7` with `N = 4` or `5` cannot exceed `deff = 4` however long it runs (measured over 3·10⁴ trials),
while `N ≥ 8` admits `deff = 7` — and `N = 10` delivered it.
The `hyper` start family excites one near-resonant mode per pattern `k`, circularly polarised in pairs
(`2+2+1` on the five transverse directions of `d = 7`); it produced every `deff = 6` record.

## 6. `continue` — bifurcation trees in the exponent α

Pseudo-arclength continuation of a branch of critical points in `(x, α)`: tangent predictor, corrector =
Newton–LM on the bordered system `[H, ∂∇A/∂α; Gᵀ, 0; tᵀ]` (G = analytic gauge generators, so the corrector is
gauge-fixed), folds are traversed. Along a branch, a change of Morse index means a non-gauge eigenvalue crossed
zero; the crossing is located by bisection, the new branch is entered with the crossing eigenvector as the
arclength tangent (which keeps the corrector off the parent), and followed in both directions. Every crossing
of α = 1 is polished, certified and catalogued. Roots: catalogue records and the k-fold covers of the circle
(`--covers`), because the transverse Lyapunov resonances live on the covers (on the 5-fold circle the (5,6)
hip-hop bifurcates at α ≈ 0.725 with Morse index 10→8; its own secondary bifurcation at α ≈ 1.369 leads back
to the circle).

## 7. Symmetry DSL (`--sym`)

A generator `(ε, θ = 2πp/q, R ∈ O(d))` imposes `q(εt + θ) = R q(t)`; it acts linearly and orthogonally on
each mode's `(c_m, s_m) ∈ R^{2d}`; the generated group is closed per mode (rejected if > 8192 elements),
averaged into the projector onto the fixed subspace, whose eigenvectors are the search basis.

    t+p/q | t-p/q   time map (sign = ε)      s[±i,±j,…]  signed permutation      r(i,j,p/q)  rotation by 2πp/q in plane (i,j)

Figure-eight class: `"t+1/2 s[-1,2]; t-0 s[1,-2]"` — the eight's half-period map is `diag(−1,+1)` (its `x`
carries the odd harmonics and its `y` the even ones), not `−I`. This README used to say `s[-1,-2]`, which is
the *circle's* class: 400 trials under it return `A = 6.534776057` every time, while `s[-1,2]` returns the
eight `A = 8.123975492` in 296 of 300. The automatic detection below found the error.

Two named classes expand to generators: **`cyc:p`** cycles the `⌊d/2⌋` coordinate planes, `q(t + 2πp/n) =
σ q(t)` with `σ` of order `n` (`cyc:1` at `d = 6` is `t+1/3 s[3,4,5,6,1,2]`), and **`fano:p`** is the `d = 7`
Fano 7-cycle of §14. Both are cyclic with trivial core, so unlike a group containing a fixed rotation they do
not cap `deff` (§5.1) — their value is cost. At `d = 6, N = 5, K = 16` the fixed subspace of `cyc:1` has
`r = 52` of `n = 156` parameters, so the Hessian eigendecomposition is ~27× cheaper.

## 8. Arbitrary precision (`refine`)

The shooting Newton runs on the Taylor integrator instantiated with `mpreal` (minimal MPFR wrapper). All
series arithmetic goes through in-place primitives (`fma_add`, `mul`, `div_ui`, …) on preallocated storage —
no temporaries, no allocation inside the integrator. The Jacobian uses central differences with step
`2^{−bits/3}`, so each iteration gains ≈ 2/3 of the working precision. Output: initial conditions, energy,
action and the Fourier coefficients (DFT of a dense-output period) to the requested digits.

## 9. Files, resumability, scaling out

* `catalog.bin` — binary, little-endian, bit-exact doubles: `"HYPCHOR1"` then records (`RecHdr`, 152 bytes:
  ids, N, d, deff, K, M, cover, Morse, nullity, hits, action, energy, r.m.s., minsep, |L|, gradient norm, ODE
  residual, seed, trial, …; then `int32 modes[nm]`, `double coef[nm·2·d]` (per mode `c` then `s`),
  `double Lsv[d]`, `double pca[d]`, `char sym[]`). `hyperchoreography show` dumps a record as JSON.
* `catalog.bin.state` — seed and next trial index; trials are deterministic functions of `(seed, trial)`, so
  a run resumed after Ctrl-C/kill continues exactly, whatever the thread scheduling. Checkpoints every 30 s
  (atomic rename). Use a different `--seed` per machine and `merge` the results. When a duplicate is found
  (in `search`, `continue` or `merge`) the better-resolved record is kept (smaller ODE residual, then fewer
  modes) and the hit counts are added.
* `catalog/` — consolidated results of the development runs (`d2-3_n3.bin`, `d2-4_n4.bin`, `d2-4_n5.bin`,
  `d3-4_n6.bin`, and `d7_n8.bin`, `d7_n9.bin`, `d7_n10.bin` — twenty records, nineteen of them `deff ≥ 6`,
  in seven dimensions) and `d2n3_eight_60digits.txt`, the figure eight refined
  to 60 digits by `refine`. `list` shows deff, d, N, K, action, energy, Morse, nullity, minsep, residual,
  twist and its scale-free form, hits, cover, symmetry/provenance (`--sort twist` orders by the last). A record's `d` is the dimension it was found in (its Morse index is taken there:
  the eight has index 0 in 2-D but 1 in 3-D); `deff` is the dimension the loop actually spans.
* Records carry a trailing `double extra[]` (`nextra` in the header, 0 in legacy files, so old catalogues
  load unchanged): `[0]` the calibration twist `χ*` (§14), `[1]` `χ*` over the jet scale, `[2]` the rung's
  jet order `k`, `[3]` `‖A_k‖` over the jet scale, `[4..7]` reserved, then `d²` entries of `Ω` when the record
  came from a rotating frame (rotated into the record's canonical frame). `extras` recomputes all of these,
  and the Morse index, over an existing catalogue.
* `src/` — `action.hpp` (basis, kernels, symmetry), `calib.hpp` (the calibration ladder: `Λ^k`, jet moments,
  the `SO(d)` ascent), `g2.hpp` (the `d = 7` Fano frame and torus), `optim.hpp` (L-BFGS, Newton–LM, inertia),
  `taylor.hpp` (generic Taylor integrator, shooting Newton), `invariants.hpp` (canonical form, Procrustes,
  symmetry detection), `catalog.hpp`, `search.hpp` (starts, trial pipeline, certification), `continue.hpp`,
  `mpreal.hpp`, `linalg.hpp`, `main.cpp` (CLI), `tests.cpp`.

## 10. Reference values (period 2π, unit masses, G = 1, action per body; `E_total = −N·A/(6π)` by the virial theorem)

| deff | N | solution | action | energy | Morse | modes | provenance |
|---|---|----------|--------|--------|-------|-------|---|
| 2 | 3 | Lagrange circle | 6.534776057 | −1.0400419 | 0 | 1 | any |
| 2 | 3 | figure eight | 8.123975492 | −1.29297085712209404 | 0 | 59 | search |
| 2 | 3 | third planar N=3 | 11.152080126 | −1.7749087 | 1 | 65 | search (rare, 0.7 % of trials) |
| 2 | 4 | square RE | 9.153307580 | −1.9423922 | 0 | 1 | any |
| 2 | 4 | planar N=4 | 11.109471463 | −2.3575030 | 0 | 98 | search |
| 2 | 4 | planar N=4 | 15.047956193 / 17.447668440 | −3.193 / −3.703 | 3 / 6 (in 3-D) | ≥ 500 | search |
| 4 | 4 | Clifford-torus RE, modes (2,3) | 16.579572165 | −3.5182940 | 6 | 3 | torus start |
| 3 | 4 | reflection–reversal symmetric | 16.973317577 | −3.6018500 | 4 | ~130 | `--sym random` |
| 3 | 4 | hip-hop resonance (5,6) | 26.761760441 | −5.6790220 | 8 | 66 | vertical start / continue |
| 2 | 5 | pentagon RE | 11.661751074 | −3.0933760 | 0 | 1 | any |
| 2 | 5 | planar N=5 (several) | 13.77 … 21.98 | | 0–2 | 64–512 | search |
| 4 | 5 | Clifford-torus RE, modes (6,7) | 37.945568979 | −10.065375 | 18 | 7 | torus start |
| 6 | 8 | eight, `d = 7` | 43.015354764 … 123.239443967 | −18.26 … −52.30 | 32–202 | 10–203 | `hyper` start |
| 6 | 9 | four, `d = 7` | 82.775484633 … 153.710633495 | −39.52 … −73.39 | 88–260 | 13–68 | `hyper` start |
| 6 | 10 | six, `d = 7` | 47.139791794 … 122.853934835 | −25.01 … −65.18 | 28–142 | 5–142 | `hyper` start |
| **7** | 10 | **genuinely 7-dimensional**, twist 0.639413 | 61.724155215 | −32.7456814 | 48 | 248 | `hyper` start |
| **7** | 10 | **genuinely 7-dimensional**, twist 0.646759 | 72.684195412 | −38.5601633 | 63 | 327 | `hyper` start |
| 3 | 5 | hip-hop resonance (7,9) | 42.644345524 | −11.311764 | 12 | 99 | continue |

## 11. Known issues and loose ends

* **Some `continue` branches stall before reaching α = 1.** The transverse (vertical) families work end to
  end. In-plane bifurcations of the circle covers (e.g. N=3: ×4 at α≈1.437, ×5 at 1.640 and 0.560; N=4 in 3-D:
  ×7 at 1.674 and 0.694) are detected and switched, and the children travel a finite distance (×4: from 1.44
  down to 1.19; ×7: 1.67 → 1.19) before the corrector fails with the step at `hmin` ("stuck"). Two things
  are mixed here: branches that genuinely end in a collision before α = 1 (Simó reports strong-force
  solutions that "seem not to exist for the Newtonian potential"), which the code should recognise by the
  shrinking minimum separation and report as such; and switching/conditioning problems for branches that
  stall immediately (×5 at 0.56, ×7 at 0.69), where a degenerate *pair* of eigenvalues crosses (rotational
  symmetry makes the bifurcating set a circle of solutions) and the bordered system needs an explicit phase
  condition or continuation in the symmetry-reduced subspace.
* ~~**Morse index at truncation.**~~ Fixed (§4), and the catalogues have now been recomputed with
  `hyperchoreography extras`. Many stored values predated the gauge fix — a stored `nullity = 0` is
  impossible, the time shift alone forces 1 — and are now the gauge-consistent ones (`d=2 → 2`,
  `d=3 → 4`, `d=4 → 6` or `7`), unchanged over a 4× larger index basis.
* **Small-amplitude branches near a parent** are identified with the parent by the 1e-3 distance tolerance.
  Lower `--tol-dist` (records are accurate to ~1e-13) when hunting near bifurcations.
* **`--sym random` is wasteful and, in `d ≥ 5`, useless**: 60–75 % of draws are rejected outright and the
  survivors land on the circle (§5). The named classes `cyc:p` and `fano:p` (§7) are the useful ones —
  cyclic, trivial core, large cost saving — but a curated table per `(d, N)` is still not wired in.
* **De-duplication had a false-negative** on many-mode loops: `d3-4_n6.bin` held three copies of one orbit
  (identical mode-power spectra to 4.9e-15) that the Procrustes distance separated by 7e-3…1e-2. Now
  guarded by a mode-power fingerprint, which is invariant under the whole equivalence group; `merge` folds
  the extra copies. Raising `max_modes` from 64 to 256 does *not* help — the failure is in the Procrustes
  alignment on many-mode loops, not the truncation.
* **Trial cost in d ≥ 4 is dominated by the dense Fourier Hessian and its eigendecomposition** (n = 2·nm·d);
  Accelerate's `dsyevd` covers the latter. A matrix-free Newton–Krylov step was measured and dropped: it only
  breaks even near n ≈ 400 and is erratic beyond (0.88×–1.57×) once force evaluation dominates.
* Double precision suffices: the action matches 40-digit MPFR to 1–2 ulp at 40 k quadrature nodes, and
  re-shooting rejected candidates in MPFR recovers nothing that double misses. The ~1e-11 floor on the T/N
  shift residual is Lyapunov amplification, not summation roundoff.
* **Calibrated starts do not work below `N ≈ 12`, for a structural reason.** §14 says a twisted loop needs
  modes in *arithmetic* resonance, `Σn = 0`; a start that converges needs the *dynamical* resonance
  `m₂/m₁ ≈ ω_k` of §5. The two are incompatible at the `N` the search actually runs at — the closest
  `1 + ω_j = ω_k` is off by 0.87 at `N=6`, **0.64 at `N=8`**, 0.37 at `N=10`, and only 0.087 at `N=12` and
  0.021 at `N=14` (transverse-only triples `ω_i + ω_j = ω_k` are worse everywhere), and deep covers make it
  *harder*, since the rounding slack is `N/(2m₁)`. Measured A/B at `d=7, N=8`, same seed, 90 s: forcing the
  resonance gives **4 records against `hyper`'s 12**, and putting it on the leftover axis gives **2 against
  14**, none of them twisted — while the trial *rate* rises, i.e. the starts fail earlier. So twist in the
  found orbits comes from the nonlinear mode cascade, not from the start. The family was written, measured
  and reverted; roadmap item 5 is the experiment that would change the verdict.
* **The catalogue predates the rigidity test and is contaminated by relative equilibria** — 31 of 170
  records, concentrated exactly where it hurts: 4 of 4 in `d7_n9.bin`, 5 of 8 in `d7_n10.bin`, 3 of 8 in
  `d7_n8.bin`, against 6–34 % in `d ≤ 4`. They are not wrong (a rigid configuration whose bodies lie on one
  curve *is* a choreography) but they are the trivial family, and they were being counted as high-`deff`
  finds. `extras` now stores the defect on every record so they can be listed (`--sort rigid`) or dropped;
  the files themselves are unchanged. The two `deff = 7` headline orbits are not among them.
* Close-approach loops (minsep < 0.05) get 500+ modes; the shooting certification is fine, but the Fourier
  representation is then a poor basis — parallel shooting (Simó) would suit them better.
* The `hits` counter is only saved at checkpoints; a kill between checkpoints loses ≤ 30 s of counts.

## 12. Roadmap — what remains to be explored

1. **Continuation in `Ω`.** The rotating frame is in (§13) but is only ever used at a fixed `Ω`. Continuing
   in it turns isolated choreographies into families — Chenciner–Féjoz's "unchained polygons" mechanism —
   and in `d ≥ 4` along the calibrated line `Σw = 0` of §14, which is where the torus relative equilibria and
   their Lyapunov families live. `continue` already has the pseudo-arclength machinery; it needs `Ω` as a
   second continuation parameter.
2. **Mine the symmetry classes instead of guessing them.** `hyperchoreography symmetry` emits exactly the
   text `--sym` consumes, so the loop closes: detect the groups that actually occur in the catalogue, keep
   the ones with a small fixed subspace, and re-search inside them. That is the curated table §11 wants, and
   it is the one part of the stabiliser programme (§14) still done by hand — the continuous subgroups are in
   (`--omega su:`), the invariant is in (`χ*`), the finite subgroups are not.
3. **Linear stability** (Floquet multipliers from the reduced monodromy `S⁻¹ DΦ_{T/N}`, whose eigenvalues are
   the N-th roots of the full multipliers): the shooting Jacobian is already computed, but the codebase has
   no unsymmetric eigensolver — a Hessenberg + Francis QR (~130 lines) is the missing piece. **The largest
   remaining gap**: the catalogue reports a Morse index, which is not linear stability. Measured on
   `d2-3_n3 #2`: the reduced monodromy has symplectic defect 1.3e-07 against 9.4e+02 for `DΦ_T` taken
   directly, which loses five digits on the largest multiplier and all of them on the small ones.
4. ~~**The `d = 8` Spin(7) frontier.**~~ Opened — see §1. A rotating frame on the calibrated line `Σw = 0`
   reaches `deff = 8` at `d = 8, N = 10` and `deff = 10` at `d = 10, N = 10`, where the inertial search finds
   nothing, and **the frame is essential**: relaxing `Ω = sΩ₀` on five orbits, two branches reach `s = 0` on
   the *rigid* family (`deff` 2 and 4) and three stall, one already rigid to 2e-15 with `deff` still reading
   8. That is natural-parameter continuation, so a stall may be a fold; the pseudo-arclength test wants `Ω`
   in `Continuer` (hard-wired to `α` by `grad_alpha`), cheap since the kinetic operator is exactly quadratic
   in `Ω`: `∂(∇A)/∂s = C + 2sD`. The original argument, still open, was:
   at `k = 4` time reversal does not silence the twist, and every one of the first 13 records carries it,
   against 24 of 98 at `k = 3`. A 20-minute run opened it (§1) and stopped at `deff = 6`; the budget of §5.1
   says `deff = 8` is reachable at `N ≥ 8`. `d = 5` and `d = 6` are equally unexplored — 30 s per `N` at
   `d = 6` produced 43 records where the catalogue had none. The same law says the `d = 7` twist is really a
   *chirality* invariant: at `k = 3` it is non-zero only on orbits that are not reversible with an `R` fixing
   `A₃`, so "maximise twist" there means "find the time-irreversible orbits".
5. **Test the `N ≳ 12` prediction.** §11 shows the calibrated start fails because the N-gon's out-of-plane
   frequencies are never arithmetically resonant at small `N` — but the gap closes as `N` grows (off by 0.64
   at `N = 8`, 0.087 at `N = 12`, 0.021 at `N = 14`). If that is the real obstruction, twisted orbits should
   become common there and a calibrated start should start paying. This is a sharp, falsifiable prediction.
6. **The in-plane half of the N-gon spectrum.** Only the out-of-plane frequencies `ω_k² = C₀ − C_k` are in
   the code (`ngon_vertical_freq`); the in-plane block is a quartic in `ν` whose real roots give a second,
   entirely unexploited family of resonances `m₂/m₁ = 1 ± ν̂_k` with `k ≡ m₂m₁⁻¹ − 1 (mod N)`. ~30 lines for
   `ngon_C(N,k,α)` plus a quartic solve; it generalises `ngon_vertical_freq` to any `α` for free.
7. **Better branch switching for degenerate (rotational) bifurcations** — see §11 — and multi-threaded
   `continue` (branches are independent).
8. **FFT-based synthesis** for K ≳ 200.
9. **Strong-force homotopy classes.** At α ≥ 2 minimisers exist in every homotopy class in the plane (Moore,
   Montgomery); systematic enumeration by braid type followed by `--alpha-start` continuation.
10. **Mountain-pass / string method** between catalogued minimisers to harvest index-1 saddles systematically.
11. **Rigorous verification** (interval Newton–Kantorovich on the shooting map in MPFR) to turn certified
    orbits into proofs, à la Kapela–Simó.
12. **Visualisation**: an HTML/WebGL viewer of the catalogue (projections for d ≥ 4, animated bodies).
13. Beyond equal masses and choreographies proper: relative choreographies with several curves
    (Ferrario–Terracini core symmetries) reuse almost all of this machinery.

## 13. The rotating frame (`--omega`)

    q_j(t) = exp(Ωt) q(t + 2πj/N),   Ω ∈ so(d)   ⇔   q_j(t) = exp(−2πΩj/N) q_0(t + 2πj/N)

— all bodies on one curve *up to a rotation*. Integer rotation rates give inertial choreographies;
non-integer rates give relative periodic orbits, closing only up to `exp(2πΩ)`. The reduction is
remarkably cheap: the potential term is bit-identically unchanged and only the kinetic term moves from
`½∫|q̇|²` to `½∫|q̇ + Ωq|²`, i.e. from the scalar `kin[i] = π m²` to one `2d × 2d` block per mode,

    K_m = π [ m²I − Ω²   −2mΩ ;   2mΩ   m²I − Ω² ] ,

so the gradient gains a per-mode matrix–vector product and the Hessian a per-mode block, both `O(nm·d²)`
— negligible against the potential. Certification changes by one rotation: the shift residual becomes
`Φ_{T/N}(Z) − G S Z`, initial velocities pick up `Ωq`, and the re-extracted loop is un-rotated by
`exp(−Ωt)` before the DFT. Gauge: the time shift survives; of the rotations, only the centraliser of `Ω`.

Syntax: `--omega "w1,w2,…"` sets rates in the successive coordinate planes `(1,2), (3,4), …`;
`--omega su:w1,…` closes the list with `w_n = −Σw` (the calibrated family of §14) and `--omega g2:p,q` is
its `d = 7` octonionic form. `--sym` and `--omega` combine only when each
generator satisfies `RᵀΩR = εΩ` (otherwise the twisted action is not invariant); trials that violate it are
rejected with a named reason rather than silently minimising the wrong functional.

Validated exactly: at `Ω = −1` the Lagrange circle must reappear as the **mode-2** curve with an unchanged
action, and it does — `A = 6.534776057` against the inertial `6.534776057`, `|∇A| = 1.3e-15`, twisted shift
residual `1.7e-15`.

## 14. The calibration ladder

Let `G ⊂ SO(d)` be a proper subgroup stabilising a `k`-form `ψ`. Then

    A_k[q] = (1/2π) ∮ q ∧ q̇ ∧ … ∧ q^(k−1) dt ∈ Λ^k(R^d)        the **jet moment** of the loop
    χ*_ψ[q] = max over R ∈ O(d) of ⟨A_k, R·ψ⟩                    the **twist**

`A_k` is `O(d)`-equivariant, so `χ*` is invariant under the whole equivalence group of a record — `O(d)`,
time shift, time reversal, relabelling — while *seeing* `G`: it is not a function of any `O(d)` invariant,
because the `G`-orbit of `ψ` is a proper subvariety of its `O(d)` orbit. `src/calib.hpp` implements the
whole ladder uniformly (`Λ^k` on increasing multi-indices, exact quadrature, Riemannian ascent on `SO(d)`).

**The rungs, with their stabiliser dimensions computed, not looked up** (null space of `A ↦ A·ψ` on `so(d)`;
these are pinned by `make test`):

| `d` | `k` | `ψ` | `G` | `dim G` | `dim so(d)` |
|-----|-----|-----|-----|---------|-------------|
| 3   | 3   | volume form            | SO(3)    | 3  | 3  |
| 4   | 2   | `Re(dz₁∧dz₂)`          | SU(2)·U(1)| 4 | 6  |
| 6   | 3   | `Re(dz₁∧dz₂∧dz₃)`      | SU(3)    | 8  | 15 |
| 7   | 3   | associative `φ`        | **G₂**   | 14 | 21 |
| 8   | 4   | Cayley `dx₀∧φ + ⋆φ`    | **Spin(7)**| 21 | 28 |
| 10  | 5   | `Re(dz₁∧…∧dz₅)`        | SU(5)    | 24 | 45 |

`R³`'s cross product buys nothing — its 3-form is invariant under all of `SO(3)` — and **7 is the only odd
dimension with any proper reduction**, since the subgroups acting transitively on `S^{n−1}` are `SO`, `U`,
`SU`, `Sp`, `Sp·Sp(1)`, `Sp·U(1)`, `G₂ (n=7)`, `Spin(7) (n=8)`, `Spin(9) (n=16)`. The rung actually used per
dimension is `Spin(7)` at `d = 8`, `G₂` at `d = 7` and `d ≥ 9` (`φ` on a 7-dimensional subspace, stabiliser
`14 + dim so(d−7)` — measured 15 at `d = 9`), `SU(3)` at `d = 6`, and the simple form `e₁∧e₂∧e₃` below that,
whose stabiliser `S(O(3)×O(d−3))` is still proper.

**Why the twist is a resonance invariant.** In the complex Fourier basis `q = Σ_{m≠0} z_m e^{imt}`,
`z_{−m} = z̄_m`, the jet moment is exactly

    A_k = i^{k(k−1)/2} · Σ over sets {n₁<…<n_k} ⊂ Z∖{0} with Σ n = 0  of  V(n) · z_{n₁} ∧ … ∧ z_{n_k},
    V(n) = Π_{r<s} (n_s − n_r)   (Vandermonde)

— verified against the quadrature to 1e-15 relative in `make test`. Three consequences, all of them
observed in the catalogue before they were derived:

* only **resonant** tuples of modes summing to zero contribute, so a non-zero twist requires the loop to
  carry a resonance `m₁ + m₂ = m₃` among its significant modes;
* `V` vanishes on a repeat, so the `k` modes must be **distinct**, and the Vandermonde weight makes
  *spread-out* mode sets dominate;
* if `deff < k` every wedge vanishes, so `χ* ≠ 0` is a **certificate of effective dimension**: measured
  exactly `0` on all 70 planar records in the catalogue, and `0` to 1e-12 for a loop placed in a
  `(k−1)`-plane by construction.

Measured properties of the implementation: `χ*` is unchanged to **3e-15 relative** under a random `O(d)`
frame change for `d = 6…10` (the frame-dependent `⟨A_k, ψ⟩` moves by `O(1)` over the same change), and costs
8 ms per record at `d = 7`, 155 ms at `d = 8`. Getting there needed two things — the ascent is restarted
from signed permutations of the loop's *own* jet frame (every `ψ` in the ladder is `±1` on coordinate
multi-indices, so the maximiser is a signed permutation of that frame plus a local rotation) and the
direction is conjugate-gradient; plain gradient ascent from random frames plateaus at 1e-3 relative in
`d ≥ 9` however many restarts it is given.

**Distinguished rotating frames: the rates must sum to zero.** `exp(Ωt)` multiplies `Θ = (e₀+ie₁)∧…` by
`exp(i Σw_j t)`, so it preserves `Re Θ` **iff `Σ w_j = 0`** — verified as an exact zero of the calibration
defect for `SU(2), SU(3), SU(4), SU(5)`, against `0.30` for a perturbed `Ω`. `--omega su:w₁,…,w_{n−1}`
closes the list with `w_n = −Σw`. The `d = 7` relation `w₃ = w₁ + w₂` of `g2_omega` is the *same* condition
seen through `R⁷ = R ⊕ C³`: the `g₂` maximal torus **is** the `su(3)` maximal torus.

**`d = 7` specifics.** The Fano convention (lines `{i, i+1, i+3} mod 7`) is fixed once, in `calib_psi`; the
7-dimensional cross product used to have its own implementation and axiom test, and both are gone, because
`χ*` depends on `ψ` only through its `O(d)` orbit and `dim stab φ = 14` pins that orbit far more sharply than
the axioms do. `src/g2.hpp` keeps only what the *frame* structure needs: `g2_omega(w₁,w₂)` (the maximal
torus — calibration defect `6e-16` on `w₃ = w₁+w₂`, `O(1)` off it) and `fano_sym(p)` — the twisted class
`q(t + 2πp/7) = σ q(t)` with `σ` the Fano 7-cycle, reachable as `--sym fano:p`. Exactly **1344** of the
645120 signed permutations preserve `φ` (the group `2³·PSL(3,2)`, absolutely irreducible on
`R⁷`), but an irreducible group has an *empty* fixed subspace (§5.1), so the 7-cycle alone is the usable
class: `r = 22…42` of `n = 154…294`. Under it mode `m` lives in the `σ`-eigenplane `k = mp mod 7`, so
reaching `deff = 7` needs a mode `≡ 0 (mod 7)` on the fixed axis — impossible at `N = 7`.

**What it separates.** Over the whole catalogue (`hyperchoreography list --sort twist`):

| `deff` | records | `χ* > 0` | max `χ*` |
|--------|---------|----------|----------|
| 2 | 70 | **0** | 0 |
| 3 | 50 | 9 | **1.217** (`d=3 N=6`, `A = 24.238353047`) |
| 4 | 28 | 8 | 1.216 |
| 5 | 1 | 0 | 0 |
| 6 | 17 | 5 | 0.913 (`d=7 N=8`) |
| 7 | 2  | **2** | 0.647 |

It is a genuine discriminator, not a restatement of `deff`: it cleanly separates the two `d = 7, N = 8`
orbits whose actions agree to four digits (`75.973642` with `χ* = 0`, `75.979900` with `χ* = 0.913`), and
both `deff = 7` orbits carry non-zero twist. The scale-free companion `χ*/Π_r ‖q^(r)‖_{L²}` (column
`tw_rel`) makes records of different size comparable, and its maximum over the catalogue is the same
`d = 3, N = 6` orbit — the most calibrated loop found so far is three-dimensional, not seven.

**What a zero twist means — the causes are now identified, and they are not the obvious one.** Classifying
every record with `deff ≥ k` by whether its significant modes admit a resonant tuple at all:

| | records | continuous `S¹` | mean finite `\|G\|` |
|---|---|---|---|
| `χ* > 0`, resonant | 23 | 0 | 6.7 |
| `χ* = 0`, **no** resonance | 48 | **30** | 21.4 |
| `χ* = 0`, resonant | 26 | 0 | 3.8 |
| `χ* > 0`, non-resonant | 1 | 0 | 2.0 |

(the last row is an artefact of the 1e-3 power cut used to call a mode "significant" — its resonance runs
through a weaker mode). The 48 are explained by the identity itself: a relative equilibrium sits at the
rotation rates of `Λ`, which are generically non-resonant, so every tuple in the sum is empty — and 30 of
them *are* relative equilibria. The 26 are not, and they are not more symmetric either (mean `|G|` 3.8
against 6.7 for the twisted ones), so "heavier symmetry group" is the wrong explanation. Term by term,
**21 of the 26 carry individual terms of size `1e-2` relative to the jet scale that cancel to `1e-16`**
(the other 5 have individually degenerate wedges) — fourteen orders of magnitude, which is forced,
not accidental. The forcing is time reversal. If `q(−t + θ) = R q(t)` then `z̄_n e^{−inθ} = R z_n`, so each
resonant wedge obeys `R·W = e^{−iθΣn} W̄ = W̄`, and therefore

    R · A_k = (−1)^{k(k−1)/2} A_k .

At `k = 3` that sign is **−1**: a reversal whose `R` fixes `A₃` — a palindromic loop being the extreme case —
annihilates the twist outright, however resonant the modes are. At `k = 4` the sign is **+1** and the same
reversal costs nothing. Since essentially every catalogued orbit is reversible, this predicts that the
Spin(7) rung should be far less often silenced than the G₂ one, and the data agrees sharply: **13 of 13**
records at `d = 8` carry twist, against **24 of 98** at `k = 3`. `make test` pins the law directly — a
cosine-only loop on the modes `{1,2,3,4,7}` (which carry both `1+2=3` and `1+2+4=7`) has `‖A₃‖ = 7e-17`
against `0.042` for the same modes with the sines restored, while `‖A₄‖ = 0.0022` is untouched.

So `χ*` is silenced by exactly two things — **no resonance** (the relative equilibria) and **time reversal at
the wrong parity of `k`** — and both are now sharp enough to design against rather than merely observe.

## 15. References

* C. Simó, *New families of solutions in N-body problems*, 3rd European Congress of Mathematics (2000) —
  variational Fourier representation, trapezoid rule, shooting Newton `Φ(2π/N, Z) − S(Z) = 0`.
* G. Minton, *Choreographies* (gminton.org/choreo.html) — Fourier/FFT representation with L-SR1.
* A. Chenciner, R. Montgomery, *A remarkable periodic solution of the three-body problem* (2000).
* A. Chenciner, J. Féjoz, *Unchained polygons and the N-body problem* (2009) — vertical Lyapunov families
  of the rotating N-gon and spatial choreographies at rational rotation numbers.
* D. Ferrario, S. Terracini, *On the existence of collisionless equivariant minimizers* (2004) — symmetry
  groups and Palais' principle for the N-body action.
