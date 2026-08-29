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
Dependencies: a C++20 compiler; MPFR + GMP for `refine` (`make NOMPFR=1` drops them). No other libraries —
the linear algebra (symmetric eigensolver, pivoted LU), the optimisers, the Taylor integrator and the MPFR
wrapper are all in `src/` (~2 500 lines).

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
* Symmetry groups in any dimension via a small DSL, including *irreducible* groups that force `deff = d`.

## 2. The variational formulation

The loop is a real Fourier series with period 2π (fixed by the Kepler scaling `q → λq, t → λ^{(α+2)/2} t`):

    q(t) = Σ_m  c_m cos(mt) + s_m sin(mt),      c_m, s_m ∈ R^d,   1 ≤ m ≤ K,   m ≢ 0 (mod N).

Modes `m ≡ 0 (mod N)` only move the centre of mass (the potential is blind to them and the kinetic term forces
them to zero at a critical point), so they are excluded. Using the Z_N symmetry the action of the whole
configuration reduces to a functional of the single curve (overall factor N dropped):

    A[q] = ½ ∫ |q̇|² dt  +  ½ Σ_{k=1}^{N-1} ∫ |q(t) − q(t + 2πk/N)|^{−α} dt .

By Palais' symmetric criticality, critical points of `A` are choreographic solutions; the Euler–Lagrange
equation is literally `q̈ = Σ_k (q_k − q)/|q_k − q|³`. The kinetic term is analytic in the coefficients
(`½ Σ π m² (|c_m|²+|s_m|²)`); the potential is a uniform-sample trapezoid quadrature with `M` nodes (multiple
of N and of 8), spectrally accurate for analytic (collision-free) loops. Only `⌊(N−1)/2⌋` shifts (plus `N/2`
with weight ½) are needed. Sampled data are structure-of-arrays with a doubled, wrap-free layout so every inner
loop is a contiguous, gather-free, vectorisable loop over the sample index (`#pragma omp simd`).
Closed-form: value, gradient, Hessian–vector products (`O(M·N·d)`), full Hessian (`O(nb²·M·d²)`) and
`∂∇A/∂α` (for continuation). Kernel timings (N=3, d=2, K=16): action+gradient 1.1 µs, Hessian 16 µs,
Hessian eigendecomposition 80 µs, ODE shift residual 120 µs.

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
  shift residual `|Φ_{T/N}(Z) − SZ|` must be ≤ 1e-5 (one mode doubling is allowed if it is > 1e-4).
* **Shooting Newton.** The shooting map is then solved to `1e-12` (finite-difference Jacobian, tiny damping
  to suppress the gauge null space, centre of mass removed). Every catalogued orbit is a solution of the true
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
* **Equivalence.** Two loops are the same choreography iff they agree up to time shift, time reversal,
  O(d) (rotations *and* reflections) and body relabelling. Candidates are filtered on invariants (action,
  energy, r.m.s. size; relative 1e-4) and then the Procrustes distance `min_{τ,ε,R} ‖A − R B(ε·+τ)‖` is
  computed: for each discrete shift the optimal rotation is given by the nuclear norm of the d×d cross-
  covariance, and the best shift is refined continuously. Below 1e-3 (relative) the loops are identified.
* **Morse index / nullity** of the full Fourier Hessian (computed only for new records, at ≤ 48 modes) are
  stored; the nullity should be 1 + the number of rotation generators acting non-trivially — a larger value
  flags a degenerate (bifurcating) solution.

## 5. Where the high-dimensional choreographies come from

Random low-mode starts collapse onto the Lagrange circle in d ≥ 3 (58 000 trials in 3-D gave circle + eight
and nothing else). The structured start families (`--starts`) and `continue` are what produce spatial
solutions:

* **torus** (even d): rotations in ⌊d/2⌋ orthogonal planes with distinct modes coprime to N. Their relative
  equilibria — N bodies on a Clifford torus — exist and are found within seconds.
* **vertical** (d ≥ 3): a rotating circle (mode m₁) plus one transverse oscillation (mode m₂) — the
  "unchained polygon" ansatz (Chenciner–Féjoz). Linearising the rotating N-gon gives its transverse mode
  frequencies in closed form, `ω_k² = Σ_{l≠0} (1 − cos 2πkl/N)/d_l³`, `d_l = 2R sin(πl/N)` (k = ±1 is the
  trivial tilt, ω = ω_N), so inertial choreographies of this type sit at rational resonances
  `m₂/m₁ ≈ ω_k/ω_N` with phase pattern `k ≡ m₂·m₁⁻¹ (mod N)`; the draw is biased to those. N=4, k=2 has
  ω₂/ω_N = 1.2156 → the (5,6) solution; N=5, k=2 has 1.3277 → the (7,9) solution. For N=3 all transverse modes
  are tilts (ω = ω_N), so no such family exists — consistent with nothing spatial being found for N=3.
* **kick** (`--seed-from`): a catalogued solution embedded in the current dimension and kicked along its
  softest Hessian directions (or randomly), then Newton. The eight has Morse index 1 in 3-D (0 in 2-D):
  such transverse negative directions are what a kick exploits.
* **irreducible symmetry groups** (`--sym random`, half of the draws in d ≥ 3): a signed d-cycle plus a
  signed transposition/flip with compatible time parts generate a group whose spatial representation is
  irreducible on R^d, so any invariant loop spans all of R^d — every solution in that subspace is guaranteed
  genuinely d-dimensional. The (5,6) solution was independently re-found this way with the group
  `t+2/3 s[2,3,1]; t-0/2 s[1,3,2]`.

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

Figure-eight class: `"t+1/2 s[-1,-2]; t-0 s[1,-2]"`.

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
* `catalog/` — consolidated results of the development runs (`d2-3_n3.bin`, `d2-4_n4.bin`, `d2-6_n5.bin`,
  `d3_n6.bin`: 24 distinct certified choreographies) and `d2n3_eight_60digits.txt`, the figure eight refined
  to 60 digits by `refine`. `list` shows deff, d, N, K, action, energy, Morse, nullity, minsep, residual, hits,
  cover, symmetry/provenance. A record's `d` is the dimension it was found in (its Morse index is taken there:
  the eight has index 0 in 2-D but 1 in 3-D); `deff` is the dimension the loop actually spans.
* `src/` — `action.hpp` (basis, kernels, symmetry), `optim.hpp` (L-BFGS, Newton–LM, inertia), `taylor.hpp`
  (generic Taylor integrator, shooting Newton), `invariants.hpp` (canonical form, Procrustes), `catalog.hpp`,
  `search.hpp` (starts, trial pipeline, certification), `continue.hpp`, `mpreal.hpp`, `linalg.hpp`,
  `main.cpp` (CLI), `tests.cpp`.

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
* **Morse index at truncation.** The index is computed from the certified series truncated to 48 modes; for
  loops with close approaches (hundreds of modes) it may be off by a degenerate direction; `nullity = 0` in
  a record means "unresolved at 48 modes", not a real solution property.
* **Small-amplitude branches near a parent** are identified with the parent by the 1e-3 distance tolerance.
  Lower `--tol-dist` (records are accurate to ~1e-13) when hunting near bifurcations.
* **`--sym random` is wasteful** (about half the draws give an empty fixed subspace or an infinite group).
  A curated list of groups per (d, N) — e.g. the crystallographic point groups and the Ferrario–Terracini
  classes — would be far more productive.
* **Trial cost in d ≥ 4 is dominated by the dense Fourier Hessian** when a mode doubling is needed
  (n = 2·nm·d). A Newton–Krylov variant (Hessian–vector products are cheap) would remove that.
* Close-approach loops (minsep < 0.05) get 500+ modes; the shooting certification is fine, but the Fourier
  representation is then a poor basis — parallel shooting (Simó) would suit them better.
* The `hits` counter is only saved at checkpoints; a kill between checkpoints loses ≤ 30 s of counts.

## 12. Roadmap — what remains to be explored

1. **Rotating-frame ("twisted") formulation.** `q_j(t) = R(2πΩj/N) q(t + 2πj/N)` in a frame rotating at rate
   Ω (integer Ω ⇔ inertial choreographies; non-integer Ω ⇔ relative choreographies). Continuation in Ω turns
   isolated choreographies into families (this is Chenciner–Féjoz's "unchained polygons" mechanism) and in
   d ≥ 4 extends to several rotation planes: the natural home of the torus relative equilibria and their
   Lyapunov families. Cost: a 2×2 rotation per shift in the pair kernel, a per-mode quadratic form for the
   Coriolis/centrifugal terms, and a mode-0 term.
2. **Better branch switching for degenerate (rotational) bifurcations** — see §11 — and multi-threaded
   `continue` (branches are independent).
3. **Curated symmetry groups** for each (d, N) and an automatic symmetry *detection* for found solutions
   (test each candidate generator on the certified series), so records carry their full symmetry group.
4. **Linear stability** (Floquet multipliers from the reduced monodromy `(S⁻¹ DΦ_{T/N})`, N-th roots of the
   full multipliers); the shooting Jacobian is already computed.
5. **Newton–Krylov / preconditioned MINRES** for large K and d ≥ 4; FFT-based synthesis for K ≳ 200.
6. **Strong-force homotopy classes.** At α ≥ 2 minimisers exist in every homotopy class in the plane (Moore,
   Montgomery); systematic enumeration by braid type followed by `--alpha-start` continuation.
7. **Mountain-pass / string method** between catalogued minimisers to harvest index-1 saddles systematically.
8. **Rigorous verification** (interval Newton–Kantorovich on the shooting map in MPFR) to turn certified
   orbits into proofs, à la Kapela–Simó.
9. **Visualisation**: an HTML/WebGL viewer of the catalogue (projections for d ≥ 4, animated bodies).
10. Beyond equal masses and choreographies proper: relative choreographies with several curves
    (Ferrario–Terracini core symmetries) reuse almost all of this machinery.

## 13. References

* C. Simó, *New families of solutions in N-body problems*, 3rd European Congress of Mathematics (2000) —
  variational Fourier representation, trapezoid rule, shooting Newton `Φ(2π/N, Z) − S(Z) = 0`.
* G. Minton, *Choreographies* (gminton.org/choreo.html) — Fourier/FFT representation with L-SR1.
* A. Chenciner, R. Montgomery, *A remarkable periodic solution of the three-body problem* (2000).
* A. Chenciner, J. Féjoz, *Unchained polygons and the N-body problem* (2009) — vertical Lyapunov families
  of the rotating N-gon and spatial choreographies at rational rotation numbers.
* D. Ferrario, S. Terracini, *On the existence of collisionless equivariant minimizers* (2004) — symmetry
  groups and Palais' principle for the N-body action.
