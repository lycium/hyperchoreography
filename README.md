# hyperchoreography

A search engine for **N-body choreographies in arbitrary dimension**: N equal masses chasing each other along
one closed curve, `q_k(t) = q(t + kT/N)`, under the Newtonian (or any homogeneous `1/r^α`) potential. Built
to run for weeks on many cores, to be resumable, to never double-count, and to certify every catalogued orbit
against the true equations of motion — in double, and in MPFR to any number of digits.

```
make                  # native build (clang/gcc; -mcpu=native on arm64, -march=native elsewhere)
make test             # self-checks of every kernel (derivatives, symmetry, integrator, dedup, I/O)

./hyperchoreography search   --d 3 --N 4 --K 24 --threads 16 --minutes 600 --out d3n4.bin   # Ctrl-C any time; rerun to resume
./hyperchoreography continue --root circle --d 3 --N 4 --K 32 --covers 9 --depth 2 --out d3n4.bin
./hyperchoreography list     d3n4.bin                 ./hyperchoreography show   d3n4.bin --id 7
./hyperchoreography verify   d3n4.bin [--id 7]        ./hyperchoreography export d3n4.bin --id 7 --out orbit.csv
./hyperchoreography refine   d3n4.bin --id 7 --digits 100 --out d3n4_7.txt
./hyperchoreography prove    d3n4.bin --id 7 --digits 40 --write
./hyperchoreography merge    all.bin d3n4.bin other.bin [--min-rigid r --min-deff k]
./hyperchoreography extras   d3n4.bin                 ./hyperchoreography symmetry d3n4.bin
make gallery                                          # docs/index.html, served at lycium.github.io/hyperchoreography
```

Dependencies: a C++20 compiler; MPFR + GMP for `refine` and `prove` (`make NOMPFR=1` drops them); on macOS,
Accelerate for LAPACK `dsyevd` above n = 64 (`make NOACCEL=1` drops it). Everything else — symmetric
eigensolver, pivoted LU, optimisers, Taylor integrator, MPFR wrapper — is in `src/`.

Naming: **dimension first, then bodies** (`d3n4.bin`); `list` sorts by (deff, N, action).

The paper is `paper/hyperchoreography.tex` (`tectonic hyperchoreography.tex`), with the
existence certificates it cites in `paper/certificates/`. The 35-minute exposition film is
`exposition/` — sixteen manim scenes, narrated and mastered into one lossless upload;
`exposition/comp/README.md` is its build.

---

## 1. What is in the catalogue

`deff` is the dimension the motion actually occupies; `d` is the ambient dimension it was found in. The
budget of §5, `2⌊N/2⌋`, bounds every record with `N ≥ 4` and is reached exactly at `N = 6`, `7` and `10`; it is
a count for one start family, not a theorem, and `N = 3` has three non-planar records (`d2-3_n3:2–4`, proven).

Every record stores the **certified initial state** as well as the Fourier series rendering it: `h.ret_err`
is the state's residual (median 2.2e-15, worst 6.0e-11 over 1800 records), `extra[6]` the series', which the
`m²` weighting always makes looser. Most of these orbits are violently unstable — only 76 are action minima,
median Morse index 20 — which costs precision, not existence (§4).

| `d` | `N` | max `deff` | file | notes |
|---|---|---|---|---|
| 2–4 | 3–6 | 4 | `d2-3_n3`, `d2-4_n4`, `d2-4_n5`, `d3-4_n6` | the classical planar and spatial families |
| 2–4 | 7–11 | 4 | `d2_n7`, `d3-4_n7`, `d3-4_n8`, `d3-4_n9`, `d3_n10`, `d3_n11` | 256 records, all inertial; `d = 3` is spatial up to `N = 11` |
| 5–6 | 6–8 | 6 | `d5-6_n6`, `d5-6_n7`, `d5-6_n8` | 551 records, 128 at `deff = 6`; frames `su:1,2` (`N = 6, 7`) and `su:1,3` (`N = 8`) at `d = 6`, `"1,2"` at `d = 5` |
| 7 | 8, 10 | 7 | `d7_n8`, `d7_n10` | first `deff = 7`, inertial |
| 7 | 10 | 7 | `d7_n10_g2`, `d7_n10_g2_16` | 154 records, 76 at `deff = 7`, frame `g2:1,6` |
| 7 | 11, 12 | 7 | `d7_n11_g2`, `d7_n12_g2` | 336 and 219 records, frame `g2:2,3` |
| 9 | 10 | 9 | `d9_n10_g2` | 31 records, 13 at `deff = 9`, all isolated |
| 11 | 12 | **11** | `d11_n12_g2` | 122 records, 38 at `deff = 11`, all isolated; frame `g2:1,2,4,5` (§6.3) |

Everything above `d = 6` lives in a **rotating frame** (§6): with relative equilibria filtered out, the
inertial search finds nothing at `d = 8` in 60 331 trials while a calibrated frame finds six in 8 799.
`d = 6` is where that stops — **eight of the `deff = 6` records are inertial** (2, 1 and 5 at `N = 6, 7, 8`),
as are 13 of the `deff = 5` ones, so genuine choreographies, not merely relative ones, reach `deff = 6`. No
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

A Fourier critical point is never accepted on its own merits. `recertify` re-derives a whole file's states
and loops; `verify <file>` with no `--id` audits one. The Newton runs in double, but where an orbit's
monodromy stalls it above `--mpfr-gate` (1e-12) the solve is repeated in MPFR and rounded back
(`d7_n11_g2:127`, 4.7e-12 → 1.3e-15); where rounding is amplified back up (`:142`, 9e-37 in MPFR, 3.7e-12
rounded) the double answer is kept.

* **ODE validation** — initial conditions are read off the series and integrated with a 22nd-order Taylor
  method over `T/N`; the shift residual must be ≤ `--ret-reject`. Its default 1e-1 is deliberately loose: the
  shooting Newton is the real gate, and tightening it to 1e-5 costs 2.5–4.5× in unique orbits per second.
* **Shooting Newton** — `Φ_{T/N}(Z) = G S Z` with `G = exp(2πΩ/N)`, solved to 1e-12 with Levenberg damping
  adapted on the residual (fixed damping abandoned 27 % of admissible candidates).
* **Fourier re-extraction** — the loop is re-sampled from the certified orbit over one `T/N` segment (exact,
  N× cheaper, immune to the error an unstable orbit accumulates over a full period).
* **k-fold covers** — a loop traversed k times is the same choreography, detected from the gcd of the
  significant modes and unwound.
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
* **Equivalence** — two loops are the same choreography iff they agree up to time shift, time reversal,
  `O(d)` and relabelling; candidates are filtered on invariants, then compared by Procrustes distance. One
  case defeats that on purpose, a **continuous family**: the action is exactly constant along it, so matching
  to round-off folds its members into one record with a hit count, and `nullity` above the gauge dimension is
  the independent signature.
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
`Φ_{T/N}(Z) = G S Z` exactly*, hence a choreography of period 2π; on the same run it verifies that the bodies
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
  remainder). With `--write` the proven radius goes into `extra[7]`; `list` marks such records `proven`.

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
(below), **fano** (`d = 7`) and **kick** (`--seed-from`, a catalogued solution pushed along its softest
Hessian directions) complete the set.

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

**Why `deff` saturates.** Two exact facts bound it. *The reflection stratum is attracting*: with `σ_r ∈ O(d)`
fixing `e_1…e_r` and `A∘σ_r = A`, the loops in `R^r` are invariant under L-BFGS *and* under the Newton–LM
step, so a `deff = r` critical point can never be left, only entered — the Hessian there is exactly block
diagonal and `∇_⊥A` odd in `ξ_⊥`, so Newton contracts transversally like `‖ξ_⊥‖³` **whatever the sign of
`H_⊥`**. Measured: phase-1 L-BFGS keeps ~50 % of `d`-dimensional starts `d`-dimensional and the Newton step
destroys 99 % of those; kicking along the softest, softest-negative or random directions raises `deff` in
**0 of 816** trials. **Starts are everything**, and symmetry cannot substitute — forcing `deff = d` needs a
trivial core, making the group cyclic or dihedral, whose real irreps have dimension ≤ 2, so `Fix(G)` always
contains planar loops (500 `--sym random` draws at `d = 5, N = 4` and `d = 7, N = 6` gave none with
`deff > 2`).

*The budget.* A vertical-family loop has `deff = 2 + Σ_j rank[c_j, s_j]`: two directions oscillating at the
same mode span only a plane, so extra *modes*, not extra directions, buy dimension. With `⌊N/2⌋ − 1` usable
transverse patterns,

    deff_max(N) = 2⌊N/2⌋ :   N=3→2, 4→4, 5→4, 6→6, 7→6, 8→8, 9→8, 10→10, 12→12 .

This reproduces the whole catalogue at `N ≥ 4` and predicted the rest before the runs: `d = 7` with `N = 4` or
`5` cannot exceed `deff = 4` however long it runs (measured over 3·10⁴ trials), while `N ≥ 8` admits
`deff = 7`, and `deff = 11` needs `N ≥ 12`. It counts the vertical family's resonances and is not a theorem:
the three non-planar `N = 3` records are not of that type.

## 6. The rotating frame (`--omega`)

`--omega "w1,w2,…"` puts rates in the successive coordinate planes; `su:w1,…` closes the list with
`w_n = −Σw`; `g2:p,q[,r,…]` is the maximal torus of `g₂` — rates `(p, q, p+q)` on the Fourier planes of the
Fano 7-cycle, plus one free rate per leftover coordinate plane. Validated exactly: at `Ω = −1` the Lagrange
circle reappears as the mode-2 curve with an identical action and `|∇A| = 1.3e-15`.

The frame must be **calibrated** (§7) or the stored `χ*` is meaningless: `su:` puts its rates on coordinate
planes rather than Fano planes, and its calibration defect reads 3.0 at `d = 7`, 6.0 at `d = 9` and 4.0 at
`d = 11`, against 1e-15 for the `g₂` frame.

The associative 3-form is **not** a `d = 7` object. It lives on `span(e₀…e₆)` in every `d ≥ 7`, and any
rotation of the remaining coordinates fixes it pointwise, so `stab(φ) = g₂ ⊕ so(d−7)` — dimension **14, 15,
20** at `d = 7, 9, 11`, computed and pinned by `make test`. That is what makes `d ≥ 9` reachable.

`continue --param omega` walks `Ω = sΩ₀` under the same pseudo-arclength corrector as `α`, so folds are
traversed and `s = 0` is certified against the genuine inertial problem. Because the kinetic operator is
exactly quadratic in `Ω`, the frame derivative is one kernel rather than a new solver.

### 6.1 Choosing a frame

The search problem depends on `Ω` only through the **multiset of `|rates|`**, which collapses the parameter
box. Score cells by **distinct families**, not raw records: a frame parked on a continuous family resamples it
indefinitely (`g2:3,4` at `d = 7, N = 10` returns 21 records that are all one family). Degenerate cells
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
seeds 1, 2 and 11 merged) gives **122 records, 38 of them `deff = 11`**, every one isolated — `nullity = 6`,
exactly the gauge dimension of the 5-torus plus the time shift — and validating clean: no missing Morse data,
no near-collision, every return error below 1e-9, nothing inside the rigidity gate. The `χ*` champion is id
13, `A = 22.745582770`, `χ* = 81.3` (`tw_rel = 0.32`), certified at shift residual 1.8e-12, full-period
return 1.6e-12 and energy drift 1.8e-15.

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

`R³`'s cross product buys nothing (its 3-form is `SO(3)`-invariant), and **7 is the only odd dimension with a
proper reduction**. In the complex Fourier basis the jet moment is exactly

    A_k = i^{k(k−1)/2} · Σ over {n₁<…<n_k} ⊂ Z∖{0} with Σn = 0  of  V(n) · z_{n₁} ∧ … ∧ z_{n_k},
    V(n) = Π_{r<s} (n_s − n_r)   (Vandermonde)

— verified against the quadrature to 1e-15 in `make test`. So only **resonant** mode tuples contribute; the
modes must be **distinct**, with spread-out sets dominating; and `χ* ≠ 0` certifies `deff ≥ k`.

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

## 8. Files, resumability, tooling

* `catalog.bin` — little-endian, bit-exact doubles: `"HYPCHOR1"` then records (`RecHdr`, 152 bytes, then
  `int32 modes[nm]`, `double coef[nm·2·d]`, `double Lsv[d]`, `double pca[d]`, `char sym[]`). Records carry a
  trailing `double extra[]`: `[0]` `χ*`, `[1]` `χ*` over the jet scale, `[2]` the rung's jet order, `[3]`
  `‖A_k‖` over the jet scale, `[4]` the rigidity defect, `[5]` the layout, `[6]` the coefficients' residual,
  `[7]` the proven box radius (`prove --write`, 0 = none), then `d²` entries of `Ω` and, in layout 1, the
  certified state. Legacy files load unchanged. Accessors must read `extra.size()`, not the header's `nextra`.
* `catalog.bin.state` — seed and next trial index. Trials are deterministic functions of `(seed, trial)`, so a
  run resumed after Ctrl-C continues exactly, whatever the thread scheduling. Checkpoints every 30 s. Use a
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
* **Two records have a floor MPFR does not lift** — the same residual in double and at 30 digits, so there is
  no solution near that state at that `Ω`: `d7_n12_g2:61` at 6.0e-11 and `d7_n10_g2_16:111` at 2.8e-12, the
  only records above 3.7e-12. `:61`'s canonical frame has principal values 0.0984/0.0981, a relative gap of
  9e-7, leaving the stored `Ω = R Ω Rᵀ` wrong at ~1e-10 — the frame the record names is not the one its orbit
  solves. Solving for `Ω` alongside `Z` would settle it.
* **Calibrated starts do not work below `N ≈ 12`, for a structural reason.** A twisted loop needs modes in
  *arithmetic* resonance (`Σn = 0`); a start that converges needs the *dynamical* resonance `m₂/m₁ ≈ ω_k`.
  The two are incompatible at small `N` — the closest `1 + ω_j = ω_k` is off by 0.64 at `N = 8` and 0.37 at
  `N = 10`, but only 0.087 at `N = 12`; forcing the resonance at `d = 7, N = 8` gives 4 records against
  `hyper`'s 12. Twist comes from the nonlinear mode cascade, not the start — a falsifiable prediction for
  `N ≳ 12`.
* Double precision suffices: the action matches 40-digit MPFR to 1–2 ulp, and re-shooting rejected candidates
  in MPFR recovers nothing. The ~1e-11 floor on the shift residual is Lyapunov amplification, not roundoff.
  Newton–Krylov was measured and dropped (0.88×–1.57×, no predictor). The `hits` counter is saved only at
  checkpoints, so a kill between them loses ≤ 30 s of counts.

## 11. Roadmap

1. **Branch-switch off the `Ω` isola.** Of 18 branches followed from the `deff = 7` records of
   `d7_n10_g2.bin`, 5 reach `s = 0` and **every one arrives on a relative equilibrium**; 1 is a **closed
   isola** (the `χ* = 16.2` champion) never leaving `s ∈ [0.98, 1.06]`. The Morse index does cross (6↔7) at
   `s ≈ 0.9829` and `s ≈ 0.9986` — what `--depth` is for, and the one untried escape.
2. **Linear stability (Floquet)** — the largest remaining gap, since a Morse index is not linear stability.
   Multipliers come from the reduced monodromy `S⁻¹ DΦ_{T/N}`, whose symplectic defect is 1.3e-07 against
   9.4e+02 for `DΦ_T` taken directly. The shooting Jacobian is already computed; the missing piece is an
   unsymmetric eigensolver (Hessenberg + Francis QR, ~130 lines).
3. **The continuous families.** `g2:3,4` sits on a 2-parameter family at `A = 19.827310867` (nullity 6 = 4
   gauge + 2); `d = 9` and `d = 11` have none, so families are so far a `d = 7` phenomenon worth explaining.
4. **Mine the symmetry classes instead of guessing them.** `symmetry` emits exactly the text `--sym` consumes,
   so the loop closes: detect the groups that occur, keep those with a small fixed subspace, re-search inside.
5. **More cells.** `d = 13` needs `--omega g2:p,q,r₁,r₂,r₃` and `N ≥ 14`; `d = 12, N = 12` gave 0 records in
   8 537 trials at a quarter of the trial rate, which is the warning.
6. Better branch switching for degenerate bifurcations; multi-threaded `continue`; FFT synthesis for
   `K ≳ 200`; strong-force homotopy classes; mountain-pass between catalogued minimisers; relative
   choreographies with several curves. For `prove`: a Lohner-type wrapping control would cut the precision the
   large-`N` orbits need, and `α ≠ 1` only needs the potential's interval power wired through.

## 12. References

* C. Simó, *New families of solutions in N-body problems*, 3rd European Congress of Mathematics (2000).
* G. Minton, *Choreographies* (gminton.org/choreo.html).
* A. Chenciner, R. Montgomery, *A remarkable periodic solution of the three-body problem* (2000).
* A. Chenciner, J. Féjoz, *Unchained polygons and the N-body problem* (2009).
* D. Ferrario, S. Terracini, *On the existence of collisionless equivariant minimizers* (2004).
