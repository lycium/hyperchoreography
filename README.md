# hyperchoreography

A search engine for **N-body choreographies in arbitrary dimension**: N equal masses chasing each other along
one closed curve, `q_k(t) = q(t + kT/N)`, under the Newtonian (or any homogeneous `1/r^α`) potential.
Built to run for weeks on many cores, to be resumable, to never double-count, and to certify every
catalogued orbit against the true equations of motion — in double, and in MPFR to any number of digits.

```
make                  # native build (clang/gcc; -mcpu=native on arm64, -march=native elsewhere)
make test             # self-checks of every kernel (derivatives, symmetry, integrator, dedup, I/O)

./hyperchoreography search   --d 3 --N 4 --K 24 --threads 16 --minutes 600 --out d3n4.bin   # Ctrl-C any time; rerun to resume
./hyperchoreography continue --root circle --d 3 --N 4 --K 32 --covers 9 --depth 2 --out d3n4.bin
./hyperchoreography list     d3n4.bin                 ./hyperchoreography show   d3n4.bin --id 7
./hyperchoreography verify   d3n4.bin [--id 7]        ./hyperchoreography export d3n4.bin --id 7 --out orbit.csv
./hyperchoreography refine   d3n4.bin --id 7 --digits 100 --out d3n4_7.txt
./hyperchoreography merge    all.bin d3n4.bin other.bin [--min-rigid r --min-deff k]
./hyperchoreography extras   d3n4.bin                 ./hyperchoreography symmetry d3n4.bin
make gallery                                          # docs/index.html, served at lycium.github.io/hyperchoreography
```

Dependencies: a C++20 compiler; MPFR + GMP for `refine` (`make NOMPFR=1` drops them); on macOS, Accelerate
for LAPACK `dsyevd` above n = 64 (`make NOACCEL=1` drops it). Everything else — symmetric eigensolver,
pivoted LU, optimisers, Taylor integrator, MPFR wrapper — is in `src/`.

Naming: **dimension first, then bodies** (`d3n4.bin`); `list` sorts by (deff, N, action).

---

## 1. What is in the catalogue

`deff` is the dimension the motion actually occupies; `d` is the ambient dimension it was found in. The
budget of §5 caps `deff` at `2⌊N/2⌋`, and that ceiling has been reached exactly at `N = 6`, `7` and `10`.

Every record stores the **certified initial state** as well as the Fourier series rendering it: `h.ret_err`
is the state's residual (median 2.2e-15, worst 6.0e-11 over 1800 records), `extra[6]` the series', which the
`m²` weighting of the equations always makes looser. Most of these orbits are violently unstable — only 76
are action minima, median Morse index 20 — which costs precision, not existence (§4).

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
as are 13 of the `deff = 5` ones, so genuine choreographies, not merely relative ones, reach `deff = 6`.
No high-`deff` relative choreography found so far deforms back to an inertial one (§11, item 1).

Arbitrary-precision refinement works in **any** frame: the MPFR Newton carries `G = exp(2πΩ/N)` and threads
its Jacobian one column per task. The figure eight reaches **1e-71**, `d11_n12_g2:106` 4.8e-4 → **1.4e-48**.

## 2. The variational formulation

The loop is a real Fourier series with period 2π (fixed by the Kepler scaling `q → λq, t → λ^{(α+2)/2} t`):

    q(t) = Σ_m  c_m cos(mt) + s_m sin(mt),      c_m, s_m ∈ R^d,   1 ≤ m ≤ K,   m ≢ 0 (mod N).

Modes `m ≡ 0 (mod N)` only move the centre of mass, so they are excluded. Using the `Z_N` symmetry the
action of the whole configuration reduces to a functional of the single curve (overall factor N dropped):

    A[q] = ½ ∫ |q̇|² dt  +  ½ Σ_{k=1}^{N-1} ∫ |q(t) − q(t + 2πk/N)|^{−α} dt .

By Palais' symmetric criticality, critical points are choreographic solutions. The kinetic term is analytic
in the coefficients; the potential is a trapezoid quadrature with `M` nodes, spectrally accurate for
collision-free loops, needing only `⌊(N−1)/2⌋` shifts. Sampled data are structure-of-arrays with a doubled,
wrap-free layout, so every inner loop is contiguous and vectorisable. Value, gradient, Hessian–vector
products, the full Hessian (assembled from two DFT bins per mode pair) and `∂∇A/∂α` are closed-form.

In a frame rotating with `Ω ∈ so(d)` the ansatz is `q_j(t) = exp(Ωt) q(t + 2πj/N)`. The potential is
**unchanged** (rotations are orthogonal); only the kinetic term becomes `½∫|q̇ + Ωq|²`, a per-mode `2d × 2d`
block `π[[m²I − Ω², −2mΩ], [2mΩ, m²I − Ω²]]` in place of the scalar `π m²`.

At `(d, N, K) = (7, 7, 32)` the Hessian costs 355 µs against 5.6 ms for its eigendecomposition, which is the
dominant cost of a high-dimensional trial.

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

A Fourier critical point is never accepted on its own merits.

`recertify` re-derives a whole file's states and loops in three minutes; `verify <file>` with no `--id`
audits one. The Newton runs in double, but where an orbit's monodromy stalls it above `--mpfr-gate` (1e-12)
the solve is repeated in MPFR and rounded back (`d7_n11_g2:127`, 4.7e-12 → 1.3e-15). Where rounding is
amplified back up (`:142`, 9e-37 in MPFR, 3.7e-12 rounded) the double answer is kept.

* **ODE validation** — initial conditions are read off the series and integrated with a 22nd-order Taylor
  method over `T/N`; the shift residual must be ≤ `--ret-reject`. Its default 1e-1 is deliberately loose:
  the shooting Newton is the real gate, and tightening it to 1e-5 costs 2.5–4.5× in unique orbits per second.
* **Shooting Newton** — `Φ_{T/N}(Z) = G S Z` with `G = exp(2πΩ/N)`, solved to 1e-12 with Levenberg damping
  adapted on the residual (fixed damping abandoned 27 % of admissible candidates).
* **Fourier re-extraction** — the loop is re-sampled from the certified orbit over one `T/N` segment
  (exact, N× cheaper, immune to the error an unstable orbit accumulates over a full period).
* **k-fold covers** — a loop traversed k times is the same choreography, detected from the gcd of the
  significant modes and unwound.
* **Canonical frame** — coefficients rotated into the principal axes of `XᵀX`; `deff` is the number of
  non-negligible principal values, so a planar loop found in a 4-D search compares equal to its 2-D twin.
* **Rigidity defect** — a relative equilibrium moves rigidly, so every mutual distance is constant; `⌊N/2⌋`
  separations exhaust the test, and a common rotation cancels, so the answer is independent of `Ω`. These
  are trivial choreographies however large their `deff`, and they **dominate** the high-`d` search. The
  stored defect is trimodal: a cluster at `≤ 2e-11`, a second at `1.3e-6 … 4.0e-6`, then a **4300× empty
  gap** to the first genuine orbit at `1.7e-2`. `--min-rigid` defaults to **1e-4**, above the second
  cluster; its members give themselves away with `tw_rel = 0` and `deff = d − 1`. `merge --min-rigid r
  --min-deff k` re-gates a file on disk, so a threshold change costs seconds rather than a re-harvest.
* **`deff` in a rotating frame** is not the loop's own rank: `exp(Ωt)` raises it when a fixed point sweeps a
  circle and lowers it when a circularly polarised mode at rate `−w` becomes a linear oscillation. At
  `Ω = 0` the two agree exactly, so no inertial record changes.
* **Equivalence** — two loops are the same choreography iff they agree up to time shift, time reversal,
  `O(d)` and relabelling. Candidates are filtered on invariants, then compared by Procrustes distance. One
  case defeats that on purpose: a **continuous family**, whose points are genuinely different loops but one
  object. The action is exactly constant along a family, so matching it to round-off folds them into a
  single record with a hit count; `nullity` above the gauge dimension is the independent signature.
* **Morse index / nullity** of the full Fourier Hessian. The exact gauge directions — the time shift and the
  rotations commuting with `Ω` — are lifted out by `H + σGGᵀ` rather than classified by magnitude. That
  centraliser is the kernel of `ad_Ω` on `so(d)`, **not** the coordinate generators `E_ab` tested one at a
  time: when the frame's planes are not coordinate planes (the `g2` torus) no single `E_ab` commutes, and
  the per-generator test leaves the torus directions in the spectrum as numerical zeros. A stored
  `nullity = 0` is impossible — the time shift alone forces 1 — and always meant this bug.
* **The calibration twist `χ*`** (§7) and **the symmetry group of the loop**, detected rather than assumed
  (`symmetry`): for fixed `(ε, θ)` the best `R` is the Procrustes fit, so the group is the zero set of the
  loop's *self*-distance. The relative squared residual bottoms out near 1e-15, which sets the default
  tolerance at 1e-6 on the r.m.s.; a loop for which every shift works is reported as `S¹`.

A generator of `--sym` is `(ε, θ = 2πp/q, R ∈ O(d))` imposing `q(εt + θ) = R q(t)`, written
`t+p/q s[±i,±j,…]` (signed permutation) or `r(i,j,p/q)` (plane rotation). The figure-eight class is
`"t+1/2 s[-1,2]; t-0 s[1,-2]"` — its half-period map is `diag(−1,+1)`, not `−I`; `s[-1,-2]` is the
*circle's* class. Named classes: `cyc:p` cycles the `⌊d/2⌋` coordinate planes, `fano:p` is the `d = 7` Fano
7-cycle. Both are cyclic with trivial core, so they do not cap `deff` (§5) — their value is cost, the fixed
subspace of `cyc:1` at `d = 6, N = 5` having `r = 52` of `n = 156` parameters.

## 5. Where the high-dimensional choreographies come from

Random low-mode starts collapse onto the Lagrange circle in `d ≥ 3` (58 000 trials in 3-D gave circle and
eight, nothing else). The structured starts (`--starts`) and `continue` are what produce spatial solutions:
**torus** (rotations in ⌊d/2⌋ orthogonal planes), **vertical** (a rotating circle plus one transverse
oscillation — the unchained-polygon ansatz), **hyper** (one near-resonant transverse mode per pattern
`k = 2…⌊N/2⌋`, circularly polarised in pairs; this produced every `deff ≥ 5` record), **inplane** (the
in-plane resonances below), **fano** (`d = 7`), and **kick** (`--seed-from`, a catalogued solution embedded and pushed along its softest Hessian directions).

Linearising the rotating N-gon gives the transverse frequencies in closed form,
`ω_k² = Σ_{l≠0} (1 − cos 2πkl/N)/d_l³`, `d_l = 2R sin(πl/N)`, so inertial choreographies of vertical type sit
at rational resonances `m₂/m₁ ≈ ω_k/ω_N` with pattern `k ≡ m₂·m₁⁻¹ (mod N)`. N=4, k=2 gives 1.2156 → the
(5,6) solution; N=5, k=2 gives 1.3277 → the (7,9) solution. For N=3 every transverse mode is a tilt, which is
why nothing spatial exists there.

The **in-plane** block is the other half of the same linearisation (`ngon_inplane_freq`): a quartic in the
same sums `C_k = Σ_{p≠0} g_p cos(2πkp/N)`, of which `C_0 − C_1 = 1` *is* the radius equation. Its real roots,
which exist for every `k ≥ 2` at `N ≥ 7`, give a second resonance family `m₂/m₁ = 1 ± ν̂_k` with
`k ≡ m₂·m₁⁻¹ − 1 (mod N)`, seeded by `--starts inplane`. A/B against the `vertical,hyper,torus` mixture at 20
min per arm: fewer records (21 against 25 at `d = 3, N = 9`) but roughly half of them ones the mixture never
found, so it belongs **in** the mix, not in place of it.

**Why `deff` saturates.** Two exact facts bound it.

*The reflection stratum is attracting.* Let `σ_r ∈ O(d)` fix `e_1…e_r`. Since `A∘σ_r = A`, the set of loops
in `R^r` is invariant under L-BFGS *and* under the Newton–LM step: a `deff = r` critical point can never be
left, only entered. At such a point the Hessian is exactly block diagonal and `∇_⊥A` is odd in `ξ_⊥`, so
Newton contracts transversally like `‖ξ_⊥‖³` **whatever the sign of `H_⊥`**. Measured: phase-1 L-BFGS keeps
~50 % of `d`-dimensional starts `d`-dimensional and the Newton step destroys 99 % of those; kicking along
the softest, softest-negative, or random directions raises `deff` in **0 of 816** trials. **Starts are
everything**, and symmetry cannot substitute — forcing `deff = d` needs a trivial core, which makes the
group cyclic or dihedral, whose real irreps have dimension ≤ 2, so `Fix(G)` always contains planar loops
(500 `--sym random` draws at `d = 5, N = 4` and `d = 7, N = 6` produced not one loop with `deff > 2`).

*The budget.* A vertical-family loop has `deff = 2 + Σ_j rank[c_j, s_j]`: two directions oscillating at the
same mode span only a plane, so extra *modes*, not extra directions, buy dimension. With `⌊N/2⌋ − 1` usable
transverse patterns,

    deff_max(N) = 2⌊N/2⌋ :   N=3→2, 4→4, 5→4, 6→6, 7→6, 8→8, 9→8, 10→10, 12→12 .

This reproduces the whole earlier catalogue and predicted the rest before the runs: `d = 7` with `N = 4` or
`5` cannot exceed `deff = 4` however long it runs (measured over 3·10⁴ trials), while `N ≥ 8` admits
`deff = 7`, and `deff = 11` needs `N ≥ 12`.

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
box; cells are indexed by that multiset. Sweep small distinct rates at ~45–90 s per cell, always including
an all-`r = 0` control, score **distinct families** rather than raw records (a frame parked on a continuous
family resamples it indefinitely — `g2:3,4` at `d = 7, N = 10` returns 21 records in 45 s that are all one
family), then re-run the top few at 3–4× the budget on a second seed.

Three regularities hold at every `d` tested: **degenerate** cells (a repeated rate) collapse and sometimes
cap `deff` outright; the cell carrying the two largest rates returns few or no records; and the smallest
distinct rates win. The ranking **does not transfer between `N`** — `{2,3,5}` is seventh of eight at
`N = 10` and first at both `N = 11` and `N = 12`, while `N = 10`'s winner `{1,6,7}` drops to seventh — so
every new cell needs its own sweep. It is stable in budget and seed *within* a cell, so it is the `N` axis
that destabilises it. Short sweeps overstate long-harvest rates: over 25–30 min the three `d = 7` cells land
within 30 % of each other, though a 150 s sweep made `N = 11` look 1.3× better.

Harvest frames: `"1,2"` at `d = 5`; `su:1,2` at `d = 6, N = 6, 7` and `su:1,3` at `N = 8`; `g2:1,6` at
`d = 7, N = 10`; `g2:2,3` at `N = 11, 12`; `g2:1,3,2` at `d = 9`; `g2:1,2,4,5` at `d = 11`.

### 6.2 The leftover rates buy the last dimension

Setting the leftover rate to zero leaves the `d = 7` frame sitting inside `R^d` with a fixed subspace. At
`d = 9` that reaches `deff = 8` but **never 9** — at 45 s it appeared to cap at 7, which invited the tidy
conclusion that the frame must rotate every direction the loop occupies, and 4× the budget refuted it (the
same frames reach 8). At `d = 8`, where `so(1) = 0` and a `g₂` frame has no leftover rate at all, `deff = 8`
is still reached. So a direction the frame does not rotate **can** be occupied.

What survives at both budgets and every rate set is narrower: **no `r = 0` run has ever reached the top
dimension.** At `d = 11` the same control, at 4× the sweep budget:

| frame | `deff` 11 / 10 / 9 @ 90 s | @ 6 min (4×) |
|---|---|---|
| `g2:1,2,4,0`, one rate zeroed | 0 / 1 / 6 | 0 / 4 / 12 |
| `g2:1,2,0,0`, both zeroed | 0 / 0 / 0 | no records at all |

The extra budget lifts the histogram exactly as at `d = 9` and still produces no `deff = 11` orbit. With
both rates zeroed the frame is the bare `d = 7` torus in `R¹¹` with a fixed 4-space, and it returns nothing
while running at the **highest** trial rate of any cell — it dies cheaply on the `deff` pre-filter. The
likely mechanism is mode-naming: the frame zeroes the kinetic eigenvalue `π(m − w)²` of the mode `m = w`, so
a rate pays only by naming a mode that exists, and at `N = 10` the modes are `m ≢ 0 (mod 10)`, so `r = 0`
names nothing. That is a statement about *rate*, not a proven structural gate.

### 6.3 `d = 11`

`--omega g2:p,q,r₁,r₂` with `N = 12`, the cheapest `N` the budget allows for `deff = 11`; `K = 24`,
`--starts hyper`, ten cells at 90 s scored by distinct `deff = 11` records:

| rates | `deff` 11 | 10 | 9 | best `χ*` | trials/s |
|---|---|---|---|---|---|
| {1,2,3,4,5} (`g2:1,2,4,5`) | **8** | 1 | 4 | 81.3 | 16.8 |
| {1,2,3,5,6} (`g2:1,2,5,6`) | 2 | 2 | 3 | 133.0 | 16.5 |
| {1,2,3,4,7} (`g2:3,4,1,2`) | 2 | 0 | 0 | 18.0 | 14.9 |
| {1,2,3,6,7} (`g2:1,6,2,3`) | 1 | 0 | 3 | 45.3 | 13.9 |
| {1,2,3,4,6} (`g2:1,2,4,6`) | 0 | 3 | 6 | — | 14.2 |
| {1,2,3,5,7} (`g2:2,5,1,3`) | 0 | 0 | 1 | — | 12.9 |
| {2,3,4,5,6} (`g2:2,3,4,6`) | 0 | 0 | 6 | — | 14.6 |
| {1,2,3,3,4} (`g2:1,2,3,4`), degenerate | 0 | 0 | 1 | — | 14.0 |

Re-run on a second seed at 3.3× the budget the ranking holds — 11, 5 and 3 `deff = 11` records for the top
three — so `g2:1,2,4,5` is the frame to harvest in.

One hour on 18 cores (57 327 trials, seeds 1, 2 and 11 merged) gives **122 records, 38 of them
`deff = 11`**, every one isolated — `nullity = 6`, exactly the gauge dimension of the 5-torus plus the time
shift — and validating clean: no missing Morse data, no near-collision, every return error below 1e-9, and
nothing inside the rigidity gate. The `χ*` champion is id 13, `A = 22.745582770`, `χ* = 81.3`
(`tw_rel = 0.32`), certified at shift residual 1.8e-12, full-period return 1.6e-12 and energy drift
1.8e-15.

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

`R³`'s cross product buys nothing — its 3-form is invariant under all of `SO(3)` — and **7 is the only odd
dimension with a proper reduction**. In the complex Fourier basis the jet moment is exactly

    A_k = i^{k(k−1)/2} · Σ over {n₁<…<n_k} ⊂ Z∖{0} with Σn = 0  of  V(n) · z_{n₁} ∧ … ∧ z_{n_k},
    V(n) = Π_{r<s} (n_s − n_r)   (Vandermonde)

— verified against the quadrature to 1e-15 in `make test`. So only **resonant** mode tuples contribute; the
modes must be **distinct**, with spread-out sets dominating; and `χ* ≠ 0` certifies `deff ≥ k`.

`χ*` is silenced by exactly two things. **No resonance** — a relative equilibrium sits at the rotation rates
of `Λ`, generically non-resonant, so every tuple is empty. **Time reversal at the wrong parity of `k`** — if
`q(−t + θ) = R q(t)` then `R·A_k = (−1)^{k(k−1)/2} A_k`, which at `k = 3` is `−1`, so a reversal whose `R`
fixes `A₃` annihilates the twist however resonant the modes are, while at `k = 4` it costs nothing. The data
agrees sharply: 13 of 13 records at `d = 8` carry twist, against 24 of 98 at `k = 3`.

Two cautions, both learned the hard way:

* **`χ*` is not invariant under the rotating-frame gauge.** `exp(Ωt) q(t)` is unchanged when a rate moves by
  `N` and the modes compensate, so `(Ω, q)` is a gauge choice and `χ*`, read off `q` alone, moves with it —
  measured at a factor of **21** on one orbit whose two records agreed on every other stored scalar. Twist
  is comparable *within* a frame, never across frames, and no canonical frame is defined yet.
* **`χ*` ranks orbits only downstream of the rigidity gate.** A *resonant* relative equilibrium maximises
  exactly the pairing `χ*` measures — 20.0 against 4e-14 for a non-resonant one — so a rigid record once
  looked like the project's twist champion at 446.5. Across the gated catalogue there is no general
  twist/rigidity correlation; this is one sharp trap, not a trend.

## 8. Files, resumability, tooling

* `catalog.bin` — little-endian, bit-exact doubles: `"HYPCHOR1"` then records (`RecHdr`, 152 bytes, then
  `int32 modes[nm]`, `double coef[nm·2·d]`, `double Lsv[d]`, `double pca[d]`, `char sym[]`). Records carry a
  trailing `double extra[]`: `[0]` `χ*`, `[1]` `χ*` over the jet scale, `[2]` the rung's jet order, `[3]`
  `‖A_k‖` over the jet scale, `[4]` the rigidity defect, `[5..7]` reserved, then `d²` entries of `Ω`.
  Legacy files load unchanged. Accessors must read `extra.size()`, not the header's `nextra`.
* `catalog.bin.state` — seed and next trial index. Trials are deterministic functions of `(seed, trial)`, so
  a run resumed after Ctrl-C continues exactly, whatever the thread scheduling. Checkpoints every 30 s.
  Use a different `--seed` per machine and `merge` the results; the better-resolved record wins and hit
  counts add.
* `src/` — `action.hpp` (basis, kernels, symmetry), `calib.hpp` (the ladder), `g2.hpp` (Fano frame and
  torus), `optim.hpp`, `taylor.hpp`, `invariants.hpp`, `catalog.hpp`, `search.hpp`, `continue.hpp`,
  `mpreal.hpp`, `linalg.hpp`, `main.cpp`, `tests.cpp`.
* `tools/gallery.py` — the whole catalogue as one self-contained page, standard library only, ~20 s.
  Organised by **effective dimension**, which the page calls `d` (the ambient space is a search setting, not
  a property of the orbit). A selection sits at the top — `--heroes N`, `--hero`/`--no-hero FILE#ID` — and
  the per-dimension sections below it start closed, so nothing else is built until it is asked for.
  `--split` writes one page per dimension behind an index instead. Sampling, tile clock and the selection
  rule are each derived where they are implemented.
  `make gallery` writes `docs/index.html`, which GitHub Pages serves at
  <https://lycium.github.io/hyperchoreography/>. A `deff = 9` orbit and a circle look identical in any two
  coordinates you pick, so each tile is drawn in the record's **own principal frame** and shows a 3-D
  orthographic shadow, one panel per principal plane in three states (curve, line, grey rule — so
  `deff = 2×curves + lines` is countable), a bar meter of `√(λₖ/λ₁)`, and a mutual-distance ribbon that is
  flat exactly for a relative equilibrium. Measured against plotting the first two coordinates: **2 of 316**
  tiles read as a near-circle, against **80 of 333**. Only body 0 is stored; the rest are rebuilt from
  `body_j(f) = Gʲ·body_0((f + jS/N) mod S)`.

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

* **Some `continue` branches stall before α = 1.** Two things are mixed: branches that genuinely end in a
  collision (which the code should recognise from the shrinking minimum separation and report as such), and
  switching problems where a degenerate *pair* of eigenvalues crosses — rotational symmetry makes the
  bifurcating set a circle of solutions, and the bordered system needs an explicit phase condition.
* **Two records have a floor MPFR does not lift** — the same residual in double and at 30 digits, so there is
  no solution near that state at that `Ω`: `d7_n12_g2:61` at 6.0e-11 and `d7_n10_g2_16:111` at 2.8e-12, the
  only records above 3.7e-12. `:61`'s canonical frame has principal values 0.0984/0.0981, a relative gap of
  9e-7, which leaves the stored `Ω = R Ω Rᵀ` wrong at ~1e-10 — the frame the record names is not the one its
  orbit solves. Solving for `Ω` alongside `Z` would settle it.
* **`--sym random` is wasteful and, in `d ≥ 5`, useless**: 60–75 % of draws are rejected and the survivors
  land on the circle. The named classes `cyc:p` and `fano:p` are the useful ones, but a curated table per
  `(d, N)` is still not wired in.
* **Calibrated starts do not work below `N ≈ 12`, for a structural reason.** A twisted loop needs modes in
  *arithmetic* resonance (`Σn = 0`); a start that converges needs the *dynamical* resonance `m₂/m₁ ≈ ω_k`.
  The two are incompatible at small `N` — the closest `1 + ω_j = ω_k` is off by 0.64 at `N = 8`, 0.37 at
  `N = 10`, and only 0.087 at `N = 12`. Measured A/B at `d = 7, N = 8`: forcing the resonance gives 4
  records against `hyper`'s 12. Twist in the found orbits comes from the nonlinear mode cascade, not the
  start. This is a sharp, falsifiable prediction for `N ≳ 12`.
* **Close-approach loops** (minsep < 0.05) get 500+ modes; certification is fine but the Fourier
  representation is then a poor basis — parallel shooting would suit them better.
* **Six records' stored `deff` disagrees with the gallery's pooled PCA**, which `make gallery` reports —
  five `d5-6_n6` members of one near-rigid family whose third principal value at 2e-4 falls inside
  `inertial_deff`'s absolute cut and outside the gallery's relative one, plus `d11_n12_g2#0`. All certify to
  1e-13 or better; the two measures differ only for rotating records, and neither is established as right.
* Double precision suffices: the action matches 40-digit MPFR to 1–2 ulp, and re-shooting rejected
  candidates in MPFR recovers nothing. The ~1e-11 floor on the shift residual is Lyapunov amplification,
  not roundoff. Newton–Krylov was measured and dropped (0.88×–1.57×, no predictor).
* The `hits` counter is only saved at checkpoints; a kill between them loses ≤ 30 s of counts.

## 11. Roadmap

1. **Branch-switch off the `Ω` isola.** Following both directions from each `deff = 7` record of
   `d7_n10_g2.bin` gives 18 branches: 5 reach `s = 0` and **every one arrives on a relative equilibrium**;
   1 is a **closed isola** (the `χ* = 16.2` champion) that never leaves `s ∈ [0.98, 1.06]` and returns to
   its own starting orbit; the rest run to the `s` bounds or stick. So no high-`deff` relative choreography
   deforms to an inertial one. The Morse index genuinely crosses (6↔7) at `s ≈ 0.9829` and `s ≈ 0.9986`,
   which is what `--depth` is for — the one untried escape.
2. **Linear stability (Floquet).** The largest remaining gap: the catalogue reports a Morse index, which is
   not linear stability. Multipliers come from the reduced monodromy `S⁻¹ DΦ_{T/N}`, whose symplectic defect
   is 1.3e-07 against 9.4e+02 for `DΦ_T` taken directly. The shooting Jacobian is already computed; the
   missing piece is an unsymmetric eigensolver (Hessenberg + Francis QR, ~130 lines).
3. **The continuous families.** `g2:3,4` sits on a 2-parameter family at `A = 19.827310867` (nullity 6 =
   4 gauge + 2). `d = 9` and `d = 11` have none — every record is pure gauge — so families are so far a
   `d = 7` phenomenon, which is itself worth explaining. Folding is lossy, so re-search with looser
   tolerances rather than expecting to recover members.
4. **Mine the symmetry classes instead of guessing them.** `symmetry` emits exactly the text `--sym`
   consumes, so the loop closes: detect the groups that occur, keep those with a small fixed subspace, and
   re-search inside them.
5. **The in-plane resonances in `continue`.** `ngon_inplane_freq` drives `--starts inplane` (§5), but the
   deterministic bifurcation list it gives `continue` is still not wired in.
6. **More cells.** `d = 13` needs `--omega g2:p,q,r₁,r₂,r₃` and `N ≥ 14`; `d = 12, N = 12` gave 0 records in
   8 537 trials at a quarter of the trial rate, which is the warning.
7. Better branch switching for degenerate bifurcations; multi-threaded `continue`; FFT synthesis for
   `K ≳ 200`; strong-force homotopy classes; mountain-pass between catalogued minimisers; rigorous interval
   Newton–Kantorovich verification; relative choreographies with several curves.

## 12. References

* C. Simó, *New families of solutions in N-body problems*, 3rd European Congress of Mathematics (2000).
* G. Minton, *Choreographies* (gminton.org/choreo.html).
* A. Chenciner, R. Montgomery, *A remarkable periodic solution of the three-body problem* (2000).
* A. Chenciner, J. Féjoz, *Unchained polygons and the N-body problem* (2009).
* D. Ferrario, S. Terracini, *On the existence of collisionless equivariant minimizers* (2004).
