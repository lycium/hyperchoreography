# Mathematics, certification, and performance audit

Audit date: 2026-09-05. Branch: `audit/math-physics-performance`, based on `02a558f`.

The engine has a strong reusable core: one Taylor recurrence supports ordinary, arbitrary-precision,
and interval computations, and the reduced shooting proof uses conservation laws to close its missing
equations. Several implementation details nevertheless invalidated claimed interval bounds or discarded
potentially distinct orbits. Those came before optimization in this audit.

This is a substantial tested correction and research foundation, not a claim that every line is now optimal,
every record is new, or the entire catalogue has been re-proved. The generated `docs/` gallery and its
source generator are included; `exposition/` is outside the audit's scope. Historical paper/PDF
claims and proof transcripts have not been comprehensively revised.

The sibling `../hyperchoreography-math-review` is a separate Git worktree on branch `math-review` at
`9b743bb`, containing older proof/paper work. The audit's base also includes later portability changes.
The sibling worktree was inspected read-only and left intact.

## 1. Main findings and changes

| Area | Finding | Change / evidence |
|---|---|---|
| Interval trigonometry | Reordering already-rounded endpoints can turn lower bounds into upper bounds; binary64 angle tests miss extrema. | Outward evaluation at both endpoints plus an interval `π/2` critical-point test. Endpoint/interior regressions. |
| Interval magnitudes | Converting a negative lower endpoint toward `+∞` underestimates its absolute value, including underflow to zero. | Away-from-zero conversion; explicit between-double and subnormal regressions. |
| Proof norm bounds | Ordinary double sums/products and an empirical safety factor around `exp` were not outward bounds. | Upper-rounded positive arithmetic and MPFR `expm1` for the Grönwall enclosure, inclusion and closure norms. |
| Action enclosure | The remainder used pair coefficients of an order the primal recurrence had not computed. | A consistent order-p primitive remainder; tested against an independent analytic radial Kepler action at orders 2, 3, 5 and 8. |
| Frame conservation laws | Rates declared congruent by floating subtraction need not be exactly congruent as real numbers. | Canonical logarithms modulo N and exactly identical/opposite snapped dyadic rates. A `0.1/5.1`, N=5 regression checks generator commutation. |
| Printed certificates | Nearest-rounded decimal endpoints need not enclose the internal interval. | Directed decimal formatting. Old proof markers require recomputation. |
| Orbit identity | Equal action/energy/separation or equal mode powers were sufficient to discard a record. Neither establishes equivalence. | Require a full-coefficient alignment and compatible frames; retain possible family members and unmatched high modes. |
| Ordinary integration | Stalls, collisions or step-limit exhaustion could appear as partial success. | Explicit failure propagation, finite/shape/tolerance checks, and correct zero-time dense output. |
| Rotating searches | The early rank test used the rotating loop rather than the inertial motion; cover unwinding retained the old frame. | Early rank filter and Fourier-gcd unwinding now apply only in the inertial case. |
| Continuation | A zero action gradient could be accepted while gauge/arclength constraints remained unsatisfied. | Check all corrector residuals, including after the final iteration. |
| Matrix exponential | Scaling used maximum entry rather than a matrix norm; high precision had a fixed series-length ceiling. | Row-sum scaling and a precision-dependent term count; a 64-dimensional analytic rotation regression. |
| Linear algebra | Inversion refactored the same matrix for every right-hand side: O(n⁴). | Reusable pivoted LU gives O(n³); MPFR updates use in-place arithmetic. Nonfinite eigensolver inputs/failure are no longer silent. |
| Allocation / concurrency | Shooting workspaces were rebuilt per iteration; proof matrix products constructed temporary intervals; expensive metadata held the catalogue lock. | Reuse worker scratch, mixed double–interval operations, metadata outside the lock with a second duplicate check, and synchronized progress snapshots. Pool exceptions now propagate after all workers finish. |
| Persistence / attribution | Truncation could be mistaken for EOF; failed writes or missing residuals could damage record selection. Prior-work provenance was absent. | Checked record lengths/writes, preserve the original on write failure, rank unknown residuals conservatively, and retain arXiv orbit citations through serialization and replacement. |

Relevant tests are in `src/tests.cpp` and `src/audit_tests.cpp`. The initial adversarial tests failed on
the original implementation while its existing test suite passed. The new cases therefore test gaps in
coverage, not just the same numerical identities at a different tolerance.

## 2. What physical problem is implemented?

For unit masses and positive α, the inertial equations and energy are

```text
Q̈_i = −α Σ_(j≠i) (Q_i−Q_j) / |Q_i−Q_j|^(α+2)
E = ½ Σ_i |V_i|² − Σ_(i<j) |Q_i−Q_j|^(−α).
```

The signs, α factor, pair accumulation, energy and Fourier action derivatives agree with this model in
the checked kernels. In arbitrary dimension, α=1 is the mathematical extension of the three-dimensional
Newton potential; it is not the Gauss-law gravitational potential for every d. Unequal masses are not
implemented. Importing the unequal-mass portion of an external dataset into the present equal-mass engine
would change the equations, not just the input format.

The scaling is `Q → λQ`, `t → λ^((α+2)/2)t`, `V → λ^(−α/2)V`. For a collision-free periodic or relative
periodic orbit whose virial boundary term cancels, `2<K> = α<U>`. The action stored here is per body over
the normalized `2π` window. For α=1, `E_total = −N A/(6π)`. These normalizations matter when comparing
published orbits with different periods.

### A real dimension bound, and a seed-family heuristic

In the center-of-mass frame let

```text
W = span{Q_i(0), V_i(0) : i=1,…,N}.
```

Acceleration is a linear combination of position differences. The ODE therefore preserves W: the
components normal to W remain zero by uniqueness of the collision-free flow. The two center-of-mass
relations give `dim W ≤ 2(N−1)`. Hence

```text
effective inertial dimension ≤ min(ambient dimension, 2(N−1)).
```

By contrast, `2 floor(N/2)` counts one vertical/hyper seed family's transverse patterns. It is not an
upper bound for all choreographies, not a proof that d=11 needs N≥12, and not a proof that d=13 needs
N≥14. Those are sensible parameters for the current recipe, not necessity theorems. Embedding an orbit
in a larger ambient space alone cannot raise its true dimension.

## 3. The twisting frame: three conditions that must not be conflated

Write `τ=2π/N`, `G=exp(τΩ)`, and let S send body j to body j+1. The ansatz is

```text
Q_j(t) = exp(Ωt) q(t+jτ),       q(t+2π)=q(t).
F(Z) = Φ_τ(Z) − GSZ.
```

Flow equivariance under rotations and equal-mass permutations gives, at a zero of F,

```text
Φ_(kτ)(Z) = (GS)^k Z,       Φ_(2π)(Z) = G^N Z,
Q_j(t) = G^(−j) Q_0(t+jτ).
```

Consequently:

1. `F=0` establishes relative choreography, with a common curve in the rotating description.
2. `G^N=I` establishes inertial `2π` periodicity. Integer frame rates are sufficient. Noninteger rates
   do not justify that assertion; a longer return or identity on the occupied span can still occur.
3. A single common inertial curve is an additional property. `G=I` is sufficient, but not necessary if
   extra spatiotemporal symmetry identifies the rotated paths. Integer rates alone do not prove it.

The certified frame is stated in its own axes, with exact dyadic rates printed in hexadecimal. Its
snapping is a choice of a nearby exact problem, not a proof that a numerically diagonalized input Ω has
exactly those rates. Bounding the discrepancy from the original input is separate from proving existence
for the explicitly chosen frame. Likewise, numerical `exp(2πΩ)≈I` in the gallery is only a rendering test.

### Why rotation changes which seeds work

With the covariant derivative `D_Ω = d/dt + Ω`, the potential is unchanged and the kinetic term is
`½∫|D_Ω q|²`. Its mode-m Hessian block is

```text
π [ m²I−Ω²    −2mΩ ]
  [  2mΩ      m²I−Ω² ].
```

On a plane of rate w the circular polarizations have kinetic eigenvalues `π(m±w)²`. A twisting frame
can make previously costly oscillations cheap. This is a concrete mechanism for different attraction
basins and accessible resonances, without asserting that rotation is physically necessary for high dimension.
At exact `m±w=0`, allowed modes can also lose kinetic coercivity: expanding a zero-kinetic shape sends its
positive potential contribution toward zero. Search must distinguish genuine critical points from large-scale,
small-gradient escape. Relative error and rigidity gates help, but an interval root proof is stronger.

Any skew Ω is admissible for the isotropic dynamics. A G₂-preserving choice restricts the geometry of a
chosen 3-form, not the force law. Embedding that form in larger d supplies a stabilizer algebra
`g₂ ⊕ so(d−7)` and useful free frame rates; it does not certify motion in those extra directions.

The current `χ*` uses ordinary jets of q and is frame-gauge dependent. A promising alternative is a moment
of the actual inertial curve, or a covariant-jet pairing `ψ(q,D_Ωq,…,D_Ω^(k−1)q)` when Ω preserves ψ.
The latter follows directly from differentiating `exp(Ωt)q(t)` and pairing with the preserved form. Any
proposed ranking still needs checks for gauge, time reversal, rigid motion and period normalization.
`calib_max` is numerical multistart optimization, not a verified global maximum; its L²-product scale is
not a general Hölder/Hadamard bound for k>2. A nonzero k-form moment can indicate rank at least k, not d.

## 4. Integration, shooting, and Newton

There is no Bulirsch–Stoer implementation in `src/`. The code does not use modified-midpoint sequences
and extrapolation. It is a direct high-order Taylor method. This is consistent with the established
[Taylor-method literature](https://www.maia.ub.es/~angel/taylor/taylor-1.4/taylor.pdf), rather than a reason
to replace the integrator by name alone.

For a pair, put `s=|Q_i−Q_j|²` and `w=s^β`, `β=−(α+2)/2`. In Taylor coefficients,

```text
s_k = Σ_(j=0)^k D_j · D_(k−j),
w_0 = s_0^β,
k s_0 w_k = Σ_(j=0)^(k−1) [β(k−j)−j] s_(k−j) w_j,
X_(k+1) = V_k/(k+1),       V_(k+1) = acceleration_k/(k+1).
```

The symmetry in the s convolution and equal/opposite pair force updates are appropriate. `Tangent`
differentiates these recurrences; it must index the primal coefficient table using the primal's order,
not its own potentially shorter order. The recurrence and tangent tests cover this structure.

The ordinary step-size selector is a high-order coefficient heuristic, not a certified error bound and
not a symplectic method. Small energy drift is useful diagnostic evidence, not an existence proof or a
guarantee of long-time phase accuracy. For searching unstable periodic orbits, accurate short return maps
are more relevant than an attractive full-period animation.

There are two distinct Newton calculations:

- Fourier search solves `∇A=0` with an analytic Hessian and eigenbasis LM filtering. This finds saddles
  as well as minima, but within a finite Fourier/quadrature approximation. Near zero eigenvalues do not
  automatically mean physical bifurcations; gauge directions and truncation errors must be separated.
- Shooting solves `F=0` by central-difference Jacobians and damped normal equations. Each Jacobian needs
  roughly twice the state dimension in integrations. Step-size noise and instability can contaminate
  finite differences, and `JᵀJ` squares the conditioning before damping.

The next solver experiment should share primal Taylor coefficients across tangent columns and compare
analytic J/Jv with central differences over several perturbation sizes. Solve the augmented LM problem
`min ||[J; sqrt(μ)I]δ + [F;0]||` by QR, with an SVD/rank-revealing option near bifurcations. LAPACK's
[DGELS documentation](https://netlib.org/lapack/lapack_routine/dgels.f) describes QR/LQ least squares and
explicitly warns that exact-rank detection is insufficient near rank deficiency. This is a proposed
experiment, not a speedup measured in this branch.

For very unstable/long-period solutions, use multiple shooting over shorter segments, including the GS
boundary condition and explicit phase/gauge constraints. Simply raising precision or replacing Taylor by
B-S does not remove an ill-conditioned long return map.

## 5. What the revised existence proof establishes

The approximate center and slice Q may be numerical constructions. Treating their stored values as exact
dyadic constants is legitimate: subsequent interval inequalities, not the accuracy of the approximate
inverse alone, do the certification work.

The proof selects a slice `Z=Z0+Qu`, drops k equations, and validates the reduced map on `u∈[-r,r]^m`.
It encloses a point residual and the derivative over the whole box. With an approximate inverse Y it checks

```text
||Y F_reduced(0)||∞ + ||I−Y[DF_reduced](B)||∞ r < r,
||I−Y[DF_reduced](B)||∞ < 1.
```

The dropped equations are not assumed away. Momentum components fixed by G, energy, and angular momentum
components commuting with exact G are conserved by both Φ and GS. If the retained equations agree, the
remaining state difference is supported on the dropped coordinates. An interval mean-value matrix of
those invariants on the hull of the two states is verified nonsingular, forcing the missing difference
to vanish. Exact frame commutation is therefore essential, not a cosmetic symmetry optimization.

The flow enclosure uses a Taylor polynomial from the tight initial box and an order-p remainder over a
rough enclosure W. A strict high-order inclusion supplies a collision-free segment enclosure; failure
shrinks the step. Tangents use a Grönwall bound with outward-rounded growth. The accumulated action is
another primitive with its own consistent remainder, not merely a quadrature of the approximate path.

The span check proves positive definiteness of a Gram matrix assembled from enclosed centered positions
at observed times. The non-rigidity check finds a pair distance with disjoint enclosures at two times.
Both are reported separately: successful shooting existence need not automatically establish either.
Uniqueness is inside the reported slice box, not global uniqueness of a family or the entire orbit space.

### Recomputed examples

These runs used `--force` without `--write`; existing catalogues were not rewritten to claim new certificates.

| Record | Result | Slice / gauges | Radius | Contraction | Closure | Validated steps |
|---|---|---|---|---|---|---|
| `d2-3_n3.bin:0` (figure eight) | Existence, span 2, non-rigidity | 8 / 4 | `1e-15` | `1.02e-8` | `8.23e-11` | 59, 1 rejected |
| `d11_n12_g2.bin:13` | Relative existence, span 11, non-rigidity | 257 / 7 | `1e-15` | `2.03e-9` | `1.99e-12` | 7, 1 rejected |

The eleven-dimensional run's proof phase took 424.1 seconds with four threads. Its rates are exactly
`1,2,3,4,5` in the certified axes, so inertial `2π` periodicity also follows. A shared inertial curve is
not claimed. Its outward decimal bounds were

```text
energy: [−1.448028772531032234e+01, −1.448028772530817353e+01]
action: [ 2.274558276984544620e+01,  2.274558276985254506e+01].
```

Reproduction commands:

```sh
./hyperchoreography prove catalog/d2-3_n3.bin --id 0 --digits 30 --threads 2 --force
./hyperchoreography prove catalog/d11_n12_g2.bin --id 13 --digits 30 --threads 4 --force
```

These are checked computations with the revised implementation, not independently replayable proof
packages. A self-contained witness still needs the high-precision center, exact frame, Q, dropped equations,
preconditioners, enclosures, arithmetic/compiler assumptions, and version/hash metadata. The stored radius
alone cannot reconstruct the theorem. `--write` now marks proof revision 2; legacy markers are not trusted.

## 6. Efficiency and memory

Tests ran on an Apple M5 Pro, 18 logical CPUs, 48 GiB RAM, Apple clang 21, Accelerate, GMP 6.3.0 and MPFR
4.2.2. The installed GMP is an armv8-a build; MPFR reports thread-local storage and default tuning.
Compiler vectorization diagnostics confirmed two-double vectors in the new double LU kernels. This is
evidence for this ARM build, not a measured AVX-512 result or proof of optimal instruction selection.

The double action kernels already use structure-of-arrays, contiguous sample loops and SIMD directives.
Retain those strengths. Do not introduce fast-math into proof translation units: reassociation, finite-only
assumptions or flushing subnormals can invalidate explicit bounds. The interval header now rejects the
compiler's fast-math/finite-only flags. MPFR's directed rounding is specified in its
[primary manual](https://www.mpfr.org/mpfr-current/mpfr.html).

### Measured kernel correction

`inverse_d` formerly performed n independent LU factorizations. A deterministic n=257 test went from
226.8 ms to 8.35 ms in one before/after pair, about 27×, with the same maximum inverse residual
`3.685e-15`. Another post-change run was 15.1 ms: these are short, scheduling-sensitive measurements,
not a controlled throughput study. `make audit-bench` reproduces the matrix and checks the residual.

In the whole-trial benchmark at N=12, d=11, K=24 and `g2:1,2,4,5`, the observed means were approximately
750 ms before and 726 ms after. That difference is too small and noisy to claim a general search speedup.
The roughly 9 ms Hessian eigensolve remains much more expensive than the roughly 1 ms Hessian assembly
in that kernel benchmark. Proof time is dominated by validated flow/tangents, not just inversion.

### Allocation and parallel work

Shooting now keeps one scalar workspace per worker and reuses it across Newton iterations. LU no longer
copies/refactors for every right-hand side. Mixed binary64–interval products avoid creating an MPFR
interval for every entry of a preconditioned matrix product. The pool is persistent and exception-safe
for job execution and partial construction. Expensive metadata no longer blocks every search worker at
the catalogue mutex; publication rechecks duplicates to handle concurrent discoveries correctly.

High precision has a real memory cost: each interval owns two MPFR endpoints and two limb buffers.
Coefficient storage per Taylor object grows roughly as
`(p+1)[2Nd + N(N−1)d/2 + N(N−1)]`, before scratch and tangents; each tangent worker repeats much of this.
The full derivative and reduced Jacobian also cost O((2Nd)²) intervals. More threads can therefore lose
to allocation/cache pressure. Benchmark thread counts and precision together, including peak resident
memory, rather than defaulting every nested task to all hardware threads.

Next candidates are batched tangent columns with shared primal coefficients, blocked mixed-precision
preconditioning, better interval wrapping, and an explicit memory budget for workers. Wider x86 vectors,
custom allocators, FFT synthesis, or tuned GMP should be justified with representative workloads on that
machine. This audit did not benchmark an x86 CPU or rebuild system dependencies.

## 7. The newly supplied Li–Liao paper changes the reference set

The paper reports a large set of general spatial periodic orbits, including 21 equal-mass choreographies;
its unequal-mass orbits and reported Floquet stability are different categories from this engine's
equal-mass single-cycle shooting and interval existence tests.
[Li & Liao, 2025, v1](https://arxiv.org/html/2508.08568v1).

The practical response is implemented, not just bibliographic: a pinned 21-row source table, a bounded
reference importer, record-level arXiv citations preserved by merge, and gallery attribution. See
[data/README.md](data/README.md) and [CREDITS.md](CREDITS.md).

All 21 published states passed the independent `T/3` cyclic shift check after center-of-mass translation,
period scaling and body-order selection. Eighteen also met the Fourier seed gate; IDs 1062, 1414 and 1488
remain flagged as under-resolved by this representation. None of the admitted seeds was identified with
an existing record in `d2-3_n3.bin`. Their normalized virial actions range from approximately 23.70 to
90.20, versus 16.98, 17.30 and 18.43 for the existing three spatial records. This comparison uses the
published periods and current numerical equivalence test; it is not an exhaustive minimal-period,
cover, or literature classification, and does not establish novelty of our records.

No source orbit is relabelled as this project's discovery. Source S/U labels remain source metadata;
they are not proof flags, and this audit did not re-prove their Floquet stability. The seed workflow
already benefits from independently published states without reproducing a supercomputer-scale grid.

The admitted references require K≈733–2048 in the current extraction, despite having only 18 state
variables. This is a concrete representation bottleneck: making state-space continuation/branch search
first-class can avoid a thousands-dimensional Fourier Hessian, while retaining Fourier output for plotting
and comparison. Simply truncating these seeds to the default search K does not faithfully preserve them.

A two-thread attempt at `prove catalog/li-liao-2025-n3.bin --id 62 --digits 30 --threads 2 --force --write`
was stopped after 809.25 seconds without a certificate. Sampling at about ten minutes showed active
interval tangent/MPFR multiplication work and approximately 11 MB resident memory, not a deadlock or
memory-capacity limit. No proof flag was written. This is an inconclusive bounded experiment, not a
failed existence theorem or a stability classification. A four-trial `kick` smoke test at K=64–128 also
completed without a new candidate (two Newton failures, two relative equilibria); it tests the seed-loading
workflow, not its discovery yield. The importer and all eighteen stored-state checks remain successful.

## 8. Remaining risks and highest-value next work

1. **Replayable certificates and revalidation.** Serialize exact witness data, build a small checker,
   re-prove the catalogue progressively, and distinguish existence/span/non-rigidity/inertial periodicity/
   common-curve evidence in the schema. Test interval primitives at multiple precisions and adversarial
   values, not only random ordinary-size operands. Numerical agreement to many digits, including CNS
   agreement, is not a replacement for validated inequalities.
2. **Conditioning before brute force.** Introduce analytic tangent shooting, augmented QR/SVD LM and
   multiple shooting. Then evaluate Lohner-type QR/doubleton representations for validated flow, with
   independent regression oracles. The shared recurrence architecture is worth preserving.
3. **Bifurcation-aware discovery.** Project out gauge directions; follow small transverse eigenvalues and
   singular subspaces, not raw eigenvectors selected by magnitude. The current `kick` chooses the first
   eigenvectors in ascending eigenvalue order, i.e. most negative, despite calling them “soft”. The
   continuation switch is also heuristic near degenerate pairs. Continuation success must not mean only
   small action gradient; the corrected arclength/gauge check is a first repair, not a full bifurcation solver.
4. **Search beyond one frame ray.** Vary independent frame rates, not only `Ω=sΩ0`. Catalog spectra modulo
   permutation/sign and account for logarithm shifts modulo N. Separate an integer/rational closing-frame
   search from a general relative search. Explore non-calibrated spectra too. For a target high dimension,
   deliberately excite transverse modes near an actual branch crossing; do not interpret failed random
   starts or the vertical-pattern budget as impossibility.
5. **Stability and distinctness.** At a relative solution use `M=(GS)^(-1)DΦ_τ`. Equivariance gives
   `DΦ_(Nτ)=G^N M^N`; when `G^N=I`, full-period multipliers are those of `M^N`. Quotient neutral directions,
   check symplectic defects, resolve reciprocal/conjugate pairs and near-unit-circle conditioning. Unit
   modulus alone does not exclude non-semisimple growth. Finite Fourier Morse data do not establish
   Floquet stability or a rigorous infinite-dimensional Morse index. Close-orbit identity/distinctness
   ultimately needs compatible proof boxes or sharper independently checked invariants.
6. **Runtime/persistence limits.** The issued-trial checkpoint can skip in-flight work after a crash; a
   durable completion frontier or replay queue is needed. `SearchState::save` still lacks robust error
   reporting/fsync/transaction pairing with the catalogue. Time limits and signals do not cooperatively
   cancel a long inner integration. In the measured d=13, N=14, 128-trial run with a two-minute launch
   budget, one trial extended the run to 614.8 seconds. No new candidate was found. That is a scheduling
   and sampling observation, not evidence of nonexistence.
7. **Other explicit limits.** `Ctx` caches by K and assumes a fixed configuration; arbitrary reuse across
   N/d/frame changes is unsafe. `inertia_gauge` adds the gauge count to nullity without removing lifted
   positive eigenvalues, so its total inertia counts are inconsistent; its max-entry residual tolerance
   is not an operator-norm theorem. CLI sizes and raw binary payloads need more exhaustive overflow,
   semantic and fuzz tests. Calibration storage has d≤16 and fixed-size minor limits; “arbitrary dimension”
   does not describe every ancillary routine. Proof currently supports α=1, and double-valued tolerance/
   step controls limit the “any digits” claim. These are not marked solved by the new regressions.

The near-term scientific target should be a few genuinely distinct, attributed, fully replay-certified
high-dimensional or unusual orbits, together with a measured improvement in difficult return-map
conditioning. Raw record count, a large calibration score, or a fast unconstrained trial rate is not an
adequate substitute.

## 9. Validation and reproducibility

- Final `make test`: original kernel/physics tests and all 25 adversarial audit regressions passed.
- Portable builds without MPFR/Accelerate, and AddressSanitizer/UndefinedBehaviorSanitizer runs of the
  adversarial tests, were exercised during the audit.
- The ordinary checks of all five original `d2-3_n3.bin` records passed the shift gate. Some full-window
  errors are larger, consistent with why a short shooting segment is used.
- The two proof commands above rechecked existence, span and non-rigidity with corrected bounds.
- The attributed import independently checked the published states; all eighteen admitted records passed
  `verify` with maximum state residual `3.90e-13`. Merging the reference file with itself retained eighteen
  records and their citations. The bounded reference-proof attempt makes no new proof/stability claim.
- `make gallery` regenerated `docs/index.html` from the existing generator, with 1,818 records and eighteen
  attributed references. Static checks parsed its JSON and JavaScript and verified record keys, attribution,
  finite metadata, and packed-curve lengths. `make gallery-check` now repeats these packaging checks and
  compares every record's dimensions, resolution, source citations and proof status with its catalogue.
  These are static checks, not browser/visual QA. Six old records triggered
  sampled-PCA/stored-dimension disagreements; these are reported diagnostics, not silently changed
  scientific dimensions. Geometry remains a numerical rendering.
- `git diff --check` and final builds/tests passed; rerun them after further changes. Publishing follows
  the existing GitHub Pages configuration, `main:/docs`. No original catalogue or sibling worktree was replaced.
