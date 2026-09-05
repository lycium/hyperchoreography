# Attributed reference orbits

[li-liao-2025-choreographies.tsv](li-liao-2025-choreographies.tsv) contains the 21 equal-mass spatial
choreographies identified in Table S2 of Li & Liao (2025). See [CREDITS.md](../CREDITS.md) for the source,
snapshot identifier, and attribution. Values are numerical approximations, not interval enclosures.

Build and reproduce the comparison/import with:

```sh
make hyperchoreography_reference
./hyperchoreography_reference data/li-liao-2025-choreographies.tsv catalog/d2-3_n3.bin reference-seeds.bin
./hyperchoreography verify reference-seeds.bin
```

The third path must not already exist. Omit it to run the comparison without writing a catalogue. The
importer never overwrites or merges the original catalogue and never sets a proof or stability flag.
It accepts this equal-mass table only; it is not a general unequal-mass importer.

For original period `T`, normalization uses

```text
lambda = (2π/T)^(2/3)
q_new = lambda * (q_old − center_of_mass)
v_new = lambda^(−1/2) * v_old
```

Both orientations of the cyclic body order are checked. An initial relative shift residual below `1e-4`
is required; numerical shooting must finish below `1e-10`, changing the normalized state by less than
`1e-5` relative to its scale. A retained Fourier representation must have residual at most `1e-6`,
numerical rank three, and non-rigidity above `1e-4`. These are engineering gates, not existence proofs or
proofs that Newton stayed on the identical mathematical orbit. The tool reports the state change explicitly.

On the audit machine all 21 source states passed the direct shift check (maximum `6.63e-13` before
canonicalization). Eighteen passed all seed gates. IDs 1062, 1414 and 1488 failed the Fourier gate, with
coefficient residuals approximately `5.63e-6`, `3.22e-5`, and `1.62e-6`; their source states are retained in
the table. This is a representation issue, not evidence against the published orbits. None of the eighteen
admitted entries matched `catalog/d2-3_n3.bin` under the engine's current numerical comparison.

The branch includes the generated, attributed subset as `catalog/li-liao-2025-n3.bin`. It can seed the
existing search without re-running a billion-point grid:

```sh
./hyperchoreography search --N 3 --d 3 --K 64 --Kmax 256 --starts kick \
  --seed-from catalog/li-liao-2025-n3.bin --threads 4 --trials 128 --out li-liao-descendants.bin
```

This is an example workflow, not a completed harvest or a prediction of new solutions. Long-period
reference loops may require substantially more modes than this initial search resolution. Higher-dimensional
embedding and branch switching need transverse/gauge-aware continuation; an arbitrary kick is not guaranteed
to leave an invariant lower-dimensional family. The current `--minutes` setting does not interrupt an
already-running trial.

For rigorous existence use `prove` separately. Neither the source's `S/U` classification nor successful
import is a proof of stability. Do not describe these imported references as new discoveries of this project.
The audit's two-thread, 30-digit proof attempt for source ID 62 was stopped after 809 seconds without a
certificate; no proof flag was written. This is inconclusive, not evidence against the source orbit.

## Binary citation suffix

After the layout ≥ 1 frame and state, an optional sequence of binary64 values is
`[0x53524331, count, yymm, arxiv_number, version, orbit_id, ...]`. All fields are small exact integers;
each citation occupies four values. Older readers can ignore the suffix. The current reader retains it
through serialization, proof-marker updates and duplicate replacement. Author names and scientific context
belong to the bibliographic record, not the symmetry-generator string.
