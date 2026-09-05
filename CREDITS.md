# Scientific credit and provenance

Discovery, numerical reproduction, rigorous existence, and novelty are different claims. A record in
this repository is not, merely by being present, a new discovery attributable to this project. A small
shooting residual is not an interval existence certificate; a numerical Morse index is not Floquet stability.

## Li & Liao's spatial three-body orbits

Xiaoming Li and Shijun Liao, *Discovery of 10,059 new three-dimensional periodic orbits of general
three-body problem* (2025), [arXiv:2508.08568v1](https://arxiv.org/html/2508.08568v1).

The authors report 21 equal-mass spatial choreographies within their larger numerical periodic-orbit
dataset. Their search combines DOP853 candidate screening, Newton correction and high-precision CNS;
they obtain stability information from variational equations and Floquet analysis. These are prior work,
not discoveries of this repository. The reported stability labels are not re-certified here.
[Paper, sections 2–3 and Table S2](https://arxiv.org/html/2508.08568v1#Sx1.T2).

The [reference table](data/li-liao-2025-choreographies.tsv) uses the Table S2 orbit IDs and the more precise
numeric initial conditions from the
[authors' dataset](https://github.com/sjtu-liao/three-body/blob/main/initial-condition-of-3D-periodic-orbits.txt).
Retrieved 2026-09-05; the complete source file's Git blob SHA-1 is
`e242b0e60df08d4992c6063594e105da16108149`. The arXiv HTML labels the paper CC BY 4.0. No external solver
code or figures have been copied into the engine.

Our transformations are center-of-mass translation, Newtonian period normalization to `2π`, body-order
selection, optional numerical shooting refinement, Fourier extraction, and orthogonal canonicalization.
The generated `catalog/li-liao-2025-n3.bin` is a separately attributed numerical reference/seed catalogue.
Its record IDs retain the authors' `O_n(1.0)` indices. It contains only entries meeting the stated seed
gates, not all 10,059 general periodic orbits. The source table retains all 21 reference entries.

Records carry citations in an optional backwards-compatible `extra[]` suffix. `show` exposes `sources`;
merge/replacement preserves their union, and the gallery links each attributed record to the paper.
Source citations mean prior identification, not rigorous identity certification. They are not automatically
assigned to existing unmatched records or to distinct descendants of a seeded search. Proper seed-ancestry
tracking remains separate work.

## Established mathematical and computational foundations

The README's existing references to Simó, Minton, Chenciner–Montgomery, Chenciner–Féjoz, and
Ferrario–Terracini remain relevant. The present audit does not establish historical priority over that
literature or perform an exhaustive literature/dataset comparison.

High-order Taylor integration has a substantial prior literature, including
[Jorba & Zou's Taylor package](https://www.maia.ub.es/~angel/taylor/taylor-1.4/taylor.pdf) and
[Biscani & Izzo's celestial-mechanics work](https://academic.oup.com/mnras/article/504/2/2614/6226683).
The validated solver also relies on the established interval Newton/Krawczyk framework. The software
dependencies include [MPFR](https://www.mpfr.org/mpfr-current/mpfr.html), GMP, and, on the tested macOS
build, Accelerate/LAPACK. These dependencies and mathematical methods are not project inventions.

See [AUDIT.md](AUDIT.md) for what was checked, what was corrected, and what remains unproved.
