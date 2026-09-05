// Invariants, canonical form (cover unwinding, principal axes) and the Procrustes equivalence test.
#pragma once
#include "action.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>

struct CurveStats {
  double action = 0, kinE = 0, potE = 0, energy = 0, energy_std = 0, minsep = 0, rms = 0, maxr = 0, Lnorm = 0;
  std::vector<double> Lsv;   // singular values of the angular momentum 2-form
};

inline CurveStats curve_stats(const Problem& P, const double* x, Work& w) {
  CurveStats S; const int M = P.M, d = P.d, N = P.N, MN = M / N;
  S.action = action_grad(P, x, nullptr, w, &S.kinE, &S.potE);
  std::vector<double> xd; deriv_coeffs(P, x, xd);
  synth(P, x, w.Q.data()); synth(P, xd.data(), w.W.data());
  const double* Q = w.Q.data(); const double* V = w.W.data();
  auto q = [&](int j, int a) { return Q[(size_t)a * 2 * M + j]; };
  auto v = [&](int j, int a) { return V[(size_t)a * 2 * M + j]; };
  double Esum = 0, E2sum = 0, mn2 = INF;
  for (int j = 0; j < MN; j++) {
    double KE = 0, PE = 0;
    for (int i = 0; i < N; i++) for (int a = 0; a < d; a++) KE += 0.5 * v(j + i * MN, a) * v(j + i * MN, a);
    for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++) { double r2 = 0; for (int a = 0; a < d; a++) { double df = q(j + i * MN, a) - q(j + l * MN, a); r2 += df * df; }
      mn2 = std::min(mn2, r2); PE -= P.alpha == 1.0 ? 1.0 / std::sqrt(r2) : std::pow(r2, -0.5 * P.alpha); }
    double E = KE + PE; Esum += E; E2sum += E * E;
  }
  S.energy = Esum / MN; S.energy_std = std::sqrt(std::max(0.0, E2sum / MN - S.energy * S.energy));
  S.minsep = std::sqrt(mn2);
  double r2sum = 0; for (int i = 0; i < P.nb; i++) for (int a = 0; a < d; a++) r2sum += x[i * d + a] * x[i * d + a];
  S.rms = std::sqrt(0.5 * r2sum);
  for (int j = 0; j < M; j++) { double r2 = 0; for (int a = 0; a < d; a++) r2 += q(j, a) * q(j, a); S.maxr = std::max(S.maxr, std::sqrt(r2)); }
  std::vector<double> L((size_t)d * d, 0.0);
  for (int i = 0; i < N; i++) for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) L[(size_t)a * d + b] += q(i * MN, a) * v(i * MN, b) - q(i * MN, b) * v(i * MN, a);
  for (int a = 0; a < d; a++) for (int b = a + 1; b < d; b++) S.Lnorm += L[(size_t)a * d + b] * L[(size_t)a * d + b];
  S.Lnorm = std::sqrt(S.Lnorm);
  S.Lsv = la::singular_values(d, L);
  return S;
}

// Rigidity defect: a relative equilibrium moves rigidly, so every mutual distance is constant. The
// choreography shift sends the pair (j,k) to (j−k,0), so ⌊N/2⌋ separations exhaust the test, and a common
// rotation cancels, so the answer does not depend on Ω. Zero ⇔ rigid — the high-d analogue of the N-gon.
inline double rigid_defect(const Problem& P, const double* x, Work& w) {
  const int M = P.M, d = P.d, MN = M / P.N; double worst = 0;
  synth(P, x, w.Q.data()); const double* Q = w.Q.data();
  for (int p = 1; p <= P.N / 2; p++) {
    double lo = INF, hi = 0;
    for (int j = 0; j < M; j++) { double r2 = 0;
      for (int a = 0; a < d; a++) { double e = Q[(size_t)a * 2 * M + j + p * MN] - Q[(size_t)a * 2 * M + j]; r2 += e * e; }
      lo = std::min(lo, r2); hi = std::max(hi, r2); }
    if (hi > 0) worst = std::max(worst, 1.0 - std::sqrt(lo / hi));
  }
  return worst;
}

// The same test on a record's own (modes, coef): the final series after shooting, not the pre-certification iterate.
inline double rigid_defect_coef(int N, int d, const std::vector<int>& modes, const double* coef) {
  const int mmax = modes.empty() ? 1 : modes.back();                 // ascending
  const int M = N * std::max(8, 4 * mmax / N + 2), MN = M / N;       // above Nyquist
  std::vector<double> Q((size_t)M * d, 0.0);
  for (int j = 0; j < M; j++) { const double t = 2 * PI * j / M; double* q = &Q[(size_t)j * d];
    for (size_t mu = 0; mu < modes.size(); mu++) { const double c = std::cos(modes[mu] * t), sn = std::sin(modes[mu] * t);
      for (int a = 0; a < d; a++) q[a] += c * coef[mu * 2 * d + a] + sn * coef[mu * 2 * d + d + a]; } }
  double worst = 0;
  for (int p = 1; p <= N / 2; p++) { double lo = INF, hi = 0;
    for (int j = 0; j < M; j++) { double r2 = 0; const double* a0 = &Q[(size_t)j * d]; const double* a1 = &Q[(size_t)((j + p * MN) % M) * d];
      for (int a = 0; a < d; a++) { const double e = a1[a] - a0[a]; r2 += e * e; }
      lo = std::min(lo, r2); hi = std::max(hi, r2); }
    if (hi > 0) worst = std::max(worst, 1.0 - std::sqrt(lo / hi)); }
  return worst;
}

// per-mode power P_m = |c_m|² + |s_m|²
inline std::vector<double> mode_power(const Problem& P, const double* x) {
  std::vector<double> pw(P.nm, 0.0);
  for (int mu = 0; mu < P.nm; mu++) for (int a = 0; a < P.d; a++) { double c = x[(2 * mu) * P.d + a], s = x[(2 * mu + 1) * P.d + a]; pw[mu] += c * c + s * s; }
  return pw;
}
// gcd of the significant modes
inline int cover_multiplicity(const Problem& P, const double* x, double rel_tol = 1e-9) {
  auto pw = mode_power(P, x); double pmax = 0; for (double v : pw) pmax = std::max(pmax, v);
  int g = 0; for (int mu = 0; mu < P.nm; mu++) if (pw[mu] > rel_tol * pmax) g = la::gcd(g, P.modes[mu]);
  return std::max(g, 1);
}
// mode km → m, amplitude × k^{2/(α+2)}
inline void unwind_cover(const Problem& P, const double* x, int k, Problem& P2, std::vector<double>& x2, int Kmin = 8) {
  P2.init(P.N, P.d, std::max(P.K / k, Kmin), 0, P.alpha); P2.minsep = P.minsep;
  x2.assign(P2.n, 0.0); double lam = std::pow((double)k, 2.0 / (P.alpha + 2.0));
  for (int mu = 0; mu < P.nm; mu++) { int m = P.modes[mu]; if (m % k) continue;
    auto it = std::lower_bound(P2.modes.begin(), P2.modes.end(), m / k); if (it == P2.modes.end() || *it != m / k) continue;
    int m2 = (int)(it - P2.modes.begin());
    for (int a = 0; a < P.d; a++) { x2[(2 * m2) * P.d + a] = lam * x[(2 * mu) * P.d + a]; x2[(2 * m2 + 1) * P.d + a] = lam * x[(2 * mu + 1) * P.d + a]; } }
}
// principal axes with fixed column signs; returns the effective dimension
inline int canonical_frame(int nb, int d, std::vector<double>& x, std::vector<double>& sv, double rel_tol = 1e-8, std::vector<double>* Rout = nullptr) {
  std::vector<double> S((size_t)d * d, 0.0), w;
  for (int i = 0; i < nb; i++) for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) S[(size_t)a * d + b] += x[i * d + a] * x[i * d + b];
  la::sym_eig(d, S, w);
  std::vector<double> y((size_t)nb * d, 0.0);
  for (int i = 0; i < nb; i++) for (int c = 0; c < d; c++) { int k = d - 1 - c; double s = 0; for (int a = 0; a < d; a++) s += x[i * d + a] * S[(size_t)a * d + k]; y[i * d + c] = s; }
  std::vector<double> sg(d, 1.0);
  for (int c = 0; c < d; c++) { int ib = 0; double best = -1; for (int i = 0; i < nb; i++) if (std::fabs(y[i * d + c]) > best) { best = std::fabs(y[i * d + c]); ib = i; }
    if (y[ib * d + c] < 0) { sg[c] = -1; for (int i = 0; i < nb; i++) y[i * d + c] = -y[i * d + c]; } }
  if (Rout) { Rout->assign((size_t)d * d, 0.0);                 // y = R x, so any Ω must become R Ω Rᵀ
    for (int c = 0; c < d; c++) for (int a = 0; a < d; a++) (*Rout)[(size_t)c * d + a] = sg[c] * S[(size_t)a * d + (d - 1 - c)]; }
  x.swap(y); sv.resize(d); int deff = 0;
  for (int c = 0; c < d; c++) { sv[c] = std::sqrt(std::max(0.0, w[d - 1 - c])); if (w[d - 1 - c] > rel_tol * w[d - 1]) deff++; }
  return deff;
}

// loop at times eps·t + tau, zero-padded to D dimensions
inline void sample_curve(const std::vector<int>& modes, int d, const double* x, int Mc, double eps, double tau, int D, std::vector<double>& out) {
  out.assign((size_t)Mc * D, 0.0); const int nm = (int)modes.size();
  for (int j = 0; j < Mc; j++) { double t = eps * 2.0 * PI * j / Mc + tau;
    for (int mu = 0; mu < nm; mu++) { double c = std::cos(modes[mu] * t), s = std::sin(modes[mu] * t);
      for (int a = 0; a < d; a++) out[(size_t)j * D + a] += c * x[(2 * mu) * d + a] + s * x[(2 * mu + 1) * d + a]; } }
}
// deff of the physical motion. A rotating frame can raise the loop's rank (a fixed point sweeps a circle) or
// lower it (a circularly polarised mode at rate −w becomes a linear oscillation), so the loop's own rank is
// not the dimension the orbit occupies. Body k traces exp(Ωt)q(t + 2πk/N) = exp(−2πΩk/N)·(body 0's path), so
// the occupied subspace is the sum of N rotated copies. Invariant under a common rotation of q and Ω.
inline int inertial_deff(int N, const std::vector<int>& modes, int d, const double* x, const double* om, std::vector<double>& sv, double rel_tol = 1e-8) {
  int mx = 0; for (int m : modes) mx = std::max(mx, m);
  double wmx = 0; for (int i = 0; i < d * d; i++) wmx = std::max(wmx, std::fabs(om[i]));
  const int Mc = std::max(64, 4 * (mx + (int)std::ceil(wmx) + 1));
  std::vector<double> q, S((size_t)d * d, 0.0), Ai, R, u(d), w;
  sample_curve(modes, d, x, Mc, 1.0, 0.0, d, q);
  for (int j = 0; j < Mc; j++) {
    Ai.assign(om, om + (size_t)d * d); const double t = 2 * PI * j / Mc; for (double& e : Ai) e *= t;
    la::expm_skew(d, Ai, R);
    for (int a = 0; a < d; a++) { double t2 = 0; for (int b = 0; b < d; b++) t2 += R[(size_t)a * d + b] * q[(size_t)j * d + b]; u[a] = t2; }
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) S[(size_t)a * d + b] += u[a] * u[b];
  }
  std::vector<double> G((size_t)d * d, 0.0), Rk, R2;
  for (int k = 0; k < N; k++) {                                     // Σ_k exp(−2πΩk/N) S exp(+2πΩk/N)
    Ai.assign(om, om + (size_t)d * d); const double t = -2 * PI * k / N; for (double& e : Ai) e *= t;
    la::expm_skew(d, Ai, Rk); R2.assign((size_t)d * d, 0.0);
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double t2 = 0;
      for (int i = 0; i < d; i++) t2 += Rk[(size_t)a * d + i] * S[(size_t)i * d + b]; R2[(size_t)a * d + b] = t2; }
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double t2 = 0;
      for (int i = 0; i < d; i++) t2 += R2[(size_t)a * d + i] * Rk[(size_t)b * d + i]; G[(size_t)a * d + b] += t2; }
  }
  la::sym_eig(d, G, w); sv.resize(d); int deff = 0;
  for (int c = 0; c < d; c++) { sv[c] = std::sqrt(std::max(0.0, w[d - 1 - c])); if (w[d - 1 - c] > rel_tol * w[d - 1]) deff++; }
  return deff;
}

// ||A||² + ||B||² − 2 ||AᵀB||_*
inline double procrustes_res2(int Mc, int D, const std::vector<double>& A, const std::vector<double>& B, double nA, double nB) {
  std::vector<double> C((size_t)D * D, 0.0);
  for (int j = 0; j < Mc; j++) for (int a = 0; a < D; a++) { double qa = A[(size_t)j * D + a]; if (qa == 0) continue; for (int b = 0; b < D; b++) C[(size_t)a * D + b] += qa * B[(size_t)j * D + b]; }
  auto sv = la::singular_values(D, C); double nuc = 0; for (double s : sv) nuc += s;
  return nA + nB - 2 * nuc;
}
// relative distance modulo time shift, reversal and O(d)
struct LoopAlignment { int eps = 1; double tau = 0; };
inline double loop_distance(int N, const std::vector<int>& mA_, int dA, const double* xa, const std::vector<int>& mB_, int dB, const double* xb, int max_modes = 64, LoopAlignment* alignment = nullptr) {
  const int D = std::max(dA, dB); const int Mc = N * std::max(1, (256 + N - 1) / N);
  std::vector<int> mA(mA_.begin(), mA_.begin() + std::min<size_t>(mA_.size(), max_modes)), mB(mB_.begin(), mB_.begin() + std::min<size_t>(mB_.size(), max_modes));
  std::vector<double> A, B, Bs; sample_curve(mA, dA, xa, Mc, 1, 0, D, A); sample_curve(mB, dB, xb, Mc, 1, 0, D, B);
  double nA = 0, nB = 0; for (double v : A) nA += v * v; for (double v : B) nB += v * v;
  double best = INF; int beps = 1, bj = 0;
  for (int eps = 1; eps >= -1; eps -= 2) for (int j0 = 0; j0 < Mc; j0++) {
    Bs.assign((size_t)Mc * D, 0.0);
    for (int j = 0; j < Mc; j++) { int jj = ((eps * j + j0) % Mc + Mc) % Mc; for (int a = 0; a < D; a++) Bs[(size_t)j * D + a] = B[(size_t)jj * D + a]; }
    double r2 = procrustes_res2(Mc, D, A, Bs, nA, nB);
    if (r2 < best) { best = r2; beps = eps; bj = j0; }
  }
  // golden-section refinement of the shift
  auto res = [&](double tau) { sample_curve(mB, dB, xb, Mc, beps, tau, D, Bs); return procrustes_res2(Mc, D, A, Bs, nA, nB); };
  double h = 2.0 * PI / Mc, lo = 2.0 * PI * bj / Mc - h, hi = 2.0 * PI * bj / Mc + h, gr = 0.6180339887498949;
  double c = hi - gr * (hi - lo), dd = lo + gr * (hi - lo), fc = res(c), fd = res(dd);
  for (int it = 0; it < 40; it++) { if (fc < fd) { hi = dd; dd = c; fd = fc; c = hi - gr * (hi - lo); fc = res(c); } else { lo = c; c = dd; fc = fd; dd = lo + gr * (hi - lo); fd = res(dd); } }
  double tau = 2.0 * PI * bj / Mc;
  if (fc < best) { best = fc; tau = c; }
  if (fd < best) { best = fd; tau = dd; }
  if (alignment) { alignment->eps = beps; alignment->tau = tau; }
  return std::sqrt(std::max(0.0, best) / (0.5 * (nA + nB)));
}
inline double loop_distance(const Problem& PA, const double* xa, const Problem& PB, const double* xb) { return loop_distance(PA.N, PA.modes, PA.d, xa, PB.modes, PB.d, xb); }

// Symmetries of the loop itself: all (ε, θ, R) ∈ O(2) × O(d) with q(εt + θ) = R q(t). For a fixed (ε, θ) the
// best R is the Procrustes fit, so the group is exactly the set of zeros of the self-distance residual — the
// same function loop_distance minimises over. In the canonical frame R commutes with the inertia tensor, so
// it is diag(±1) wherever the principal values are simple, which is what the --sym DSL expresses.
struct LoopSym { int eps = 1, p = 0, q = 1; double theta = 0, res = 0; std::string text; };

// polar factor of C: the orthogonal Q maximising ⟨Q, C⟩. Undetermined off the loop's span, completed there
// by Gram–Schmidt so the answer is deterministic.
inline void procrustes_rot(int d, const std::vector<double>& C, std::vector<double>& Q) {
  std::vector<double> G((size_t)d * d, 0.0), w, U((size_t)d * d, 0.0), V((size_t)d * d), u(d);
  for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) { double s = 0;
    for (int k = 0; k < d; k++) s += C[(size_t)k * d + i] * C[(size_t)k * d + j]; G[(size_t)i * d + j] = s; }
  la::sym_eig(d, G, w);
  for (int a = 0; a < d; a++) for (int c = 0; c < d; c++) V[(size_t)a * d + c] = G[(size_t)a * d + (d - 1 - c)];
  const double top = std::sqrt(std::max(0.0, w[d - 1]));
  int nu = 0;
  for (; nu < d; nu++) { double sv = std::sqrt(std::max(0.0, w[d - 1 - nu]));
    if (sv <= 1e-8 * std::max(top, 1e-300)) break;                  // descending, so the rest are null too
    for (int a = 0; a < d; a++) { double s = 0;
      for (int b = 0; b < d; b++) s += C[(size_t)a * d + b] * V[(size_t)b * d + nu]; U[(size_t)a * d + nu] = s / sv; } }
  for (int e = 0; e < d && nu < d; e++) {                           // off the loop's span Q is arbitrary
    std::fill(u.begin(), u.end(), 0.0); u[e] = 1;
    for (int p = 0; p < nu; p++) { double dp = 0;
      for (int a = 0; a < d; a++) dp += U[(size_t)a * d + p] * u[a];
      for (int a = 0; a < d; a++) u[a] -= dp * U[(size_t)a * d + p]; }
    double n2 = 0; for (double v : u) n2 += v * v;
    if (n2 < 1e-12) continue;
    n2 = std::sqrt(n2);
    for (int a = 0; a < d; a++) U[(size_t)a * d + nu] = u[a] / n2;
    nu++; }
  Q.assign((size_t)d * d, 0.0);
  for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double s = 0;
    for (int c = 0; c < d; c++) s += U[(size_t)a * d + c] * V[(size_t)b * d + c]; Q[(size_t)a * d + b] = s; }
}

// Validate a proposed alignment on ALL Fourier coefficients, using an explicit squared residual
// instead of subtracting nearly equal norms. The same rotation must intertwine the two frames.
// A failed fit may retain a duplicate; scalar invariants must never delete a distinct orbit.
inline double aligned_loop_distance(const std::vector<int>& ma, int da, const double* xa, const double* oa,
                                    const std::vector<int>& mb, int db, const double* xb, const double* ob,
                                    const LoopAlignment& fit, double frame_tol = 1e-7) {
  const int d = std::max(da, db);
  std::vector<int> modes; std::set_union(ma.begin(), ma.end(), mb.begin(), mb.end(), std::back_inserter(modes));
  const size_t rows = 2 * modes.size();
  std::vector<double> a(rows * d, 0), b(rows * d, 0), cov((size_t)d * d, 0), rot;
  for (size_t k = 0; k < modes.size(); k++) {
    auto ia = std::lower_bound(ma.begin(), ma.end(), modes[k]), ib = std::lower_bound(mb.begin(), mb.end(), modes[k]);
    if (ia != ma.end() && *ia == modes[k]) for (int h = 0; h < 2; h++) for (int c = 0; c < da; c++)
      a[(2 * k + h) * d + c] = xa[(2 * (ia - ma.begin()) + h) * da + c];
    if (ib != mb.end() && *ib == modes[k]) {
      const double co = std::cos(modes[k] * fit.tau), si = std::sin(modes[k] * fit.tau);
      for (int c = 0; c < db; c++) {
        const double cm = xb[2 * (ib - mb.begin()) * db + c], sm = xb[(2 * (ib - mb.begin()) + 1) * db + c];
        b[2 * k * d + c] = co * cm + si * sm;
        b[(2 * k + 1) * d + c] = fit.eps * (-si * cm + co * sm);
      }
    }
  }
  for (size_t k = 0; k < rows; k++) for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) cov[(size_t)i * d + j] += b[k * d + i] * a[k * d + j];
  procrustes_rot(d, cov, rot);
  double orth = 0, err = 0, scale = 0;
  for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) { double v = -(i == j);
    for (int k = 0; k < d; k++) v += rot[(size_t)k * d + i] * rot[(size_t)k * d + j]; orth = std::max(orth, std::fabs(v)); }
  if (!(orth < 1e-6)) return INF;
  for (size_t k = 0; k < rows; k++) for (int i = 0; i < d; i++) {
    double v = -b[k * d + i]; for (int j = 0; j < d; j++) v += rot[(size_t)i * d + j] * a[k * d + j];
    err += v * v; scale += 0.5 * (a[k * d + i] * a[k * d + i] + b[k * d + i] * b[k * d + i]);
  }
  if (oa || ob) {
    double mismatch = 0, norm = 1;
    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) { double v = 0;
      for (int k = 0; k < d; k++) {
        if (ob && i < db && k < db) v += ob[(size_t)i * db + k] * rot[(size_t)k * d + j];
        if (oa && k < da && j < da) v -= fit.eps * rot[(size_t)i * d + k] * oa[(size_t)k * da + j];
      }
      mismatch = std::max(mismatch, std::fabs(v));
      if (oa && i < da && j < da) norm = std::max(norm, std::fabs(oa[(size_t)i * da + j]));
      if (ob && i < db && j < db) norm = std::max(norm, std::fabs(ob[(size_t)i * db + j]));
    }
    if (!(mismatch <= frame_tol * norm)) return INF;
  }
  return scale > 0 ? std::sqrt(err / scale) : INF;
}
inline void best_frac(double v, int qmax, int& p, int& q) {
  double best = INF; p = 0; q = 1;
  for (int qq = 1; qq <= qmax; qq++) { int pp = (int)std::lround(v * qq);
    double e = std::fabs(v - (double)pp / qq); if (e < best - 1e-14) { best = e; p = pp; q = qq; } }
}
// All symmetries of one loop. Work in the deff-dimensional span (the canonical frame puts it first): off the
// span R is undetermined, and fitting it there turns an exact signed permutation into a generic rotation.
// cont = true means every shift works, i.e. the loop is a relative equilibrium with a continuous S¹ symmetry.
inline std::vector<LoopSym> detect_symmetry(int N, const std::vector<int>& modes_, int d, int deff, const double* x,
                                            int* cont = nullptr, double tol = 1e-6, int max_modes = 256) {
  const int Mc = N * std::max(1, (512 + N - 1) / N), r = std::max(1, std::min(deff, d));
  std::vector<int> md(modes_.begin(), modes_.begin() + std::min<size_t>(modes_.size(), max_modes));
  std::vector<double> xr(md.size() * 2 * r, 0.0);                          // restrict to the span
  for (size_t mu = 0; mu < md.size(); mu++) for (int h = 0; h < 2; h++) for (int a = 0; a < r; a++)
    xr[(2 * mu + h) * r + a] = x[(2 * mu + h) * d + a];
  std::vector<double> A, Bs, C, Q; sample_curve(md, r, xr.data(), Mc, 1, 0, r, A);
  double nA = 0; for (double v : A) nA += v * v;
  if (cont) *cont = 0;
  if (nA <= 0) return {};
  auto res = [&](int eps, double tau) { sample_curve(md, r, xr.data(), Mc, eps, tau, r, Bs); return procrustes_res2(Mc, r, A, Bs, nA, nA) / nA; };
  const double tol2 = tol * tol;                                           // nA + nB − 2‖C‖_* cancels, so the
  // relative squared residual bottoms out near 1e-15 however exact the loop is
  std::vector<LoopSym> out;
  for (int eps = 1; eps >= -1; eps -= 2) {
    std::vector<double> f(Mc); int below = 0;
    for (int j = 0; j < Mc; j++) { f[j] = res(eps, 2.0 * PI * j / Mc); if (f[j] < tol2) below++; }
    if (below > Mc / 2) { if (cont) *cont |= eps > 0 ? 1 : 2; continue; }      // continuous: every shift works
    for (int j = 0; j < Mc; j++) {
      if (f[j] > 1e-3 || f[j] > f[(j + Mc - 1) % Mc] || f[j] > f[(j + 1) % Mc]) continue;
      double h = 2.0 * PI / Mc, lo = 2.0 * PI * j / Mc - h, hi = lo + 2 * h, gr = 0.6180339887498949;
      double c = hi - gr * (hi - lo), e = lo + gr * (hi - lo), fc = res(eps, c), fe = res(eps, e);
      for (int it = 0; it < 60; it++) { if (fc < fe) { hi = e; e = c; fe = fc; c = hi - gr * (hi - lo); fc = res(eps, c); }
        else { lo = c; c = e; fc = fe; e = lo + gr * (hi - lo); fe = res(eps, e); } }
      double th = fc < fe ? c : e, rr = std::min(fc, fe);
      if (rr > tol2) continue;
      LoopSym S; S.eps = eps; S.res = std::sqrt(std::max(0.0, rr)); S.theta = th - 2 * PI * std::floor(th / (2 * PI));
      best_frac(S.theta / (2 * PI), 4 * N + 8, S.p, S.q);
      if (S.p == S.q) { S.p = 0; S.q = 1; }
      bool dup = false;
      for (auto& o : out) if (o.eps == eps && std::fabs(o.theta - S.theta) < 1e-6) dup = true;
      if (dup) continue;
      sample_curve(md, r, xr.data(), Mc, eps, S.theta, r, Bs);
      C.assign((size_t)r * r, 0.0);
      for (int jj = 0; jj < Mc; jj++) for (int a = 0; a < r; a++) { double qa = Bs[(size_t)jj * r + a]; if (qa == 0) continue;
        for (int b = 0; b < r; b++) C[(size_t)a * r + b] += qa * A[(size_t)jj * r + b]; }
      procrustes_rot(r, C, Q);                                                // q(eps t + theta) = Q q(t)
      char tb[32]; std::snprintf(tb, sizeof tb, "t%c%d/%d", eps > 0 ? '+' : '-', S.p, S.q);
      S.text = tb; std::string perm = " s["; bool sp = true;
      for (int b = 0; b < r && sp; b++) { int hit = -1; double sg = 0;
        for (int a = 0; a < r; a++) { double v = Q[(size_t)a * r + b];
          if (std::fabs(v) > 1e-5 && std::fabs(std::fabs(v) - 1) > 1e-5) { sp = false; break; }
          if (std::fabs(std::fabs(v) - 1) <= 1e-5) { if (hit >= 0) { sp = false; break; } hit = a; sg = v; } }
        if (!sp || hit < 0) { sp = false; break; }
        perm += (b ? "," : "") + std::string(sg < 0 ? "-" : "") + std::to_string(hit + 1); }
      for (int b = r; b < d && sp; b++) perm += "," + std::to_string(b + 1);   // the span's complement is free
      S.text += sp ? perm + "]" : " R";
      out.push_back(std::move(S));
    }
  }
  return out;
}
