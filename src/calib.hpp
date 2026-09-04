// The calibration ladder. For a proper subgroup G ⊂ SO(d) stabilising a k-form ψ, the pairing of ψ with the
// loop's jet moment  A_k = (1/2π)∮ q ∧ q̇ ∧ … ∧ q^(k−1) dt ∈ Λ^k(R^d)  is an invariant that no O(d) invariant
// of the loop can see. d = 7 (G₂, k = 3) is one rung of a tower that reaches every dimension.
#pragma once
#include "g2.hpp"
#include <bit>
#include <map>
#include <stdexcept>

// Λ^k(R^d) in the basis of increasing multi-indices
struct Wedge {
  int d = 0, k = 0, n = 0;
  std::vector<int> idx;                        // n rows of k increasing indices
  std::vector<int> at;                         // index-set bitmask → row
  void init(int d_, int k_) {
    if (d_ < 1 || d_ > 16 || k_ < 1 || k_ > d_) throw std::runtime_error("wedge: bad (d,k)");
    d = d_; k = k_; idx.clear(); at.assign((size_t)1 << d, -1);
    std::vector<int> c(k); for (int i = 0; i < k; i++) c[i] = i;
    for (;;) {
      int m = 0; for (int i = 0; i < k; i++) m |= 1 << c[i];
      at[m] = (int)(idx.size() / k); idx.insert(idx.end(), c.begin(), c.end());
      int i = k - 1; while (i >= 0 && c[i] == d - k + i) i--;
      if (i < 0) break;
      c[i]++; for (int j = i + 1; j < k; j++) c[j] = c[j - 1] + 1;
    }
    n = (int)(idx.size() / k);
  }
};
inline const Wedge& wedge_basis(int d, int k) {
  static thread_local std::map<std::pair<int, int>, Wedge> cache;
  auto key = std::make_pair(d, k); auto it = cache.find(key);
  if (it == cache.end()) { Wedge W; W.init(d, k); it = cache.emplace(key, std::move(W)).first; }
  return it->second;
}
// permutation sign of an index tuple; 0 on a repeat
inline double wedge_sort(int k, int* a) {
  double sg = 1;
  for (int i = 1; i < k; i++) for (int j = i; j > 0 && a[j] < a[j - 1]; j--) { std::swap(a[j], a[j - 1]); sg = -sg; }
  for (int i = 1; i < k; i++) if (a[i] == a[i - 1]) return 0;
  return sg;
}
inline void wedge_add(const Wedge& W, const int* I, double v, double* psi) {
  int a[16]; for (int i = 0; i < W.k; i++) a[i] = I[i];
  double sg = wedge_sort(W.k, a); if (sg == 0) return;
  int m = 0; for (int i = 0; i < W.k; i++) m |= 1 << a[i];
  psi[W.at[m]] += sg * v;
}
inline double wedge_get(const Wedge& W, const double* psi, const int* I) {
  int a[16]; for (int i = 0; i < W.k; i++) a[i] = I[i];
  double sg = wedge_sort(W.k, a); if (sg == 0) return 0;
  int m = 0; for (int i = 0; i < W.k; i++) m |= 1 << a[i];
  return sg * psi[W.at[m]];
}

inline double det_small(int k, const double* S) {
  if (k <= 0) return 1;
  if (k == 1) return S[0];
  if (k == 2) return S[0] * S[3] - S[1] * S[2];
  if (k == 3) return S[0] * (S[4] * S[8] - S[5] * S[7]) - S[1] * (S[3] * S[8] - S[5] * S[6]) + S[2] * (S[3] * S[7] - S[4] * S[6]);
  double a[64], det = 1;
  for (int i = 0; i < k * k; i++) a[i] = S[i];
  for (int c = 0; c < k; c++) {
    int p = c; for (int r = c + 1; r < k; r++) if (std::fabs(a[r * k + c]) > std::fabs(a[p * k + c])) p = r;
    if (a[p * k + c] == 0) return 0;
    if (p != c) { for (int j = 0; j < k; j++) std::swap(a[c * k + j], a[p * k + j]); det = -det; }
    det *= a[c * k + c];
    for (int r = c + 1; r < k; r++) { double f = a[r * k + c] / a[c * k + c];
      for (int j = c; j < k; j++) a[r * k + j] -= f * a[c * k + j]; }
  }
  return det;
}
// determinant plus the cofactor matrix ∂det/∂S
inline double det_cof(int k, const double* S, double* cof) {
  double m[64];
  for (int s = 0; s < k; s++) for (int t = 0; t < k; t++) {
    int p = 0;
    for (int i = 0; i < k; i++) { if (i == s) continue; for (int j = 0; j < k; j++) if (j != t) m[p++] = S[i * k + j]; }
    cof[s * k + t] = ((s + t) & 1 ? -1.0 : 1.0) * det_small(k - 1, m);
  }
  double det = 0; for (int t = 0; t < k; t++) det += S[t] * cof[t];
  return det;
}

// A_I = (1/2π) ∮ det[ q^(r)_{I_s} ]_{r,s} dt.  The integrand is a trigonometric polynomial of degree k·m_max,
// so the trapezoid rule with Ms > k·m_max samples is exact. coef is [c_0..c_{d−1}, s_0..s_{d−1}] per mode.
inline void jet_form(const std::vector<int>& modes, const double* coef, int d, const Wedge& W, double* A, int Ms = 0) {
  const int k = W.k, nm = (int)modes.size();
  int mmax = 1; for (int m : modes) mmax = std::max(mmax, m);
  if (Ms <= 0) Ms = k * mmax + 2;
  std::vector<double> D((size_t)k * d);
  for (int i = 0; i < W.n; i++) A[i] = 0.0;
  double S[64];
  for (int j = 0; j < Ms; j++) {
    const double t = 6.283185307179586 * j / Ms;
    std::fill(D.begin(), D.end(), 0.0);
    for (int mu = 0; mu < nm; mu++) {
      const double m = modes[mu], *cm = coef + (size_t)mu * 2 * d, *sm = cm + d;
      double pc = std::cos(m * t), ps = std::sin(m * t), mr = 1;
      for (int r = 0; r < k; r++) {                                  // q^(r) = m^r (pc·c + ps·s)
        double* Dr = &D[(size_t)r * d];
        for (int a = 0; a < d; a++) Dr[a] += mr * (pc * cm[a] + ps * sm[a]);
        double t2 = pc; pc = -ps; ps = t2; mr *= m;
      }
    }
    for (int i = 0; i < W.n; i++) { const int* I = &W.idx[(size_t)i * k];
      for (int r = 0; r < k; r++) for (int s = 0; s < k; s++) S[r * k + s] = D[(size_t)r * d + I[s]];
      A[i] += det_small(k, S); }
  }
  for (int i = 0; i < W.n; i++) A[i] /= Ms;
}

// Π_{r<k} ‖q^(r)‖_{L²} from the coefficients — the Hadamard bound on the k-jet wedge, so χ*/jet_scale is a
// scale-free twist comparable across records of different size.
inline double jet_scale(const std::vector<int>& modes, const double* coef, int d, int k) {
  double p = 1;
  for (int r = 0; r < k; r++) { double s = 0;
    for (size_t mu = 0; mu < modes.size(); mu++) { double e = 0;
      for (int a = 0; a < 2 * d; a++) { double v = coef[mu * 2 * d + a]; e += v * v; }
      s += 0.5 * e * std::pow((double)modes[mu], 2 * r); }
    p *= std::sqrt(s); }
  return p;
}

enum Calib { CAL_NONE = 0, CAL_SIMPLE, CAL_SL, CAL_G2, CAL_SPIN7 };
inline const char* calib_group(int id) {
  static const char* n[] = {"-", "S(OxO)", "SU(n)", "G2", "Spin(7)"};
  return n[id >= 0 && id <= CAL_SPIN7 ? id : 0];
}
inline std::vector<double> calib_psi(const Wedge& W, int id) {
  std::vector<double> psi(W.n, 0.0); const int d = W.d, k = W.k;
  if (id == CAL_SIMPLE) { std::vector<int> I(k); for (int i = 0; i < k; i++) I[i] = i; wedge_add(W, I.data(), 1.0, psi.data()); }
  else if (id == CAL_SL) {                                  // Re[(e0+ie1) ∧ (e2+ie3) ∧ … ∧ (e_{2k−2}+ie_{2k−1})]
    if (d < 2 * k) throw std::runtime_error("special-Lagrangian form needs d >= 2k");
    std::vector<int> I(k);
    for (int S = 0; S < (1 << k); S++) { int pc = (int)std::popcount((unsigned)S); if (pc & 1) continue;
      for (int p = 0; p < k; p++) I[p] = 2 * p + ((S >> p) & 1);
      wedge_add(W, I.data(), (pc / 2) & 1 ? -1.0 : 1.0, psi.data()); } }
  else if (id == CAL_G2) {                                  // associative 3-form: the seven Fano lines
    if (k != 3 || d < 7) throw std::runtime_error("associative form needs k = 3, d >= 7");
    for (int i = 0; i < 7; i++) { int I[3] = {i, (i + 1) % 7, (i + 3) % 7}; wedge_add(W, I, 1.0, psi.data()); } }
  else if (id == CAL_SPIN7) {                               // Cayley 4-form e_0 ∧ φ + ⋆φ, φ on span(e_1..e_7)
    if (k != 4 || d < 8) throw std::runtime_error("Cayley form needs k = 4, d >= 8");
    for (int i = 0; i < 7; i++) {
      const int l[3] = {i, (i + 1) % 7, (i + 3) % 7};
      int I[4] = {0, l[0] + 1, l[1] + 1, l[2] + 1}; wedge_add(W, I, 1.0, psi.data());
      int perm[7], c[4], p = 3;                             // ⋆(e_a∧e_b∧e_c) = ε e_p∧e_q∧e_r∧e_s
      for (int j = 0; j < 3; j++) perm[j] = l[j];
      for (int j = 0; j < 7; j++) if (j != l[0] && j != l[1] && j != l[2]) perm[p++] = j;
      double eps = 1;                                       // sign of (a,b,c,p,q,r,s) as a permutation of 0..6
      for (int a = 0; a < 7; a++) for (int b = a + 1; b < 7; b++) if (perm[a] > perm[b]) eps = -eps;
      for (int j = 0; j < 4; j++) c[j] = perm[3 + j] + 1;
      wedge_add(W, c, eps, psi.data()); } }
  return psi;
}
// The rung used for a record in dimension d: the largest proper stabiliser available, and the smallest jet
// order that can see it. deff ≥ k is necessary for a non-zero pairing, so a low k reaches more records.
inline void calib_pick(int d, int& k, int& id) {
  if (d < 3)       { k = 0; id = CAL_NONE; }
  else if (d == 8) { k = 4; id = CAL_SPIN7; }
  else if (d >= 7) { k = 3; id = CAL_G2; }
  else if (d == 6) { k = 3; id = CAL_SL; }
  else             { k = 3; id = CAL_SIMPLE; }
}

// dim { A ∈ so(d) : A·ψ = 0 },  (A·ψ)_{i_1…i_k} = Σ_s Σ_j A_{i_s j} ψ_{i_1…j…i_k}
inline int calib_stab_dim(const Wedge& W, const std::vector<double>& psi, double rel_tol = 1e-9) {
  const int d = W.d, k = W.k, ng = d * (d - 1) / 2;
  std::vector<double> M((size_t)W.n * ng, 0.0);
  for (int a = 0, g = 0; a < d; a++) for (int b = a + 1; b < d; b++, g++)
    for (int i = 0; i < W.n; i++) { const int* I = &W.idx[(size_t)i * k]; double v = 0; int J[16];
      for (int s = 0; s < k; s++) { for (int u = 0; u < k; u++) J[u] = I[u];
        if (I[s] == a) { J[s] = b; v += wedge_get(W, psi.data(), J); }
        else if (I[s] == b) { J[s] = a; v -= wedge_get(W, psi.data(), J); } }
      M[(size_t)i * ng + g] = v; }
  std::vector<double> G((size_t)ng * ng, 0.0), w;
  for (int p = 0; p < ng; p++) for (int q = 0; q < ng; q++) { double s = 0;
    for (int i = 0; i < W.n; i++) s += M[(size_t)i * ng + p] * M[(size_t)i * ng + q]; G[(size_t)p * ng + q] = s; }
  la::sym_eig(ng, G, w);
  double top = w.empty() ? 0 : w.back(); int z = 0;
  for (double v : w) if (v <= rel_tol * std::max(top, 1e-300)) z++;
  return z;
}
// ‖Ω·ψ‖∞ — zero exactly when exp(Ωt) preserves the calibration
inline double calib_defect(const Wedge& W, const std::vector<double>& psi, const std::vector<double>& Om) {
  const int d = W.d, k = W.k; double worst = 0; int J[16];
  for (int i = 0; i < W.n; i++) { const int* I = &W.idx[(size_t)i * k]; double v = 0;
    for (int s = 0; s < k; s++) for (int j = 0; j < d; j++) { double a = Om[(size_t)I[s] * d + j]; if (a == 0) continue;
      for (int u = 0; u < k; u++) J[u] = I[u]; J[s] = j; v += a * wedge_get(W, psi.data(), J); }
    worst = std::max(worst, std::fabs(v)); }
  return worst;
}

// Principal frame of A: eigenvectors of Q_ab = Σ_L A_{aL} A_{bL}, most energetic first. Every ψ in the
// ladder is ±1 on coordinate multi-indices, so the maximiser is a signed permutation of this frame plus a
// local rotation — which is what makes the ascent below converge from so few restarts.
inline void jet_frame(const Wedge& W, const double* A, std::vector<double>& E) {
  const int d = W.d, k = W.k; const Wedge& L = wedge_basis(d, k > 1 ? k - 1 : 1);
  std::vector<double> Q((size_t)d * d, 0.0), w; int J[16];
  for (int a = 0; a < d; a++) for (int b = a; b < d; b++) { double s = 0;
    for (int l = 0; l < L.n; l++) { const int* Li = &L.idx[(size_t)l * L.k];
      J[0] = a; for (int i = 0; i < k - 1; i++) J[i + 1] = Li[i]; double va = wedge_get(W, A, J);
      if (va == 0) continue;
      J[0] = b; s += va * wedge_get(W, A, J); }
    Q[(size_t)a * d + b] = Q[(size_t)b * d + a] = s; }
  la::sym_eig(d, Q, w);
  E.assign((size_t)d * d, 0.0);
  for (int a = 0; a < d; a++) for (int t = 0; t < d; t++) E[(size_t)a * d + t] = Q[(size_t)a * d + (d - 1 - t)];
}

// Riemannian ascent of f(R) = ⟨A, R·ψ⟩ = Σ_{I,J} A_I ψ_J det(R[I,J]) over SO(d)
inline double calib_ascend(const Wedge& W, const double* A, const std::vector<double>& psi, int iters, int restarts) {
  const int d = W.d, k = W.k;
  std::vector<int> nzJ; for (int j = 0; j < W.n; j++) if (psi[j] != 0) nzJ.push_back(j);
  std::vector<int> nzI; double amax = 0;
  for (int i = 0; i < W.n; i++) amax = std::max(amax, std::fabs(A[i]));
  for (int i = 0; i < W.n; i++) if (std::fabs(A[i]) > 1e-14 * amax) nzI.push_back(i);
  if (nzI.empty() || nzJ.empty()) return 0.0;
  std::vector<double> R, R2, Rn, Om((size_t)d * d), G((size_t)d * d), Dir, E;
  jet_frame(W, A, E);
  std::vector<int> perm(d); for (int i = 0; i < d; i++) perm[i] = i;
  double S[64], C[64], best = -1e300; la::Rng rng(0x9E3779B97F4A7C15ULL);
  auto eval = [&](const std::vector<double>& Rm, double* Gr) {
    double f = 0; if (Gr) std::fill(Gr, Gr + (size_t)d * d, 0.0);
    for (int ii : nzI) { const int* I = &W.idx[(size_t)ii * k]; const double a = A[ii];
      for (int jj : nzJ) { const int* J = &W.idx[(size_t)jj * k]; const double p = a * psi[jj];
        for (int s = 0; s < k; s++) for (int t = 0; t < k; t++) S[s * k + t] = Rm[(size_t)I[s] * d + J[t]];
        if (!Gr) { f += p * det_small(k, S); continue; }
        f += p * det_cof(k, S, C);
        for (int s = 0; s < k; s++) for (int t = 0; t < k; t++) Gr[(size_t)I[s] * d + J[t]] += p * C[s * k + t]; } }
    return f;
  };
  for (int r0 = 0; r0 < restarts; r0++) {
    R.assign((size_t)d * d, 0.0);
    if (r0 == 0) { for (int i = 0; i < d; i++) R[(size_t)i * d + i] = 1; }        // the record's own frame
    else {
      if (r0 > 1) { for (int i = d - 1; i > 0; i--) std::swap(perm[i], perm[rng.below(i + 1)]); }
      double det = 1;                                                             // signed permutation of E
      for (int i = 0; i < d; i++) for (int j = i + 1; j < d; j++) if (perm[i] > perm[j]) det = -det;
      for (int t = 0; t < d; t++) { double sg = (r0 > 1 && rng.uniform() < 0.5) ? -1.0 : 1.0; det *= sg;
        for (int a = 0; a < d; a++) R[(size_t)a * d + t] = sg * E[(size_t)a * d + perm[t]]; }
      if (det < 0) for (int a = 0; a < d; a++) R[(size_t)a * d + d - 1] = -R[(size_t)a * d + d - 1];
    }
    double f = eval(R, G.data()), eta = 0.5, gprev = 0;
    Dir.assign((size_t)d * d, 0.0);
    for (int it = 0; it < iters; it++) {
      for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) { double m = 0, mt = 0;
        for (int l = 0; l < d; l++) { m += R[(size_t)l * d + i] * G[(size_t)l * d + j]; mt += R[(size_t)l * d + j] * G[(size_t)l * d + i]; }
        Om[(size_t)i * d + j] = 0.5 * (m - mt); }
      double gn2 = 0, gmax = 0;
      for (double v : Om) { gn2 += v * v; gmax = std::max(gmax, std::fabs(v)); }
      if (gmax < 1e-14) break;
      double beta = gprev > 0 ? gn2 / gprev : 0, dg = 0;      // Fletcher–Reeves, restarted if not ascending
      for (size_t i = 0; i < Dir.size(); i++) { Dir[i] = Om[i] + beta * Dir[i]; dg += Dir[i] * Om[i]; }
      if (dg <= 0) Dir = Om;
      gprev = gn2;
      double dmax = 0; for (double v : Dir) dmax = std::max(dmax, std::fabs(v));
      bool up = false;
      for (int ls = 0; ls < 4 && eta > 1e-13; ls++) {
        Rn.assign((size_t)d * d, 0.0);
        for (size_t i = 0; i < Dir.size(); i++) Om[i] = Dir[i] * (eta / dmax);
        la::expm_skew(d, Om, Rn);
        R2.assign((size_t)d * d, 0.0);
        for (int i = 0; i < d; i++) for (int l = 0; l < d; l++) { double rl = R[(size_t)i * d + l]; if (rl == 0) continue;
          for (int j = 0; j < d; j++) R2[(size_t)i * d + j] += rl * Rn[(size_t)l * d + j]; }
        double f2 = eval(R2, nullptr);
        if (f2 > f) { R.swap(R2); f = f2; eta = std::min(2 * eta, 1.0); up = true; break; }
        eta *= 0.25; gprev = 0;
      }
      if (!up) break;
      f = eval(R, G.data());
    }
    best = std::max(best, f);
  }
  return best;
}
// χ* = max over O(d) of ⟨A, R·ψ⟩ — how close the loop's jet moment comes to being ψ-calibrated. A reflection
// is covered by negating the ψ components that contain e_0, so χ* is invariant under the full equivalence
// group of a record (O(d), time shift, reversal, relabelling) and is exactly 0 whenever deff < k.
inline double calib_max(const Wedge& W, const double* A, const std::vector<double>& psi, int iters = 200, int restarts = 3) {
  double na = 0, np = 0;
  for (int i = 0; i < W.n; i++) na += A[i] * A[i];
  for (double v : psi) np += v * v;
  if (W.k >= W.d - 1) return std::sqrt(na * np);            // Λ^{d−1} ≅ Λ^1: every ψ is simple, SO(d) transitive
  std::vector<double> ref(psi);
  for (int i = 0; i < W.n; i++) if (W.idx[(size_t)i * W.k] == 0) ref[i] = -ref[i];
  return std::max(calib_ascend(W, A, psi, iters, restarts), calib_ascend(W, A, ref, iters, restarts));
}

// The distinguished rotating frames of the SU(n) rung: rates in the coordinate planes summing to zero, so
// exp(Ωt) multiplies Θ = (e_0+ie_1)∧… by exp(i Σw t) and preserves Re Θ. The d = 7 g₂ torus w₃ = w₁+w₂ is
// the same condition on the Fano planes (R⁷ = R ⊕ C³ under SU(3) ⊂ G₂).
inline std::vector<double> su_omega(int d, std::vector<double> w) {
  const int n = d / 2;
  if ((int)w.size() > n) throw std::runtime_error("--omega su: more rates than rotation planes");
  w.resize(n, 0.0); double s = 0; for (int i = 0; i + 1 < n; i++) s += w[i];
  w[n - 1] = -s;                                            // the last rate closes Σ w = 0
  std::vector<double> Om((size_t)d * d, 0.0);
  for (int p = 0; p < n; p++) { Om[(size_t)(2 * p) * d + 2 * p + 1] = -w[p]; Om[(size_t)(2 * p + 1) * d + 2 * p] = w[p]; }
  return Om;
}
