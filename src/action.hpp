// Choreography action in a Fourier basis: A = ½∫|q̇|² + ½ Σ_k ∫ |q(t) − q(t+2πk/N)|^{-α}, modes m ≢ 0 (mod N).
// Sampled data are structure-of-arrays with a doubled (wrap-free) layout so the inner loops vectorize.
#pragma once
#include "linalg.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <algorithm>

static constexpr double PI = 3.14159265358979323846;
static constexpr double INF = std::numeric_limits<double>::infinity();

struct Problem {
  int N = 3, d = 2, K = 16, M = 0;
  double alpha = 1.0;
  std::vector<int> modes;
  int nm = 0, nb = 0, n = 0;              // modes, basis functions (2 nm), parameters (nb d)
  std::vector<double> C;                  // nb × M basis samples
  std::vector<double> kin;                // π m²
  std::vector<int> shifts; std::vector<double> shw;
  double wq = 0;                          // 2π/M
  double minsep = 1e-3;                   // collision guard

  static int auto_M(int N, int K) {
    int l = std::lcm(N, 8), target = std::max(64, 8 * K);
    return l * ((target + l - 1) / l);
  }
  void init(int N_, int d_, int K_, int M_ = 0, double alpha_ = 1.0) {
    N = N_; d = d_; K = K_; alpha = alpha_;
    M = M_ > 0 ? M_ : auto_M(N, K);
    if (M % N) throw std::runtime_error("M must be a multiple of N");
    modes.clear();
    for (int m = 1; m <= K; m++) if (m % N) modes.push_back(m);
    nm = (int)modes.size(); nb = 2 * nm; n = nb * d;
    C.assign((size_t)nb * M, 0.0); kin.assign(nb, 0.0);
    for (int mu = 0; mu < nm; mu++) {
      int m = modes[mu];
      for (int j = 0; j < M; j++) {
        double t = 2.0 * PI * (double)(((long)m * j) % M) / M;
        C[(size_t)(2 * mu) * M + j] = std::cos(t);
        C[(size_t)(2 * mu + 1) * M + j] = std::sin(t);
      }
      kin[2 * mu] = kin[2 * mu + 1] = PI * (double)m * m;
    }
    shifts.clear(); shw.clear();
    for (int k = 1; k <= (N - 1) / 2; k++) { shifts.push_back(k * M / N); shw.push_back(1.0); }
    if (N % 2 == 0) { shifts.push_back(M / 2); shw.push_back(0.5); }
    wq = 2.0 * PI / M;
  }
  int idx(int i, int a) const { return i * d + a; }
  void transfer(const Problem& from, const double* xf, double* xt) const {
    std::fill(xt, xt + n, 0.0);
    for (int mu = 0; mu < nm; mu++) {
      auto it = std::lower_bound(from.modes.begin(), from.modes.end(), modes[mu]);
      if (it == from.modes.end() || *it != modes[mu]) continue;
      int mf = (int)(it - from.modes.begin());
      for (int a = 0; a < d; a++) { xt[(2 * mu) * d + a] = xf[(2 * mf) * from.d + a]; xt[(2 * mu + 1) * d + a] = xf[(2 * mf + 1) * from.d + a]; }
    }
  }
};

struct Work {                      // per-thread scratch
  std::vector<double> Q, G, W, HG, rho2, f, fA, fI, dot, U, B, p;
  void resize(const Problem& P) {
    size_t d2M = (size_t)P.d * 2 * P.M;
    Q.resize(d2M); G.resize(d2M); W.resize(d2M); HG.resize(d2M);
    rho2.resize(P.M); f.resize(P.M); fA.resize(P.M); fI.resize(P.M); dot.resize(P.M); p.resize(P.M);
    U.resize((size_t)P.nb * P.M); B.resize((size_t)P.d * (P.d + 1) / 2 * P.M);
  }
};

// Q[a][j], j ∈ [0,2M): doubled so j+shift never wraps
inline void synth(const Problem& P, const double* x, double* Q) {
  const int M = P.M, d = P.d, nb = P.nb;
  for (int a = 0; a < d; a++) {
    double* q = Q + (size_t)a * 2 * M;
    std::fill(q, q + M, 0.0);
    for (int i = 0; i < nb; i++) {
      const double xi = x[i * d + a]; if (xi == 0.0) continue;
      const double* c = &P.C[(size_t)i * M];
#pragma omp simd
      for (int j = 0; j < M; j++) q[j] += xi * c[j];
    }
    std::copy(q, q + M, q + M);
  }
}
// derivative coefficients: d/dt [c cos mt + s sin mt] = (m s) cos mt + (−m c) sin mt
inline void deriv_coeffs(const Problem& P, const double* x, std::vector<double>& xd) {
  xd.assign(P.n, 0.0);
  for (int mu = 0; mu < P.nm; mu++) { double m = P.modes[mu];
    for (int a = 0; a < P.d; a++) { xd[(2 * mu) * P.d + a] = m * x[(2 * mu + 1) * P.d + a]; xd[(2 * mu + 1) * P.d + a] = -m * x[(2 * mu) * P.d + a]; } }
}

// per shift: rho2 and f = w α ρ^{-α-2}; returns Σ_j w ρ^{-α}, or -1 on collision
inline double pair_kernel(const Problem& P, const double* Q, int s, double wgt, double* rho2, double* f, double* fA, double* fI) {
  const int M = P.M, d = P.d;
  std::fill(rho2, rho2 + M, 0.0);
  for (int a = 0; a < d; a++) { const double* q = Q + (size_t)a * 2 * M;
#pragma omp simd
    for (int j = 0; j < M; j++) { double r = q[j] - q[j + s]; rho2[j] += r * r; } }
  double mn = INF;
  for (int j = 0; j < M; j++) mn = std::min(mn, rho2[j]);
  if (mn < P.minsep * P.minsep) return -1.0;
  double U = 0.0;
  if (P.alpha == 1.0) {
#pragma omp simd reduction(+:U)
    for (int j = 0; j < M; j++) { double inv = 1.0 / std::sqrt(rho2[j]), i3 = inv * inv * inv; U += inv; f[j] = wgt * i3;
      if (fA) { fA[j] = 3.0 * wgt * i3 * inv * inv; fI[j] = -wgt * i3; } }
    U *= wgt;
  } else {
    const double al = P.alpha;
    for (int j = 0; j < M; j++) { double pa = std::pow(rho2[j], -0.5 * al); U += pa; double g = al * pa / rho2[j]; f[j] = wgt * g;
      if (fA) { fA[j] = wgt * (al + 2.0) * g / rho2[j]; fI[j] = -wgt * g; } }
    U *= wgt;
  }
  return U;
}

// A(x) and optionally ∇A; +inf on collision
inline double action_grad(const Problem& P, const double* x, double* grad, Work& w, double* kin_out = nullptr, double* pot_out = nullptr) {
  const int M = P.M, d = P.d, nb = P.nb;
  synth(P, x, w.Q.data());
  double* Q = w.Q.data(); double* G = w.G.data();
  if (grad) std::fill(G, G + (size_t)d * 2 * M, 0.0);
  double U = 0.0;
  for (size_t si = 0; si < P.shifts.size(); si++) {
    const int s = P.shifts[si];
    double u = pair_kernel(P, Q, s, P.shw[si], w.rho2.data(), w.f.data(), nullptr, nullptr);
    if (u < 0) return INF;
    U += u;
    if (grad) { const double* f = w.f.data();
      for (int a = 0; a < d; a++) { const double* q = Q + (size_t)a * 2 * M; double* g = G + (size_t)a * 2 * M;
#pragma omp simd
        for (int j = 0; j < M; j++) g[j] -= f[j] * (q[j] - q[j + s]);
#pragma omp simd
        for (int j = 0; j < M; j++) g[j + s] += f[j] * (q[j] - q[j + s]); } }
  }
  double kinE = 0.0;
  for (int i = 0; i < nb; i++) for (int a = 0; a < d; a++) kinE += P.kin[i] * x[i * d + a] * x[i * d + a];
  kinE *= 0.5;
  if (kin_out) *kin_out = kinE;
  if (pot_out) *pot_out = P.wq * U;
  if (grad) {
    for (int a = 0; a < d; a++) { double* g = G + (size_t)a * 2 * M;
#pragma omp simd
      for (int j = 0; j < M; j++) g[j] += g[j + M]; }
    for (int i = 0; i < nb; i++) { const double* c = &P.C[(size_t)i * M];
      for (int a = 0; a < d; a++) { const double* g = G + (size_t)a * 2 * M; double sum = 0.0;
#pragma omp simd reduction(+:sum)
        for (int j = 0; j < M; j++) sum += c[j] * g[j];
        grad[i * d + a] = P.kin[i] * x[i * d + a] + P.wq * sum; } }
  }
  return kinE + P.wq * U;
}

// ∂∇A/∂α: per pair ∂/∂α[−w α ρ^{−α−2} r] = −w ρ^{−α−2}(1 − α ln ρ) r
inline bool grad_alpha(const Problem& P, const double* x, double* ga, Work& w) {
  const int M = P.M, d = P.d, nb = P.nb;
  synth(P, x, w.Q.data()); double* Q = w.Q.data(); double* G = w.G.data();
  std::fill(G, G + (size_t)d * 2 * M, 0.0);
  for (size_t si = 0; si < P.shifts.size(); si++) {
    const int s = P.shifts[si]; const double wgt = P.shw[si];
    if (pair_kernel(P, Q, s, wgt, w.rho2.data(), w.f.data(), nullptr, nullptr) < 0) return false;
    double* f = w.f.data();
    for (int j = 0; j < M; j++) { double rho = std::sqrt(w.rho2[j]); f[j] = wgt * std::pow(rho, -P.alpha - 2.0) * (1.0 - P.alpha * std::log(rho)); }
    for (int a = 0; a < d; a++) { const double* q = Q + (size_t)a * 2 * M; double* g = G + (size_t)a * 2 * M;
      for (int j = 0; j < M; j++) { double fr = f[j] * (q[j] - q[j + s]); g[j] -= fr; g[j + s] += fr; } }
  }
  for (int a = 0; a < d; a++) { double* g = G + (size_t)a * 2 * M; for (int j = 0; j < M; j++) g[j] += g[j + M]; }
  for (int i = 0; i < nb; i++) { const double* c = &P.C[(size_t)i * M];
    for (int a = 0; a < d; a++) { const double* g = G + (size_t)a * 2 * M; double sum = 0.0; for (int j = 0; j < M; j++) sum += c[j] * g[j]; ga[i * d + a] = P.wq * sum; } }
  return true;
}

// out = H(x) v
inline bool hessvec(const Problem& P, const double* x, const double* v, double* out, Work& w) {
  const int M = P.M, d = P.d, nb = P.nb;
  synth(P, x, w.Q.data()); synth(P, v, w.W.data());
  double* Q = w.Q.data(); double* W = w.W.data(); double* HG = w.HG.data();
  std::fill(HG, HG + (size_t)d * 2 * M, 0.0);
  for (size_t si = 0; si < P.shifts.size(); si++) {
    const int s = P.shifts[si];
    if (pair_kernel(P, Q, s, P.shw[si], w.rho2.data(), w.f.data(), w.fA.data(), w.fI.data()) < 0) return false;
    double* dot = w.dot.data(); std::fill(dot, dot + M, 0.0);
    for (int a = 0; a < d; a++) { const double* q = Q + (size_t)a * 2 * M; const double* ww = W + (size_t)a * 2 * M;
#pragma omp simd
      for (int j = 0; j < M; j++) dot[j] += (q[j] - q[j + s]) * (ww[j] - ww[j + s]); }
    const double* fA = w.fA.data(); const double* fI = w.fI.data();
    for (int a = 0; a < d; a++) { const double* q = Q + (size_t)a * 2 * M; const double* ww = W + (size_t)a * 2 * M; double* h = HG + (size_t)a * 2 * M;
#pragma omp simd
      for (int j = 0; j < M; j++) h[j] += fA[j] * (q[j] - q[j + s]) * dot[j] + fI[j] * (ww[j] - ww[j + s]);
#pragma omp simd
      for (int j = 0; j < M; j++) h[j + s] -= fA[j] * (q[j] - q[j + s]) * dot[j] + fI[j] * (ww[j] - ww[j + s]); }
  }
  for (int a = 0; a < d; a++) { double* h = HG + (size_t)a * 2 * M;
#pragma omp simd
    for (int j = 0; j < M; j++) h[j] += h[j + M]; }
  for (int i = 0; i < nb; i++) { const double* c = &P.C[(size_t)i * M];
    for (int a = 0; a < d; a++) { const double* h = HG + (size_t)a * 2 * M; double sum = 0.0;
#pragma omp simd reduction(+:sum)
      for (int j = 0; j < M; j++) sum += c[j] * h[j];
      out[i * d + a] = P.kin[i] * v[i * d + a] + P.wq * sum; } }
  return true;
}

// H(x), n×n row-major
inline bool hessian(const Problem& P, const double* x, std::vector<double>& H, Work& w) {
  const int M = P.M, d = P.d, nb = P.nb, n = P.n;
  H.assign((size_t)n * n, 0.0);
  synth(P, x, w.Q.data());
  double* Q = w.Q.data(); double* U = w.U.data(); double* B = w.B.data(); double* p = w.p.data();
  for (size_t si = 0; si < P.shifts.size(); si++) {
    const int s = P.shifts[si];
    if (pair_kernel(P, Q, s, P.shw[si], w.rho2.data(), w.f.data(), w.fA.data(), w.fI.data()) < 0) return false;
    for (int i = 0; i < nb; i++) { const double* c = &P.C[(size_t)i * M]; double* u = U + (size_t)i * M;
      for (int j = 0; j < M; j++) u[j] = c[j] - c[(j + s) % M]; }
    int ab = 0;
    for (int a = 0; a < d; a++) for (int b = a; b < d; b++, ab++) {
      const double* qa = Q + (size_t)a * 2 * M; const double* qb = Q + (size_t)b * 2 * M; double* Bab = B + (size_t)ab * M;
      const double* fA = w.fA.data(); const double* fI = w.fI.data();
      if (a == b) { for (int j = 0; j < M; j++) { double ra = qa[j] - qa[j + s]; Bab[j] = fA[j] * ra * ra + fI[j]; } }
      else { for (int j = 0; j < M; j++) Bab[j] = fA[j] * (qa[j] - qa[j + s]) * (qb[j] - qb[j + s]); }
    }
    for (int i = 0; i < nb; i++) for (int i2 = i; i2 < nb; i2++) {
      const double* u1 = U + (size_t)i * M; const double* u2 = U + (size_t)i2 * M;
#pragma omp simd
      for (int j = 0; j < M; j++) p[j] = u1[j] * u2[j];
      int ab2 = 0;
      for (int a = 0; a < d; a++) for (int b = a; b < d; b++, ab2++) {
        const double* Bab = B + (size_t)ab2 * M; double h = 0.0;
#pragma omp simd reduction(+:h)
        for (int j = 0; j < M; j++) h += p[j] * Bab[j];
        H[(size_t)(i * d + a) * n + i2 * d + b] += h;
        if (a != b) H[(size_t)(i * d + b) * n + i2 * d + a] += h;
      }
    }
  }
  for (int i = 0; i < nb; i++) for (int i2 = i; i2 < nb; i2++) for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) {
    size_t r = (size_t)(i * d + a) * n + i2 * d + b;
    H[r] *= P.wq;
    if (i2 > i) H[(size_t)(i2 * d + b) * n + i * d + a] = H[r];
  }
  for (int i = 0; i < nb; i++) for (int a = 0; a < d; a++) H[(size_t)(i * d + a) * n + i * d + a] += P.kin[i];
  return true;
}

// Symmetry generator (ε, θ=2πp/q, R): q(εt+θ) = R q(t). The group is closed per mode and averaged into the
// projector onto the fixed subspace. DSL: generators ';'-separated; tokens t±p/q, s[±i,±j,...] (signed
// permutation), r(i,j,p/q) (rotation in plane (i,j), 1-based); spatial tokens compose left to right.
struct SymGen { int eps = 1, p = 0, q = 1; std::vector<double> R; };

struct Symmetry {
  int d = 0; std::vector<SymGen> gens; std::string text;
  bool empty() const { return gens.empty(); }

  static std::vector<double> ident(int d) { std::vector<double> I((size_t)d * d, 0.0); for (int i = 0; i < d; i++) I[(size_t)i * d + i] = 1; return I; }
  static std::vector<double> matmul(int d, const std::vector<double>& A, const std::vector<double>& B) {
    std::vector<double> C((size_t)d * d, 0.0);
    for (int i = 0; i < d; i++) for (int k = 0; k < d; k++) { double a = A[(size_t)i * d + k]; if (a == 0) continue; for (int j = 0; j < d; j++) C[(size_t)i * d + j] += a * B[(size_t)k * d + j]; }
    return C;
  }
  static Symmetry parse(const std::string& s, int d) {
    Symmetry S; S.d = d; S.text = s;
    if (s.empty() || s == "none") return S;
    std::stringstream gs(s); std::string gstr;
    while (std::getline(gs, gstr, ';')) {
      std::stringstream ts(gstr); std::string tok; SymGen g; g.R = ident(d); bool any = false;
      while (ts >> tok) {
        any = true;
        if (tok[0] == 't') {
          if (tok.size() < 2 || (tok[1] != '+' && tok[1] != '-')) throw std::runtime_error("bad time token: " + tok);
          g.eps = tok[1] == '+' ? 1 : -1;
          std::string frac = tok.substr(2); size_t sl = frac.find('/');
          g.p = std::stoi(frac.substr(0, sl)); g.q = sl == std::string::npos ? 1 : std::stoi(frac.substr(sl + 1));
          if (g.q <= 0) throw std::runtime_error("bad time token: " + tok);
        } else if (tok[0] == 's') {
          std::vector<double> R((size_t)d * d, 0.0);
          std::string body = tok.substr(tok.find('[') + 1); body = body.substr(0, body.find(']'));
          std::stringstream bs(body); std::string e; int k = 0;
          while (std::getline(bs, e, ',')) { int v = std::stoi(e); if (v == 0 || std::abs(v) > d || k >= d) throw std::runtime_error("bad perm: " + tok);
            R[(size_t)(std::abs(v) - 1) * d + k] = v > 0 ? 1 : -1; k++; }
          if (k != d) throw std::runtime_error("perm needs d entries: " + tok);
          g.R = matmul(d, R, g.R);
        } else if (tok[0] == 'r') {
          std::string body = tok.substr(tok.find('(') + 1); body = body.substr(0, body.find(')'));
          std::stringstream bs(body); std::string e; std::vector<std::string> parts;
          while (std::getline(bs, e, ',')) parts.push_back(e);
          if (parts.size() != 3) throw std::runtime_error("bad rotation: " + tok);
          int i = std::stoi(parts[0]) - 1, j = std::stoi(parts[1]) - 1; size_t sl = parts[2].find('/');
          int p = std::stoi(parts[2].substr(0, sl)), q = sl == std::string::npos ? 1 : std::stoi(parts[2].substr(sl + 1));
          if (i < 0 || j < 0 || i >= d || j >= d || i == j) throw std::runtime_error("bad rotation plane: " + tok);
          double ph = 2.0 * PI * p / q; std::vector<double> R = ident(d);
          R[(size_t)i * d + i] = std::cos(ph); R[(size_t)i * d + j] = -std::sin(ph); R[(size_t)j * d + i] = std::sin(ph); R[(size_t)j * d + j] = std::cos(ph);
          g.R = matmul(d, R, g.R);
        } else throw std::runtime_error("unknown symmetry token: " + tok);
      }
      if (any) S.gens.push_back(g);
    }
    return S;
  }
  // 2d×2d block of g on (c_m, s_m)
  std::vector<double> block(const SymGen& g, int m) const {
    int D = 2 * d; std::vector<double> Bm((size_t)D * D, 0.0);
    double th = 2.0 * PI * (double)(((long)m * g.p) % g.q) / g.q, c = std::cos(th), s = std::sin(th);
    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) {
      double Rt = g.R[(size_t)j * d + i];
      Bm[(size_t)i * D + j] = c * Rt;               Bm[(size_t)i * D + d + j] = s * Rt;
      Bm[(size_t)(d + i) * D + j] = -g.eps * s * Rt; Bm[(size_t)(d + i) * D + d + j] = g.eps * c * Rt;
    }
    return Bm;
  }
  // orthonormal basis (n × r) of the fixed subspace; r = -1 if the group is too large
  std::vector<double> basis(const Problem& P, int& r) const {
    r = P.n; if (gens.empty()) return {};
    int D = 2 * d; std::vector<std::vector<double>> cols; std::vector<int> col_mode;
    for (int mu = 0; mu < P.nm; mu++) {
      int m = P.modes[mu];
      std::vector<std::vector<double>> G; G.push_back(ident(D));
      std::vector<std::vector<double>> gb; for (auto& g : gens) gb.push_back(block(g, m));
      for (size_t e = 0; e < G.size(); e++) for (auto& g : gb) {
        if (G.size() >= 8192) { r = -1; return {}; }
        std::vector<double> F = matmul(D, g, G[e]); bool seen = false;
        for (auto& E : G) { double md = 0; for (size_t k = 0; k < F.size(); k++) md = std::max(md, std::fabs(F[k] - E[k])); if (md < 1e-9) { seen = true; break; } }
        if (!seen) G.push_back(F);
      }
      std::vector<double> Pm((size_t)D * D, 0.0), w;
      for (auto& E : G) for (size_t k = 0; k < Pm.size(); k++) Pm[k] += E[k] / G.size();
      for (int i = 0; i < D; i++) for (int j = i + 1; j < D; j++) { double v = 0.5 * (Pm[(size_t)i * D + j] + Pm[(size_t)j * D + i]); Pm[(size_t)i * D + j] = Pm[(size_t)j * D + i] = v; }
      la::sym_eig(D, Pm, w);
      for (int k = 0; k < D; k++) if (w[k] > 0.5) { std::vector<double> v(D); for (int i = 0; i < D; i++) v[i] = Pm[(size_t)i * D + k]; cols.push_back(v); col_mode.push_back(mu); }
    }
    r = (int)cols.size();
    std::vector<double> Bmat((size_t)P.n * r, 0.0);
    for (int c = 0; c < r; c++) { int base = 2 * col_mode[c] * d; for (int k = 0; k < D; k++) Bmat[(size_t)(base + k) * r + c] = cols[c][k]; }
    return Bmat;
  }
  // random small group; in d ≥ 3 half the draws are irreducible-type (signed d-cycle + signed transposition/flip),
  // whose fixed loops necessarily span R^d
  static std::string random_text(int d, la::Rng& rng) {
    auto frac = [&](int qmax) { int q = 1 + rng.below(qmax); int p = rng.below(q); return std::to_string(p) + "/" + std::to_string(q); };
    // θ = 2πp/ord so that g^ord acts trivially
    auto tpart_ord = [&](int ord) { int p = rng.below(ord); return std::string("t") + (rng.uniform() < 0.5 ? "+" : "-") + std::to_string(p) + "/" + std::to_string(ord); };
    if (d >= 3 && rng.uniform() < 0.5) {
      int nneg = 0; std::string c1 = " s["; for (int i = 0; i < d; i++) { bool neg = rng.uniform() < 0.5; nneg += neg; c1 += (neg ? "-" : "") + std::to_string((i + 1) % d + 1) + (i + 1 < d ? "," : ""); } c1 += "]";
      std::string g1 = tpart_ord(nneg % 2 ? 2 * d : d) + c1;
      std::vector<int> perm(d); std::iota(perm.begin(), perm.end(), 1); int i = rng.below(d), j = rng.below(d - 1); if (j >= i) j++;
      if (rng.uniform() < 0.5) std::swap(perm[i], perm[j]); else perm[i] = -perm[i];
      std::string g2 = tpart_ord(2) + " s["; for (int k = 0; k < d; k++) g2 += std::to_string(perm[k]) + (k + 1 < d ? "," : ""); g2 += "]";
      return g1 + "; " + g2;
    }
    int ngen = 1 + (rng.uniform() < 0.35 ? 1 : 0); std::string out;
    for (int g = 0; g < ngen; g++) {
      std::string s = std::string("t") + (rng.uniform() < 0.5 ? "+" : "-") + frac(6);
      int kind = rng.below(3);
      if (kind == 0 || d < 2) {          // signed permutation
        std::vector<int> perm(d); std::iota(perm.begin(), perm.end(), 1);
        for (int i = d - 1; i > 0; i--) std::swap(perm[i], perm[rng.below(i + 1)]);
        s += " s["; for (int i = 0; i < d; i++) { s += (rng.uniform() < 0.5 ? "-" : "") + std::to_string(perm[i]) + (i + 1 < d ? "," : ""); } s += "]";
      } else {                            // plane rotation, optionally with a reflection
        int i = rng.below(d), j = rng.below(d - 1); if (j >= i) j++;
        static const int crys[] = {2, 3, 4, 6}; int q = d >= 3 ? crys[rng.below(4)] : 1 + rng.below(6); int p = 1 + rng.below(std::max(1, q - 1));
        s += " r(" + std::to_string(i + 1) + "," + std::to_string(j + 1) + "," + std::to_string(p) + "/" + std::to_string(q) + ")";
        if (kind == 2) { s += " s["; for (int k = 0; k < d; k++) s += (rng.uniform() < 0.3 ? "-" : "") + std::to_string(k + 1) + (k + 1 < d ? "," : ""); s += "]"; }
      }
      out += (g ? "; " : "") + s;
    }
    return out;
  }
};

// search coordinates y, x = B y (identity without symmetry)
struct Reduced {
  const Problem* P = nullptr; const std::vector<double>* B = nullptr; int r = 0;
  std::vector<double> xf, gf, Hf, T;
  Reduced() {}
  Reduced(const Problem& P_, const std::vector<double>* B_, int r_) : P(&P_), B(B_ && !B_->empty() ? B_ : nullptr), r(B_ && !B_->empty() ? r_ : P_.n) {}
  int dim() const { return r; }
  void expand(const double* y, double* x) const {
    if (!B) { std::copy(y, y + P->n, x); return; }
    for (int i = 0; i < P->n; i++) { double s = 0; const double* b = &(*B)[(size_t)i * r]; for (int c = 0; c < r; c++) s += b[c] * y[c]; x[i] = s; }
  }
  void reduce(const double* g, double* gy) const {
    if (!B) { std::copy(g, g + P->n, gy); return; }
    std::fill(gy, gy + r, 0.0);
    for (int i = 0; i < P->n; i++) { const double* b = &(*B)[(size_t)i * r]; double gi = g[i]; for (int c = 0; c < r; c++) gy[c] += b[c] * gi; }
  }
  double f(const double* y, double* gy, Work& w, double* kin = nullptr, double* pot = nullptr) {
    xf.resize(P->n); gf.resize(P->n); expand(y, xf.data());
    double A = action_grad(*P, xf.data(), gy ? gf.data() : nullptr, w, kin, pot);
    if (gy && std::isfinite(A)) reduce(gf.data(), gy);
    return A;
  }
  std::vector<double> vf, hf;
  bool hv(const double* y, const double* v, double* out, Work& w) {
    xf.resize(P->n); vf.resize(P->n); hf.resize(P->n); expand(y, xf.data()); expand(v, vf.data());
    if (!hessvec(*P, xf.data(), vf.data(), hf.data(), w)) return false;
    reduce(hf.data(), out); return true;
  }
  bool hess(const double* y, std::vector<double>& Hy, Work& w) {
    xf.resize(P->n); expand(y, xf.data());
    if (!hessian(*P, xf.data(), Hf, w)) return false;
    if (!B) { Hy = Hf; return true; }
    int n = P->n; T.assign((size_t)n * r, 0.0);
    for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) { double h = Hf[(size_t)i * n + k]; if (h == 0) continue; const double* b = &(*B)[(size_t)k * r]; double* t = &T[(size_t)i * r]; for (int c = 0; c < r; c++) t[c] += h * b[c]; }
    Hy.assign((size_t)r * r, 0.0);
    for (int i = 0; i < n; i++) { const double* b = &(*B)[(size_t)i * r]; const double* t = &T[(size_t)i * r]; for (int c = 0; c < r; c++) { double bc = b[c]; if (bc == 0) continue; double* h = &Hy[(size_t)c * r]; for (int c2 = 0; c2 < r; c2++) h[c2] += bc * t[c2]; } }
    return true;
  }
};

// random loop: modes ≤ K0, amplitude ∝ m^{-γ}
inline void random_guess(const Problem& P, la::Rng& rng, int K0, double gamma, std::vector<double>& x) {
  x.assign(P.n, 0.0);
  for (int mu = 0; mu < P.nm; mu++) { int m = P.modes[mu]; if (m > K0) break; double amp = std::pow((double)m, -gamma);
    for (int a = 0; a < P.d; a++) { x[(2 * mu) * P.d + a] = amp * rng.normal(); x[(2 * mu + 1) * P.d + a] = amp * rng.normal(); } }
}
// λ minimising A(λx) = λ² K + λ^{-α} U  ⇒  λ^{α+2} = αU/(2K)
inline double optimal_scale(double kinE, double potE, double alpha) { return std::pow(alpha * potE / (2.0 * kinE), 1.0 / (alpha + 2.0)); }
