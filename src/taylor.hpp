// Taylor-series N-body integrator, generic over the scalar type (double / mpreal); in-place arithmetic only.
#pragma once
#include "mpreal.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <atomic>
#include <thread>
#include <chrono>
#include "linalg.hpp"

template <class T> T PI_T();
template <> inline double PI_T<double>() { return 3.14159265358979323846; }
#ifdef HAVE_MPFR
template <> inline mpreal PI_T<mpreal>() { return mpreal::pi(); }
#endif

template <class T>
struct NBody {
  int N, d, order, np; double alpha, aexp;
  std::vector<T> X, V;          // Taylor coefficients, (order+1) × N·d
  std::vector<T> D, S, W;       // per pair: differences, |D|², |D|^{2 aexp}
  std::vector<T> Acc, conv, cw;
  T sk, acc, tmp, invs0, h, t, dt, px, pv;
  std::vector<T> p2, v2;
  NBody(int N_, int d_, double alpha_, int order_) : N(N_), d(d_), order(order_), np(N_ * (N_ - 1) / 2), alpha(alpha_), aexp(-(alpha_ + 2.0) / 2.0) {
    size_t nd = (size_t)N * d, K = order + 1;
    X.assign(K * nd, T(0)); V.assign(K * nd, T(0)); D.assign(K * np * d, T(0)); S.assign(K * np, T(0)); W.assign(K * np, T(0));
    Acc.assign(nd, T(0)); conv.assign(d, T(0)); p2.assign(nd, T(0)); v2.assign(nd, T(0));
    cw.assign(K * K, T(0));
    for (int k = 1; k <= order; k++) for (int j = 0; j < k; j++) set_d(cw[(size_t)k * K + j], aexp * (k - j) - j);   // w_k = Σ_j cw[k][j] s_{k−j} w_j / (k s_0)
  }
  T& x(int k, int i, int a) { return X[((size_t)k * N + i) * d + a]; }
  T& v(int k, int i, int a) { return V[((size_t)k * N + i) * d + a]; }
  T& dd(int k, int p, int a) { return D[((size_t)k * np + p) * d + a]; }
  T& s(int k, int p) { return S[(size_t)k * np + p]; }
  T& w(int k, int p) { return W[(size_t)k * np + p]; }

  void series(const T* pos, const T* vel) {
    const size_t nd = (size_t)N * d, K = order + 1;
    for (size_t i = 0; i < nd; i++) { set(X[i], pos[i]); set(V[i], vel[i]); }
    for (int k = 0; k < order; k++) {
      int p = 0;
      for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++, p++) {
        for (int a = 0; a < d; a++) sub(dd(k, p, a), x(k, i, a), x(k, l, a));
        set_zero(sk);                                                         // s_k = Σ_j D_j·D_{k−j}
        for (int j = 0; j < (k + 1) / 2; j++) for (int a = 0; a < d; a++) fma_add(sk, dd(j, p, a), dd(k - j, p, a));
        mul_d(sk, sk, 2.0);
        if (k % 2 == 0) for (int a = 0; a < d; a++) fma_add(sk, dd(k / 2, p, a), dd(k / 2, p, a));
        set(s(k, p), sk);
        if (k == 0) {
          if (alpha == 1.0) { sqrt_(tmp, sk); mul(tmp, tmp, sk); div(w(0, p), T(1), tmp); }       // s^{-3/2}
          else pow_d(w(0, p), sk, aexp);
        } else {
          set_zero(acc);
          for (int j = 0; j < k; j++) { mul(tmp, s(k - j, p), w(j, p)); fma_add(acc, tmp, cw[(size_t)k * K + j]); }
          div(tmp, acc, s(0, p)); div_ui(w(k, p), tmp, (unsigned)k);
        }
      }
      for (size_t i = 0; i < nd; i++) set_zero(Acc[i]);
      p = 0;
      for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++, p++) {
        for (int a = 0; a < d; a++) { set_zero(conv[a]); for (int j = 0; j <= k; j++) fma_add(conv[a], w(j, p), dd(k - j, p, a)); if (alpha != 1.0) mul_d(conv[a], conv[a], alpha); }
        for (int a = 0; a < d; a++) { sub_inplace(Acc[(size_t)i * d + a], conv[a]); add_inplace(Acc[(size_t)l * d + a], conv[a]); }
      }
      for (int i = 0; i < N; i++) for (int a = 0; a < d; a++) { div_ui(x(k + 1, i, a), v(k, i, a), (unsigned)(k + 1)); div_ui(v(k + 1, i, a), Acc[(size_t)i * d + a], (unsigned)(k + 1)); }
    }
  }
  // log2tol: log2 of the absolute tolerance
  double stepsize(double log2tol) {
    size_t nd = (size_t)N * d; double lp = -1e300, lp1 = -1e300;
    for (size_t i = 0; i < nd; i++) {
      lp = std::max(lp, std::max(log2abs(X[(size_t)order * nd + i]), log2abs(V[(size_t)order * nd + i])));
      lp1 = std::max(lp1, std::max(log2abs(X[(size_t)(order - 1) * nd + i]), log2abs(V[(size_t)(order - 1) * nd + i])));
    }
    double l2h = std::min((log2tol - lp) / order, (log2tol - lp1) / (order - 1));
    return 0.9 * std::exp2(std::min(l2h, 8.0));
  }
  void eval(const T& hh, T* pos, T* vel) {
    size_t nd = (size_t)N * d;
    for (size_t i = 0; i < nd; i++) {
      set(px, X[(size_t)order * nd + i]); set(pv, V[(size_t)order * nd + i]);
      for (int k = order - 1; k >= 0; k--) { mul(px, px, hh); add_inplace(px, X[(size_t)k * nd + i]); mul(pv, pv, hh); add_inplace(pv, V[(size_t)k * nd + i]); }
      set(pos[i], px); set(vel[i], pv);
    }
  }
  // optional dense output at sorted times ts
  int integrate(std::vector<T>& pos, std::vector<T>& vel, const T& tend, double tol, const std::vector<T>* ts = nullptr, std::vector<T>* samp = nullptr) {
    size_t nd = (size_t)N * d;
    double scale = 0; for (size_t i = 0; i < nd; i++) scale = std::max(scale, std::max(std::fabs(to_double(pos[i])), std::fabs(to_double(vel[i]))));
    double log2tol = std::log2(tol * std::max(1.0, scale));
    set_zero(t); int steps = 0; size_t si = 0;
    if (samp) samp->assign(ts->size() * 2 * nd, T(0));
    while (t < tend) {
      series(pos.data(), vel.data());
      set_d(h, stepsize(log2tol));
      add(tmp, t, h); if (tmp > tend) sub(h, tend, t);
      if (ts) { add(tmp, t, h); while (si < ts->size() && (*ts)[si] <= tmp) { sub(dt, (*ts)[si], t); eval(dt, &(*samp)[si * 2 * nd], &(*samp)[si * 2 * nd + nd]); si++; } }
      eval(h, p2.data(), v2.data());
      for (size_t i = 0; i < nd; i++) { set(pos[i], p2[i]); set(vel[i], v2[i]); }
      add_inplace(t, h); steps++;
      if (steps > 2000000) break;
    }
    return steps;
  }
};

template <class T>
T nbody_energy(int N, int d, double alpha, const T* pos, const T* vel) {
  T E(0);
  for (int i = 0; i < N; i++) for (int a = 0; a < d; a++) E += vel[i * d + a] * vel[i * d + a] * 0.5;
  for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++) { T r2(0); for (int a = 0; a < d; a++) { T df = pos[i * d + a] - pos[l * d + a]; r2 += df * df; }
    if (alpha == 1.0) E -= T(1) / sqrt(r2); else E -= pow(r2, T(-alpha / 2)); }
  return E;
}

// shooting residual F = Φ_{2π/N}(Z) − G S Z, (SZ)_j = Z_{j+1}, G = exp(2πΩ/N) (identity when Ω = 0);
// returns max|F|
template <class T>
double chore_residual(NBody<T>& nb, const std::vector<T>& pos, const std::vector<T>& vel, double tol, const std::vector<T>* G = nullptr, std::vector<T>* F = nullptr) {
  int N = nb.N, d = nb.d; size_t nd = (size_t)N * d;
  std::vector<T> p = pos, v = vel;
  T tend = T(2) * T(PI_T<T>()) / N;
  nb.integrate(p, v, tend, tol);
  double mx = 0; if (F) F->assign(2 * nd, T(0));
  for (int j = 0; j < N; j++) for (int a = 0; a < d; a++) { int j1 = (j + 1) % N;
    T tp = pos[j1 * d + a], tv = vel[j1 * d + a];
    if (G) { tp = T(0); tv = T(0);
      for (int b = 0; b < d; b++) { tp += (*G)[(size_t)a * d + b] * pos[j1 * d + b]; tv += (*G)[(size_t)a * d + b] * vel[j1 * d + b]; } }
    T fp = p[j * d + a] - tp, fv = v[j * d + a] - tv;
    if (F) { (*F)[j * d + a] = fp; (*F)[nd + j * d + a] = fv; }
    double afp = std::fabs(to_double(fp)), afv = std::fabs(to_double(fv));
    if (!(afp < INFINITY) || !(afv < INFINITY)) return INFINITY;          // collision
    mx = std::max(mx, std::max(afp, afv)); }
  return mx;
}

template <class T> T two_pow_T(long e);
template <> inline double two_pow_T<double>(long e) { return std::ldexp(1.0, (int)e); }
#ifdef HAVE_MPFR
template <> inline mpreal two_pow_T<mpreal>(long e) { return mpreal::two_pow(e); }
#endif

template <class T> struct ShootWork { std::vector<T> F, Fn, Fp, Fm, J, JtJ, A, rhs, rhs0, Zp, Zm, Z0; };

// Newton–LM on F(Z) = Φ_{2π/N}(Z) − G S Z; central-difference Jacobian (step 2^hstep_log2), damping
// μ·max diag(JᵀJ) from 2^mu_log2, raised ×8 and retried against the same Jacobian until ‖F‖ falls (fixed
// damping abandoned 27 % of admissible candidates). Returns final max|F|, INFINITY on blow-up.
template <class T>
double shoot_newton(NBody<T>& nb, std::vector<T>& Z, double itol, int max_iter, double target, long hstep_log2, long mu_log2, ShootWork<T>& W, bool verbose = false, const std::vector<T>* G = nullptr, int threads = 1, double stall = 0, int lm_retries = 8) {
  const int N = nb.N, d = nb.d, nd = N * d, n2 = 2 * nd;
  auto residual = [&](const std::vector<T>& Zc, std::vector<T>& F) { std::vector<T> p(Zc.begin(), Zc.begin() + nd), v(Zc.begin() + nd, Zc.end()); return chore_residual(nb, p, v, itol, G, &F); };
  auto remove_cm = [&](std::vector<T>& Zc) { for (int half = 0; half < 2; half++) for (int c = 0; c < d; c++) { T m(0); for (int k = 0; k < N; k++) m += Zc[half * nd + k * d + c]; m /= N; for (int k = 0; k < N; k++) Zc[half * nd + k * d + c] -= m; } };
  W.J.assign((size_t)n2 * n2, T(0)); W.JtJ.assign((size_t)n2 * n2, T(0)); W.rhs.assign(n2, T(0));
  remove_cm(Z);
  auto t0 = std::chrono::steady_clock::now();
  auto secs = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
  double maxF = residual(Z, W.F); if (!(maxF < INFINITY)) return INFINITY;
  if (verbose) { std::printf("  iter 0: residual %.3e  [%.1fs]\n", maxF, secs()); std::fflush(stdout); }
  T hstep = two_pow_T<T>(hstep_log2), mu = two_pow_T<T>(mu_log2), mu0 = two_pow_T<T>(mu_log2), tmp(0);
  int flat = 0;                    // consecutive iterations that bought less than a factor `stall`
  for (int it = 1; it <= max_iter && maxF > target; it++) {
    if (threads <= 1) {
      for (int c = 0; c < n2; c++) {
        W.Zp = Z; W.Zp[c] += hstep; if (!(residual(W.Zp, W.Fp) < INFINITY)) return INFINITY;
        W.Zm = Z; W.Zm[c] -= hstep; if (!(residual(W.Zm, W.Fm) < INFINITY)) return INFINITY;
        for (int r = 0; r < n2; r++) { sub(tmp, W.Fp[r], W.Fm[r]); div(tmp, tmp, hstep); div_ui(W.J[(size_t)r * n2 + c], tmp, 2u); }
      }
    } else {
      // one column per task; W.J is preallocated at the working precision, so no element is resized
      std::atomic<int> nextc{0}; std::atomic<bool> blew{false};
      std::vector<std::thread> pool; pool.reserve(threads);
      for (int t = 0; t < threads; t++) pool.emplace_back([&] {
        NBody<T> nbl(nb.N, nb.d, nb.alpha, nb.order);
        std::vector<T> Zc, p, v, Fp, Fm; T tl(0);
        for (int c = nextc.fetch_add(1); c < n2 && !blew.load(); c = nextc.fetch_add(1)) {
          for (int sgn = 0; sgn < 2; sgn++) {
            Zc = Z; if (sgn) Zc[c] -= hstep; else Zc[c] += hstep;
            p.assign(Zc.begin(), Zc.begin() + nd); v.assign(Zc.begin() + nd, Zc.end());
            if (!(chore_residual(nbl, p, v, itol, G, sgn ? &Fm : &Fp) < INFINITY)) { blew = true; break; }
          }
          if (blew.load()) break;
          for (int r = 0; r < n2; r++) { sub(tl, Fp[r], Fm[r]); div(tl, tl, hstep); div_ui(W.J[(size_t)r * n2 + c], tl, 2u); }
        }
      });
      for (auto& t : pool) t.join();
      if (blew.load()) return INFINITY;
    }
    T dmax(0);
    for (int i = 0; i < n2; i++) { for (int j = i; j < n2; j++) { T& acc = W.JtJ[(size_t)i * n2 + j]; set_zero(acc); for (int k = 0; k < n2; k++) fma_add(acc, W.J[(size_t)k * n2 + i], W.J[(size_t)k * n2 + j]); set(W.JtJ[(size_t)j * n2 + i], acc); }
      if (W.JtJ[(size_t)i * n2 + i] > dmax) set(dmax, W.JtJ[(size_t)i * n2 + i]); }
    W.rhs0.resize(n2);
    for (int i = 0; i < n2; i++) { T& r = W.rhs0[i]; set_zero(r); for (int k = 0; k < n2; k++) fma_sub(r, W.J[(size_t)k * n2 + i], W.F[k]); }
    W.Z0 = Z; bool ok = false; double newF = INFINITY, dn = 0;
    for (int t = 0; t < lm_retries && !ok; t++) {                      // LM retries reuse the Jacobian
      W.A = W.JtJ; mul(tmp, mu, dmax); for (int i = 0; i < n2; i++) add_inplace(W.A[(size_t)i * n2 + i], tmp);
      W.rhs = W.rhs0; if (!la::lu_solve(n2, W.A, W.rhs)) break;
      Z = W.Z0; dn = 0; for (int i = 0; i < n2; i++) { Z[i] += W.rhs[i]; dn = std::max(dn, std::fabs(to_double(W.rhs[i]))); }
      remove_cm(Z); newF = residual(Z, W.Fn);
      if (newF < maxF) { ok = true; W.F.swap(W.Fn); if (mu > mu0) mul_d(mu, mu, 0.2); } else mul_d(mu, mu, 8.0);
    }
    if (verbose) { std::printf("  iter %d: residual %.3e  step %.3e  [%.1fs]\n", it, newF, dn, secs()); std::fflush(stdout); }
    if (!ok) { if (newF < maxF) maxF = newF; else Z = W.Z0; break; }
    // a residual that has stopped falling is at this problem's floor, not on its way down
    flat = (stall > 0 && newF > stall * maxF) ? flat + 1 : 0;
    maxF = newF;
    if (flat >= 2) { if (verbose) std::printf("  residual has stopped falling — stopping at %.3e\n", maxF); break; }
  }
  return maxF;
}
