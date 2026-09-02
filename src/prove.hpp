// Existence proof by interval arithmetic. A choreography is a zero of F(Z) = Φ_{2π/N}(Z) − S Z on the 2Nd
// state (S shifts the bodies). The proof is Krawczyk's test on a box around the refined state: the flow and its
// derivative over the whole box are enclosed by an interval Taylor method (the series recurrence of taylor.hpp
// run on intervals, its linearisation for the variational flow, Lagrange remainders from a Picard enclosure of
// the step). Time shift, translation and rotation make DF singular, so the test runs on a slice transverse to
// the group orbit with as many equations dropped as generators; the dropped ones follow from the conservation of
// energy, momentum and angular momentum, closed by a second, tiny interval Newton argument.
#pragma once
#ifdef HAVE_MPFR
#include "interval.hpp"
#include "taylor.hpp"
#include <thread>
#include <atomic>
#include <functional>
#include <chrono>
#include <memory>

template <> inline ival PI_T<ival>() { return ival::pi(); }

// Linearisation of NBody<T>::series along one direction: the same recurrences, differentiated. δ-quantities
// share the coefficient layout of the primal series they are read from.
template <class T>
struct Tangent {
  int N, d, order, np, nd; double alpha, aexp;
  std::vector<T> X, V, D, S, W, Acc, conv; T tmp, acc;
  Tangent(int N_, int d_, double alpha_, int order_) : N(N_), d(d_), order(order_), np(N_ * (N_ - 1) / 2), nd(N_ * d_), alpha(alpha_), aexp(-(alpha_ + 2.0) / 2.0) {
    const size_t K = order + 1;
    X.assign(K * nd, T(0)); V.assign(K * nd, T(0)); D.assign(K * np * d, T(0)); S.assign(K * np, T(0)); W.assign(K * np, T(0)); Acc.assign(nd, T(0)); conv.assign(d, T(0));
  }
  T& x(int k, int i, int a) { return X[((size_t)k * N + i) * d + a]; }
  T& v(int k, int i, int a) { return V[((size_t)k * N + i) * d + a]; }
  T& dd(int k, int p, int a) { return D[((size_t)k * np + p) * d + a]; }
  T& s(int k, int p) { return S[(size_t)k * np + p]; }
  T& w(int k, int p) { return W[(size_t)k * np + p]; }
  // nb holds a computed series; (dpos, dvel) is the direction at t = 0
  void series(NBody<T>& nb, const T* dpos, const T* dvel) {
    const size_t K = order + 1;
    for (int i = 0; i < nd; i++) { set(X[i], dpos[i]); set(V[i], dvel[i]); }
    for (int k = 0; k < order; k++) {
      int p = 0;
      for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++, p++) {
        for (int a = 0; a < d; a++) sub(dd(k, p, a), x(k, i, a), x(k, l, a));
        set_zero(acc);                                                          // δs_k = 2 Σ_j D_j·δD_{k−j}
        for (int j = 0; j <= k; j++) for (int a = 0; a < d; a++) fma_add(acc, nb.dd(j, p, a), dd(k - j, p, a));
        mul_d(s(k, p), acc, 2.0);
        if (k == 0) { div(tmp, nb.w(0, p), nb.s(0, p)); mul_d(tmp, tmp, aexp); mul(w(0, p), tmp, s(0, p)); }   // δw₀ = a w₀/s₀ δs₀
        else {                                                                  // δw_k = (δA_k − k w_k δs₀)/(k s₀)
          set_zero(acc);
          for (int j = 0; j < k; j++) { mul(tmp, s(k - j, p), nb.w(j, p)); fma_add(tmp, nb.s(k - j, p), w(j, p)); fma_add(acc, tmp, nb.cw[(size_t)k * K + j]); }
          mul(tmp, nb.w(k, p), s(0, p)); mul_d(tmp, tmp, (double)k); sub_inplace(acc, tmp);
          div(tmp, acc, nb.s(0, p)); div_ui(w(k, p), tmp, (unsigned)k);
        }
      }
      for (int i = 0; i < nd; i++) set_zero(Acc[i]);
      p = 0;
      for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++, p++) {
        for (int a = 0; a < d; a++) { set_zero(conv[a]);
          for (int j = 0; j <= k; j++) { fma_add(conv[a], w(j, p), nb.dd(k - j, p, a)); fma_add(conv[a], nb.w(j, p), dd(k - j, p, a)); }
          if (alpha != 1.0) mul_d(conv[a], conv[a], alpha); }
        for (int a = 0; a < d; a++) { sub_inplace(Acc[(size_t)i * d + a], conv[a]); add_inplace(Acc[(size_t)l * d + a], conv[a]); }
      }
      for (int i = 0; i < N; i++) for (int a = 0; a < d; a++) { div_ui(x(k + 1, i, a), v(k, i, a), (unsigned)(k + 1)); div_ui(v(k + 1, i, a), Acc[(size_t)i * d + a], (unsigned)(k + 1)); }
    }
  }
};

// Σ_{k<p} h^k c_k by Horner; c_k at coef[k·stride + i]
template <class T>
inline void horner(const std::vector<T>& coef, size_t stride, int p, const T& h, T* out, T& px) {
  for (size_t i = 0; i < stride; i++) { set(px, coef[(size_t)(p - 1) * stride + i]);
    for (int k = p - 2; k >= 0; k--) { mul(px, px, h); add_inplace(px, coef[(size_t)k * stride + i]); } set(out[i], px); }
}

// Validated flow of the N-body problem: state box Z (2Nd), m tangent columns Psi (directions of the box),
// and the action integral. Each step: order-p Taylor polynomial on the tight box, Lagrange remainder h^p c_p
// on a rough enclosure W of the step (first-order Picard for the state, Gronwall for the tangents).
struct Verified {
  int N, d, nd, n2, order, m, threads; double alpha;
  NBody<ival> nb, nbw, nbq; Tangent<ival> tg1; std::vector<std::unique_ptr<Tangent<ival>>> tg;
  std::vector<ival> Z, Psi, Zn, Psin, Wz, Wpsi, tmpv; ival A, t, px, hI, hp, coefp;
  double log2tol, hprev = 0; long steps = 0, rejects = 0; double maxwid = 0, hmin = INFINITY, hmax = 0; bool debug = false, took_last = false;
  Verified(int N_, int d_, double alpha_, int order_, int m_, int threads_, double tol)
    : N(N_), d(d_), nd(N_ * d_), n2(2 * N_ * d_), order(order_), m(m_), threads(std::max(1, threads_)), alpha(alpha_),
      nb(N_, d_, alpha_, order_), nbw(N_, d_, alpha_, order_), nbq(N_, d_, alpha_, std::min(order_, 6)), tg1(N_, d_, alpha_, 1), log2tol(std::log2(tol)) {
    for (int i = 0; i < threads; i++) tg.emplace_back(std::make_unique<Tangent<ival>>(N_, d_, alpha_, order_));
    Z.assign(n2, ival(0)); Psi.assign((size_t)m * n2, ival(0)); Zn = Z; Psin = Psi; Wz = Z; Wpsi = Psi; tmpv.assign(n2, ival(0));
  }
  ival* psi(int j) { return &Psi[(size_t)j * n2]; }
  // rough enclosure W of the step: z(τ) ∈ Σ_{k<q} [0,h]^k Z_k(Z) + [0,h]^q Z_q(W) ⊂ int W for τ ∈ [0,h] (Taylor with
  // Lagrange remainder, closed by continuation). The operator contracts like (h/ρ)^q, so the step is set by the
  // tolerance, not by h‖Df‖ < 1. Leaves nbq holding the order-q series on the final W.
  bool picard(double h) {
    const ival span(0.0, h); const int q = nbq.order; ival sq(ipow(span, q)), c;
    horner(nb.X, nd, q, span, tmpv.data(), px); horner(nb.V, nd, q, span, tmpv.data() + nd, px);      // the fixed part
    for (int i = 0; i < n2; i++) { mul(c, sq, i < nd ? nb.X[(size_t)q * nd + i] : nb.V[(size_t)q * nd + i - nd]); add(Wz[i], tmpv[i], c); Wz[i].inflate(0.0, 0.1 * Wz[i].wid() + 1e-3 * h * (1 + nb.V[i].mag())); }
    for (int it = 0; it < 20; it++) {
      nbq.series(Wz.data(), Wz.data() + nd);
      bool inside = true;
      for (int i = 0; i < n2; i++) { mul(c, sq, i < nd ? nbq.X[(size_t)q * nd + i] : nbq.V[(size_t)q * nd + i - nd]); add_inplace(c, tmpv[i]);
        if (!c.finite()) return false;
        if (!Wz[i].interior(c)) { inside = false; set(Wz[i], c); Wz[i].inflate(0.0, 0.1 * c.wid() + 1e-3 * h * (1 + nbq.V[i].mag())); } }
      if (inside) return true;
    }
    return false;
  }
  // ‖Df(W)‖_∞ from the order-1 tangent series on W with unit directions
  double lipschitz() {
    std::vector<double> row(nd, 0.0); std::vector<ival> e(n2, ival(0));
    for (int j = 0; j < nd; j++) { set_d(e[j], 1.0); tg1.series(nbq, e.data(), e.data() + nd); set_zero(e[j]);
      for (int i = 0; i < nd; i++) row[i] += tg1.v(1, i / d, i % d).mag(); }
    double a = 1.0; for (double r : row) a = std::max(a, r); return a;
  }
  // one step of size ≤ h; returns the step taken (0 on failure)
  double step(double h, bool last, const ival& hlast) {   // nb holds the series on Z
    const int p = order;
    if (hprev > 0 && h > hprev) { h = hprev; last = false; }      // hprev: the controller's prediction
    double scale = 1.0; for (int i = 0; i < n2; i++) scale = std::max(scale, Z[i].mag());
    for (int attempt = 0; attempt < 40; attempt++, rejects++) {
      if (!last) set_d(hI, h); else set(hI, hlast);
      if (!picard(h)) { if (debug) std::printf("    step %ld: picard failed at h=%.3e\n", steps, h); h *= 0.5; last = false; continue; }
      nbw.series(Wz.data(), Wz.data() + nd);
      // the step is set by the remainder as the intervals see it on W, which dependency inflates far beyond
      // the tight-box estimate: shrink to land it at the tolerance
      set(hp, ipow(hI, p)); double rem = 0;
      for (int i = 0; i < n2; i++) { mul(coefp, hp, i < nd ? nbw.X[(size_t)p * nd + i] : nbw.V[(size_t)p * nd + i - nd]); rem = std::max(rem, coefp.mag()); }
      if (!(rem <= std::exp2(log2tol) * scale)) { const double f = std::isfinite(rem) ? std::max(0.3, 0.8 * std::pow(std::exp2(log2tol) * scale / rem, 1.0 / p)) : 0.5;
        if (debug) std::printf("    step %ld: remainder %.2e at h=%.3e, shrinking by %.2f\n", steps, rem, h, f); h *= f; last = false; continue; }
      const double a = lipschitz(), gron = 1.001 * std::exp(a * h);                 // Gronwall: |ψ(τ) − ψ(0)| ≤ τ a ‖ψ(0)‖ e^{aτ}
      if (a * h > 4.0) { if (debug) std::printf("    step %ld: a h = %.3e at h=%.3e\n", steps, a * h, h); h *= 0.5; last = false; continue; }
      // state: polynomial on Z, remainder on W
      bool ok = true;
      horner(nb.X, nd, p, hI, Zn.data(), px); horner(nb.V, nd, p, hI, Zn.data() + nd, px);
      for (int i = 0; i < n2 && ok; i++) { mul(coefp, hp, i < nd ? nbw.X[(size_t)p * nd + i] : nbw.V[(size_t)p * nd + i - nd]); add_inplace(Zn[i], coefp); ok = Zn[i].finite(); }
      if (!ok) { if (debug) std::printf("    step %ld: state remainder not finite at h=%.3e\n", steps, h); h *= 0.5; last = false; continue; }
      // tangents: Gronwall box for the step, then the same polynomial + remainder per column
      std::atomic<int> next{0}; std::atomic<bool> bad{false};
      auto work = [&](int tid) { Tangent<ival>& tg_ = *tg[tid]; ival q, hq(hp), c;
        for (int j = next.fetch_add(1); j < m && !bad.load(); j = next.fetch_add(1)) {
          const ival* pj = psi(j); ival* wj = &Wpsi[(size_t)j * n2]; ival* nj = &Psin[(size_t)j * n2];
          double mx = 0; for (int i = 0; i < n2; i++) mx = std::max(mx, pj[i].mag());
          const double rho = gron * h * a * mx;
          for (int i = 0; i < n2; i++) { set(wj[i], pj[i]); wj[i].inflate(0.0, rho); }
          tg_.series(nb, pj, pj + nd); horner(tg_.X, nd, p, hI, nj, q); horner(tg_.V, nd, p, hI, nj + nd, q);
          tg_.series(nbw, wj, wj + nd);
          for (int i = 0; i < n2; i++) { mul(c, hq, i < nd ? tg_.X[(size_t)p * nd + i] : tg_.V[(size_t)p * nd + i - nd]); add_inplace(nj[i], c); if (!nj[i].finite()) bad = true; }
        } };
      if (m > 0) { std::vector<std::thread> pool; for (int tid = 1; tid < threads; tid++) pool.emplace_back(work, tid); work(0); for (auto& th : pool) th.join(); }
      if (bad.load()) { if (debug) std::printf("    step %ld: tangent remainder not finite at h=%.3e\n", steps, h); h *= 0.5; last = false; continue; }
      // action: ∫L = Σ_{k<p} h^{k+1} L_k/(k+1) + h^{p+1} L_p(W)/(p+1), L = ½|v|² + Σ s·w
      for (int k = 0; k <= p; k++) { NBody<ival>& src = k < p ? nb : nbw; set_zero(coefp);
        for (int i = 0; i < nd; i++) for (int j = 0; j <= k; j++) fma_add(coefp, src.V[(size_t)j * nd + i], src.V[(size_t)(k - j) * nd + i]);
        mul_d(coefp, coefp, 0.5);
        for (int q = 0; q < nb.np; q++) for (int j = 0; j <= k; j++) fma_add(coefp, src.s(j, q), src.w(k - j, q));
        mul(coefp, coefp, ipow(hI, k + 1)); div_ui(coefp, coefp, (unsigned)(k + 1)); add_inplace(A, coefp); }
      Z.swap(Zn); Psi.swap(Psin); add_inplace(t, hI); steps++;
      for (int i = 0; i < n2; i++) maxwid = std::max(maxwid, Z[i].wid());
      hmin = std::min(hmin, h); hmax = std::max(hmax, h);
      took_last = last; hprev = h * std::min(1.25, 0.9 * std::pow(std::exp2(log2tol) * scale / std::max(rem, 1e-300), 1.0 / p));
      if (debug) std::printf("    step %ld: h=%.3e t=%.6f width %.2e\n", steps, h, t.mid(), maxwid);
      // the step's rough enclosure of body 0's first two coordinates, for anyone drawing the corridor
      if (debug && d >= 2) std::printf("    wbox %.9f %.6e %.9f %.9f %.9f %.9f\n", t.mid() - h, h, mpfr_get_d(Wz[0].lo, MPFR_RNDD), mpfr_get_d(Wz[0].hi, MPFR_RNDU), mpfr_get_d(Wz[1].lo, MPFR_RNDD), mpfr_get_d(Wz[1].hi, MPFR_RNDU));
      return h;
    }
    return 0;
  }
  // flow to tend; false when a step could not be validated (collision reached, or the box blew up)
  bool integrate(const ival& tend) {
    set_zero(t); set_zero(A); ival rem;
    for (;;) {
      sub(rem, tend, t); if (mpfr_sgn(rem.lo) <= 0) return true;
      nb.series(Z.data(), Z.data() + nd);
      double h = nb.stepsize(log2tol);
      const bool last = h >= mpfr_get_d(rem.lo, MPFR_RNDD);
      if (step(last ? mpfr_get_d(rem.hi, MPFR_RNDU) : h, last, rem) <= 0) return false;
      if (took_last) return true;
      if (steps > 200000) return false;
    }
  }
};

// cos and sin of a thin interval: the endpoint values, outward, and the extreme value hulled in when an
// extremum (kπ for cos, π/2 + kπ for sin) lies inside or within a hair of the interval — a zero crossing is
// monotone and needs nothing, so a quarter turn stays sharp
inline void icos_isin(const ival& th, ival& c, ival& s) {
  mpfr_cos(c.lo, th.lo, MPFR_RNDD); mpfr_cos(c.hi, th.hi, MPFR_RNDU); if (mpfr_cmp(c.lo, c.hi) > 0) mpfr_swap(c.lo, c.hi);
  mpfr_sin(s.lo, th.lo, MPFR_RNDD); mpfr_sin(s.hi, th.hi, MPFR_RNDU); if (mpfr_cmp(s.lo, s.hi) > 0) mpfr_swap(s.lo, s.hi);
  const double lo = mpfr_get_d(th.lo, MPFR_RNDD), hi = mpfr_get_d(th.hi, MPFR_RNDU), tol = 1e-6 + (hi - lo);
  const double kc = std::round(0.5 * (lo + hi) / PI), ks = std::round((0.5 * (lo + hi) - PI / 2) / PI);
  if (std::fabs(0.5 * (lo + hi) - kc * PI) < tol) c.hull(ival(std::fmod(kc, 2.0) == 0.0 ? 1.0 : -1.0));
  if (std::fabs(0.5 * (lo + hi) - (PI / 2 + ks * PI)) < tol) s.hull(ival(std::fmod(ks, 2.0) == 0.0 ? 1.0 : -1.0));
}

// A rotating frame in its own coordinates: planes (2i, 2i+1) turning at rate[i], the remaining axes fixed. The
// angle classes of G = exp(2πΩ/N) decide the gauge group: translations along axes G fixes, the time shift, and
// the rotations commuting with G — each paired with a conserved quantity (momentum, energy, angular momentum).
// Rates are snapped so that coinciding angles coincide exactly; the frame proved is this one.
struct Frame {
  int d = 0, N = 0; std::vector<double> R;                 // rows: frame axes in the record's coordinates
  std::vector<double> rate; std::vector<int> cls;          // per plane; class 0 = angle 0, 1 = π, 2 = generic
  std::vector<std::vector<double>> trans, rots;            // gauge generators c ∈ R^d and ξ ∈ so(d), exact entries
  bool rotating() const { return !rate.empty(); }
  int nplanes() const { return (int)rate.size(); }
  // G at the working type: block rotations
  template <class T, class Cos> void G(std::vector<T>& G, Cos cs) const {
    G.assign((size_t)d * d, T(0)); for (int a = 0; a < d; a++) set_d(G[(size_t)a * d + a], 1.0);
    for (int i = 0; i < nplanes(); i++) { const int u = 2 * i, v = 2 * i + 1;
      if (cls[i] == 0) continue;
      if (cls[i] == 1) { set_d(G[(size_t)u * d + u], -1.0); set_d(G[(size_t)v * d + v], -1.0); continue; }
      T c, sn; cs(rate[i], c, sn);
      set(G[(size_t)u * d + u], c); set(G[(size_t)v * d + v], c); mul_d(G[(size_t)u * d + v], sn, -1.0); set(G[(size_t)v * d + u], sn); }
  }
  void G_ival(std::vector<ival>& Gi) const { G<ival>(Gi, [&](double w, ival& c, ival& s) { ival th = ival(w) * ival::pi() * 2.0 / ival(N); icos_isin(th, c, s); }); }
  void G_mpreal(std::vector<mpreal>& Gm) const { G<mpreal>(Gm, [&](double w, mpreal& c, mpreal& s) { mpreal th = mpreal(w) * mpreal::pi() * 2 / N; c = cos(th); s = sin(th); }); }
  // rotate a 2Nd state into the frame's axes
  void apply(int N_, const double* Z, std::vector<double>& out) const {
    out.assign(2 * (size_t)N_ * d, 0.0);
    for (int h = 0; h < 2; h++) for (int k = 0; k < N_; k++) for (int a = 0; a < d; a++) { double t = 0; for (int b = 0; b < d; b++) t += R[(size_t)a * d + b] * Z[h * N_ * d + k * d + b]; out[h * N_ * d + k * d + a] = t; }
  }
};

inline Frame frame_of(int N, int d, const double* Om, double snap = 1e-9) {
  Frame fr; fr.d = d; fr.N = N; fr.R.assign((size_t)d * d, 0.0); for (int a = 0; a < d; a++) fr.R[(size_t)a * d + a] = 1.0;
  std::vector<std::vector<double>> axes;                                       // frame axes in order: u1 v1 u2 v2 … fixed
  if (Om) {
    std::vector<double> A((size_t)d * d), S((size_t)d * d, 0.0), lam;
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) A[(size_t)a * d + b] = 0.5 * (Om[(size_t)a * d + b] - Om[(size_t)b * d + a]);
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double t = 0; for (int k = 0; k < d; k++) t -= A[(size_t)a * d + k] * A[(size_t)k * d + b]; S[(size_t)a * d + b] = t; }
    la::sym_eig(d, S, lam); double lmax = 0; for (double l : lam) lmax = std::max(lmax, l);
    std::vector<std::pair<double, std::vector<double>>> rem; std::vector<std::vector<double>> fixed;
    for (int j = 0; j < d; j++) { std::vector<double> u(d); for (int a = 0; a < d; a++) u[a] = S[(size_t)a * d + j]; if (lam[j] > 1e-10 * lmax) rem.push_back({lam[j], u}); else fixed.push_back(u); }
    auto dot = [&](const std::vector<double>& x, const std::vector<double>& y) { double t = 0; for (int a = 0; a < d; a++) t += x[a] * y[a]; return t; };
    while (!rem.empty()) {
      std::vector<double> u = rem[0].second, v(d); const double w = std::sqrt(rem[0].first);
      for (int a = 0; a < d; a++) { double t = 0; for (int b = 0; b < d; b++) t += A[(size_t)a * d + b] * u[b]; v[a] = t / w; }
      { double c = dot(v, u); for (int a = 0; a < d; a++) v[a] -= c * u[a]; double n = std::sqrt(dot(v, v)); for (double& x : v) x /= n; }
      axes.push_back(u); axes.push_back(v); fr.rate.push_back(w);
      std::vector<std::pair<double, std::vector<double>>> next;
      for (size_t r = 1; r < rem.size(); r++) { std::vector<double> x = rem[r].second;
        for (const std::vector<double>* p : {&u, &v}) { double c = dot(x, *p); for (int a = 0; a < d; a++) x[a] -= c * (*p)[a]; }
        for (auto& q : next) { double c = dot(x, q.second); for (int a = 0; a < d; a++) x[a] -= c * q.second[a]; }
        double n = std::sqrt(dot(x, x)); if (n < 1e-6) continue; for (double& y : x) y /= n; next.push_back({rem[r].first, x}); }
      rem.swap(next);
    }
    for (auto& f : fixed) axes.push_back(f);
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) fr.R[(size_t)a * d + b] = axes[a][b];
    // angle classes, with rates snapped to make coincidences exact
    const int np = fr.nplanes(); fr.cls.assign(np, 2);
    for (int i = 0; i < np; i++) { const double k = std::round(fr.rate[i] / N), h = std::round(fr.rate[i] / N - 0.5) + 0.5;
      if (std::fabs(fr.rate[i] - k * N) < snap) { fr.rate[i] = k * N; fr.cls[i] = 0; } else if (std::fabs(fr.rate[i] - h * N) < snap) { fr.rate[i] = h * N; fr.cls[i] = 1; } }
    for (int i = 0; i < np; i++) for (int j = i + 1; j < np; j++) if (fr.cls[i] == 2 && fr.cls[j] == 2) {
      const double ks = std::round((fr.rate[j] - fr.rate[i]) / N), ko = std::round((fr.rate[j] + fr.rate[i]) / N);
      if (std::fabs(fr.rate[j] - fr.rate[i] - ks * N) < snap) fr.rate[j] = fr.rate[i] + ks * N;
      else if (std::fabs(fr.rate[j] + fr.rate[i] - ko * N) < snap) fr.rate[j] = ko * N - fr.rate[i]; }
  }
  // gauge generators in frame coordinates
  const int np = fr.nplanes(); std::vector<int> F, P;                             // fixed axes, π axes
  for (int i = 0; i < np; i++) { if (fr.cls[i] == 0) { F.push_back(2 * i); F.push_back(2 * i + 1); } if (fr.cls[i] == 1) { P.push_back(2 * i); P.push_back(2 * i + 1); } }
  for (int a = 2 * np; a < d; a++) F.push_back(a);
  auto skew = [&](int a, int b) { std::vector<double> x((size_t)d * d, 0.0); x[(size_t)a * d + b] = 1.0; x[(size_t)b * d + a] = -1.0; return x; };
  for (int a : F) { std::vector<double> c(d, 0.0); c[a] = 1.0; fr.trans.push_back(c); }
  for (const std::vector<int>* set : {&F, &P}) for (size_t i = 0; i < set->size(); i++) for (size_t j = i + 1; j < set->size(); j++) fr.rots.push_back(skew((*set)[i], (*set)[j]));
  for (int i = 0; i < np; i++) if (fr.cls[i] == 2) fr.rots.push_back(skew(2 * i, 2 * i + 1));
  for (int i = 0; i < np; i++) for (int j = i + 1; j < np; j++) if (fr.cls[i] == 2 && fr.cls[j] == 2) {
    const int ui = 2 * i, vi = 2 * i + 1, uj = 2 * j, vj = 2 * j + 1;
    const bool same = std::fmod(fr.rate[j] - fr.rate[i], (double)N) == 0.0, opp = std::fmod(fr.rate[j] + fr.rate[i], (double)N) == 0.0;
    // the block [[0, X], [−Xᵀ, 0]] commutes with R(θ) ⊕ R(±θ) for X ∈ {I, J} (same angle) or the two reflections (opposite)
    auto blk = [&](double x00, double x01, double x10, double x11) { std::vector<double> g((size_t)d * d, 0.0);
      g[(size_t)ui * d + uj] = x00; g[(size_t)ui * d + vj] = x01; g[(size_t)vi * d + uj] = x10; g[(size_t)vi * d + vj] = x11;
      for (int a : {ui, vi}) for (int b : {uj, vj}) g[(size_t)b * d + a] = -g[(size_t)a * d + b]; return g; };
    if (same) { fr.rots.push_back(blk(1, 0, 0, 1)); fr.rots.push_back(blk(0, -1, 1, 0)); }
    else if (opp) { fr.rots.push_back(blk(1, 0, 0, -1)); fr.rots.push_back(blk(0, 1, 1, 0)); } }
  return fr;
}

// derivative rows of the conserved quantities at state Y: momentum along each c, energy, angular momentum
// L_ξ = Σ_j ⟨ξ q_j, v_j⟩ for each ξ; rows × 2Nd
template <class T>
void conserved_jac(int N, int d, double alpha, const T* Y, const std::vector<std::vector<double>>& trans, const std::vector<std::vector<double>>& rots, std::vector<T>& J) {
  const int nd = N * d, n2 = 2 * nd, rows = (int)trans.size() + 1 + (int)rots.size(); J.assign((size_t)rows * n2, T(0));
  int r = 0;
  for (auto& c : trans) { for (int k = 0; k < N; k++) for (int a = 0; a < d; a++) set_d(J[(size_t)r * n2 + nd + k * d + a], c[a]); r++; }
  { T* E = &J[(size_t)r * n2]; T df, s, w, tmp; const double aexp = -(alpha + 2.0) / 2.0;
    for (int i = 0; i < nd; i++) set(E[nd + i], Y[nd + i]);
    for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++) { set_zero(s);
      for (int a = 0; a < d; a++) { sub(df, Y[i * d + a], Y[l * d + a]); fma_add(s, df, df); }
      pow_d(w, s, aexp); mul_d(w, w, alpha);
      for (int a = 0; a < d; a++) { sub(df, Y[i * d + a], Y[l * d + a]); mul(tmp, w, df); add_inplace(E[i * d + a], tmp); sub_inplace(E[l * d + a], tmp); } }
    r++; }
  for (auto& xi : rots) { T* L = &J[(size_t)r * n2]; T tmp;                          // ∂/∂q = ξᵀv = −ξv, ∂/∂v = ξq
    for (int k = 0; k < N; k++) for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { const double x = xi[(size_t)a * d + b]; if (x == 0.0) continue;
      mul_d(tmp, Y[nd + k * d + a], -x); add_inplace(L[k * d + b], tmp); mul_d(tmp, Y[k * d + b], x); add_inplace(L[nd + k * d + a], tmp); }
    r++; }
}

// the MPFR shooting Newton of `refine`, as a function: a double Newton first, then to `digits`
inline double refine_state(int N, int d, double alpha, const std::vector<double>& Zd, const Frame& fr, int digits, int threads, std::vector<mpreal>& Z, bool verbose = false) {
  const int nd = N * d, n2 = 2 * nd;
  NBody<double> nbd(N, d, alpha, 22); std::vector<double> Zc(Zd), Gd; ShootWork<double> Wd;
  std::vector<mpreal> Gm; if (fr.rotating()) { fr.G_mpreal(Gm); Gd.resize(Gm.size()); for (size_t i = 0; i < Gm.size(); i++) Gd[i] = to_double(Gm[i]); }
  double sc = 1.0; for (double v : Zc) sc = std::max(sc, std::fabs(v));
  shoot_newton(nbd, Zc, 1e-16, 20, 1e-13 * sc, -18, -40, Wd, false, Gd.empty() ? nullptr : &Gd, threads);
  const mpfr_prec_t bits = mpreal::default_prec(); const int order = (int)(1.15 * digits) + 6; const double itol = std::pow(10.0, -(digits + 4));
  Z.assign(n2, mpreal(0)); for (int i = 0; i < n2; i++) Z[i] = Zc[i];
  NBody<mpreal> nb(N, d, alpha, order); ShootWork<mpreal> W;
  return shoot_newton(nb, Z, itol, 25, std::pow(10.0, -digits), -(long)(bits / 3), -(long)std::min<mpfr_prec_t>(bits / 2, 96), W, verbose, Gm.empty() ? nullptr : &Gm, threads, 0.9, 24);
}

struct Proof {
  bool ok = false; std::string why;
  int dp = 0, m = 0, k = 0; double radius = 0, newton = 0, kappa = 0, closure = 0, seconds = 0; long steps = 0, rejects = 0; double maxwid = 0, hmin = 0, hmax = 0;
  ival energy, action;
};

// approximate inverse of an n×n matrix given as doubles (columns solved one by one)
inline bool inverse_d(int n, const std::vector<double>& A, std::vector<double>& Y) {
  Y.assign((size_t)n * n, 0.0); std::vector<double> b(n);
  for (int c = 0; c < n; c++) { std::fill(b.begin(), b.end(), 0.0); b[c] = 1.0; if (!la::lu_solve(n, A, b)) return false; for (int r = 0; r < n; r++) Y[(size_t)r * n + c] = b[r]; }
  return true;
}

// Z0: the refined state, 2Nd mpreals with |F(Z0)| far below the box radius. Returns the certificate.
inline Proof prove_state(int N, int d, double alpha, const std::vector<mpreal>& Z0, const Frame& fr, double radius, double tol, int order, int threads, bool verbose = false) {
  Proof pr; auto t0 = std::chrono::steady_clock::now();
  const int nd = N * d, n2 = 2 * nd; pr.dp = d; pr.radius = radius;
  std::vector<double> z(n2); for (int i = 0; i < n2; i++) z[i] = to_double(Z0[i]);
  // --- generators of the gauge group at Z0: translations G fixes, the time shift, the rotations commuting
  // with G that move the orbit
  std::vector<std::vector<double>> trans = fr.trans, rots; std::vector<std::vector<double>> gen;
  for (auto& c : trans) { std::vector<double> g(n2, 0.0); for (int k = 0; k < N; k++) for (int a = 0; a < d; a++) g[k * d + a] = c[a]; gen.push_back(g); }
  { NBody<double> nbd(N, d, alpha, 1); nbd.series(z.data(), z.data() + nd); std::vector<double> g(n2);
    for (int i = 0; i < nd; i++) { g[i] = nbd.V[i]; g[nd + i] = nbd.V[(size_t)nd + i]; } gen.push_back(g); }
  double scale = 0; for (double v : z) scale = std::max(scale, std::fabs(v));
  for (auto& xi : fr.rots) { std::vector<double> g(n2, 0.0); double nrm = 0;
    for (int h = 0; h < 2; h++) for (int k = 0; k < N; k++) for (int a = 0; a < d; a++) { double t = 0; for (int b = 0; b < d; b++) t += xi[(size_t)a * d + b] * z[h * nd + k * d + b]; g[h * nd + k * d + a] = t; nrm += t * t; }
    if (std::sqrt(nrm) > 1e-9 * scale) { gen.push_back(g); rots.push_back(xi); } }
  // orthonormalise; a generator that is a combination of the others (a symmetric orbit) is dropped with its law
  std::vector<std::vector<double>> T;
  { size_t gi = 0; std::vector<int> keep;
    for (auto& g : gen) { std::vector<double> u = g; double n0 = 0; for (double v : u) n0 += v * v;
      for (auto& t : T) { double dot = 0; for (int i = 0; i < n2; i++) dot += t[i] * u[i]; for (int i = 0; i < n2; i++) u[i] -= dot * t[i]; }
      double n1 = 0; for (double v : u) n1 += v * v;
      if (n1 > 1e-16 * n0) { for (double& v : u) v /= std::sqrt(n1); T.push_back(u); keep.push_back((int)gi); } gi++; }
    std::vector<std::vector<double>> tr2, ro2; const int nt = (int)trans.size();
    for (int i : keep) { if (i < nt) tr2.push_back(trans[i]); else if (i > nt) ro2.push_back(rots[i - nt - 1]); }
    trans = tr2; rots = ro2; }
  const int k = (int)T.size(), m = n2 - k; pr.k = k; pr.m = m;
  if ((int)trans.size() + 1 + (int)rots.size() != k) { pr.why = "generators and conservation laws do not pair"; return pr; }
  // the slice: eigenvectors of I − TTᵀ with eigenvalue 1
  std::vector<double> Q((size_t)n2 * m), Pm((size_t)n2 * n2, 0.0), ev;
  for (int i = 0; i < n2; i++) { Pm[(size_t)i * n2 + i] = 1.0; for (auto& t : T) for (int j = 0; j < n2; j++) Pm[(size_t)i * n2 + j] -= t[i] * t[j]; }
  la::sym_eig(n2, Pm, ev);
  { int c = 0; for (int j = 0; j < n2; j++) if (ev[j] > 0.5) { if (c < m) for (int i = 0; i < n2; i++) Q[(size_t)i * m + c] = Pm[(size_t)i * n2 + j]; c++; }
    if (c != m) { pr.why = "slice dimension mismatch"; return pr; } }
  auto shifted = [&](int i) { int h = i / nd, j = (i % nd) / d, a = i % d; return h * nd + ((j + 1) % N) * d + a; };
  // (G S Z)_i as an interval combination of Z: G S Z_{j,a} = Σ_b G_ab Z_{j+1,b}
  std::vector<ival> Gi; if (fr.rotating()) fr.G_ival(Gi);
  auto gs = [&](const std::vector<ival>& Zv, int i, ival& out) { if (Gi.empty()) { set(out, Zv[shifted(i)]); return; }
    const int h = i / nd, j = (i % nd) / d, a = i % d; set_zero(out); for (int b = 0; b < d; b++) fma_add(out, Gi[(size_t)a * d + b], Zv[h * nd + ((j + 1) % N) * d + b]); };
  auto gs_d = [&](const std::vector<double>& Zv, int i) { if (Gi.empty()) return Zv[shifted(i)];
    const int h = i / nd, j = (i % nd) / d, a = i % d; double t = 0; for (int b = 0; b < d; b++) t += Gi[(size_t)a * d + b].mid() * Zv[h * nd + ((j + 1) % N) * d + b]; return t; };
  // equations to drop: the columns on which the conserved quantities are most sensitive, greedily
  std::vector<int> D; std::vector<char> dropped(n2, 0);
  { std::vector<double> SZ(n2); for (int i = 0; i < n2; i++) SZ[i] = gs_d(z, i);
    std::vector<double> Jq; conserved_jac<double>(N, d, alpha, SZ.data(), trans, rots, Jq);
    if (verbose) { std::printf("  conservation rows (%zu translations, energy, %zu rotations), max entries:", trans.size(), rots.size());
      for (int r = 0; r < k; r++) { double mx = 0; for (int c = 0; c < n2; c++) mx = std::max(mx, std::fabs(Jq[(size_t)r * n2 + c])); std::printf(" %.3g", mx); } std::printf("\n"); }
    for (int r = 0; r < k; r++) { int best = -1; double bv = 0;
      for (int c = 0; c < n2; c++) if (!dropped[c] && std::fabs(Jq[(size_t)r * n2 + c]) > bv) { bv = std::fabs(Jq[(size_t)r * n2 + c]); best = c; }
      if (best < 0) { double mx = 0; for (int c = 0; c < n2; c++) mx = std::max(mx, std::fabs(Jq[(size_t)r * n2 + c]));
        pr.why = "conservation row " + std::to_string(r) + " of " + std::to_string(k) + " dependent (max entry " + std::to_string(mx) + ", " + std::to_string(trans.size()) + " translations)"; return pr; }
      D.push_back(best); dropped[best] = 1;
      for (int r2 = 0; r2 < k; r2++) if (r2 != r) { double f = Jq[(size_t)r2 * n2 + best] / Jq[(size_t)r * n2 + best]; for (int c = 0; c < n2; c++) Jq[(size_t)r2 * n2 + c] -= f * Jq[(size_t)r * n2 + c]; } } }
  std::vector<int> keepEq; for (int c = 0; c < n2; c++) if (!dropped[c]) keepEq.push_back(c);
  const ival tend = ival::pi() * 2.0 / ival(N);
  // --- point run: F(Z0)
  std::vector<ival> F0(n2);
  { Verified V(N, d, alpha, order, 0, 1, tol); V.debug = verbose;
    for (int i = 0; i < n2; i++) V.Z[i] = ival(Z0[i].v);
    if (!V.integrate(tend)) { pr.why = "point integration failed"; return pr; }
    std::vector<ival> Zi(n2); for (int i = 0; i < n2; i++) Zi[i] = ival(Z0[i].v);
    ival g; for (int i = 0; i < n2; i++) { gs(Zi, i, g); sub(F0[i], V.Z[i], g); }
    if (verbose) std::printf("  point run: %ld steps, %ld rejected, h in [%.2e, %.2e], |F(Z0)| <= %.2e\n", V.steps, V.rejects, V.hmin, V.hmax, [&] { double mx = 0; for (auto& f : F0) mx = std::max(mx, f.mag()); return mx; }()); }
  // --- box run: Φ and DΦ·Q over B = Z0 + Q[−r, r]^m, then Krawczyk: K = −Y F0 + (I − Y[J]) [−r, r]^m ⊂ int [−r, r]^m.
  // The width of [J] is linear in the radius, so a failed contraction with a Newton term far below the radius is
  // retried on a smaller box; only a Newton term too close to the radius needs more digits.
  std::vector<ival> B(n2), J((size_t)m * m); std::vector<double> Jm((size_t)m * m), Y;
  std::unique_ptr<Verified> Vp;
  for (int attempt = 0; attempt < 4; attempt++) {
    Vp = std::make_unique<Verified>(N, d, alpha, order, m, threads, tol); Verified& V = *Vp; V.debug = verbose;
    for (int i = 0; i < n2; i++) { double w = 0; for (int c = 0; c < m; c++) w += std::fabs(Q[(size_t)i * m + c]); V.Z[i] = ival(Z0[i].v); V.Z[i].inflate(0.0, radius * w); }
    for (int c = 0; c < m; c++) for (int i = 0; i < n2; i++) set_d(V.Psi[(size_t)c * n2 + i], Q[(size_t)i * m + c]);
    B = V.Z;
    if (!V.integrate(tend)) { pr.why = "box integration failed"; return pr; }
    pr.steps = V.steps; pr.rejects = V.rejects; pr.maxwid = V.maxwid; pr.hmin = V.hmin; pr.hmax = V.hmax; pr.radius = radius;
    if (verbose) std::printf("  box run:   %ld steps, %ld rejected, h in [%.2e, %.2e], max state width %.2e (radius %.1e)\n", V.steps, V.rejects, V.hmin, V.hmax, V.maxwid, radius);
    { std::vector<ival> Qc(n2); ival g;
      for (int c = 0; c < m; c++) { for (int i = 0; i < n2; i++) set_d(Qc[i], Q[(size_t)i * m + c]);
        for (int r = 0; r < m; r++) { const int i = keepEq[r]; ival& e = J[(size_t)r * m + c]; gs(Qc, i, g); sub(e, V.Psi[(size_t)c * n2 + i], g); Jm[(size_t)r * m + c] = e.mid(); } } }
    if (!inverse_d(m, Jm, Y)) { pr.why = "reduced Jacobian singular"; return pr; }
    double kmax = 0, kappa = 0, nmax = 0; ival g, acc;
    for (int r = 0; r < m; r++) {
      set_zero(acc); for (int c = 0; c < m; c++) fma_add(acc, ival(Y[(size_t)r * m + c]), F0[keepEq[c]]);
      nmax = std::max(nmax, acc.mag());
      double rowsum = 0;
      for (int c = 0; c < m; c++) { set_zero(g); if (r == c) set_d(g, 1.0); for (int l = 0; l < m; l++) fma_sub(g, ival(Y[(size_t)r * m + l]), J[(size_t)l * m + c]); rowsum += g.mag(); }
      kappa = std::max(kappa, rowsum); kmax = std::max(kmax, acc.mag() + rowsum * radius); }
    pr.newton = nmax; pr.kappa = kappa;
    if (kmax < radius && kappa < 1.0) break;
    const double rnew = std::max(radius / (20 * std::max(kappa, 1.0)), 1e3 * nmax);
    if (verbose) std::printf("  Krawczyk failed (|Y F| %.1e, contraction %.2e at radius %.1e)%s\n", nmax, kappa, radius, rnew < 0.5 * radius ? ", retrying" : "");
    if (!(rnew < 0.5 * radius) || attempt == 3) { pr.why = nmax * 1e3 > radius ? "Krawczyk inclusion failed: refine with more digits" : "Krawczyk inclusion failed"; return pr; }
    radius = rnew;
  }
  Verified& V = *Vp;
  // closure: the dropped components agree because the conserved quantities do; ∂(P,E,L)/∂Y_D must be
  // nonsingular on the hull of Φ(B) and S B
  { std::vector<ival> H(n2), Jc; ival gh; for (int i = 0; i < n2; i++) { H[i] = V.Z[i]; gs(B, i, gh); H[i].hull(gh); }
    conserved_jac<ival>(N, d, alpha, H.data(), trans, rots, Jc);
    std::vector<double> Mm((size_t)k * k), Ym; std::vector<ival> M((size_t)k * k);
    for (int r = 0; r < k; r++) for (int c = 0; c < k; c++) { M[(size_t)r * k + c] = Jc[(size_t)r * n2 + D[c]]; Mm[(size_t)r * k + c] = M[(size_t)r * k + c].mid(); }
    if (!inverse_d(k, Mm, Ym)) { pr.why = "closure matrix singular"; return pr; }
    double cl = 0; ival g;
    for (int r = 0; r < k; r++) { double rowsum = 0;
      for (int c = 0; c < k; c++) { set_zero(g); if (r == c) set_d(g, 1.0); for (int l = 0; l < k; l++) fma_sub(g, ival(Ym[(size_t)r * k + l]), M[(size_t)l * k + c]); rowsum += g.mag(); }
      cl = std::max(cl, rowsum); }
    pr.closure = cl; if (!(cl < 1.0)) { pr.why = "closure matrix not verifiably nonsingular"; return pr; } }
  // the invariants of the proven orbit: energy on the box, the action from the run (N bodies over T/N)
  { ival E(0), r2, df, w;
    for (int i = 0; i < nd; i++) fma_add(E, B[nd + i], B[nd + i]); mul_d(E, E, 0.5);
    for (int i = 0; i < N; i++) for (int l = i + 1; l < N; l++) { set_zero(r2); for (int a = 0; a < d; a++) { sub(df, B[i * d + a], B[l * d + a]); fma_add(r2, df, df); }
      pow_d(w, r2, -alpha / 2); sub_inplace(E, w); }
    pr.energy = E; pr.action = V.A; }
  pr.ok = true; pr.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return pr;
}
#endif
