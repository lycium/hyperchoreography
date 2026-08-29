// Pseudo-arclength continuation in the potential exponent α, with branch switching at Morse index changes.
#pragma once
#include "search.hpp"
#include <functional>

struct ContCfg { double alpha_lo = 0.6, alpha_hi = 2.4, h0 = 0.02, hmax = 0.1, hmin = 1e-4, kick = 0.05; int depth = 2, max_branches = 200, max_steps = 1500; bool verbose = true; };

struct Continuer {
  const Config& cfg; const ContCfg& cc; Ctx& ctx; const Problem& P0;
  std::function<void(std::vector<double>&, const std::string&)> on_solution;
  int branches = 0; const int n; std::vector<double> H, w, Gv, J, JtJ, rhs, ga;
  Continuer(const Config& c, const ContCfg& cc_, Ctx& ctx_, const Problem& P) : cfg(c), cc(cc_), ctx(ctx_), P0(P), n(P.n) {}

  // orthonormal gauge generators (time shift, rotations) as columns of Gv
  int gauge(const std::vector<double>& x) {
    std::vector<std::vector<double>> gens; std::vector<double> xd; deriv_coeffs(P0, x.data(), xd); gens.push_back(xd);
    for (int a = 0; a < P0.d; a++) for (int b = a + 1; b < P0.d; b++) { std::vector<double> v(n, 0.0); for (int i = 0; i < P0.nb; i++) { v[i * P0.d + a] = -x[i * P0.d + b]; v[i * P0.d + b] = x[i * P0.d + a]; } gens.push_back(v); }
    std::vector<std::vector<double>> on; double xn = 0; for (double v : x) xn += v * v; xn = std::sqrt(xn);
    for (auto& g : gens) { for (auto& u : on) { double dp = 0; for (int i = 0; i < n; i++) dp += u[i] * g[i]; for (int i = 0; i < n; i++) g[i] -= dp * u[i]; }
      double nn = 0; for (double v : g) nn += v * v; nn = std::sqrt(nn); if (nn < 1e-8 * xn) continue; for (double& v : g) v /= nn; on.push_back(g); }
    int ng = (int)on.size(); Gv.assign((size_t)n * ng, 0.0); for (int k = 0; k < ng; k++) for (int i = 0; i < n; i++) Gv[(size_t)i * ng + k] = on[k][i];
    return ng;
  }
  bool hess_at(Problem& Pa, double alpha, const std::vector<double>& x) { Pa.alpha = alpha; return hessian(Pa, x.data(), H, ctx.w); }
  int index_at(Problem& Pa, double alpha, const std::vector<double>& x, int* nullity = nullptr) {
    if (!hess_at(Pa, alpha, x)) return -1; Inertia I = inertia(n, H); if (nullity) *nullity = I.zero; return I.neg;
  }
  // least squares for the bordered system [H g_α; Gᵀ 0; tᵀ] (δx, δα) = −r
  bool bordered_solve(int ng, const std::vector<double>& t, const std::vector<double>& r_top, const std::vector<double>& r_g, double r_t, std::vector<double>& sol, double mu) {
    const int m = n + 1, rows = n + ng + 1; J.assign((size_t)rows * m, 0.0);
    for (int i = 0; i < n; i++) { for (int k = 0; k < n; k++) J[(size_t)i * m + k] = H[(size_t)i * n + k]; J[(size_t)i * m + n] = ga[i]; }
    for (int g = 0; g < ng; g++) for (int i = 0; i < n; i++) J[(size_t)(n + g) * m + i] = Gv[(size_t)i * ng + g];
    for (int k = 0; k < m; k++) J[(size_t)(n + ng) * m + k] = t[k];
    JtJ.assign((size_t)m * m, 0.0); rhs.assign(m, 0.0); double dmax = 0;
    for (int i = 0; i < m; i++) for (int j = i; j < m; j++) { double s = 0; for (int r = 0; r < rows; r++) s += J[(size_t)r * m + i] * J[(size_t)r * m + j]; JtJ[(size_t)i * m + j] = JtJ[(size_t)j * m + i] = s; if (i == j) dmax = std::max(dmax, s); }
    for (int i = 0; i < m; i++) { JtJ[(size_t)i * m + i] += mu * dmax; double s = 0;
      for (int r = 0; r < n; r++) s += J[(size_t)r * m + i] * r_top[r]; for (int g = 0; g < ng; g++) s += J[(size_t)(n + g) * m + i] * r_g[g]; s += J[(size_t)(n + ng) * m + i] * r_t; rhs[i] = -s; }
    if (!la::lu_solve(m, JtJ, rhs)) return false; sol = rhs; return true;
  }
  // unit branch tangent: kernel of [H g_α] ⟂ gauge, oriented along tprev (or +α)
  bool tangent(Problem& Pa, const std::vector<double>& x, double alpha, const std::vector<double>* tprev, std::vector<double>& t) {
    if (!hess_at(Pa, alpha, x)) return false; ga.resize(n); if (!grad_alpha(Pa, x.data(), ga.data(), ctx.w)) return false;
    int ng = gauge(x); std::vector<double> tref(n + 1, 0.0); if (tprev) tref = *tprev; else tref[n] = 1.0;
    std::vector<double> z(n, 0.0), zg(ng, 0.0), sol;
    if (!bordered_solve(ng, tref, z, zg, -1.0, sol, 1e-12)) return false;      // tref·t = 1
    double nn = 0; for (double v : sol) nn += v * v; nn = std::sqrt(nn); if (nn == 0) return false;
    t.resize(n + 1); for (int k = 0; k <= n; k++) t[k] = sol[k] / nn;
    return true;
  }
  // returns iterations or -1
  int corrector(Problem& Pa, std::vector<double>& x, double& alpha, const std::vector<double>& t, const std::vector<double>& xp, double ap, int max_it = 20) {
    std::vector<double> g(n), sol, rg; double mu = 1e-12;
    for (int it = 0; it < max_it; it++) {
      Pa.alpha = alpha; double A = action_grad(Pa, x.data(), g.data(), ctx.w); if (!std::isfinite(A)) return -1;
      double gn = 0; for (double v : g) gn += v * v; gn = std::sqrt(gn);
      if (gn < cfg.gtol) return it;
      if (!hess_at(Pa, alpha, x)) return -1; ga.resize(n); grad_alpha(Pa, x.data(), ga.data(), ctx.w);
      int ng = gauge(x); rg.assign(ng, 0.0); for (int k = 0; k < ng; k++) for (int i = 0; i < n; i++) rg[k] += Gv[(size_t)i * ng + k] * (x[i] - xp[i]);
      double rt = 0; for (int i = 0; i < n; i++) rt += t[i] * (x[i] - xp[i]); rt += t[n] * (alpha - ap);
      if (!bordered_solve(ng, t, g, rg, rt, sol, mu)) return -1;
      for (int i = 0; i < n; i++) x[i] += sol[i]; alpha += sol[n];
      if (alpha < 0.05 || alpha > 4.0) return -1;
    }
    Pa.alpha = alpha; action_grad(Pa, x.data(), g.data(), ctx.w); double gn = 0; for (double v : g) gn += v * v;
    return std::sqrt(gn) < 100 * cfg.gtol ? max_it : -1;
  }
  bool polish_at(Problem& Pa, double alpha, std::vector<double>& x, int iters = 40) { Pa.alpha = alpha; Symmetry none; return polish(Pa, none, ctx, x, iters, cfg.gtol); }

  void report(std::vector<double> x, double a_prev, double a_new, const std::string& tag) {
    Problem Pa = P0; if (!polish_at(Pa, 1.0, x)) { if (cc.verbose) std::printf("  [%s] crossing of α=1 near %.4f→%.4f lost in polish\n", tag.c_str(), a_prev, a_new); return; }
    std::string why; Record rec; const Problem* Pp = &P0; std::vector<double> xs = x;
    if (certify(cfg, Pp, xs, Symmetry(), ctx, rec, why)) on_solution(x, tag); else if (cc.verbose) std::printf("  [%s] α=1 crossing not certified: %s\n", tag.c_str(), why.c_str());
  }

  // tref: initial direction (roots ±e_α, switched branches ±crossing eigenvector)
  void follow(std::vector<double> x, double alpha, const std::vector<double>& tref, int depth, const std::string& tag) {
    if (branches++ > cc.max_branches) return;
    Problem Pa = P0; if (!polish_at(Pa, alpha, x)) return;
    std::vector<double> t; if (!tangent(Pa, x, alpha, &tref, t)) return;
    double dp = 0; for (int k = 0; k <= n; k++) dp += t[k] * tref[k]; if (dp < 0) for (double& v : t) v = -v;
    int idx = index_at(Pa, alpha, x); double h = cc.h0; const std::vector<double> x0 = x; const double a0 = alpha;
    if (cc.verbose) std::printf("  [%s] start α=%.4f morse=%d\n", tag.c_str(), alpha, idx);
    for (int step = 0; step < cc.max_steps; step++) {
      std::vector<double> xp(n), xn; double ap = alpha + h * t[n]; for (int i = 0; i < n; i++) xp[i] = x[i] + h * t[i];
      xn = xp; double an = ap; int its = corrector(Pa, xn, an, t, xp, ap);
      if (its < 0) { h *= 0.5; if (h < cc.hmin) { if (cc.verbose) std::printf("  [%s] stuck at α=%.4f\n", tag.c_str(), alpha); return; } continue; }
      std::vector<double> tn; if (!tangent(Pa, xn, an, &t, tn)) return;
      int idn = index_at(Pa, an, xn);
      if (idn != idx && idn >= 0 && depth < cc.depth) branch_from(Pa, x, alpha, xn, an, idx, idn, depth, tag);
      if ((alpha - 1.0) * (an - 1.0) <= 0 && step > 0) report(std::fabs(alpha - 1) < std::fabs(an - 1) ? x : xn, alpha, an, tag);
      x.swap(xn); alpha = an; t.swap(tn); idx = idn;
      h = its <= 3 ? std::min(h * 1.5, cc.hmax) : its >= 8 ? h * 0.6 : h;
      if (alpha < cc.alpha_lo || alpha > cc.alpha_hi) break;
      if (step > 20 && std::fabs(alpha - a0) < 0.5 * h && loop_distance(P0, x.data(), P0, x0.data()) < 1e-4) { if (cc.verbose) std::printf("  [%s] closed loop\n", tag.c_str()); break; }
    }
    if (cc.verbose) std::printf("  [%s] end α=%.4f morse=%d\n", tag.c_str(), alpha, idx);
  }
  // an eigenvalue crossed zero between (xa, aa) and (xb, ab): bisect, kick along its eigenvector, follow new branches
  void branch_from(Problem& Pa, const std::vector<double>& xa, double aa, const std::vector<double>& xb, double ab, int ia, int ib, int depth, const std::string& tag) {
    double lo = aa, hi = ab; std::vector<double> xlo = xa, xhi = xb;
    for (int b = 0; b < 10; b++) { double mid = 0.5 * (lo + hi); std::vector<double> xm = xlo; if (!polish_at(Pa, mid, xm)) break; int im = index_at(Pa, mid, xm); if (im == ia) { lo = mid; xlo = xm; } else { hi = mid; xhi = xm; } }
    hess_at(Pa, hi, xhi); std::vector<double> Hc = H, wv; la::sym_eig(n, Hc, wv);
    double lmax = 0; for (double v : wv) lmax = std::max(lmax, std::fabs(v)); int kc = -1; double best = INF;
    for (int k = 0; k < n; k++) if (std::fabs(wv[k]) > 1e-8 * lmax && std::fabs(wv[k]) < best) { best = std::fabs(wv[k]); kc = k; }
    if (kc < 0) return;
    double astar = 0.5 * (lo + hi);
    if (cc.verbose) std::printf("  [%s] bifurcation at α≈%.5f (morse %d→%d), branching\n", tag.c_str(), astar, ia, ib);
    // arclength constraint along v* keeps the corrector off the parent
    double rms = 0; for (double v : xhi) rms += 0.5 * v * v; rms = std::sqrt(rms); int nb = 0;
    std::vector<double> t(n + 1, 0.0); for (int i = 0; i < n; i++) t[i] = Hc[(size_t)i * n + kc];
    for (double eps : {1.0, 2.0, 4.0, 8.0}) { bool got = false;
      for (double sgn : {1.0, -1.0}) {
        double h = sgn * eps * cc.kick * rms; std::vector<double> xp = xhi, xn; for (int i = 0; i < n; i++) xp[i] += h * t[i];
        xn = xp; double an = astar; std::vector<double> ts = t; for (double& v : ts) v *= sgn;
        if (corrector(Pa, xn, an, ts, xp, astar, 20) < 0) continue;
        if (an < cc.alpha_lo || an > cc.alpha_hi) continue;
        std::vector<double> parent = xhi; if (!polish_at(Pa, an, parent) || loop_distance(P0, xn.data(), P0, parent.data()) < 1e-3) continue;
        nb++; got = true; std::string ct = tag + "." + std::to_string(nb);
        std::vector<double> tm = ts; for (double& v : tm) v = -v;
        follow(xn, an, ts, depth + 1, ct + "+"); follow(xn, an, tm, depth + 1, ct + "-");
      }
      if (got) break;
    }
  }
  void explore(const std::vector<double>& x, const std::string& tag) {
    std::vector<double> ea(n + 1, 0.0); ea[n] = 1.0; follow(x, 1.0, ea, 0, tag + "+"); ea[n] = -1.0; follow(x, 1.0, ea, 0, tag + "-"); }
};
