// L-BFGS and Levenberg–Marquardt Newton on ∇f = 0 in the Hessian eigenbasis (gauge null space gets no step).
#pragma once
#include "linalg.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <limits>

using FnGrad = std::function<double(const double*, double*)>;
using FnHess = std::function<bool(const double*, std::vector<double>&)>;

struct OptResult { int iters = 0; double f = 0, gnorm = 0; bool converged = false; };

inline double vnorm(const std::vector<double>& v) { double s = 0; for (double a : v) s += a * a; return std::sqrt(s); }
inline double vdot(const std::vector<double>& a, const std::vector<double>& b) { double s = 0; for (size_t i = 0; i < a.size(); i++) s += a[i] * b[i]; return s; }

inline OptResult lbfgs(int n, std::vector<double>& x, const FnGrad& fn, int max_iter, double gtol, int mem = 10) {
  std::vector<double> g(n), gnew(n), xnew(n), dir(n), al(mem), rho(mem);
  std::vector<double> S((size_t)mem * n), Y((size_t)mem * n);      // ring buffers
  int m = 0, head = 0;
  OptResult R;
  double f = fn(x.data(), g.data());
  if (!std::isfinite(f)) { R.f = f; R.gnorm = INFINITY; return R; }
  double gn = vnorm(g);
  auto Si = [&](int i) { return &S[(size_t)((head - m + i + mem) % mem) * n]; };
  auto Yi = [&](int i) { return &Y[(size_t)((head - m + i + mem) % mem) * n]; };
  auto rhoi = [&](int i) -> double& { return rho[(head - m + i + mem) % mem]; };
  auto dotp = [&](const double* a, const double* b) { double s = 0; for (int k = 0; k < n; k++) s += a[k] * b[k]; return s; };
  int it = 0;
  for (; it < max_iter && gn > gtol; it++) {
    dir = g;
    for (int i = m - 1; i >= 0; i--) { double a = rhoi(i) * dotp(Si(i), dir.data()); al[i] = a; const double* y = Yi(i); for (int k = 0; k < n; k++) dir[k] -= a * y[k]; }
    double h0 = m ? dotp(Si(m - 1), Yi(m - 1)) / dotp(Yi(m - 1), Yi(m - 1)) : 1.0 / std::max(gn, 1e-30);
    for (double& v : dir) v *= h0;
    for (int i = 0; i < m; i++) { double b = rhoi(i) * dotp(Yi(i), dir.data()); const double* s = Si(i); for (int k = 0; k < n; k++) dir[k] += (al[i] - b) * s[k]; }
    for (double& v : dir) v = -v;
    double gd = vdot(g, dir);
    if (!(gd < 0)) { for (int k = 0; k < n; k++) dir[k] = -g[k] / std::max(gn, 1e-30); gd = -gn; m = 0; }
    double step = 1.0, fnew = INFINITY; bool ok = false;
    for (int ls = 0; ls < 40; ls++) {
      for (int k = 0; k < n; k++) xnew[k] = x[k] + step * dir[k];
      fnew = fn(xnew.data(), gnew.data());
      if (std::isfinite(fnew) && fnew <= f + 1e-4 * step * gd) { ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;
    double* s = &S[(size_t)head * n]; double* y = &Y[(size_t)head * n]; double sy = 0, ss = 0, yy = 0;
    for (int k = 0; k < n; k++) { s[k] = xnew[k] - x[k]; y[k] = gnew[k] - g[k]; sy += s[k] * y[k]; ss += s[k] * s[k]; yy += y[k] * y[k]; }
    if (sy > 1e-12 * std::sqrt(ss * yy)) { rho[head] = 1.0 / sy; head = (head + 1) % mem; m = std::min(m + 1, mem); }
    x.swap(xnew); g.swap(gnew); f = fnew; gn = vnorm(g);
  }
  R.iters = it; R.f = f; R.gnorm = gn; R.converged = gn <= gtol;
  return R;
}

// δ = −Σ_k q_k λ_k/(λ_k²+μ) (q_kᵀg); μ adapted on ||∇f||
inline OptResult newton_lm(int n, std::vector<double>& x, const FnGrad& fn, const FnHess& fnH, int max_iter, double gtol, double* mu_io = nullptr) {
  std::vector<double> g(n), gnew(n), xnew(n), H, w, gq(n), delta(n);
  OptResult R;
  double f = fn(x.data(), g.data());
  if (!std::isfinite(f)) { R.f = f; R.gnorm = INFINITY; return R; }
  double gn = vnorm(g), mu = mu_io ? *mu_io : -1.0;
  int it = 0;
  for (; it < max_iter && gn > gtol; it++) {
    if (!fnH(x.data(), H)) break;
    la::sym_eig(n, H, w);
    double lmax = 0; for (double v : w) lmax = std::max(lmax, std::fabs(v));
    if (mu < 0) mu = 1e-4 * lmax * lmax;
    for (int k = 0; k < n; k++) { double s = 0; for (int i = 0; i < n; i++) s += H[(size_t)i * n + k] * g[i]; gq[k] = s; }
    bool accepted = false;
    for (int tries = 0; tries < 25 && !accepted; tries++) {
      for (int i = 0; i < n; i++) { double s = 0; for (int k = 0; k < n; k++) s += H[(size_t)i * n + k] * (w[k] / (w[k] * w[k] + mu)) * gq[k]; delta[i] = -s; xnew[i] = x[i] + delta[i]; }
      double fnew = fn(xnew.data(), gnew.data());
      double gnn = std::isfinite(fnew) ? vnorm(gnew) : INFINITY;
      if (gnn < gn) { accepted = true; x.swap(xnew); g.swap(gnew); f = fnew; gn = gnn; mu = std::max(mu / 5.0, 1e-18 * lmax * lmax); }
      else mu *= 8.0;
    }
    if (!accepted) break;
  }
  if (mu_io) *mu_io = mu;
  R.iters = it; R.f = f; R.gnorm = gn; R.converged = gn <= gtol;
  return R;
}

// negative / zero / positive eigenvalue counts
struct Inertia { int neg = 0, zero = 0, pos = 0; std::vector<double> eig; };
inline Inertia inertia(int n, std::vector<double> H, double rel_tol = 1e-8) {
  Inertia I; la::sym_eig(n, H, I.eig);
  double lmax = 0; for (double v : I.eig) lmax = std::max(lmax, std::fabs(v));
  for (double v : I.eig) { if (std::fabs(v) <= rel_tol * lmax) I.zero++; else if (v < 0) I.neg++; else I.pos++; }
  return I;
}
