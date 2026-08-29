// Numerical self-checks (make test).
#include "catalog.hpp"
#include "optim.hpp"
#include "taylor.hpp"
#include <cstdio>
#include <chrono>

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; std::printf("  FAIL: " __VA_ARGS__); std::printf("\n"); } else { std::printf("  ok:   " __VA_ARGS__); std::printf("\n"); } } while (0)

static double circle_radius(int N, double alpha = 1.0) {   // N-gon with ω = 1
  double s = 0; for (int k = 1; k < N; k++) { double ck = 2 * std::sin(PI * k / N); s += alpha * std::sin(PI * k / N) / std::pow(ck, alpha + 1); }
  return std::pow(s, 1.0 / (alpha + 2.0));
}

int main() {
  std::printf("[sym_eig]\n");
  { int n = 40; la::Rng rng(1); std::vector<double> A(n * n), A0; for (int i = 0; i < n; i++) for (int j = i; j < n; j++) A[i * n + j] = A[j * n + i] = rng.normal(); A0 = A;
    std::vector<double> w; la::sym_eig(n, A, w); double err = 0, orth = 0;
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) { double s = 0; for (int j = 0; j < n; j++) s += A0[i * n + j] * A[j * n + k]; err = std::max(err, std::fabs(s - w[k] * A[i * n + k])); }
    for (int k = 0; k < n; k++) for (int l = 0; l < n; l++) { double s = 0; for (int i = 0; i < n; i++) s += A[i * n + k] * A[i * n + l]; orth = std::max(orth, std::fabs(s - (k == l))); }
    CHECK(err < 1e-10 && orth < 1e-10, "eigen residual %.2e orthonormality %.2e", err, orth); }

  std::printf("[action derivatives]\n");
  for (int cfg = 0; cfg < 3; cfg++) {
    int N = cfg == 0 ? 3 : cfg == 1 ? 4 : 5, d = cfg == 0 ? 2 : cfg == 1 ? 3 : 4; double alpha = cfg == 2 ? 1.5 : 1.0;
    Problem P; P.init(N, d, 10, 0, alpha); Work w; w.resize(P); la::Rng rng(7 + cfg);
    std::vector<double> x, g(P.n), gp(P.n), gm(P.n), xp, v(P.n), Hv(P.n), H;
    random_guess(P, rng, 4, 1.0, x); for (double& xi : x) xi *= 1.5;
    double kin, pot; action_grad(P, x.data(), nullptr, w, &kin, &pot); double lam = optimal_scale(kin, pot, alpha); for (double& xi : x) xi *= lam;
    double A = action_grad(P, x.data(), g.data(), w); double h = 1e-5, gerr = 0, hverr = 0, herr = 0;
    for (double& vi : v) vi = rng.normal();
    hessvec(P, x.data(), v.data(), Hv.data(), w); hessian(P, x.data(), H, w);
    for (int i = 0; i < P.n; i++) {
      xp = x; xp[i] += h; double fp = action_grad(P, xp.data(), gp.data(), w); xp[i] -= 2 * h; double fm = action_grad(P, xp.data(), gm.data(), w);
      gerr = std::max(gerr, std::fabs((fp - fm) / (2 * h) - g[i]) / (1 + std::fabs(g[i])));
      double hv = 0; for (int k = 0; k < P.n; k++) hv += H[(size_t)i * P.n + k] * v[k]; herr = std::max(herr, std::fabs(hv - Hv[i]) / (1 + std::fabs(Hv[i])));
    }
    xp = x; for (int i = 0; i < P.n; i++) xp[i] += h * v[i]; action_grad(P, xp.data(), gp.data(), w); for (int i = 0; i < P.n; i++) xp[i] -= 2 * h * v[i]; action_grad(P, xp.data(), gm.data(), w);
    for (int i = 0; i < P.n; i++) hverr = std::max(hverr, std::fabs((gp[i] - gm[i]) / (2 * h) - Hv[i]) / (1 + std::fabs(Hv[i])));
    CHECK(gerr < 1e-6 && hverr < 1e-5 && herr < 1e-9, "N=%d d=%d alpha=%.1f  A=%.6f  grad %.1e  Hv %.1e  H-vs-Hv %.1e", N, d, alpha, A, gerr, hverr, herr);
  }

  std::printf("[symmetry basis]\n");
  { Problem P; P.init(3, 2, 12); Symmetry S = Symmetry::parse("t+1/2 s[-1,-2]; t-0 s[1,-2]", 2); int r; auto B = S.basis(P, r);
    Reduced R(P, &B, r); la::Rng rng(3); std::vector<double> y(r), x(P.n), qa, qb; for (double& yi : y) yi = rng.normal(); R.expand(y.data(), x.data());
    sample_curve(P.modes, 2, x.data(), 64, 1, 0, 2, qa); sample_curve(P.modes, 2, x.data(), 64, 1, PI, 2, qb); double e1 = 0;
    for (size_t k = 0; k < qa.size(); k++) e1 = std::max(e1, std::fabs(qa[k] + qb[k]));
    sample_curve(P.modes, 2, x.data(), 64, -1, 0, 2, qb); double e2 = 0;
    for (int j = 0; j < 64; j++) { e2 = std::max(e2, std::fabs(qa[j * 2] - qb[j * 2])); e2 = std::max(e2, std::fabs(qa[j * 2 + 1] + qb[j * 2 + 1])); }
    double orth = 0; for (int c = 0; c < r; c++) for (int c2 = 0; c2 < r; c2++) { double s = 0; for (int i = 0; i < P.n; i++) s += B[(size_t)i * r + c] * B[(size_t)i * r + c2]; orth = std::max(orth, std::fabs(s - (c == c2))); }
    CHECK(e1 < 1e-12 && e2 < 1e-12 && orth < 1e-12 && r == 8, "eight group: r=%d (expect 8), q(t+π)=-q(t) err %.1e, q(-t)=σq(t) err %.1e, orthonormal %.1e", r, e1, e2, orth);
    Work w; w.resize(P); std::vector<double> gy(r), gp(r), yp; R.f(y.data(), gy.data(), w); double err = 0, h = 1e-6;
    for (int c = 0; c < r; c++) { yp = y; yp[c] += h; double fp = R.f(yp.data(), nullptr, w); yp[c] -= 2 * h; double fm = R.f(yp.data(), nullptr, w); err = std::max(err, std::fabs((fp - fm) / (2 * h) - gy[c]) / (1 + std::fabs(gy[c]))); }
    std::vector<double> Hy; R.hess(y.data(), Hy, w); std::vector<double> v(r), yv, g1(r), g2(r); for (double& vi : v) vi = rng.normal(); double herr = 0;
    yv = y; for (int c = 0; c < r; c++) yv[c] += h * v[c]; R.f(yv.data(), g1.data(), w); for (int c = 0; c < r; c++) yv[c] -= 2 * h * v[c]; R.f(yv.data(), g2.data(), w);
    for (int c = 0; c < r; c++) { double hv = 0; for (int c2 = 0; c2 < r; c2++) hv += Hy[(size_t)c * r + c2] * v[c2]; herr = std::max(herr, std::fabs(hv - (g1[c] - g2[c]) / (2 * h)) / (1 + std::fabs(hv))); }
    CHECK(err < 1e-6 && herr < 1e-5, "reduced gradient %.1e, reduced Hessian %.1e", err, herr);
    Symmetry S3 = Symmetry::parse("t+1/3 r(1,2,1/3)", 3); Problem P3; P3.init(4, 3, 9); int r3; auto B3 = S3.basis(P3, r3);
    CHECK(r3 > 0, "3D rotational generator parses, r=%d of %d", r3, P3.n); }

  std::printf("[Taylor integrator]\n");
  for (int N : {3, 5}) { int d = 2; double R = circle_radius(N);
    std::vector<double> pos(N * d), vel(N * d); for (int j = 0; j < N; j++) { double th = 2 * PI * j / N; pos[j * 2] = R * std::cos(th); pos[j * 2 + 1] = R * std::sin(th); vel[j * 2] = -R * std::sin(th); vel[j * 2 + 1] = R * std::cos(th); }
    NBody<double> nb(N, d, 1.0, 22); double E0 = nbody_energy(N, d, 1.0, pos.data(), vel.data());
    auto t0 = std::chrono::steady_clock::now(); double res = chore_residual(nb, pos, vel, 1e-16); double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::vector<double> p2 = pos, v2 = vel; int steps = nb.integrate(p2, v2, 2 * PI, 1e-16); double E1 = nbody_energy(N, d, 1.0, p2.data(), v2.data()); double ret = 0;
    for (int i = 0; i < N * d; i++) ret = std::max(ret, std::max(std::fabs(p2[i] - pos[i]), std::fabs(v2[i] - vel[i])));
    CHECK(res < 1e-12 && ret < 1e-11 && std::fabs(E1 - E0) < 1e-12, "circle N=%d: shift residual %.1e, full-period return %.1e, dE %.1e, %d steps/period, %.2f ms per T/N", N, res, ret, E1 - E0, steps, ms); }
#ifdef HAVE_MPFR
  { mpreal::set_default_prec(400); int N = 3, d = 2; mpreal R = pow(mpreal(1) / sqrt(mpreal(3)), mpreal(1) / 3);   // R³ = 1/√3
    std::vector<mpreal> pos(N * d), vel(N * d); mpreal pi = mpreal::pi();
    for (int j = 0; j < N; j++) { mpreal th = pi * 2 * j / N; pos[j * 2] = R * cos(th); pos[j * 2 + 1] = R * sin(th); vel[j * 2] = -R * sin(th); vel[j * 2 + 1] = R * cos(th); }
    NBody<mpreal> nb(N, d, 1.0, 140); auto t0 = std::chrono::steady_clock::now(); double res = chore_residual(nb, pos, vel, 1e-110);
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    CHECK(res < 1e-100, "MPFR (400 bits) circle: shift residual %.1e, %.0f ms per T/N", res, ms); }
#endif

  std::printf("[loop distance / canonicalisation]\n");
  { Problem P; P.init(3, 3, 8); la::Rng rng(11); std::vector<double> x, y; random_guess(P, rng, 5, 1.0, x);
    std::vector<double> Rm(9); for (double& v : Rm) v = rng.normal();
    for (int c = 0; c < 3; c++) { for (int p = 0; p < c; p++) { double s = 0; for (int i = 0; i < 3; i++) s += Rm[i * 3 + c] * Rm[i * 3 + p]; for (int i = 0; i < 3; i++) Rm[i * 3 + c] -= s * Rm[i * 3 + p]; } double nn = 0; for (int i = 0; i < 3; i++) nn += Rm[i * 3 + c] * Rm[i * 3 + c]; for (int i = 0; i < 3; i++) Rm[i * 3 + c] /= std::sqrt(nn); }
    double tau = 0.7123; y.assign(P.n, 0.0);
    for (int mu = 0; mu < P.nm; mu++) { int m = P.modes[mu]; double c = std::cos(m * tau), s = std::sin(m * tau);
      for (int a = 0; a < 3; a++) { double cc = 0, ss = 0; for (int b = 0; b < 3; b++) { cc += Rm[a * 3 + b] * x[(2 * mu) * 3 + b]; ss += Rm[a * 3 + b] * x[(2 * mu + 1) * 3 + b]; }
        y[(2 * mu) * 3 + a] = c * cc + s * ss; y[(2 * mu + 1) * 3 + a] = -(-s * cc + c * ss); } }
    double dist = loop_distance(P, x.data(), P, y.data()); std::vector<double> z; random_guess(P, rng, 5, 1.0, z); double dist2 = loop_distance(P, x.data(), P, z.data());
    CHECK(dist < 1e-9 && dist2 > 0.05, "equivalent loops distance %.1e, unrelated loops %.3f", dist, dist2);
    std::vector<double> xc = x, sv; int deff = canonical_frame(P.nb, P.d, xc, sv); double dist3 = loop_distance(P, x.data(), P, xc.data());
    CHECK(deff == 3 && dist3 < 1e-9, "canonical frame deff=%d distance to original %.1e", deff, dist3);
    Problem P2; P2.init(3, 3, 16); std::vector<double> x2(P2.n, 0.0);          // 2-fold cover
    for (int mu = 0; mu < P.nm; mu++) { int m2 = 2 * P.modes[mu]; auto it = std::lower_bound(P2.modes.begin(), P2.modes.end(), m2); int k = (int)(it - P2.modes.begin());
      for (int a = 0; a < 3; a++) { x2[(2 * k) * 3 + a] = x[(2 * mu) * 3 + a]; x2[(2 * k + 1) * 3 + a] = x[(2 * mu + 1) * 3 + a]; } }
    int cov = cover_multiplicity(P2, x2.data()); Problem P3; std::vector<double> x3; unwind_cover(P2, x2.data(), cov, P3, x3);
    for (double& v : x3) v /= std::pow(2.0, 2.0 / 3.0); double dist4 = loop_distance(P, x.data(), P3, x3.data());
    CHECK(cov == 2 && dist4 < 1e-9, "cover multiplicity %d, unwound distance %.1e", cov, dist4); }

  std::printf("[Newton on the circular choreography]\n");
  { int N = 4, d = 3; Problem P; P.init(N, d, 12); Work w; w.resize(P); double R = circle_radius(N);
    std::vector<double> x(P.n, 0.0); x[0] = R * 1.05; x[d + 1] = R * 0.97; x[2] = 0.05;
    Reduced Rd(P, nullptr, 0);
    auto fn = [&](const double* y, double* g) { return Rd.f(y, g, w); }; auto fh = [&](const double* y, std::vector<double>& H) { return Rd.hess(y, H, w); };
    OptResult res = newton_lm(P.n, x, fn, fh, 50, 1e-12);
    CurveStats S = curve_stats(P, x.data(), w); std::vector<double> H; hessian(P, x.data(), H, w); Inertia I = inertia(P.n, H);
    std::vector<double> sv; std::vector<double> xc = x; int deff = canonical_frame(P.nb, P.d, xc, sv);
    CHECK(res.converged && std::fabs(S.rms - R) < 1e-10 && deff == 2, "N=%d d=%d converged in %d its, |grad| %.1e, rms %.10f (R=%.10f), deff %d, energy %.8f, morse %d nullity %d", N, d, res.iters, res.gnorm, S.rms, R, deff, S.energy, I.neg, I.zero); }

  std::printf("[random search smoke test N=3 d=2]\n");
  { Problem P; P.init(3, 2, 16); Work w; w.resize(P); Reduced Rd(P, nullptr, 0); int found = 0; std::vector<double> actions;
    auto t0 = std::chrono::steady_clock::now();
    for (int trial = 0; trial < 12; trial++) {
      la::Rng rng(1000 + trial); std::vector<double> x; random_guess(P, rng, 3, 1.0, x); double kin, pot; double A = action_grad(P, x.data(), nullptr, w, &kin, &pot); if (!std::isfinite(A)) continue;
      double lam = optimal_scale(kin, pot, 1.0); for (double& v : x) v *= lam;
      auto fn = [&](const double* y, double* g) { return Rd.f(y, g, w); }; auto fh = [&](const double* y, std::vector<double>& H) { return Rd.hess(y, H, w); };
      lbfgs(P.n, x, fn, 150, 1e-9); OptResult r2 = newton_lm(P.n, x, fn, fh, 40, 1e-11);
      if (r2.converged) { found++; actions.push_back(r2.f); }
    }
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::string s; for (double a : actions) { char b[32]; std::snprintf(b, sizeof b, "%.6f ", a); s += b; }
    CHECK(found >= 3, "%d/12 trials converged (%.0f ms total); actions: %s", found, ms, s.c_str()); }

  std::printf("[binary catalog round trip]\n");
  { Catalog cat; Problem P; P.init(3, 2, 8); la::Rng rng(5);
    for (int k = 0; k < 3; k++) { Record r; std::vector<double> x; random_guess(P, rng, 3, 1.0, x); r.set_solution(P, x.data()); r.h.action = 1.5 + k; r.h.energy = -1; r.h.rms = 1; r.sym = k ? "t+1/2 s[-1,-2]" : ""; r.Lsv = {1, 0.5}; r.pca = {1, 1}; cat.push(r); }
    std::string path = "/tmp/hypchor_test.cat"; cat.save(path); Catalog c2; bool ok = c2.load(path); std::remove(path.c_str());
    double err = 0; bool same = ok && c2.recs.size() == 3;
    if (same) for (int k = 0; k < 3; k++) { same = same && c2.recs[k].sym == cat.recs[k].sym && c2.recs[k].modes == cat.recs[k].modes && c2.recs[k].h.id == k; for (size_t i = 0; i < cat.recs[k].coef.size(); i++) err = std::max(err, std::fabs(cat.recs[k].coef[i] - c2.recs[k].coef[i])); }
    long dup = same ? c2.find_duplicate(cat.recs[1]) : -2;
    CHECK(same && err == 0 && dup == 1, "3 records written/read bit-exact (%zu coefficients each), duplicate lookup -> %ld", cat.recs[0].coef.size(), dup); }

  std::printf("\n%s (%d failures)\n", fails ? "SOME TESTS FAILED" : "ALL TESTS PASSED", fails);
  return fails ? 1 : 0;
}
