// Numerical self-checks (make test).
#include "catalog.hpp"
#include "search.hpp"
#include "optim.hpp"
#include "taylor.hpp"
#include <cstdio>
#include <chrono>
#include <complex>

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
#ifdef HAVE_ACCELERATE
  { int n = 130; la::Rng rng(5); std::vector<double> A((size_t)n * n), A0;      // above the threshold: dsyevd path
    for (int i = 0; i < n; i++) for (int j = i; j < n; j++) A[i * n + j] = A[j * n + i] = rng.normal();
    A0 = A; std::vector<double> w; la::sym_eig(n, A, w); double err = 0, orth = 0;
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) { double s = 0; for (int j = 0; j < n; j++) s += A0[(size_t)i * n + j] * A[(size_t)j * n + k]; err = std::max(err, std::fabs(s - w[k] * A[(size_t)i * n + k])); }
    for (int k = 0; k < n; k++) { double s = 0; for (int i = 0; i < n; i++) s += A[(size_t)i * n + k] * A[(size_t)i * n + k]; orth = std::max(orth, std::fabs(s - 1)); }
    CHECK(err < 1e-9 && orth < 1e-12, "dsyevd residual %.2e orthonormality %.2e", err, orth); }
#endif

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
    double sym = 0; for (int i = 0; i < P.n; i++) for (int j = 0; j < P.n; j++) sym = std::max(sym, std::fabs(H[(size_t)i * P.n + j] - H[(size_t)j * P.n + i]));
    CHECK(gerr < 1e-6 && hverr < 1e-5 && herr < 1e-9 && sym == 0.0, "N=%d d=%d alpha=%.1f  A=%.6f  grad %.1e  Hv %.1e  H-vs-Hv %.1e  H exactly symmetric (%.0e)", N, d, alpha, A, gerr, hverr, herr, sym);
  }

  std::printf("[symmetry basis]\n");
  { Problem P; P.init(3, 2, 12); Symmetry S = Symmetry::parse("t+1/2 s[-1,-2]; t-0 s[1,-2]", 2); int r; auto B = S.basis(P, r);
    Reduced R(P, &B, r); la::Rng rng(3); std::vector<double> y(r), x(P.n), qa, qb; for (double& yi : y) yi = rng.normal(); R.expand(y.data(), x.data());
    sample_curve(P.modes, 2, x.data(), 64, 1, 0, 2, qa); sample_curve(P.modes, 2, x.data(), 64, 1, PI, 2, qb); double e1 = 0;
    for (size_t k = 0; k < qa.size(); k++) e1 = std::max(e1, std::fabs(qa[k] + qb[k]));
    sample_curve(P.modes, 2, x.data(), 64, -1, 0, 2, qb); double e2 = 0;
    for (int j = 0; j < 64; j++) { e2 = std::max(e2, std::fabs(qa[j * 2] - qb[j * 2])); e2 = std::max(e2, std::fabs(qa[j * 2 + 1] + qb[j * 2 + 1])); }
    double orth = 0; for (int c = 0; c < r; c++) for (int c2 = 0; c2 < r; c2++) { double s = 0; for (int i = 0; i < P.n; i++) s += B[(size_t)i * r + c] * B[(size_t)i * r + c2]; orth = std::max(orth, std::fabs(s - (c == c2))); }
    CHECK(e1 < 1e-12 && e2 < 1e-12 && orth < 1e-12 && r == 8, "antipodal-reversal group r=%d (expect 8), q(t+π)=-q(t) err %.1e, q(-t)=σq(t) err %.1e, orthonormal %.1e", r, e1, e2, orth);
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

  std::printf("[rigidity defect]\n");
  { int N = 5, d = 4; Problem P; P.init(N, d, 12); Work w; w.resize(P); double R = circle_radius(N);
    std::vector<double> x(P.n, 0.0);
    x[0] = R; x[d + 1] = R;                                        // the N-gon: rigid, in any dimension
    double rg = rigid_defect(P, x.data(), w);
    x[2 * d + 2] = 0.2 * R;                                        // one transverse mode breaks rigidity
    double rb = rigid_defect(P, x.data(), w);
    CHECK(rg < 1e-14 && rb > 1e-2, "N-gon %.1e vs N-gon + one transverse mode %.4f", rg, rb); }

  std::printf("[inertial effective dimension]\n");
  { const int N = 4, d = 4; std::vector<int> md = {1};                          // a planar circle, rank 2
    std::vector<double> c(2 * d, 0.0); c[0] = 1.0; c[d + 1] = 1.0;
    std::vector<double> Om((size_t)d * d, 0.0), sv;
    int flat = inertial_deff(N, md, d, c.data(), Om.data(), sv);                // Ω = 0: the loop's own rank
    Om[0 * d + 3] = 1.0; Om[3 * d + 0] = -1.0;                                  // a plane the loop does not span
    int up = inertial_deff(N, md, d, c.data(), Om.data(), sv);
    // a rate that counter-rotates the loop freezes each body, but the N rotated copies still span the plane
    Om.assign((size_t)d * d, 0.0); Om[0 * d + 1] = 1.0; Om[1 * d + 0] = -1.0;
    int keep = inertial_deff(N, md, d, c.data(), Om.data(), sv);
    md[0] = N; Om[0 * d + 1] = N; Om[1 * d + 0] = -N;                           // w ≡ 0 (mod N): the copies coincide
    int down = inertial_deff(N, md, d, c.data(), Om.data(), sv);
    CHECK(flat == 2 && up == 3 && keep == 2 && down == 1, "circle in R^4: Omega=0 -> %d, transverse rate -> %d, counter-rotating -> %d, counter-rotating with w = 0 mod N -> %d (expect 2, 3, 2, 1)", flat, up, keep, down); }

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


  std::printf("[rotating frame]\n");
  { const int N = 3, d = 2, K = 8; double R3 = 0; for (int k = 1; k < N; k++) R3 += 1.0 / (4.0 * std::sin(PI * k / N));
    double R = std::cbrt(R3);                                  // at Ω = −1 the circle is the mode-2 curve, same action
    Problem P; P.init(N, d, K); P.set_omega(parse_omega("-1", d)); Work w; w.resize(P);
    std::vector<double> x(P.n, 0.0), g(P.n); int mu = 1;        // modes are 1,2,4,5,7,8 -> index 1 is m=2
    x[(2 * mu) * d + 0] = R; x[(2 * mu + 1) * d + 1] = R;
    double A = action_grad(P, x.data(), g.data(), w), gn = vnorm(g), ret = return_error(P, x.data());
    CHECK(std::fabs(A - 6.534776057) < 1e-8 && gn < 1e-12 && ret < 1e-12, "twisted circle: A=%.9f (inertial 6.534776057), |grad| %.1e, shift residual %.1e", A, gn, ret); }
  { const int N = 5, d = 7, K = 10; Problem P; P.init(N, d, K); la::Rng rng(31);   // derivatives with a general Ω
    P.set_omega(g2_omega(1, 2)); Work w; w.resize(P);
    std::vector<double> x, g(P.n), gp(P.n), gm(P.n), v(P.n), Hv(P.n), H;
    random_guess(P, rng, 5, 1.0, x); for (double& e : v) e = rng.normal();
    action_grad(P, x.data(), g.data(), w); hessvec(P, x.data(), v.data(), Hv.data(), w); hessian(P, x.data(), H, w);
    double h = 1e-6, gerr = 0, herr = 0;
    for (int i = 0; i < P.n; i++) { auto xp = x; xp[i] += h; double fp = action_grad(P, xp.data(), gp.data(), w);
      xp[i] -= 2 * h; double fm = action_grad(P, xp.data(), gm.data(), w);
      gerr = std::max(gerr, std::fabs((fp - fm) / (2 * h) - g[i]) / (1 + std::fabs(g[i])));
      double hv = 0; for (int k = 0; k < P.n; k++) hv += H[(size_t)i * P.n + k] * v[k];
      herr = std::max(herr, std::fabs(hv - Hv[i]) / (1 + std::fabs(Hv[i]))); }
    CHECK(gerr < 1e-6 && herr < 1e-9, "d=7 g2 frame: grad %.1e, H vs Hv %.1e", gerr, herr); }

  std::printf("[calibration ladder]\n");
  { struct { int d, k, id, dim; } rung[] = { {3,3,CAL_SIMPLE,3}, {4,2,CAL_SL,4}, {6,3,CAL_SL,8}, {7,3,CAL_G2,14}, {8,4,CAL_SPIN7,21}, {10,5,CAL_SL,24} };
    std::string got, want; bool ok = true;
    for (auto& r : rung) { const Wedge& W = wedge_basis(r.d, r.k); int dim = calib_stab_dim(W, calib_psi(W, r.id));
      ok = ok && dim == r.dim; got += " " + std::to_string(dim); want += " " + std::to_string(r.dim); }
    CHECK(ok, "stabiliser dims (vol R3, SU(2), SU(3), G2, Spin(7), SU(5)) =%s of so(d) (expect%s)", got.c_str(), want.c_str()); }

  { la::Rng rng(43); const int d = 7, k = 3; std::vector<int> modes = {1, 2, 3, 5};   // A_k from resonant mode
    std::vector<double> coef(modes.size() * 2 * d);                                   // triples, Vandermonde weights
    for (double& v : coef) v = rng.normal();
    const Wedge& W = wedge_basis(d, k); std::vector<double> A(W.n); jet_form(modes, coef.data(), d, W, A.data());
    const int nm = (int)modes.size(), ns = 2 * nm; std::vector<int> sm(ns);
    std::vector<std::vector<std::complex<double>>> z(ns, std::vector<std::complex<double>>(d));
    for (int i = 0; i < ns; i++) { int mu = i < nm ? nm - 1 - i : i - nm, sn = i < nm ? -1 : 1;   // ascending signed modes
      sm[i] = sn * modes[mu];
      for (int a = 0; a < d; a++) z[i][a] = std::complex<double>(0.5 * coef[mu * 2 * d + a], -0.5 * sn * coef[mu * 2 * d + d + a]); }
    std::vector<std::complex<double>> Ac(W.n, {0, 0});
    for (int p = 0; p < ns; p++) for (int q = p + 1; q < ns; q++) for (int r = q + 1; r < ns; r++) {
      if (sm[p] + sm[q] + sm[r]) continue;
      double V = (double)(sm[q] - sm[p]) * (sm[r] - sm[p]) * (sm[r] - sm[q]);
      for (int i = 0; i < W.n; i++) { const int* I = &W.idx[(size_t)i * k];
        const auto &u = z[p], &v = z[q], &w = z[r];
        Ac[i] += V * (u[I[0]] * (v[I[1]] * w[I[2]] - v[I[2]] * w[I[1]])
                    - u[I[1]] * (v[I[0]] * w[I[2]] - v[I[2]] * w[I[0]])
                    + u[I[2]] * (v[I[0]] * w[I[1]] - v[I[1]] * w[I[0]])); } }
    double err = 0, im = 0, sc = 0;
    for (int i = 0; i < W.n; i++) { std::complex<double> v = std::complex<double>(0, -1) * Ac[i];   // i^{k(k-1)/2} = -i
      err = std::max(err, std::fabs(v.real() - A[i])); im = std::max(im, std::fabs(v.imag())); sc = std::max(sc, std::fabs(A[i])); }
    CHECK(err < 1e-12 * sc && im < 1e-12 * sc, "A_3 = -i sum_{n1+n2+n3=0} V(n) z^z^z reproduces the quadrature to %.1e (scale %.1f, |Im| %.1e)", err, sc, im); }

  { la::Rng rng(11); bool ok = true; std::string line;
    for (int d : {6, 7, 8}) { int k, id; calib_pick(d, k, id);
      const Wedge& W = wedge_basis(d, k); auto psi = calib_psi(W, id);
      std::vector<int> modes = {1, 2, 3, 5, 7}; std::vector<double> coef(modes.size() * 2 * d);
      for (double& v : coef) v = rng.normal();
      std::vector<double> A(W.n); jet_form(modes, coef.data(), d, W, A.data());
      double x0 = calib_max(W, A.data(), psi);
      std::vector<double> S((size_t)d * d, 0.0), R, c2(coef.size(), 0.0);        // random O(d): rotation x reflection
      for (int i = 0; i < d; i++) for (int j = i + 1; j < d; j++) { double g = rng.normal(); S[(size_t)i*d+j] = g; S[(size_t)j*d+i] = -g; }
      la::expm_skew(d, S, R); for (int a = 0; a < d; a++) R[(size_t)a * d] = -R[(size_t)a * d];
      for (size_t mu = 0; mu < modes.size(); mu++) for (int h = 0; h < 2; h++) for (int a = 0; a < d; a++) { double t = 0;
        for (int b = 0; b < d; b++) t += R[(size_t)a*d+b] * coef[mu*2*d + h*d + b]; c2[mu*2*d + h*d + a] = t; }
      std::vector<double> A2(W.n); jet_form(modes, c2.data(), d, W, A2.data());
      double x1 = calib_max(W, A2.data(), psi);
      std::vector<double> flat(coef.size(), 0.0);                                 // loop inside a (k-1)-plane
      for (size_t mu = 0; mu < modes.size(); mu++) for (int h = 0; h < 2; h++) for (int a = 0; a < k - 1; a++) flat[mu*2*d + h*d + a] = rng.normal();
      std::vector<double> Af(W.n); jet_form(modes, flat.data(), d, W, Af.data());
      double lo = calib_max(W, Af.data(), psi);
      ok = ok && std::fabs(x1 - x0) < 1e-9 * x0 && std::fabs(lo) < 1e-12 * x0;
      char b[128]; std::snprintf(b, sizeof b, " %s(d=%d,k=%d) %.1e", calib_group(id), d, k, std::fabs(x1 - x0) / x0); line += b; }
    CHECK(ok, "chi* is O(d)-invariant and vanishes below the rung:%s", line.c_str()); }

  { bool ok = true; std::string line;                                             // rotation rates must sum to zero
    for (int d : {4, 6, 8, 10}) { const int n = d / 2; const Wedge& W = wedge_basis(d, n); auto psi = calib_psi(W, CAL_SL);
      std::vector<double> w(n - 1); for (int i = 0; i < n - 1; i++) w[i] = 0.5 + i;
      auto Om = su_omega(d, w), bad = Om; bad[1] += 0.3; bad[(size_t)d] -= 0.3;
      double good = calib_defect(W, psi, Om), off = calib_defect(W, psi, bad);
      ok = ok && good < 1e-12 && off > 0.1;
      char b[64]; std::snprintf(b, sizeof b, " SU(%d) %.0e/%.2f", n, good, off); line += b; }
    const Wedge& W7 = wedge_basis(7, 3); auto phi = calib_psi(W7, CAL_G2);
    std::vector<double> bad7(49, 0.0); bad7[1] = -1; bad7[7] = 1;
    double g7 = calib_defect(W7, phi, g2_omega(1, 2)), b7 = calib_defect(W7, phi, bad7);
    ok = ok && g7 < 1e-12 && b7 > 0.1;
    CHECK(ok, "sum of rates zero preserves the calibration:%s, g2 torus %.0e/%.2f", line.c_str(), g7, b7); }

  { std::string t6 = named_sym("cyc:1", 6), t7 = named_sym("fano:1", 7);          // named twisted-choreography classes
    Problem P; P.init(5, 6, 16); Symmetry S = Symmetry::parse(t6, 6); int r = 0; auto B = S.basis(P, r);
    CHECK(t6 == "t+1/3 s[3,4,5,6,1,2]" && t7 == "t+1/7 s[2,3,4,5,6,7,1]" && r > 0 && r * 3 <= P.n + 6,
      "named classes: cyc:1 -> \"%s\" (d=6 N=5 K=16: r=%d of n=%d), fano:1 -> \"%s\"", t6.c_str(), r, P.n, t7.c_str()); }

  std::printf("[symmetry detection]\n");
  { const int N = 3, d = 2; std::vector<int> modes = {1, 2, 4, 5, 7, 8};   // cos(odd t) on e_0, sin(even t) on e_1
    std::vector<double> coef(modes.size() * 2 * d, 0.0);                   // gives q(t+π) = diag(-1,1) q, q(-t) = diag(1,-1) q
    la::Rng rng(9);
    for (size_t mu = 0; mu < modes.size(); mu++) { double a = 0.3 + rng.uniform();
      if (modes[mu] % 2) coef[mu * 2 * d + 0] = a; else coef[mu * 2 * d + d + 1] = a; }
    int cont = -1; auto g = detect_symmetry(N, modes, d, d, coef.data(), &cont);
    std::string t; double worst = 0;
    for (auto& e : g) { t += (t.empty() ? "" : " ") + e.text; worst = std::max(worst, e.res); }
    CHECK(g.size() == 4 && cont == 0 && worst < 1e-7 && t == "t+0/1 s[1,2] t+1/2 s[-1,2] t-0/1 s[1,-2] t-1/2 s[-1,-2]",
      "recovers the figure-eight group of a designed loop: |G|=%zu, \"%s\", residual %.1e", g.size(), t.c_str(), worst); }

  { la::Rng rng(17); const int d = 8; std::vector<int> modes = {1, 2, 3, 4, 7};   // 1+2=3 and 1+2+4=7 are both
    std::vector<double> cos_only(modes.size() * 2 * d, 0.0), full(modes.size() * 2 * d);   // resonant
    for (double& v : full) v = rng.normal();
    for (size_t mu = 0; mu < modes.size(); mu++) for (int a = 0; a < d; a++) cos_only[mu * 2 * d + a] = full[mu * 2 * d + a];
    auto rel = [&](const std::vector<double>& cf, int k) {                        // ‖A_k‖ / jet scale
      const Wedge& W = wedge_basis(d, k); std::vector<double> A(W.n); double n = 0;
      jet_form(modes, cf.data(), d, W, A.data());
      for (double v : A) n += v * v;
      return std::sqrt(n) / jet_scale(modes, cf.data(), d, k); };
    double z3 = rel(cos_only, 3), z4 = rel(cos_only, 4), g3 = rel(full, 3);
    CHECK(z3 < 1e-14 && z4 > 1e-3 && g3 > 1e-3,
      "reversal forces R A_k = (-1)^{k(k-1)/2} A_k: with q(-t) = q(t) it kills the twist at k=3 (%.1e, against "
      "%.4f without) but not at k=4 (%.4f)", z3, g3, z4); }

  { const int d = 7, k = 3; const Wedge& W = wedge_basis(d, k); auto psi = calib_psi(W, CAL_G2);
    double chi[2]; int cont[2]; size_t ng[2];
    const int mm[2][3] = { {1, 2, 5}, {1, 2, 3} };                         // 3-plane torus REs: rates non-resonant, then 1+2=3
    for (int c = 0; c < 2; c++) {
      std::vector<int> modes(mm[c], mm[c] + 3); std::vector<double> coef(3 * 2 * d, 0.0);
      for (int p = 0; p < 3; p++) { coef[p * 2 * d + 2 * p] = 1.0; coef[p * 2 * d + d + 2 * p + 1] = 1.0; }
      std::vector<double> A(W.n); jet_form(modes, coef.data(), d, W, A.data());
      chi[c] = calib_max(W, A.data(), psi);
      cont[c] = -1; ng[c] = detect_symmetry(3, modes, d, 6, coef.data(), &cont[c]).size();
    }
    CHECK(chi[0] < 1e-12 && chi[1] > 0.1 && (cont[0] & 1) && (cont[1] & 1) && ng[0] == 0 && ng[1] == 0,
      "both 3-plane relative equilibria have a continuous S1 symmetry, but only the resonant one (1+2=3) has twist: %.1e vs %.4f",
      chi[0], chi[1]); }

  std::printf("[binary catalog round trip]\n");
  { Catalog cat; Problem P; P.init(3, 2, 8); la::Rng rng(5);
    for (int k = 0; k < 3; k++) { Record r; std::vector<double> x; random_guess(P, rng, 3, 1.0, x); r.set_solution(P, x.data()); r.h.action = 1.5 + k; r.h.energy = -1; r.h.rms = 1; r.sym = k ? "t+1/2 s[-1,-2]" : ""; r.Lsv = {1, 0.5}; r.pca = {1, 1};
      if (k) { r.extra.assign(Record::NEX + 4, 0.0); r.extra[0] = 1.25 + k; r.extra[1] = -0.5; r.extra[Record::NEX + 1] = 3.5; }
      cat.push(r); }
    std::string path = "/tmp/hypchor_test.cat"; cat.save(path); Catalog c2; bool ok = c2.load(path); std::remove(path.c_str());
    double err = 0; bool same = ok && c2.recs.size() == 3;
    if (same) for (int k = 0; k < 3; k++) { same = same && c2.recs[k].sym == cat.recs[k].sym && c2.recs[k].modes == cat.recs[k].modes && c2.recs[k].h.id == k && c2.recs[k].extra == cat.recs[k].extra; for (size_t i = 0; i < cat.recs[k].coef.size(); i++) err = std::max(err, std::fabs(cat.recs[k].coef[i] - c2.recs[k].coef[i])); }
    same = same && cat.recs[1].omega() != nullptr && c2.recs[1].omega() != nullptr && cat.recs[0].omega() == nullptr;
    long dup = same ? c2.find_duplicate(cat.recs[1]) : -2;
    CHECK(same && err == 0 && dup == 1, "3 records written/read bit-exact (%zu coefficients, extra[] round trips), duplicate lookup -> %ld", cat.recs[0].coef.size(), dup); }

  std::printf("\n%s (%d failures)\n", fails ? "SOME TESTS FAILED" : "ALL TESTS PASSED", fails);
  return fails ? 1 : 0;
}
