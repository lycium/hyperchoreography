// Numerical self-checks (make test).
#include "catalog.hpp"
#include "search.hpp"
#include "continue.hpp"
#include "optim.hpp"
#include "taylor.hpp"
#ifdef HAVE_MPFR
#include "prove.hpp"
#endif
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

  std::printf("[N-gon spectrum]\n");
  { double e1 = 0, e2 = 0;
    for (int N = 3; N <= 12; N++) for (double al : {1.0, 1.5}) { double C[128]; ngon_C(N, al, C); e1 = std::max(e1, std::fabs(C[0] - C[1] - 1)); }
    for (int N = 4; N <= 10; N++) { double C[128]; ngon_C(N, 1.0, C);
      for (int k = 0; k < N; k++) e2 = std::max(e2, std::fabs(std::sqrt(std::max(0.0, C[0] - C[k])) - ngon_vertical_freq(N, k))); }
    CHECK(e1 < 1e-13 && e2 < 1e-13, "radius equation C_0 - C_1 = 1 to %.1e, omega_k^2 = C_0 - C_k to %.1e", e1, e2); }
  { double a[4], b[4], c[4]; const int na = ngon_inplane_freq(7, 2, 1.0, a), nb = ngon_inplane_freq(6, 3, 1.0, b), nc = ngon_inplane_freq(8, 0, 1.0, c);
    CHECK(na == 2 && std::fabs(a[0] - 0.501584) < 1e-6 && std::fabs(a[1] - 1.014133) < 1e-6 &&
          nb == 2 && std::fabs(b[1] - 0.468728) < 1e-6 &&
          nc == 2 && std::fabs(c[0] + 1) < 1e-12 && std::fabs(c[1] - 1) < 1e-12,     // k=0 is nu^2 (nu^2 + alpha - 2)
          "in-plane roots N=7 k=2 (%.6f, %.6f), N=6 k=3 (%.6f), Kepler pair at k=0 (%.3f, %.3f)", a[0], a[1], b[1], c[0], c[1]); }
  { Problem P; P.init(8, 3, 24); la::Rng rng(3); std::vector<double> x; int got = 0, least = 99;
    for (int t = 0; t < 32; t++) if (inplane_guess(P, rng, x)) { got++; int nmod = 0;
      for (int mu = 0; mu < P.nm; mu++) for (int c = 0; c < 2; c++) if (std::fabs(x[2 * mu * 3 + c]) + std::fabs(x[(2 * mu + 1) * 3 + c]) > 1e-12) { nmod++; break; }
      least = std::min(least, nmod); }
    CHECK(got >= 28 && least >= 2, "inplane start: %d/32 resonant, at least %d in-plane modes each", got, least); }

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
  { const int N = 5, d = 7, K = 10; Problem P; P.init(N, d, K); la::Rng rng(7);   // ∂(∇A)/∂s along Ω = sΩ₀
    ContParam pp; pp.frame(g2_omega(1, 2), d); Work w; std::vector<double> x, ga(P.n), gp(P.n), gm(P.n);
    const double s = 0.7, h = 1e-6;
    pp.set(P, s); w.resize(P); random_guess(P, rng, 5, 1.0, x);
    pp.dgrad(P, x.data(), ga.data(), s, w);
    pp.set(P, s + h); action_grad(P, x.data(), gp.data(), w);
    pp.set(P, s - h); action_grad(P, x.data(), gm.data(), w);
    double e = 0, sc = 0; for (int i = 0; i < P.n; i++) { e = std::max(e, std::fabs((gp[i] - gm[i]) / (2 * h) - ga[i])); sc = std::max(sc, std::fabs(ga[i])); }
    CHECK(e < 1e-6 * (1 + sc), "d(grad A)/ds at Ω = sΩ₀ vs central differences: %.1e (scale %.1f)", e, sc); }
  { const int N = 5, d = 7, K = 10; Problem P; P.init(N, d, K); la::Rng rng(11);  // gauge = time shift + centraliser
    std::vector<double> x, G; random_guess(P, rng, 5, 1.0, x);
    int g0 = gauge_basis(P, x.data(), G);                                    // inertial: all of so(7)
    P.set_omega(g2_omega(1, 2)); int g1 = gauge_basis(P, x.data(), G);       // distinct rates: the 3-torus
    P.set_omega(g2_omega(0, 1)); int g2r = gauge_basis(P, x.data(), G);      // rates (0,1,1): so(3) ⊕ u(2)
    CHECK(g0 == 22 && g1 == 4 && g2r == 8, "gauge dim: inertial %d, g2 torus %d, degenerate rates %d (expect 22 4 8)", g0, g1, g2r); }

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

  { bool ok = true; std::string line;                                             // the g2 frame extends past d = 7
    for (int d : {7, 9, 11}) {
      const Wedge& W = wedge_basis(d, 3); auto phi = calib_psi(W, CAL_G2);
      std::vector<double> rest((d - 7) / 2); for (size_t i = 0; i < rest.size(); i++) rest[i] = 1.5 + (double)i;
      double g = calib_defect(W, phi, g2_omega(1, 6, d, rest));
      std::vector<double> sw(d / 2 - 1); for (size_t i = 0; i < sw.size(); i++) sw[i] = 1.0 + (double)i;
      double su = calib_defect(W, phi, su_omega(d, sw));
      int stab = calib_stab_dim(W, phi), want = 14 + (d - 7) * (d - 8) / 2;
      ok = ok && g < 1e-12 && su > 0.1 && stab == want;
      char b[96]; std::snprintf(b, sizeof b, " d=%d %.0e/%.2f dim %d", d, g, su, stab); line += b; }
    CHECK(ok, "g2 frame is phi-calibrated in every d >= 7 where su: is not (defect g2/su, stab dim 14+so(d-7)):%s", line.c_str()); }

  { const int N = 5, d = 3, K = 12; Problem P; P.init(N, d, K); Work w; w.resize(P);   // rigidity from coefficients
    std::vector<double> x(P.n, 0.0); const double R = 0.9;                             // must match the Problem form
    x[0 * d + 0] = R; x[1 * d + 1] = R;
    std::vector<int> modes(P.modes.begin(), P.modes.end());
    std::vector<double> coef(modes.size() * 2 * d, 0.0);
    for (size_t mu = 0; mu < modes.size(); mu++) for (int a = 0; a < d; a++)
      { coef[mu * 2 * d + a] = x[(2 * mu) * d + a]; coef[mu * 2 * d + d + a] = x[(2 * mu + 1) * d + a]; }
    double c0 = rigid_defect_coef(N, d, modes, coef.data()), p0 = rigid_defect(P, x.data(), w);
    x[2 * d + 2] = 0.35; coef[1 * 2 * d + 2] = 0.35;                                   // one transverse mode breaks it
    double c1 = rigid_defect_coef(N, d, modes, coef.data()), p1 = rigid_defect(P, x.data(), w);
    CHECK(c0 < 1e-12 && p0 < 1e-12 && c1 > 0.1 && std::fabs(c1 - p1) < 1e-2 * p1,   // sampled on two grids
      "rigid_defect from coefficients matches the Problem form: N-gon %.0e/%.0e, tilted %.4f/%.4f", c0, p0, c1, p1); }

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
  { // layout 1: the Ω block is always present, so an inertial record carries d² zeros and must still read
    // as inertial, and the 2Nd certified state sits after it
    Catalog cat; Problem P; P.init(3, 2, 8); la::Rng rng(11);
    const int N = 3, d = 2, nst = 2 * N * d;
    for (int k = 0; k < 2; k++) { Record r; std::vector<double> x; random_guess(P, rng, 3, 1.0, x); r.set_solution(P, x.data());
      r.h.action = 4.0 + k; r.h.energy = -1; r.h.rms = 1; r.Lsv = {1, 0.5}; r.pca = {1, 1};
      r.extra.assign(Record::NEX + (size_t)d * d + nst, 0.0); r.extra[5] = 1; r.extra[6] = 3.25e-12;
      if (k) { r.extra[Record::NEX + 1] = 0.75; r.extra[Record::NEX + 2] = -0.75; }   // record 1 rotates
      for (int i = 0; i < nst; i++) r.extra[Record::NEX + (size_t)d * d + i] = 0.5 + i;
      cat.push(r); }
    std::string path = "/tmp/hypchor_test2.cat"; cat.save(path); Catalog c2; bool ok = c2.load(path); std::remove(path.c_str());
    bool same = ok && c2.recs.size() == 2 && c2.recs[0].extra == cat.recs[0].extra && c2.recs[1].extra == cat.recs[1].extra;
    same = same && c2.recs[0].layout() == 1 && c2.recs[0].omega() == nullptr && c2.recs[1].omega() != nullptr;
    const double* z = ok ? c2.recs[0].state() : nullptr;
    same = same && z != nullptr && c2.recs[0].coef_err() == 3.25e-12;
    if (z) for (int i = 0; i < nst; i++) same = same && z[i] == 0.5 + i;
    // and the state block must not be mistaken for part of Ω
    same = same && c2.recs[1].omega()[1] == 0.75 && c2.recs[1].state()[0] == 0.5;
    CHECK(same, "layout 1 round trips: zeroed Omega block reads as inertial, %d-entry state and coef_err recovered", nst); }
  { Catalog cat; Problem P; P.init(3, 2, 8); la::Rng rng(9);        // two points of one continuous family
    Record a, b; std::vector<double> x, y;                          // different loops, identical action
    random_guess(P, rng, 3, 1.0, x); random_guess(P, rng, 3, 1.0, y);
    a.set_solution(P, x.data()); b.set_solution(P, y.data());
    for (Record* r : {&a, &b}) { r->h.action = 19.827310867; r->h.energy = -10.5187151; r->h.rms = 1.0; r->h.minsep = 1.8944; }
    cat.push(a);
    double raw = loop_distance(a.h.N, a.mode_list(), a.h.d, a.coef.data(), b.mode_list(), b.h.d, b.coef.data());
    long dup = cat.find_duplicate(b);                               // same action, genuinely different loops
    Record c = b; c.h.action *= 1 + 1e-6;    long fa = cat.find_duplicate(c);   // a real action gap: distinct orbit
    Record e = b; e.h.minsep *= 1 + 1e-2;    long fm = cat.find_duplicate(e);   // same action, different geometry
    CHECK(dup == 0 && fa < 0 && fm < 0 && raw > 1e-2, "continuous family folds on action+energy though the loops are %.2f apart; a 1e-6 action gap or a 1%% minsep gap does not fold", raw); }


#ifdef HAVE_MPFR
  std::printf("[interval arithmetic]\n");
  { ival::prec() = 128; la::Rng rng(11); double worst = 0; int bad = 0;
    for (int t = 0; t < 2000; t++) {
      double a0 = rng.normal(), b0 = rng.normal(), wa = std::fabs(rng.normal()) * 0.1, wb = std::fabs(rng.normal()) * 0.1;
      ival A(a0 - wa, a0 + wa), B(b0 - wb, b0 + wb), R;
      for (int s = 0; s < 8; s++) { double x = a0 + wa * (2.0 * rng.uniform() - 1), y = b0 + wb * (2.0 * rng.uniform() - 1);
        auto in = [&](const ival& I, double v) { return mpfr_cmp_d(I.lo, v) <= 0 && mpfr_cmp_d(I.hi, v) >= 0; };
        add(R, A, B); if (!in(R, x + y)) bad++;   sub(R, A, B); if (!in(R, x - y)) bad++;   mul(R, A, B); if (!in(R, x * y)) bad++;
        if (!(mpfr_sgn(B.lo) <= 0 && mpfr_sgn(B.hi) >= 0)) { div(R, A, B); if (!in(R, x / y)) bad++; }
        R = A; mul(R, R, B); if (!in(R, x * y)) bad++;                                 // aliasing
        R = A; fma_sub(R, A, B); if (!in(R, x - x * y)) bad++;
        ival S2; mul(S2, A, A); add_inplace(S2, ival(1)); sqrt_(R, S2); if (!in(R, std::sqrt(x * x + 1))) bad++;
        pow_d(R, S2, -1.5); if (!in(R, std::pow(x * x + 1, -1.5))) bad++;
        worst = std::max(worst, R.wid()); }
      R = A; sub_inplace(R, A); if (!(mpfr_sgn(R.lo) <= 0 && mpfr_sgn(R.hi) >= 0)) bad++; }
    CHECK(bad == 0, "%d containment failures in 2000 random operand pairs (aliasing included)", bad); }

  std::printf("[tangent series]\n");
  { int N = 4, d = 3, order = 10, nd = N * d, n2 = 2 * nd; NBody<double> nb(N, d, 1.0, order), nbp(N, d, 1.0, order), nbm(N, d, 1.0, order); Tangent<double> tg(N, d, 1.0, order);
    la::Rng rng(5); std::vector<double> Z(n2), dir(n2), Zp, Zm; for (double& z : Z) z = rng.normal(); for (double& v : dir) v = rng.normal();
    nb.series(Z.data(), Z.data() + nd); tg.series(nb, dir.data(), dir.data() + nd);
    const double eps = 1e-6; Zp = Z; Zm = Z; for (int i = 0; i < n2; i++) { Zp[i] += eps * dir[i]; Zm[i] -= eps * dir[i]; }
    nbp.series(Zp.data(), Zp.data() + nd); nbm.series(Zm.data(), Zm.data() + nd); double err = 0;
    for (int k = 0; k <= order; k++) for (int i = 0; i < nd; i++) {
      double fx = (nbp.X[(size_t)k * nd + i] - nbm.X[(size_t)k * nd + i]) / (2 * eps), fv = (nbp.V[(size_t)k * nd + i] - nbm.V[(size_t)k * nd + i]) / (2 * eps);
      err = std::max(err, std::fabs(fx - tg.X[(size_t)k * nd + i]) / (1 + std::fabs(fx))); err = std::max(err, std::fabs(fv - tg.V[(size_t)k * nd + i]) / (1 + std::fabs(fv))); }
    CHECK(err < 1e-6, "linearised Taylor coefficients match central differences to %.1e over orders 0..%d", err, order); }

  std::printf("[validated flow]\n");
  { // the figure eight (Simo's initial conditions), rescaled to period 2pi and polished by the shooting Newton
    const int N = 3, d = 2, nd = 6, n2 = 12; const double T0 = 6.32591398, lam = T0 / (2 * PI), sq = std::pow(lam, -2.0 / 3), sv = std::pow(lam, 1.0 / 3);
    std::vector<double> Z = { 0.97000436, -0.24308753, -0.97000436, 0.24308753, 0, 0, 0.93240737 / 2 * sv, 0.86473146 / 2 * sv, 0.93240737 / 2 * sv, 0.86473146 / 2 * sv, -0.93240737 * sv, -0.86473146 * sv };
    for (int i = 0; i < nd; i++) Z[i] *= sq;
    NBody<double> nbd(N, d, 1.0, 22); ShootWork<double> W;
    { std::vector<double> p0(Z.begin(), Z.begin() + nd), v0(Z.begin() + nd, Z.end());        // the bodies chase in one of two orders
      if (chore_residual(nbd, p0, v0, 1e-16) > 0.1) for (int h = 0; h < 2; h++) for (int c = 0; c < d; c++) std::swap(Z[h * nd + d + c], Z[h * nd + 2 * d + c]); } double res = shoot_newton(nbd, Z, 1e-16, 30, 1e-14, -18, -40, W, false, (const std::vector<double>*)nullptr);
    ival::prec() = 200; mpreal::set_default_prec(200); const double tol = 1e-30; const int order = 40;
    Verified V(N, d, 1.0, order, 2, 2, tol); for (int i = 0; i < n2; i++) V.Z[i] = ival(Z[i]);
    set_d(V.Psi[0], 1.0); set_d(V.Psi[(size_t)n2 + 7], 1.0);                         // e_0 and e_7
    ival tend = ival::pi() * 2.0 / ival(N); bool ok = V.integrate(tend);
    NBody<mpreal> nbm(N, d, 1.0, 50); std::vector<mpreal> p(nd), v(nd); for (int i = 0; i < nd; i++) { p[i] = Z[i]; v[i] = Z[nd + i]; }
    nbm.integrate(p, v, mpreal::pi() * 2 / N, 1e-40); int inside = 0; double wid = 0;
    for (int i = 0; i < n2; i++) { const mpreal& x = i < nd ? p[i] : v[i - nd]; if (mpfr_cmp(V.Z[i].lo, x.v) <= 0 && mpfr_cmp(x.v, V.Z[i].hi) <= 0) inside++; wid = std::max(wid, V.Z[i].wid()); }
    // tangent columns against central differences of the double flow
    double jerr = 0; for (int c = 0; c < 2; c++) { const int col = c == 0 ? 0 : 7; const double eps = 1e-6;
      std::vector<double> pp(Z.begin(), Z.begin() + nd), vp(Z.begin() + nd, Z.end()), pm = pp, vm = vp;
      (col < nd ? pp[col] : vp[col - nd]) += eps; (col < nd ? pm[col] : vm[col - nd]) -= eps;
      nbd.integrate(pp, vp, 2 * PI / N, 1e-16); nbd.integrate(pm, vm, 2 * PI / N, 1e-16);
      for (int i = 0; i < n2; i++) { double fd = ((i < nd ? pp[i] : vp[i - nd]) - (i < nd ? pm[i] : vm[i - nd])) / (2 * eps); jerr = std::max(jerr, std::fabs(fd - V.Psi[(size_t)c * n2 + i].mid())); } }
    CHECK(res < 1e-12 && ok && inside == n2 && wid < 1e-20 && jerr < 1e-6, "eight over T/3: %ld steps, MPFR flow inside the enclosure (%d/%d, width %.1e), tangents vs finite differences %.1e (Newton %.1e)", V.steps, inside, n2, wid, jerr, res);
    // a frame with a quarter turn (rate N/4) and a half turn: cos(π/2) is a zero crossing, sin(π/2) an
    // extremum, and neither may widen the block beyond rounding
    { const int Nf = 12, df = 5; std::vector<double> Om((size_t)df * df, 0.0);
      Om[0 * df + 1] = -3.0; Om[1 * df + 0] = 3.0; Om[2 * df + 3] = -6.0; Om[3 * df + 2] = 6.0;
      Frame fq = frame_of(Nf, df, Om.data()); std::vector<ival> G; fq.G_ival(G);
      double w = 0; for (auto& g : G) w = std::max(w, g.wid());
      bool blk = std::fabs(G[0 * df + 1].mid() + 1) < 1e-12 && std::fabs(G[1 * df + 0].mid() - 1) < 1e-12 && std::fabs(G[0].mid()) < 1e-12
              && std::fabs(G[2 * df + 2].mid() + 1) < 1e-12 && fq.cls[0] == 2 && fq.cls[1] == 1 && fq.trans.size() == 1 && fq.rots.size() == 2;
      CHECK(blk && w < 1e-30, "frame at a quarter turn and a half turn: G blocks exact to %.1e, classes %d/%d, %zu translation, %zu commuting rotations", w, fq.cls[0], fq.cls[1], fq.trans.size(), fq.rots.size()); }
    std::vector<mpreal> Zm; Frame fr = frame_of(N, d, nullptr); double r2 = refine_state(N, d, 1.0, Z, fr, 20, 2, Zm);
    Proof P = prove_state(N, d, 1.0, Zm, fr, 1e-10, 1e-24, 32, 2);
    CHECK(r2 < 1e-20 && P.ok && P.kappa < 0.1 && mpfr_cmp_d(P.action.lo, 8.12397549) <= 0 && mpfr_cmp_d(P.action.hi, 8.12397550) >= 0,
          "existence proof of the eight: %s, contraction %.1e, closure %.1e, action %s", P.ok ? "proven" : P.why.c_str(), P.kappa, P.closure, P.action.str(12).c_str()); }
#endif

  std::printf("\n%s (%d failures)\n", fails ? "SOME TESTS FAILED" : "ALL TESTS PASSED", fails);
  return fails ? 1 : 0;
}
