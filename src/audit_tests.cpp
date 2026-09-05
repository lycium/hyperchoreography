// Adversarial numerical regressions for the mathematics/performance audit.
#include "search.hpp"
#include "continue.hpp"
#ifdef HAVE_MPFR
#include "prove.hpp"
#endif
#include <cstdio>

static int failures = 0;
#define CHECK(c, msg) do { bool audit_check_ok = (c); std::printf("%s: %s\n", audit_check_ok ? "ok" : "FAIL", msg); failures += !audit_check_ok; } while (0)

int main() {
  {
    Pool pool(3); bool ok = true;
    for (int throwing : {0, 1, 2}) {
      std::atomic<int> finished{0}; bool caught = false;
      try { pool.run([&](int k) { ++finished; if (k == throwing) throw std::runtime_error("audit pool error"); }); }
      catch (const std::runtime_error&) { caught = true; }
      ok &= caught && finished == 3;
      finished = 0; pool.run([&](int) { ++finished; }); ok &= finished == 3;
    }
    CHECK(ok, "worker and caller exceptions finish the job, propagate, and leave the pool reusable");
  }
  {
    Problem p; p.init(3, 2, 4); p.set_omega({0, -1, 1, 0}); p.init(3, 3, 4);
    CHECK(p.Om.empty() && p.Gsh.empty(), "reinitializing a problem clears the previous frame");
    std::vector<double> a = {0, 0, 0, NAN}, w; bool rejected = false;
    try { la::sym_eig(2, a, w); } catch (const std::invalid_argument&) { rejected = true; }
    CHECK(rejected, "eigensolver rejects nonfinite matrices");
  }
  {
    std::vector<double> a = {0, 2, 1, 1, 0, 3, 4, 1, 2}, lu = a; std::vector<int> piv;
    bool ok = la::lu_factor(3, lu, piv);
    for (auto x : {std::vector<double>{2, -1, 3}, std::vector<double>{-3, 2, 1}}) {
      std::vector<double> b(3, 0); for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) b[i] += a[i * 3 + j] * x[j];
      la::lu_substitute(3, lu, piv, b); for (int i = 0; i < 3; i++) ok &= std::fabs(b[i] - x[i]) < 1e-13;
    }
    CHECK(ok, "one pivoted factorization solves multiple independent right-hand sides");
  }
  {
    NBody<double> nb(2, 1, 1, 12); std::vector<double> p = {0.5, -0.5}, v = {-1, 1}, ts = {0}, samples;
    CHECK(nb.integrate(p, v, 0.0, 1e-12, &ts, &samples) == 0 && samples == std::vector<double>({0.5, -0.5, -1, 1}),
          "zero-length dense output returns the initial state");
    CHECK(nb.integrate(p, v, 0.1, 0.0) < 0, "invalid tolerance fails without a stalled time loop");
    p = {0, 0}; v = {0, 0};
    CHECK(!std::isfinite(chore_residual(nb, p, v, 1e-12)), "collision cannot be reported as a successful return");
  }
  {
    Config cfg; cfg.N = 3; cfg.d = 2; cfg.K = 4; cfg.gtol = 1e-10; Ctx ctx;
    const Problem& p = ctx.problem(cfg, cfg.K); ctx.w.resize(p); ContCfg cc;
    Continuer cont(cfg, cc, ctx, p); Problem pa = p;
    std::vector<double> x(p.n, 0), tangent(p.n + 1, 0);
    x[0] = x[p.d + 1] = std::cbrt(1 / std::sqrt(3.0)); const auto xp = x;
    tangent[p.n] = 1; double alpha = 1;
    int it = cont.corrector(pa, x, alpha, tangent, xp, 1.2);
    CHECK(it > 0 && std::fabs(alpha - 1.2) < 1e-9, "continuation enforces arclength even when the initial action gradient is zero");
  }
  {
    const int d = 64, h = d / 2; const double theta = 8;
    std::vector<double> a((size_t)d * d, 0), r;
    for (int i = 0; i < h; i++) for (int j = h; j < d; j++) { a[(size_t)i * d + j] = -theta / h; a[(size_t)j * d + i] = theta / h; }
    la::expm_skew(d, a, r); double err = 0;
    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) {
      double expected = (i == j) + ((i < h) == (j < h) ? (std::cos(theta) - 1) / h : (i < h ? -1 : 1) * std::sin(theta) / h);
      err = std::max(err, std::fabs(r[(size_t)i * d + j] - expected));
    }
    CHECK(err < 1e-12, "matrix exponential scales by a matrix norm in high dimension");
  }
#ifdef HAVE_MPFR
  ival::prec() = 128;
  mpreal::set_default_prec(256);
  {
    bool ok = true;
    for (double a : {-1.1, 0.0, 3.0}) {
      const ival b(-2.4, 0.7); ival add_ref(-0.3, 0.2), sub_ref(add_ref), added(add_ref), subbed(add_ref);
      fma_add(add_ref, ival(a), b); fma_sub(sub_ref, ival(a), b);
      fma_add_d(added, a, b); fma_sub_d(subbed, a, b);
      ok &= added.contains(add_ref) && subbed.contains(sub_ref);
      ival alias(b), ref(b); fma_sub(ref, ival(a), b); fma_sub_d(alias, a, alias); ok &= alias.contains(ref);
    }
    CHECK(ok, "allocation-free mixed interval products retain outward bounds and alias safety");
  }
  {
    // Exact zero-energy radial Kepler fall: q=+/-0.5(1-3t)^(2/3), v=mp(1-3t)^(-1/3).
    // Its action is 2[1-(1-3t)^(1/3)], independent of the Taylor implementation.
    bool enclosed = true; const double h = std::ldexp(1.0, -10);
    mpreal u = pow(mpreal(1) - mpreal(h) * 3, mpreal(1) / 3), exact = (mpreal(1) - u) * 2;
    for (int order : {2, 3, 5, 8}) {
      Verified flow(2, 1, 1, order, 0, 1, 1e-8);
      flow.Z = {ival(0.5), ival(-0.5), ival(-1), ival(1)};
      bool ok = flow.integrate(ival(h));
      enclosed &= ok && mpfr_cmp(flow.A.lo, exact.v) <= 0 && mpfr_cmp(flow.A.hi, exact.v) >= 0;
    }
    CHECK(enclosed, "validated action encloses an analytic nonconstant Kepler action at several orders");
  }
  {
    std::vector<double> om(16, 0); om[1] = -0.1; om[4] = 0.1; om[11] = -5.1; om[14] = 5.1;
    Frame fr = frame_of(5, 4, om.data()); std::vector<ival> g; fr.G_ival(g);
    bool commute = true;
    for (const auto& xi : fr.rots) for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
      ival entry;
      for (int k = 0; k < 4; k++) { fma_add(entry, ival(xi[i * 4 + k]), g[k * 4 + j]); fma_sub(entry, g[i * 4 + k], ival(xi[k * 4 + j])); }
      commute &= entry.contains(ival(0));
    }
    CHECK(commute, "snapped frame generators commute with the exact certified rotation");
  }
  {
    mpreal exact = -mpreal(1) - mpreal::two_pow(-60);
    ival x(exact.v);
    CHECK(mpfr_cmp_d(exact.v, -x.mag()) >= 0,
          "magnitude bounds a negative endpoint between adjacent doubles");
    mpreal tiny = -mpreal::two_pow(-1100);
    ival y(tiny.v);
    CHECK(y.mag() > 0, "magnitude does not erase a negative subnormal bound");
  }
  {
    int bad = 0;
    for (auto bounds : {std::pair<double, double>{0.25, 0.5}, {2.0, 2.5}, {-2.5, -2.0}, {-8.0, 8.0}}) {
      ival th(bounds.first, bounds.second), c, s;
      icos_isin(th, c, s);
      for (int j = 0; j <= 64; j++) {
        mpreal t = mpreal(bounds.first) + (mpreal(bounds.second) - mpreal(bounds.first)) * j / 64;
        mpreal cv = cos(t), sv = sin(t);
        bad += mpfr_cmp(c.lo, cv.v) > 0 || mpfr_cmp(c.hi, cv.v) < 0;
        bad += mpfr_cmp(s.lo, sv.v) > 0 || mpfr_cmp(s.hi, sv.v) < 0;
      }
    }
    CHECK(bad == 0, "trigonometric intervals enclose endpoints and interior extrema");
  }
  {
    // Serialize at fewer decimal digits than the MPFR precision: printing is part of the certificate.
    ival x(1.23456789, 2.34567891);
    std::string s = x.str(4);
    auto comma = s.find(',');
    mpreal lo(s.substr(1, comma - 1)), hi(s.substr(comma + 2, s.size() - comma - 3));
    CHECK(mpfr_cmp(lo.v, x.lo) <= 0 && mpfr_cmp(hi.v, x.hi) >= 0,
          "printed decimal interval endpoints remain outward rounded");
  }
#endif
  {
    Catalog cat; Problem p; p.init(3, 2, 8); la::Rng rng(9);
    std::vector<double> x, y; random_guess(p, rng, 3, 1.0, x); random_guess(p, rng, 3, 1.0, y);
    Record a, b; a.set_solution(p, x.data()); b.set_solution(p, y.data());
    for (Record* r : {&a, &b}) { r->h.action = 19.827310867; r->h.energy = -10.5187151; r->h.rms = 1; r->h.minsep = 1.8944; }
    cat.push(a);
    CHECK(cat.find_duplicate(b) < 0, "equal action, energy and separation do not prove orbit equivalence");
    for (int mu = 0; mu < p.nm; mu++) {
      double na = 0, nb = 0;
      for (int c = 0; c < 2 * p.d; c++) { na += x[mu * 2 * p.d + c] * x[mu * 2 * p.d + c]; nb += y[mu * 2 * p.d + c] * y[mu * 2 * p.d + c]; }
      if (nb > 0) for (int c = 0; c < 2 * p.d; c++) y[mu * 2 * p.d + c] *= std::sqrt(na / nb);
    }
    b.set_solution(p, y.data()); b.h.minsep *= 1.1;
    CHECK(cat.find_duplicate(b) < 0, "equal per-mode powers do not prove orbit equivalence");
  }
  {
    Problem p; p.init(3, 2, 120); la::Rng rng(3); std::vector<double> x; random_guess(p, rng, 4, 1, x);
    Record a; a.set_solution(p, x.data()); a.h.action = 10; a.h.energy = -2; a.h.rms = 1;
    Catalog cat; cat.push(a);
    bool equivalent = true;
    for (int eps : {-1, 1}) {
      Record b = a;
      for (int mu = 0; mu < p.nm; mu++) {
        double co = std::cos(p.modes[mu] * 0.37), si = std::sin(p.modes[mu] * 0.37);
        for (int h = 0; h < 2; h++) {
          double u[2]; for (int c = 0; c < 2; c++) u[c] = h ? eps * (-si * x[4 * mu + c] + co * x[4 * mu + 2 + c]) : co * x[4 * mu + c] + si * x[4 * mu + 2 + c];
          b.coef[4 * mu + 2 * h] = 0.8 * u[0] - 0.6 * u[1]; b.coef[4 * mu + 2 * h + 1] = 0.6 * u[0] + 0.8 * u[1];
        }
      }
      equivalent &= cat.find_duplicate(b) == 0;
    }
    CHECK(equivalent, "full-coefficient matching retains rotation, phase and reversal equivalence");
    Record b = a; b.coef.back() = 0.1;
    CHECK(cat.find_duplicate(b) < 0, "unmatched high modes survive the low-mode alignment filter");
    a.extra.assign(Record::NEX + 4, 0); b = a; a.extra[Record::NEX + 1] = -0.25; a.extra[Record::NEX + 2] = 0.25;
    b.extra[Record::NEX + 1] = -0.5; b.extra[Record::NEX + 2] = 0.5; Catalog rotating; rotating.push(a);
    CHECK(rotating.find_duplicate(b) < 0, "identical loop coefficients in different frames are not automatically equivalent");
    a.extra[7] = 1e-20; a.extra[5] = 1;
    CHECK(a.proven() == 0 && a.recorded_proof() > 0, "legacy proof flags require recomputation");
  }
  {
    FILE* f = std::tmpfile(); bool rejected = false;
    if (f) {
      RecHdr h; h.N = 3; h.d = 2; std::fwrite(&h, 1, sizeof h - 1, f); std::rewind(f);
      try { Record r; r.read(f); } catch (const std::runtime_error&) { rejected = true; }
      std::fclose(f);
    }
    CHECK(rejected, "a truncated catalogue header is an error, not end of catalogue");
    Record missing, measured; missing.h.ret_err = -1; measured.h.ret_err = 1e-14;
    CHECK(Catalog::better(measured, missing), "missing residuals cannot outrank measured residuals");
  }
  {
    Record r; r.h.N = 3; r.h.d = 3; r.Lsv.resize(3); r.pca.resize(3);
    r.extra.assign(Record::NEX + 9 + 18, 0); r.extra[5] = 1;
    const Record::Citation source{2508, 8568, 1, 62}; r.add_citation(source);
    Record restored; bool ok = restored.from_bytes(r.bytes()) && restored.citations() == std::vector<Record::Citation>{source};
    Catalog cat; r.h.ret_err = 1e-12; cat.push(r); Record better = r;
    better.extra.resize(better.source_offset()); better.h.ret_err = 1e-15;
    cat.absorb(0, better); ok &= cat.recs[0].citations() == std::vector<Record::Citation>{source};
    Catalog second; second.push(better); second.absorb(0, r);
    ok &= second.recs[0].citations() == std::vector<Record::Citation>{source};
    CHECK(ok, "prior-work citations survive serialization and both directions of catalogue replacement");
  }
  std::printf("%d audit regression failures\n", failures);
  return failures ? 1 : 0;
}
