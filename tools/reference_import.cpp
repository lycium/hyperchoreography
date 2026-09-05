// Compare/import the equal-mass Table S2 reference set of Li & Liao (2025).
// No network, no overwrite, no automatic merge and no interval-proof or stability claim.
#include "../src/search.hpp"
#include <fstream>
#include <set>

int main(int argc, char** argv) {
  if (argc < 3 || argc > 4) {
    std::fprintf(stderr, "usage: %s data/li-liao-2025-choreographies.tsv existing.bin [new-reference-seeds.bin]\n", argv[0]);
    return 2;
  }
  try {
    std::ifstream input(argv[1]); if (!input) throw std::runtime_error("cannot read reference table");
    Catalog existing, seeds; if (!existing.load(argv[2])) throw std::runtime_error("cannot read comparison catalogue");
    if (argc == 4 && std::filesystem::exists(argv[3])) throw std::runtime_error("refusing to overwrite an existing seed catalogue");
#ifdef HAVE_MPFR
    mpreal::set_default_prec(256);
#endif
    std::puts("# Li & Liao, arXiv:2508.08568v1, equal masses. S/U labels are reported, not re-proved.");
    std::puts("source_id\treported_stability\tvirial_action\traw_shift\trefined_shift\tstate_change\tcoefficient_shift\tmatched_id\tstatus");
    std::string line; std::set<int> seen; int accepted = 0;
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream row(line); int id; double z0, vx, vy, vz, period; char stable; std::string trailing;
      if (!(row >> id >> z0 >> vx >> vy >> vz >> period >> stable) || (row >> trailing) || id < 1 || !seen.insert(id).second ||
          (stable != 'S' && stable != 'U') || !(period > 0) || !std::isfinite(period) ||
          !std::isfinite(z0) || !std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz)) throw std::runtime_error("invalid reference row");
      // Kepler scaling: q_new=lambda(q_old-CM), v_new=lambda^(-1/2)v_old, T_new=2pi.
      const double lambda = std::pow(2 * PI / period, 2.0 / 3.0), vs = 1 / std::sqrt(lambda);
      std::vector<double> p = {-lambda, 0, -lambda*z0/3, lambda, 0, -lambda*z0/3, 0, 0, 2*lambda*z0/3};
      std::vector<double> v = {vs*vx, vs*vy, vs*vz, vs*vx, vs*vy, -vs*vz, -2*vs*vx, -2*vs*vy, 0};
      const double action = -2 * PI * nbody_energy(3, 3, 1, p.data(), v.data()); // periodic-orbit virial diagnostic
      NBody<double> nb(3, 3, 1, 22); double scale = 1;
      for (double x : p) scale = std::max(scale, std::fabs(x)); for (double x : v) scale = std::max(scale, std::fabs(x));
      const double fwd = chore_residual(nb, p, v, 1e-16);
      for (int a = 0; a < 3; a++) { std::swap(p[3+a], p[6+a]); std::swap(v[3+a], v[6+a]); }
      const double rev = chore_residual(nb, p, v, 1e-16);
      if (fwd < rev) for (int a = 0; a < 3; a++) { std::swap(p[3+a], p[6+a]); std::swap(v[3+a], v[6+a]); }
      const double raw = std::min(fwd, rev) / scale;
      Orbit orbit; orbit.Z = p; orbit.Z.insert(orbit.Z.end(), v.begin(), v.end()); const auto before = orbit.Z;
      Ctx ctx; double residual = INFINITY, change = 0, coef_error = INFINITY; long match = -1; std::string status = "raw shift gate failed";
      if (raw < 1e-4) {
        // Keep the benchmark bounded: no arbitrary-precision rescue of a distant, unstable seed.
        residual = certify_state(3, 3, 1, orbit.Z, nullptr, nullptr, nb, ctx.sw, 1e-12, 10, 0, 1);
        for (size_t i = 0; i < before.size(); i++) change = std::max(change, std::fabs(orbit.Z[i] - before[i]) / scale);
        status = "refinement/identity gate failed";
        if (residual <= 1e-10 && change < 1e-5) {
          status = "Fourier fit failed";
          if (fit_loop(3, 3, 1, orbit, ctx, 6144, 512, 2048, nullptr, nb, &coef_error)) {
            Record r; std::vector<double> sv, rotation;
            const int deff = canonical_frame(2 * (int)orbit.modes.size(), 3, orbit.coef, sv, 1e-8, &rotation);
            r.h.id = id; r.h.N = 3; r.h.d = 3; r.h.deff = deff; r.h.K = orbit.K; r.h.M = orbit.Ms;
            r.h.action = orbit.action; r.h.energy = orbit.energy; r.h.rms = orbit.rms; r.h.maxr = orbit.maxr;
            r.h.minsep = orbit.minsep; r.h.Lnorm = orbit.Lnorm; r.Lsv = orbit.Lsv; r.pca = sv;
            r.modes.assign(orbit.modes.begin(), orbit.modes.end()); r.h.nm = (int)r.modes.size(); r.coef = orbit.coef;
            r.h.morse = r.h.nullity = -1; r.h.grad_norm = -1; r.h.hits = 0;
            r.extra.assign(Record::NEX + 9 + 18, 0); r.extra[5] = 1;
            for (int h = 0; h < 2; h++) for (int j = 0; j < 3; j++) for (int a = 0; a < 3; a++)
              for (int b = 0; b < 3; b++) r.extra[Record::NEX + 9 + h*9 + j*3 + a] += rotation[a*3+b] * orbit.Z[h*9+j*3+b];
            record_residuals(r, r.h.ret_err, r.extra[6]); coef_error = r.extra[6];
            r.extra[4] = rigid_defect_coef(3, 3, orbit.modes, orbit.coef.data());
            status = "post-fit gate failed";
            if (r.h.ret_err <= 1e-10 && coef_error <= 1e-6 && deff == 3 && r.rigid() > 1e-4) {
              r.add_citation({2508, 8568, 1, id});
              match = existing.find_duplicate(r); if (match >= 0) match = existing.recs[match].h.id;
              seeds.push(r); ++accepted; status = "numerical reference seed";
            }
          }
        }
      }
      std::printf("%d\t%c\t%.12g\t%.3e\t%.3e\t%.3e\t%.3e\t%ld\t%s\n", id, stable, action, raw, residual, change, coef_error, match, status.c_str()); std::fflush(stdout);
    }
    if (seen.size() != 21) throw std::runtime_error("expected the 21 Table S2 references; no output written");
    if (argc == 4) seeds.save(argv[3]);
    std::fprintf(stderr, "%d of 21 references passed the numerical seed gates; no existence/stability proof asserted.\n", accepted);
    return 0;
  } catch (const std::exception& e) { std::fprintf(stderr, "reference import: %s\n", e.what()); return 1; }
}
