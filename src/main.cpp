#include "search.hpp"
#include "continue.hpp"
#ifdef HAVE_MPFR
#include "prove.hpp"
#endif
#include <thread>
#include <atomic>
#include <mutex>
#include <csignal>
#include <cstring>
#include <cstdio>
#include <iostream>

static std::atomic<bool> g_stop{false};
static void on_signal(int) { g_stop = true; }

struct Args {
  std::vector<std::string> pos; std::map<std::string, std::string> opt;
  Args(int argc, char** argv) { for (int i = 2; i < argc; i++) { std::string a = argv[i];
    if (a.rfind("--", 0) == 0) { std::string k = a.substr(2); size_t eq = k.find('='); if (eq != std::string::npos) { opt[k.substr(0, eq)] = k.substr(eq + 1); continue; }
      if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) opt[k] = argv[++i]; else opt[k] = "1"; } else pos.push_back(a); } }
  bool has(const std::string& k) const { return opt.count(k) > 0; }
  std::string get(const std::string& k, const std::string& def = "") const { auto it = opt.find(k); return it == opt.end() ? def : it->second; }
  double num(const std::string& k, double def) const { auto it = opt.find(k); return it == opt.end() ? def : std::stod(it->second); }
};

static void usage() {
  std::puts("hyperchoreography — N-body choreography search in any dimension\n"
    "  hyperchoreography search  --N 3 --d 2 [--K 16 --Kmax 64] [--sym none|random|cyc:p|fano:p|\"t+1/2 s[-1,-2]; ...\"] [--threads T] [--seed S]\n"
    "                [--trials n] [--minutes m] [--out catalog.bin] [--alpha-start 2 --alpha-steps 8] [--min-deff k --min-rigid 1e-4]\n"
    "                [--lbfgs-min 20 --lbfgs-max 400 --newton 60 --gtol 1e-10 --ret-tol 1e-8 --K0 2 --K0max 6 --minsep 2e-3]\n"
    "                [--phase1 action|gradnorm|mixed] [--seed-from other.bin --kick-min 0.02 --kick-max 0.5]\n"
    "                [--starts random,torus,vertical,hyper,inplane,fano,kick] [--K-index 48] [--Ms 2048 --Kout-max 512 --shoot-tol 1e-12 --ret-double 1e-4]\n"
    "                [--omega \"w1,w2,...\" | --omega su:w1,... | --omega g2:p,q[,r,...]]   rotating frame: q_j(t) = exp(Omega t) q(t + 2 pi j/N)\n"
    "                [--tol-inv 1e-4 --tol-dist 1e-3 --checkpoint 30 --ret-reject 1e-1]   (resumes if catalog/state exist)\n"
    "  hyperchoreography list    catalog.bin [--N n] [--deff k] [--min-deff k] [--sort action|id|hits|twist|rigid]\n"
    "  hyperchoreography show    catalog.bin --id i                (JSON dump of one record)\n"
    "  hyperchoreography export  catalog.bin --id i [--samples 720] [--out curve.csv]   (body positions over one period)\n"
    "  hyperchoreography verify  catalog.bin [--id i] [--gate 1e-9]  (Taylor checks; no --id sweeps the whole file)\n"
    "  hyperchoreography refine  catalog.bin --id i --digits 60 [--K 64] [--threads T] [--out refined.txt]   (MPFR shooting Newton)\n"
    "  hyperchoreography prove   catalog.bin [--id i] [--digits 40] [--order 55] [--radius 1e-20] [--threads T] [--write] [--force] [--verbose]\n"
    "                                                  (interval-arithmetic existence proof: Krawczyk on the shooting map; --write stores the radius, --force redoes proven records)\n"
    "  hyperchoreography merge   out.bin in1.bin in2.bin ... [--min-rigid r --min-deff k]   (union, de-duplicate, re-gate)\n"
    "  hyperchoreography extras  catalog.bin [--K-index 48]        (recompute Morse index, nullity, twist and rigidity)\n"

    "  hyperchoreography symmetry catalog.bin [--id i] [--tol 1e-6] (detect the symmetry group of each stored loop)\n"
    "  hyperchoreography continue catalog.bin [--id i | --root circle] --N n --d d [--K 24] [--covers 7] [--depth 2] [--out found.bin]\n"
    "                [--param alpha --alpha-lo 0.6 --alpha-hi 2.4 | --param omega --s-lo -0.25 --s-hi 1.5]   (Omega = s Omega_0, s = 0 is inertial)\n"
    "                                                  (bifurcation-tree exploration by continuation in the potential exponent)\n"
    "  hyperchoreography bench   [--N 3 --d 2 --K 16]              (kernel timings)");
}

static Catalog load_cat(const std::string& p) { Catalog c; if (!c.load(p)) throw std::runtime_error("cannot open catalog " + p); return c; }
static const Record& rec_by_id(const Catalog& c, long id) { for (auto& r : c.recs) if (r.h.id == id) return r; throw std::runtime_error("no record with id " + std::to_string(id)); }

static int cmd_search(const Args& a) {
  Config cfg;
  cfg.N = (int)a.num("N", 3); cfg.d = (int)a.num("d", 2); cfg.K = (int)a.num("K", 16); cfg.Kmax = (int)a.num("Kmax", std::max(4 * cfg.K, 64));
  cfg.sym = named_sym(a.get("sym", "none"), cfg.d); cfg.threads = (int)a.num("threads", std::thread::hardware_concurrency()); cfg.seed = std::stoull(a.get("seed", "1"));
  cfg.trials = a.has("trials") ? std::stol(a.get("trials")) : LONG_MAX; cfg.minutes = a.num("minutes", 1e30); cfg.out = a.get("out", "catalog.bin");
  if (cfg.trials < 0 || !(cfg.minutes >= 0) || !std::isfinite(cfg.minutes)) throw std::invalid_argument("search: invalid trial or time budget");
  cfg.alpha_start = a.num("alpha-start", 1.0); cfg.alpha_steps = (int)a.num("alpha-steps", 8); cfg.min_deff = (int)a.num("min-deff", 1); cfg.min_rigid = a.num("min-rigid", cfg.min_rigid);
  cfg.lbfgs_min = (int)a.num("lbfgs-min", 20); cfg.lbfgs_max = (int)a.num("lbfgs-max", 400); cfg.newton_iters = (int)a.num("newton", 60);
  cfg.gtol = a.num("gtol", 1e-10); cfg.ret_tol = a.num("ret-tol", 1e-8); cfg.K0min = (int)a.num("K0", 2); cfg.K0max = (int)a.num("K0max", 6);
  cfg.minsep = a.num("minsep", 2e-3); cfg.verbose = a.has("verbose"); cfg.checkpoint_secs = (int)a.num("checkpoint", 30);
  cfg.tol_inv = a.num("tol-inv", 1e-4); cfg.tol_dist = a.num("tol-dist", 1e-3);
  cfg.phase1 = a.get("phase1", "mixed"); cfg.kick_min = a.num("kick-min", 0.02); cfg.kick_max = a.num("kick-max", 0.5);
  cfg.Ms = (int)a.num("Ms", 2048); cfg.Kout_max = (int)a.num("Kout-max", 512); cfg.shoot_tol = a.num("shoot-tol", 1e-12); cfg.ret_reject = a.num("ret-reject", 1e-1); cfg.ret_double = a.num("ret-double", 1e-4);
  cfg.starts = a.get("starts", a.has("seed-from") ? "random,torus,vertical,kick" : "random,torus,vertical"); cfg.K_index = (int)a.num("K-index", 48);
  cfg.omega_text = a.get("omega", ""); cfg.omega = parse_omega(cfg.omega_text, cfg.d);
  if (cfg.threads < 1) cfg.threads = 1;
  std::vector<Record> seeds;
  if (a.has("seed-from")) { Catalog sc = load_cat(a.get("seed-from")); seeds = sc.recs; cfg.seeds = &seeds; std::printf("seeding from %zu catalogue solutions (%s)\n", seeds.size(), a.get("seed-from").c_str()); }
  if (cfg.sym != "none" && cfg.sym != "random") Symmetry::parse(cfg.sym, cfg.d);

  Catalog cat; bool resumed = false; try { resumed = cat.load(cfg.out); } catch (std::exception& e) { std::fprintf(stderr, "%s\n", e.what()); return 1; }
  SearchState st; if (!(st.load(cfg.out + ".state") && st.seed == (int64_t)cfg.seed)) { st = SearchState(); st.seed = (int64_t)cfg.seed; }
  std::printf("hyperchoreography search: d=%d N=%d K=%d..%d sym=\"%s\" omega=\"%s\" starts=%s phase1=%s alpha_start=%g threads=%d seed=%llu out=%s\n", cfg.d, cfg.N, cfg.K, cfg.Kmax, cfg.sym.c_str(), cfg.omega_text.c_str(), cfg.starts.c_str(), cfg.phase1.c_str(), cfg.alpha_start, cfg.threads, (unsigned long long)cfg.seed, cfg.out.c_str());
  if (resumed) std::printf("resuming: %zu records in catalog, trial counter at %lld\n", cat.recs.size(), (long long)st.next_trial);
  std::signal(SIGINT, on_signal); std::signal(SIGTERM, on_signal);

  std::atomic<uint64_t> counter{(uint64_t)st.next_trial}; std::atomic<long> done{0}, found{0}, dups{0}, failed{0}; std::atomic<int> active{cfg.threads};
  std::mutex mu; std::map<std::string, long> reasons; auto tstart = std::chrono::steady_clock::now();
  auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - tstart).count(); };
  std::vector<std::thread> th;
  for (int t = 0; t < cfg.threads; t++) th.emplace_back([&] {
    Ctx ctx;
    while (!g_stop) {
      uint64_t tr = counter.fetch_add(1);
      if (tr >= (uint64_t)cfg.trials || elapsed() > cfg.minutes * 60.0) { g_stop = true; break; }
      TrialOut o; try { o = run_trial(cfg, tr, ctx); } catch (std::exception& e) { o.ok = false; o.why = std::string("exception: ") + e.what(); }
      done++;
      if (o.ok) {
        { std::lock_guard<std::mutex> lk(mu);
          long dup = cat.find_duplicate(o.rec, cfg.tol_inv, cfg.tol_dist);
          if (dup >= 0 && !Catalog::better(o.rec, cat.recs[dup])) { cat.absorb((size_t)dup, o.rec); dups++; continue; } }
        // The Hessian/index and calibration can take minutes. Keep them off the catalog lock;
        // another worker may insert a duplicate meanwhile, so recheck before publishing.
        try { record_extras(cfg, o.rec, ctx); }
        catch (std::exception& e) { o.ok = false; o.why = std::string("metadata exception: ") + e.what(); }
      }
      std::lock_guard<std::mutex> lk(mu);
      if (!o.ok) { failed++; reasons[o.why]++; continue; }
      double dist; long dup = cat.find_duplicate(o.rec, cfg.tol_inv, cfg.tol_dist, &dist);
      if (dup >= 0) { cat.absorb((size_t)dup, o.rec); dups++; continue; }
      o.rec.h.id = -1; size_t idx = cat.push(o.rec); found++; const Record& r = cat.recs[idx];
      std::printf("+ id=%lld d=%d/%d N=%d K=%d A=%.9f E=%.6f morse=%d null=%d minsep=%.3f ret=%.1e cover=%d sym=\"%s\" (trial %llu, %.2fs)\n", (long long)r.h.id, r.h.deff, r.h.d, r.h.N, r.h.K, r.h.action, r.h.energy, r.h.morse, r.h.nullity, r.h.minsep, r.h.ret_err, r.h.cover, r.sym.c_str(), (unsigned long long)tr, r.h.secs);
      std::fflush(stdout);
    }
    active--;
  });
  double last_ck = 0;
  const int64_t base_done = st.trials_done; const double base_elapsed = st.elapsed;
  auto checkpoint = [&] { std::lock_guard<std::mutex> lk(mu); cat.save(cfg.out); st.next_trial = (int64_t)counter.load(); st.trials_done = base_done + done.load(); st.found = (int64_t)cat.recs.size(); st.elapsed = base_elapsed + elapsed(); st.save(cfg.out + ".state"); return cat.recs.size(); };
  while (active.load() > 0) {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    if (elapsed() - last_ck >= cfg.checkpoint_secs) { last_ck = elapsed(); const size_t records = checkpoint();
      std::printf("[%6.0fs] trials %ld (%.1f/s) | unique %zu (+%ld) | duplicates %ld | failed %ld\n", elapsed(), done.load(), done.load() / std::max(1e-9, elapsed()), records, found.load(), dups.load(), failed.load()); std::fflush(stdout); }
  }
  for (auto& t : th) t.join();
  checkpoint();
  std::printf("done: %ld trials in %.1fs (%.1f/s): %ld new unique, %ld duplicates, %ld failed; catalog now %zu records -> %s\n", done.load(), elapsed(), done.load() / std::max(1e-9, elapsed()), found.load(), dups.load(), failed.load(), cat.recs.size(), cfg.out.c_str());
  for (auto& kv : reasons) std::printf("   %-40s %ld\n", kv.first.c_str(), kv.second);
  return 0;
}

static int cmd_list(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); int Nf = (int)a.num("N", 0), deff = (int)a.num("deff", 0), mind = (int)a.num("min-deff", 0); std::string sort = a.get("sort", "action");
  std::vector<const Record*> rs; for (auto& r : cat.recs) if ((!Nf || r.h.N == Nf) && (!deff || r.h.deff == deff) && r.h.deff >= mind) rs.push_back(&r);
  std::sort(rs.begin(), rs.end(), [&](const Record* x, const Record* y) { if (sort == "id") return x->h.id < y->h.id; if (sort == "hits") return x->h.hits > y->h.hits; if (sort == "twist") return x->twist() > y->twist(); if (sort == "rigid") return x->rigid() < y->rigid();
    if (x->h.deff != y->h.deff) return x->h.deff < y->h.deff; if (x->h.N != y->h.N) return x->h.N < y->h.N; return x->h.action < y->h.action; });
  std::printf("%5s %2s/%-2s %2s %4s %14s %12s %5s %4s %7s %8s %10s %6s %8s %5s %3s %s\n", "id", "de", "d", "N", "K", "action", "energy", "morse", "null", "minsep", "ret_err", "twist", "tw_rel", "rigid", "hits", "cov", "sym");
  for (auto r : rs) std::printf("%5lld %2d/%-2d %2d %4d %14.9f %12.7f %5d %4d %7.4f %8.1e %10.4g %6.3f %8.2e %5d %3d %s%s\n", (long long)r->h.id, r->h.deff, r->h.d, r->h.N, r->h.K, r->h.action, r->h.energy, r->h.morse, r->h.nullity, r->h.minsep, r->h.ret_err, r->twist(), r->extra.size() > 1 ? r->extra[1] : 0.0, r->rigid(), r->h.hits, r->h.cover, r->sym.c_str(), r->proven() > 0 ? " proven" : r->recorded_proof() > 0 ? " legacy-proof (rerun)" : "");
  std::printf("%zu records\n", rs.size());
  return 0;
}

static int cmd_show(const Args& a) { Catalog cat = load_cat(a.pos.at(0)); std::printf("%s\n", rec_by_id(cat, (long)a.num("id", 0)).to_json().c_str()); return 0; }

static int cmd_export(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); const Record& r = rec_by_id(cat, (long)a.num("id", 0)); int S = (int)a.num("samples", 720);
  std::string out = a.get("out", "curve_" + std::to_string(r.h.id) + ".csv"); FILE* f = std::fopen(out.c_str(), "w"); if (!f) throw std::runtime_error("cannot write " + out);
  const int N = r.h.N, d = r.h.d; std::fprintf(f, "t");
  for (int k = 0; k < N; k++) for (int c = 0; c < d; c++)
    if (c < 4) std::fprintf(f, ",q%d_%c", k, "xyzw"[c]); else std::fprintf(f, ",q%d_x%d", k, c);
  std::fprintf(f, "\n");
  const double* om = r.omega(); std::vector<double> Rt, Ai;
  for (int j = 0; j < S; j++) { double t = 2 * PI * j / S; std::fprintf(f, "%.10f", t);
    if (om) { Ai.assign(om, om + (size_t)d * d); for (double& e : Ai) e *= t; la::expm_skew(d, Ai, Rt); }
    for (int k = 0; k < N; k++) { double tk = t + 2 * PI * k / N; std::vector<double> q(d, 0.0);
      for (size_t mu = 0; mu < r.modes.size(); mu++) { double c = std::cos(r.modes[mu] * tk), s = std::sin(r.modes[mu] * tk); for (int c2 = 0; c2 < d; c2++) q[c2] += c * r.coef[mu * 2 * d + c2] + s * r.coef[mu * 2 * d + d + c2]; }
      if (om) { std::vector<double> qi(d, 0.0);
        for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) qi[a] += Rt[(size_t)a * d + b] * q[b]; q.swap(qi); }
      for (int c2 = 0; c2 < d; c2++) std::fprintf(f, ",%.12f", q[c2]); }
    std::fprintf(f, "\n"); }
  std::fclose(f); std::printf("wrote %s (%d samples, %d bodies, %dD)\n", out.c_str(), S, N, d); return 0;
}

static int cmd_verify_all(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); double gate = a.num("gate", 1e-9);
  std::printf("%5s %2s/%-2s %2s %4s %10s %10s %10s %10s\n", "id", "de", "d", "N", "K", "ret_err", "state", "coef", "period");
  int bad = 0, nolayout = 0; double worst = 0; long worst_id = -1;
  std::vector<std::pair<double, long>> rank;
  for (const Record& r : cat.recs) { double sres = 0, cres = 0, per = 0; record_residuals(r, sres, cres, &per);
    if (r.layout() < 1) nolayout++;
    rank.emplace_back(sres, (long)r.h.id);
    std::printf("%5lld %2d/%-2d %2d %4d %10.1e %10.1e %10.1e %10.1e%s\n", (long long)r.h.id, r.h.deff, r.h.d, r.h.N, r.h.K, r.h.ret_err, sres, cres, per, sres > gate ? "  **" : "");
    if (!(sres <= gate)) bad++;
    if (sres > worst) { worst = sres; worst_id = (long)r.h.id; } }
  std::printf("%zu records; %d above %.0e (marked **); worst is id %ld at %.2e\n", cat.recs.size(), bad, gate, worst_id, worst);
  if (nolayout) std::printf("%d records predate stored numerical states and carry only coefficients\n", nolayout);
  if (bad) { std::sort(rank.begin(), rank.end(), [](auto& x, auto& y) { return x.first > y.first; });
    std::printf("worst ids:  "); for (size_t i = 0; i < rank.size() && i < 8 && rank[i].first > gate; i++) std::printf("%ld ", rank[i].second); std::printf("\n"); }
  return bad ? 1 : 0;
}

static int cmd_verify(const Args& a) {
  if (!a.has("id")) return cmd_verify_all(a);
  Catalog cat = load_cat(a.pos.at(0)); const Record& r = rec_by_id(cat, (long)a.num("id", 0));
  const int N = r.h.N, d = r.h.d, nd = N * d; const double* om = r.omega();
  double sres = 0, cres = 0, per = 0; record_residuals(r, sres, cres, &per);
  std::vector<double> pos, vel; initial_state(N, d, r.mode_list(), r.coef.data(), om, pos, vel);
  if (const double* z = r.state()) { pos.assign(z, z + nd); vel.assign(z + nd, z + 2 * nd); }
  NBody<double> nb(N, d, r.h.alpha, 22);
  double E0 = nbody_energy(N, d, r.h.alpha, pos.data(), vel.data());
  std::vector<double> p = pos, v = vel; int steps = nb.integrate(p, v, 2 * PI, 1e-16);
  if (steps < 0) throw std::runtime_error("full-period verification integration failed");
  double E1 = nbody_energy(N, d, r.h.alpha, p.data(), v.data());
  std::printf("record %lld: d=%d (deff=%d) N=%d K=%d action=%.12f energy=%.12f\n", (long long)r.h.id, r.h.d, r.h.deff, r.h.N, r.h.K, r.h.action, r.h.energy);
  if (om) std::printf("  rotating frame: Omega stored with the record\n");
  if (!r.state()) std::printf("  ** no numerical state stored (layout 0) — the numbers below are the coefficients' **\n");
  std::printf("  numerical state, shift residual        = %.3e   (stored ret_err %.1e)\n", sres, r.h.ret_err);
  std::printf("  stored coefficients, shift residual    = %.3e   (stored extra[6] %.1e)\n", cres, r.coef_err());
  std::printf("  full-period return |Phi_T(Z) - G^N Z|  = %.3e  (%d Taylor steps)\n  energy drift over one period          = %.3e\n  energy from initial state             = %.12f\n", per, steps, E1 - E0, E0);
  if (std::fabs(E0 - r.h.energy) > 1e-6 * std::max(1.0, std::fabs(r.h.energy))) std::printf("  ** stale record: reconstructed energy %.9g disagrees with stored %.9g **\n", E0, r.h.energy);
  std::printf("  initial conditions (t=0):\n"); for (int k = 0; k < N; k++) { std::printf("   body %d  q =", k); for (int c = 0; c < d; c++) std::printf(" %+.15f", pos[k * d + c]); std::printf("   v ="); for (int c = 0; c < d; c++) std::printf(" %+.15f", vel[k * d + c]); std::printf("\n"); }
  return sres <= a.num("gate", 1e-9) ? 0 : 1;
}

static int cmd_merge(const Args& a) {
  if (a.pos.size() < 2) { usage(); return 1; }
  Catalog out; out.load(a.pos[0]);  // an absent destination is fine; a corrupt one must not be overwritten
  double tol_inv = a.num("tol-inv", 1e-4), tol_dist = a.num("tol-dist", 1e-3);
  double min_rigid = a.num("min-rigid", 0.0); int min_deff = (int)a.num("min-deff", 0);
  long added = 0, merged = 0, dropped = 0;
  auto keep = [&](const Record& r) {
    if (r.h.deff < min_deff) return false;
    if (min_rigid > 0 && (int)r.extra.size() > 4 && r.extra[4] < min_rigid) return false;
    return true; };
  if (min_rigid > 0 || min_deff > 0) {
    std::vector<Record> kept; for (auto& r : out.recs) { if (keep(r)) kept.push_back(r); else dropped++; }
    if (dropped) { Catalog re; for (auto& r : kept) { Record c = r; c.h.id = -1; re.push(c); } out = std::move(re); } }
  for (size_t i = 1; i < a.pos.size(); i++) { Catalog in = load_cat(a.pos[i]);
    for (auto& r : in.recs) { if (!keep(r)) { dropped++; continue; }
      long dup = out.find_duplicate(r, tol_inv, tol_dist); if (dup >= 0) { out.absorb((size_t)dup, r, r.h.hits); merged++; } else { Record c = r; c.h.id = -1; out.push(c); added++; } } }
  out.save(a.pos[0]);
  std::printf("merged: %ld added, %ld duplicates folded, %ld dropped by gate; %s has %zu records\n", added, merged, dropped, a.pos[0].c_str(), out.recs.size()); return 0;
}

static int cmd_extras(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); const int Ki = (int)a.num("K-index", 48); long n = 0;
  for (auto& r : cat.recs) {
    if (r.extra.empty()) r.extra.assign(Record::NEX, 0.0);
    Config cfg; cfg.N = r.h.N; cfg.d = r.h.d; cfg.K = r.h.K; cfg.K_index = Ki; cfg.minsep = 0;
    Ctx ctx;
    record_extras(cfg, r, ctx); n++;
  }
  cat.save(a.pos[0]); std::printf("recomputed extras for %ld records -> %s\n", n, a.pos[0].c_str()); return 0;
}

static int cmd_symmetry(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); double tol = a.num("tol", 1e-6); long id = a.has("id") ? (long)a.num("id", 0) : -1;
  std::printf("%5s %2s/%-2s %2s %14s %10s %5s %5s %4s  %s\n", "id", "de", "d", "N", "action", "twist", "|G|", "shift", "rev", "generators of  q(eps t + theta) = R q(t)");
  for (auto& r : cat.recs) {
    if (id >= 0 && r.h.id != id) continue;
    int cont = 0;
    auto g = detect_symmetry(r.h.N, r.mode_list(), r.h.d, r.h.deff, r.coef.data(), &cont, tol);
    const LoopSym *gen = nullptr, *rev = nullptr;
    int nsh = (cont & 1) ? 0 : 1, nrv = 0;
    for (auto& e : g) {
      if (e.eps > 0) { if (e.p == 0) continue; nsh++; if (!gen || (double)e.p / e.q < (double)gen->p / gen->q) gen = &e; }
      else { nrv++; if (!rev) rev = &e; }
    }
    std::string t = gen ? gen->text : "";
    if (rev) { size_t sp = rev->text.find(' ');
      t += (t.empty() ? "" : "; ") + std::string("t-0/1") + (sp == std::string::npos ? "" : rev->text.substr(sp)); }
    if (cont & 2) t += (t.empty() ? "" : "; ") + std::string("t-* (a circle of reversals)");
    std::string sh = (cont & 1) ? "S1" : std::to_string(nsh), go = cont ? "inf" : std::to_string(nsh + nrv);
    std::printf("%5lld %2d/%-2d %2d %14.9f %10.4g %5s %5s %4s  %s\n", (long long)r.h.id, r.h.deff, r.h.d, r.h.N, r.h.action, r.twist(),
                go.c_str(), sh.c_str(), (nrv || (cont & 2)) ? "yes" : "no", t.empty() ? "(trivial)" : t.c_str());
  }
  return 0;
}

static int cmd_continue(const Args& a) {
  Config cfg; ContCfg cc; cfg.N = (int)a.num("N", 3); cfg.d = (int)a.num("d", 2); cfg.K = (int)a.num("K", 24); cfg.Kmax = (int)a.num("Kmax", 2 * cfg.K);
  const bool omode = a.get("param", "alpha") == "omega";
  cc.lo = omode ? a.num("s-lo", -0.25) : a.num("alpha-lo", 0.6); cc.hi = omode ? a.num("s-hi", 1.5) : a.num("alpha-hi", 2.4); cc.h0 = a.num("h0", 0.02); cc.hmax = a.num("hmax", 0.1); cc.depth = (int)a.num("depth", 2); cc.kick = a.num("kick", 0.05); cc.max_steps = (int)a.num("max-steps", 1500); cc.verbose = !a.has("quiet");
  std::string outp = a.get("out", a.pos.empty() ? "continued.bin" : a.pos[0]);
  Catalog cat; if (!a.pos.empty()) { try { cat.load(a.pos[0]); } catch (...) {} }
  Catalog outc; if (outp != (a.pos.empty() ? "" : a.pos[0])) { try { outc.load(outp); } catch (...) {} } else outc = cat;
  std::vector<std::pair<std::vector<double>, std::string>> roots; std::vector<std::vector<double>> roms; Ctx ctx; long skipped = 0;
  if (!omode && (a.get("root") == "circle" || a.pos.empty())) {
    const Problem& P = ctx.problem(cfg, cfg.K); std::vector<double> x(P.n, 0.0); double R3 = 0; for (int k = 1; k < cfg.N; k++) R3 += 1.0 / (4.0 * std::sin(PI * k / cfg.N)); double R = std::cbrt(R3);
    x[0] = R; x[P.d + 1] = R; roots.emplace_back(x, "circle"); roms.emplace_back();
  } else {
    for (auto& r : cat.recs) { if (a.has("id") && r.h.id != (long)a.num("id", 0)) continue; if (r.h.N != cfg.N || r.h.d > cfg.d) continue;
      const double* om = r.omega();
      if (omode ? !om : om != nullptr) { skipped++; continue; }
      const Problem& P = ctx.problem(cfg, cfg.K); std::vector<double> x(P.n, 0.0);
      for (size_t k = 0; k < r.modes.size(); k++) { int m = r.modes[k]; if (m > P.K) break; auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) continue; int mu = (int)(it - P.modes.begin());
        for (int c = 0; c < r.h.d; c++) { x[(2 * mu) * P.d + c] = r.coef[k * 2 * r.h.d + c]; x[(2 * mu + 1) * P.d + c] = r.coef[k * 2 * r.h.d + r.h.d + c]; } }
      std::vector<double> O; if (omode) { O.assign((size_t)cfg.d * cfg.d, 0.0);
        for (int i = 0; i < r.h.d; i++) for (int j = 0; j < r.h.d; j++) O[(size_t)i * cfg.d + j] = om[(size_t)i * r.h.d + j]; }
      roots.emplace_back(x, "id" + std::to_string(r.h.id)); roms.push_back(O); }
  }
  int covers = (int)a.num("covers", 1); const Problem& P0 = ctx.problem(cfg, cfg.K); size_t nroots = roots.size();
  for (int k = 2; k <= covers; k++) { if (la::gcd(k, cfg.N) != 1) continue;
    for (size_t r0 = 0; r0 < nroots; r0++) { std::vector<double> x(P0.n, 0.0); double lam = std::pow((double)k, -2.0 / 3.0); bool ok = true;
      for (int mu = 0; mu < P0.nm; mu++) { int m = P0.modes[mu]; bool has = false; for (int c = 0; c < P0.d && !has; c++) has = roots[r0].first[(2 * mu) * P0.d + c] != 0 || roots[r0].first[(2 * mu + 1) * P0.d + c] != 0;
        if (!has) continue; auto it = std::lower_bound(P0.modes.begin(), P0.modes.end(), k * m); if (it == P0.modes.end() || *it != k * m) { ok = false; break; } int mu2 = (int)(it - P0.modes.begin());
        for (int c = 0; c < P0.d; c++) { x[(2 * mu2) * P0.d + c] = lam * roots[r0].first[(2 * mu) * P0.d + c]; x[(2 * mu2 + 1) * P0.d + c] = lam * roots[r0].first[(2 * mu + 1) * P0.d + c]; } }
      if (ok) { roots.emplace_back(x, roots[r0].second + "x" + std::to_string(k)); roms.push_back(roms[r0]); } } }
  ctx.w.resize(ctx.problem(cfg, cfg.Kmax)); long added = 0, dups = 0;
  Continuer C(cfg, cc, ctx, P0);
  C.on_solution = [&](std::vector<double>& x, const std::string& tag) {
    std::string why; Record rec; const Problem* Pp = &P0; std::vector<double> xs = x;
    if (!certify(cfg, Pp, xs, Symmetry(), ctx, rec, why)) return;
    long dup = outc.find_duplicate(rec, cfg.tol_inv, cfg.tol_dist);
    if (dup >= 0) { outc.absorb((size_t)dup, rec); dups++; std::printf("  [%s] → known solution id=%lld (A=%.9f)\n", tag.c_str(), (long long)outc.recs[dup].h.id, rec.h.action); return; }
    record_extras(cfg, rec, ctx); rec.sym = "continue:" + tag; outc.push(rec); added++; const Record& r = outc.recs.back();
    std::printf("+ id=%lld d=%d/%d N=%d K=%d A=%.9f E=%.6f morse=%d null=%d minsep=%.3f ret=%.1e cover=%d  (%s)\n", (long long)r.h.id, r.h.deff, r.h.d, r.h.N, r.h.K, r.h.action, r.h.energy, r.h.morse, r.h.nullity, r.h.minsep, r.h.ret_err, r.h.cover, tag.c_str());
    outc.save(outp);
  };
  if (skipped) std::printf("skipped %ld record(s) whose frame does not match --param %s\n", skipped, omode ? "omega" : "alpha");
  std::printf("continuation in %s ∈ [%.2f, %.2f] from %zu root(s), d=%d N=%d K=%d depth=%d → %s\n", omode ? "s (Ω = sΩ₀)" : "α", cc.lo, cc.hi, roots.size(), cfg.d, cfg.N, cfg.K, cc.depth, outp.c_str());
  for (size_t i = 0; i < roots.size(); i++) { C.branches = 0;
    if (omode) { C.pp = ContParam(); C.pp.frame(roms[i], cfg.d); C.pp.guard_lo = cc.lo - 1.0; C.pp.guard_hi = cc.hi + 1.0; }
    C.explore(roots[i].first, roots[i].second); }
  outc.save(outp); std::printf("continuation done: %ld new solutions, %ld already known; %s has %zu records\n", added, dups, outp.c_str(), outc.recs.size());
  return 0;
}

static int cmd_bench(const Args& a) {
  Config cfg; cfg.N = (int)a.num("N", 3); cfg.d = (int)a.num("d", 2); cfg.K = (int)a.num("K", 16); cfg.Kmax = (int)a.num("Kmax", 4 * cfg.K);
  cfg.omega = parse_omega(a.get("omega", ""), cfg.d);
  Ctx ctx; const Problem& P = ctx.problem(cfg, cfg.K); ctx.w.resize(P); la::Rng rng(1); std::vector<double> x, g(P.n), v(P.n), Hv(P.n), H; random_guess(P, rng, 4, 1.0, x);
  double kin, pot; action_grad(P, x.data(), nullptr, ctx.w, &kin, &pot); double lam = optimal_scale(kin, pot, 1.0); for (double& e : x) e *= lam; for (double& e : v) e = rng.normal();
  auto time = [&](const char* name, int reps, auto fn) { auto t0 = std::chrono::steady_clock::now(); for (int i = 0; i < reps; i++) fn(); double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / reps; std::printf("  %-28s %10.1f us\n", name, us); };
  std::printf("bench N=%d d=%d K=%d (modes %d, params %d, sample M=%d)\n", P.N, P.d, P.K, P.nm, P.n, P.M);
  time("action+gradient", 200, [&] { action_grad(P, x.data(), g.data(), ctx.w); });
  time("hessian-vector", 200, [&] { hessvec(P, x.data(), v.data(), Hv.data(), ctx.w); });
  time("full hessian", 20, [&] { hessian(P, x.data(), H, ctx.w); });
  time("sym_eig of hessian", 20, [&] { std::vector<double> Hc = H, w; la::sym_eig(P.n, Hc, w); });
  time("ODE return error (Taylor)", 20, [&] { return_error(P, x.data()); });
  time("full trial", 10, [&] { static uint64_t t = 0; run_trial(cfg, t++, ctx); });
  return 0;
}

#ifdef HAVE_MPFR

static int cmd_prove(const Args& a) {
  const std::string path = a.pos.at(0); Catalog cat = load_cat(path);
  const int digits = (int)a.num("digits", 40); int threads = (int)a.num("threads", (double)std::thread::hardware_concurrency()); if (threads < 1) threads = 1;
  const int order = (int)a.num("order", 1.15 * (digits + 4) + 4); const bool verbose = a.has("verbose"), write = a.has("write");
  const double radius = a.num("radius", std::pow(10.0, -digits / 2.0)), tol = std::pow(10.0, -(digits + 4));
  const mpfr_prec_t bits = (mpfr_prec_t)(digits * 3.3219280948873626 + 256); mpreal::set_default_prec(bits); ival::prec() = bits;
  std::vector<Record*> recs; for (auto& r : cat.recs) if (!a.has("id") || r.h.id == (long)a.num("id", 0)) recs.push_back(&r);
  if (recs.empty()) throw std::runtime_error("no such record");
  int proven = 0;
  for (Record* rp : recs) { Record& r = *rp; const int N = r.h.N, d = r.h.d, nd = N * d;
    std::printf("record %lld (N=%d d=%d deff=%d action=%.12f):", (long long)r.h.id, N, d, r.h.deff, r.h.action);
    if (r.h.alpha != 1.0) { std::printf(" alpha != 1 — not supported yet\n"); continue; }
    if (r.proven() > 0 && !a.has("force")) { std::printf(" already proven (radius %.1e)\n", r.proven()); proven++; continue; }
    std::vector<double> pos, vel; initial_state(N, d, r.mode_list(), r.coef.data(), r.omega(), pos, vel);
    if (const double* z = r.state()) { pos.assign(z, z + nd); vel.assign(z + nd, z + 2 * nd); }
    std::vector<double> Zd(pos); Zd.insert(Zd.end(), vel.begin(), vel.end()); int dp = d; Frame fr;
    if (r.omega()) { fr = frame_of(N, d, r.omega()); std::vector<double> Zr; fr.apply(N, Zd.data(), Zr); Zd.swap(Zr);
      if (verbose) { std::printf("\n  frame: %d planes, rates", fr.nplanes()); for (int i = 0; i < fr.nplanes(); i++) std::printf(" %.12g%s", fr.rate[i], fr.cls[i] == 0 ? "(fixed)" : fr.cls[i] == 1 ? "(pi)" : ""); std::printf("; %zu translations, %zu commuting rotations\n", fr.trans.size(), fr.rots.size()); } }
    else {
      double sc = 0; for (double v : Zd) sc = std::max(sc, std::fabs(v));
      std::vector<int> axes; for (int c = 0; c < d; c++) { double amp = 0; for (int k = 0; k < N; k++) amp = std::max(amp, std::max(std::fabs(pos[k * d + c]), std::fabs(vel[k * d + c]))); if (amp > 1e-7 * sc) axes.push_back(c); }
      dp = (int)axes.size(); Zd.assign(2 * (size_t)N * dp, 0.0);
      for (int k = 0; k < N; k++) for (int c = 0; c < dp; c++) { Zd[k * dp + c] = pos[k * d + axes[c]]; Zd[(size_t)N * dp + k * dp + c] = vel[k * d + axes[c]]; }
      fr = frame_of(N, dp, nullptr); }
    std::fflush(stdout); auto t0 = std::chrono::steady_clock::now();
    std::vector<mpreal> Z; double res = refine_state(N, dp, 1.0, Zd, fr, digits, threads, Z, verbose);
    if (verbose) std::printf("  refined to residual %.2e in %.1fs (%ld bits, order %d, radius %.1e, tol %.0e)\n", res, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), (long)bits, order, radius, tol);
    if (!(res < radius * 1e-3)) { std::printf(" refinement stalled at %.2e — not proven\n", res); continue; }
    Proof P = prove_state(N, dp, 1.0, Z, fr, radius, tol, order, threads, verbose);
    if (!P.ok) { std::printf(" NOT PROVEN (%s; Newton %.1e/%.1e, contraction %.2e, closure %.2e)\n", P.why.c_str(), P.newton, radius, P.kappa, P.closure); continue; }
    proven++;
    std::printf(" PROVEN in %.1fs\n  %s of N=%d bodies in R^%d exists with initial state within %.1e (max norm) of the refined state,\n"
                "  the only zero of the shooting map in the slice box Z0 + Q[-r, r]^m, r = %.1e; the bodies span R^%d: %s; not a relative equilibrium: %s.\n"
                "  Krawczyk: |Y F| = %.1e, contraction %.2e, closure %.2e; slice dim %d, %d gauge generators;\n"
                "  %ld validated Taylor steps (%ld rejected), h in [%.1e, %.1e], state width <= %.1e\n"
                "  energy in %s\n  action in %s\n", P.seconds, fr.rotating() ? "a relative choreography, rotating-frame period 2pi," : "an inertial choreography, period 2pi,", N, dp, P.hull, P.radius, dp, P.span ? "verified" : "NOT verified",
                P.nonrigid ? "verified" : "NOT verified", P.newton, P.kappa, P.closure, P.m, P.k, P.steps, P.rejects, P.hmin, P.hmax, P.maxwid, P.energy.str(digits / 2 + 4).c_str(), P.action.str(digits / 2 + 4).c_str());
    if (fr.rotating()) {
      bool periodic = true; std::printf("  certified frame rates (modulo N):");
      for (double w : fr.rate) { std::printf(" %a", w); periodic &= w == std::round(w); }
      std::printf("\n  Phi_(2pi)(Z) = G^N Z; G^N = I: %s. A common inertial curve is an additional condition.\n", periodic ? "verified from integer rates" : "not asserted");
    }
    if (write) {
      // Upgrade the payload without changing its initial state or frame. Old proof markers are
      // never trusted by this revision; the layout tag records that this run recomputed the proof.
      const auto sources = r.citations();
      const double* oldom = r.omega(); std::vector<double> om;
      if (oldom) om.assign(oldom, oldom + (size_t)d * d);
      r.extra.resize(Record::NEX + (size_t)d * d + 2 * (size_t)nd, 0.0);
      if (!om.empty()) std::copy(om.begin(), om.end(), r.extra.begin() + Record::NEX);
      std::copy(pos.begin(), pos.end(), r.extra.begin() + Record::NEX + (size_t)d * d);
      std::copy(vel.begin(), vel.end(), r.extra.begin() + Record::NEX + (size_t)d * d + nd);
      r.extra[5] = Record::PROOF_LAYOUT; r.extra[7] = P.radius;
      for (const auto& source : sources) r.add_citation(source);
      cat.save(path);
    }
  }
  if (write) std::printf("%d proven, catalogue saved\n", proven); else std::printf("%d of %zu proven\n", proven, recs.size());
  return proven == (int)recs.size() ? 0 : 1;
}

static int cmd_refine(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); const Record& r = rec_by_id(cat, (long)a.num("id", 0)); Problem P; std::vector<double> x; r.to_problem(P, x);
  int digits = (int)a.num("digits", 50), Kout = (int)a.num("K", P.K); std::string outp = a.get("out", "");
  int threads = (int)a.num("threads", (double)std::thread::hardware_concurrency()); if (threads < 1) threads = 1;
  const int N = P.N, d = P.d, nd = N * d, n2 = 2 * nd; const double* om = r.omega();

  std::vector<double> pos0, vel0; initial_state(P, x.data(), pos0, vel0);
  NBody<double> nbd(N, d, P.alpha, 22);
  if (const double* z = r.state()) { pos0.assign(z, z + nd); vel0.assign(z + nd, z + 2 * nd); }
  double sc0 = 1.0; for (double v : pos0) sc0 = std::max(sc0, std::fabs(v)); for (double v : vel0) sc0 = std::max(sc0, std::fabs(v));
  double stored = chore_residual(nbd, pos0, vel0, 1e-16, P.gshift());

  std::vector<double> Zd(pos0); Zd.insert(Zd.end(), vel0.begin(), vel0.end());
  ShootWork<double> Wd;
  double dres = shoot_newton(nbd, Zd, 1e-16, 20, 1e-13 * sc0, -18, -40, Wd, false, P.gshift(), threads);

  mpfr_prec_t bits = (mpfr_prec_t)(digits * 3.3219280948873626 + 96); mpreal::set_default_prec(bits);
  int order = (int)(1.15 * digits) + 6; double itol = std::pow(10.0, -(digits + 4));
  mpreal twopi = mpreal::pi() * 2;
  std::vector<mpreal> Z(n2); for (int i = 0; i < n2; i++) Z[i] = Zd[i];
  std::vector<mpreal> G, GN; const std::vector<mpreal>* Gp = nullptr;
  if (om) { omega_exp(d, om, twopi / N, G); omega_exp(d, om, twopi, GN); Gp = &G; }
  NBody<mpreal> nb(N, d, P.alpha, order);
  std::printf("refine record %lld (N=%d d=%d%s) to %d digits: %ld bits, Taylor order %d, %d threads\n",
              (long long)r.h.id, N, d, om ? ", rotating frame" : "", digits, (long)bits, order, threads);
  std::printf("  stored state %.3e (record claims ret_err %.1e), coefficients %.1e\n  after the double Newton: %.3e\n", stored, r.h.ret_err, r.coef_err(), dres);
  ShootWork<mpreal> W;
  double maxF = shoot_newton(nb, Z, itol, 25, std::pow(10.0, -digits), -(long)(bits / 3),
                             -(long)std::min<mpfr_prec_t>(bits / 2, 96), W, true, Gp, threads, 0.9, 24);

  mpreal E = nbody_energy(N, d, P.alpha, Z.data(), Z.data() + nd);
  double fullret = 0;
  { std::vector<mpreal> pf(Z.begin(), Z.begin() + nd), vf(Z.begin() + nd, Z.end());
    if (nb.integrate(pf, vf, twopi, itol) < 0) throw std::runtime_error("refinement full-period integration failed");
    for (int k = 0; k < N; k++) for (int c = 0; c < d; c++) {
      mpreal tp = Z[k * d + c], tv = Z[nd + k * d + c];
      if (!GN.empty()) { tp = mpreal(0); tv = mpreal(0);
        for (int b = 0; b < d; b++) { tp += GN[(size_t)c * d + b] * Z[k * d + b]; tv += GN[(size_t)c * d + b] * Z[nd + k * d + b]; } }
      fullret = std::max(fullret, std::max(std::fabs(to_double(pf[k * d + c] - tp)), std::fabs(to_double(vf[k * d + c] - tv)))); } }

  const int Mseg = std::max(8, (std::max(64, 16 * Kout) + N - 1) / N), Ms = Mseg * N;
  std::vector<mpreal> ts(Mseg); for (int i = 0; i < Mseg; i++) ts[i] = twopi * i / Ms;
  std::vector<mpreal> ps(Z.begin(), Z.begin() + nd), vs(Z.begin() + nd, Z.end()), seg;
  if (nb.integrate(ps, vs, twopi / N, itol, &ts, &seg) < 0) throw std::runtime_error("refinement sampling integration failed");
  mpreal A(0);
  for (int i = 0; i < Mseg; i++) { const mpreal* sp = &seg[(size_t)i * 2 * nd];
    for (int b = 0; b < N; b++) { mpreal ke(0); for (int c = 0; c < d; c++) ke += sp[nd + b * d + c] * sp[nd + b * d + c]; A += ke * 0.5; }
    for (int b = 0; b < N; b++) for (int l = b + 1; l < N; l++) { mpreal r2(0);
      for (int c = 0; c < d; c++) { mpreal df = sp[b * d + c] - sp[l * d + c]; r2 += df * df; }
      A += P.alpha == 1.0 ? 1 / sqrt(r2) : pow(r2, mpreal(-P.alpha / 2)); } }
  A = A * twopi / Ms;
  double virial = P.alpha == 1.0 ? std::fabs(to_double(A) * N + 6 * PI * to_double(E)) : 0.0;
  std::vector<mpreal> qs((size_t)Ms * d), Rj, Rs, Rn;
  if (om) { omega_exp(d, om, -twopi / Ms, Rs); Rj.assign((size_t)d * d, mpreal(0)); Rn.assign((size_t)d * d, mpreal(0));
            for (int i = 0; i < d; i++) Rj[(size_t)i * d + i] = mpreal(1); }
  for (int i = 0; i < Mseg; i++) {
    for (int j = 0; j < N; j++) { const mpreal* Q = &seg[(size_t)i * 2 * nd + j * d];
      for (int c = 0; c < d; c++) { if (!om) { qs[((size_t)j * Mseg + i) * d + c] = Q[c]; continue; }
        mpreal t(0); for (int b = 0; b < d; b++) t += Rj[(size_t)c * d + b] * Q[b]; qs[((size_t)j * Mseg + i) * d + c] = t; } }
    if (om) { for (int u = 0; u < d; u++) for (int l = 0; l < d; l++) { mpreal t(0);
                for (int b = 0; b < d; b++) t += Rj[(size_t)u * d + b] * Rs[(size_t)b * d + l]; Rn[(size_t)u * d + l] = t; }
              Rj.swap(Rn); } }
  seg.clear(); seg.shrink_to_fit();
  std::vector<mpreal> cosj(Ms), sinj(Ms); for (int j = 0; j < Ms; j++) { mpreal th = twopi * j / Ms; cosj[j] = cos(th); sinj[j] = sin(th); }
  std::vector<mpreal> cm((size_t)(Kout + 1) * d), sm((size_t)(Kout + 1) * d);
  for (int m = 0; m <= Kout; m++) for (int c = 0; c < d; c++) { mpreal sc(0), ss(0);
    for (int j = 0; j < Ms; j++) { int idx = (int)(((long)m * j) % Ms); sc += qs[(size_t)j * d + c] * cosj[idx]; ss += qs[(size_t)j * d + c] * sinj[idx]; }
    cm[(size_t)m * d + c] = sc * (m ? 2 : 1) / Ms; sm[(size_t)m * d + c] = ss * (m ? 2 : 1) / Ms; }
  double amax = 0, atail = 0;
  for (int m = 1; m <= Kout; m++) { double mg = 0;
    for (int c = 0; c < d; c++) mg = std::max(mg, std::max(std::fabs(to_double(cm[(size_t)m * d + c])), std::fabs(to_double(sm[(size_t)m * d + c]))));
    amax = std::max(amax, mg); if (m > Kout - 8) atail = std::max(atail, mg); }
  Problem Pt; Pt.init(N, d, Kout, 0, P.alpha); if (om) Pt.set_omega(std::vector<double>(om, om + (size_t)d * d));
  std::vector<double> xt(Pt.n, 0.0);
  for (int mu = 0; mu < Pt.nm; mu++) { int m = Pt.modes[mu]; if (m > Kout) break;
    for (int c = 0; c < d; c++) { xt[(2 * mu) * d + c] = to_double(cm[(size_t)m * d + c]); xt[(2 * mu + 1) * d + c] = to_double(sm[(size_t)m * d + c]); } }
  std::vector<double> pt, vt; initial_state(Pt, xt.data(), pt, vt);
  double trunc = chore_residual(nbd, pt, vt, 1e-16, Pt.gshift());

  FILE* fo = outp.empty() ? stdout : std::fopen(outp.c_str(), "w"); if (!fo) throw std::runtime_error("cannot write " + outp);
  std::fprintf(fo, "# hyperchoreography refined record %lld: N=%d d=%d period=2*pi alpha=%g digits=%d%s\n", (long long)r.h.id, N, d, P.alpha, digits, om ? " rotating_frame=1" : "");
  std::fprintf(fo, "stored_shift_residual %.3e\nshift_residual %.3e\nfull_period_return %.3e\nvirial_defect %.3e\nfourier_tail %.3e\nfourier_truncation_residual %.3e\nenergy %s\naction_per_body %s\n",
               stored, maxF, fullret, virial, amax > 0 ? atail / amax : 0.0, trunc, E.str(digits).c_str(), A.str(digits).c_str());
  if (om) { std::fprintf(fo, "# Omega (rotating frame, canonical axes), row-major %dx%d\n", d, d);
    for (int i = 0; i < d; i++) { std::fprintf(fo, "omega"); for (int j = 0; j < d; j++) std::fprintf(fo, " %+.17g", om[(size_t)i * d + j]); std::fprintf(fo, "\n"); }
    std::fprintf(fo, "# q_j(t) = exp(Omega t) q(t + 2 pi j/N); the states below are inertial (Q, Q-dot) at t = 0\n"); }
  for (int k = 0; k < N; k++) { std::fprintf(fo, "body %d q", k); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", Z[k * d + c].str(digits).c_str()); std::fprintf(fo, "\nbody %d v", k); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", Z[nd + k * d + c].str(digits).c_str()); std::fprintf(fo, "\n"); }
  std::fprintf(fo, "# Fourier coefficients of %s: mode m, c_m[0..d-1], s_m[0..d-1]   (from %d dense samples)\n", om ? "the rotating-frame loop q(t) = exp(-Omega t) Q_0(t)" : "body 0", Ms);
  for (int m = 0; m <= Kout; m++) { double mag = 0; for (int c = 0; c < d; c++) mag = std::max(mag, std::max(std::fabs(to_double(cm[(size_t)m * d + c])), std::fabs(to_double(sm[(size_t)m * d + c]))));
    if (m && mag < 1e-300) continue; std::fprintf(fo, "mode %d", m); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", cm[(size_t)m * d + c].str(digits).c_str()); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", sm[(size_t)m * d + c].str(digits).c_str()); std::fprintf(fo, "\n"); }
  if (fo != stdout) { std::fclose(fo); std::printf("  residual %.3e (was %.3e), full-period %.3e, virial %.1e, K=%d tail %.1e, truncation %.1e\n  energy %s\n  wrote %s\n",
                                                   maxF, stored, fullret, virial, Kout, amax > 0 ? atail / amax : 0.0, trunc, E.str(20).c_str(), outp.c_str()); }
  return 0;
}
#endif

int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 1; }
  std::string cmd = argv[1]; Args a(argc, argv);
#ifdef HAVE_MPFR
  mpreal::set_default_prec((mpfr_prec_t)(MP_DIGITS * 3.3219280948873626 + 96));
#endif
  try {
    if (cmd == "search") return cmd_search(a);
    if (cmd == "list") return cmd_list(a);
    if (cmd == "show") return cmd_show(a);
    if (cmd == "export") return cmd_export(a);
    if (cmd == "verify") return cmd_verify(a);
    if (cmd == "merge") return cmd_merge(a);
    if (cmd == "extras") return cmd_extras(a);
    if (cmd == "symmetry") return cmd_symmetry(a);
    if (cmd == "bench") return cmd_bench(a);
    if (cmd == "continue") return cmd_continue(a);
#ifdef HAVE_MPFR
    if (cmd == "refine") return cmd_refine(a);
    if (cmd == "prove") return cmd_prove(a);
#else
    if (cmd == "prove") { std::fprintf(stderr, "built without MPFR (make without NOMPFR=1)\n"); return 1; }
    if (cmd == "refine") { std::fprintf(stderr, "built without MPFR (make without NOMPFR=1)\n"); return 1; }
#endif
  } catch (std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 1; }
  usage(); return 1;
}
