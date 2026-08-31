// CLI: search | continue | list | show | export | verify | refine | merge | extras | symmetry | bench
#include "search.hpp"
#include "continue.hpp"
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
    "                [--starts random,torus,vertical,hyper,fano,kick] [--K-index 48] [--Ms 2048 --Kout-max 512 --shoot-tol 1e-12 --ret-double 1e-4]\n"
    "                [--omega \"w1,w2,...\" | --omega su:w1,... | --omega g2:p,q[,r,...]]   rotating frame: q_j(t) = exp(Omega t) q(t + 2 pi j/N)\n"
    "                [--tol-inv 1e-4 --tol-dist 1e-3 --checkpoint 30 --ret-reject 1e-1]   (resumes if catalog/state exist)\n"
    "  hyperchoreography list    catalog.bin [--N n] [--deff k] [--min-deff k] [--sort action|id|hits|twist|rigid]\n"
    "  hyperchoreography show    catalog.bin --id i                (JSON dump of one record)\n"
    "  hyperchoreography export  catalog.bin --id i [--samples 720] [--out curve.csv]   (body positions over one period)\n"
    "  hyperchoreography verify  catalog.bin --id i                (double-precision Taylor integration checks)\n"
    "  hyperchoreography refine  catalog.bin --id i --digits 60 [--K 64] [--out refined.txt]   (MPFR shooting Newton)\n"
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
  cfg.sym = named_sym(a.get("sym", "none"), cfg.d); cfg.threads = (int)a.num("threads", std::thread::hardware_concurrency()); cfg.seed = (uint64_t)a.num("seed", 1);
  cfg.trials = (long)a.num("trials", (double)LONG_MAX); cfg.minutes = a.num("minutes", 1e30); cfg.out = a.get("out", "catalog.bin");
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
      std::lock_guard<std::mutex> lk(mu);
      if (!o.ok) { failed++; reasons[o.why]++; continue; }
      double dist; long dup = cat.find_duplicate(o.rec, cfg.tol_inv, cfg.tol_dist, &dist);
      if (dup >= 0) { if (cat.absorb((size_t)dup, o.rec)) record_extras(cfg, cat.recs[dup], ctx); dups++; continue; }
      record_extras(cfg, o.rec, ctx);                      // new records only
      o.rec.h.id = -1; size_t idx = cat.push(o.rec); found++; const Record& r = cat.recs[idx];
      std::printf("+ id=%lld d=%d/%d N=%d K=%d A=%.9f E=%.6f morse=%d null=%d minsep=%.3f ret=%.1e cover=%d sym=\"%s\" (trial %llu, %.2fs)\n", (long long)r.h.id, r.h.deff, r.h.d, r.h.N, r.h.K, r.h.action, r.h.energy, r.h.morse, r.h.nullity, r.h.minsep, r.h.ret_err, r.h.cover, r.sym.c_str(), (unsigned long long)tr, r.h.secs);
      std::fflush(stdout);
    }
    active--;
  });
  double last_ck = 0;
  const int64_t base_done = st.trials_done; const double base_elapsed = st.elapsed;
  auto checkpoint = [&] { std::lock_guard<std::mutex> lk(mu); cat.save(cfg.out); st.next_trial = (int64_t)counter.load(); st.trials_done = base_done + done.load(); st.found = (int64_t)cat.recs.size(); st.elapsed = base_elapsed + elapsed(); st.save(cfg.out + ".state"); };
  while (active.load() > 0) {
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    if (elapsed() - last_ck >= cfg.checkpoint_secs) { last_ck = elapsed(); checkpoint();
      std::printf("[%6.0fs] trials %ld (%.1f/s) | unique %zu (+%ld) | duplicates %ld | failed %ld\n", elapsed(), done.load(), done.load() / std::max(1e-9, elapsed()), cat.recs.size(), found.load(), dups.load(), failed.load()); std::fflush(stdout); }
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
  for (auto r : rs) std::printf("%5lld %2d/%-2d %2d %4d %14.9f %12.7f %5d %4d %7.4f %8.1e %10.4g %6.3f %8.2e %5d %3d %s\n", (long long)r->h.id, r->h.deff, r->h.d, r->h.N, r->h.K, r->h.action, r->h.energy, r->h.morse, r->h.nullity, r->h.minsep, r->h.ret_err, r->twist(), r->extra.size() > 1 ? r->extra[1] : 0.0, r->rigid(), r->h.hits, r->h.cover, r->sym.c_str());
  std::printf("%zu records\n", rs.size());
  return 0;
}

static int cmd_show(const Args& a) { Catalog cat = load_cat(a.pos.at(0)); std::printf("%s\n", rec_by_id(cat, (long)a.num("id", 0)).to_json().c_str()); return 0; }

static int cmd_export(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); const Record& r = rec_by_id(cat, (long)a.num("id", 0)); int S = (int)a.num("samples", 720);
  std::string out = a.get("out", "curve_" + std::to_string(r.h.id) + ".csv"); FILE* f = std::fopen(out.c_str(), "w"); if (!f) throw std::runtime_error("cannot write " + out);
  const int N = r.h.N, d = r.h.d; std::fprintf(f, "t");
  for (int k = 0; k < N; k++) for (int c = 0; c < d; c++)                       // x,y,z,w then x4…
    if (c < 4) std::fprintf(f, ",q%d_%c", k, "xyzw"[c]); else std::fprintf(f, ",q%d_x%d", k, c);
  std::fprintf(f, "\n");
  const double* om = r.omega(); std::vector<double> Rt, Ai;
  for (int j = 0; j < S; j++) { double t = 2 * PI * j / S; std::fprintf(f, "%.10f", t);
    if (om) { Ai.assign(om, om + (size_t)d * d); for (double& e : Ai) e *= t; la::expm_skew(d, Ai, Rt); }
    for (int k = 0; k < N; k++) { double tk = t + 2 * PI * k / N; std::vector<double> q(d, 0.0);
      for (size_t mu = 0; mu < r.modes.size(); mu++) { double c = std::cos(r.modes[mu] * tk), s = std::sin(r.modes[mu] * tk); for (int c2 = 0; c2 < d; c2++) q[c2] += c * r.coef[mu * 2 * d + c2] + s * r.coef[mu * 2 * d + d + c2]; }
      if (om) { std::vector<double> qi(d, 0.0);                                   // inertial: exp(Ωt) q(t + 2πk/N)
        for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) qi[a] += Rt[(size_t)a * d + b] * q[b]; q.swap(qi); }
      for (int c2 = 0; c2 < d; c2++) std::fprintf(f, ",%.12f", q[c2]); }
    std::fprintf(f, "\n"); }
  std::fclose(f); std::printf("wrote %s (%d samples, %d bodies, %dD)\n", out.c_str(), S, N, d); return 0;
}

static int cmd_verify(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); const Record& r = rec_by_id(cat, (long)a.num("id", 0)); Problem P; std::vector<double> x; r.to_problem(P, x);
  std::vector<double> pos, vel; initial_state(P, x.data(), pos, vel); NBody<double> nb(P.N, P.d, P.alpha, 22);
  double E0 = nbody_energy(P.N, P.d, P.alpha, pos.data(), vel.data());
  double res = chore_residual(nb, pos, vel, 1e-16, P.gshift());
  std::vector<double> p = pos, v = vel; int steps = nb.integrate(p, v, 2 * PI, 1e-16); double E1 = nbody_energy(P.N, P.d, P.alpha, p.data(), v.data()), ret = 0;
  std::vector<double> GN;                                    // after one period a twisted orbit closes up to exp(2πΩ)
  if (!P.Om.empty()) { std::vector<double> A(P.Om); for (double& e : A) e *= 2 * PI; la::expm_skew(P.d, A, GN); }
  for (int k = 0; k < P.N; k++) for (int c = 0; c < P.d; c++) { double tp = pos[k * P.d + c], tv = vel[k * P.d + c];
    if (!GN.empty()) { tp = tv = 0; for (int b = 0; b < P.d; b++) { tp += GN[(size_t)c * P.d + b] * pos[k * P.d + b]; tv += GN[(size_t)c * P.d + b] * vel[k * P.d + b]; } }
    ret = std::max(ret, std::max(std::fabs(p[k * P.d + c] - tp), std::fabs(v[k * P.d + c] - tv))); }
  std::printf("record %lld: d=%d (deff=%d) N=%d K=%d action=%.12f energy=%.12f\n", (long long)r.h.id, r.h.d, r.h.deff, r.h.N, r.h.K, r.h.action, r.h.energy);
  if (!P.Om.empty()) std::printf("  rotating frame: Omega stored with the record\n");
  std::printf("  shift residual |Phi_{T/N}(Z) - G S Z|  = %.3e\n  full-period return |Phi_T(Z) - G^N Z| = %.3e  (%d Taylor steps)\n  energy drift over one period          = %.3e\n  energy from initial state             = %.12f\n", res, ret, steps, E1 - E0, E0);
  if (std::fabs(E0 - r.h.energy) > 1e-6 * std::max(1.0, std::fabs(r.h.energy))) std::printf("  ** stale record: reconstructed energy %.9g disagrees with stored %.9g **\n", E0, r.h.energy);
  std::printf("  initial conditions (t=0):\n"); for (int k = 0; k < P.N; k++) { std::printf("   body %d  q =", k); for (int c = 0; c < P.d; c++) std::printf(" %+.15f", pos[k * P.d + c]); std::printf("   v ="); for (int c = 0; c < P.d; c++) std::printf(" %+.15f", vel[k * P.d + c]); std::printf("\n"); }
  return 0;
}

// --min-rigid/--min-deff re-apply a gate to catalogues already on disk; extra[4] is the stored defect.
static int cmd_merge(const Args& a) {
  if (a.pos.size() < 2) { usage(); return 1; }
  Catalog out; try { out.load(a.pos[0]); } catch (...) {}
  double tol_inv = a.num("tol-inv", 1e-4), tol_dist = a.num("tol-dist", 1e-3);
  double min_rigid = a.num("min-rigid", 0.0); int min_deff = (int)a.num("min-deff", 0);
  long added = 0, merged = 0, dropped = 0;
  auto keep = [&](const Record& r) {
    if (r.h.deff < min_deff) return false;
    if (min_rigid > 0 && (int)r.extra.size() > 4 && r.extra[4] < min_rigid) return false;
    return true; };
  if (min_rigid > 0 || min_deff > 0) {                       // re-gate what is already in the destination
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
    if (r.extra.empty()) r.extra.assign(Record::NEX, 0.0);                 // legacy record, no extras yet
    Config cfg; cfg.N = r.h.N; cfg.d = r.h.d; cfg.K = r.h.K; cfg.K_index = Ki; cfg.minsep = 0;
    Ctx ctx;                                  // fresh: Ctx caches problems by K alone, and (N, d) vary here
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
    int nsh = (cont & 1) ? 0 : 1, nrv = 0;                        // the identity is always among the shifts
    for (auto& e : g) {
      if (e.eps > 0) { if (e.p == 0) continue; nsh++; if (!gen || (double)e.p / e.q < (double)gen->p / gen->q) gen = &e; }
      else { nrv++; if (!rev) rev = &e; }
    }
    std::string t = gen ? gen->text : "";
    if (rev) { size_t sp = rev->text.find(' ');                   // theta is a time origin; re-origin it to 0
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
      const double* om = r.omega();                                             // Ω-mode relaxes a frame; α-mode must not silently drop one
      if (omode ? !om : om != nullptr) { skipped++; continue; }
      const Problem& P = ctx.problem(cfg, cfg.K); std::vector<double> x(P.n, 0.0);
      for (size_t k = 0; k < r.modes.size(); k++) { int m = r.modes[k]; if (m > P.K) break; auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) continue; int mu = (int)(it - P.modes.begin());
        for (int c = 0; c < r.h.d; c++) { x[(2 * mu) * P.d + c] = r.coef[k * 2 * r.h.d + c]; x[(2 * mu + 1) * P.d + c] = r.coef[k * 2 * r.h.d + r.h.d + c]; } }
      std::vector<double> O; if (omode) { O.assign((size_t)cfg.d * cfg.d, 0.0);      // the record's frame, embedded in d
        for (int i = 0; i < r.h.d; i++) for (int j = 0; j < r.h.d; j++) O[(size_t)i * cfg.d + j] = om[(size_t)i * r.h.d + j]; }
      roots.emplace_back(x, "id" + std::to_string(r.h.id)); roms.push_back(O); }
  }
  // k-fold covers of the roots (the transverse resonances live on the covers)
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
static int cmd_refine(const Args& a) {
  Catalog cat = load_cat(a.pos.at(0)); const Record& r = rec_by_id(cat, (long)a.num("id", 0)); Problem P; std::vector<double> x; r.to_problem(P, x);
  if (r.omega()) throw std::runtime_error("refine does not support rotating-frame records yet (record has Omega)");
  int digits = (int)a.num("digits", 50), Kout = (int)a.num("K", P.K); std::string outp = a.get("out", "");
  mpfr_prec_t bits = (mpfr_prec_t)(digits * 3.3219280948873626 + 96); mpreal::set_default_prec(bits);
  const int N = P.N, d = P.d, nd = N * d, n2 = 2 * nd; int order = (int)(1.15 * digits) + 6; double itol = std::pow(10.0, -(digits + 4));
  std::vector<double> pos0, vel0; initial_state(P, x.data(), pos0, vel0);
  std::vector<mpreal> Z(n2); for (int i = 0; i < nd; i++) { Z[i] = pos0[i]; Z[nd + i] = vel0[i]; }
  NBody<mpreal> nb(N, d, P.alpha, order);
  std::printf("refine record %lld (N=%d d=%d) to %d digits: %ld bits, Taylor order %d\n", (long long)r.h.id, N, d, digits, (long)bits, order);
  ShootWork<mpreal> W; double maxF = shoot_newton(nb, Z, itol, 25, std::pow(10.0, -digits), -(long)(bits / 3), -(long)(bits / 2), W, true);
  (void)n2;
  // diagnostics and Fourier coefficients from a dense-output period
  mpreal E = nbody_energy(N, d, P.alpha, Z.data(), Z.data() + nd);
  int Ms = std::max(64, 16 * Kout); std::vector<mpreal> ts(Ms), samp; mpreal twopi = mpreal::pi() * 2;
  for (int j = 0; j < Ms; j++) ts[j] = twopi * j / Ms;
  std::vector<mpreal> p(Z.begin(), Z.begin() + nd), v(Z.begin() + nd, Z.end()); nb.integrate(p, v, twopi, itol, &ts, &samp);
  double fullret = 0; for (int i = 0; i < nd; i++) fullret = std::max(fullret, std::max(std::fabs(to_double(p[i] - Z[i])), std::fabs(to_double(v[i] - Z[nd + i]))));
  std::vector<mpreal> cosj(Ms), sinj(Ms); for (int j = 0; j < Ms; j++) { cosj[j] = cos(ts[j]); sinj[j] = sin(ts[j]); }
  std::vector<mpreal> cm((size_t)(Kout + 1) * d), sm((size_t)(Kout + 1) * d);
  for (int m = 0; m <= Kout; m++) for (int c = 0; c < d; c++) { mpreal sc(0), ss(0);
    for (int j = 0; j < Ms; j++) { int idx = (int)(((long)m * j) % Ms); sc += samp[(size_t)j * 2 * nd + c] * cosj[idx]; ss += samp[(size_t)j * 2 * nd + c] * sinj[idx]; }
    cm[(size_t)m * d + c] = sc * (m ? 2 : 1) / Ms; sm[(size_t)m * d + c] = ss * (m ? 2 : 1) / Ms; }
  mpreal A(0);   // action per body
  for (int j = 0; j < Ms; j++) { const mpreal* s = &samp[(size_t)j * 2 * nd]; mpreal ke(0); for (int c = 0; c < d; c++) ke += s[nd + c] * s[nd + c];
    mpreal pe(0); for (int k = 1; k < N; k++) { mpreal r2(0); for (int c = 0; c < d; c++) { mpreal df = s[c] - s[k * d + c]; r2 += df * df; } pe += P.alpha == 1.0 ? 1 / sqrt(r2) : pow(r2, mpreal(-P.alpha / 2)); }
    A += (ke + pe) * 0.5; }
  A = A * twopi / Ms;
  FILE* fo = outp.empty() ? stdout : std::fopen(outp.c_str(), "w"); if (!fo) throw std::runtime_error("cannot write " + outp);
  std::fprintf(fo, "# hyperchoreography refined record %lld: N=%d d=%d period=2*pi alpha=%g digits=%d\n", (long long)r.h.id, N, d, P.alpha, digits);
  std::fprintf(fo, "shift_residual %.3e\nfull_period_return %.3e\nenergy %s\naction_per_body %s\n", maxF, fullret, E.str(digits).c_str(), A.str(digits).c_str());
  for (int k = 0; k < N; k++) { std::fprintf(fo, "body %d q", k); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", Z[k * d + c].str(digits).c_str()); std::fprintf(fo, "\nbody %d v", k); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", Z[nd + k * d + c].str(digits).c_str()); std::fprintf(fo, "\n"); }
  std::fprintf(fo, "# Fourier coefficients of body 0: mode m, c_m[0..d-1], s_m[0..d-1]   (from %d dense samples)\n", Ms);
  for (int m = 0; m <= Kout; m++) { double mag = 0; for (int c = 0; c < d; c++) mag = std::max(mag, std::max(std::fabs(to_double(cm[(size_t)m * d + c])), std::fabs(to_double(sm[(size_t)m * d + c]))));
    if (m && mag < 1e-300) continue; std::fprintf(fo, "mode %d", m); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", cm[(size_t)m * d + c].str(digits).c_str()); for (int c = 0; c < d; c++) std::fprintf(fo, " %s", sm[(size_t)m * d + c].str(digits).c_str()); std::fprintf(fo, "\n"); }
  if (fo != stdout) { std::fclose(fo); std::printf("residual %.3e, full-period return %.3e, energy %s\nwrote %s\n", maxF, fullret, E.str(20).c_str(), outp.c_str()); }
  return 0;
}
#endif

int main(int argc, char** argv) {
  if (argc < 2) { usage(); return 1; }
  std::string cmd = argv[1]; Args a(argc, argv);
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
#else
    if (cmd == "refine") { std::fprintf(stderr, "built without MPFR (make without NOMPFR=1)\n"); return 1; }
#endif
  } catch (std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 1; }
  usage(); return 1;
}
