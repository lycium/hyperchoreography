// One search trial: start → L-BFGS → Newton–LM → ODE validation → shooting certification → Fourier re-extraction → record.
#pragma once
#include "catalog.hpp"
#include "optim.hpp"
#include "taylor.hpp"
#include "calib.hpp"
#include <map>
#include <chrono>
#include <climits>
#include <memory>
#include <sstream>
#include <algorithm>

struct Config {
  int N = 3, d = 2, K = 16, Kmax = 64;
  double alpha_start = 1.0; int alpha_steps = 8;
  std::string sym = "none";                    // none | random | generator DSL
  std::string phase1 = "mixed";                // action | gradnorm | mixed
  int threads = 0; uint64_t seed = 1; long trials = LONG_MAX; double minutes = 1e30;
  std::string out = "catalog.bin";
  int lbfgs_min = 20, lbfgs_max = 400, newton_iters = 60; double gtol = 1e-10, ret_tol = 1e-8, ret_reject = 1e-1;
  int K0min = 2, K0max = 6; double minsep = 2e-3; int min_deff = 1; double min_rigid = 1e-4;   // rigid clusters end at 4e-6, the first real orbit is at 1.7e-2
  const std::vector<Record>* seeds = nullptr; double kick_min = 0.02, kick_max = 0.5;
  int Ms = 2048, Kout_max = 512, Kout_grow = 4096; double shoot_tol = 1e-12, ret_double = 1e-4, mpfr_gate = 1e-12;
  std::string starts = "random,torus,vertical";                                        // comma list; kick needs seeds
  int K_index = 48;                                                                    // modes for the Morse index
  std::vector<double> omega; std::string omega_text;                                   // rotating frame Ω ∈ so(d)
  bool verbose = false; int checkpoint_secs = 30; double tol_inv = 1e-4, tol_dist = 1e-3;
};

// --omega "w1,w2,…": rates in the successive coordinate planes (1,2),(3,4),…;  "su:w1,…,w_{n−1}" closes the
// list with w_n = −Σw, the condition for exp(Ωt) to preserve the special-Lagrangian form of the SU(n) rung;
// "g2:p,q" (d = 7) is the maximal torus of g2 — rates (p, q, p+q) on the Fourier planes of the Fano 7-cycle,
// which is the same "rates sum to zero" condition seen through R⁷ = R ⊕ C³.
inline std::vector<double> parse_omega(const std::string& t, int d) {
  if (t.empty() || t == "none") return {};
  if (t.rfind("su:", 0) == 0) {
    std::vector<double> w; std::stringstream ws(t.substr(3)); std::string e;
    while (std::getline(ws, e, ',')) w.push_back(std::stod(e));
    return su_omega(d, w);
  }
  if (t.rfind("g2:", 0) == 0) {
    if (d < 7) throw std::runtime_error("--omega g2: needs d >= 7");
    std::vector<double> w; std::stringstream ss(t.substr(3)); std::string e;
    while (std::getline(ss, e, ',')) w.push_back(e.empty() ? 0.0 : std::stod(e));
    w.resize(std::max<size_t>(w.size(), 2), 0.0);            // "g2:p" is (p, 0, p)
    return g2_omega(w[0], w[1], d, std::vector<double>(w.begin() + 2, w.end()));
  }
  std::vector<double> Om((size_t)d * d, 0.0); std::stringstream ss(t); std::string e;
  for (int p = 0; std::getline(ss, e, ','); p++) {
    if (2 * p + 1 >= d) throw std::runtime_error("--omega: more rates than rotation planes");
    double w = std::stod(e); Om[(size_t)(2 * p) * d + 2 * p + 1] = -w; Om[(size_t)(2 * p + 1) * d + 2 * p] = w;
  }
  return Om;
}

// Named symmetry classes: twisted choreographies q(t + 2πp/q) = σ q(t) with σ of order q. "fano:p" is the
// d = 7 Fano 7-cycle (which also preserves φ); "cyc:p" cycles the ⌊d/2⌋ coordinate planes in
// any dimension. Both have trivial core, so unlike a group with a fixed rotation they do not cap deff — their
// value is cost: the fixed subspace has r ≈ n/q parameters and the Hessian eigendecomposition is O(q³) cheaper.
inline std::string named_sym(const std::string& t, int d) {
  if (t.rfind("fano:", 0) == 0) {
    if (d != 7) throw std::runtime_error("--sym fano: needs d = 7");
    return fano_sym(std::stoi(t.substr(5)));
  }
  if (t.rfind("cyc:", 0) == 0) {
    const int n = d / 2;
    if (n < 2) throw std::runtime_error("--sym cyc: needs d >= 4");
    int p = std::stoi(t.substr(4));
    std::string s = "t+" + std::to_string(((p % n) + n) % n) + "/" + std::to_string(n) + " s[";
    for (int k = 0; k < d; k++) s += (k ? "," : "") + std::to_string((k < 2 * n ? (k + 2) % (2 * n) : k) + 1);
    return s + "]";
  }
  return t;
}

struct Ctx {                                  // per-thread state
  Work w;
  std::map<int, Problem> probs;               // by K
  std::map<std::pair<int, std::string>, std::pair<std::vector<double>, int>> bases;
  std::unique_ptr<NBody<double>> nb;
  ShootWork<double> sw; std::vector<double> cosT, sinT; int Ms = 0;
  std::vector<double> H, g, eig, pos, vel, Z, samp, xd, tmp;
  const Problem& problem(const Config& cfg, int K) {
    auto it = probs.find(K);
    if (it == probs.end()) { Problem P; P.init(cfg.N, cfg.d, K, 0, 1.0); P.minsep = cfg.minsep; P.set_omega(cfg.omega); it = probs.emplace(K, std::move(P)).first; }
    return it->second;
  }
  const std::pair<std::vector<double>, int>& basis(const Problem& P, const Symmetry& S) {
    auto key = std::make_pair(P.K, S.text);
    auto it = bases.find(key);
    if (it == bases.end()) { int r; auto B = S.basis(P, r); it = bases.emplace(key, std::make_pair(std::move(B), r)).first; }
    return it->second;
  }
  NBody<double>& integrator(const Config& cfg) { if (!nb) nb = std::make_unique<NBody<double>>(cfg.N, cfg.d, 1.0, 22); return *nb; }
  void tables(int Ms_) { if (Ms == Ms_) return; Ms = Ms_; cosT.resize(Ms); sinT.resize(Ms); for (int j = 0; j < Ms; j++) { cosT[j] = std::cos(2 * PI * j / Ms); sinT[j] = std::sin(2 * PI * j / Ms); } }
};

// state of body k at t = 0, from (modes, coefficients) — no Problem, whose tables are tens of MB at K = 512
inline void initial_state(int N, int d, const std::vector<int>& modes, const double* coef, const double* Om,
                          std::vector<double>& pos, std::vector<double>& vel) {
  const int nm = (int)modes.size(); pos.assign((size_t)N * d, 0.0); vel.assign((size_t)N * d, 0.0);
  for (int k = 0; k < N; k++) { const double t = 2.0 * PI * k / N;
    for (int mu = 0; mu < nm; mu++) { const double m = modes[mu], c = std::cos(m * t), s = std::sin(m * t);
      const double* cm = coef + (size_t)mu * 2 * d; const double* sm = cm + d;
      for (int a = 0; a < d; a++) { pos[k * d + a] += c * cm[a] + s * sm[a]; vel[k * d + a] += m * (c * sm[a] - s * cm[a]); } } }
  if (Om) for (int k = 0; k < N; k++) for (int a = 0; a < d; a++)                          // Q̇ = q̇ + Ωq
    for (int b = 0; b < d; b++) vel[k * d + a] += Om[(size_t)a * d + b] * pos[k * d + b];
}
inline void initial_state(const Problem& P, const double* x, std::vector<double>& pos, std::vector<double>& vel) {
  initial_state(P.N, P.d, P.modes, x, P.Om.empty() ? nullptr : P.Om.data(), pos, vel);
}
// relative shift residual of the Fourier loop under the ODE
inline double return_error(const Problem& P, const double* x, NBody<double>& nb) {
  std::vector<double> pos, vel; initial_state(P, x, pos, vel);
  double sc = 1.0; for (double v : pos) sc = std::max(sc, std::fabs(v)); for (double v : vel) sc = std::max(sc, std::fabs(v));
  return chore_residual(nb, pos, vel, 1e-16, P.gshift()) / sc;
}
inline double return_error(const Problem& P, const double* x) { NBody<double> nb(P.N, P.d, P.alpha, 22); return return_error(P, x, nb); }

// A stored record re-checked against the ODE: the T/N shift residual of the certified state
// and of the coefficients, and optionally the full-period return, all relative to the loop's scale.
inline void record_residuals(const Record& r, double& state_res, double& coef_res, double* period_ret = nullptr) {
  const int N = r.h.N, d = r.h.d, nd = N * d;
  const double* om = r.omega();
  NBody<double> nb(N, d, r.h.alpha, 22);
  std::vector<double> A, G, GN; const std::vector<double>* Gp = nullptr;
  if (om) { A.assign(om, om + (size_t)d * d); for (double& e : A) e *= 2 * PI / N; la::expm_skew(d, A, G); Gp = &G; }
  std::vector<double> cp, cv; initial_state(N, d, r.mode_list(), r.coef.data(), om, cp, cv);
  double sc = 1.0; for (double v : cp) sc = std::max(sc, std::fabs(v)); for (double v : cv) sc = std::max(sc, std::fabs(v));
  coef_res = chore_residual(nb, cp, cv, 1e-16, Gp) / sc;
  std::vector<double> sp = cp, sv = cv;
  const double* z = r.state();
  if (z) { sp.assign(z, z + nd); sv.assign(z + nd, z + 2 * nd); }
  state_res = z ? chore_residual(nb, sp, sv, 1e-16, Gp) / sc : coef_res;
  if (period_ret) {
    std::vector<double> p = sp, v = sv; nb.integrate(p, v, 2 * PI, 1e-16);
    if (om) { A.assign(om, om + (size_t)d * d); for (double& e : A) e *= 2 * PI; la::expm_skew(d, A, GN); }
    double ret = 0;
    for (int k = 0; k < N; k++) for (int c = 0; c < d; c++) { double tp = sp[k * d + c], tv = sv[k * d + c];
      if (!GN.empty()) { tp = tv = 0; for (int b = 0; b < d; b++) { tp += GN[(size_t)c * d + b] * sp[k * d + b]; tv += GN[(size_t)c * d + b] * sv[k * d + b]; } }
      ret = std::max(ret, std::max(std::fabs(p[k * d + c] - tp), std::fabs(v[k * d + c] - tv))); }
    *period_ret = ret / sc; }
}

struct Orbit {
  std::vector<double> Z; int Ms = 0; double residual = INF, energy = 0, action = 0, minsep = INF, maxr = 0, Lnorm = 0, rms = 0;
  std::vector<double> Lsv, coef; std::vector<int> modes; int K = 0;
};
// integrate one T/N segment, assemble the loop from all N bodies (q(t+2πj/N) = q_j(t)), invariants, Fourier series of body 0
inline bool orbit_fit(int N, int d, double alpha, Orbit& O, Ctx& ctx, int Ms_req, int Kout_max, const std::vector<double>* Om = nullptr) {
  const int nd = N * d, Mseg = std::max(8, (Ms_req + N - 1) / N), Ms = Mseg * N; ctx.tables(Ms);
  std::vector<double> ts(Mseg); for (int j = 0; j < Mseg; j++) ts[j] = 2 * PI * j / Ms;
  std::vector<double> p(O.Z.begin(), O.Z.begin() + nd), v(O.Z.begin() + nd, O.Z.end()), seg;
  NBody<double> nbi(N, d, alpha, 22);
  nbi.integrate(p, v, 2 * PI / N, 1e-16, &ts, &seg);
  for (double x : seg) if (!std::isfinite(x)) return false;
  std::vector<double>& S = ctx.samp; S.assign((size_t)Ms * 2 * d, 0.0);
  std::vector<double> Rt, Ai;                                       // q(t) = exp(−Ωt) Q(t), q̇ = exp(−Ωt)(Q̇ − ΩQ)
  for (int i = 0; i < Mseg; i++) {
    if (Om) { Ai.assign(Om->begin(), Om->end()); for (double& z : Ai) z *= -ts[i]; la::expm_skew(d, Ai, Rt); }
    for (int j = 0; j < N; j++) { const double* qp = &seg[(size_t)i * 2 * nd + j * d]; const double* qv = &seg[(size_t)i * 2 * nd + nd + j * d];
      for (int c = 0; c < d; c++) { double pc = qp[c], vc = qv[c];
        if (Om) { pc = vc = 0;
          for (int b = 0; b < d; b++) { double wb = qv[b]; for (int e = 0; e < d; e++) wb -= (*Om)[(size_t)b * d + e] * qp[e];
            pc += Rt[(size_t)c * d + b] * qp[b]; vc += Rt[(size_t)c * d + b] * wb; } }
        S[((size_t)j * Mseg + i) * 2 * d + c] = pc; S[((size_t)j * Mseg + i) * 2 * d + d + c] = vc; } } }
  O.Ms = Ms; O.energy = nbody_energy(N, d, alpha, O.Z.data(), O.Z.data() + nd);
  double act = 0; O.minsep = INF; O.maxr = 0;
  for (int i = 0; i < Mseg; i++) { const double* s = &seg[(size_t)i * 2 * nd];
    for (int b = 0; b < N; b++) { double ke = 0, r2b = 0; for (int c = 0; c < d; c++) { ke += s[nd + b * d + c] * s[nd + b * d + c]; r2b += s[b * d + c] * s[b * d + c]; } act += 0.5 * ke; O.maxr = std::max(O.maxr, std::sqrt(r2b)); }
    for (int b = 0; b < N; b++) for (int l = b + 1; l < N; l++) { double r2 = 0; for (int c = 0; c < d; c++) { double df = s[b * d + c] - s[l * d + c]; r2 += df * df; }
      O.minsep = std::min(O.minsep, std::sqrt(r2)); act += alpha == 1.0 ? 1.0 / std::sqrt(r2) : std::pow(r2, -0.5 * alpha); } }
  O.action = act * (2 * PI / Ms);            // per body
  if (alpha == 1.0 && std::fabs(O.action * N + 6 * PI * O.energy) > 1e-6 * std::fabs(O.action * N)) return false;   // virial
  std::vector<double> L((size_t)d * d, 0.0);
  for (int i = 0; i < N; i++) for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) L[(size_t)a * d + b] += O.Z[i * d + a] * O.Z[nd + i * d + b] - O.Z[i * d + b] * O.Z[nd + i * d + a];
  O.Lnorm = 0; for (int a = 0; a < d; a++) for (int b = a + 1; b < d; b++) O.Lnorm += L[(size_t)a * d + b] * L[(size_t)a * d + b]; O.Lnorm = std::sqrt(O.Lnorm); O.Lsv = la::singular_values(d, L);
  int Kmax = std::min(Kout_max, Ms / 4); std::vector<double> cs((size_t)(Kmax + 1) * 2 * d, 0.0); double amax = 0;
  for (int m = 1; m <= Kmax; m++) for (int c = 0; c < d; c++) { double sc = 0, ss = 0;
    for (int j = 0; j < Ms; j++) { int idx = (int)(((long)m * j) % Ms); double q = S[(size_t)j * 2 * d + c]; sc += q * ctx.cosT[idx]; ss += q * ctx.sinT[idx]; }
    cs[(size_t)m * 2 * d + c] = 2 * sc / Ms; cs[(size_t)m * 2 * d + d + c] = 2 * ss / Ms; amax = std::max(amax, std::max(std::fabs(sc), std::fabs(ss)) * 2 / Ms); }
  int K = 4; for (int m = Kmax; m >= 1; m--) { double a = 0; for (int c = 0; c < 2 * d; c++) a = std::max(a, std::fabs(cs[(size_t)m * 2 * d + c])); if (a > 1e-13 * amax) { K = std::max(K, m); break; } }
  O.K = K; O.modes.clear(); O.coef.clear(); O.rms = 0;
  for (int m = 1; m <= K; m++) { if (m % N == 0) continue; O.modes.push_back(m); for (int c = 0; c < 2 * d; c++) { double v2 = cs[(size_t)m * 2 * d + c]; O.coef.push_back(v2); O.rms += 0.5 * v2 * v2; } }
  O.rms = std::sqrt(O.rms);
  return true;
}

struct TrialOut { bool ok = false; Record rec; std::string why; };

inline uint64_t trial_seed(uint64_t seed, uint64_t trial) { la::Rng r(seed * 0x9E3779B97F4A7C15ULL + trial * 0xD1B54A32D192ED03ULL + 0x632BE59BD9B4E019ULL); r.next(); return r.next(); }

// Newton polish inside symmetry S (or the full space)
inline bool polish(const Problem& P, const Symmetry& S, Ctx& ctx, std::vector<double>& x, int iters, double gtol, double* gnorm = nullptr) {
  const auto& BR = S.empty() ? std::pair<std::vector<double>, int>{} : ctx.basis(P, S);
  Reduced R(P, S.empty() ? nullptr : &BR.first, BR.second);
  std::vector<double> y(R.dim()); R.reduce(x.data(), y.data());
  auto fn = [&](const double* yy, double* g) { return R.f(yy, g, ctx.w); };
  auto fh = [&](const double* yy, std::vector<double>& H) { return R.hess(yy, H, ctx.w); };
  OptResult res = newton_lm(R.dim(), y, fn, fh, iters, gtol);
  R.expand(y.data(), x.data());
  if (gnorm) *gnorm = res.gnorm;
  return res.converged;
}

// catalogue solution embedded in the current space, kicked along soft Hessian directions or randomly
inline bool seeded_start(const Config& cfg, const Problem& P, Ctx& ctx, la::Rng& rng, std::vector<double>& x) {
  const auto& seeds = *cfg.seeds; if (seeds.empty()) return false;
  for (int attempt = 0; attempt < 12; attempt++) {
    const Record& r = seeds[rng.below((int)seeds.size())];
    if (r.h.N != P.N || r.h.d > P.d) continue;
    const int dr = r.h.d; x.assign(P.n, 0.0);
    for (size_t k = 0; k < r.modes.size(); k++) { int m = r.modes[k]; if (m > P.K) break; auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) continue;
      int mu = (int)(it - P.modes.begin()); for (int a = 0; a < dr; a++) { x[(2 * mu) * P.d + a] = r.coef[k * 2 * dr + a]; x[(2 * mu + 1) * P.d + a] = r.coef[k * 2 * dr + dr + a]; } }
    std::vector<double>& v = ctx.tmp; v.assign(P.n, 0.0);
    if (rng.uniform() < 0.5 && hessian(P, x.data(), ctx.H, ctx.w)) {
      la::sym_eig(P.n, ctx.H, ctx.eig); int nsoft = std::min(6, P.n);
      for (int k = 0; k < nsoft; k++) { double gk = rng.normal(); for (int i = 0; i < P.n; i++) v[i] += gk * ctx.H[(size_t)i * P.n + k]; }
    } else random_guess(P, rng, P.K, 1.0, v);
    double vn = 0; for (double e : v) vn += 0.5 * e * e; vn = std::sqrt(vn); if (vn == 0) continue;
    double eps = cfg.kick_min * std::pow(cfg.kick_max / cfg.kick_min, rng.uniform()) * r.h.rms / vn;
    for (int i = 0; i < P.n; i++) x[i] += eps * v[i];
    if (std::isfinite(action_grad(P, x.data(), nullptr, ctx.w))) return true;
  }
  return false;
}

// torus start: rotations in ⌊d/2⌋ planes with distinct modes coprime to N
inline int coprime_mode(int N, int K, la::Rng& rng) { for (int t = 0; t < 64; t++) { int m = 1 + rng.below(std::max(1, std::min(K, 7))); if (la::gcd(m, N) == 1) return m; } return 1; }
inline void torus_guess(const Problem& P, la::Rng& rng, std::vector<double>& x) {
  x.assign(P.n, 0.0); int planes = P.d / 2; std::vector<int> used;
  for (int pl = 0; pl < planes; pl++) {
    int m = coprime_mode(P.N, P.K, rng); for (int t = 0; t < 8 && std::find(used.begin(), used.end(), m) != used.end(); t++) m = coprime_mode(P.N, P.K, rng); used.push_back(m);
    auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) continue; int mu = (int)(it - P.modes.begin());
    double a = 0.3 + rng.uniform(), sgn = rng.uniform() < 0.5 ? 1 : -1;
    x[(2 * mu) * P.d + 2 * pl] = a; x[(2 * mu + 1) * P.d + 2 * pl + 1] = sgn * a;
  }
  if (P.d % 2) { int c = P.d - 1; for (int mu = 0; mu < P.nm && P.modes[mu] <= 5; mu++) { double a = 0.2 * std::pow((double)P.modes[mu], -1.0); x[(2 * mu) * P.d + c] = a * rng.normal(); x[(2 * mu + 1) * P.d + c] = a * rng.normal(); } }
}
// transverse mode frequencies of the rotating N-gon (ω_N = 1): ω_k² = Σ_{l≠0} (1 − cos(2πkl/N)) / d_l³, d_l = 2R sin(πl/N)
inline double ngon_vertical_freq(int N, int k) {
  double R3 = 0; for (int l = 1; l < N; l++) R3 += 1.0 / (4.0 * std::sin(PI * l / N)); double R = std::cbrt(R3), w2 = 0;
  for (int l = 1; l < N; l++) { double dl = 2 * R * std::sin(PI * l / N); w2 += (1 - std::cos(2 * PI * k * l / N)) / (dl * dl * dl); }
  return std::sqrt(w2);
}
// N-gon sums C_k = Σ_{p≠0} g_p cos(2πkp/N), g_p = α/d_p^{α+2}; C_0 − C_1 = 1 *is* the radius equation.
inline void ngon_C(int N, double alpha, double* C) {
  double S = 0; for (int p = 1; p < N; p++) S += alpha * std::sin(PI * p / N) / std::pow(2 * std::sin(PI * p / N), alpha + 1);
  const double R = std::pow(S, 1.0 / (alpha + 2));
  for (int k = 0; k < N; k++) { double c = 0;
    for (int p = 1; p < N; p++) { const double dp = 2 * R * std::sin(PI * p / N); c += alpha * std::cos(2 * PI * k * p / N) / std::pow(dp, alpha + 2); }
    C[k] = c; }
}
// real roots of t³ + p t + q, ascending
inline int cubic_depressed(double p, double q, double* r) {
  const double D = q * q / 4 + p * p * p / 27;
  if (D >= 0) { const double s = std::sqrt(D); r[0] = std::cbrt(-q / 2 + s) + std::cbrt(-q / 2 - s); return 1; }
  const double m = 2 * std::sqrt(-p / 3), th = std::acos(std::clamp(3 * q / (p * m), -1.0, 1.0)) / 3;
  for (int i = 0; i < 3; i++) r[i] = m * std::cos(th - 2 * PI * i / 3);
  std::sort(r, r + 3); return 3;
}
// real roots of ν⁴ + a ν² + b ν + c, bracketed between the critical points (the cubic 4ν³ + 2aν + b)
inline int quartic_real(double a, double b, double c, double* r) {
  auto f = [&](double v) { return ((v * v + a) * v + b) * v + c; };
  double cr[3]; const int nc = cubic_depressed(a / 2, b / 4, cr);
  const double M = 1 + std::fabs(a) + std::fabs(b) + std::fabs(c);
  double xs[5]; int nx = 0; xs[nx++] = -M;
  for (int i = 0; i < nc; i++) if (cr[i] > -M && cr[i] < M) xs[nx++] = cr[i];
  xs[nx++] = M; int n = 0;
  for (int i = 0; i + 1 < nx; i++) { double lo = xs[i], hi = xs[i + 1]; const double fl = f(lo);
    if (fl * f(hi) > 0) continue;
    for (int it = 0; it < 100; it++) { const double mid = 0.5 * (lo + hi); if ((f(mid) < 0) == (fl < 0)) lo = mid; else hi = mid; }
    r[n++] = 0.5 * (lo + hi); }
  return n;
}
// In-plane block of the rotating N-gon (Maxwell's ring): bodies in their own radial/tangential frames,
// η_j = v e^{2πikj/N}. Real roots give the resonances m2/m1 = 1 ± ν, pattern k ≡ m2·m1⁻¹ − 1 (mod N).
inline int ngon_inplane_freq(int N, int k, double alpha, double* nu) {
  double C[128]; if (N > 128) return 0;
  ngon_C(N, alpha, C);
  auto Cj = [&](int j) { return C[((j % N) + N) % N]; };
  const double ck = Cj(k), cs = Cj(k + 1) + Cj(k - 1);
  const double Pk = 0.5 * ((alpha + 2) * ck - 0.5 * alpha * cs - 2 * C[0] + alpha + 4);
  const double Qk = 0.5 * (-(alpha + 2) * ck - 0.5 * alpha * cs + (2 * alpha + 2) * C[0] - alpha);
  const double bk = 0.25 * alpha * (Cj(k - 1) - Cj(k + 1));
  return quartic_real(Pk + Qk - 4, -4 * bk, Pk * Qk - bk * bk, nu);
}
inline void set_mode(const Problem& P, std::vector<double>& x, int m, int c, double cc, double ss) {
  auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) return;
  const size_t mu = (size_t)(it - P.modes.begin()); x[2 * mu * P.d + c] = cc; x[(2 * mu + 1) * P.d + c] = ss;
}
// vertical start: circle (mode m1) + one transverse mode m2, biased to resonances m2/m1 ≈ ω_k/ω_N, k ≡ m2·m1⁻¹ (mod N)
inline void vertical_guess(const Problem& P, la::Rng& rng, std::vector<double>& x) {
  x.assign(P.n, 0.0); const int N = P.N, d = P.d; int m1 = 1, m2 = 2;
  for (int t = 0; t < 200; t++) {
    int a = 1 + rng.below(std::min(P.K, 12)), b = 1 + rng.below(std::min(P.K, 12));
    if (la::gcd(a, N) != 1 || b == a || b % N == 0 || b > P.K) continue;
    int ainv = 1; for (int u = 1; u < N; u++) if ((a * u) % N == 1) ainv = u;
    if (rng.uniform() < 0.7) { double target = ngon_vertical_freq(N, (b * ainv) % N); if (std::fabs((double)b / a - target) > 0.2 * target) continue; }
    m1 = a; m2 = b; break;
  }
  double sgn = rng.uniform() < 0.5 ? 1 : -1, amp = 0.08 * std::pow(15.0, rng.uniform());
  set_mode(P, x, m1, 0, 1.0, 0.0); set_mode(P, x, m1, 1, 0.0, sgn);
  double ph = 2 * PI * rng.uniform(); set_mode(P, x, m2, 2, amp * std::cos(ph), amp * std::sin(ph));
  for (int c = 3; c < d; c++) if (rng.uniform() < 0.5) { int m3 = 1 + rng.below(std::min(P.K, 8)); if (m3 % N) set_mode(P, x, m3, c, 0.3 * amp * rng.normal(), 0.3 * amp * rng.normal()); }
}

// transverse resonances m ≈ m1·ω_k, m ≡ k·m1 (mod N), distinct patterns k, circularly polarised in pairs
// on axes c0…d−1. Each pair adds 2 to deff, so deff ≤ 2 + 2(⌊N/2⌋−1).
inline void add_vertical_modes(const Problem& P, la::Rng& rng, std::vector<double>& x, int m1, int c0) {
  const int N = P.N, d = P.d, nk = N / 2 - 1; if (nk <= 0) return;
  std::vector<int> ks(nk); for (int i = 0; i < nk; i++) ks[i] = i + 2;
  for (int i = nk - 1; i > 0; i--) std::swap(ks[i], ks[rng.below(i + 1)]);
  for (int c = c0, i = 0; c < d && i < 4 * nk; i++) {
    const int k = ks[i % nk]; const double tgt = m1 * ngon_vertical_freq(N, k);
    const int base = (k * m1) % N, m = base + N * (int)std::lround((tgt - base) / N);
    if (m < 1 || m % N == 0 || m > P.K || std::fabs(m - tgt) > 0.25 * tgt) continue;
    const double amp = 0.03 * std::pow(20.0, rng.uniform()), ph = 2 * PI * rng.uniform();
    set_mode(P, x, m, c, amp * std::cos(ph), amp * std::sin(ph));
    if (c + 1 < d) { const double s2 = rng.uniform() < 0.5 ? 1 : -1;
      set_mode(P, x, m, c + 1, -s2 * amp * std::sin(ph), s2 * amp * std::cos(ph)); c += 2; } else c++;
  }
}
// m1-fold circle in plane (0,1) plus the transverse resonances
inline void hyper_guess(const Problem& P, la::Rng& rng, std::vector<double>& x) {
  x.assign(P.n, 0.0);
  int m1 = 1; for (int t = 0; t < 64; t++) { const int a = 1 + rng.below(std::max(1, P.K / 2)); if (la::gcd(a, P.N) == 1) { m1 = a; break; } }
  set_mode(P, x, m1, 0, 1.0, 0.0); set_mode(P, x, m1, 1, 0.0, rng.uniform() < 0.5 ? 1 : -1);
  add_vertical_modes(P, rng, x, m1, 2);
}

// in-plane start: m1-fold circle plus a resonance m2 ≈ m1(1 ± ν_k), m2 ≡ (k+1)·m1 (mod N). The in-plane
// Hessian couples m2 only to |2m1 − m2|; transverse modes go on top, else the loop stays planar.
inline bool inplane_guess(const Problem& P, la::Rng& rng, std::vector<double>& x) {
  x.assign(P.n, 0.0); const int N = P.N, K = P.K; if (N < 4) return false;
  int m1 = 1; for (int t = 0; t < 64; t++) { const int a = 1 + rng.below(std::max(1, std::min(K / 2, 12))); if (la::gcd(a, N) == 1) { m1 = a; break; } }
  int cand[32], nc = 0; double nu[4];
  for (int k = 2; k <= N / 2 && nc < 32; k++) {
    const int nr = ngon_inplane_freq(N, k, P.alpha, nu);
    for (int i = 0; i < nr && nc < 32; i++) for (int s = 0; s < 2; s++) {
      const double ratio = 1 + (s ? -nu[i] : nu[i]); if (ratio < 0.1) continue;
      const double tgt = ratio * m1; const int base = ((k + 1) * m1) % N;
      const int m = base + N * (int)std::lround((tgt - base) / N);
      if (m < 1 || m > K || m % N == 0 || m == m1 || std::fabs(m - tgt) > 0.25 * tgt) continue;
      cand[nc++] = m;
    }
  }
  if (!nc) return false;
  const int m2 = cand[rng.below(nc)], m3 = std::abs(2 * m1 - m2);
  const double sgn = rng.uniform() < 0.5 ? 1 : -1, amp = 0.05 * std::pow(10.0, rng.uniform()), ph = 2 * PI * rng.uniform();
  set_mode(P, x, m1, 0, 1.0, 0.0); set_mode(P, x, m1, 1, 0.0, sgn);
  set_mode(P, x, m2, 0, amp * std::cos(ph), amp * std::sin(ph));
  set_mode(P, x, m2, 1, -sgn * amp * std::sin(ph), sgn * amp * std::cos(ph));
  if (m3 % N) { const double a3 = 0.4 * amp;
    set_mode(P, x, m3, 0, a3 * std::sin(ph), -a3 * std::cos(ph));
    set_mode(P, x, m3, 1, sgn * a3 * std::cos(ph), sgn * a3 * std::sin(ph)); }
  if (P.d > 2 && rng.uniform() < 0.75) add_vertical_modes(P, rng, x, m1, 2);
  return true;
}

// d = 7 start adapted to q(t + 2πp/7) = σ q(t): mode m sits entirely in the σ-eigenplane k = mp mod 7, so the
// loop spans three planes (deff ≥ 6, and only then is the G2 twist non-zero) only if its modes hit three
// distinct k, and reaches deff = 7 only by also carrying a mode ≡ 0 (mod 7) on the fixed axis. Modes are drawn
// from a resonant triple m1 + m2 = m3, which is what makes the twist integral non-zero.
inline bool fano_guess(const Problem& P, la::Rng& rng, int p, std::vector<double>& x) {
  const int K = P.K, N = P.N; x.assign(P.n, 0.0);
  std::vector<int> tri;
  for (int t = 0; t < 400 && tri.empty(); t++) {
    int m1 = 1 + rng.below(std::min(K, 9)), m2 = 1 + rng.below(std::min(K, 9)), m3 = m1 + m2;
    if (m3 > K || m1 % N == 0 || m2 % N == 0 || m3 % N == 0) continue;
    int a = fano_plane(m1, p), b = fano_plane(m2, p), c = fano_plane(m3, p);
    if (a && b && c && a != b && b != c && a != c) tri = {m1, m2, m3};
  }
  if (tri.empty()) return false;
  for (int ax = 7; ax <= K; ax += 7) if (ax % N) { tri.push_back(ax); break; }      // the fixed axis
  double u[7], v[7];
  for (int m : tri) {
    auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) continue;
    const int mu = (int)(it - P.modes.begin()), k = (((m * p) % 7) + 7) % 7;
    fano_basis(k > 3 ? 7 - k : k, u, v);
    const double sg = k > 3 ? -1.0 : 1.0, A = (0.3 + rng.uniform()) * (k ? 1.0 : 0.5), ps = 2 * PI * rng.uniform();
    for (int a = 0; a < 7; a++) { x[(2 * mu) * 7 + a] = A * (std::cos(ps) * u[a] + sg * std::sin(ps) * v[a]);
      x[(2 * mu + 1) * 7 + a] = A * (-std::sin(ps) * u[a] + sg * std::cos(ps) * v[a]); }
  }
  return true;
}

// orthonormal gauge directions at x: the time shift, plus the rotations that survive the rotating frame.
// Those are the centraliser {L in so(d) : [L, Omega] = 0}, taken as the kernel of ad_Omega in the E_ab basis:
// when the frame's planes are not coordinate planes (the g2 torus) no single E_ab commutes, so testing them
// one at a time finds nothing and leaves the true null directions in the spectrum.
inline int gauge_basis(const Problem& P, const double* x, std::vector<double>& G) {
  const int d = P.d, D = d * (d - 1) / 2;
  std::vector<int> ia(D), ib(D);
  for (int a = 0, k = 0; a < d; a++) for (int b = a + 1; b < d; b++, k++) { ia[k] = a; ib[k] = b; }
  std::vector<double> L((size_t)D * D, 0.0); int nl = 0;                  // centraliser basis in E_ab coordinates
  if (P.Om.empty()) { for (int k = 0; k < D; k++) L[(size_t)k * D + k] = 1.0; nl = D; }
  else {
    std::vector<double> M((size_t)D * D, 0.0), S((size_t)D * D, 0.0), ev; const double* O = P.Om.data();
    for (int r = 0; r < D; r++) for (int c = 0; c < D; c++) { const int a = ia[r], b = ib[r], e = ia[c], f = ib[c];
      double v = 0;                                                       // [Omega, E_ef]_{ab}, E_ef = e_e f^T - e_f e^T
      if (b == f) v += O[(size_t)a * d + e];
      if (b == e) v -= O[(size_t)a * d + f];
      if (a == e) v -= O[(size_t)f * d + b];
      if (a == f) v += O[(size_t)e * d + b];
      M[(size_t)r * D + c] = v; }
    for (int i = 0; i < D; i++) for (int j = i; j < D; j++) { double s = 0;
      for (int k = 0; k < D; k++) s += M[(size_t)k * D + i] * M[(size_t)k * D + j];
      S[(size_t)i * D + j] = S[(size_t)j * D + i] = s; }
    la::sym_eig(D, S, ev);
    double emax = 0; for (double v : ev) emax = std::max(emax, std::fabs(v));
    for (int k = 0; k < D; k++) if (std::fabs(ev[k]) <= 1e-12 * emax)
      { for (int i = 0; i < D; i++) L[(size_t)nl * D + i] = S[(size_t)i * D + k]; nl++; }
  }
  std::vector<std::vector<double>> gs(1); deriv_coeffs(P, x, gs[0]);
  for (int k = 0; k < nl; k++) { std::vector<double> v(P.n, 0.0);
    for (int r = 0; r < D; r++) { const double c = L[(size_t)k * D + r]; if (c == 0.0) continue;
      const int a = ia[r], b = ib[r];
      for (int i = 0; i < P.nb; i++) { v[i * d + a] -= c * x[i * d + b]; v[i * d + b] += c * x[i * d + a]; } }
    gs.push_back(v); }
  double xn = 0; for (int i = 0; i < P.n; i++) xn += x[i] * x[i]; xn = std::sqrt(xn);
  std::vector<std::vector<double>> on;
  for (auto& g : gs) { for (auto& u : on) { double dp = 0; for (int i = 0; i < P.n; i++) dp += u[i] * g[i];
      for (int i = 0; i < P.n; i++) g[i] -= dp * u[i]; }
    double nn = 0; for (double v : g) nn += v * v; nn = std::sqrt(nn);
    if (nn < 1e-8 * xn) continue; for (double& v : g) v /= nn; on.push_back(g); }
  int ng = (int)on.size(); G.assign((size_t)P.n * ng, 0.0);
  for (int k = 0; k < ng; k++) for (int i = 0; i < P.n; i++) G[(size_t)i * ng + k] = on[k][i];
  return ng;
}
// H + σ G Gᵀ lifts the gauge directions out of the spectrum, so the index no longer depends on a magnitude
// threshold; the tolerance for the rest comes from ‖H G‖ by Weyl's inequality
inline Inertia inertia_gauge(const Problem& P, const double* x, std::vector<double> H, int* ng_out = nullptr) {
  std::vector<double> G; int ng = gauge_basis(P, x, G); if (ng_out) *ng_out = ng;
  double sg = 0, tol = 0; for (int i = 0; i < P.n; i++) sg = std::max(sg, std::fabs(H[(size_t)i * P.n + i]));
  for (int k = 0; k < ng; k++) for (int i = 0; i < P.n; i++) { double s = 0;
    for (int j = 0; j < P.n; j++) s += H[(size_t)i * P.n + j] * G[(size_t)j * ng + k]; tol = std::max(tol, std::fabs(s)); }
  for (int i = 0; i < P.n; i++) for (int j = 0; j < P.n; j++) { double s = 0;
    for (int k = 0; k < ng; k++) s += G[(size_t)i * ng + k] * G[(size_t)j * ng + k]; H[(size_t)i * P.n + j] += sg * s; }
  Inertia I; la::sym_eig(P.n, H, I.eig); tol = std::max(tol, 16 * 2.220446049250313e-16 * sg);
  for (double v : I.eig) { if (std::fabs(v) <= tol) I.zero++; else if (v < 0) I.neg++; else I.pos++; }
  I.zero += ng; return I;
}

// Morse index, nullity and the calibration twist — computed only for records that are new or replaced.
// `index` off stops before the index, which is the expensive half (177 s a record at d = 11) and the
// only half a caller deciding whether it wants the record at all does not need.
inline void record_extras(const Config& cfg, Record& rec, Ctx& ctx, bool index = true) {
  int ck = 0, cid = 0; calib_pick(rec.h.d, ck, cid);
  if (ck && (int)rec.extra.size() >= Record::NEX) {
    const Wedge& W = wedge_basis(rec.h.d, ck); const std::vector<double> psi = calib_psi(W, cid);
    std::vector<double> A(W.n);
    jet_form(rec.mode_list(), rec.coef.data(), rec.h.d, W, A.data());
    double na = 0; for (double v : A) na += v * v;
    const double sc = jet_scale(rec.mode_list(), rec.coef.data(), rec.h.d, ck);
    rec.extra[0] = calib_max(W, A.data(), psi);
    rec.extra[1] = sc > 0 ? rec.extra[0] / sc : 0.0;              // scale-free twist
    rec.extra[2] = ck; rec.extra[3] = sc > 0 ? std::sqrt(na) / sc : 0.0;
  }
  if (rec.h.alpha != 1.0) return;                        // the index below is taken at the Newtonian problem
  const Problem& Pc = ctx.problem(cfg, std::max(cfg.K_index, cfg.K));
  Problem Pl; const Problem* Pp = &Pc;                     // the record's Ω is in its own canonical frame
  if (const double* om = rec.omega()) { Pl = Pc; Pl.set_omega(std::vector<double>(om, om + (size_t)rec.h.d * rec.h.d)); Pp = &Pl; }
  const Problem& P = *Pp; ctx.w.resize(P);
  std::vector<double> x(P.n, 0.0);
  for (size_t k = 0; k < rec.modes.size(); k++) { int m = rec.modes[k]; if (m > P.K) break; auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) continue;
    int mu = (int)(it - P.modes.begin()); for (int a = 0; a < P.d; a++) { x[(2 * mu) * P.d + a] = rec.coef[k * 2 * P.d + a]; x[(2 * mu + 1) * P.d + a] = rec.coef[k * 2 * P.d + P.d + a]; } }
  if ((int)rec.extra.size() >= Record::NEX) rec.extra[4] = rigid_defect(P, x.data(), ctx.w);
  if (const double* om = rec.omega()) { std::vector<double> isv; rec.h.deff = inertial_deff(rec.h.N, rec.mode_list(), rec.h.d, rec.coef.data(), om, isv); }
  double gn = 0; ctx.g.resize(P.n); double A0 = action_grad(P, x.data(), ctx.g.data(), ctx.w);
  if (!std::isfinite(A0)) return; for (double v : ctx.g) gn += v * v;
  rec.h.grad_norm = std::sqrt(gn);                       // residual of the truncation, before the polish
  if (!index) return;
  Symmetry none; polish(P, none, ctx, x, cfg.newton_iters, cfg.gtol);   // the index is only defined at a critical point
  double A1 = action_grad(P, x.data(), nullptr, ctx.w);
  if (std::fabs(A1 - A0) < 1e-3 * std::fabs(A0) && hessian(P, x.data(), ctx.H, ctx.w))
    { Inertia I = inertia_gauge(P, x.data(), ctx.H); rec.h.morse = I.neg; rec.h.nullity = I.zero; }
}

// residual of the stored coefficients; h.ret_err is the state's, and the gap is the K-mode truncation
inline double coef_residual(int N, int d, const std::vector<int>& modes, const double* coef, const double* om, NBody<double>& nb) {
  std::vector<double> pos, vel; initial_state(N, d, modes, coef, om, pos, vel);
  double sc = 1.0; for (double v : pos) sc = std::max(sc, std::fabs(v)); for (double v : vel) sc = std::max(sc, std::fabs(v));
  std::vector<double> G; const std::vector<double>* Gp = nullptr;
  if (om) { std::vector<double> A(om, om + (size_t)d * d); for (double& e : A) e *= 2 * PI / N; la::expm_skew(d, A, G); Gp = &G; }
  return chore_residual(nb, pos, vel, 1e-16, Gp) / sc;
}
inline double coef_residual(const Record& rec, NBody<double>& nb) {
  return coef_residual(rec.h.N, rec.h.d, rec.mode_list(), rec.coef.data(), rec.omega(), nb);
}

#ifdef HAVE_MPFR
// exp(sΩ) at the working precision; Ω is antisymmetrised first, or exp is not orthogonal and G S Z is not
// the symmetry the orbit satisfies
inline void omega_exp(int d, const double* Om, const mpreal& s, std::vector<mpreal>& R) {
  std::vector<mpreal> A((size_t)d * d);
  for (int a = 0; a < d; a++) for (int b = 0; b < d; b++)
    A[(size_t)a * d + b] = mpreal(0.5 * (Om[(size_t)a * d + b] - Om[(size_t)b * d + a])) * s;
  la::expm_skew(d, A, R);
}
#endif
inline constexpr int MP_DIGITS = 30;                  // the fallback precision; main() sets it once, globally

// Certify a state in place: the double shooting Newton, then MPFR where the monodromy stalls it above
// mp_gate. Rounding the MPFR answer back is amplified by that same monodromy, but lands where double cannot.
// Returns the relative residual, and in mp_out the MPFR one (-1 if it was not needed).
inline double certify_state(int N, int d, double alpha, std::vector<double>& Z, const double* Om,
                            const std::vector<double>* G, NBody<double>& nb, ShootWork<double>& sw,
                            double tol, int iters, double mp_gate, int threads,
                            double* mp_out = nullptr, bool* mp_used = nullptr) {
  const int nd = N * d;
  double sc = 1.0; for (double v : Z) sc = std::max(sc, std::fabs(v));
  double res = shoot_newton(nb, Z, 1e-16, iters, tol * sc, -18, -40, sw, false, G, 1, 0.9, 24) / sc;
  if (mp_out) *mp_out = -1;
  if (mp_used) *mp_used = false;
#ifdef HAVE_MPFR
  if (mp_gate > 0 && res > mp_gate) {
    const int order = (int)(1.15 * MP_DIGITS) + 6; const double itol = std::pow(10.0, -(MP_DIGITS + 4));
    const mpfr_prec_t bits = mpreal::default_prec();
    std::vector<mpreal> Zm(2 * (size_t)nd); for (int c = 0; c < 2 * nd; c++) Zm[c] = Z[c];
    std::vector<mpreal> Gm; const std::vector<mpreal>* Gmp = nullptr;
    if (Om) { omega_exp(d, Om, mpreal::pi() * 2 / N, Gm); Gmp = &Gm; }
    NBody<mpreal> nbm(N, d, alpha, order); ShootWork<mpreal> swm;
    double mp = shoot_newton(nbm, Zm, itol, 8, std::pow(10.0, -MP_DIGITS), -(long)(bits / 3),
                             -(long)std::min<mpfr_prec_t>(bits / 2, 96), swm, false, Gmp, threads, 0.9, 24) / sc;
    if (mp_out) *mp_out = mp;
    std::vector<double> zp(nd), zv(nd);
    for (int c = 0; c < nd; c++) { zp[c] = to_double(Zm[c]); zv[c] = to_double(Zm[nd + c]); }
    double s2 = 1.0; for (double v : zp) s2 = std::max(s2, std::fabs(v)); for (double v : zv) s2 = std::max(s2, std::fabs(v));
    const double rr = chore_residual(nb, zp, zv, 1e-16, G) / s2;
    if (rr < res) { res = rr; if (mp_used) *mp_used = true; for (int c = 0; c < nd; c++) { Z[c] = zp[c]; Z[nd + c] = zv[c]; } }
  }
#else
  (void)alpha; (void)Om; (void)mp_gate; (void)threads; (void)nd; (void)mp_used;
#endif
  return res;
}

// The loop a certified state defines, with enough modes and samples to represent it: more modes when the
// spectrum is still above threshold at the cap, and more samples because Ms = 2048 is four per mode at
// K ≈ 500, where the tail aliases back down. Only large-K loops pay for either.
inline bool fit_loop(int N, int d, double alpha, Orbit& O, Ctx& ctx, int Ms, int Kout_max, int Kgrow,
                     const std::vector<double>* Om, NBody<double>& nb, double* coef_err = nullptr) {
  int kmax = Kout_max, ms = Ms;
  bool fit = orbit_fit(N, d, alpha, O, ctx, ms, kmax, Om);
  auto cres = [&](const Orbit& q) { return coef_residual(N, d, q.modes, q.coef.data(), Om ? Om->data() : nullptr, nb); };
  double best = fit ? cres(O) : INF;
  for (int pass = 0; fit && pass < 4; pass++) {
    const int want_k = (O.K >= kmax - 1 && kmax < Kgrow) ? std::min(Kgrow, kmax * 4) : kmax;
    const int want_ms = std::max(ms, 16 * O.K);
    if (want_k == kmax && want_ms <= ms) break;
    Orbit keep = O; const int k0 = kmax, m0 = ms;
    kmax = want_k; ms = want_ms;
    fit = orbit_fit(N, d, alpha, O, ctx, ms, kmax, Om);
    const double got = fit ? cres(O) : INF;
    // a spectrum that will not decay, because the state is not quite a solution, would grow to the cap
    if (!(got < 0.5 * best)) { O = keep; kmax = k0; ms = m0; fit = true; break; }
    best = got;
  }
  if (coef_err) *coef_err = fit ? best : INF;
  return fit;
}


// ODE validation, cover unwinding, shooting certification, Fourier re-extraction, canonical frame → record
inline bool certify(const Config& cfg, const Problem*& P, std::vector<double>& x, const Symmetry& S, Ctx& ctx, Record& rec, std::string& why) {
  // at Ω = 0 deff can only be lost downstream, so gate before the ODE work; 1e-12 keeps the gate permissive
  if (cfg.min_deff > 1) { std::vector<double> xc(x), sv; if (canonical_frame(P->nb, P->d, xc, sv, 1e-12) < cfg.min_deff) { why = "effective dimension below filter"; return false; } }
  // rigid loops are the N-gon's high-d analogue and dominate the deff filter in a rotating frame
  if (rigid_defect(*P, x.data(), ctx.w) < cfg.min_rigid) { why = "relative equilibrium"; return false; }
  // one mode doubling if badly under-resolved
  NBody<double>& nb = ctx.integrator(cfg);
  double ret = return_error(*P, x.data(), nb);
  std::vector<double> x2;
  if (ret > cfg.ret_double && P->K < cfg.Kmax) {
    const Problem& P2 = ctx.problem(cfg, std::min(2 * P->K, cfg.Kmax));
    x2.assign(P2.n, 0.0); P2.transfer(*P, x.data(), x2.data());
    if (!polish(P2, S, ctx, x2, cfg.newton_iters, cfg.gtol)) { why = "refinement lost"; return false; }
    x.swap(x2); P = &P2; ret = return_error(*P, x.data(), nb);
  }
  if (!(ret <= cfg.ret_reject)) { why = "spurious (return error)"; return false; }
  int cover = cover_multiplicity(*P, x.data());
  if (cover > 1) {
    Problem Pu; unwind_cover(*P, x.data(), cover, Pu, x2); Pu.minsep = cfg.minsep;
    const Problem& P3 = ctx.problem(cfg, std::max(Pu.K, cfg.K)); std::vector<double> x3(P3.n, 0.0); P3.transfer(Pu, x2.data(), x3.data());
    Symmetry none; if (!polish(P3, none, ctx, x3, cfg.newton_iters, cfg.gtol)) { why = "unwound polish lost"; return false; }
    x.swap(x3); P = &P3; ret = return_error(*P, x.data(), nb);
    if (!(ret <= cfg.ret_reject)) { why = "spurious (return error)"; return false; }
  }
  Orbit O; initial_state(*P, x.data(), ctx.pos, ctx.vel); O.Z = ctx.pos; O.Z.insert(O.Z.end(), ctx.vel.begin(), ctx.vel.end());
  double res = certify_state(P->N, P->d, 1.0, O.Z, P->Om.empty() ? nullptr : P->Om.data(), P->gshift(),
                             nb, ctx.sw, cfg.shoot_tol, 10, cfg.mpfr_gate, 1);
  if (!(res <= cfg.shoot_tol * 100)) { why = "shooting did not certify"; return false; }
  if (!fit_loop(P->N, P->d, 1.0, O, ctx, cfg.Ms, cfg.Kout_max, cfg.Kout_grow, P->Om.empty() ? nullptr : &P->Om, nb))
    { why = "orbit sampling failed"; return false; }
  O.residual = res;
  std::vector<double> sv, Rc; int deff = canonical_frame((int)O.modes.size() * 2, P->d, O.coef, sv, 1e-8, &Rc);
  const int d = P->d; std::vector<double> Omc;
  if (!P->Om.empty()) { Omc.assign((size_t)d * d, 0.0);                                    // Ω → R Ω Rᵀ
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double t = 0;
      for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) t += Rc[(size_t)a * d + i] * P->Om[(size_t)i * d + j] * Rc[(size_t)b * d + j];
      Omc[(size_t)a * d + b] = t; }
    std::vector<double> isv; deff = inertial_deff(P->N, O.modes, d, O.coef.data(), Omc.data(), isv);   // the dimension actually occupied
  }
  if (deff < cfg.min_deff) { why = "effective dimension below filter"; return false; }
  // the gate above ran before shooting; a trial can still land on a relative equilibrium in between
  if (rigid_defect_coef(P->N, d, O.modes, O.coef.data()) < cfg.min_rigid) { why = "relative equilibrium"; return false; }
  rec.h.N = P->N; rec.h.d = P->d; rec.h.K = O.K; rec.h.M = O.Ms; rec.h.alpha = 1.0; rec.modes.assign(O.modes.begin(), O.modes.end()); rec.coef = O.coef; rec.h.nm = (int32_t)O.modes.size();
  rec.h.deff = deff; rec.h.cover = cover; rec.h.action = O.action; rec.h.energy = O.energy; rec.h.energy_std = 0; rec.h.rms = O.rms; rec.h.maxr = O.maxr;
  rec.h.minsep = O.minsep; rec.h.Lnorm = O.Lnorm; rec.Lsv = O.Lsv; rec.pca = sv; rec.h.morse = -1; rec.h.nullity = -1; rec.h.grad_norm = -1; rec.h.ret_err = res;
  const size_t nom = (size_t)d * d, nst = 2 * (size_t)P->N * d; const int nd = P->N * d;
  rec.extra.assign(Record::NEX + nom + nst, 0.0);
  rec.extra[5] = 1;                                          // layout 1: Ω block always present, then the state
  if (!Omc.empty()) std::copy(Omc.begin(), Omc.end(), rec.extra.begin() + Record::NEX);
  // the certified state, rotated into the coefficients' canonical axes; h.ret_err describes this, not them
  double* Zc = &rec.extra[Record::NEX + nom];
  for (int half = 0; half < 2; half++) for (int k = 0; k < P->N; k++) for (int a = 0; a < d; a++) {
    double t = 0; for (int b = 0; b < d; b++) t += Rc[(size_t)a * d + b] * O.Z[half * nd + k * d + b];
    Zc[half * nd + k * d + a] = t; }
  rec.extra[6] = coef_residual(rec, nb);
  return true;
}

inline TrialOut run_trial(const Config& cfg, uint64_t trial, Ctx& ctx) {
  TrialOut out; auto t0 = std::chrono::steady_clock::now();
  la::Rng rng(trial_seed(cfg.seed, trial));
  std::string symtext = cfg.sym == "random" ? Symmetry::random_text(cfg.d, rng) : (cfg.sym == "none" ? "" : cfg.sym);
  Symmetry S = Symmetry::parse(symtext, cfg.d);
  const Problem* P = &ctx.problem(cfg, cfg.K);
  // q(εt+θ) = R q(t) leaves ½∫|q̇ + Ωq|² invariant only when RᵀΩR = εΩ
  if (!P->Om.empty()) for (const auto& g : S.gens) { const int d = cfg.d; double e = 0;
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double t = 0;
      for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) t += g.R[(size_t)i * d + a] * P->Om[(size_t)i * d + j] * g.R[(size_t)j * d + b];
      e = std::max(e, std::fabs(t - g.eps * P->Om[(size_t)a * d + b])); }
    if (e > 1e-9) { out.why = "symmetry incompatible with the rotating frame"; return out; } }
  Problem Pa; if (cfg.alpha_start != 1.0) { Pa = *P; Pa.alpha = cfg.alpha_start; P = &Pa; }
  ctx.w.resize(ctx.problem(cfg, cfg.Kmax));
  const auto& BR = S.empty() ? std::pair<std::vector<double>, int>{} : ctx.basis(*P, S);
  if (!S.empty() && BR.second < 0) { out.why = "symmetry group too large"; return out; }
  if (!S.empty() && BR.second < cfg.d) { out.why = "symmetry subspace too small"; return out; }
  Reduced R(*P, S.empty() ? nullptr : &BR.first, BR.second);
  const int r = R.dim();
  std::vector<double> x(P->n), y(r);
  bool started = false;
  std::vector<std::string> fam; { std::stringstream ss(cfg.starts); std::string f; while (std::getline(ss, f, ',')) if (!f.empty()) fam.push_back(f); }
  std::string family = fam.empty() ? "random" : fam[rng.below((int)fam.size())];
  if (family == "kick" && cfg.seeds && S.empty() && cfg.alpha_start == 1.0 && seeded_start(cfg, *P, ctx, rng, x)) { R.reduce(x.data(), y.data()); started = true; }
  for (int attempt = 0; attempt < 8 && !started; attempt++) {
    int K0 = cfg.K0min + rng.below(cfg.K0max - cfg.K0min + 1); double gamma = 0.7 + 1.3 * rng.uniform();
    if (family == "hyper" && P->d >= 3) hyper_guess(*P, rng, x);
    else if (family == "inplane") { if (!inplane_guess(*P, rng, x)) random_guess(*P, rng, K0, gamma, x); }
    else if (family == "fano" && P->d == 7) { if (!fano_guess(*P, rng, 1 + rng.below(6), x)) random_guess(*P, rng, K0, gamma, x); }
    else if (family == "torus" && P->d >= 4) torus_guess(*P, rng, x); else if (family == "vertical" && P->d >= 3) vertical_guess(*P, rng, x); else random_guess(*P, rng, K0, gamma, x);
    R.reduce(x.data(), y.data());
    double kin, pot; double A = R.f(y.data(), nullptr, ctx.w, &kin, &pot);
    if (!std::isfinite(A)) continue;
    double lam = kin > 0 ? optimal_scale(kin, pot, P->alpha) : 1.0; for (double& v : y) v *= lam;
    if (std::isfinite(R.f(y.data(), nullptr, ctx.w))) started = true;
  }
  if (!started) { out.why = "no collision-free start"; return out; }
  auto fn = [&](const double* yy, double* g) { return R.f(yy, g, ctx.w); };
  auto fh = [&](const double* yy, std::vector<double>& H) { return R.hess(yy, H, ctx.w); };
  // phase 1: descent on A or on ½|∇A|², randomised length
  int n1 = (int)std::lround(cfg.lbfgs_min * std::pow((double)cfg.lbfgs_max / std::max(1, cfg.lbfgs_min), rng.uniform()));
  bool gradnorm = cfg.phase1 == "gradnorm" || (cfg.phase1 == "mixed" && rng.uniform() < 0.5);
  if (gradnorm) {
    std::vector<double> gA(r), hv(r);
    auto fg = [&](const double* yy, double* g) { double A = R.f(yy, gA.data(), ctx.w); if (!std::isfinite(A)) return A; R.hv(yy, gA.data(), hv.data(), ctx.w); double s = 0; for (int i = 0; i < r; i++) { s += gA[i] * gA[i]; g[i] = hv[i]; } return 0.5 * s; };
    lbfgs(r, y, fg, n1, cfg.gtol * cfg.gtol);
  } else lbfgs(r, y, fn, n1, cfg.gtol);
  OptResult nr = newton_lm(r, y, fn, fh, cfg.newton_iters, cfg.gtol);
  if (!nr.converged) { out.why = "newton did not converge"; return out; }
  R.expand(y.data(), x.data());
  if (cfg.alpha_start != 1.0) {
    Problem Pk = *P;
    for (int s = 1; s <= cfg.alpha_steps; s++) {
      Pk.alpha = cfg.alpha_start + (1.0 - cfg.alpha_start) * s / cfg.alpha_steps;
      if (!polish(Pk, S, ctx, x, cfg.newton_iters, cfg.gtol)) { out.why = "alpha continuation lost"; return out; }
    }
    P = &ctx.problem(cfg, cfg.K);
  }
  if (!certify(cfg, P, x, S, ctx, out.rec, out.why)) return out;
  out.rec.sym = symtext; out.rec.h.seed = (int64_t)cfg.seed; out.rec.h.trial = (int64_t)trial;
  out.rec.h.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  out.ok = true; return out;
}
