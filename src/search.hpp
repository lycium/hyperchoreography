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

struct Config {
  int N = 3, d = 2, K = 16, Kmax = 64;
  double alpha_start = 1.0; int alpha_steps = 8;
  std::string sym = "none";                    // none | random | generator DSL
  std::string phase1 = "mixed";                // action | gradnorm | mixed
  int threads = 0; uint64_t seed = 1; long trials = LONG_MAX; double minutes = 1e30;
  std::string out = "catalog.bin";
  int lbfgs_min = 20, lbfgs_max = 400, newton_iters = 60; double gtol = 1e-10, ret_tol = 1e-8, ret_reject = 1e-1;
  int K0min = 2, K0max = 6; double minsep = 2e-3; int min_deff = 1;
  const std::vector<Record>* seeds = nullptr; double kick_min = 0.02, kick_max = 0.5;
  int Ms = 2048, Kout_max = 512; double shoot_tol = 1e-12, ret_double = 1e-4;
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
    if (d != 7) throw std::runtime_error("--omega g2: needs d = 7");
    std::stringstream ss(t.substr(3)); std::string a, b; std::getline(ss, a, ','); std::getline(ss, b, ',');
    return g2_omega(std::stod(a), b.empty() ? 0.0 : std::stod(b));
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

// state of body k at t = 0
inline void initial_state(const Problem& P, const double* x, std::vector<double>& pos, std::vector<double>& vel) {
  const int N = P.N, d = P.d; pos.assign((size_t)N * d, 0.0); vel.assign((size_t)N * d, 0.0);
  std::vector<double> xd; deriv_coeffs(P, x, xd);
  for (int k = 0; k < N; k++) { double t = 2.0 * PI * k / N;
    for (int mu = 0; mu < P.nm; mu++) { double c = std::cos(P.modes[mu] * t), s = std::sin(P.modes[mu] * t);
      for (int a = 0; a < d; a++) { pos[k * d + a] += c * x[(2 * mu) * d + a] + s * x[(2 * mu + 1) * d + a]; vel[k * d + a] += c * xd[(2 * mu) * d + a] + s * xd[(2 * mu + 1) * d + a]; } } }
  if (!P.Om.empty()) for (int k = 0; k < N; k++) for (int a = 0; a < d; a++)              // Q̇ = q̇ + Ωq
    for (int b = 0; b < d; b++) vel[k * d + a] += P.Om[(size_t)a * d + b] * pos[k * d + b];
}
// relative shift residual of the Fourier loop under the ODE
inline double return_error(const Problem& P, const double* x, NBody<double>& nb) {
  std::vector<double> pos, vel; initial_state(P, x, pos, vel);
  double sc = 1.0; for (double v : pos) sc = std::max(sc, std::fabs(v)); for (double v : vel) sc = std::max(sc, std::fabs(v));
  return chore_residual(nb, pos, vel, 1e-16, P.gshift()) / sc;
}
inline double return_error(const Problem& P, const double* x) { NBody<double> nb(P.N, P.d, P.alpha, 22); return return_error(P, x, nb); }

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
  auto set_mode = [&](int m, int c, double cc, double ss) { auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) return; int mu = (int)(it - P.modes.begin()); x[(2 * mu) * d + c] = cc; x[(2 * mu + 1) * d + c] = ss; };
  double sgn = rng.uniform() < 0.5 ? 1 : -1, amp = 0.08 * std::pow(15.0, rng.uniform());
  set_mode(m1, 0, 1.0, 0.0); set_mode(m1, 1, 0.0, sgn);
  double ph = 2 * PI * rng.uniform(); set_mode(m2, 2, amp * std::cos(ph), amp * std::sin(ph));
  for (int c = 3; c < d; c++) if (rng.uniform() < 0.5) { int m3 = 1 + rng.below(std::min(P.K, 8)); if (m3 % N) set_mode(m3, c, 0.3 * amp * rng.normal(), 0.3 * amp * rng.normal()); }
}

// Resonant hyper-start: m1-fold circle in plane (0,1) plus transverse modes at the vertical resonances
// m ≈ m1·ω_k with m ≡ k·m1 (mod N), distinct patterns k drawn without replacement. The transverse block is
// (d−2)-fold degenerate, so each circularly polarised pair adds 2 to deff: deff ≤ 2 + 2(⌊N/2⌋−1) caps the
// reachable dimension, which is why d = 7 needs N ≥ 8.
inline void hyper_guess(const Problem& P, la::Rng& rng, std::vector<double>& x) {
  x.assign(P.n, 0.0); const int N = P.N, d = P.d, nk = N / 2 - 1;
  auto set_mode = [&](int m, int c, double cc, double ss) {
    auto it = std::lower_bound(P.modes.begin(), P.modes.end(), m); if (it == P.modes.end() || *it != m) return;
    int mu = (int)(it - P.modes.begin()); x[(2 * mu) * d + c] = cc; x[(2 * mu + 1) * d + c] = ss; };
  int m1 = 1; for (int t = 0; t < 64; t++) { int a = 1 + rng.below(std::max(1, P.K / 2)); if (la::gcd(a, N) == 1) { m1 = a; break; } }
  set_mode(m1, 0, 1.0, 0.0); set_mode(m1, 1, 0.0, rng.uniform() < 0.5 ? 1 : -1);
  if (nk <= 0) return;
  std::vector<int> ks(nk); for (int i = 0; i < nk; i++) ks[i] = i + 2;
  for (int i = nk - 1; i > 0; i--) std::swap(ks[i], ks[rng.below(i + 1)]);
  for (int c = 2, i = 0; c < d && i < 4 * nk; i++) {
    int k = ks[i % nk]; double tgt = m1 * ngon_vertical_freq(N, k);
    int base = (k * m1) % N, m = base + N * (int)std::lround((tgt - base) / N);
    if (m < 1 || m % N == 0 || m > P.K || std::fabs(m - tgt) > 0.25 * tgt) continue;
    double amp = 0.03 * std::pow(20.0, rng.uniform()), ph = 2 * PI * rng.uniform();
    set_mode(m, c, amp * std::cos(ph), amp * std::sin(ph));
    if (c + 1 < d) { double s2 = rng.uniform() < 0.5 ? 1 : -1;
      set_mode(m, c + 1, -s2 * amp * std::sin(ph), s2 * amp * std::cos(ph)); c += 2; } else c++;
  }
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

// orthonormal gauge directions at x: time shift plus the rotations that survive the rotating frame
inline int gauge_basis(const Problem& P, const double* x, std::vector<double>& G) {
  std::vector<std::vector<double>> gs(1); deriv_coeffs(P, x, gs[0]);
  for (int a = 0; a < P.d; a++) for (int b = a + 1; b < P.d; b++) {
    if (!P.Om.empty()) { double cm = 0;                       // in a rotating frame only [L, Ω] = 0 stays gauge
      for (int i = 0; i < P.d; i++) for (int j = 0; j < P.d; j++) { double lo = 0, ol = 0;
        if (i == a) lo -= P.Om[(size_t)b * P.d + j];
        if (i == b) lo += P.Om[(size_t)a * P.d + j];
        if (j == b) ol -= P.Om[(size_t)i * P.d + a];
        if (j == a) ol += P.Om[(size_t)i * P.d + b];
        cm = std::max(cm, std::fabs(lo - ol)); }
      if (cm > 1e-12) continue; }
    std::vector<double> v(P.n, 0.0);
    for (int i = 0; i < P.nb; i++) { v[i * P.d + a] = -x[i * P.d + b]; v[i * P.d + b] = x[i * P.d + a]; } gs.push_back(v); }
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

// Morse index, nullity and the calibration twist — computed only for records that are new or replaced
inline void record_extras(const Config& cfg, Record& rec, Ctx& ctx) {
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
  double gn = 0; ctx.g.resize(P.n); double A0 = action_grad(P, x.data(), ctx.g.data(), ctx.w);
  if (!std::isfinite(A0)) return; for (double v : ctx.g) gn += v * v;
  rec.h.grad_norm = std::sqrt(gn);                       // residual of the truncation, before the polish
  Symmetry none; polish(P, none, ctx, x, cfg.newton_iters, cfg.gtol);   // the index is only defined at a critical point
  double A1 = action_grad(P, x.data(), nullptr, ctx.w);
  if (std::fabs(A1 - A0) < 1e-3 * std::fabs(A0) && hessian(P, x.data(), ctx.H, ctx.w))
    { Inertia I = inertia_gauge(P, x.data(), ctx.H); rec.h.morse = I.neg; rec.h.nullity = I.zero; }
}

// ODE validation, cover unwinding, shooting certification, Fourier re-extraction, canonical frame → record
inline bool certify(const Config& cfg, const Problem*& P, std::vector<double>& x, const Symmetry& S, Ctx& ctx, Record& rec, std::string& why) {
  // deff can only be lost downstream, so gate before the ODE work; 1e-12 keeps the gate permissive
  if (cfg.min_deff > 1) { std::vector<double> xc(x), sv; if (canonical_frame(P->nb, P->d, xc, sv, 1e-12) < cfg.min_deff) { why = "effective dimension below filter"; return false; } }
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
  double sc = 1.0; for (double v : O.Z) sc = std::max(sc, std::fabs(v));
  double res = shoot_newton(nb, O.Z, 1e-16, 10, cfg.shoot_tol * sc, -18, -40, ctx.sw, false, P->gshift()) / sc;
  if (!(res <= cfg.shoot_tol * 100)) { why = "shooting did not certify"; return false; }
  if (!orbit_fit(P->N, P->d, 1.0, O, ctx, cfg.Ms, cfg.Kout_max, P->Om.empty() ? nullptr : &P->Om)) { why = "orbit sampling failed"; return false; }
  O.residual = res;
  std::vector<double> sv, Rc; int deff = canonical_frame((int)O.modes.size() * 2, P->d, O.coef, sv, 1e-8, P->Om.empty() ? nullptr : &Rc);
  if (deff < cfg.min_deff) { why = "effective dimension below filter"; return false; }
  rec.h.N = P->N; rec.h.d = P->d; rec.h.K = O.K; rec.h.M = O.Ms; rec.h.alpha = 1.0; rec.modes.assign(O.modes.begin(), O.modes.end()); rec.coef = O.coef; rec.h.nm = (int32_t)O.modes.size();
  rec.h.deff = deff; rec.h.cover = cover; rec.h.action = O.action; rec.h.energy = O.energy; rec.h.energy_std = 0; rec.h.rms = O.rms; rec.h.maxr = O.maxr;
  rec.h.minsep = O.minsep; rec.h.Lnorm = O.Lnorm; rec.Lsv = O.Lsv; rec.pca = sv; rec.h.morse = -1; rec.h.nullity = -1; rec.h.grad_norm = -1; rec.h.ret_err = res;
  const int d = P->d; rec.extra.assign(Record::NEX + (P->Om.empty() ? 0 : (size_t)d * d), 0.0);
  if (!P->Om.empty()) for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { double t = 0;      // Ω → R Ω Rᵀ
    for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) t += Rc[(size_t)a * d + i] * P->Om[(size_t)i * d + j] * Rc[(size_t)b * d + j];
    rec.extra[Record::NEX + (size_t)a * d + b] = t; }
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
