// Binary catalog: "HYPCHOR1", then records = RecHdr | int32 modes[nm] | double coef[nm*2*d] | double Lsv[d] | double pca[d] | char sym[symlen]
#pragma once
#include "invariants.hpp"
#include <string>
#include <map>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <stdexcept>

struct RecHdr {                     // fixed-size record header, 152 bytes
  uint32_t magic = 0x31434552u;     // "REC1"
  uint32_t nbytes = 0;              // total record size including header
  int64_t id = -1, seed = 0, trial = -1;
  int32_t N = 0, d = 0, deff = 0, K = 0, M = 0, cover = 1, morse = -1, nullity = -1, hits = 1, nm = 0, symlen = 0, nextra = 0;
  double alpha = 1, action = 0, energy = 0, energy_std = 0, rms = 0, maxr = 0, minsep = 0, Lnorm = 0, grad_norm = 0, ret_err = -1, secs = 0;
};
static_assert(sizeof(RecHdr) == 8 + 24 + 48 + 88, "RecHdr layout");

// extra[]: 0 = calibration twist χ*, 1 = χ* / jet scale, 2 = the rung's jet order k, 3 = ‖A_k‖ / jet scale,
// 4 = rigidity defect, 5..7 reserved, then d² entries of Ω for a rotating frame.
// Legacy files have nextra = 0 and read back unchanged.
struct Record {
  RecHdr h;
  std::vector<int32_t> modes; std::vector<double> coef, Lsv, pca, extra; std::string sym;
  static constexpr int NEX = 8;
  double twist() const { return extra.empty() ? 0.0 : extra[0]; }
  double rigid() const { return extra.size() > 4 ? extra[4] : -1.0; }   // 0 = relative equilibrium, −1 = not computed
  const double* omega() const { return extra.size() >= (size_t)(NEX + h.d * h.d) ? &extra[NEX] : nullptr; }
  int64_t& id() { return h.id; }  int N() const { return h.N; }  int d() const { return h.d; }

  void set_solution(const Problem& P, const double* x) {
    h.N = P.N; h.d = P.d; h.K = P.K; h.M = P.M; h.alpha = P.alpha; h.nm = P.nm; modes.assign(P.modes.begin(), P.modes.end());
    coef.assign((size_t)P.nm * 2 * P.d, 0.0);
    for (int mu = 0; mu < P.nm; mu++) for (int a = 0; a < P.d; a++) { coef[(size_t)mu * 2 * P.d + a] = x[(2 * mu) * P.d + a]; coef[(size_t)mu * 2 * P.d + P.d + a] = x[(2 * mu + 1) * P.d + a]; }
  }
  void to_problem(Problem& P, std::vector<double>& x) const {
    int M = h.M % h.N ? 0 : h.M;                          // legacy records stored an M not divisible by N
    P.init(h.N, h.d, h.K, M, h.alpha); x.assign(P.n, 0.0); const int d = h.d;
    for (size_t k = 0; k < modes.size(); k++) { auto it = std::lower_bound(P.modes.begin(), P.modes.end(), modes[k]); if (it == P.modes.end() || *it != modes[k]) continue;
      int mu = (int)(it - P.modes.begin()); for (int a = 0; a < d; a++) { x[(2 * mu) * d + a] = coef[k * 2 * d + a]; x[(2 * mu + 1) * d + a] = coef[k * 2 * d + d + a]; } }
    if (const double* om = omega()) P.set_omega(std::vector<double>(om, om + (size_t)d * d));
  }
  std::vector<int> mode_list() const { return std::vector<int>(modes.begin(), modes.end()); }

  bool write(FILE* f) const {
    RecHdr hh = h; hh.nm = (int32_t)modes.size(); hh.symlen = (int32_t)sym.size(); hh.nextra = (int32_t)extra.size();
    hh.nbytes = (uint32_t)(sizeof(RecHdr) + modes.size() * 4 + (coef.size() + Lsv.size() + pca.size() + extra.size()) * 8 + sym.size());
    if (Lsv.size() != (size_t)h.d || pca.size() != (size_t)h.d) throw std::runtime_error("record: Lsv/pca must have d entries");
    return std::fwrite(&hh, sizeof hh, 1, f) == 1 && (modes.empty() || std::fwrite(modes.data(), 4, modes.size(), f) == modes.size())
        && std::fwrite(coef.data(), 8, coef.size(), f) == coef.size() && std::fwrite(Lsv.data(), 8, Lsv.size(), f) == Lsv.size()
        && std::fwrite(pca.data(), 8, pca.size(), f) == pca.size() && (extra.empty() || std::fwrite(extra.data(), 8, extra.size(), f) == extra.size())
        && (sym.empty() || std::fwrite(sym.data(), 1, sym.size(), f) == sym.size());
  }
  bool read(FILE* f) {
    if (std::fread(&h, sizeof h, 1, f) != 1) return false;
    if (h.magic != 0x31434552u || h.nm < 0 || h.d <= 0 || h.symlen < 0 || h.nextra < 0) return false;
    size_t nb = sizeof(RecHdr) + (size_t)h.nm * 4 + ((size_t)h.nm * 2 * h.d + 2 * (size_t)h.d + h.nextra) * 8 + h.symlen;
    if (nb != h.nbytes) return false;
    modes.resize(h.nm); coef.resize((size_t)h.nm * 2 * h.d); Lsv.resize(h.d); pca.resize(h.d); extra.resize(h.nextra); sym.resize(h.symlen);
    return (h.nm == 0 || std::fread(modes.data(), 4, h.nm, f) == (size_t)h.nm) && std::fread(coef.data(), 8, coef.size(), f) == coef.size()
        && std::fread(Lsv.data(), 8, Lsv.size(), f) == Lsv.size() && std::fread(pca.data(), 8, pca.size(), f) == pca.size()
        && (extra.empty() || std::fread(extra.data(), 8, extra.size(), f) == extra.size())
        && (h.symlen == 0 || std::fread(&sym[0], 1, h.symlen, f) == (size_t)h.symlen);
  }
  // JSON dump
  std::string to_json(bool with_coef = true) const {
    auto num = [](double v) { char b[32]; std::snprintf(b, sizeof b, "%.17g", v); return std::string(b); };
    auto arr = [&](const std::vector<double>& v) { std::string o = "["; for (size_t i = 0; i < v.size(); i++) o += (i ? "," : "") + num(v[i]); return o + "]"; };
    std::string o = "{\"id\":" + std::to_string(h.id) + ",\"N\":" + std::to_string(h.N) + ",\"d\":" + std::to_string(h.d) + ",\"deff\":" + std::to_string(h.deff) + ",\"K\":" + std::to_string(h.K) + ",\"M\":" + std::to_string(h.M);
    o += ",\"cover\":" + std::to_string(h.cover) + ",\"alpha\":" + num(h.alpha) + ",\"action\":" + num(h.action) + ",\"energy\":" + num(h.energy) + ",\"energy_std\":" + num(h.energy_std) + ",\"rms\":" + num(h.rms) + ",\"maxr\":" + num(h.maxr);
    o += ",\"minsep\":" + num(h.minsep) + ",\"Lnorm\":" + num(h.Lnorm) + ",\"Lsv\":" + arr(Lsv) + ",\"pca\":" + arr(pca) + ",\"morse\":" + std::to_string(h.morse) + ",\"nullity\":" + std::to_string(h.nullity);
    if (!extra.empty()) { o += ",\"twist\":" + num(extra[0]) + ",\"twist_rel\":" + num(extra[1]) + ",\"calib_k\":" + num(extra[2]) + ",\"jet_rel\":" + num(extra[3]);
      if (omega()) o += ",\"omega\":" + arr(std::vector<double>(extra.begin() + NEX, extra.end())); }
    o += ",\"grad_norm\":" + num(h.grad_norm) + ",\"ret_err\":" + num(h.ret_err) + ",\"sym\":\"" + sym + "\",\"seed\":" + std::to_string(h.seed) + ",\"trial\":" + std::to_string(h.trial) + ",\"hits\":" + std::to_string(h.hits) + ",\"secs\":" + num(h.secs);
    o += ",\"modes\":["; for (size_t i = 0; i < modes.size(); i++) o += (i ? "," : "") + std::to_string(modes[i]); o += "]";
    if (with_coef) o += ",\"coef\":" + arr(coef);
    return o + "}";
  }
};

struct Catalog {
  std::vector<Record> recs;
  std::multimap<double, size_t> index;      // action → record position
  int64_t next_id = 0;
  static constexpr const char* MAGIC = "HYPCHOR1";

  bool load(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) return false;
    char m[8]; if (std::fread(m, 1, 8, f) != 8 || std::memcmp(m, MAGIC, 8) != 0) { std::fclose(f); throw std::runtime_error("not a catalog file: " + path); }
    Record r; while (r.read(f)) { push(r); r = Record(); }
    std::fclose(f); return true;
  }
  size_t push(Record r) { if (r.h.id < 0) r.h.id = next_id; next_id = std::max(next_id, r.h.id + 1); index.emplace(r.h.action, recs.size()); recs.push_back(std::move(r)); return recs.size() - 1; }
  // per-mode power is invariant under O(d), time shift, reversal and relabelling — the whole equivalence
  // group — and catches many-mode loops on which the Procrustes optimisation fails to find the alignment
  static bool same_spectrum(const Record& a, const Record& b, double tol = 1e-9) {
    if (a.h.d != b.h.d || a.modes != b.modes) return false;
    double s = 0, e = 0; const int d2 = 2 * a.h.d;
    for (size_t k = 0; k < a.modes.size(); k++) { double pa = 0, pb = 0;
      for (int c = 0; c < d2; c++) { double va = a.coef[k * d2 + c], vb = b.coef[k * d2 + c]; pa += va * va; pb += vb * vb; }
      s += pa; e = std::max(e, std::fabs(pa - pb)); }
    return e <= tol * s;
  }
  // invariants within tol_inv, then Procrustes distance below tol_dist
  long find_duplicate(const Record& r, double tol_inv = 1e-4, double tol_dist = 1e-3, double* dist_out = nullptr) const {
    auto lo = index.lower_bound(r.h.action * (1 - tol_inv)), hi = index.upper_bound(r.h.action * (1 + tol_inv));
    long best = -1; double bestd = INF;
    for (auto it = lo; it != hi; ++it) {
      const Record& c = recs[it->second];
      if (c.h.N != r.h.N || std::fabs(c.h.energy - r.h.energy) > tol_inv * std::fabs(r.h.energy) || std::fabs(c.h.rms - r.h.rms) > tol_inv * r.h.rms) continue;
      double dist = loop_distance(r.h.N, r.mode_list(), r.h.d, r.coef.data(), c.mode_list(), c.h.d, c.coef.data());
      if (same_spectrum(r, c)) dist = 0;
      // Points of a continuous family are genuinely different loops — loop_distance is right to separate them
      // and χ* really does vary along the family — but the catalogue wants one entry with a hit count, not a
      // sampling of the family. Action and energy are exactly constant along a family, so match them to
      // round-off; minsep is not, but it separates the case this must not fold — distinct orbits whose
      // actions happen to agree closely, which at d=7 N=10 differ in minsep by percents (and in Morse index,
      // not yet computed at this point). Heuristic, and deliberately narrow.
      if (std::fabs(c.h.action - r.h.action) <= 1e-8 * std::fabs(r.h.action) &&
          std::fabs(c.h.energy - r.h.energy) <= 1e-8 * std::fabs(r.h.energy) &&
          std::fabs(c.h.minsep - r.h.minsep) <= 1e-3 * std::fabs(r.h.minsep)) dist = 0;
      if (dist < bestd) { bestd = dist; best = (long)it->second; }
    }
    if (dist_out) *dist_out = bestd;
    return bestd < tol_dist ? best : -1;
  }
  // keep the better-resolved duplicate (smaller residual, then fewer modes)
  static bool better(const Record& cand, const Record& cur) {
    double rc = cand.h.ret_err > 0 ? cand.h.ret_err : 1e-300, ro = cur.h.ret_err > 0 ? cur.h.ret_err : 1e-300;
    if (rc < 0.1 * ro) return true; if (ro < 0.1 * rc) return false; return cand.h.K < cur.h.K;
  }
  bool absorb(size_t idx, const Record& cand, int add_hits = 1) {          // true when cand replaced the record
    Record& cur = recs[idx]; cur.h.hits += add_hits;
    if (better(cand, cur)) { int64_t id = cur.h.id; int hits = cur.h.hits; double a = cur.h.action; cur = cand; cur.h.id = id; cur.h.hits = hits;
      if (a != cur.h.action) { for (auto it = index.lower_bound(a); it != index.end() && it->first == a; ++it) if (it->second == idx) { index.erase(it); break; } index.emplace(cur.h.action, idx); }
      return true; }
    return false;
  }
  void save(const std::string& path) const {
    std::string tmp = path + ".tmp"; FILE* f = std::fopen(tmp.c_str(), "wb"); if (!f) throw std::runtime_error("cannot write " + tmp);
    std::fwrite(MAGIC, 1, 8, f); for (auto& r : recs) r.write(f); std::fclose(f);
    std::rename(tmp.c_str(), path.c_str());
  }
};

struct SearchState {                 // <catalog>.state — resumable trial counter per seed
  uint64_t magic = 0x3154415453ull;  // "STAT1"
  int64_t seed = 0, next_trial = 0, trials_done = 0, found = 0; double elapsed = 0;
  bool load(const std::string& p) { FILE* f = std::fopen(p.c_str(), "rb"); if (!f) return false; bool ok = std::fread(this, sizeof *this, 1, f) == 1 && magic == 0x3154415453ull; std::fclose(f); return ok; }
  void save(const std::string& p) const { std::string t = p + ".tmp"; FILE* f = std::fopen(t.c_str(), "wb"); if (!f) return; std::fwrite(this, sizeof *this, 1, f); std::fclose(f); std::rename(t.c_str(), p.c_str()); }
};
