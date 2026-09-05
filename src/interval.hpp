// Interval scalar on MPFR with outward rounding: the same in-place primitives as double/mpreal, so the Taylor
// recurrences run unchanged on intervals. Every result encloses the exact result of the same operation on any
// members of the operands; an operation that could touch a singularity (division by an interval holding 0,
// sqrt or pow of one reaching ≤ 0) yields the whole line, which no inclusion test can pass.
#pragma once
#ifdef HAVE_MPFR
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
#error "Validated arithmetic requires finite checks and outward bounds; disable fast-math."
#endif
#include <mpfr.h>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdio>

// Upper bounds on nonnegative binary64 arithmetic, including underflow. These bounds are used by
// inclusion tests, not just step-size heuristics. Do not compile this code with -ffast-math.
inline double add_up(double a, double b) {
  if (a == 0) return b; if (b == 0) return a;
  return std::nextafter(a + b, INFINITY);
}
inline double mul_up(double a, double b) {
  if (a == 0 || b == 0) return 0;
  return std::nextafter(a * b, INFINITY);
}

struct ival {
  mpfr_t lo, hi;
  static mpfr_prec_t& prec() { static mpfr_prec_t p = 256; return p; }
  ival() { mpfr_init2(lo, prec()); mpfr_init2(hi, prec()); mpfr_set_zero(lo, 1); mpfr_set_zero(hi, 1); }
  ival(double x) : ival() { mpfr_set_d(lo, x, MPFR_RNDD); mpfr_set_d(hi, x, MPFR_RNDU); }
  ival(int x) : ival() { mpfr_set_si(lo, x, MPFR_RNDD); mpfr_set_si(hi, x, MPFR_RNDU); }
  ival(const mpfr_t x) : ival() { mpfr_set(lo, x, MPFR_RNDD); mpfr_set(hi, x, MPFR_RNDU); }
  ival(double a, double b) : ival() { mpfr_set_d(lo, a, MPFR_RNDD); mpfr_set_d(hi, b, MPFR_RNDU); }
  ival(const ival& o) { mpfr_init2(lo, mpfr_get_prec(o.lo)); mpfr_init2(hi, mpfr_get_prec(o.hi)); mpfr_set(lo, o.lo, MPFR_RNDD); mpfr_set(hi, o.hi, MPFR_RNDU); }
  ival(ival&& o) noexcept : ival() { mpfr_swap(lo, o.lo); mpfr_swap(hi, o.hi); }
  ~ival() { mpfr_clear(lo); mpfr_clear(hi); }
  ival& operator=(const ival& o) { if (this != &o) { mpfr_set(lo, o.lo, MPFR_RNDD); mpfr_set(hi, o.hi, MPFR_RNDU); } return *this; }
  ival& operator=(ival&& o) noexcept { mpfr_swap(lo, o.lo); mpfr_swap(hi, o.hi); return *this; }
  static ival pi() { ival r; mpfr_const_pi(r.lo, MPFR_RNDD); mpfr_const_pi(r.hi, MPFR_RNDU); return r; }
  static ival entire() { ival r; mpfr_set_inf(r.lo, -1); mpfr_set_inf(r.hi, 1); return r; }
  static ival sym(double m) { return ival(-m, m); }                 // [−m, m]
  double mid() const { return 0.5 * (mpfr_get_d(lo, MPFR_RNDN) + mpfr_get_d(hi, MPFR_RNDN)); }
  double mag() const { return finite() ? std::max(std::fabs(mpfr_get_d(lo, MPFR_RNDA)), std::fabs(mpfr_get_d(hi, MPFR_RNDA))) : INFINITY; }   // ≥ sup|x|
  double wid() const { mpfr_t t; mpfr_init2(t, 64); mpfr_sub(t, hi, lo, MPFR_RNDU); double w = mpfr_get_d(t, MPFR_RNDU); mpfr_clear(t); return w; }
  bool finite() const { return mpfr_number_p(lo) && mpfr_number_p(hi) && mpfr_cmp(lo, hi) <= 0; }
  bool contains(const ival& o) const { return mpfr_cmp(lo, o.lo) <= 0 && mpfr_cmp(o.hi, hi) <= 0; }
  bool interior(const ival& o) const { return mpfr_cmp(lo, o.lo) < 0 && mpfr_cmp(o.hi, hi) < 0; }   // o ⊂ int(this)
  bool positive() const { return mpfr_sgn(lo) > 0; }
  void hull(const ival& o) { if (mpfr_cmp(o.lo, lo) < 0) mpfr_set(lo, o.lo, MPFR_RNDD); if (mpfr_cmp(o.hi, hi) > 0) mpfr_set(hi, o.hi, MPFR_RNDU); }
  void inflate(double rel, double abs_) {                           // widen by rel·mag + abs on each side
    ival m(add_up(mul_up(rel, mag()), abs_)); mpfr_sub(lo, lo, m.hi, MPFR_RNDD); mpfr_add(hi, hi, m.hi, MPFR_RNDU); }
  std::string str(int digits = 20) const {
    char buf[32]; std::snprintf(buf, sizeof buf, "%%.%dRDe", std::max(1, digits) - 1);
    char* a = nullptr; char* b = nullptr; mpfr_asprintf(&a, buf, lo);
    std::snprintf(buf, sizeof buf, "%%.%dRUe", std::max(1, digits) - 1); mpfr_asprintf(&b, buf, hi);
    std::string s = std::string("[") + a + ", " + b + "]"; mpfr_free_str(a); mpfr_free_str(b); return s; }
};

// per-thread scratch so every operation is alias-safe and allocation-free
struct ival_scratch { mpfr_t a, b, c, d; ival inv; ival_scratch() { for (mpfr_ptr p : {a, b, c, d}) mpfr_init2(p, ival::prec()); }
  ~ival_scratch() { for (mpfr_ptr p : {a, b, c, d}) mpfr_clear(p); } };
inline ival_scratch& iscratch() { static thread_local ival_scratch s;
  if (mpfr_get_prec(s.a) != ival::prec()) for (mpfr_ptr p : {s.a, s.b, s.c, s.d, s.inv.lo, s.inv.hi}) mpfr_set_prec(p, ival::prec()); return s; }

inline double expm1_up(double x) {
  ival_scratch& s = iscratch(); mpfr_set_d(s.a, x, MPFR_RNDN);
  mpfr_expm1(s.a, s.a, MPFR_RNDU); return mpfr_get_d(s.a, MPFR_RNDU);
}

inline double to_double(const ival& x) { return x.mid(); }
inline double log2abs(const ival& x) { double m = x.mag(); return m == 0 ? -1e300 : std::log2(m); }
inline long prec_bits(const ival&) { return (long)ival::prec(); }
inline void set_zero(ival& r) { mpfr_set_zero(r.lo, 1); mpfr_set_zero(r.hi, 1); }
inline void set_d(ival& r, double v) { mpfr_set_d(r.lo, v, MPFR_RNDD); mpfr_set_d(r.hi, v, MPFR_RNDU); }
inline void set(ival& r, const ival& v) { r = v; }
inline void add(ival& r, const ival& a, const ival& b) { ival_scratch& s = iscratch(); mpfr_add(s.a, a.lo, b.lo, MPFR_RNDD); mpfr_add(s.b, a.hi, b.hi, MPFR_RNDU); mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
inline void sub(ival& r, const ival& a, const ival& b) { ival_scratch& s = iscratch(); mpfr_sub(s.a, a.lo, b.hi, MPFR_RNDD); mpfr_sub(s.b, a.hi, b.lo, MPFR_RNDU); mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
// products by sign case: two multiplications unless both operands straddle zero
inline void mul_into(mpfr_ptr lo, mpfr_ptr hi, const ival& a, const ival& b, ival_scratch& s) {
  const int sa = mpfr_sgn(a.lo) >= 0 ? 1 : mpfr_sgn(a.hi) <= 0 ? -1 : 0, sb = mpfr_sgn(b.lo) >= 0 ? 1 : mpfr_sgn(b.hi) <= 0 ? -1 : 0;
  if (sa == 0 && sb == 0) {
    mpfr_mul(s.c, a.lo, b.hi, MPFR_RNDD); mpfr_mul(s.d, a.hi, b.lo, MPFR_RNDD); if (mpfr_cmp(s.d, s.c) < 0) mpfr_swap(s.c, s.d); mpfr_swap(lo, s.c);
    mpfr_mul(s.c, a.lo, b.lo, MPFR_RNDU); mpfr_mul(s.d, a.hi, b.hi, MPFR_RNDU); if (mpfr_cmp(s.d, s.c) > 0) mpfr_swap(s.c, s.d); mpfr_swap(hi, s.c); return; }
  mpfr_srcptr l1, l2, h1, h2;
  if (sa > 0)      { if (sb > 0) { l1 = a.lo; l2 = b.lo; h1 = a.hi; h2 = b.hi; } else if (sb < 0) { l1 = a.hi; l2 = b.lo; h1 = a.lo; h2 = b.hi; } else { l1 = a.hi; l2 = b.lo; h1 = a.hi; h2 = b.hi; } }
  else if (sa < 0) { if (sb > 0) { l1 = a.lo; l2 = b.hi; h1 = a.hi; h2 = b.lo; } else if (sb < 0) { l1 = a.hi; l2 = b.hi; h1 = a.lo; h2 = b.lo; } else { l1 = a.lo; l2 = b.hi; h1 = a.lo; h2 = b.lo; } }
  else             { if (sb > 0) { l1 = a.lo; l2 = b.hi; h1 = a.hi; h2 = b.hi; } else { l1 = a.hi; l2 = b.lo; h1 = a.lo; h2 = b.lo; } }
  mpfr_mul(s.c, l1, l2, MPFR_RNDD); mpfr_mul(s.d, h1, h2, MPFR_RNDU); mpfr_swap(lo, s.c); mpfr_swap(hi, s.d);
}
inline void mul(ival& r, const ival& a, const ival& b) { ival_scratch& s = iscratch(); mul_into(s.a, s.b, a, b, s); mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
inline void mul_d(ival& r, const ival& a, double c) { ival_scratch& s = iscratch();
  if (c >= 0) { mpfr_mul_d(s.a, a.lo, c, MPFR_RNDD); mpfr_mul_d(s.b, a.hi, c, MPFR_RNDU); } else { mpfr_mul_d(s.a, a.hi, c, MPFR_RNDD); mpfr_mul_d(s.b, a.lo, c, MPFR_RNDU); }
  mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
inline void div(ival& r, const ival& a, const ival& b) {
  if (mpfr_sgn(b.lo) <= 0 && mpfr_sgn(b.hi) >= 0) { r = ival::entire(); return; }
  ival_scratch& s = iscratch(); mpfr_ui_div(s.inv.lo, 1, b.hi, MPFR_RNDD); mpfr_ui_div(s.inv.hi, 1, b.lo, MPFR_RNDU);
  mul_into(s.a, s.b, a, s.inv, s); mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
inline void div_ui(ival& r, const ival& a, unsigned k) { mpfr_div_ui(r.lo, a.lo, k, MPFR_RNDD); mpfr_div_ui(r.hi, a.hi, k, MPFR_RNDU); }
inline void fma_add(ival& acc, const ival& a, const ival& b) { ival_scratch& s = iscratch(); mul_into(s.a, s.b, a, b, s); mpfr_add(acc.lo, acc.lo, s.a, MPFR_RNDD); mpfr_add(acc.hi, acc.hi, s.b, MPFR_RNDU); }
inline void fma_sub(ival& acc, const ival& a, const ival& b) { ival_scratch& s = iscratch(); mul_into(s.a, s.b, a, b, s); mpfr_sub(acc.lo, acc.lo, s.b, MPFR_RNDD); mpfr_sub(acc.hi, acc.hi, s.a, MPFR_RNDU); }
// Exact binary64 preconditioners need no temporary MPFR interval per matrix product.
inline void fma_add_d(ival& acc, double a, const ival& b) { ival_scratch& s = iscratch();
  mpfr_mul_d(s.a, a >= 0 ? b.lo : b.hi, a, MPFR_RNDD); mpfr_mul_d(s.b, a >= 0 ? b.hi : b.lo, a, MPFR_RNDU);
  mpfr_add(acc.lo, acc.lo, s.a, MPFR_RNDD); mpfr_add(acc.hi, acc.hi, s.b, MPFR_RNDU); }
inline void fma_sub_d(ival& acc, double a, const ival& b) { ival_scratch& s = iscratch();
  mpfr_mul_d(s.a, a >= 0 ? b.lo : b.hi, a, MPFR_RNDD); mpfr_mul_d(s.b, a >= 0 ? b.hi : b.lo, a, MPFR_RNDU);
  mpfr_sub(acc.lo, acc.lo, s.b, MPFR_RNDD); mpfr_sub(acc.hi, acc.hi, s.a, MPFR_RNDU); }
inline void sqrt_(ival& r, const ival& a) { if (mpfr_sgn(a.lo) <= 0) { r = ival::entire(); return; } mpfr_sqrt(r.lo, a.lo, MPFR_RNDD); mpfr_sqrt(r.hi, a.hi, MPFR_RNDU); }
inline void pow_d(ival& r, const ival& a, double e) {                // a > 0; monotone in a, direction by the sign of e
  if (mpfr_sgn(a.lo) <= 0) { r = ival::entire(); return; }
  ival_scratch& s = iscratch(); mpfr_set_d(s.c, e, MPFR_RNDN);
  if (e >= 0) { mpfr_pow(s.a, a.lo, s.c, MPFR_RNDD); mpfr_pow(s.b, a.hi, s.c, MPFR_RNDU); } else { mpfr_pow(s.a, a.hi, s.c, MPFR_RNDD); mpfr_pow(s.b, a.lo, s.c, MPFR_RNDU); }
  mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
inline void add_inplace(ival& r, const ival& a) { mpfr_add(r.lo, r.lo, a.lo, MPFR_RNDD); mpfr_add(r.hi, r.hi, a.hi, MPFR_RNDU); }
inline void sub_inplace(ival& r, const ival& a) { ival_scratch& s = iscratch(); mpfr_sub(s.a, r.lo, a.hi, MPFR_RNDD); mpfr_sub(s.b, r.hi, a.lo, MPFR_RNDU); mpfr_swap(r.lo, s.a); mpfr_swap(r.hi, s.b); }
inline ival operator+(const ival& a, const ival& b) { ival r; add(r, a, b); return r; }
inline ival operator-(const ival& a, const ival& b) { ival r; sub(r, a, b); return r; }
inline ival operator*(const ival& a, const ival& b) { ival r; mul(r, a, b); return r; }
inline ival operator/(const ival& a, const ival& b) { ival r; div(r, a, b); return r; }
inline ival& operator+=(ival& a, const ival& b) { add_inplace(a, b); return a; }
inline ival& operator-=(ival& a, const ival& b) { sub_inplace(a, b); return a; }
inline ival& operator*=(ival& a, const ival& b) { mul(a, a, b); return a; }
inline ival& operator/=(ival& a, const ival& b) { div(a, a, b); return a; }
inline ival operator*(const ival& a, double c) { ival r; mul_d(r, a, c); return r; }
inline ival operator-(const ival& a) { ival r; mul_d(r, a, -1.0); return r; }
inline ival sqrt(const ival& a) { ival r; sqrt_(r, a); return r; }
inline ival pow(const ival& a, double e) { ival r; pow_d(r, a, e); return r; }
// h^k for an interval step h ≥ 0 and the k-th power as a plain interval product
inline ival ipow(const ival& h, int k) { ival r(1); for (int i = 0; i < k; i++) mul(r, r, h); return r; }
#endif
