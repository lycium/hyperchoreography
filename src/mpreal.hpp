// Minimal MPFR wrapper and the in-place scalar primitives used by the hot loops.
#pragma once
#include <cmath>
#ifdef HAVE_MPFR
#include <mpfr.h>
#include <string>
#include <cmath>
#include <utility>

struct mpreal {
  mpfr_t v;
  static mpfr_prec_t& default_prec() { static mpfr_prec_t p = 256; return p; }
  static void set_default_prec(mpfr_prec_t p) { default_prec() = p; }
  mpreal() { mpfr_init2(v, default_prec()); mpfr_set_zero(v, 1); }
  mpreal(double x) { mpfr_init2(v, default_prec()); mpfr_set_d(v, x, MPFR_RNDN); }
  mpreal(int x) { mpfr_init2(v, default_prec()); mpfr_set_si(v, x, MPFR_RNDN); }
  mpreal(long x) { mpfr_init2(v, default_prec()); mpfr_set_si(v, x, MPFR_RNDN); }
  explicit mpreal(const char* s) { mpfr_init2(v, default_prec()); mpfr_set_str(v, s, 10, MPFR_RNDN); }
  explicit mpreal(const std::string& s) : mpreal(s.c_str()) {}
  mpreal(const mpreal& o) { mpfr_init2(v, mpfr_get_prec(o.v)); mpfr_set(v, o.v, MPFR_RNDN); }
  mpreal(mpreal&& o) noexcept { mpfr_init2(v, default_prec()); mpfr_swap(v, o.v); }
  ~mpreal() { mpfr_clear(v); }
  mpreal& operator=(const mpreal& o) { if (this != &o) { mpfr_set_prec(v, mpfr_get_prec(o.v)); mpfr_set(v, o.v, MPFR_RNDN); } return *this; }
  mpreal& operator=(mpreal&& o) noexcept { mpfr_swap(v, o.v); return *this; }
  mpreal& operator=(double x) { mpfr_set_d(v, x, MPFR_RNDN); return *this; }

  mpreal& operator+=(const mpreal& o) { mpfr_add(v, v, o.v, MPFR_RNDN); return *this; }
  mpreal& operator-=(const mpreal& o) { mpfr_sub(v, v, o.v, MPFR_RNDN); return *this; }
  mpreal& operator*=(const mpreal& o) { mpfr_mul(v, v, o.v, MPFR_RNDN); return *this; }
  mpreal& operator/=(const mpreal& o) { mpfr_div(v, v, o.v, MPFR_RNDN); return *this; }
  mpreal& operator*=(double x) { mpfr_mul_d(v, v, x, MPFR_RNDN); return *this; }
  mpreal& operator/=(double x) { mpfr_div_d(v, v, x, MPFR_RNDN); return *this; }
  mpreal& operator*=(int x) { mpfr_mul_si(v, v, x, MPFR_RNDN); return *this; }
  mpreal& operator/=(int x) { mpfr_div_si(v, v, x, MPFR_RNDN); return *this; }
  mpreal operator-() const { mpreal r(*this); mpfr_neg(r.v, r.v, MPFR_RNDN); return r; }

  friend mpreal operator+(mpreal a, const mpreal& b) { a += b; return a; }
  friend mpreal operator-(mpreal a, const mpreal& b) { a -= b; return a; }
  friend mpreal operator*(mpreal a, const mpreal& b) { a *= b; return a; }
  friend mpreal operator/(mpreal a, const mpreal& b) { a /= b; return a; }
  friend mpreal operator*(mpreal a, double b) { a *= b; return a; }
  friend mpreal operator*(double b, mpreal a) { a *= b; return a; }
  friend mpreal operator/(mpreal a, double b) { a /= b; return a; }
  friend mpreal operator*(mpreal a, int b) { a *= b; return a; }
  friend mpreal operator*(int b, mpreal a) { a *= b; return a; }
  friend mpreal operator/(mpreal a, int b) { a /= b; return a; }
  friend mpreal operator+(mpreal a, double b) { mpfr_add_d(a.v, a.v, b, MPFR_RNDN); return a; }
  friend mpreal operator-(mpreal a, double b) { mpfr_sub_d(a.v, a.v, b, MPFR_RNDN); return a; }
  friend mpreal operator/(double b, const mpreal& a) { mpreal r(b); r /= a; return r; }

  friend bool operator<(const mpreal& a, const mpreal& b) { return mpfr_cmp(a.v, b.v) < 0; }
  friend bool operator>(const mpreal& a, const mpreal& b) { return mpfr_cmp(a.v, b.v) > 0; }
  friend bool operator<=(const mpreal& a, const mpreal& b) { return mpfr_cmp(a.v, b.v) <= 0; }
  friend bool operator>=(const mpreal& a, const mpreal& b) { return mpfr_cmp(a.v, b.v) >= 0; }
  friend bool operator==(const mpreal& a, const mpreal& b) { return mpfr_cmp(a.v, b.v) == 0; }
  friend bool operator!=(const mpreal& a, const mpreal& b) { return mpfr_cmp(a.v, b.v) != 0; }
  friend bool operator<(const mpreal& a, double b) { return mpfr_cmp_d(a.v, b) < 0; }
  friend bool operator>(const mpreal& a, double b) { return mpfr_cmp_d(a.v, b) > 0; }

  friend mpreal sqrt(const mpreal& a) { mpreal r; mpfr_sqrt(r.v, a.v, MPFR_RNDN); return r; }
  friend mpreal abs(const mpreal& a) { mpreal r; mpfr_abs(r.v, a.v, MPFR_RNDN); return r; }
  friend mpreal fabs(const mpreal& a) { return abs(a); }
  friend mpreal exp(const mpreal& a) { mpreal r; mpfr_exp(r.v, a.v, MPFR_RNDN); return r; }
  friend mpreal log(const mpreal& a) { mpreal r; mpfr_log(r.v, a.v, MPFR_RNDN); return r; }
  friend mpreal sin(const mpreal& a) { mpreal r; mpfr_sin(r.v, a.v, MPFR_RNDN); return r; }
  friend mpreal cos(const mpreal& a) { mpreal r; mpfr_cos(r.v, a.v, MPFR_RNDN); return r; }
  friend mpreal pow(const mpreal& a, const mpreal& b) { mpreal r; mpfr_pow(r.v, a.v, b.v, MPFR_RNDN); return r; }
  friend mpreal pow(const mpreal& a, double b) { return pow(a, mpreal(b)); }
  friend mpreal ldexp2(const mpreal& a, long e) { mpreal r; mpfr_mul_2si(r.v, a.v, e, MPFR_RNDN); return r; }
  static mpreal pi() { mpreal r; mpfr_const_pi(r.v, MPFR_RNDN); return r; }
  static mpreal two_pow(long e) { mpreal r(1.0); mpfr_mul_2si(r.v, r.v, e, MPFR_RNDN); return r; }

  double to_double() const { return mpfr_get_d(v, MPFR_RNDN); }
  explicit operator double() const { return to_double(); }
  bool is_zero() const { return mpfr_zero_p(v); }
  double log2abs() const {
    if (mpfr_zero_p(v)) return -1e300;
    long e; double m = mpfr_get_d_2exp(&e, v, MPFR_RNDN);
    return std::log2(std::fabs(m)) + (double)e;
  }
  std::string str(int digits = 40) const {
    char buf[64]; std::snprintf(buf, sizeof buf, "%%.%dRe", digits - 1);
    char* s = nullptr; mpfr_asprintf(&s, buf, v);
    std::string out(s); mpfr_free_str(s); return out;
  }
};
inline double to_double(const mpreal& x) { return x.to_double(); }
inline double log2abs(const mpreal& x) { return x.log2abs(); }
inline long prec_bits(const mpreal& x) { return (long)mpfr_get_prec(x.v); }
#endif
inline double to_double(double x) { return x; }
inline double log2abs(double x) { return x == 0 ? -1e300 : std::log2(std::fabs(x)); }
inline double ldexp2(double x, long e) { return std::ldexp(x, (int)e); }
inline long prec_bits(double) { return 53; }

// in-place primitives: no temporaries, no allocation
inline void set_zero(double& r) { r = 0.0; }
inline void set_d(double& r, double v) { r = v; }
inline void set(double& r, double v) { r = v; }
inline void add(double& r, double a, double b) { r = a + b; }
inline void sub(double& r, double a, double b) { r = a - b; }
inline void mul(double& r, double a, double b) { r = a * b; }
inline void mul_d(double& r, double a, double b) { r = a * b; }
inline void div(double& r, double a, double b) { r = a / b; }
inline void div_ui(double& r, double a, unsigned k) { r = a / k; }
inline void fma_add(double& acc, double a, double b) { acc += a * b; }
inline void fma_sub(double& acc, double a, double b) { acc -= a * b; }
inline void sqrt_(double& r, double a) { r = std::sqrt(a); }
inline void pow_d(double& r, double a, double e) { r = std::pow(a, e); }
inline void add_inplace(double& r, double a) { r += a; }
inline void sub_inplace(double& r, double a) { r -= a; }
#ifdef HAVE_MPFR
inline void set_zero(mpreal& r) { mpfr_set_zero(r.v, 1); }
inline void set_d(mpreal& r, double v) { mpfr_set_d(r.v, v, MPFR_RNDN); }
inline void set(mpreal& r, const mpreal& v) { mpfr_set(r.v, v.v, MPFR_RNDN); }
inline void add(mpreal& r, const mpreal& a, const mpreal& b) { mpfr_add(r.v, a.v, b.v, MPFR_RNDN); }
inline void sub(mpreal& r, const mpreal& a, const mpreal& b) { mpfr_sub(r.v, a.v, b.v, MPFR_RNDN); }
inline void mul(mpreal& r, const mpreal& a, const mpreal& b) { mpfr_mul(r.v, a.v, b.v, MPFR_RNDN); }
inline void mul_d(mpreal& r, const mpreal& a, double b) { mpfr_mul_d(r.v, a.v, b, MPFR_RNDN); }
inline void div(mpreal& r, const mpreal& a, const mpreal& b) { mpfr_div(r.v, a.v, b.v, MPFR_RNDN); }
inline void div_ui(mpreal& r, const mpreal& a, unsigned k) { mpfr_div_ui(r.v, a.v, k, MPFR_RNDN); }
inline void fma_add(mpreal& acc, const mpreal& a, const mpreal& b) { mpfr_fma(acc.v, a.v, b.v, acc.v, MPFR_RNDN); }
inline void fma_sub(mpreal& acc, const mpreal& a, const mpreal& b) { mpfr_fms(acc.v, a.v, b.v, acc.v, MPFR_RNDN); mpfr_neg(acc.v, acc.v, MPFR_RNDN); }
inline void sqrt_(mpreal& r, const mpreal& a) { mpfr_sqrt(r.v, a.v, MPFR_RNDN); }
inline void pow_d(mpreal& r, const mpreal& a, double e) { mpreal ee(e); mpfr_pow(r.v, a.v, ee.v, MPFR_RNDN); }
inline void add_inplace(mpreal& r, const mpreal& a) { mpfr_add(r.v, r.v, a.v, MPFR_RNDN); }
inline void sub_inplace(mpreal& r, const mpreal& a) { mpfr_sub(r.v, r.v, a.v, MPFR_RNDN); }
#endif
