// Dense linear algebra: symmetric eigensolver (tred2/tql2, LAPACK dsyevd when large), pivoted LU, RNG.
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#ifdef HAVE_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif

namespace la {

// eigenvectors in the columns of V, eigenvalues ascending in d
inline void sym_eig_ref(int n, std::vector<double>& V, std::vector<double>& d) {
  static thread_local std::vector<double> e; e.assign(n, 0.0);
  d.assign(n, 0.0);
  auto A = [&](int i, int j) -> double& { return V[(size_t)i * n + j]; };
  for (int j = 0; j < n; j++) d[j] = A(n - 1, j);
  for (int i = n - 1; i > 0; i--) {
    double scale = 0.0, h = 0.0;
    for (int k = 0; k < i; k++) scale += std::fabs(d[k]);
    if (scale == 0.0) {
      e[i] = d[i - 1];
      for (int j = 0; j < i; j++) { d[j] = A(i - 1, j); A(i, j) = 0.0; A(j, i) = 0.0; }
    } else {
      for (int k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
      double f = d[i - 1], g = std::sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g; h -= f * g; d[i - 1] = f - g;
      for (int j = 0; j < i; j++) e[j] = 0.0;
      for (int j = 0; j < i; j++) {
        f = d[j]; A(j, i) = f; g = e[j] + A(j, j) * f;
        for (int k = j + 1; k <= i - 1; k++) { g += A(k, j) * d[k]; e[k] += A(k, j) * f; }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) e[j] -= hh * d[j];
      for (int j = 0; j < i; j++) {
        f = d[j]; g = e[j];
        for (int k = j; k <= i - 1; k++) A(k, j) -= (f * e[k] + g * d[k]);
        d[j] = A(i - 1, j); A(i, j) = 0.0;
      }
    }
    d[i] = h;
  }
  for (int i = 0; i < n - 1; i++) {
    A(n - 1, i) = A(i, i); A(i, i) = 1.0;
    double h = d[i + 1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) d[k] = A(k, i + 1) / h;
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) g += A(k, i + 1) * A(k, j);
        for (int k = 0; k <= i; k++) A(k, j) -= g * d[k];
      }
    }
    for (int k = 0; k <= i; k++) A(k, i + 1) = 0.0;
  }
  for (int j = 0; j < n; j++) { d[j] = A(n - 1, j); A(n - 1, j) = 0.0; }
  A(n - 1, n - 1) = 1.0; e[0] = 0.0;
  for (int i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0.0;
  double f = 0.0, tst1 = 0.0;
  const double eps = 2.220446049250313e-16;
  for (int l = 0; l < n; l++) {
    tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
    int m = l;
    while (m < n) { if (std::fabs(e[m]) <= eps * tst1) break; m++; }
    if (m > l) {
      int iter = 0;
      do {
        iter++;
        double g = d[l], p = (d[l + 1] - g) / (2.0 * e[l]), r = std::hypot(p, 1.0);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r); d[l + 1] = e[l] * (p + r);
        double dl1 = d[l + 1], h = g - d[l];
        for (int i = l + 2; i < n; i++) d[i] -= h;
        f += h; p = d[m];
        double c = 1.0, c2 = c, c3 = c, el1 = e[l + 1], s = 0.0, s2 = 0.0;
        for (int i = m - 1; i >= l; i--) {
          c3 = c2; c2 = c; s2 = s; g = c * e[i]; h = c * p; r = std::hypot(p, e[i]);
          e[i + 1] = s * r; s = e[i] / r; c = p / r; p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          for (int k = 0; k < n; k++) { h = A(k, i + 1); A(k, i + 1) = s * A(k, i) + c * h; A(k, i) = c * A(k, i) - s * h; }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1; e[l] = s * p; d[l] = c * p;
      } while (std::fabs(e[l]) > eps * tst1 && iter < 80);
    }
    d[l] += f; e[l] = 0.0;
  }
  for (int i = 0; i < n - 1; i++) {
    int k = i; double p = d[i];
    for (int j = i + 1; j < n; j++) if (d[j] < p) { k = j; p = d[j]; }
    if (k != i) { d[k] = d[i]; d[i] = p; for (int j = 0; j < n; j++) std::swap(A(j, i), A(j, k)); }
  }
}

#ifdef HAVE_ACCELERATE
// dsyevd returns eigenvectors as columns of a column-major array; transpose into our row-major layout
inline bool sym_eig_lapack(int n, std::vector<double>& V, std::vector<double>& d) {
  static thread_local std::vector<double> work; static thread_local std::vector<__LAPACK_int> iwork;
  __LAPACK_int N = n, lda = n, lwork = 1 + 6 * n + 2 * n * n, liwork = 3 + 5 * n, info = 0;
  d.resize(n);
  if (work.size() < (size_t)lwork) work.resize(lwork);
  if (iwork.size() < (size_t)liwork) iwork.resize(liwork);
  dsyevd_("V", "L", &N, V.data(), &lda, d.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info) return false;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) std::swap(V[(size_t)i * n + j], V[(size_t)j * n + i]);
  return true;
}
#endif

inline void sym_eig(int n, std::vector<double>& V, std::vector<double>& d) {
#ifdef HAVE_ACCELERATE
  if (n >= 64) {                                   // measured break-even against tred2/tql2
    static thread_local std::vector<double> A0; A0.assign(V.begin(), V.end());
    if (sym_eig_lapack(n, V, d)) return;
    V.assign(A0.begin(), A0.end());                // dsyevd leaves V undefined when it fails
  }
#endif
  sym_eig_ref(n, V, d);
}

// A x = b by partial-pivot elimination; false if singular
template <class T>
bool lu_solve(int n, std::vector<T> A, std::vector<T>& b) {
  using std::abs;
  for (int c = 0; c < n; c++) {
    int piv = c; T best = abs(A[(size_t)c * n + c]);
    for (int r = c + 1; r < n; r++) { T a = abs(A[(size_t)r * n + c]); if (a > best) { best = a; piv = r; } }
    if (!(best > T(0))) return false;
    if (piv != c) { for (int j = 0; j < n; j++) std::swap(A[(size_t)c * n + j], A[(size_t)piv * n + j]); std::swap(b[c], b[piv]); }
    for (int r = c + 1; r < n; r++) {
      T f = A[(size_t)r * n + c] / A[(size_t)c * n + c];
      if (f == T(0)) continue;
      for (int j = c; j < n; j++) A[(size_t)r * n + j] -= f * A[(size_t)c * n + j];
      b[r] -= f * b[c];
    }
  }
  for (int r = n - 1; r >= 0; r--) {
    T s = b[r];
    for (int j = r + 1; j < n; j++) s -= A[(size_t)r * n + j] * b[j];
    b[r] = s / A[(size_t)r * n + r];
  }
  return true;
}

// singular values of a d×d matrix, descending
inline std::vector<double> singular_values(int d, const std::vector<double>& C) {
  std::vector<double> G((size_t)d * d, 0.0), w;
  for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) { double s = 0; for (int k = 0; k < d; k++) s += C[(size_t)k * d + i] * C[(size_t)k * d + j]; G[(size_t)i * d + j] = s; }
  sym_eig(d, G, w);
  std::vector<double> sv(d);
  for (int i = 0; i < d; i++) sv[i] = std::sqrt(std::max(0.0, w[d - 1 - i]));
  return sv;
}

inline int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a < 0 ? -a : a; }

// splitmix64
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() { uint64_t z = (s += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
  double uniform() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
  double normal() { double u = uniform(), v = uniform(); return std::sqrt(-2.0 * std::log(u + 1e-300)) * std::cos(6.283185307179586 * v); }
  int below(int n) { return (int)(next() % (uint64_t)n); }
};

} // namespace la
