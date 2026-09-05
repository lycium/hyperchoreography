// Deterministic, bounded timings; compare the same compiler, flags and machine.
#include "../src/search.hpp"
#include "../src/prove.hpp"
#include <chrono>
#include <cstdio>

int main() {
  int failures = 0;
#ifdef HAVE_MPFR
  for (int n : {64, 128, 257}) {
    la::Rng rng(42); std::vector<double> a((size_t)n * n), y;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[(size_t)i * n + j] = rng.normal() + (i == j ? n : 0);
    auto start = std::chrono::steady_clock::now(); bool ok = inverse_d(n, a, y);
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    if (!ok) { std::printf("inverse n=%d: FAILED\n", n); failures++; continue; }
    double err = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
      double s = -(i == j); for (int k = 0; k < n; k++) s += a[(size_t)i * n + k] * y[(size_t)k * n + j];
      if (!std::isfinite(s)) err = INFINITY;
      else err = std::max(err, std::fabs(s));
    }
    std::printf("inverse n=%d: %.3f ms, residual %.3e, success %d\n", n, ms, err, ok);
    failures += !(err < 1e-11);
  }
#else
  std::printf("inverse benchmark requires HAVE_MPFR (the proof header owns inverse_d)\n");
#endif
  return failures != 0;
}
