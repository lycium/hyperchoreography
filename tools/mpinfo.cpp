// Which mpn assembly GMP is actually running, and what MPFR costs at prove/refine precisions.
// GMP's CPU table stops at AMD family 0x19 in 6.3.0, so a Zen 4 or 5 host -- and every MSYS2
// build, configured for a generic host -- silently falls back to the baseline k8 path.

#include <mpfr.h>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#if defined(__unix__) || defined(__APPLE__)
#include <dlfcn.h>
#define MPINFO_DLSYM 1
#endif

// Longest first: __gmpn_mod_1s_4p_cps_x86_64 does not end at its last underscore.
static const char* const CPUS[] = {
  "silvermont", "coreisbr", "coreihwl", "coreibwl", "coreinhm", "pentium4", "goldmont", "skylake",
  "x86_64", "core2", "zen3", "zen2", "atom", "nano", "zen", "k10", "bd1", "bd2", "bd3", "bd4",
  "bt1", "bt2", "k8", "fat" };

static bool in(const char* set, const std::string& s) {
  return std::strstr(set, (" " + s + " ").c_str()) != nullptr;
}
static const char* cpu_suffix(const char* sym) {
  size_t n = std::strlen(sym);
  for (const char* c : CPUS) { size_t m = std::strlen(c);
    if (n > m + 1 && sym[n - m - 1] == '_' && !std::strcmp(sym + n - m, c)) return c; }
  return nullptr;
}
static std::string gmp_march() {
  const char* f = std::strstr(__GMP_CFLAGS, "-march=");
  if (!f) return "";
  f += 7; size_t n = 0; while (f[n] && f[n] != ' ' && f[n] != '\t') n++;
  return std::string(f, n);
}

// Best of three: a boost clock that wanders is worth more here than a mean over it.
template <class F> static double ns(F&& f, long it) {
  double best = 1e300;
  for (int rep = 0; rep < 3; rep++) {
    for (long i = 0; i < it / 8; i++) f();
    auto t0 = std::chrono::steady_clock::now(); for (long i = 0; i < it; i++) f();
    double d = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / it;
    if (d < best) best = d;
  }
  return best;
}

int main() {
  std::printf("GMP   %-8s  cc: %s\n            flags: %s\n", gmp_version, __GMP_CC, __GMP_CFLAGS);
  bool generic = false, decided = false;

  // Not a fat binary: the assembly was fixed at configure time and -march is the only trace left.
  auto from_march = [&] {
    std::string m = gmp_march();
    std::printf("      mpn path: fixed at build time, not a fat binary; -march=%s\n",
                m.empty() ? "(none recorded)" : m.c_str());
    if (m.empty() || in(" x86-64 x86-64-v2 x86-64-v3 x86-64-v4 nocona k8 athlon64 core2 generic ", m))
      generic = decided = true;
    else if (m == "native")
      std::printf("                -march=native is the compiler's idea of the CPU, not GMP's.\n");
    else decided = true;
  };
#ifdef MPINFO_DLSYM
  { mp_limb_t a[8] = {1,2,3,4,5,6,7,8}, b[8] = {8,7,6,5,4,3,2,1}, r[16];
    mpn_add_n(r, a, b, 8); mpn_mul_n(r, a, b, 8); mpn_sqr(r, a, 8);   // force the fat dispatch
    void** cpuvec = (void**)dlsym(RTLD_DEFAULT, "__gmpn_cpuvec");
    if (!cpuvec) from_march();
    else {
      std::string seen = " ", list;                // cpuvec_t holds 42 in 6.3.0; 32 is inside it
      bool nonbase = false;
      for (int i = 0; i < 32; i++) {
        Dl_info d;
        if (!dladdr(cpuvec[i], &d) || !d.dli_sname || std::strncmp(d.dli_sname, "__gmpn_", 7)) break;
        const char* c = cpu_suffix(d.dli_sname); if (!c) continue;
        if (!in(" x86_64 k8 fat ", c)) nonbase = true;
        if (!in(seen.c_str(), c)) { seen += std::string(c) + " "; list += (list.empty() ? "" : ", ") + std::string(c); }
      }
      decided = !list.empty(); generic = decided && !nonbase;
      std::printf("      mpn path: fat binary, dispatched to %s\n",
                  decided ? list.c_str() : "names this does not recognise");
    }
  }
#else
  from_march();
#endif

  std::printf("MPFR  %-8s  tune: %s   gmp-internals: %s   tls: %s\n\n", mpfr_get_version(),
              mpfr_buildopt_tune_case(), mpfr_buildopt_gmpinternals_p() ? "yes" : "no",
              mpfr_buildopt_tls_p() ? "yes" : "no");
  std::printf("%-28s %9s %9s %9s %9s\n", "ns/op", "mul", "add", "div", "sqrt");
  for (mpfr_prec_t p : {(mpfr_prec_t)256, (mpfr_prec_t)389, (mpfr_prec_t)1024}) {
    mpfr_t x, y, z; mpfr_inits2(p, x, y, z, (mpfr_ptr)0);
    mpfr_const_pi(x, MPFR_RNDN); mpfr_sqrt_ui(y, 2, MPFR_RNDN);   // every limb in use
    long it = 4000000 / (1 + p / 64);
    char tag[48]; std::snprintf(tag, sizeof tag, "%ld bits (%ld limbs)", (long)p, (long)((p + 63) / 64));
    std::printf("%-28s %9.2f %9.2f %9.2f %9.2f\n", tag,
                ns([&]{ mpfr_mul(z, x, y, MPFR_RNDN); }, it), ns([&]{ mpfr_add(z, x, y, MPFR_RNDN); }, it),
                ns([&]{ mpfr_div(z, x, y, MPFR_RNDN); }, it), ns([&]{ mpfr_sqrt(z, x, MPFR_RNDN); }, it));
    mpfr_clears(x, y, z, (mpfr_ptr)0);
  }
  mpfr_free_cache();

  if (generic) std::printf("\nSLOW: this GMP is on its baseline path, the k8 code.  tools/deps.sh"
                           " builds one that knows what it is running on.\n");
  else std::printf("\n%s\n", decided ? "OK: GMP is on a CPU-specific path."
                                     : "Undecided: build with tools/deps.sh and compare the numbers above.");
  return generic ? 1 : 0;
}
