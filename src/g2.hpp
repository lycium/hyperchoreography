// The d = 7 Fano structure: the maximal torus of g2 ⊂ so(7) and the Fano-twisted symmetry class. The
// associative 3-form itself lives in calib.hpp, which is where the cross product's content actually is.
#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include "linalg.hpp"

// Ω in the maximal torus of g2: rates (w1, w2, w1+w2) on the Fourier planes of the Fano 7-cycle
// (u_k[j] = √(2/7) cos 2πkj/7, v_k[j] = √(2/7) sin 2πkj/7, k = 1,2,3; the axis (1,…,1)/√7 is fixed).
// The relation w3 = w1 + w2 is exactly the condition that exp(Ωt) preserves the associative 3-form.
// For d > 7 the same φ lives on span(e_0…e_6) and stab(φ) grows to g2 ⊕ so(d−7), so `rest` carries free
// rates on the leftover planes (7,8), (9,10), …. su_omega is not an alternative: its planes are not Fano.
inline std::vector<double> g2_omega(double w1, double w2, int d = 7, const std::vector<double>& rest = {}) {
  if (d < 7) throw std::runtime_error("g2 frame needs d >= 7");
  if ((int)rest.size() > (d - 7) / 2) throw std::runtime_error("g2 frame: more extra rates than leftover planes");
  const double PI7 = 6.283185307179586 / 7, w[4] = {0, w1, w2, w1 + w2};
  std::vector<double> Om((size_t)d * d, 0.0), u(7), v(7);
  for (int k = 1; k <= 3; k++) {
    for (int j = 0; j < 7; j++) { u[j] = std::sqrt(2.0 / 7) * std::cos(PI7 * k * j); v[j] = std::sqrt(2.0 / 7) * std::sin(PI7 * k * j); }
    for (int i = 0; i < 7; i++) for (int j = 0; j < 7; j++) Om[(size_t)i * d + j] += w[k] * (v[i] * u[j] - u[i] * v[j]);
  }
  for (size_t p = 0; p < rest.size(); p++) {                 // these planes act trivially on φ
    const int a = 7 + 2 * (int)p, b = a + 1;
    Om[(size_t)a * d + b] = -rest[p]; Om[(size_t)b * d + a] = rest[p];
  }
  return Om;
}
// q(t + 2πp/7) = σ q(t) with σ the Fano 7-cycle. Mode m then lives in the σ-eigenplane k = mp mod 7,
// so modes ≡ 0 (mod 7) carry the fixed axis: reaching deff = 7 needs one of those, impossible at N = 7.
inline std::string fano_sym(int p) { return "t+" + std::to_string(((p % 7) + 7) % 7) + "/7 s[2,3,4,5,6,7,1]"; }
inline int fano_plane(int m, int p) { int k = (((m * p) % 7) + 7) % 7; return k > 3 ? 7 - k : k; }   // 0 = the fixed axis
// σ-eigenplane k (1..3): u_k[j] = √(2/7) cos 2πkj/7, v_k[j] = √(2/7) sin 2πkj/7; k = 0 gives the fixed axis
inline void fano_basis(int k, double* u, double* v) {
  const double c = std::sqrt(2.0 / 7), PI7 = 6.283185307179586 / 7;
  if (k == 0) { for (int j = 0; j < 7; j++) { u[j] = std::sqrt(1.0 / 7); v[j] = 0; } return; }
  for (int j = 0; j < 7; j++) { u[j] = c * std::cos(PI7 * k * j); v[j] = c * std::sin(PI7 * k * j); }
}
