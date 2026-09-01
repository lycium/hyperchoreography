"""The action functional, in numpy, mirroring src/action.hpp.

Everything the presentation animates is computed here rather than mocked up, so a
frame that shows a gradient dropping is showing this module's gradient dropping.

Loop:      q(t) = sum_m  c_m cos(mt) + s_m sin(mt),   m in 1..K, m % N != 0
Action:    A[q] = 1/2 int |q'|^2 dt  +  1/2 sum_{k=1}^{N-1} int |q(t) - q(t+2pi k/N)|^-alpha dt

The kinetic half is closed form (pi m^2 per coefficient); the potential half is a
trapezoid sum over M nodes, which is spectrally accurate for a collision-free loop.
Coefficients are stored as x[mu, 0, :] = c_m and x[mu, 1, :] = s_m, matching the
record layout in src/catalog.hpp.
"""

from __future__ import annotations

import numpy as np

PI = np.pi


def mode_list(N: int, K: int) -> np.ndarray:
    """1..K with the multiples of N dropped: those only move the centre of mass."""
    m = np.arange(1, K + 1)
    return m[m % N != 0]


def auto_M(N: int, K: int) -> int:
    """Quadrature nodes: a multiple of N (so the shifts are exact) and of 8, >= 8K."""
    l = np.lcm(N, 8)
    target = max(64, 8 * K)
    return int(l * ((target + l - 1) // l))


class Action:
    """A (N, d, K) choreography problem: the action, its gradient and its Hessian."""

    def __init__(self, N: int, d: int, K: int, M: int | None = None,
                 alpha: float = 1.0, modes=None, minsep: float = 2e-3, omega=None):
        self.N, self.d, self.alpha, self.minsep = N, d, float(alpha), minsep
        self.modes = np.asarray(mode_list(N, K) if modes is None else modes, dtype=int)
        self.K = int(self.modes.max()) if len(self.modes) else 0
        self.nm = len(self.modes)
        self.n = self.nm * 2 * d                       # free parameters
        self.M = int(M) if M else auto_M(N, self.K)
        if self.M % N:
            raise ValueError("M must be a multiple of N")

        j = np.arange(self.M)
        t = 2 * PI * j / self.M
        # basis[mu, 0, j] = cos(m t_j), basis[mu, 1, j] = sin(m t_j)
        self.basis = np.empty((self.nm, 2, self.M))
        self.basis[:, 0, :] = np.cos(np.outer(self.modes, t))
        self.basis[:, 1, :] = np.sin(np.outer(self.modes, t))
        self.kin = PI * self.modes.astype(float) ** 2   # per coefficient pair

        # the shifts k -> k + 2pi k/N; k and N-k give the same integral, so only half are needed
        self.shifts = [(k * self.M // N, 1.0) for k in range(1, (N - 1) // 2 + 1)]
        if N % 2 == 0:
            self.shifts.append((self.M // 2, 0.5))
        self.wq = 2 * PI / self.M
        self.set_omega(omega)

    # -- the rotating frame -----------------------------------------------
    def set_omega(self, omega):
        """q_j(t) = exp(Omega t) q(t + 2 pi j/N).

        Only the kinetic term changes, from |q'|^2 to |q' + Omega q|^2, and it stays
        quadratic: per mode it becomes the 2d x 2d block

            pi [[ m^2 I - Omega^2 ,  -2 m Omega ],
                [   2 m Omega     , m^2 I - Omega^2 ]]

        in place of the scalar pi m^2. The potential is untouched, because a rotation
        preserves every distance between bodies.
        """
        self.Om = None if omega is None else np.asarray(omega, dtype=float).reshape(self.d, self.d)
        self.Om2 = None if self.Om is None else self.Om @ self.Om
        return self

    # -- shape helpers ----------------------------------------------------
    def zeros(self) -> np.ndarray:
        return np.zeros((self.nm, 2, self.d))

    def flat(self, x) -> np.ndarray:
        return np.asarray(x).reshape(-1)

    def shaped(self, x) -> np.ndarray:
        return np.asarray(x).reshape(self.nm, 2, self.d)

    # -- the curve --------------------------------------------------------
    def samples(self, x) -> np.ndarray:
        """q(t_j) on the quadrature grid, shape (M, d)."""
        x = self.shaped(x)
        return np.einsum("mcj,mca->ja", self.basis, x)

    def curve(self, x, ts) -> np.ndarray:
        """q(t) at arbitrary times, shape (len(ts), d)."""
        x = self.shaped(x)
        ts = np.atleast_1d(np.asarray(ts, dtype=float))
        ph = np.outer(ts, self.modes)
        B = np.stack([np.cos(ph), np.sin(ph)], axis=2)      # (T, nm, 2)
        return np.einsum("tmc,mca->ta", B, x)

    def bodies(self, x, ts) -> np.ndarray:
        """All N bodies in the inertial frame. Shape (N, T, d)."""
        ts = np.atleast_1d(np.asarray(ts, dtype=float))
        out = np.stack([self.curve(x, ts + 2 * PI * k / self.N) for k in range(self.N)])
        if self.Om is not None:
            from .catalog import expm_skew
            R = np.stack([expm_skew(self.Om * t) for t in ts])
            out = np.einsum("tab,ktb->kta", R, out)
        return out

    def velocity(self, x, ts) -> np.ndarray:
        x = self.shaped(x)
        xd = np.empty_like(x)
        m = self.modes[:, None].astype(float)
        xd[:, 0, :] = m * x[:, 1, :]
        xd[:, 1, :] = -m * x[:, 0, :]
        return self.curve(xd, ts)

    # -- action and gradient ---------------------------------------------
    def value(self, x, grad: bool = False):
        """A[q], and optionally dA/dx with the same shape as x.

        Returns inf (and a zero gradient) on a colliding loop, exactly as the
        C++ kernel does: the line search then simply rejects the step.
        """
        x = self.shaped(x)
        Q = self.samples(x)                                   # (M, d)
        kinE = 0.5 * float(np.sum(self.kin[:, None, None] * x * x))
        kw = self._kin_omega(x)
        if kw is not None:
            kinE += 0.5 * float(np.sum(x * kw))

        U = 0.0
        G = np.zeros_like(Q) if grad else None
        for s, wgt in self.shifts:
            D = Q - np.roll(Q, -s, axis=0)                    # q(t_j) - q(t_j + shift)
            r2 = np.einsum("ja,ja->j", D, D)
            if r2.min() < self.minsep ** 2:
                return (np.inf, np.zeros_like(x)) if grad else np.inf
            if self.alpha == 1.0:
                inv = 1.0 / np.sqrt(r2)
                U += wgt * float(inv.sum())
                f = wgt * inv ** 3
            else:
                pa = r2 ** (-0.5 * self.alpha)
                U += wgt * float(pa.sum())
                f = wgt * self.alpha * pa / r2
            if grad:
                fD = f[:, None] * D
                G -= fD
                G += np.roll(fD, s, axis=0)

        A = kinE + self.wq * U
        if not grad:
            return A
        g = self.kin[:, None, None] * x + self.wq * np.einsum("mcj,ja->mca", self.basis, G)
        if kw is not None:
            g = g + kw
        return A, g

    def _kin_omega(self, x):
        """The frame's contribution to d(kinetic)/dx, or None with no frame."""
        if self.Om is None:
            return None
        m = self.modes[:, None].astype(float)
        c, s_ = x[:, 0, :], x[:, 1, :]
        out = np.empty_like(x)
        out[:, 0, :] = -c @ self.Om2.T - 2 * m * (s_ @ self.Om.T)
        out[:, 1, :] = -s_ @ self.Om2.T + 2 * m * (c @ self.Om.T)
        return PI * out

    def gradient(self, x) -> np.ndarray:
        return self.value(x, grad=True)[1]

    def hessian(self, x) -> np.ndarray:
        """The exact Hessian, (n, n), in the flattened coefficient order.

        d2A = kinetic diagonal  +  sum_shifts  B^T [ fA (D D^T) - f I ] B, assembled
        directly rather than by finite differences.
        """
        x = self.shaped(x)
        n, d, M = self.n, self.d, self.M
        Q = self.samples(x)
        H = np.zeros((n, n))

        # basis matrix E[j, i]: value of basis function i at node j, i indexing (mu, c)
        E = self.basis.reshape(self.nm * 2, M).T              # (M, 2*nm)
        nb = self.nm * 2

        for s, wgt in self.shifts:
            D = Q - np.roll(Q, -s, axis=0)
            r2 = np.einsum("ja,ja->j", D, D)
            if self.alpha == 1.0:
                inv = 1.0 / np.sqrt(r2)
                i3 = inv ** 3
                f = wgt * i3
                fA = 3.0 * wgt * i3 * inv * inv
            else:
                al = self.alpha
                pa = r2 ** (-0.5 * al)
                gk = al * pa / r2
                f = wgt * gk
                fA = wgt * (al + 2.0) * gk / r2
            # difference operator on the basis: Ds[j, i] = E[j, i] - E[j+s, i]
            Ds = E - np.roll(E, -s, axis=0)
            # W[j, a, b] = fA[j] D[j,a] D[j,b] - f[j] delta_ab
            for a in range(d):
                for b in range(a, d):
                    w = fA * D[:, a] * D[:, b] - (f if a == b else 0.0)
                    blk = self.wq * (Ds.T * w) @ Ds          # (nb, nb)
                    ia = np.arange(nb) * d + a
                    ib = np.arange(nb) * d + b
                    H[np.ix_(ia, ib)] += blk
                    if a != b:
                        H[np.ix_(ib, ia)] += blk.T
        # kinetic part: pi m^2 on the diagonal, plus the frame's block if there is one
        kd = np.repeat(np.repeat(self.kin, 2), d)
        H[np.diag_indices(n)] += kd
        if self.Om is not None:
            for mu, mm in enumerate(self.modes.astype(float)):
                i0 = (2 * mu) * d
                i1 = (2 * mu + 1) * d
                H[i0:i0 + d, i0:i0 + d] += -PI * self.Om2
                H[i1:i1 + d, i1:i1 + d] += -PI * self.Om2
                H[i0:i0 + d, i1:i1 + d] += -2 * PI * mm * self.Om
                H[i1:i1 + d, i0:i0 + d] += 2 * PI * mm * self.Om
        return 0.5 * (H + H.T)

    # -- scaling ----------------------------------------------------------
    def optimal_scale(self, x) -> float:
        """Kepler rescaling q -> lam q that makes the action stationary in size."""
        A, _ = self.value(x, grad=True)
        if not np.isfinite(A):
            return 1.0
        x = self.shaped(x)
        kinE = 0.5 * float(np.sum(self.kin[:, None, None] * x * x))
        potE = A - kinE
        if kinE <= 0 or potE <= 0:
            return 1.0
        return float((self.alpha * potE / (2.0 * kinE)) ** (1.0 / (self.alpha + 2.0)))

    # -- diagnostics ------------------------------------------------------
    def min_separation(self, x) -> float:
        Q = self.samples(x)
        best = np.inf
        for k in range(1, self.N):
            s = k * self.M // self.N
            D = Q - np.roll(Q, -s, axis=0)
            best = min(best, float(np.sqrt(np.einsum("ja,ja->j", D, D)).min()))
        return best

    def energy(self, x) -> float:
        """Total energy per body, from the virial theorem E = -N A / (6 pi) at alpha = 1."""
        A = self.value(x)
        return -self.N * A / (6 * PI)

    def effective_dimension(self, x, rel_tol: float = 1e-8):
        """Principal values of X^T X and the count above tolerance -- the record's deff."""
        x = self.shaped(x).reshape(-1, self.d)
        w = np.linalg.eigvalsh(x.T @ x)[::-1]
        sv = np.sqrt(np.maximum(w, 0.0))
        return int(np.sum(w > rel_tol * w[0])), sv
