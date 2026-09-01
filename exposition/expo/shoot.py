"""Integrating the real equations, and the shooting Newton that certifies a record.

A critical point of the truncated Fourier problem is not yet a solution of the
N-body problem: it is a curve whose action is stationary among curves with that
many modes. What settles the question is leaving the Fourier world entirely --
read the positions and velocities off the series, integrate Newton's equations,
and ask whether the state really does come back to itself shifted by one place
after a period over N.

The solver uses a 22nd-order Taylor method in C++; this is a plain fourth-order
Runge-Kutta with a small step, which reaches about 1e-12 and is enough to show what
the Newton does. Results are cached, since a Jacobian costs 2*N*d integrations.
"""

from __future__ import annotations

import os

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def accelerations(pos: np.ndarray) -> np.ndarray:
    """Newtonian acceleration of N unit masses at `pos`, shape (N, d)."""
    diff = pos[None, :, :] - pos[:, None, :]                # (N, N, d)
    r2 = np.einsum("ija,ija->ij", diff, diff)
    np.fill_diagonal(r2, np.inf)
    inv3 = r2 ** -1.5
    return np.einsum("ij,ija->ia", inv3, diff)


def rk4(pos, vel, T: float, steps: int):
    """Fixed-step classical Runge-Kutta on (q, v)."""
    h = T / steps
    p, v = np.array(pos, dtype=float), np.array(vel, dtype=float)
    for _ in range(steps):
        k1p, k1v = v, accelerations(p)
        k2p, k2v = v + 0.5 * h * k1v, accelerations(p + 0.5 * h * k1p)
        k3p, k3v = v + 0.5 * h * k2v, accelerations(p + 0.5 * h * k2p)
        k4p, k4v = v + h * k3v, accelerations(p + h * k3p)
        p = p + (h / 6.0) * (k1p + 2 * k2p + 2 * k3p + k4p)
        v = v + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
    return p, v


def energy(pos, vel) -> float:
    ke = 0.5 * float(np.einsum("ia,ia->", vel, vel))
    d = pos[None, :, :] - pos[:, None, :]
    r = np.sqrt(np.einsum("ija,ija->ij", d, d))
    iu = np.triu_indices(len(pos), 1)
    return ke - float(np.sum(1.0 / r[iu]))


def state_from_loop(P, x) -> np.ndarray:
    """The 2*N*d numbers the record is really about: every body's place and speed."""
    ts = np.array([2 * np.pi * k / P.N for k in range(P.N)])
    pos = P.curve(x, ts)
    vel = P.velocity(x, ts)
    return np.concatenate([pos.reshape(-1), vel.reshape(-1)])


def shift_residual(Z, N: int, d: int, T: float, steps: int, detail=False):
    """How far the state is from coming back to itself with the bodies relabelled.

    After a period over N, body k should be exactly where body k+1 started. That is
    the whole choreography condition, and it is one integration to test.
    """
    nd = N * d
    p0 = Z[:nd].reshape(N, d)
    v0 = Z[nd:].reshape(N, d)
    p1, v1 = rk4(p0, v0, T / N, steps)
    want_p = np.roll(p0, -1, axis=0)
    want_v = np.roll(v0, -1, axis=0)
    F = np.concatenate([(p1 - want_p).reshape(-1), (v1 - want_v).reshape(-1)])
    if detail:
        return F, (p1, v1)
    return F


def shoot_newton(Z, N: int, d: int, T: float = 2 * np.pi, steps: int = 6000,
                 iters: int = 6, damping: float = 1e-12):
    """Correct the state until the shift condition holds.

    A Gauss-Newton on 2*N*d equations in 2*N*d unknowns, with the Jacobian taken by
    central differences -- one integration per column, which is why this is cached.
    """
    Z = np.array(Z, dtype=float)
    n = Z.size
    scale = max(1.0, float(np.abs(Z).max()))
    hist = []
    for _ in range(iters + 1):
        F = shift_residual(Z, N, d, T, steps)
        hist.append((Z.copy(), float(np.abs(F).max()) / scale))
        if len(hist) > iters or hist[-1][1] < 1e-13:
            break
        J = np.empty((n, n))
        h = 1e-6
        for c in range(n):
            Zp, Zm = Z.copy(), Z.copy()
            Zp[c] += h
            Zm[c] -= h
            J[:, c] = (shift_residual(Zp, N, d, T, steps)
                       - shift_residual(Zm, N, d, T, steps)) / (2 * h)
        JtJ = J.T @ J + damping * np.trace(J.T @ J) / n * np.eye(n)
        Z = Z - np.linalg.solve(JtJ, J.T @ F)
    return Z, hist


def cached_run(name: str, build):
    path = os.path.join(DATA, name + ".npz")
    if os.path.exists(path):
        with np.load(path) as f:
            return {k: f[k] for k in f.files}
    out = build()
    os.makedirs(DATA, exist_ok=True)
    np.savez_compressed(path, **out)
    return out
