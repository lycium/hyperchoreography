"""Reading catalogued orbits.

`./hyperchoreography show <file> --id <i>` prints one record as JSON; that is the
only interface used here, so a record drawn in the presentation is byte-for-byte
the record in the catalogue. make_data.py caches the JSON under data/ so the
scenes render without the binary present.
"""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # the repository root
DATA = os.path.join(os.path.dirname(HERE), "data")
BINARY = os.path.join(ROOT, "hyperchoreography")
CATALOG = os.path.join(ROOT, "catalog")


@dataclass
class Orbit:
    """One catalogued choreography: the loop, its scalars and its frame."""
    name: str
    id: int
    N: int
    d: int
    deff: int
    K: int
    action: float
    energy: float
    morse: int
    nullity: int
    minsep: float
    ret_err: float
    twist: float
    modes: np.ndarray                # (nm,) integers
    coef: np.ndarray                 # (nm, 2, d): [:, 0] = cos, [:, 1] = sin
    omega: np.ndarray | None         # (d, d) skew, or None when inertial
    pca: np.ndarray
    sym: str

    # -- geometry ---------------------------------------------------------
    def loop(self, ts) -> np.ndarray:
        """The single curve q(t) in the rotating frame, shape (len(ts), d)."""
        ts = np.atleast_1d(np.asarray(ts, dtype=float))
        ph = np.outer(ts, self.modes)
        B = np.stack([np.cos(ph), np.sin(ph)], axis=2)
        return np.einsum("tmc,mca->ta", B, self.coef)

    def bodies(self, ts) -> np.ndarray:
        """All N bodies in the inertial frame: q_k(t) = exp(Omega t) q(t + 2 pi k/N).

        Shape (N, len(ts), d). With no frame this is just the shifted loop.
        """
        ts = np.atleast_1d(np.asarray(ts, dtype=float))
        out = np.stack([self.loop(ts + 2 * np.pi * k / self.N) for k in range(self.N)])
        if self.omega is not None:
            R = np.stack([expm_skew(self.omega * t) for t in ts])       # (T, d, d)
            out = np.einsum("tab,ktb->kta", R, out)
        return out

    def frame_bodies(self, ts) -> np.ndarray:
        """The bodies as seen from the rotating frame: no exp(Omega t) applied.

        In that frame the motion is just the one loop with N time shifts, which is
        usually a simple figure; the complication people see is the frame turning.
        """
        ts = np.atleast_1d(np.asarray(ts, dtype=float))
        return np.stack([self.loop(ts + 2 * np.pi * k / self.N) for k in range(self.N)])

    def principal_frame(self) -> np.ndarray:
        """Rows are the principal axes of the motion, most-occupied first.

        A deff = 9 orbit and a circle look identical in whichever two coordinates
        you happen to pick, so every picture is drawn in this frame.
        """
        ts = np.linspace(0, 2 * np.pi, 512, endpoint=False)
        X = self.bodies(ts).reshape(-1, self.d)
        w, V = np.linalg.eigh(X.T @ X)
        V = V[:, ::-1]
        for c in range(V.shape[1]):                        # deterministic signs
            k = np.argmax(np.abs(V[:, c]))
            if V[k, c] < 0:
                V[:, c] = -V[:, c]
        return V.T

    def principal_values(self) -> np.ndarray:
        ts = np.linspace(0, 2 * np.pi, 512, endpoint=False)
        X = self.bodies(ts).reshape(-1, self.d)
        w = np.linalg.eigvalsh(X.T @ X)[::-1]
        return np.sqrt(np.maximum(w, 0.0))

    def mode_power(self) -> np.ndarray:
        return np.einsum("mca,mca->m", self.coef, self.coef)


def expm_skew(A: np.ndarray) -> np.ndarray:
    """exp of a skew matrix by its eigenvalues -- always orthogonal to round-off."""
    w, V = np.linalg.eigh(1j * A)
    return np.real(V @ np.diag(np.exp(-1j * w)) @ V.conj().T)


def _from_json(js: dict, name: str) -> Orbit:
    d, nm = js["d"], len(js["modes"])
    coef = np.asarray(js["coef"], dtype=float).reshape(nm, 2, d)
    om = js.get("omega")
    return Orbit(
        name=name, id=js["id"], N=js["N"], d=d, deff=js["deff"], K=js["K"],
        action=js["action"], energy=js["energy"], morse=js.get("morse", -1),
        nullity=js.get("nullity", -1), minsep=js.get("minsep", 0.0),
        ret_err=js.get("ret_err", -1.0), twist=js.get("twist", 0.0),
        modes=np.asarray(js["modes"], dtype=int), coef=coef,
        omega=np.asarray(om, dtype=float).reshape(d, d) if om else None,
        pca=np.asarray(js.get("pca", []), dtype=float), sym=js.get("sym", ""))


def load(name: str) -> Orbit:
    """A cached orbit from data/<name>.json."""
    path = os.path.join(DATA, name + ".json")
    with open(path) as f:
        return _from_json(json.load(f), name)


def fetch(catalog_file: str, rid: int, name: str) -> Orbit:
    """Pull a record straight out of a catalogue with the solver binary."""
    out = subprocess.run([BINARY, "show", os.path.join(CATALOG, catalog_file),
                          "--id", str(rid)], capture_output=True, text=True, check=True)
    js = json.loads(out.stdout.splitlines()[0])
    if js["id"] != rid:
        raise RuntimeError("show returned record %d, asked for %d" % (js["id"], rid))
    return _from_json(js, name)
