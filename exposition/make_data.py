#!/usr/bin/env python3
"""Cache the orbits the presentation draws.

Runs `./hyperchoreography show` once per record and writes data/<name>.json, so
rendering needs neither the solver binary nor the catalogue. Re-run after a
harvest renumbers ids.
"""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BINARY = os.path.join(ROOT, "hyperchoreography")
DATA = os.path.join(HERE, "data")

# name -> (catalogue file, record id, what it is)
WANTED = {
    "eight":        ("d2-3_n3.bin",       0,  "the figure eight, N = 3"),
    "n3_spatial":   ("d2-3_n3.bin",       4,  "N = 3 out of the plane, deff = 3"),
    "n4_planar":    ("d2-4_n4.bin",       0,  "the planar N = 4 minimiser"),
    "hiphop":       ("d2-4_n4.bin",       9,  "the N = 4 hip-hop resonance, modes (5,6)"),
    "n4_deff4":     ("d2-4_n4.bin",      15,  "N = 4 filling all four dimensions"),
    "n6_deff6":     ("d5-6_n6.bin",     139,  "N = 6 at deff = 6, frame su:1,2"),
    "d7_inertial":  ("d7_n10.bin",        0,  "deff = 7 with no rotating frame"),
    "d7_twist":     ("d7_n10_g2_16.bin", 13,  "the G2 twist champion, frame g2:1,6"),
    "d11":          ("d11_n12_g2.bin",   13,  "deff = 11, frame g2:1,2,4,5"),
}


def main() -> int:
    if not os.path.exists(BINARY):
        print("no solver binary at %s -- run `make` in the repository root" % BINARY)
        return 1
    os.makedirs(DATA, exist_ok=True)
    index = {}
    for name, (cat, rid, what) in WANTED.items():
        path = os.path.join(ROOT, "catalog", cat)
        if not os.path.exists(path):
            print("  skip %-12s (no %s)" % (name, cat))
            continue
        out = subprocess.run([BINARY, "show", path, "--id", str(rid)],
                             capture_output=True, text=True)
        if out.returncode != 0 or not out.stdout.strip():
            print("  fail %-12s %s" % (name, out.stderr.strip()[:120]))
            continue
        js = json.loads(out.stdout.splitlines()[0])
        if js["id"] != rid:                       # ids are positions; a rewrite renumbers them
            print("  fail %-12s: asked for id %d, got %d" % (name, rid, js["id"]))
            continue
        with open(os.path.join(DATA, name + ".json"), "w") as f:
            json.dump(js, f)
        index[name] = {"catalog": cat, "id": rid, "what": what, "N": js["N"],
                       "d": js["d"], "deff": js["deff"], "action": js["action"]}
        print("  %-12s %-18s #%-4d N=%-3d d=%d deff=%-2d A=%.9f"
              % (name, cat, rid, js["N"], js["d"], js["deff"], js["action"]))
    with open(os.path.join(DATA, "index.json"), "w") as f:
        json.dump(index, f, indent=1)
    print("%d records -> %s" % (len(index), DATA))
    return 0


if __name__ == "__main__":
    sys.exit(main())
