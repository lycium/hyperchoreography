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


# name -> (catalogue file, record id, digits): the existence proofs the film shows, with their step traces
PROOFS = {
    "prove_eight": ("d2-3_n3.bin", 0, 40),
    "prove_g2":    ("d7_n10_g2.bin", 0, 40),
    "prove_d11":   ("d11_n12_g2.bin", 13, 40),
}


def prove(cat, rid, digits):
    """Run `prove --verbose --write` and keep everything it measured: the certificate, the step traces
    of the point run and of the (last) box run, and the rough enclosure of body 0 at every step. The
    proof is recorded in the catalogue as well, as `prove --write` would."""
    import re
    path = os.path.join(ROOT, "catalog", cat)
    r = subprocess.run([BINARY, "prove", path, "--id", str(rid), "--digits", str(digits), "--verbose",
                        "--force", "--write", "--threads", str(os.cpu_count() or 4)], capture_output=True, text=True)
    log = r.stdout
    if "PROVEN in" not in log:
        return None
    out = {"catalog": cat, "id": rid, "digits": digits, "point": [], "box": [], "wbox": [], "frame": None}
    section, box, wbox = "point", [], []
    for line in log.splitlines():
        m = re.search(r"step \d+: h=([0-9.e+-]+) t=([0-9.]+) width ([0-9.e+-]+)", line)
        if m:
            (out["point"] if section == "point" else box).append([float(m.group(2)), float(m.group(1)), float(m.group(3))])
            continue
        m = re.match(r"\s*wbox ([0-9.e+-]+) ([0-9.e+-]+) ([0-9.e+-]+) ([0-9.e+-]+) ([0-9.e+-]+) ([0-9.e+-]+)", line)
        if m:
            if section != "point":
                wbox.append([float(m.group(i)) for i in range(1, 7)])
            continue
        if line.strip().startswith("point run:"):
            section = "box"
        elif line.strip().startswith("box run:"):
            out["box"], out["wbox"] = box, wbox
            box, wbox = [], []
            out["maxwid"] = float(re.search(r"max state width ([0-9.e+-]+)", line).group(1))
        m = re.search(r"frame: (\d+) planes, rates ([^;]+); (\d+) translations, (\d+) commuting rotations", line)
        if m:
            out["frame"] = {"planes": int(m.group(1)), "rates": m.group(2).split(),
                            "trans": int(m.group(3)), "rots": int(m.group(4))}
    g = lambda pat, conv=float: conv(re.search(pat, log).group(1))
    out.update({
        "N": g(r"choreography of N=(\d+) bodies", int), "d": g(r"bodies in R\^(\d+)", int),
        "seconds": g(r"PROVEN in ([0-9.]+)s"), "radius": g(r"within ([0-9.e+-]+) \(max norm\)"),
        "newton": g(r"\|Y F\| = ([0-9.e+-]+)"), "kappa": g(r"contraction ([0-9.e+-]+)"), "closure": g(r"closure ([0-9.e+-]+)"),
        "slice": g(r"slice dim (\d+)", int), "gauge": g(r"(\d+) gauge generators", int),
        "steps": g(r"(\d+) validated Taylor steps", int),
        "energy": list(re.search(r"energy in \[([^,]+), ([^\]]+)\]", log).groups()),
        "action": list(re.search(r"action in \[([^,]+), ([^\]]+)\]", log).groups()),
    })
    return out


def main() -> int:
    """With no arguments everything is refreshed; names (`eight`, `prove_d11`, ...) select."""
    if not os.path.exists(BINARY):
        print("no solver binary at %s -- run `make` in the repository root" % BINARY)
        return 1
    only = set(sys.argv[1:])
    os.makedirs(DATA, exist_ok=True)
    index = {}
    if os.path.exists(os.path.join(DATA, "index.json")):
        with open(os.path.join(DATA, "index.json")) as f:
            index = json.load(f)
    for name, (cat, rid, what) in WANTED.items():
        if only and name not in only:
            continue
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
    for name, (cat, rid, digits) in PROOFS.items():
        if only and name not in only:
            continue
        if not os.path.exists(os.path.join(ROOT, "catalog", cat)):
            print("  skip %-12s (no %s)" % (name, cat))
            continue
        pr = prove(cat, rid, digits)
        if pr is None:
            print("  fail %-12s: not proven" % name)
            continue
        with open(os.path.join(DATA, name + ".json"), "w") as f:
            json.dump(pr, f)
        print("  %-12s %-18s #%-4d proven in %.1f s: radius %.0e, contraction %.1e, %d steps"
              % (name, cat, rid, pr["seconds"], pr["radius"], pr["kappa"], pr["steps"]))
    print("%d records -> %s" % (len(index), DATA))
    return 0


if __name__ == "__main__":
    sys.exit(main())
