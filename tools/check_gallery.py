#!/usr/bin/env python3
"""Check the generated gallery against its catalogues, without executing browser code.

Python standard library only. If Node is installed, also parse every inline script.
This checks packaging and metadata consistency, not the mathematics or visual rendering.
"""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

from gallery import cat_list, cat_show, fingerprint


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Scripts(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.scripts = []
        self.current = None

    def handle_starttag(self, tag, attrs):
        if tag == "script":
            self.current = [dict(attrs), ""]

    def handle_data(self, text):
        if self.current is not None:
            self.current[1] += text

    def handle_endtag(self, tag):
        if tag == "script" and self.current is not None:
            self.scripts.append(self.current)
            self.current = None


def finite_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite gallery metadata")
    elif isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


def validate_payload(data):
    finite_tree(data)
    records = data["recs"]
    paths = data["meta"]["paths"]
    require(bool(records), "empty gallery")
    keys = set()
    for r in records:
        key = (r["f"], r["i"])
        require(key not in keys, "duplicate record key: %s" % (key,))
        keys.add(key)
        require(r["f"] in paths, "missing catalogue path: %s" % r["f"])
        for field in ("N", "d", "K", "M", "S"):
            require(type(r[field]) is int and r[field] > 0, "invalid %s: %s" % (field, key))
        require(r["N"] >= 2 and r["S"] % r["N"] == 0, "invalid cyclic sampling: %s" % (key,))
        require(0 <= r["view_rank"] <= r["d"], "invalid display rank: %s" % (key,))
        require(r["proof"] in ("none", "legacy", "revision2"), "unknown proof status: %s" % (key,))
        require(len(r["sc"]) == r["d"], "invalid curve scale: %s" % (key,))
        require("G" not in r or len(r["G"]) == r["d"] ** 2, "invalid reconstruction matrix: %s" % (key,))
        curve = base64.b64decode(r["q"], validate=True)
        expected = 2 * r["S"] * r["d"] * (r["N"] if r.get("all") else 1)
        require(len(curve) == expected, "truncated packed curve: %s" % (key,))
        for source in r["sources"]:
            require(re.fullmatch(r"\d{4}\.\d{4,5}", source["arxiv"]), "invalid source identifier")
            require(source["version"] >= 1 and source["orbit_id"] >= 1, "invalid source version/orbit")
    return keys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", nargs="?", default="docs/index.html")
    parser.add_argument("--binary", default="./hyperchoreography")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    html = (root / args.html).read_text(encoding="utf-8")
    require(not re.search(r"__(?:DATA|JS|CSS|CARDS)__", html), "unexpanded gallery placeholder")
    scripts = Scripts()
    scripts.feed(html)
    payloads = [body for attrs, body in scripts.scripts if attrs.get("id") == "DATA"]
    require(len(payloads) == 1, "expected exactly one embedded DATA payload")
    data = json.loads(payloads[0])
    keys = validate_payload(data)
    paths = {label: root / path for label, path in data["meta"]["paths"].items()}
    before = {label: fingerprint(path) for label, path in paths.items()}
    binary = str((root / args.binary).resolve())
    expected = {(label, row["id"]) for label, path in paths.items() for row in cat_list(binary, str(path))}
    require(keys == expected, "gallery records differ from the catalogues; regenerate with make gallery")

    def check_record(r):
        actual = cat_show(binary, str(paths[r["f"]]), r["i"])
        key = "%s#%s" % (r["f"], r["i"])
        for field, stored in (("N", "N"), ("d", "d"), ("de", "deff"), ("K", "K"), ("M", "M")):
            require(r[field] == actual[stored], "stale %s: %s" % (field, key))
        require(r["sources"] == actual["sources"], "lost or stale attribution: " + key)
        proof = "revision2" if actual.get("proven", 0) > 0 else "legacy" if actual.get("legacy_proof", 0) > 0 else "none"
        require(r["proof"] == proof, "stale proof status: " + key)

    with ThreadPoolExecutor(max_workers=4) as pool:
        for _ in pool.map(check_record, data["recs"]):
            pass
    require(before == {label: fingerprint(path) for label, path in paths.items()}, "catalogue changed during validation")
    node = shutil.which("node")
    if node:
        bodies = [body for attrs, body in scripts.scripts if attrs.get("type", "") != "application/json"]
        subprocess.run([node, "-e", "const fs=require('fs'),vm=require('vm');"
                        "for(const body of JSON.parse(fs.readFileSync(0,'utf8')))new vm.Script(body);"],
                       input=json.dumps(bodies), text=True, check=True, timeout=30)
    attributed = sum(bool(r["sources"]) for r in data["recs"])
    print("gallery check: %d records, %d attributed; catalogue metadata and packed curves agree" % (len(keys), attributed))
    print("inline JavaScript: parsed" if node else "inline JavaScript: not checked (Node not installed)")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError, subprocess.SubprocessError) as error:
        print("gallery check failed: %s" % error, file=sys.stderr)
        sys.exit(1)
