"""Synthesise every narration line into the cache, in parallel, before a render."""
import ast
import glob
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

os.environ.setdefault("EXPO_VOICE", "kokoro:bm_daniel")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from expo import narrate

TAG = re.compile(r"<[^>]+>")


def spoken_lines(path):
    tree = ast.parse(open(path).read())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr not in ("say", "say_with") or not node.args:
            continue
        first = node.args[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            yield TAG.sub("", first.value)


texts, seen = [], set()
for f in sorted(glob.glob(os.path.join(ROOT, "scenes", "s*.py"))):
    for t in spoken_lines(f):
        if t not in seen:
            seen.add(t)
            texts.append(t)
todo = [t for t in texts if not narrate.measured(t)]
print("%d literal lines, %d not yet synthesised" % (len(texts), len(todo)))


def one(t):
    return t if narrate.duration(t) is None else None


with ThreadPoolExecutor(max_workers=6) as ex:
    failed = [t for t in ex.map(one, todo) if t]
for t in failed:
    print("FAILED:", t[:60])
narrate._save_index()
print("done,", len(failed), "failures")
