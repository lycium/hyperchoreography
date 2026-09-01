"""Throwaway: synthesise every narration line into the cache, in parallel."""
import glob, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
os.environ["EXPO_VOICE"] = "kokoro:bm_daniel"
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from expo import narrate

texts = []
seen = set()
for f in sorted(glob.glob("out/*.cues.json")):
    for c in json.load(open(f))["cues"]:
        if c["text"] not in seen:
            seen.add(c["text"])
            texts.append(c["text"])
print(len(texts), "unique lines")

def one(t):
    return t if narrate.clip(t) is None else None

with ThreadPoolExecutor(max_workers=6) as ex:
    failed = [t for t in ex.map(one, texts) if t]
for t in failed:
    print("FAILED:", t[:60])
print("done,", len(failed), "failures")
