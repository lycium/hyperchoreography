"""Throwaway: one film line in each British-male Kokoro voice, for choosing."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

LINE = ("A choreography is the boldest possible guess about the answer: "
        "that all N bodies trace out one and the same curve.")

os.makedirs("comp/audition", exist_ok=True)
for v in ("bm_daniel", "bm_fable", "bm_george", "bm_lewis"):
    os.environ["EXPO_VOICE"] = "kokoro:" + v
    # a fresh import per voice, since narrate binds VOICE at import time
    for m in [m for m in list(sys.modules) if m.startswith("expo")]:
        del sys.modules[m]
    from expo import narrate
    p = narrate.clip(LINE)
    if p:
        import shutil
        shutil.copy(p, f"comp/audition/{v}.wav")
        print(v, "ok", f"{narrate.duration(LINE):.2f}s")
    else:
        print(v, "FAILED")
