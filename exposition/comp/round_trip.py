#!/usr/bin/env python3
"""RGB in, RGB back out of the upload: it has to be exact. Called by verify_pipe.sh."""

import sys

import numpy as np


def main() -> int:
    a = np.fromfile(sys.argv[1], np.uint8).astype(np.int16)
    b = np.fromfile(sys.argv[2], np.uint8).astype(np.int16)
    if a.size != b.size or a.size == 0:
        print("  FAIL the two dumps are %d and %d samples" % (a.size, b.size))
        return 1
    d = np.abs(a - b)
    if d.max() == 0:
        print("  ok   RGB round trip: bit-exact over %d samples" % a.size)
        return 0
    print("  FAIL RGB round trip: off by %d code values, on %.4f %% of %d samples"
          % (d.max(), 100.0 * (d > 0).mean(), a.size))
    return 1


if __name__ == "__main__":
    sys.exit(main())
