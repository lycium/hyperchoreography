#!/bin/sh
# Every plate in plates.py to comp/plates/<name>.png (alpha kept), the backdrop as
# an exact-valued opaque plate, and the end card's turning orbit as a movie.
# The backdrop goes through a raw pipe because ffmpeg's color source lands one code
# value off through yuv.
set -e
cd "$(dirname "$0")/.."
W=${W:-1920}
H=${H:-1080}
RES=${RES:-1080p}
for s in TitleOrbit TitleName EndOrbit EndSubscribe EndCatalogue EndCode EndCredits; do
  .venv/bin/manim render -v WARNING -r "$W,$H" -s -t --media_dir comp/.media comp/plates.py "$s"
  cp "comp/.media/images/plates/${s}_ManimCE_v0.21.0.png" "comp/plates/$s.png"
done
rm -rf comp/.media
W=$W H=$H python3 -c "
import subprocess
import os
w, h = int(os.environ['W']), int(os.environ['H'])
raw = bytes([10, 12, 17]) * (w * h)
subprocess.run(['ffmpeg', '-v', 'error', '-y', '-f', 'rawvideo', '-pix_fmt', 'rgb24',
                '-s', '%dx%d' % (w, h), '-i', '-', '-frames:v', '1',
                'comp/plates/Backdrop.png'], input=raw, check=True)"
.venv/bin/python comp/make_plate_movie.py EndOrbitTurning plates/EndOrbit.mkv --res "$RES"

ls -la comp/plates/
