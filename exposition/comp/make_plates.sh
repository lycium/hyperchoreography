#!/bin/sh
# Render every plate in plates.py to comp/plates/<name>.png — 1920x1080, alpha kept —
# plus the backdrop, which is the film's own background colour as an OPAQUE plate.
# It is written through a raw rgb24 pipe on purpose: ffmpeg's lavfi color source
# renders via yuv420p and lands one code value off (10,11,16 for 0x0A0C11), a step
# the master's dips to black would show against the scenes' exact (10,12,17).
set -e
cd "$(dirname "$0")/.."
for s in TitleOrbit TitleName EndOrbit EndSubscribe EndCatalogue EndCode EndAllura; do
  .venv/bin/manim render -v WARNING -r 1920,1080 -s -t --media_dir comp/.media comp/plates.py "$s"
  cp "comp/.media/images/plates/${s}_ManimCE_v0.21.0.png" "comp/plates/$s.png"
done
rm -rf comp/.media
python3 -c "
import subprocess
raw = bytes([10, 12, 17]) * (1920 * 1080)
subprocess.run(['ffmpeg', '-v', 'error', '-y', '-f', 'rawvideo', '-pix_fmt', 'rgb24',
                '-s', '1920x1080', '-i', '-', '-frames:v', '1',
                'comp/plates/Backdrop.png'], input=raw, check=True)"
ls -la comp/plates/
