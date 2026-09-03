#!/bin/sh
# Render the comp straight to the lossless upload, through a fifo.
# Usage: ./make_upload.sh [output.mp4]
#
# Studio range, not full: the video_full_range_flag does not reliably reach the
# file, and full-range samples read as studio lose everything below code 16 --
# which is this film's background. Do not change tv back to pc.
set -eu
cd "$(dirname "$0")"

ALLURA=${ALLURA:-../../../allura_studio/build/allura}
OUT=${1:-../out/hyperchoreography-upload.mp4}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT INT TERM
pipe="$work/master.mkv"
mkfifo "$pipe"

echo "rendering -> $OUT (lossless HEVC 4:4:4 10-bit studio range)"
"$ALLURA" --open hyperchoreography.allura --comp hyperchoreography \
    --render "$pipe" --preset ffv1 --exact > "$work/allura.log" 2>&1 &
allura=$!

ffmpeg -hide_banner -loglevel error -y -i "$pipe" \
    -c:v libx265 -pix_fmt yuv444p10le \
    -color_range tv -colorspace bt709 -color_primaries bt709 -color_trc bt709 \
    -x265-params lossless=1:colorprim=bt709:transfer=bt709:colormatrix=bt709:range=limited:log-level=error \
    -tag:v hvc1 \
    -c:a aac -b:a 384k -movflags +faststart \
    "$OUT"

if ! wait "$allura"; then
    echo "allura failed:" >&2
    cat "$work/allura.log" >&2
    exit 1
fi
tail -n 1 "$work/allura.log"

want=$(python3 -c 'import json;print(json.load(open("timeline.json"))["master_frames"])')
got=$(ffprobe -v error -select_streams v:0 -count_packets \
      -show_entries stream=nb_read_packets -of csv=p=0 "$OUT")
fps=$(ffprobe -v error -select_streams v:0 \
      -show_entries stream=r_frame_rate -of csv=p=0 "$OUT")
echo "  $got frames at $fps fps  (timeline says $want)"
[ "$got" = "$want" ] || { echo "FAIL: the upload is not the whole film" >&2; exit 1; }
echo "  $(du -h "$OUT" | cut -f1)  $OUT"
