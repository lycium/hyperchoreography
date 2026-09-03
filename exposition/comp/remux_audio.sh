#!/bin/sh
# Rebuild the upload's soundtrack without re-encoding a frame. A narration change
# moves no picture, so this is minutes rather than hours.
# Usage: ./remux_audio.sh [upload.mp4]
set -eu
cd "$(dirname "$0")"

ALLURA=${ALLURA:-../../../allura_studio/build/allura}
OUT=${1:-../out/hyperchoreography-upload.mp4}
[ -f "$OUT" ] || { echo "no upload to remux at $OUT -- run make upload first" >&2; exit 1; }

work=$(mktemp -d)
allura=
trap 'test -n "$allura" && kill $allura 2>/dev/null; rm -rf "$work"; true' EXIT INT TERM

echo "bouncing the comp's audio"
"$ALLURA" --open hyperchoreography.allura --comp hyperchoreography \
    --render "$work/bounce.mkv" --preset ffv1 --exact > "$work/allura.log" 2>&1 &
allura=$!

while ! grep -q "^bounced " "$work/allura.log" 2>/dev/null; do
    kill -0 $allura 2>/dev/null || { cat "$work/allura.log" >&2; exit 1; }
    sleep 1
done
kill $allura 2>/dev/null || true
wait $allura 2>/dev/null || true
grep "^bounced " "$work/allura.log"

wav="$work/bounce.audio.wav"
[ -s "$wav" ] || { echo "the bounce wrote nothing" >&2; exit 1; }

echo "remuxing against the picture already in $(basename "$OUT")"
ffmpeg -hide_banner -loglevel error -y -i "$OUT" -i "$wav" \
    -map 0:v:0 -map 1:a:0 -c:v copy -c:a aac -b:a 384k -movflags +faststart \
    "$work/out.mp4"

before=$(ffmpeg -v error -i "$OUT" -map 0:v:0 -c copy -f md5 -)
after=$(ffmpeg -v error -i "$work/out.mp4" -map 0:v:0 -c copy -f md5 -)
[ "$before" = "$after" ] || { echo "FAIL: the picture changed ($before vs $after)" >&2; exit 1; }
echo "  ok   the video stream is bit-identical ($after)"

want=$(python3 -c 'import json;print(json.load(open("timeline.json"))["master_frames"])')
mv "$work/out.mp4" "$OUT"
got=$(ffprobe -v error -select_streams v:0 -count_packets \
      -show_entries stream=nb_read_packets -of csv=p=0 "$OUT")
echo "  ok   $got frames (timeline says $want)"
[ "$got" = "$want" ] || { echo "FAIL: wrong length" >&2; exit 1; }
echo "  $(du -h "$OUT" | cut -f1)  $OUT"
exit 0
