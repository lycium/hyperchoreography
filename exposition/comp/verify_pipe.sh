#!/bin/sh
# Prove the fifo carries the file's own bytes, the encode says studio range,
# and the RGB round trip is exact.
set -eu
cd "$(dirname "$0")"

ALLURA=${ALLURA:-../../../allura_studio/build/allura}
FRAMES=${FRAMES:-90}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT INT TERM

hevc() {                    # $1 in, $2 out -- the upload's own encoder settings
    ffmpeg -hide_banner -loglevel error -y -i "$1" \
        -c:v libx265 -pix_fmt yuv444p10le \
        -color_range tv -colorspace bt709 -color_primaries bt709 -color_trc bt709 \
        -x265-params lossless=1:colorprim=bt709:transfer=bt709:colormatrix=bt709:range=limited:log-level=error \
        -tag:v hvc1 -an "$2"
}
hashes() {
    ffmpeg -v error -i "$1" -map 0:v:0 -f framemd5 - | grep -v '^#'
}

mkfifo "$work/pipe.mkv"
"$ALLURA" --open hyperchoreography.allura --comp hyperchoreography \
    --render "$work/pipe.mkv" --preset ffv1 --exact --range 0 "$FRAMES" \
    > "$work/a.log" 2>&1 &
hevc "$work/pipe.mkv" "$work/piped.mp4"
wait

"$ALLURA" --open hyperchoreography.allura --comp hyperchoreography \
    --render "$work/file.mkv" --preset ffv1 --exact --range 0 "$FRAMES" \
    > "$work/b.log" 2>&1
hevc "$work/file.mkv" "$work/filed.mp4"

if [ "$(hashes "$work/piped.mp4" | md5)" = "$(hashes "$work/filed.mp4" | md5)" ]; then
    echo "  ok   $FRAMES frames through the fifo are the file's own bytes"
else
    echo "  FAIL the fifo changed the picture" >&2
    exit 1
fi

range=$(ffprobe -v error -select_streams v:0 -show_entries stream=color_range \
        -of csv=p=0 "$work/filed.mp4")
if [ "$range" = "tv" ]; then
    echo "  ok   the encode says studio range, which is what its samples are"
else
    echo "  FAIL the encode is tagged '$range', not tv" >&2
    exit 1
fi

ffmpeg -v error -i "$work/file.mkv" -frames:v 6 -f rawvideo -pix_fmt rgb24 "$work/src.raw"
ffmpeg -v error -i "$work/filed.mp4" -frames:v 6 -f rawvideo -pix_fmt rgb24 "$work/back.raw"
../.venv/bin/python round_trip.py "$work/src.raw" "$work/back.raw"
