#!/bin/sh
# Build the environment the presentation renders in.
#
#   ./setup.sh          make .venv and install manim
#
# On macOS pkgconf and pango are needed to build pycairo and ManimPango; ffmpeg
# does the time filtering and the encoding. If Homebrew is not present, install
# the three of them however your system prefers.
set -e
cd "$(dirname "$0")"

if command -v brew >/dev/null 2>&1; then
  brew list --formula 2>/dev/null | grep -qx pkgconf || brew install pkgconf
  brew list --formula 2>/dev/null | grep -qx pango   || brew install pango
  brew list --formula 2>/dev/null | grep -qx ffmpeg  || brew install ffmpeg
fi

python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
PKG_CONFIG_PATH="${PKG_CONFIG_PATH:-/opt/homebrew/lib/pkgconfig}" \
  .venv/bin/pip install -r requirements.txt

if [ -x ../hyperchoreography ]; then
  .venv/bin/python make_data.py
else
  echo "note: no solver binary at ../hyperchoreography, so data/ was not refreshed"
  echo "      (the presentation renders from the files already in data/)"
fi

echo
echo "ready:  .venv/bin/python render.py --preview s00_open"
