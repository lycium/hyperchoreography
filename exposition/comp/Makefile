# The Allura Studio production of the exposition film.
#
#   make upload      the deliverable: ../out/hyperchoreography-upload.mp4, lossless
#   make verify-pipe prove the one-pass route carries the comp's own bytes
#   make audio       re-lay the soundtrack only, when just the narration changed
#   make master      the archival FFV1 + FLAC master, when you want one on disk
#   make verify      prove every scene's picture and sound bit-exact in that master
#
# `upload` is one pass with nothing between the comp and the encoder: Allura writes
# FFV1 into a fifo and x265 reads it out, so the film is delivered without a
# twenty-two gigabyte intermediate in a container no player wants. The master is
# still there for the archive and for `verify`, which needs to seek.
#
# The committed inputs are hyperchoreography.allura (the comp), spec.json,
# chapters.txt and the plates' PNGs; ../out/ carries the rendered scenes and
# everything this Makefile produces, and so does the one plate that is a render
# rather than a still (plates/EndOrbit.mkv, the end card's turning orbit).
# Regenerating the comp itself (plates -> spec -> conform) needs the film's venv
# and the allura binary:
#
#   make plates spec conform
#
ALLURA ?= ../../../allura_studio/build/allura
export ALLURA

UPLOAD = ../out/hyperchoreography-upload.mp4

upload: $(UPLOAD)
$(UPLOAD): hyperchoreography.allura make_upload.sh
	./make_upload.sh $@

verify-pipe:
	./verify_pipe.sh

# Only the narration changed: rebuild the upload's soundtrack against the picture it
# already carries. A pronunciation fix or a real recording moves no frame, so this is
# a quarter of an hour rather than four.
audio:
	./remux_audio.sh $(UPLOAD)

master: ../out/master.mkv
../out/master.mkv: hyperchoreography.allura
	$(ALLURA) --open hyperchoreography.allura --comp hyperchoreography \
	    --render ../out/master.mkv --preset ffv1 --exact

verify: ../out/master.mkv
	./verify_master.py

# The film itself, re-timed and re-rendered against its narration. The production
# voice is Kokoro's bm_daniel (neural, local; narration/.venv-tts); the film is
# PACED by the measured readings, so changing the voice re-renders everything.
# One process per scene and mostly one core each, longest first.
film:
	cd .. && EXPO_VOICE=kokoro:bm_daniel .venv/bin/python render.py \
	    --res 4k --ss 16 --jobs 8

plates:
	./make_plates.sh
spec:
	python3 build_spec.py
conform:
	$(ALLURA) --conform spec.json --save-project hyperchoreography.allura

.PHONY: upload verify-pipe audio master verify film plates spec conform
