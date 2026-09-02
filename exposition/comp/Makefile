# The Allura Studio production of the exposition film.
#
#   make master      the lossless master (FFV1 + FLAC), ../out/master.mkv
#   make mp4         the lossless YouTube upload copy beside it
#   make verify      prove every scene's picture and sound bit-exact in the master
#
# The committed inputs are hyperchoreography.allura (the comp), spec.json,
# plates/ and chapters.txt; ../out/ carries the rendered scenes and everything
# this Makefile produces. Regenerating the comp itself (plates -> spec ->
# conform) needs the film's venv and the allura binary:
#
#   make plates spec conform
#
ALLURA ?= ../../../allura_studio/build/allura

master: ../out/master.mkv
../out/master.mkv: hyperchoreography.allura
	$(ALLURA) --open hyperchoreography.allura --comp hyperchoreography \
	    --render ../out/master.mkv --preset ffv1 --exact

# A lossless HEVC encode of the master in 10-bit 4:4:4. Not RGB: HEVC can carry the
# master's GBR planes bit-exactly, but must then signal matrix_coefficients = identity,
# and players that ignore that tag read G-B-R as Y-Cb-Cr and paint the film green. The
# BT.709 matrix is irrational and so never exactly invertible, but 10 bits of headroom
# cost at most one code value on 0.03 % of samples (8-bit costs two on 99.6 %; 12-bit is
# no better and twice the size). 4:2:0 would discard three quarters of the chroma before
# the encoder ran, so 4:4:4 stays. YouTube re-encodes whatever it is given; the upload's
# job is only to hand it the master's picture intact.
mp4: ../out/hyperchoreography-upload.mp4
../out/hyperchoreography-upload.mp4: ../out/master.mkv
	ffmpeg -y -i ../out/master.mkv \
	    -c:v libx265 -pix_fmt yuv444p10le \
	    -color_range pc -colorspace bt709 -color_primaries bt709 -color_trc bt709 \
	    -x265-params lossless=1:colorprim=bt709:transfer=bt709:colormatrix=bt709:range=full \
	    -tag:v hvc1 \
	    -c:a aac -b:a 384k -movflags +faststart \
	    ../out/hyperchoreography-upload.mp4

verify: ../out/master.mkv
	./verify_master.py

# The film itself, re-timed and re-rendered against its narration. The production
# voice is Kokoro's bm_daniel (neural, local; narration/.venv-tts); the film is
# PACED by the measured readings, so changing the voice re-renders everything.
film:
	cd .. && EXPO_VOICE=kokoro:bm_daniel .venv/bin/python render.py

plates:
	./make_plates.sh
spec:
	python3 build_spec.py
conform:
	$(ALLURA) --conform spec.json --save-project hyperchoreography.allura

.PHONY: master mp4 verify film plates spec conform
