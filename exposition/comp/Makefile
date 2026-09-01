# The Allura Studio production of the exposition film.
#
#   make master      the lossless master (FFV1 + FLAC), ../out/master.mkv
#   make mp4         the YouTube upload copy beside it
#   make verify      prove every scene's picture and sound bit-exact in the master
#
# The committed inputs are hyperchoreography.allura (the comp), spec.json,
# plates/ and chapters.txt; ../out/ carries the rendered scenes and everything
# this Makefile produces. Regenerating the comp itself (plates -> spec ->
# conform) needs the film's venv and the allura binary:
#
#   make plates spec conform
#
ALLURA ?= ../../allura_studio/build/allura

master: ../out/master.mkv
../out/master.mkv: hyperchoreography.allura
	$(ALLURA) --open hyperchoreography.allura --comp hyperchoreography \
	    --render ../out/master.mkv --preset ffv1 --exact

# Full-range RGB in, tagged limited-range bt709 4:2:0 out; the master's FLAC
# becomes AAC. -crf 16 at 1080p30 for a screencast-like picture is transparent.
mp4: ../out/hyperchoreography-upload.mp4
../out/hyperchoreography-upload.mp4: ../out/master.mkv
	ffmpeg -y -i ../out/master.mkv \
	    -c:v libx264 -preset slow -crf 16 -pix_fmt yuv420p \
	    -vf scale=in_range=full:out_range=tv \
	    -color_primaries bt709 -color_trc bt709 -colorspace bt709 \
	    -x264-params colorprim=bt709:transfer=bt709:colormatrix=bt709 \
	    -c:a aac -b:a 192k -movflags +faststart \
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
