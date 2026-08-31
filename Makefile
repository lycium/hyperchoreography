# make            native-optimized build (MPFR/GMP needed for `refine`)
# make NOMPFR=1   without MPFR
# make NOACCEL=1  without Accelerate
CXX      ?= c++
ARCH     := $(shell uname -m)
CXXFLAGS ?= -std=c++20 -O3 -DNDEBUG -fno-math-errno -fno-trapping-math \
            -funroll-loops -fopenmp-simd -fomit-frame-pointer -Wall -Wextra -Wno-unused-parameter -pthread
ifeq ($(ARCH),arm64)
  CXXFLAGS += -mcpu=native                  # -march=native selects apple-m1 here
else
  CXXFLAGS += -march=native -mtune=native
endif
ifeq ($(ARCH),x86_64)
  CXXFLAGS += -mprefer-vector-width=512
endif
LDFLAGS  ?= -pthread
ifeq ($(shell uname -s),Darwin)
  ifndef NOACCEL
    CXXFLAGS += -DHAVE_ACCELERATE -DACCELERATE_NEW_LAPACK=1
    LDFLAGS  += -framework Accelerate
  endif
endif
ifndef NOMPFR
  BREW := $(shell brew --prefix 2>/dev/null)
  ifneq ($(BREW),)
    CXXFLAGS += -I$(BREW)/include
    LDFLAGS  += -L$(BREW)/lib
  endif
  CXXFLAGS += -DHAVE_MPFR
  LDFLAGS  += -lmpfr -lgmp
endif
SRC  := src/main.cpp
HDRS := $(wildcard src/*.hpp)

hyperchoreography: $(SRC) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

test: src/tests.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o hyperchoreography_test src/tests.cpp $(LDFLAGS) && ./hyperchoreography_test

# docs/ is what GitHub Pages serves.  The eight is kept out of the selection deliberately: it is the one
# orbit everybody has already seen.  The rest is the score's ranking overruled by eye.
gallery: hyperchoreography
	python3 tools/gallery.py --out docs/index.html \
	  --no-hero 'd2-3_n3.bin#0' --no-hero 'd5-6_n7.bin#100' --hero 'd2-4_n4.bin#6'

clean:
	rm -f hyperchoreography hyperchoreography_test
.PHONY: clean test gallery
