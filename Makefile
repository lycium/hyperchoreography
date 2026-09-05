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
  # A GMP that does not recognise the CPU falls back to its baseline k8 assembly and says nothing.
  # tools/deps.sh builds one that does, into MPPREFIX; `make mpinfo` says which you ended up with.
  MPPREFIX ?= $(HOME)/.local/opt/gmp-zen
  ifneq ($(wildcard $(MPPREFIX)/include/mpfr.h),)
    CXXFLAGS += -isystem $(MPPREFIX)/include
    LDFLAGS  += -L$(MPPREFIX)/lib
  else
    BREW := $(shell brew --prefix 2>/dev/null)
    ifneq ($(BREW),)
      CXXFLAGS += -I$(BREW)/include
      LDFLAGS  += -L$(BREW)/lib
    endif
  endif
  CXXFLAGS += -DHAVE_MPFR
  LDFLAGS  += -lmpfr -lgmp
endif
SRC  := src/main.cpp
HDRS := $(wildcard src/*.hpp)

hyperchoreography: $(SRC) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

test: src/tests.cpp src/audit_tests.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o hyperchoreography_test src/tests.cpp $(LDFLAGS) && ./hyperchoreography_test
	$(CXX) $(CXXFLAGS) -o hyperchoreography_audit_test src/audit_tests.cpp $(LDFLAGS) && ./hyperchoreography_audit_test

audit-test: src/audit_tests.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o hyperchoreography_audit_test src/audit_tests.cpp $(LDFLAGS) && ./hyperchoreography_audit_test

audit-bench: tools/audit_bench.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o hyperchoreography_audit_bench tools/audit_bench.cpp $(LDFLAGS) && ./hyperchoreography_audit_bench

hyperchoreography_reference: tools/reference_import.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ tools/reference_import.cpp $(LDFLAGS)

# Same flags as the binary above, so it links and reports on the same library.
mpinfo: tools/mpinfo.cpp
	$(CXX) $(CXXFLAGS) -o hyperchoreography_mpinfo tools/mpinfo.cpp $(LDFLAGS) && ./hyperchoreography_mpinfo

# docs/ is what GitHub Pages serves.  The eight is kept out of the selection deliberately: it is the one
# orbit everybody has already seen.  The rest is the score's ranking overruled by eye.
gallery: hyperchoreography
	python3 tools/gallery.py --out docs/index.html \
	  --no-hero 'd2-3_n3.bin#0' --no-hero 'd5-6_n7.bin#100' --hero 'd2-4_n4.bin#6' \
	  --hero 'li-liao-2025-n3.bin#62'

gallery-check: hyperchoreography
	python3 tools/check_gallery.py docs/index.html

clean:
	rm -f hyperchoreography hyperchoreography_test hyperchoreography_mpinfo hyperchoreography_audit_test hyperchoreography_audit_bench hyperchoreography_reference
.PHONY: clean test audit-test audit-bench gallery gallery-check mpinfo
