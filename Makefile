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

clean:
	rm -f hyperchoreography hyperchoreography_test
.PHONY: clean test
