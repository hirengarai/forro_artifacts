# Build both verification programs:  make
# Everything they need is in this folder; there is nothing to install.

CXX      ?= c++
UNAME_M  := $(shell uname -m)
ifneq (,$(filter arm64 aarch64,$(UNAME_M)))
  NATIVE := -mcpu=native
else
  NATIVE := -march=native
endif
# Only pass the native flag if this compiler accepts it. Never force -mavx2: on a
# machine without AVX2 that still compiles and then dies with SIGILL at run time,
# whereas letting -march=native decide leaves forro.cpp to pick its scalar
# fallback, which is slower but correct. Override with  make ARCHFLAGS=...
ARCHFLAGS ?= $(shell $(CXX) $(NATIVE) -E -x c++ /dev/null >/dev/null 2>&1 && echo $(NATIVE))
CXXFLAGS ?= -std=c++20 -Wall -Wextra -O3 $(ARCHFLAGS)
LDFLAGS  ?= -pthread

DEPS = forro.cpp verify.hpp syncopation.hpp
BINS = run65 run675

all: $(BINS)

run65: forro65_verification.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

run675: forro675_verification.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(BINS)

.PHONY: all clean
