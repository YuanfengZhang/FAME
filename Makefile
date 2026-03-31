# FAME Makefile
# Supports automatic detection of zstd library

OBJECTS_BASE=gzstream.o RefReader_istr.o RefGenome.o DnaBitStr.o main.o\
		ReadQueue.o Read.o CONST.o ShiftAnd.o LevenshtDP.o
PROGNAME=FAME
CXX=g++

# Compiler flags
CXXFLAGS= -std=c++14 -ggdb -Wshadow -Wall -pedantic -pipe -O3 -fopenmp -march=native \
		-I ./sparsehash/include/usr/local/include/ -I ./hopscotch-map/include/tsl/

# Always need zlib for gzip support
GZFLAGS= -lz

# Detect zstd library
# Try pkg-config first, then fallback to manual check
ZSTD_CHECK := $(shell pkg-config --exists libzstd 2>/dev/null && echo "yes")

ifeq ($(ZSTD_CHECK),yes)
	# Use pkg-config for zstd
	ZSTDFLAGS= $(shell pkg-config --libs libzstd)
	ZSTD_CFLAGS= $(shell pkg-config --cflags libzstd)
	CXXFLAGS += -DZSTD_SUPPORT $(ZSTD_CFLAGS)
	OBJECTS= zstdstream.o $(OBJECTS_BASE)
	ZSTD_STATUS=enabled
else
	# Try manual check for libzstd
	ZSTD_MANUAL := $(shell echo '\#include <zstd.h>' | $(CXX) -E -x c++ - >/dev/null 2>&1 && echo "yes")
	ifeq ($(ZSTD_MANUAL),yes)
		ZSTDFLAGS= -lzstd
		CXXFLAGS += -DZSTD_SUPPORT
		OBJECTS= zstdstream.o $(OBJECTS_BASE)
		ZSTD_STATUS=enabled
	else
		ZSTDFLAGS=
		OBJECTS= $(OBJECTS_BASE)
		ZSTD_STATUS=disabled
	endif
endif

.PHONY: all clean profile check-deps

all: check-deps ${PROGNAME}

check-deps:
	@echo "=== Checking dependencies ==="
	@echo "C++ compiler: $(CXX)"
	@echo "zlib: required (gzip support)"
	@echo "zstd: $(ZSTD_STATUS)"
ifeq ($(ZSTD_STATUS),disabled)
	@echo "  (install libzstd-dev to enable)"
endif
	@echo ""

profile: ${PROGNAME}Profile

gzstream.o: gzstream/gzstream.C gzstream/gzstream.h
	${CXX} ${CXXFLAGS} -c $<

# Only compile zstdstream if zstd is available
ifeq ($(ZSTD_STATUS),enabled)
zstdstream.o: zstdstream/zstdstream.C zstdstream/zstdstream.h
	${CXX} ${CXXFLAGS} -c $<
endif

ReadQueue.o: ReadQueue.cpp ReadQueue.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

${PROGNAME}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} ${GZFLAGS} ${ZSTDFLAGS} -o $@

${PROGNAME}Profile: ${OBJECTS}
	${CXX} ${CXXFLAGS} -pg -rdynamic ${OBJECTS} ${GZFLAGS} ${ZSTDFLAGS} -o $@

clean:
	rm -f ${OBJECTS} zstdstream.o ${PROGNAME} ${PROGNAME}Profile