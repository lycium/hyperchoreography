#!/bin/sh
# Build a GMP and an MPFR that know what CPU they are on, static, into $PREFIX, which the
# Makefile picks up if it is there.  GMP's table in 6.3.0 stops at AMD family 25, so a Zen 4 or 5
# host gets the baseline k8 assembly, as does every MSYS2 build; `make mpinfo` reports which.
#
#   tools/deps.sh                       # $HOME/.local/opt/gmp-zen, CPU guessed
#   PREFIX=/opt/gmp tools/deps.sh       # elsewhere
#   CPU=zen4 tools/deps.sh              # override the guess (gmp-6.3.0/configure lists the names)
#   NOTUNE=1 tools/deps.sh              # skip MPFR's threshold tuning, the slow half
#
# Needs a compiler, m4, make, curl, 700 MB of scratch.  Both test suites run before either is
# installed: this library decides whether a certificate in paper/ is true.

set -e
PREFIX=${PREFIX:-$HOME/.local/opt/gmp-zen}
WORK=${WORK:-${TMPDIR:-/tmp}/hc-deps-$$}
GMP=gmp-6.3.0
MPFR=mpfr-4.2.2
GMP_SHA=a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
MPFR_SHA=b67ba0383ef7e8a8563734e2e889ef5ec3c3b898a01d00fa0a6869ad81c6ce01
JOBS=${JOBS:-$( (nproc || sysctl -n hw.ncpu) 2>/dev/null || echo 4)}

# Under C23, GCC 14's default, a GMP 6.3.0 configure probe calls a `void g(){}` with six
# arguments; configure then reports finding no working compiler at all.
if [ -z "$CC" ]; then command -v cc >/dev/null 2>&1 && CC=cc || CC=gcc; fi
case $($CC -dumpversion 2>/dev/null | cut -d. -f1 || echo 0) in
  1[4-9]|[2-9][0-9]) CC="$CC -std=gnu17" ;;
esac

# .part until complete: gmplib.org drops connections halfway, and the next run would take the
# truncated file for one it already has.
fetch() { # url mirror sha name
  if [ ! -f "$4" ]; then
    curl -fLsS --retry 3 --retry-delay 2 -o "$4.part" "$1" ||
    curl -fLsS --retry 3 --retry-delay 2 -o "$4.part" "$2"
    mv "$4.part" "$4"
  fi
  have=$( (sha256sum "$4" 2>/dev/null || shasum -a 256 "$4") | cut -d' ' -f1)
  [ "$have" = "$3" ] || { rm -f "$4"; echo "$4: sha256 $have, expected $3" >&2; exit 1; }
}

mkdir -p "$WORK"; cd "$WORK"
fetch "https://gmplib.org/download/gmp/$GMP.tar.xz" "https://ftp.gnu.org/gnu/gmp/$GMP.tar.xz" \
      "$GMP_SHA" "$GMP.tar.xz"
fetch "https://www.mpfr.org/$MPFR/$MPFR.tar.xz" "https://ftp.gnu.org/gnu/mpfr/$MPFR.tar.xz" \
      "$MPFR_SHA" "$MPFR.tar.xz"
[ -d "$GMP" ]  || tar xf "$GMP.tar.xz"
[ -d "$MPFR" ] || tar xf "$MPFR.tar.xz"

# Correct config.guess only past the end of its table, by feature flag rather than family so
# this does not go stale the same way.  GMP's zen4 resolves to its zen3 assembly.
TRIPLET=$("$WORK/$GMP/config.guess")
if [ -z "$CPU" ]; then
  case $TRIPLET in x86_64-*)
    if [ -r /proc/cpuinfo ]; then
      vendor=$(grep -m1 vendor_id /proc/cpuinfo || true); flags=$(grep -m1 '^flags' /proc/cpuinfo || true)
    else
      vendor=$(sysctl -n machdep.cpu.vendor 2>/dev/null || true)
      flags=$(sysctl -n machdep.cpu.leaf7_features 2>/dev/null | tr 'A-Z' 'a-z' || true)
    fi
    has() { echo " $flags " | grep -q " $1 "; }
    if echo "$vendor" | grep -q AuthenticAMD; then
      if has avx512f; then CPU=zen4; elif has adx && has bmi2; then CPU=zen2; fi
    elif echo "$vendor" | grep -q GenuineIntel; then
      if has avx512f; then CPU=skylake; elif has adx && has bmi2; then CPU=coreibwl; fi
    fi ;;
  esac
fi
if [ -n "$CPU" ]; then TRIPLET=$CPU-$(echo "$TRIPLET" | cut -d- -f2-); fi
echo "=== building for $TRIPLET into $PREFIX ($JOBS jobs) ==="

cd "$WORK/$GMP"
[ -f Makefile ] || ./configure --prefix="$PREFIX" --build="$TRIPLET" \
                    --disable-shared --enable-static CC="$CC"
grep -m1 '^path=' config.log || true          # the mpn directories it settled on
make -j"$JOBS"
make -j"$JOBS" check
make install
make -C tune libspeed.la                        # MPFR's tuneup links against it
cp tune/.libs/libspeed.a tune/ 2>/dev/null || true

cd "$WORK/$MPFR"
# --with-gmp-build, not --with-gmp: MPFR then calls GMP's internal mpn entry points directly.
[ -f Makefile ] || ./configure --prefix="$PREFIX" --with-gmp-build="$WORK/$GMP" \
                    --disable-shared --enable-static --enable-thread-safe \
                    CC="$CC" CFLAGS="-O3 -march=native -mtune=native -fomit-frame-pointer"
make -j"$JOBS"
if [ -z "$NOTUNE" ]; then
  echo "=== tuning MPFR's thresholds (20 minutes or so; NOTUNE=1 skips it) ==="
  if (cd tune && make tune); then make -j"$JOBS"
  else echo "!! tuning failed; keeping MPFR's default thresholds" >&2; fi
fi
make -j"$JOBS" check
make install

echo
echo "=== installed into $PREFIX; now: make clean && make && make mpinfo ==="
