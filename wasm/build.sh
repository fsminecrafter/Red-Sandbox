#!/bin/sh
# Compile every WebAssembly module (wasm/*.c) and put them into index.html as base64, in the WASM_B64 table.
# Each module loads on its own: one that fails falls back to JavaScript without taking the others with it.
#   solver.c    the gas pressure solve
#   geometry.c  parts as distance shapes, and how solid each gas cell is
#   gas.c       gas mixing, burning and steam (owns the gas arrays)
# (bulk memory: memory.fill instead of a C library memset; every current browser has it)
# Needs clang with the wasm32 target and wasm-ld (on Debian/Ubuntu: apt install clang lld).
set -e
cd "$(dirname "$0")"
for m in solver geometry gas; do
  clang --target=wasm32 -O3 -mbulk-memory -ffp-contract=off -nostdlib -Wl,--no-entry -Wl,--strip-all -o "$m.wasm" "$m.c"
done
python3 - solver geometry gas <<'PY'
import base64, re, sys
p = '../index.html'
s = open(p).read()
for m in sys.argv[1:]:
    b64 = base64.b64encode(open(m + '.wasm', 'rb').read()).decode()
    s, n = re.subn(r"(\n  %s: )'[A-Za-z0-9+/=]*'," % m, lambda mo: mo.group(1) + "'" + b64 + "',", s)
    assert n == 1, m + ' not found in the WASM_B64 table in index.html'
    print('%s.wasm: %d base64 characters' % (m, len(b64)))
open(p, 'w').write(s)
PY
rm -f solver.wasm geometry.wasm gas.wasm
echo "put into index.html"
