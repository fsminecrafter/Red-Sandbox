#!/bin/sh
# Compile wasm/sim.c to WebAssembly and put it into index.html (as the base64 string SIM_WASM_B64).
# (bulk memory: memory.fill instead of a C library memset; every current browser has it)
# Needs clang with the wasm32 target and wasm-ld (on Debian/Ubuntu: apt install clang lld).
set -e
cd "$(dirname "$0")"
clang --target=wasm32 -O3 -mbulk-memory -nostdlib -Wl,--no-entry -Wl,--strip-all -o sim.wasm sim.c
B64=$(base64 -w0 sim.wasm 2>/dev/null || base64 sim.wasm | tr -d '\n')
python3 - "$B64" <<'PY'
import re, sys
p = '../index.html'
s = open(p).read()
s, n = re.subn(r"const SIM_WASM_B64 = '[A-Za-z0-9+/=]*';", "const SIM_WASM_B64 = '" + sys.argv[1] + "';", s)
assert n == 1, 'SIM_WASM_B64 not found in index.html'
open(p, 'w').write(s)
PY
rm -f sim.wasm
echo "sim.wasm: ${#B64} base64 characters, put into index.html"
