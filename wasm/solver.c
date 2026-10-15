// Red Sandbox WebAssembly module "solver": the gas pressure solve (preconditioned conjugate gradients over the
// active gas cells). The same operations in the same order as the JavaScript fallback, so the same results.
//
// index.html carries every compiled module as base64 (WASM_B64). After changing this file run wasm/build.sh,
// which compiles each module with clang and puts them into index.html.
// No C library: only what WebAssembly does natively (sqrt, abs, floor, min, max).

#define EXPORT(name) __attribute__((export_name(name)))
#define sqrt __builtin_sqrt
#define fabs __builtin_fabs
#define floor __builtin_floor

static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double hyp2(double a, double b) { return sqrt(a * a + b * b); }

// =====================================================================
//  Pressure solve
// =====================================================================
#define NMAX 9000
#define LMAX 8192

static double diag[NMAX], rhs[NMAX], x[NMAX], r[NMAX], z[NMAX], d[NMAX], q[NMAX];
static int nb[NMAX * 6];
static double nc[NMAX * 6];
static int lia[LMAX], lib[LMAX];
static double lc[LMAX];

EXPORT("diag") double *p_diag(void) { return diag; }
EXPORT("rhs") double *p_rhs(void) { return rhs; }
EXPORT("x") double *p_x(void) { return x; }
EXPORT("nb") int *p_nb(void) { return nb; }
EXPORT("nc") double *p_nc(void) { return nc; }
EXPORT("lia") int *p_lia(void) { return lia; }
EXPORT("lib") int *p_lib(void) { return lib; }
EXPORT("lc") double *p_lc(void) { return lc; }
EXPORT("nmax") int nmax(void) { return NMAX; }
EXPORT("lmax") int lmax(void) { return LMAX; }

static void mul(int n, int nl, const double *v, double *out) {
  for (int a = 0, k = 0; a < n; a++, k += 6) {
    double s = diag[a] * v[a];
    int j;
    j = nb[k];     if (j >= 0) s -= nc[k] * v[j];
    j = nb[k + 1]; if (j >= 0) s -= nc[k + 1] * v[j];
    j = nb[k + 2]; if (j >= 0) s -= nc[k + 2] * v[j];
    j = nb[k + 3]; if (j >= 0) s -= nc[k + 3] * v[j];
    j = nb[k + 4]; if (j >= 0) s -= nc[k + 4] * v[j];
    j = nb[k + 5]; if (j >= 0) s -= nc[k + 5] * v[j];
    out[a] = s;
  }
  for (int i = 0; i < nl; i++) { int ia = lia[i], ib = lib[i]; double c = lc[i]; out[ia] -= c * v[ib]; out[ib] -= c * v[ia]; }
}

// x holds the starting pressures on entry and the solution on return; returns the passes it took
EXPORT("solve") int solve(int n, int nl, double tol, int cmax) {
  mul(n, nl, x, q);
  double rz = 0;
  for (int a = 0; a < n; a++) {
    r[a] = rhs[a] - q[a];
    z[a] = r[a] / diag[a];
    rz += r[a] * z[a];
    d[a] = z[a];
  }
  int it = 0;
  for (; it < cmax && rz > tol; it++) {
    mul(n, nl, d, q);
    double dq = 0;
    for (int a = 0; a < n; a++) dq += d[a] * q[a];
    if (dq <= 0) break;
    double alpha = rz / dq, rz2 = 0;
    for (int a = 0; a < n; a++) {
      x[a] += alpha * d[a];
      double ra = r[a] - alpha * q[a];
      r[a] = ra;
      double za = ra / diag[a];
      z[a] = za;
      rz2 += ra * za;
    }
    double beta = rz2 / rz;
    rz = rz2;
    for (int a = 0; a < n; a++) d[a] = z[a] + beta * d[a];
  }
  return it;
}
