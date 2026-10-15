// Red Sandbox WebAssembly module "gas": the per-cell chemistry of the gas grid, the parts of gasStep in
// index.html that only look at a cell and its neighbours:
//
//  * mixStep:   heat and fuel spreading between neighbouring cells, and smoke thinning out into the open air
//  * burnStep:  fuel and air burning next to a flame (or once hot enough for long enough), then cooling
//  * steamStep: steam condensing back into water as it cools
//
// The gas arrays themselves live in this module's memory: index.html points its gas.air, gas.fuel, ... at
// them when this module loads, so the JavaScript code and these functions work on the very same numbers.
// The sums are done in doubles and stored as floats, in the same order as the JavaScript fallback (which
// does exactly that with Float32Arrays), so both give the same results.
//
// index.html carries every compiled module as base64 (WASM_B64). After changing this file run wasm/build.sh.
// No C library.

#define EXPORT(name) __attribute__((export_name(name)))
#define NMAX 9000 // gas cells (GAS_MAX in index.html)
#define GH 0.1

static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dabs(double a) { return a < 0 ? -a : a; }

static float air[NMAX], fuel[NMAX], exh[NMAX], en[NMAX], stm[NMAX], ofu[NMAX], osm[NMAX];
static float burn[NMAX], ign[NMAX], f[NMAX], lw[NMAX], amb[NMAX];
static unsigned char flame[NMAX], flameNext[NMAX], copen[NMAX];
static int ci[NMAX], cj[NMAX], ck[NMAX], nb[NMAX * 6], faceNb[NMAX * 6], root[NMAX], act[NMAX];
static float faceC[NMAX * 6], faceA[NMAX * 6];
// The constants, set once by index.html (so they are only written down in one place)
//  [0] FMIN [1] T_IGN [2] T_FLAME [3] T_LIGHT [4] IGN_RATE [5] AFR [6] HEAT [7] LEAN
//  [8] T_COND [9] COND_RATE [10] STEAM_X [11] LATENT
static double K[16];
// Running totals burnStep adds to: [0] heat released [1..3] flash position sum [4] flash weight
static double acc[8];

EXPORT("nmax") int nmax(void) { return NMAX; }
EXPORT("air") float *p_air(void) { return air; }
EXPORT("fuel") float *p_fuel(void) { return fuel; }
EXPORT("exh") float *p_exh(void) { return exh; }
EXPORT("en") float *p_en(void) { return en; }
EXPORT("stm") float *p_stm(void) { return stm; }
EXPORT("ofu") float *p_ofu(void) { return ofu; }
EXPORT("osm") float *p_osm(void) { return osm; }
EXPORT("burn") float *p_burn(void) { return burn; }
EXPORT("ign") float *p_ign(void) { return ign; }
EXPORT("f") float *p_f(void) { return f; }
EXPORT("lw") float *p_lw(void) { return lw; }
EXPORT("amb") float *p_amb(void) { return amb; }
EXPORT("flame") unsigned char *p_flame(void) { return flame; }
EXPORT("flameNext") unsigned char *p_flameNext(void) { return flameNext; }
EXPORT("copen") unsigned char *p_copen(void) { return copen; }
EXPORT("ci") int *p_ci(void) { return ci; }
EXPORT("cj") int *p_cj(void) { return cj; }
EXPORT("ck") int *p_ck(void) { return ck; }
EXPORT("nb") int *p_nb(void) { return nb; }
EXPORT("faceNb") int *p_faceNb(void) { return faceNb; }
EXPORT("root") int *p_root(void) { return root; }
EXPORT("act") int *p_act(void) { return act; }
EXPORT("faceC") float *p_faceC(void) { return faceC; }
EXPORT("faceA") float *p_faceA(void) { return faceA; }
EXPORT("consts") double *p_consts(void) { return K; }
EXPORT("acc") double *p_acc(void) { return acc; }

static inline double total(int c) { return (double)air[c] + (double)fuel[c] + (double)exh[c] + (double)stm[c]; }
static int find(int c) { while (root[c] != c) { root[c] = root[root[c]]; c = root[c]; } return c; }

// the oily share goes along with fuel or exhaust that moves from a to b
static void oilShare(float *arr, const float *of, int a, int b, double ds) {
  int from = ds > 0 ? a : b;
  double base = of[from];
  if (base <= 1e-9 || arr[from] <= 0) return;
  double v = dmin(arr[from], ((double)arr[from] / base) * dabs(ds));
  arr[from] -= v; arr[ds > 0 ? b : a] += v;
}
static void mix(int c, int n, double fa, double kt, double kd) {
  double nc = total(c), nn = total(n);
  if (nc <= 1e-9 || nn <= 1e-9) return;
  double m = dmin(nc, nn) * 0.5 * fa;
  double dE = kt * (en[c] / nc - en[n] / nn) * m;
  en[c] -= dE; en[n] += dE;
  double ds = kd * (fuel[c] / nc - fuel[n] / nn) * m;
  if (ds > 0) ds = dmin(ds, air[n]); else ds = dmax(ds, -(double)air[c]);
  if ((double)ofu[c] + (double)ofu[n] > 0) oilShare(ofu, fuel, c, n, ds);
  fuel[c] -= ds; fuel[n] += ds; air[c] += ds; air[n] -= ds;
  ds = kd * (exh[c] / nc - exh[n] / nn) * m;
  if (ds > 0) ds = dmin(ds, air[n]); else ds = dmax(ds, -(double)air[c]);
  if ((double)osm[c] + (double)osm[n] > 0) oilShare(osm, exh, c, n, ds);
  exh[c] -= ds; exh[n] += ds; air[c] += ds; air[n] -= ds;
  if ((double)stm[c] + (double)stm[n] <= 1e-9) return;
  ds = kd * (stm[c] / nc - stm[n] / nn) * m;
  if (ds > 0) ds = dmin(ds, air[n]); else ds = dmax(ds, -(double)air[c]);
  stm[c] -= ds; stm[n] += ds; air[c] += ds; air[n] -= ds;
}

// Step 5 of gasStep: heat and fuel spread between neighbours; cool smoke in the open thins out
EXPORT("mixStep") void mixStep(int n, double kt, double kd, double dis, double disT) {
  const double FMIN = K[0], T_IGN = K[1], LEAN = K[7];
  for (int a = 0; a < n; a++) {
    int c = act[a], base = c * 6;
    for (int d = 0; d < 6; d += 2) {
      int o = faceNb[base + d];
      if (faceC[base + d] != 0 && o >= 0) mix(c, o, faceA[base + d], kt, kd);
    }
    int openCh = amb[c] == 0 && f[c] >= FMIN && copen[find(c)];
    if ((amb[c] != 0 || openCh) && !flame[c] && en[c] < T_IGN * total(c)) {
      double tn = total(c), dd = dis;
      double vf = fuel[c] * dd * (fuel[c] < LEAN * tn ? 1 : 0.15), vx = exh[c] * dd, vs = stm[c] * dd;
      if (ofu[c] > 0) ofu[c] *= 1 - vf / dmax(fuel[c], 1e-9);
      if (osm[c] > 0) osm[c] *= 1 - vx / dmax(exh[c], 1e-9);
      fuel[c] -= vf; exh[c] -= vx; stm[c] -= vs; air[c] += vf + vx + vs;
      en[c] += (tn - en[c]) * disT;
    }
  }
}

// Step 6 of gasStep: burning, then cooling towards the temperature of the air around
EXPORT("burnStep") void burnStep(int n, double h, double burnK, double coolK) {
  const double T_IGN = K[1], T_FLAME = K[2], T_LIGHT = K[3], IGN_RATE = K[4], AFR = K[5], HEAT = K[6];
  for (int a = 0; a < n; a++) {
    int c = act[a];
    burn[c] *= 0.9;
    double tn = total(c);
    if (tn <= 1e-9) continue;
    double T = en[c] / tn;
    if (T > T_IGN && T <= T_FLAME) ign[c] += h * (T - T_IGN) * IGN_RATE; else if (T <= T_IGN) ign[c] *= 0.9;
    int lit = T > T_FLAME || (T > T_IGN && ign[c] >= 1);
    if (!lit && T > T_LIGHT) for (int d = 0; d < 6 && !lit; d++) { int o = nb[c * 6 + d]; if (o >= 0 && flame[o]) lit = 1; }
    flameNext[c] = 0;
    if (lit && fuel[c] > 1e-6 && air[c] > 1e-6) {
      double b = dmin(fuel[c], air[c] / AFR) * burnK;
      if (ofu[c] > 0) { double bo = dmin(ofu[c], (b * ofu[c]) / fuel[c]); ofu[c] -= bo; osm[c] += bo * (1 + AFR) * 1.05; }
      fuel[c] -= b;
      air[c] -= b * AFR;
      exh[c] += b * (1 + AFR) * 1.05;
      double e = b * HEAT;
      en[c] += e;
      burn[c] += e;
      if (e > 0.02) flameNext[c] = 1;
      acc[0] += e;
      acc[1] += (ci[c] + 0.5) * GH * e; acc[2] += (cj[c] + 0.5) * GH * e; acc[3] += (ck[c] + 0.5) * GH * e;
      acc[4] += e;
    }
    en[c] += (total(c) - en[c]) * coolK;
  }
  for (int a = 0; a < n; a++) { int c = act[a]; flame[c] = flameNext[c]; }
}

// Steam that has cooled down turns back into water, giving up its heat
EXPORT("steamStep") void steamStep(int n, double h) {
  const double T_COND = K[8], COND_RATE = K[9], STEAM_X = K[10], LATENT = K[11];
  for (int a = 0; a < n; a++) {
    int c = act[a];
    if (stm[c] <= 1e-6) continue;
    double tn = total(c), T = tn > 1e-6 ? en[c] / tn : 1;
    if (T >= T_COND) continue;
    double dc = stm[c] * dmin(0.5, COND_RATE * h * (T_COND - T));
    stm[c] -= dc;
    en[c] -= dc * T;
    en[c] += (dc / STEAM_X) * LATENT * 0.5;
    lw[c] += dc / STEAM_X;
  }
}
