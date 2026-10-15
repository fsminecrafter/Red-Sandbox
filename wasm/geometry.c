// Red Sandbox WebAssembly module "geometry": the parts as signed distance shapes sorted into half-metre
// buckets, and measuring how solid each gas cell and its six faces are.
//
// index.html carries every compiled module as base64 (WASM_B64). After changing this file run wasm/build.sh.
// No C library: only what WebAssembly does natively (sqrt, abs, floor, min, max).

#define EXPORT(name) __attribute__((export_name(name)))
#define sqrt __builtin_sqrt
#define fabs __builtin_fabs
#define floor __builtin_floor

static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double hyp2(double a, double b) { return sqrt(a * a + b * b); }

#define NMAX 9000 // gas cells (GAS_MAX in index.html)
EXPORT("nmax") int nmax(void) { return NMAX; }

// =====================================================================
//  Geometry: parts as distance shapes, sorted into buckets
// =====================================================================
#define GH 0.1           // gas cell size (metres)
#define HB 2.0           // hash buckets per metre
#define OMAX 4096        // parts
#define OF 20            // numbers per part: see the layout below
#define HSIZE 32768      // hash slots (a power of two)
#define EMAX 262144      // bucket entries, all buckets together

// Each part, as index.html fills it in:
//  [0] shape   [1..3] position   [4..7] rotation (quaternion x y z w)   [8..10] half-extent of its world box
//  [11] moving (1 or 0)   [12..19] the shape's own numbers:
//    0 box: hx hy hz          1 ball: r          2 cylinder: r hh
//    3 hollow cylinder: R hh wall form(0 tube, 1 open cup, 2 closed)
//    4 tube bend: bend-radius tube-radius wall cos(half-angle) sin(half-angle)
//    5 wedge: w h half-depth  6 nothing solid
static double own[OMAX * OF];
static int hkey[HSIZE], hcount[HSIZE], hstart[HSIZE], hfill[HSIZE];
static unsigned char hlive[HSIZE];
static int entries[EMAX];
static int lastOwner = -1;

EXPORT("owners") double *p_owners(void) { return own; }
EXPORT("omax") int omax(void) { return OMAX; }
EXPORT("ofields") int ofields(void) { return OF; }
EXPORT("lastOwner") int last_owner(void) { return lastOwner; }

static inline int bk(double v) { return (int)floor(v * HB); }
static inline int bucketKey(int i, int j, int k) { return (i + 512) * 1048576 + (j + 512) * 1024 + (k + 512); }
static inline unsigned slotHash(int key) { return ((unsigned)key * 2654435761u) >> 17; } // 15 bits

static int findSlot(int key) {
  unsigned s = slotHash(key) & (HSIZE - 1);
  while (hkey[s] != -1 && hkey[s] != key) s = (s + 1) & (HSIZE - 1);
  return (int)s;
}
static inline int lookup(int key) {
  int s = findSlot(key);
  return hkey[s] == key ? s : -1;
}

// a part's bucket range, or 0 when it covers too many buckets to sort (as in index.html: then it is skipped)
static int range(const double *o, int *b) {
  double ex = o[8] + GH, ey = o[9] + GH, ez = o[10] + GH;
  b[0] = bk(o[1] - ex); b[1] = bk(o[1] + ex);
  b[2] = bk(o[2] - ey); b[3] = bk(o[2] + ey);
  b[4] = bk(o[3] - ez); b[5] = bk(o[3] + ez);
  return (b[1] - b[0] + 1) * (b[3] - b[2] + 1) * (b[5] - b[4] + 1) <= 4000;
}

// Sort n parts into buckets. Returns 1, or 0 when they do not fit (index.html then uses its own buckets).
EXPORT("hashBuild") int hashBuild(int n) {
  if (n > OMAX) return 0;
  for (int s = 0; s < HSIZE; s++) { hkey[s] = -1; hcount[s] = 0; hlive[s] = 0; }
  int used = 0, total = 0, b[6];
  for (int p = 0; p < n; p++) {
    const double *o = own + p * OF;
    if (!range(o, b)) continue;
    for (int i = b[0]; i <= b[1]; i++) for (int j = b[2]; j <= b[3]; j++) for (int k = b[4]; k <= b[5]; k++) {
      int key = bucketKey(i, j, k), s = findSlot(key);
      if (hkey[s] == -1) { if (++used > HSIZE / 2) return 0; hkey[s] = key; }
      hcount[s]++;
      if (o[11] != 0) hlive[s] = 1;
      if (++total > EMAX) return 0;
    }
  }
  for (int s = 0, at = 0; s < HSIZE; s++) { hstart[s] = at; hfill[s] = at; at += hcount[s]; }
  for (int p = 0; p < n; p++) {
    const double *o = own + p * OF;
    if (!range(o, b)) continue;
    for (int i = b[0]; i <= b[1]; i++) for (int j = b[2]; j <= b[3]; j++) for (int k = b[4]; k <= b[5]; k++)
      entries[hfill[findSlot(bucketKey(i, j, k))]++] = p;
  }
  return 1;
}

static double sdBox(double x, double y, double z, double hx, double hy, double hz) {
  double qx = fabs(x) - hx, qy = fabs(y) - hy, qz = fabs(z) - hz;
  double ox = dmax(qx, 0), oy = dmax(qy, 0), oz = dmax(qz, 0);
  return sqrt(ox * ox + oy * oy + oz * oz) + dmin(dmax(qx, dmax(qy, qz)), 0);
}
static double sdCyl(double x, double y, double z, double rad, double hh) {
  double dx = hyp2(x, z) - rad, dy = fabs(y) - hh;
  return dmin(dmax(dx, dy), 0) + hyp2(dmax(dx, 0), dmax(dy, 0));
}
static double sdShell(double x, double y, double z, double R, double hh, double w, int form) {
  double outer = sdCyl(x, y, z, R, hh), inner;
  if (form == 0) inner = sdCyl(x, y, z, R - w, hh + 1);
  else if (form == 1) inner = sdCyl(x, y - (w + 1) / 2, z, R - w, (2 * hh + 1 - w) / 2);
  else inner = sdCyl(x, y, z, R - w, dmax(0.001, hh - w));
  return dmax(outer, -inner);
}
static double sdBend(double x, double y, double z, const double *s) {
  double R = s[0], rr = s[1], w = s[2], ca = s[3], sa = s[4];
  double px = x + R, rho0 = hyp2(px, z);
  // inside the swept angle (|atan2(z, px)| <= half-angle): measured across the tube
  if (px >= ca * rho0) { double qq = rho0 - R; return fabs(hyp2(qq, y) - (rr - w / 2)) - w / 2; }
  double sg = z > 0 ? 1 : z < 0 ? -1 : (px < 0 ? 1 : 0);
  double ex = ca, ez = sg * sa, tx = -sg * ez, tz = sg * ex;
  double along = px * tx + z * tz, rho = px * ex + z * ez - R;
  double ds = dmax(fabs(hyp2(rho, y) - (rr - w / 2)) - w / 2, 0);
  return hyp2(dmax(along, 0), ds);
}
static double sdf(const double *o, double x, double y, double z) {
  const double *s = o + 12;
  switch ((int)o[0]) {
    case 0: return sdBox(x, y, z, s[0], s[1], s[2]);
    case 1: return sqrt(x * x + y * y + z * z) - s[0];
    case 2: return sdCyl(x, y, z, s[0], s[1]);
    case 3: return sdShell(x, y, z, s[0], s[1], s[2], (int)s[3]);
    case 4: return sdBend(x, y, z, s);
    case 5: { double w = s[0], h = s[1], nl = hyp2(w, h); return dmax(sdBox(x, y, z, w / 2, h / 2, s[2]), (x * h + y * w) / nl); }
  }
  return 1;
}

// How solid the point is (0 open .. 1 solid); lastOwner is the part that fills most of it, or -1
EXPORT("solid") double solid(double px0, double py0, double pz0) {
  double phi = 0.5 - py0 / GH;
  phi = phi < 0 ? 0 : phi > 1 ? 1 : phi; // the floor
  double best = 0;
  lastOwner = -1;
  int s = lookup(bucketKey(bk(px0), bk(py0), bk(pz0)));
  if (s >= 0) {
    for (int e = hstart[s], end = hstart[s] + hcount[s]; e < end; e++) {
      int idx = entries[e];
      const double *o = own + idx * OF;
      double px = px0 - o[1], py = py0 - o[2], pz = pz0 - o[3];
      if (fabs(px) > o[8] + GH || fabs(py) > o[9] + GH || fabs(pz) > o[10] + GH) continue;
      // into the part's frame with the inverse rotation
      double qx = -o[4], qy = -o[5], qz = -o[6], qw = o[7];
      double tx = 2 * (qy * pz - qz * py), ty = 2 * (qz * px - qx * pz), tz = 2 * (qx * py - qy * px);
      double lx = px + qw * tx + (qy * tz - qz * ty), ly = py + qw * ty + (qz * tx - qx * tz), lz = pz + qw * tz + (qx * ty - qy * tx);
      double c = 0.5 - sdf(o, lx, ly, lz) / GH;
      if (c <= 0) continue;
      double cc = c > 1 ? 1 : c;
      phi += cc;
      if (cc > best) { best = cc; lastOwner = idx; }
    }
  }
  return phi > 1 ? 1 : phi;
}

static inline int liveAt(double x, double y, double z) {
  int s = lookup(bucketKey(bk(x), bk(y), bk(z)));
  return s >= 0 && hlive[s];
}
EXPORT("bucketLive") int bucketLive(double x, double y, double z) { return liveAt(x, y, z); }

// Does anything near cell (i, j, k) move? (faces on a bucket's edge are measured in the next bucket)
static int cellLive(int i, int j, int k) {
  double x = (i + 0.5) * GH, y = (j + 0.5) * GH, zz = (k + 0.5) * GH, h = GH * 0.5;
  if (liveAt(x, y, zz)) return 1;
  if (bk(x + h) != bk(x) && liveAt(x + h, y, zz)) return 1;
  if (bk(x - h) != bk(x) && liveAt(x - h, y, zz)) return 1;
  if (bk(y + h) != bk(y) && liveAt(x, y + h, zz)) return 1;
  if (bk(y - h) != bk(y) && liveAt(x, y - h, zz)) return 1;
  if (bk(zz + h) != bk(zz) && liveAt(x, y, zz + h)) return 1;
  if (bk(zz - h) != bk(zz) && liveAt(x, y, zz - h)) return 1;
  return 0;
}

// Measuring cells in a batch: mi mj mk hold each cell's grid position and mflag whether it must be measured
// anyway (1) or only if something near it moves (0). On return mflag says whether it was measured, and for
// those: mph its solidity, mown the part filling it, mfo the openness of its six faces. mnb holds each cell's
// six neighbours as places in the batch (-1 for none): two neighbours share a face, and the point measured on
// it is exactly the same number for both, so a face the neighbour has already measured is copied.
static int mi[NMAX], mj[NMAX], mk[NMAX], mflag[NMAX], mown[NMAX], mnb[NMAX * 6];
static double mph[NMAX], mfo[NMAX * 6];
EXPORT("mi") int *p_mi(void) { return mi; }
EXPORT("mj") int *p_mj(void) { return mj; }
EXPORT("mk") int *p_mk(void) { return mk; }
EXPORT("mflag") int *p_mflag(void) { return mflag; }
EXPORT("mown") int *p_mown(void) { return mown; }
EXPORT("mnb") int *p_mnb(void) { return mnb; }
EXPORT("mph") double *p_mph(void) { return mph; }
EXPORT("mfo") double *p_mfo(void) { return mfo; }

static const int DXs[6] = {1, -1, 0, 0, 0, 0}, DYs[6] = {0, 0, 1, -1, 0, 0}, DZs[6] = {0, 0, 0, 0, 1, -1};
EXPORT("measure") int measure(int n) {
  int count = 0;
  for (int a = 0; a < n; a++) {
    int i = mi[a], j = mj[a], k = mk[a];
    if (!mflag[a] && !cellLive(i, j, k)) continue;
    mflag[a] = 1; count++;
    mph[a] = solid((i + 0.5) * GH, (j + 0.5) * GH, (k + 0.5) * GH);
    mown[a] = lastOwner;
    for (int dd = 0; dd < 6; dd++) {
      int b = mnb[a * 6 + dd];
      if (b >= 0 && b < a && mflag[b]) { mfo[a * 6 + dd] = mfo[b * 6 + (dd ^ 1)]; continue; }
      mfo[a * 6 + dd] = 1 - solid((i + 0.5 + DXs[dd] * 0.5) * GH, (j + 0.5 + DYs[dd] * 0.5) * GH, (k + 0.5 + DZs[dd] * 0.5) * GH);
    }
  }
  return count;
}
