#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

#include <dlfcn.h>
#include <sys/wait.h>

#include"nType.h"
#include"num_in.h"
#include"num_out.h"
#include"VandP.h"
#include"dynamic_cs.h"
#include"rootDir.h"

/* Number of random phase-space points for 2->3 tests. */
#ifndef N_POINTS
#define N_POINTS 100
#endif

/* Absolute tolerance for 4-momentum conservation / on-shell tests (GeV). */
#ifndef MOM_TOL
#define MOM_TOL 1.0e-6
#endif

/*
 * amp2 for a 2->3 process:
 *  - generates CalcHEP squared-matrix-element code for `proc`
 *  - samples `nPoints` hierarchical 2->3 phase-space points
 *  - evaluates |M|^2 at each point
 *  - writes one line per point to `outFile`:
 *      E1 px1 py1 pz1  E2 px2 py2 pz2  E3 px3 py3 pz3  E4 px4 py4 pz4  E5 px5 py5 pz5  |M|^2
 *
 * Phase space follows the recursive CalcHEP construction:
 *   CM -> particle_3 + intermediate(4+5)
 *   intermediate -> particle_4 + particle_5
 * with random polar/azimuthal angles at each splitting and a random
 * intermediate invariant mass.  Points are not unweighted events; they are
 * only required to cover a wide range of kinematic configurations.
 */
double amp2(char* proc, char* remove, char* soName, double Pin,
            char* m1Name, double m1V, double hel1,
            char* m2Name, double m2V, double hel2,
            char* m3Name, double m3V,
            char* m4Name, double m4V,
            char* m5Name, double m5V,
            char* otherParamName1, double otherParamValue1,
            char* otherParamName2, double otherParamValue2,
            char* outFile, int nPoints);

/* ---------- random numbers in (0,1) and helpers ---------- */

static double rnd_unit(void)
{
  /* Open interval (0,1) to avoid exact endpoints of singular regions. */
  return (rand() + 1.0) / (RAND_MAX + 2.0);
}

static double rnd_uniform(double a, double b)
{
  return a + (b - a) * rnd_unit();
}

/* Mix of uniform and endpoint-biased samples so soft/collinear and hard
 * regions both appear among a modest number of points. */
static double rnd_mass2(double smin, double smax)
{
  double u = rnd_unit();
  double r = rnd_unit();
  if(u < 0.34)
    return smin + (smax - smin) * r;                 /* uniform in m^2 */
  if(u < 0.67)
    return smin + (smax - smin) * r * r;             /* bias soft intermediate */
  return smin + (smax - smin) * (1.0 - r * r);       /* bias hard intermediate */
}

static double kallen_sqrt(double x, double y, double z)
{
  double a = x - y - z;
  double lam = a * a - 4.0 * y * z;
  if(lam <= 0.0) return 0.0;
  return sqrt(lam);
}

/* Isotropic two-body decay of a particle of mass M at rest. */
static void two_body_rest(double M, double ma, double mb,
                          double cos_th, double phi,
                          double pa[4], double pb[4])
{
  double p = kallen_sqrt(M * M, ma * ma, mb * mb) / (2.0 * M);
  double Ea = sqrt(p * p + ma * ma);
  double Eb = sqrt(p * p + mb * mb);
  double sin_th, cphi, sphi;
  double ct = cos_th;

  if(ct > 1.0) ct = 1.0;
  if(ct < -1.0) ct = -1.0;
  sin_th = sqrt(fmax(0.0, (1.0 - ct) * (1.0 + ct)));
  cphi = cos(phi);
  sphi = sin(phi);

  pa[0] = Ea;
  pa[1] = p * sin_th * cphi;
  pa[2] = p * sin_th * sphi;
  pa[3] = p * ct;

  pb[0] = Eb;
  pb[1] = -pa[1];
  pb[2] = -pa[2];
  pb[3] = -pa[3];
}

/* Boost p_rest from the rest frame of parent into the lab frame. */
static void boost_from_rest(const double parent[4], const double p_rest[4],
                            double p_lab[4])
{
  double M2 = parent[0] * parent[0]
            - parent[1] * parent[1]
            - parent[2] * parent[2]
            - parent[3] * parent[3];
  double M, bx, by, bz, b2, g, bp, coeff;

  if(M2 <= 0.0) {
    p_lab[0] = p_rest[0];
    p_lab[1] = p_rest[1];
    p_lab[2] = p_rest[2];
    p_lab[3] = p_rest[3];
    return;
  }
  M = sqrt(M2);
  if(parent[0] <= 0.0 || M / parent[0] > 1.0 - 1e-15) {
    /* Parent essentially at rest. */
    p_lab[0] = p_rest[0];
    p_lab[1] = p_rest[1];
    p_lab[2] = p_rest[2];
    p_lab[3] = p_rest[3];
    return;
  }

  bx = parent[1] / parent[0];
  by = parent[2] / parent[0];
  bz = parent[3] / parent[0];
  b2 = bx * bx + by * by + bz * bz;
  if(b2 < 1e-30) {
    p_lab[0] = p_rest[0];
    p_lab[1] = p_rest[1];
    p_lab[2] = p_rest[2];
    p_lab[3] = p_rest[3];
    return;
  }

  g = parent[0] / M;
  bp = bx * p_rest[1] + by * p_rest[2] + bz * p_rest[3];
  coeff = ((g - 1.0) / b2) * bp + g * p_rest[0];

  p_lab[0] = g * (p_rest[0] + bp);
  p_lab[1] = p_rest[1] + coeff * bx;
  p_lab[2] = p_rest[2] + coeff * by;
  p_lab[3] = p_rest[3] + coeff * bz;
}

static double mass2(const double p[4])
{
  return p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3];
}

/*
 * Hierarchical 2->3 kinematics (CalcHEP-style recursive splitting).
 *
 * In the CM frame of the 2->3 process:
 *   1) sample intermediate mass m45 of the (4,5) cluster
 *   2) decay CM -> particle 3 + intermediate(4+5)  with random (cos θ, φ)
 *   3) decay intermediate -> 4 + 5 in its rest frame with random (cos θ*, φ*)
 *      and boost back to the CM / lab frame
 *
 * Incoming beams are along ±z with beam momentum Pin (same convention as main_22.c).
 * Returns 0 on success, nonzero if kinematics are impossible.
 */
static int fill_momenta_23(double Pin,
                           double m1, double m2, double m3, double m4, double m5,
                           REAL *p1, REAL *p2, REAL *p3, REAL *p4, REAL *p5)
{
  double E1, E2, sqrt_S, S;
  double m45_min, m45_max, m45, m45_2;
  double cos1, phi1, cos2, phi2;
  double parent[4], p3_cm[4], pint_cm[4], p4_rf[4], p5_rf[4], p4_lab[4], p5_lab[4];
  int i;

  E1 = sqrt(Pin * Pin + m1 * m1);
  E2 = sqrt(Pin * Pin + m2 * m2);
  sqrt_S = E1 + E2;
  S = sqrt_S * sqrt_S;

  m45_min = m4 + m5;
  m45_max = sqrt_S - m3;
  if(m45_max <= m45_min) return 1;

  /* Keep a tiny buffer away from exact kinematic endpoints (singular / zero measure). */
  {
    double span = m45_max - m45_min;
    double eps = fmax(1e-8 * sqrt_S, 1e-12 * span);
    if(span > 2 * eps) {
      m45_min += eps;
      m45_max -= eps;
    }
  }

  m45_2 = rnd_mass2(m45_min * m45_min, m45_max * m45_max);
  m45 = sqrt(m45_2);

  cos1 = rnd_uniform(-1.0, 1.0);
  phi1 = rnd_uniform(0.0, 2.0 * M_PI);
  cos2 = rnd_uniform(-1.0, 1.0);
  phi2 = rnd_uniform(0.0, 2.0 * M_PI);

  /* Incoming momenta in the lab / CM frame. */
  for(i = 0; i < 4; i++) {
    p1[i] = 0;
    p2[i] = 0;
    p3[i] = 0;
    p4[i] = 0;
    p5[i] = 0;
  }
  p1[0] = E1; p1[3] = Pin;
  p2[0] = E2; p2[3] = -Pin;

  /* First splitting in the total CM (at rest): CM -> 3 + (45). */
  parent[0] = sqrt_S;
  parent[1] = 0;
  parent[2] = 0;
  parent[3] = 0;
  two_body_rest(sqrt_S, m3, m45, cos1, phi1, p3_cm, pint_cm);

  /* Second splitting in the (45) rest frame, then boost to CM. */
  two_body_rest(m45, m4, m5, cos2, phi2, p4_rf, p5_rf);
  boost_from_rest(pint_cm, p4_rf, p4_lab);
  boost_from_rest(pint_cm, p5_rf, p5_lab);

  for(i = 0; i < 4; i++) {
    p3[i] = p3_cm[i];
    p4[i] = p4_lab[i];
    p5[i] = p5_lab[i];
  }

  return 0;
}

/* Return 0 if 4-momentum is conserved and particles are on-shell within MOM_TOL. */
static int check_momentum(const REAL *p1, const REAL *p2,
                          const REAL *p3, const REAL *p4, const REAL *p5,
                          double m1, double m2, double m3, double m4, double m5,
                          int ipoint)
{
  double dE  = (double)(p1[0] + p2[0] - p3[0] - p4[0] - p5[0]);
  double dpx = (double)(p1[1] + p2[1] - p3[1] - p4[1] - p5[1]);
  double dpy = (double)(p1[2] + p2[2] - p3[2] - p4[2] - p5[2]);
  double dpz = (double)(p1[3] + p2[3] - p3[3] - p4[3] - p5[3]);
  double m2err[5];
  const REAL *ps[5];
  double ms[5];
  int i, bad = 0;

  if(fabs(dE) > MOM_TOL || fabs(dpx) > MOM_TOL ||
     fabs(dpy) > MOM_TOL || fabs(dpz) > MOM_TOL) {
    fprintf(stderr,
            "Point %d: momentum not conserved: dE=%g dpx=%g dpy=%g dpz=%g\n",
            ipoint, dE, dpx, dpy, dpz);
    bad = 1;
  }

  ps[0] = p1; ps[1] = p2; ps[2] = p3; ps[3] = p4; ps[4] = p5;
  ms[0] = m1; ms[1] = m2; ms[2] = m3; ms[3] = m4; ms[4] = m5;
  for(i = 0; i < 5; i++) {
    double p[4] = { (double)ps[i][0], (double)ps[i][1],
                    (double)ps[i][2], (double)ps[i][3] };
    m2err[i] = mass2(p) - ms[i] * ms[i];
    /* Relative-friendly absolute cut; massless legs can have tiny negative m^2 from roundoff. */
    if(fabs(m2err[i]) > fmax(MOM_TOL, 1e-8 * fmax(1.0, p[0] * p[0]))) {
      fprintf(stderr, "Point %d: particle %d off-shell: p^2 - m^2 = %g (m=%g)\n",
              ipoint, i + 1, m2err[i], ms[i]);
      bad = 1;
    }
  }
  return bad;
}

int main(void)
{
  /* e+ e- -> gamma gamma gamma  (tree-level QED). */
  amp2("e,E->A,A,A", "", "eeAAA",
       250.0,
       "Me", 0.0005, 0,
       "Me", 0.0005, 0,
       "", 0,
       "", 0,
       "", 0,
       "", 0, "", 0,
       "eeAAA_points.dat", N_POINTS);
  return 0;
}

double amp2(char* proc, char* remove, char* soName, double Pin,
            char* m1Name, double m1V, double hel1,
            char* m2Name, double m2V, double hel2,
            char* m3Name, double m3V,
            char* m4Name, double m4V,
            char* m5Name, double m5V,
            char* otherParamName1, double otherParamValue1,
            char* otherParamName2, double otherParamValue2,
            char* outFile, int nPoints)
{
  int err;
  REAL pvec[20];
  REAL *p1 = pvec, *p2 = pvec + 4, *p3 = pvec + 8, *p4 = pvec + 12, *p5 = pvec + 16;
  REAL m1, m2, m3, m4, m5, m[5];
  double GG;
  numout *cc;
  FILE *fp;
  int i, nsub, nin, nout;
  int nOK = 0, nFail = 0;
  unsigned int seed = 1;

  printf("Calculating CH |M|^2 for 2->3 process %s (remove '%s').\n", proc, remove);
  printf("Pin = %f, nPoints = %d, outFile = %s\n", Pin, nPoints, outFile);

  setModel("models", 1);

  GG = 1.238;
  printf("EE=%E\n", findValW("EE"));
  printf("GG=%E\n", GG);
  printf("MW=%E\n", findValW("MW"));
  assignValW("wW", 0);
  printf("wW=%E\n", findValW("wW"));
  printf("wZ=%E\n", findValW("wZ"));
  printf("SW=%E\n", findValW("SW"));
  printf("Mh=%E\n", findValW("Mh"));
  printf("wh=%E\n", findValW("wh"));

  if(strlen(otherParamName1) > 0) {
    assignValW(otherParamName1, otherParamValue1);
    printf("%s=%E\n", otherParamName1, findValW(otherParamName1));
  }
  if(strlen(otherParamName2) > 0) {
    assignValW(otherParamName2, otherParamValue2);
    printf("%s=%E\n", otherParamName2, findValW(otherParamName2));
  }

  if(strlen(m1Name) > 0) {
    assignValW(m1Name, m1V);
    printf("1: %s=%E\n", m1Name, findValW(m1Name));
  } else printf("1: 0\n");
  if(strlen(m2Name) > 0) {
    assignValW(m2Name, m2V);
    printf("2: %s=%E\n", m2Name, findValW(m2Name));
  } else printf("2: 0\n");
  if(strlen(m3Name) > 0) {
    assignValW(m3Name, m3V);
    printf("3: %s=%E\n", m3Name, findValW(m3Name));
  } else printf("3: 0\n");
  if(strlen(m4Name) > 0) {
    assignValW(m4Name, m4V);
    printf("4: %s=%E\n", m4Name, findValW(m4Name));
  } else printf("4: 0\n");
  if(strlen(m5Name) > 0) {
    assignValW(m5Name, m5V);
    printf("5: %s=%E\n", m5Name, findValW(m5Name));
  } else printf("5: 0\n");

  err = calcMainFunc();
  if(err) {
    printf("Can not calculate constrained parameter %s\n", varNames[err]);
    return err;
  }

  cc = getMEcode(0, 1, proc, remove, "", soName);
  if(!cc) {
    printf("getMEcode failed for process %s\n", proc);
    return 1;
  }
  err = passParameters(cc);
  if(err) {
    printf("Can not calculate constrained parameter %s\n",
           cc->interface->varName[err]);
    return err;
  }

  procInfo1(cc, &nsub, &nin, &nout);
  printf("Process library: nsub=%d nin=%d nout=%d\n", nsub, nin, nout);
  if(nin != 2 || nout != 3) {
    printf("Error: expected 2->3 process, got %d->%d\n", nin, nout);
    return 1;
  }

  procInfo2(cc, 1, NULL, m);
  m1 = m[0]; m2 = m[1]; m3 = m[2]; m4 = m[3]; m5 = m[4];
  printf("Masses: %g %g -> %g %g %g\n",
         (double)m1, (double)m2, (double)m3, (double)m4, (double)m5);

  /* Polarized massless initial states (photons / gluons) use Helicity[], as in main_22. */
  if(strlen(m1Name) == 0) {
    Helicity[0] = hel1;
    printf("Hel1 = %Lf\n", Helicity[0]);
  }
  if(strlen(m2Name) == 0) {
    Helicity[1] = hel2;
    printf("Hel2 = %Lf\n", Helicity[1]);
  }

  fp = fopen(outFile, "w");
  if(!fp) {
    perror(outFile);
    return 1;
  }
  fprintf(fp,
          "# process: %s\n"
          "# Pin=%g  seed=%u  nPoints=%d\n"
          "# columns: E1 px1 py1 pz1  E2 px2 py2 pz2  E3 px3 py3 pz3  "
          "E4 px4 py4 pz4  E5 px5 py5 pz5  |M|^2\n",
          proc, Pin, seed, nPoints);

  srand(seed);

  for(i = 0; i < nPoints; i++) {
    int sqerr = 0;
    double amp;
    int tries;

    /* Retry a few times if a rare kinematic sample fails. */
    for(tries = 0; tries < 20; tries++) {
      if(fill_momenta_23(Pin, m1, m2, m3, m4, m5, p1, p2, p3, p4, p5) == 0)
        break;
    }
    if(tries >= 20) {
      fprintf(stderr, "Point %d: failed to generate kinematics\n", i);
      nFail++;
      continue;
    }

    if(check_momentum(p1, p2, p3, p4, p5, m1, m2, m3, m4, m5, i)) {
      nFail++;
      continue;
    }

    amp = cc->interface->sqme(1, GG, pvec, NULL, &sqerr);
    if(sqerr) {
      fprintf(stderr, "Point %d: sqme err_code=%d\n", i, sqerr);
      nFail++;
      continue;
    }

    fprintf(fp,
            "%.15E %.15E %.15E %.15E "
            "%.15E %.15E %.15E %.15E "
            "%.15E %.15E %.15E %.15E "
            "%.15E %.15E %.15E %.15E "
            "%.15E %.15E %.15E %.15E "
            "%.15E\n",
            (double)p1[0], (double)p1[1], (double)p1[2], (double)p1[3],
            (double)p2[0], (double)p2[1], (double)p2[2], (double)p2[3],
            (double)p3[0], (double)p3[1], (double)p3[2], (double)p3[3],
            (double)p4[0], (double)p4[1], (double)p4[2], (double)p4[3],
            (double)p5[0], (double)p5[1], (double)p5[2], (double)p5[3],
            amp);
    nOK++;

    if((i + 1) % 20 == 0 || i == 0)
      printf("point %3d: |M|^2 = %.6E\n", i, amp);
  }

  fclose(fp);
  printf("Wrote %d / %d points to %s  (failed/skipped: %d)\n",
         nOK, nPoints, outFile, nFail);
  return (nOK > 0) ? 0 : 1;
}
