#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <error.h>
#include <errno.h>
#include <byteswap.h>
#include <glob.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

typedef struct
{
  int cell_1;
  int cell_2;
  int npts;
  int inc;
  float endtime;
  int ncfile;
  int velcnt;
} SurFileHdr;

static char *experiment;
static int trial;
static int particle1;
static int particle2;

static void
swap_hdr (SurFileHdr *h)
{
#define S(x)  __asm__ ("bswap %0" : "=r" (x) : "0" (x))
  S(h->cell_1);
  S(h->cell_2);
  S(h->npts);
  S(h->inc);
  S(h->endtime);
  S(h->ncfile);
  S(h->velcnt);
#undef S
}

static void
print_hdr (SurFileHdr *h)
{
#define P(x) printf (#x ": %d\n", h->x)
#define PF(x) printf (#x ": %g\n", h->x)
  P (cell_1);
  P (cell_2);
  P (npts);
  P (inc);
  PF (endtime);
  P (ncfile);
  P (velcnt);
#undef P  
}

static double
mixpmf (int k, double p1, double p2, int n1, int velcnt)
{
  double sum = 0;
  for (int i = 0; i <= k; i++)
    sum += gsl_ran_binomial_pdf (i, p1, n1) * gsl_ran_binomial_pdf (k - i, p2, velcnt - n1);
  return sum;
}

static void
lineidx (FILE *f, int n, int *code, int *time)
{
  static int linelen = 14;
  long offset = n * linelen;

  if (fseek (f, offset, n < 0 ? SEEK_END : SEEK_SET))
    error_at_line (1, errno, __FILE__, __LINE__, "fseek %ld", offset);
  char s[14];
  if (!fgets (s, sizeof s, f))
    error_at_line (1, errno, __FILE__, __LINE__, "fgets");
  *time = atoi (s + 5);
  s[5] = 0;
  *code = atoi (s);
}

static double
get_timadj (void)
{
  glob_t globbuf;

  if (glob("*_vertrial?.gdt", 0, NULL, &globbuf))
    error (1, errno, "glob");
  if (globbuf.gl_pathc != 1)
    error (1, 0, "error: %ld gdt files: there must be only one", (long)globbuf.gl_pathc);
  FILE *f = fopen (globbuf.gl_pathv[0], "r");
  int code, time;
  lineidx (f, 2, &code, &time);
  if (code != 21)
    error (1, 0, "third line of %s is not code 21", globbuf.gl_pathv[0]);
  return time / 2. - .5;
}

static SurFileHdr h;
static float **v;
static char *argv1;

typedef struct
{
  double mean1;
  double s_1;
  double mean2;
  double s_2;
  double p;
  double df;
  int n1;
  int n2;
} Stats;

#define TRUNC 10

static Stats
compare_means (double stop, double start)
{
  //  printf ("stop: %g, start: %g, velcnt %d\n", stop, start, h.velcnt);
  int n1 = 0;
  double mean1 = 0;
  double M2_1 = 0;

  int n2 = 0;
  double mean2 = 0;
  double M2_2 = 0;
  
  for (int i = 0; i < h.velcnt - TRUNC; i++) {
    double x = v[i][0];
    if (i + 1 <= stop) {
      n1++;
      double delta = x - mean1;
      mean1 += delta / n1;
      M2_1 += delta * (x - mean1);
    }
    else if (i >= start) {
      n2++;
      double delta = x - mean2;
      mean2 += delta / n2;
      M2_2 += delta * (x - mean2);
    }
  }
  double s2_1 =  M2_1 / (n1 - 1);
  double s2_2 =  M2_2 / (n2 - 1);
  double s12 = s2_1 / n1 + s2_2 / n2;
  double t = (mean2 - mean1) / sqrt (s12);
  //  printf ("t: %g\n", t); fflush (stdout);
#define SQ(x) ((x)*(x))
  double df = SQ(s12) / (SQ(s2_1 / n1) / (n1 - 1) + SQ(s2_2 / n2) / (n2 - 1));
  double p =  gsl_cdf_tdist_Q (fabs (t), df);
  Stats s;
  s.mean1 = mean1;
  s.s_1 = sqrt(s2_1) / sqrt (n1);
  s.mean2 = mean2;
  s.s_2 = sqrt (s2_2) / sqrt (n1);
  s.p = p;
  s.df = df;
  s.n1 = n1;
  s.n2 = n2;
  return s;
}

static void
print_stats (Stats s)
{
  printf ("%10s %d %-9s mean1: %10.7f+=%10.7f, mean2: %10.7f+-%10.7f, p: %.9f, df: %7.3f, n1: %d, n2: %d\n",
          experiment, trial, argv1, s.mean1, s.s_1, s.mean2, s.s_2, s.p, s.df, s.n1, s.n2);
}

static void
parse_dir (void)
{
  char *dir = strdup (getenv ("PWD"));
  char *tmp = 0, *try = 0;
  char *p = dir + 1;
  while ((p = strchr (p, '/')) != NULL) {
    experiment = tmp;
    tmp = try;
    try = ++p;
  }
  trial = atoi (try + strlen ("trial"));
  tmp[-1] = 0;
  particle1 = atoi (argv1);
  p = strchr (argv1, '-');
  particle2 = atoi (p + 1);
  //  printf ("experiment: %s, trial: %d, particles: %d %d\n", experiment, trial, particle1, particle2);
}

static double *
quantile (int dn0)
{
  double *q = calloc (h.velcnt, sizeof *q);
  for (int vn = 0; vn < h.velcnt; vn++) {
    int count = 0;
    for (int dn = 0; dn <= h.ncfile; dn++)
      if (v[vn][dn] < v[vn][dn0])
        count++;
    q[vn] = (double)count / h.ncfile;
  }
  return q;
}

static double *
velocity (int dn0, double *lo, double *hi)
{
  double *vel = calloc (h.velcnt, sizeof *vel);
  double min = v[0][0];
  double max = v[0][0];
  for (int vn = 0; vn < h.velcnt; vn++) {
    for (int dn = 0; dn <= h.ncfile; dn++) {
      if (v[vn][dn] < min)
        min = v[vn][dn];
      if (v[vn][dn] > max)
        max = v[vn][dn];
    }
    vel[vn] = v[vn][dn0];
  }
  *lo = min;
  *hi = max;
  return vel;
}

#define BINCNT 30

static int *
vhist (double *vel, double lo, double hi)
{
  int *hist = calloc (BINCNT, sizeof *hist);
  for (int vn = 0; vn < h.velcnt; vn++) {
    int bin = (vel[vn] - lo) / (hi - lo) * BINCNT;
    if (bin >= BINCNT)
      bin = BINCNT - 1;
    if (bin < 0)
      bin = 0;
    hist[bin]++;
  }  
  return hist;
}

static int *
qhist (double *q)
{
  int *hist = calloc (BINCNT, sizeof *hist);
  for (int vn = 0; vn < h.velcnt; vn++) {
    int bin = q[vn] * BINCNT;
    if (bin == BINCNT)
      bin = BINCNT - 1;
    if (bin < 0 || bin >= BINCNT)
      error_at_line (1, 0, __FILE__, __LINE__, "quantile out of range");
    hist[bin]++;
  }  
  return hist;
}

static double
mean_sd (float *x, int xcnt, double *sd)
{
  int n = 0;
  double mean = 0;
  double M2 = 0;

  for (int i = 0; i < xcnt; i++) {
    n++;
    double delta = x[i] - mean;
    mean += delta / n;
    M2 += delta * (x[i] - mean);
  } 
  *sd = sqrt (M2 / (n - 1));
  return mean;
}

static double *
sdvel (int dn0)
{
  double *q = calloc (h.velcnt, sizeof *q);
  for (int vn = 0; vn < h.velcnt; vn++) {
    double sd;
    double mean = mean_sd (v[vn] + 1, h.ncfile, &sd);
    q[vn] = (v[vn][dn0] - mean) / sd;
  }
  return q;
}

static int *
sdvhist (double *y)
{
  int *hist = calloc (BINCNT, sizeof *hist);
  for (int vn = 0; vn < h.velcnt; vn++) {
    int bin = (y[vn] + 6) / 12 * BINCNT;
    if (bin >= BINCNT)
      bin = BINCNT - 1;
    if (bin < 0)
      bin = 0;
    hist[bin]++;
  }  
  return hist;
}

static void
plotxvel (double *y, double y0)
{
  FILE *f = fopen ("p", "w");
  if (!f) error_at_line (1, errno, __FILE__, __LINE__, "Can't open \"p\" for write");
  double slicesec = (int)(h.endtime / (h.velcnt + 1)) / 1000.;
  for (int vn = 0; vn < h.velcnt; vn++)
    if (fprintf (f, "%g %g\n", vn * slicesec, y[vn]) < 0)
      error_at_line (1, errno, __FILE__, __LINE__, "Error writing \"p\"");
  fclose (f);
  if ((f = fopen ("g", "w")) == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't open \"b\" for write");
  if (fprintf (f, "plot \"p\" w l lt -1, %g lt -1\n", y0) < 0)
    error_at_line (1, errno, __FILE__, __LINE__, "Error writing \"g\"");
  fclose (f);
  if (system ("gnuplot -persist g"));
}

static void
plothist (int *hist, double xl, double xr)
{
  FILE *f = fopen ("p", "w");
  if (!f) error_at_line (1, errno, __FILE__, __LINE__, "Can't open \"p\" for write");
  double binw = (xr - xl) / BINCNT;
  for (int i = 0; i < BINCNT; i++)
    if (fprintf (f, "%g %d\n", xl + binw / 2 + i * binw, hist[i]) < 0)
      error_at_line (1, errno, __FILE__, __LINE__, "Error writing \"p\"");
  fclose (f);
  if ((f = fopen ("g", "w")) == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't open \"b\" for write");
  if (xl == 0 && xr == 1) {
    double mean = h.velcnt / BINCNT;
    if (fprintf (f, "plot [%g:%g] \"p\" w histeps lt -1, %g lt -1\n", xl, xr, mean) < 0)
      error_at_line (1, errno, __FILE__, __LINE__, "Error writing \"g\"");
  }
  else
    if (fprintf (f, "plot [%g:%g] \"p\" w histeps lt -1\n", xl, xr) < 0)
      error_at_line (1, errno, __FILE__, __LINE__, "Error writing \"g\"");
  fclose (f);
  if (system ("gnuplot -persist g"));
}

int
main (int argc, char **argv)
{
  if (argc < 2)
    error (0, 0, "usage: %s surfile", argv[0]);
  argv1 = argv[1];
  parse_dir ();
  FILE *f = fopen (argv[1], "rb");
  if (!f)
    error (1, errno, "fopen");

  
  if (fread (&h, sizeof h, 1, f) != 1)
    error (1, errno, "Can't read surfile header");
  unsigned short swaptest = 1;
  if (*(char *)&swaptest == 1)         /* swap if little-endian */
    swap_hdr (&h);

  if (0)
    print_hdr (&h);

  int bufcnt = h.velcnt * (h.ncfile + 1);
  float *buf = malloc (bufcnt * sizeof *buf);
  v = malloc ((h.ncfile + 1) * sizeof *v);
  for (int i = 0; i < h.ncfile + 1; i++)
    v[i] = buf + h.velcnt * i;

  if (fread (buf, sizeof *buf, bufcnt, f) != bufcnt)
    error (1, errno, "Can't read velocities\n");
  fclose (f);

  unsigned int *uibuf = (unsigned int *)buf;
  if (*(char *)&swaptest == 1)         /* swap if little-endian */
    for (int i = 0; i < bufcnt; i++)
      uibuf[i] = bswap_32 (uibuf[i]);

  float *buf2 = malloc (bufcnt * sizeof *buf);
  float **v2 = malloc (h.velcnt * sizeof *v);
  for (int i = 0; i < h.velcnt; i++)
    v2[i] = buf2 + (h.ncfile + 1) * i;

  for (int dn = 0; dn <= h.ncfile; dn++)
    for (int vn = 0; vn < h.velcnt; vn++)
      v2[vn][dn] = v[dn][vn];
  free (buf);
  free (v);
  v = v2;

  if (argc > 2 && strcmp (argv[2], "ttest") == 0) {
    double stop, start;
    if (argc == 5) {
      stop  = atof (argv[3]);
      start = atof (argv[4]);
    }
    else {
      char *fname = "Ipulse-times";
      FILE *f = fopen (fname, "r");
      if (!f)
        error (1, errno, "Can't open %s", fname);
      if (fscanf (f, "%lf %lf", &stop, &start) != 2)
        error (1, errno, "%s is invalid", fname);
    }
    double timadj = get_timadj ();
    //    printf ("timadj: %g\n", timadj);
    stop  = (stop * 1000 - timadj) / (int)(h.endtime / (h.velcnt + 1));
    start = (start * 1000 - timadj) / (int)(h.endtime / (h.velcnt + 1));
    Stats s = compare_means (stop, start);
    print_stats (s);
    return 0;
  }

  if (argc == 4 && strcmp (argv[2], "quantile") == 0) {
    plotxvel (quantile (atoi (argv[3])), .5);
    return 0;
  }

  if (argc == 4 && strcmp (argv[2], "qhist") == 0) {
    plothist (qhist (quantile (atoi (argv[3]))), 0, 1);
    return 0;
  }

  if (argc == 4 && strcmp (argv[2], "sdvel") == 0) {
    plotxvel (sdvel (atoi (argv[3])), 0);
    return 0;
  }

  if (argc == 4 && strcmp (argv[2], "sdvhist") == 0) {
    plothist (sdvhist (sdvel (atoi (argv[3]))), -6, 6);
    return 0;
  }

  if (argc == 4 && strcmp (argv[2], "vhist") == 0) {
    double lo, hi;
    double *vel = velocity (atoi (argv[3]), &lo, &hi);
    plothist (vhist (vel, lo, hi), lo, hi);
    return 0;
  }

  if (argc >= 4 && strcmp (argv[2], "mvhist") == 0) {
    double lo, hi;
    double *vel = velocity (atoi (argv[3]), &lo, &hi);
    lo = -.1;
    hi = +.1;
    if (argc > 4) lo = atof (argv[4]);
    if (argc > 5) hi = atof (argv[5]);
    plothist (vhist (vel, lo, hi), lo, hi);
    return 0;
  }

  if (argc == 3 && strcmp (argv[2], "xsplot") == 0) {
    double stop, start;
    char *fname = "Ipulse-times";
    FILE *f = fopen (fname, "r");
    if (!f)
      error (1, errno, "Can't open %s", fname);
    if (fscanf (f, "%lf %lf", &stop, &start) != 2)
      error (1, errno, "%s is invalid", fname);
    fclose (f);
    double timadj = get_timadj ();
    stop  = (stop * 1000 - timadj) / (int)(h.endtime / (h.velcnt + 1));
    start = (start * 1000 - timadj) / (int)(h.endtime / (h.velcnt + 1));

    double slicesec = (int)(h.endtime / (h.velcnt + 1)) / 1000.;
    double *pos = calloc (h.ncfile + 1, sizeof *pos);
    for (int i = 0; i < h.ncfile + 1; i++)
      pos[i] = 100;
    double x[4] = {0, floor (stop), ceil (start), h.velcnt - TRUNC};
    double y[4] = {100};
    if ((f = fopen ("p", "w")) == NULL)
      error_at_line (1, errno, __FILE__, __LINE__, "Can't open file named \"p\" for write");
    fprintf (f, "%g %g %g %g %g\n", 0., 100., 100., 100., 100.);
    double maxmax = 0;
    for (int vn = 0; vn < h.velcnt - TRUNC; vn++) {
      for (int dn = 0; dn <= h.ncfile; dn++)
        pos[dn] -= v[vn][dn];
      for (int i = 1; i < 4; i++)
        if (vn + 1 == x[i])
          y[i] = pos[0];
      double max = pos[1];
      double min = pos[1];
      double sum = pos[1];
      for (int dn = 2; dn <= h.ncfile; dn++) {
        sum += pos[dn];
        if (pos[dn] > max) max = pos[dn];
        if (pos[dn] < min) min = pos[dn];
      }
      double mean = sum / h.ncfile; 
      if (max > maxmax)
        maxmax = max;
      if (pos[0] > maxmax)
        maxmax = pos[0];
      fprintf (f, "%g %g %g %g %g\n", (vn + 1) * slicesec, min, mean, max, pos[0]);
      if (vn + 1 == x[3] && pos[0] < min)
        printf ("%s/%s-grav_P/trial%d/%d-%d.ps\n",  experiment, experiment, trial, particle1, particle2);
    }
    Stats s = compare_means (stop, start);
    assert (fabs ((y[0] - y[1]) / s.n1 - s.mean1) / s.mean1 < .001);
    if (!(fabs ((y[2] - y[3]) / s.n2 - s.mean2) / s.mean2 < .001)) {
      error_at_line (0, 0, __FILE__, __LINE__, "%s trial %d pair %d-%d", experiment, trial, particle1, particle2);
      printf ("ydif/n: %g s.mean2: %g\n", (y[2] - y[3]) / s.n2, s.mean2);
      assert (fabs ((y[2] - y[3]) / s.n2 - s.mean2) / s.mean2 < .001);
    }
    double dy1 = s.s_1 * (x[1] - x[0]) / 2;
    double dy2 = s.s_2 * (x[3] - x[2]) / 2;
    for (int i = 0; i < 4; i++)
      x[i] *= slicesec;
    fprintf (f, "\n");

    fprintf (f, "\n%g %g\n%g %g\n", x[0], y[0]      , x[1], y[1]      );
    fprintf (f, "\n%g %g\n%g %g\n", x[0], y[0] + dy1, x[1], y[1] - dy1);
    fprintf (f, "\n%g %g\n%g %g\n", x[0], y[0] - dy1, x[1], y[1] + dy1);

    fprintf (f, "\n");

    fprintf (f, "\n%g %g\n%g %g\n", x[2], y[2]      , x[3], y[3]      );
    fprintf (f, "\n%g %g\n%g %g\n", x[2], y[2] + dy2, x[3], y[3] - dy2);
    fprintf (f, "\n%g %g\n%g %g\n", x[2], y[2] - dy2, x[3], y[3] + dy2);

    double x1 = (x[0] + x[1]) / 2;
    double x2 = (x[2] + x[3]) / 2;
    double y1 = (y[0] + y[1]) / 2;
    double y2 = (y[2] + y[3]) / 2;
    double xmax = x[3];

    x[2] -= x2 - x1;
    x[3] -= x2 - x1;
    y[0] -= y[3] - dy2 + y1 - y2;
    y[1] -= y[3] - dy2 + y1 - y2;
    y[2] -= y[3] - dy2;
    y[3] -= y[3] - dy2;

    for (int i = 0; i < 4; i++) x[i] += 10, y[i] += 10;

    fprintf (f, "\n");

    fprintf (f, "\n%g %g\n%g %g\n", x[0], y[0]      , x[1], y[1]      );
    fprintf (f, "\n%g %g\n%g %g\n", x[0], y[0] + dy1, x[1], y[1] - dy1);
    fprintf (f, "\n%g %g\n%g %g\n", x[0], y[0] - dy1, x[1], y[1] + dy1);

    fprintf (f, "\n");

    fprintf (f, "\n%g %g\n%g %g\n", x[2], y[2]      , x[3], y[3]      );
    fprintf (f, "\n%g %g\n%g %g\n", x[2], y[2] + dy2, x[3], y[3] - dy2);
    fprintf (f, "\n%g %g\n%g %g\n", x[2], y[2] - dy2, x[3], y[3] + dy2);

    fclose (f);
    if ((f = fopen ("g", "w")) == NULL)
      error_at_line (1, errno, __FILE__, __LINE__, "Can't open file named \"g\" for write");
    fprintf (f, "set title \"%s Trial %d Particles: %d %d\"\n", experiment, trial, particle1, particle2);
    fprintf (f, "set xlabel \"Time (s)\"\n");
    fprintf (f, "set ylabel \"Distance\"\n");

    fprintf (f, "set terminal postscript\n");
    char *oname = strdup (argv1);
    int len = strlen (oname);
    if (strcmp (oname + len - 4, ".sur") != 0)
      error_at_line (1, 0, __FILE__, __LINE__, "input file must have a .sur extension");
    strcpy (oname + len - 4, ".ps");
    fprintf (f, "set output \"%s\"\n", oname);
    free (oname);

    fprintf (f, "plot  [0:%g] [0:%g]", xmax, maxmax);
    fprintf (f, "  \"p\" u 1:2 i 0 w l lt -1 t \"\"");
    fprintf (f, ", \"p\" u 1:3 i 0 w l lt -1 t \"\"");
    fprintf (f, ", \"p\" u 1:4 i 0 w l lt -1 t \"\"");
    fprintf (f, ", \"p\" u 1:5 i 0 w l lt -1 lw 2 t \"\"");
    fprintf (f, ", \"p\" u 1:2 i 1 w l lt 1 lc rgb \"blue\" t \"\"");
    fprintf (f, ", \"p\" u 1:2 i 2 w l lt 1 lc rgb \"red\" t \"\"");
    fprintf (f, ", \"p\" u 1:2 i 3 w l lt 1 lc rgb \"blue\" t \"\"");
    fprintf (f, ", \"p\" u 1:2 i 4 w l lt 1 lc rgb \"red\" t \"\"");
    fprintf (f, "\n");
    fclose (f);

    if (system ("gnuplot g"));
    return 0;
  }

  if (argc == 3 && (strcmp (argv[2], "minvel") == 0 || strcmp (argv[2], "plot") == 0)) {
    double *sum = calloc (h.ncfile + 1, sizeof *sum);
    for (int vn = 0; vn < h.velcnt; vn++)
      for (int dn = 1; dn <= h.ncfile; dn++)
        sum[dn] += v[vn][dn];
    int dn_at_min = 1;
    int dn_at_max = 1;
    for (int dn = 2; dn <= h.ncfile; dn++) {
      if (sum[dn] < sum[dn_at_min])
        dn_at_min = dn;
      if (sum[dn] > sum[dn_at_max])
        dn_at_max = dn;
    }
    printf ("the surrogate with the smallest velocity sum is %d\n", dn_at_min);
    printf ("the surrogate with the largest velocity sum is %d\n", dn_at_max);
    f = fopen ("surstats_plot", "w");
    if (strcmp (argv[2], "plot") == 0) {
      double minsum = 0, dsum = 0, maxsum = 0;
      for (int vn = 0; vn < h.velcnt; vn++) {
        minsum += v[vn][dn_at_min];
        maxsum += v[vn][dn_at_max];
        dsum += v[vn][0];
        fprintf (f, "%g %g %g\n", minsum, dsum, maxsum);
      }
    }
    fclose (f);
    return 0;
  }

  if (argc == 3 && strcmp (argv[2], "slices") == 0) {
    double timadj = get_timadj ();
    int slicewidth = h.endtime / (h.velcnt + 1);
    for (int vn = 0; vn < h.velcnt; vn++)
      printf ("%5d%8d\n", 51, (int)nearbyint ((vn * slicewidth + timadj) * 2));
    return 0;
  }

  int maxec = 0;
  int *dn_at_max = calloc (h.velcnt, sizeof *dn_at_max);

  for (int vn = 0; vn < h.velcnt; vn++) {
    int dnm = 0;
    int maxcnt = 1;
    for (int dn = 1; dn <= h.ncfile; dn++)
      if (v[vn][dn] > v[vn][dnm]) {
        dnm = dn;
        maxcnt = 1;
      }
      else if (v[vn][dn] == v[vn][dnm])
        maxcnt++;
    dn_at_max[vn] = -1;
    if (maxcnt == 1)
      dn_at_max[vn] = v[vn][dnm] > .01 ? dnm : dnm - 2;
    if (maxcnt == 1 && dn_at_max[vn] == 0) {
      if (vn > 0) printf ("vn %d: %11.8f", vn, v[vn - 1][0]);
      printf (" %11.8f ", v[vn][0]);
      if (vn + 1 < h.velcnt) printf (" %11.8f", v[vn + 1][0]);
      printf ("\n");
    }
  }

  if (argc == 3 && strcmp (argv[2], "flats") == 0) {
    printf ("doing flats\n");
    double timadj = get_timadj ();
    int slicewidth = h.endtime / (h.velcnt + 1);
    printf ("slicewidth: %d\n", slicewidth);
    int zcnt = 0;
    for (int vn = 0; vn < h.velcnt; vn++) {
      if (dn_at_max[vn] != -2 && dn_at_max[vn] != 0 && zcnt > 1) {
        double start = ((vn - zcnt + 1) * slicewidth + timadj) / 1000;
        double stop = ((vn + 1) * slicewidth + timadj) / 1000;
        printf ("%d-slice wide data extreme from %.4f (%.0f ticks) to %.4f seconds (%.0f ticks) (%.4f seconds) v: %g, ms: %.1f\n",
                zcnt, start, start * 2000, stop, stop * 2000, stop - start, v[vn - 1][0], start * 1000 - timadj);
      }
      zcnt = dn_at_max[vn] == 0 || dn_at_max[vn] == -2 ? zcnt + 1 : 0;
      if (zcnt)
        printf ("zcnt = %d at %d\n", zcnt, vn);
    }
    return 0;
  }

  int *extreme_count = calloc (h.ncfile + 1, sizeof *extreme_count);
  for (int vn = 0; vn < h.velcnt; vn++)
    if (dn_at_max[vn] >= 0)
      extreme_count[dn_at_max[vn]]++;

  int n = 0;
  double mean = 0;
  double M2 = 0;

  int *hist = calloc (h.velcnt, sizeof *hist);
  for (int dn = 0; dn <= h.ncfile; dn++) {
    hist[extreme_count[dn]]++;
    if (extreme_count[dn] > maxec)
      maxec = extreme_count[dn];
    {
      n++;
      double delta = extreme_count[dn] - mean;
      mean += delta / n;
      M2 += delta * (extreme_count[dn] - mean);
    }
  }

  double sd =  sqrt (M2 / (n - 1));


  fprintf (stderr, "extreme count for data: %d\n", extreme_count[0]);
  fprintf (stderr, "mean of %d surrogates: %g, sd %g\n", n, mean, sd);
  double p = mean / h.velcnt;
  for (int i = 0; i <= maxec; i++)
    if (0)
      printf ("%d %d %g %g\n", i, hist[i], n * gsl_ran_binomial_pdf (i, p, h.velcnt), n * mixpmf (i, .015, .005, h.velcnt/8, h.velcnt));
    else if (hist[i])
      printf ("%d %d\n", i, hist[i]);
  return 0;
}
