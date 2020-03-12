#define _GNU_SOURCE
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <error.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include "ogl_sb.h"

#define MDIR "movie_XXXXXX"

static char moviedir[] = MDIR;
static int frame = 0;
static bool recording_movie = false;
static uint64_t fps;

void
pace (uint64_t *interval)
{
  if (recording_movie) {
    fps = nearbyint (1e6 / *interval);
    return;
  }
  
  static uint64_t last;
  struct timeval tv;
  gettimeofday (&tv, NULL);
  uint64_t now = tv.tv_usec + (uint64_t)tv.tv_sec * 1000000;
  if (*interval == 0) {
    last = now;
    return;
  }
  uint64_t next = last + *interval;
  static int pass;
  pass++;
  if (now > next) {
    static int count = 0;
    if (0)
      if (++count < 100)
        printf ("%d behind\n", pass);
  }
  while (now < next) {
    uint64_t rem = next - now;
    struct timespec ts = {rem / 1000000, rem % 1000000 * 1000};
    nanosleep (&ts, NULL);
    gettimeofday (&tv, NULL);
    now = tv.tv_usec + (uint64_t)tv.tv_sec * 1000000;
  }
  last = next;
}

void
progress (int *pos, int *size)
{
  static time_t lasttime = 0;
  time_t now = time (0);
  if (now == lasttime)
    return;
  lasttime = now;
  printf ("%5.1f %%\r", (double)*pos / *size * 100);
  fflush (stdout);
}

/* The movie code is in this file to give the pace routine convenient
   access to the "recording_movie" flag */

void
init_movie (void)
{
  char mdir[] = MDIR;
  if (mkdtemp (mdir) == NULL) {
    char *cwd = get_current_dir_name ();
    error (0, errno, "Sorry, can't make a directory for the movie in %s", cwd);
    free (cwd);
    return;
  }
  strcpy (moviedir, mdir);
  printf ("movie will be named movie.avi\n");
  frame = 0;
  recording_movie = true;
}

void
movie_frame (long *fildes)
{
  char *filename;
  if (!recording_movie)
    return;
  if (asprintf (&filename, "%s/frame%03d.png", moviedir, frame++) == -1) exit (1);
  print_window (fildes, filename, PW_MOVIE);
}

void
end_movie (int *is_movie)
{
  *is_movie = recording_movie ? 1 : 0;
  if (!recording_movie)
    return;
  
  recording_movie = false;
  
  char *cmd;
  if (asprintf (&cmd, "ffmpeg -i %s/frame%%03d.png -r %lu  -pix_fmt bgr24 -vcodec zlib movie.avi; rm -rf %s",
                moviedir, fps, moviedir) == -1) exit (1);
  printf ("%s\n", cmd);
  if (system (cmd));
  free (cmd);
}

void
cprint (double *x)
{
  printf ("cprint: %.17f\n", *x);
}

#if 0
#include <stdio.h>

int
main (void)
{
  pace (0);
  for (int i = 0; i < 100; i++) {
    uint64_t interval = 100000;
    pace (&interval);
    printf ("tick\n");
  }
  return 0;
}
#endif
