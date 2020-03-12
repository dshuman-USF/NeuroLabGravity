#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

typedef struct
{
  double time, min, mean, max, data;
} Plot;
Plot p;

int
main (void)
{
#define S(x)  __asm__ ("bswap %0" : "=r" (x) : "0" (x))
//  printf ("%08x %g\n", *(int *)&x, x);
//  S(x);
//  printf ("%08x %g\n", *(int *)&x, x);
//  S(x);
//  printf ("%08x %g\n", *(int *)&x, x);
  
  printf ("%g\n", p.time);
  return 0;
}
