/* molecule.c: run as part of xgvis ML 2/92. */
/*
 * molecule.c: find layout with attraction/repulsion.
 * dists.c: writes out distance matrix.
 * ML: 11/91.
*/

/* Includes. */
#include <sys/types.h>
#include <math.h>
#include <X11/keysym.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgvis.h"

/* Functions. */
double distance(double *, double *, int, double);
double normalize(double *, int, double);

/* Macros. */
#define DOUBLE(x) ((double)(x))
#define FLOAT(x) ((float)(x))

/* Euclidean for now. */
double
distance(double *p1, double *p2, int dims, double dlnorm)
{
  double dsum = 0.0;
  int i;
  for (i = 0; i < dims; i++) {
    dsum += pow(fabs(p1[i]-p2[i]), dlnorm);
  }
  return(pow(dsum, 1.0/dlnorm));
}

/* Euclidean for now. */
double
normalize(double *p1, int dims, double dlnorm)
{
  double dsum = 0.0;
  double dtmp;
  int i;
  for (i = 0; i < dims; i++) {
    dtmp = pow(fabs(p1[i]), dlnorm);
    dsum += dtmp;
  }
  return(pow(dsum, 1.0/dlnorm));
}
