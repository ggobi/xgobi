/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  the xgobi code and associated documentation,  provided  *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *        Copyright (c) 1997 AT&T                           *
 *                                                          *
 *  We welcome your questions and comments, and request     *
 *  that you share any modifications with us.               *
 *                                                          *
 *    Deborah F. Swayne            Dianne Cook              *
 *   dfs@research.att.com       dicook@iastate.edu          *
 *      (973) 360-8423    www.public.iastate.edu/~dicook/   *
 *                                                          *
 *                    Andreas Buja                          *
 *                andreas@research.att.com                  *
 *              www.research.att.com/~andreas/              *
 *                                                          *
 ************************************************************/

#include <limits.h>
#include <math.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"


double
randvalue(void) {
  double drand;

#ifdef USE_DRAND48
  drand = drand48();    /* rrand on [0.0,1.0] */
#else
  /* random() returns a value on [0, (2**31)-1], or [0, INT_MAX] */
  long lrand = (long) random();
  drand = (double) lrand / (double) INT_MAX;
#endif

  return drand;
}

void
rnorm2(double *drand, double *dsave) {

#ifdef USE_DRAND48
  *drand = 2.0 * drand48() - 1.0;
  *dsave = 2.0 * drand48() - 1.0;
#else
  long lrand, lsave;
  lrand = random();
  lsave = random();
  *drand = 2.0 * (double) lrand / (double) INT_MAX - 1.0;
  *dsave = 2.0 * (double) lsave / (double) INT_MAX - 1.0;
#endif

}

int
find_selected_cols(xg, cols)
  xgobidata *xg;
  int *cols;
{
  int i, ncols = 0;

  if (xg->is_plotting1d) {
    if (xg->plot1d_vars.y != -1)
      cols[ncols++] = xg->plot1d_vars.y;
    else if (xg->plot1d_vars.x != -1)
      cols[ncols++] = xg->plot1d_vars.x;
  }
  else if (xg->is_xyplotting) {
    cols[ncols++] = xg->xy_vars.x;
    cols[ncols++] = xg->xy_vars.y;
  }
  else if (xg->is_spinning) {
    cols[ncols++] = xg->spin_vars.x;
    cols[ncols++] = xg->spin_vars.y;
    cols[ncols++] = xg->spin_vars.z;
  }
  else if (xg->is_touring) {
    for (i=0; i<xg->numvars_t; i++)
      cols[ncols++] = xg->tour_vars[i];
  }
  else if (xg->is_corr_touring) {
    for (i=0; i<xg->ncorr_xvars; i++)
      cols[ncols++] = xg->corr_xvars[i];
    for (i=0; i<xg->ncorr_yvars; i++)
      cols[ncols++] = xg->corr_yvars[i];
  }
  return(ncols);
}

void
add_vgroups(xg, cols, ncols)
/*
 * If one of the chosen columns is in a vgroup,
 * add its comrades (unless they're already present)
*/
  xgobidata *xg;
  int *cols;
  int *ncols;
{
  int nc = *ncols;
  int j, k, n;

  for (j=0; j<nc; j++) {
    int vg = xg->vgroup_ids[cols[j]];

    for (k=0; k<xg->ncols_used; k++) {
      if (xg->vgroup_ids[k] == vg && k != cols[j]) {
        /* Got one; if it isn't already in cols, add it */
        Boolean addit = True;
        for (n=0; n<nc; n++) {
          if (cols[n] == k) {
            addit = False;
            break;
          }
        }
        if (addit) cols[(*ncols)++] = k;
        if (*ncols >= xg->ncols_used)
          break;
      }
    }
  }
}

int
fcompare(const void *x1, const void *x2)
{
  int val = 0;
  float *f1 = (float *) x1;
  float *f2 = (float *) x2;

  if (*f1 < *f2)
    val = -1;
  else if (*f1 > *f2)
    val = 1;

  return(val);
}

void
resort_vgroup_ids(xgobidata *xg, int *group_ids) {
  int maxid, i, id, newid, j;
  Boolean found;

  /*
   * Find maximum vgroup id.
  */
  maxid = 0;
  for (i=1; i<xg->ncols; i++) {
    if (group_ids[i] > maxid)
      maxid = group_ids[i];
  }

  /*
   * Find minimum vgroup id, set it to 0.  Find next, set it to 1; etc.
  */
  id = 0;
  newid = -1;
  while (id <= maxid) {
    found = false;
    for (j=0; j<xg->ncols; j++) {
      if (group_ids[j] == id) {
        newid++;
        found = true;
        break;
      }
    }
    if (found)
      for (j=0; j<xg->ncols; j++)
        if (group_ids[j] == id)
          group_ids[j] = newid;
    id++;
  }
}

/* Not used anywhere yet ... */
void
fshuffle(float *x, int n) {
/*
 * Knuth, Seminumerical Algorithms, Vol2; Algorithm P.
*/
  int i, k;
  float f;

  for (i=0; i<n; i++) {
    k = (int) (randvalue() * (double) i);
    f = x[i];
    x[i] = x[k];
    x[k] = f;
  }
}

/* ---------------------------------------------------------------------*/
/* The routines below have been added for the R/S connection */

int glyphIDfromName(char *glyphName) {
  int id = -1;

  if (strcasecmp(glyphName, "plus") == 0)
    id = PLUS_GLYPH;
  else if (strcasecmp(glyphName, "x") == 0)
    id = X_GLYPH;
  else if (strcasecmp(glyphName, "point") == 0)
    id = POINT_GLYPH;
  else if ((strcasecmp(glyphName, "open rectangle") == 0) ||
           (strcasecmp(glyphName, "open_rectangle") == 0) ||
           (strcasecmp(glyphName, "openrectangle") == 0))
    id = OPEN_RECTANGLE_GLYPH;
  else if ((strcasecmp(glyphName, "filled rectangle") == 0) ||
           (strcasecmp(glyphName, "filled_rectangle") == 0) ||
           (strcasecmp(glyphName, "filledrectangle") == 0))
    id = FILLED_RECTANGLE_GLYPH;
  else if ((strcasecmp(glyphName, "open circle") == 0) ||
           (strcasecmp(glyphName, "open_circle") == 0) ||
           (strcasecmp(glyphName, "opencircle") == 0))
    id = OPEN_CIRCLE_GLYPH;
  else if ((strcasecmp(glyphName, "filled circle") == 0) ||
           (strcasecmp(glyphName, "filled_circle") == 0) ||
           (strcasecmp(glyphName, "filledcircle") == 0))
    id = FILLED_CIRCLE_GLYPH;

  return id;
}

int glyphNames(char **names) {
  int i;
  static char* glyphNames[] =
    {"plus", "x", "openrectangle", "filledrectangle", "opencircle",
    "filledcircle", "point"};
  for (i=0; i<7; i++) names[i] = glyphNames[i];
  return (7);
}

int varno_from_name(xgobidata *xg, char *name) {
  int i, varno = -1;

  for(i = 0; i < xg->ncols_used; i++) {
    if(strcmp(name, xg->collab[i]) == 0) {
      varno = i;
      break;
    }
  }
  return varno;

}
