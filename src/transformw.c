/* transformw.c */
/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  the xgobi code and associated documentation,  provided  *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *        Copyright (c) 1990, ..., 1996 Bellcore            *
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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"
#include "DrawingA.h"

#define NDOMAINBTNS 4
static int domain_ind;
#define DOMAIN_OK       0
#define RAISE_MIN_TO_0  1
#define RAISE_MIN_TO_1  2
#define NEGATE          3
static Widget domain_menu_cmd, domain_menu_btn[NDOMAINBTNS] ;
static char *domain_menu_btn_label[] = {
  "No adjustment",
  "Raise minimum to 0",
  "Raise minimum to 1",
  "Negative",
};
static Widget *varlabel;

#define NTFORMS 12
#define RESTORE      0
#define APPLY_ADJ    1
#define POWER        2
#define ABSVALUE     3
#define INVERSE      4
#define LOG10        5
#define SCALE        6
#define STANDARDIZE  7
#define DISCRETE2    8
#define PERMUTE      9
#define SORT         10
#define NORMSCORE    11

Widget tform_cmd[NTFORMS];

static char *tform_names[] = {
  "Restore",
  "Apply adj",
  "Power (Box-Cox)",
  "Absolute value",
  "Inverse",
  "Base 10 log",
  "Scale to [0,1]",
  "Standardize",
  "Discretize: 2 levels",
  "Permute",
  "Sort",
  "Normal score"
};

char message[MSGLENGTH];
#define DOMAIN_ERROR sprintf(message, "Data outside the domain of function.\n")

static Widget *var_cbox;
static Widget tpane, tpopup = NULL;
static Position popupx = -1, popupy = -1;
static Widget exponent_lbl;
static float  exponent = 1.0;
static float  domain_incr = 0.0;
static int    ntform_cols, *tform_cols = NULL;

int transform(xgobidata *, int *, int, float,
  float (*)(float), int, double);

float no_change(float x)      { return x; }
float negate(float x)         { return -x; }

float raise_min_to_0(float x) { return (x + domain_incr); }
float raise_min_to_1(float x) { return (x + domain_incr + 1.0); }
float (*domain_adj)(float x);

float inv_raise_min_to_0(float x) { return (x - domain_incr); }
float inv_raise_min_to_1(float x) { return (x - domain_incr - 1.0); }
float (*inv_domain_adj)(float x);

/*ARGSUSED*/
static void
fallback(xgobidata *xg)
{
  XtCallCallbacks(tform_cmd[RESTORE], XtNcallback, (XtPointer) xg);
}

typedef struct {
  int tform;
  float domain_incr;
  float (*domain_adj)(float x);
  float (*inv_domain_adj)(float x);
  float param;
} TFormType;
TFormType *tform_tp = NULL;  /* the transformation applied to each variable */

void
alloc_transform_tp(xgobidata *xg)
{
  /* disallow transforming the groups variable */
  tform_tp = (TFormType *) XtRealloc((char *) tform_tp,
    (Cardinal) (xg->ncols-1) * sizeof(TFormType));
}

int
which_cols(int *cols, int varno, xgobidata *xg) {
/*
 * Figure out which columns to transform.
*/
  int j, ncols = 0;
  int groupno = xg->vgroup_ids[varno];

  /* allocated before this is called */
  for (j=0; j<(xg->ncols-1); j++) {
    if (xg->vgroup_ids[j] == groupno)
      cols[ncols++] = j;
  }
  return(ncols);
}

static void
set_initial_variable(xgobidata *xg) {
  int j, jvar, gid;

  if (xg->is_plotting1d)
    jvar = (xg->plot1d_vars.y != -1) ? xg->plot1d_vars.y : xg->plot1d_vars.x;
  else if (xg->is_xyplotting)
    jvar = xg->xy_vars.x;
  else if (xg->is_spinning)
    jvar = xg->spin_vars.x;
  else if (xg->is_touring)
    jvar = xg->tour_vars[0];
  else if (xg->is_corr_touring)
    jvar = xg->corr_xvars[0];

  gid = xg->vgroup_ids[jvar];
  for (j=0; j<xg->ncols-1; j++)
    if (xg->vgroup_ids[j] == gid) {
      XtVaSetValues(var_cbox[j], XtNstate, True, NULL);
      setToggleBitmap(var_cbox[j], True);
    }
}

static void
mean_stddev(xgobidata *xg, int *cols, int ncols, float (*stage1)(float),
  float *mean, float *stddev)
/*
 * Find the minimum and maximum values of a column or variable
 * group scaling by mean and std_width standard deviations.
 * Use the function pointer to domain_adj.
*/
{
  int i, j, n;
  float x;
  double sumxi = 0.0, sumxisq = 0.0;
  double dx, dmean, dvar, dstddev;
  double dn = (double) (ncols * xg->nrows);

  for (n=0; n<ncols; n++) {
    j = cols[n];
    for (i=0; i<xg->nrows; i++) {
      dx = (double) (*stage1)(xg->raw_data[i][j]);
      sumxi = sumxi + dx;
      sumxisq = sumxisq + dx * dx;
    }
  }
  dmean = sumxi / dn;
  dvar = (sumxisq / dn) - (dmean * dmean);
  dstddev = sqrt(dvar);

  *mean = (float) dmean;
  *stddev = (float) dstddev;
}

float
median(xgobidata *xg, float **data, int *cols, int ncols)
{
/*
 * Find the minimum and maximum values of each column or variable
 * group scaling by median and largest distance
*/
  int i, j, n, np;
  float *x;
  double dmedian = 0;
  extern int fcompare(const void *, const void *);

  np = ncols * xg->nrows;
  x = (float *) XtMalloc((Cardinal) np * sizeof(float));
  for (n=0; n<ncols; n++) {
    j = cols[n];
    for (i=0; i<xg->nrows; i++) {
      x[n*xg->nrows_in_plot + i] = data[i][j];
    }
  }

  qsort((void *) x, np, sizeof(float), fcompare);
  dmedian = ((np % 2) != 0) ?  x[(np-1)/2] : (x[np/2-1] + x[np/2])/2. ;

  return (float) dmedian;
}

void
reset_tform(xgobidata *xg) {
  int j;

  for (j=0; j<xg->ncols-1; j++) {
    tform_tp[j].tform = RESTORE;
    tform_tp[j].domain_incr = 0.;
    tform_tp[j].param = 0.;
    tform_tp[j].domain_adj = no_change;
    tform_tp[j].inv_domain_adj = no_change;
  }
}

static int 
sort_compare (float *val1, float *val2)
{
  if (*val1 < *val2) 
    return (-1);
  else if (*val1 == *val2)
    return (0);
  else 
    return (1);
}

int
transform(xgobidata *xg, int *cols, int ncols, float domain_incr,
  float (*domain_adj)(float), int tfnum, double param)
{
  int i, j, n;
  float min, max, diff;
  float mean, stddev;
  float fmedian, ref;
  Boolean allequal;
  int tmpi, numperm, indx_flag;
  double dtmp;

/* This depends on the selected variables; for now, it's
   ok for a single variable
*/
  switch (domain_ind) {
    case DOMAIN_OK:
      domain_incr = 0;
      domain_adj = no_change;
      inv_domain_adj = no_change;
      break;
    case RAISE_MIN_TO_0:
      domain_incr = fabs(xg->lim_raw[ tform_cols[0] ].min);
      domain_adj = raise_min_to_0;
      inv_domain_adj = inv_raise_min_to_0;
      break;
    case RAISE_MIN_TO_1:
      domain_incr = fabs(xg->lim_raw[ tform_cols[0] ].min) + 1.0;
      domain_adj = raise_min_to_1;
      inv_domain_adj = inv_raise_min_to_1;
      break;
    case NEGATE:
      domain_incr = 0.0;
      domain_adj = negate;
      inv_domain_adj = negate;
      break;
    default:
      domain_incr = 0;
      domain_adj = no_change;
      inv_domain_adj = no_change;
  }


  switch(tfnum)
  {
    case RESTORE:    /* Restore original values -- set domain adj to null */

      XtCallCallbacks(domain_menu_btn[DOMAIN_OK], XtNcallback, (XtPointer) xg);
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = xg->raw_data[i][j];

        (void) strcpy(xg->collab_tform[j], xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case APPLY_ADJ:    /* Apply domain adj */

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = (*domain_adj)(xg->raw_data[i][j]);

        (void) strcpy(xg->collab_tform[j], xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case POWER:  /* Box-Cox power transform family */

      if (fabs(param-0) < .001) {       /* Natural log */
        for (n=0; n<ncols; n++) {
          j = cols[n];
          for (i=0; i<xg->nrows; i++) {
            if ((*domain_adj)(xg->raw_data[i][j]) <= 0) {
              printf(stderr, "%f %f\n",
                xg->raw_data[i][j], (*domain_adj)(xg->raw_data[i][j]));
              DOMAIN_ERROR;
              show_message(message, xg);
              return(0);
            }
          }
        }
        for (n=0; n<ncols; n++) {
          j = cols[n];
          for (i=0; i<xg->nrows; i++) {
            xg->tform_data[i][j] = (float)
              log((double) ((*domain_adj)(xg->raw_data[i][j])));
          }

          (void) sprintf(xg->collab_tform[j], "ln(%s)", xg->collab[j]);
          XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
        }
      }

      else {

        for (n=0; n<ncols; n++) {
          j = cols[n];
          for (i=0; i<xg->nrows; i++) {
            dtmp = pow((double) (*domain_adj)(xg->raw_data[i][j]), param);
            dtmp = (dtmp - 1.0) / param;

            /* If dtmp no good, restore and return */
            if (!finite(dtmp)) {
              printf(stderr, "%f %f %f\n",
                xg->raw_data[i][j], (*domain_adj)(xg->raw_data[i][j]), dtmp);
              DOMAIN_ERROR;
              show_message(message, xg);
              fallback(xg);
              return (0);
            }
            xg->tform_data[i][j] = (float) dtmp;
          }

          (void) sprintf(xg->collab_tform[j], "B-C(%s,%.2f)",
            xg->collab[j], param);
          XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
        }
      }
      break;

    case ABSVALUE:
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          if ((xg->raw_data[i][j] + domain_incr) < 0)
            xg->tform_data[i][j] = (float)
              fabs((double)(*domain_adj)(xg->raw_data[i][j])) ;
          else
            xg->tform_data[i][j] = (*domain_adj)(xg->raw_data[i][j]);
        }

        (void) sprintf(xg->collab_tform[j], "Abs(%s)", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case INVERSE:    /* 1/x */
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          if ( (*domain_adj)(xg->raw_data[i][j]) == 0) {
            DOMAIN_ERROR;
            show_message(message, xg);
            return(0);
          }
        }
      }

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          xg->tform_data[i][j] = (float)
            pow((double) (*domain_adj)(xg->raw_data[i][j]),
              (double) (-1.0));
        }

        (void) sprintf(xg->collab_tform[j], "1/%s", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case LOG10:    /* Base 10 log */
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          if ( (*domain_adj)(xg->raw_data[i][j]) <= 0) {
            DOMAIN_ERROR;
            show_message(message, xg);
            return(0);
          }
        }
      }
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          xg->tform_data[i][j] = (float)
            log10((double) (*domain_adj)(xg->raw_data[i][j]));
        }

        (void) sprintf(xg->collab_tform[j], "log10(%s)", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case SCALE:    /* Map onto [0,1] */
      /* First find min and max; they get updated after transformations */

      min = max = (*domain_adj)(xg->raw_data[0][cols[0]]);
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          if ( (ref = (*domain_adj)(xg->raw_data[i][j])) < min)
            min = ref;
          else if (ref > max) max = ref;
        }
      }

      adjust_limits(&min, &max);
      diff = max - min;

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = 
            ((*domain_adj)(xg->raw_data[i][j]) - min)/diff;

        (void) sprintf(xg->collab_tform[j], "%s [0,1]", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case STANDARDIZE:    /* (x-mean)/sigma */

      mean_stddev(xg, cols, ncols, domain_adj, &mean, &stddev);
      /* DOMAIN_ERROR if stddev == 0 */

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] =
            ((*domain_adj)(xg->raw_data[i][j]) - mean)/stddev;

        (void) sprintf(xg->collab_tform[j], "(%s-m)/s", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case DISCRETE2:    /* x>median */
      /* refuse to discretize if all values are the same */
      allequal = True;
      ref = xg->raw_data[0][cols[0]];
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          if (xg->raw_data[i][j] != ref) {
            allequal = False;
            break;
          }
        }
      }
      if (allequal) {
        DOMAIN_ERROR;
        show_message(message, xg);
        return(0);
      }

      /* First find median */

      fmedian = median(xg, xg->raw_data, cols, ncols);
      fmedian = (*domain_adj)(fmedian);

      /* Then find the true min and max */
      min = max = (*domain_adj)(xg->raw_data[0][cols[0]]);
      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          if ( (ref = (*domain_adj)(xg->raw_data[i][j])) < min)
            min = ref;
          else if (ref > max) max = ref;
        }
      }

      /* This prevents the collapse of the data in a special case */
      if (max == fmedian)
        fmedian = (min + max)/2.0;

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] =
            ( (*domain_adj)(xg->raw_data[i][j]) > fmedian ) ? 1.0 : 0.0;

        (void) sprintf(xg->collab_tform[j], "%s:0,1", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      break;

    case PERMUTE:    /* (x-mean)/sigma */
    { /* This curly brace allows me to define variables for this case */

      /* First generate a random vector */
      /* Allocate array for permutation indices */

      int *permindx;
      permindx = (int *) XtMalloc((Cardinal) xg->nrows * sizeof(int));
      numperm = 0;
      while (numperm < xg->nrows) {

#ifdef USE_DRAND48
        tmpi = (int) (drand48() * (double)xg->nrows);
#else
        tmpi =  (int) ((double) random()/ (double) INT_MAX * xg->nrows);
#endif

        indx_flag = 0;
        for (i=0; i<numperm; i++) {
          if (tmpi == permindx[i]) {
            indx_flag = 1;
	      }
	    }
        if (!indx_flag) {
          permindx[numperm] = tmpi;
          numperm++;
        }
      }

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = (*domain_adj)(xg->raw_data[permindx[i]][j]);

        (void) sprintf(xg->collab_tform[j], "perm(%s)", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }

      XtFree((XtPointer) permindx);
    }
      break;

    case SORT: 
    {
      /*  create sort index  here - mallika */
      float *sort_data; /* mallika */
      /* Allocate array for sorted columns - mallika */
      sort_data = (float *) XtMalloc((Cardinal) xg->nrows * sizeof(float));

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          sort_data[i] = (*domain_adj)(xg->raw_data[i][j]);

        qsort((char *) sort_data, xg->nrows, sizeof(float), sort_compare);
   
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = sort_data[i]; 

        (void) sprintf(xg->collab_tform[j], "sort(%s)", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      XtFree((XtPointer) sort_data);
    }
      break;

    case NORMSCORE:
    {
      float *norm_score_data;
      float ftmp;

      /* Allocate array for normalized scores */
      norm_score_data = (float *)
        XtMalloc((Cardinal) xg->nrows * sizeof(float));

     for (n=0; n<ncols; n++) {
        float normmean=0, normvar=0;
        j = cols[n];
        for (i=0; i<xg->nrows; i++) {
          ftmp = (*domain_adj)(xg->raw_data[i][j]);
          norm_score_data[i] = ftmp;
          normmean += ftmp;
          normvar += (ftmp * ftmp);
        }
        normmean /= xg->nrows;
        normvar = (float)sqrt((float)(normvar/xg->nrows - normmean*normmean));
        for (i=0; i<xg->nrows; i++)
          norm_score_data[i] = (norm_score_data[i]-normmean)/normvar;

        for (i=0; i<xg->nrows; i++) {
          if (norm_score_data[i]>0)
            norm_score_data[i] = erf(norm_score_data[i]/sqrt(2.))/
              2.8284271+0.5;
          else if (norm_score_data[i]<0)
            norm_score_data[i] = 0.5 - erf((float) fabs((double) 
              norm_score_data[i])/sqrt(2.))/2.8284271;
          else 
            norm_score_data[i]=0.5;
        }
        
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = norm_score_data[i]; 

        (void) sprintf(xg->collab_tform[j], "normsc(%s)", xg->collab[j]);
        XtVaSetValues(varlabel[j], XtNlabel, xg->collab_tform[j], NULL);
      }
      XtFree((XtPointer) norm_score_data);/* mallika */
    }
    break;
  }
  return(1);
}

static void
tform_response(xgobidata *xg, int *cols, int ncols,
  float domain_incr, float (*stage1)(float), float (*inv_stage1)(float),
  int tfno, double param)
{
  int j, n, gid;
  Boolean reset_vgroups = false;

  /*
   * Reset vgroups?
  */
  gid = xg->vgroup_ids[cols[0]];
  for (n=1; n<ncols; n++) {
    if (xg->vgroup_ids[cols[n]] != gid) {
      reset_vgroups = true;
      break;
    }
  }
  if (reset_vgroups) {
    gid = numvargroups(xg)+1;
    for (n=0; n<ncols; n++) {
      xg->vgroup_ids[cols[n]] = gid;
    }

    resort_vgroup_ids(xg, xg->vgroup_ids);
  }

  if (xg->ncols_used > 2)
    update_sphered(xg, cols, ncols);
  update_lims(xg);
  update_world(xg);

  /* Set tform_tp[] for transformed columns */
  for (n=0; n<ncols; n++) {
    tform_tp[cols[n]].tform = tfno;
    tform_tp[cols[n]].domain_incr = domain_incr;
    tform_tp[cols[n]].param = param;
    tform_tp[cols[n]].domain_adj = stage1;
    tform_tp[cols[n]].inv_domain_adj = inv_stage1;
  }

  world_to_plane(xg);
  plane_to_screen(xg);

  /*
    This bit of init_axes() is needed.
  */
  for (n=0; n<ncols; n++) {
    j = cols[n];
    xg->nicelim[j].min = xg->lim0[j].min;
    xg->nicelim[j].max = xg->lim0[j].max;
    SetNiceRange(j, xg);
    xg->deci[j] = set_deci(xg->tickdelta[j]);
  }

  if (xg->is_xyplotting) {
    init_ticks(&xg->xy_vars, xg);
  }
  else if (xg->is_plotting1d)
    init_ticks(&xg->plot1d_vars, xg);

  if (xg->is_brushing) {
    assign_points_to_bins(xg);
    if (xg->brush_mode == transient)
      reinit_transient_brushing(xg);
  }

  plot_once(xg);

  if (xg->is_cprof_plotting)
    update_cprof_plot(xg);
}

/* ARGSUSED */
static XtCallbackProc
tform_cback(Widget w, xgobidata *xg, XtPointer cbdata)
{
  int tfno, ntform_cols, j, k, groupno;
  Boolean state;

  for (tfno=0; tfno<NTFORMS; tfno++)
    if (w == tform_cmd[tfno])
      break;

  ntform_cols = 0;
  for (j=0; j<xg->ncols-1; j++) {
    XtVaGetValues(var_cbox[j], XtNstate, &state, NULL);
    if (state) {
      tform_cols[ntform_cols++] = j;
      groupno = xg->vgroup_ids[j];
      for (k=j+1; k<xg->ncols-1; k++) {
        if (xg->vgroup_ids[k] == groupno)
          tform_cols[ntform_cols++] = k;
      }
    }
  }

  if (ntform_cols > 0) {
    if (transform(xg, tform_cols, ntform_cols, domain_incr,
      domain_adj, tfno, (double) exponent))
    {
      tform_response(xg, tform_cols, ntform_cols, domain_incr,
        domain_adj, inv_domain_adj, tfno, (double) exponent);
    }
  }
}

/* power transform */
static void
reset_exp(xgobidata *xg) {
  char lbl[32];
  Boolean state;
  int j, k, groupno;
  int tfno = POWER;

  /* reset the power label */
  sprintf(lbl, "%4.2f", exponent);
  XtVaSetValues(exponent_lbl, XtNlabel, lbl, NULL);

  ntform_cols = 0;
  for (j=0; j<xg->ncols-1; j++) {
    XtVaGetValues(var_cbox[j], XtNstate, &state, NULL);
    if (state) {
      tform_cols[ntform_cols++] = j;
      groupno = xg->vgroup_ids[j];
      for (k=j+1; k<xg->ncols-1; k++) {
        if (xg->vgroup_ids[k] == groupno)
          tform_cols[ntform_cols++] = k;
      }
    }
  }

  if (ntform_cols > 0) {
    if (transform(xg, tform_cols, ntform_cols, domain_incr, domain_adj,
      tfno, (double) exponent))
    {
      tform_response(xg, tform_cols, ntform_cols,
        domain_incr, domain_adj, inv_domain_adj,
        tfno, (double) exponent);
    }
  }
}

/* ARGSUSED */
static XtCallbackProc
reduce_exp_cback(Widget w, xgobidata *xg, XtPointer slideposp)
{
  exponent -= .01;
  reset_exp(xg);
}
/* ARGSUSED */
static XtCallbackProc
increase_exp_cback(Widget w, xgobidata *xg, XtPointer slideposp)
{
  exponent += .01;
  reset_exp(xg);
}

/* ARGSUSED */
static XtCallbackProc
exp_sbar_cback(Widget w, xgobidata *xg, XtPointer slideposp)
{
  int iexp, third_dec_place;
  float ftmp;
  float fslidepos = * (float *) slideposp;

  /* rescale from [0,1] to [-4, 4] */
  ftmp = 8*fslidepos - 4;

/*
 * restrict the exponent to two decimal places
*/

  third_dec_place = ((int) (ftmp * 1000)) % 10;
  iexp = (int) (ftmp * 100.);
  if (third_dec_place < 5)
    exponent = ((float) iexp) / 100.;
  else
    exponent = ((float) (iexp+1)) / 100.;

  reset_exp(xg);
}

/* ARGSUSED */
static XtCallbackProc
var_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  Boolean state;
  XtVaGetValues(w, XtNstate, &state, NULL);
  setToggleBitmap(w, state);
}

/* ARGSUSED */
static XtCallbackProc
set_domain_incr_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  int k;

  domain_ind = DOMAIN_OK;
  for (k=0; k<NDOMAINBTNS; k++)
    if (domain_menu_btn[k] == w) {
      domain_ind = k;
      break;
    }

  XtVaSetValues(domain_menu_cmd,
    XtNlabel, domain_menu_btn_label[domain_ind],
    NULL);
}

/* ARGSUSED */
static XtCallbackProc
reset_vgroups_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  int j;
  for (j=0; j<xg->ncols-1; j++)
    xg->vgroup_ids[j] = xg->vgroup_ids_ori[j];

  /*
   * I don't know what has to happen at this point ...
  */

}

/* ARGSUSED */
static XtCallbackProc
close_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  XtDestroyWidget(tpopup);
  tpopup = NULL;
}


#include <X11/Xaw/Viewport.h>
/* ARGSUSED */
XtCallbackProc
open_tform_popup_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  Widget close;
  Widget form0;
  Dimension width, height;
  register int j, k;
  Widget box_variables, box_varlabels, box_tforms;
  Widget power_form, power_lbl, reduce_exp_arr, exp_sbar;
  Widget increase_exp_arr;
  char str[32];
  static char *domain_menu_str = "Domain adjustment:";
  Widget domain_menu_box, domain_menu_lab, domain_menu;
  Widget vport, vport_child, reset_vgroups_cmd;
  Widget box_stage[4];
  Boolean doit = false;

  Dimension maxwidth = 0;

  if (tpopup == NULL) {
    tform_cols = (int *) XtMalloc((Cardinal) (xg->ncols-1) * sizeof(int));
    var_cbox = (Widget *) XtMalloc((Cardinal) (xg->ncols-1) * sizeof(Widget));

    alloc_transform_tp(xg);

    if (popupx == -1 && popupy == -1) {
      XtVaGetValues(xg->workspace,
        XtNwidth, &width,
        XtNheight, &height, NULL);
      XtTranslateCoords(xg->workspace,
        (Position) width, (Position) (height/2), &popupx, &popupy);
    }

    tpopup = XtVaCreatePopupShell("Variable Transformation",
      topLevelShellWidgetClass, xg->shell,
      XtNx,        popupx,
      XtNy,        popupy,
      XtNinput,    True,
      XtNtitle,    "Transform variables",
      XtNiconName, "Tform",
      NULL);
    if (mono) set_mono(tpopup);

    /*
     * Create a paned widget so the 'Click here ...'
     * can be all across the bottom.
    */
    tpane = XtVaCreateManagedWidget("Form",
      panedWidgetClass, tpopup,
      XtNorientation, (XtOrientation) XtorientVertical,
      XtNresizable, False,
      NULL);

    form0 = XtVaCreateManagedWidget("Form",
      formWidgetClass, tpane,
      XtNresizable, False,
      NULL);
    if (mono) set_mono(form0);

    box_tforms = XtVaCreateManagedWidget("Close",
      boxWidgetClass, form0,
      XtNorientation, (XtOrientation) XtorientVertical,
      XtNleft, (XtEdgeType) XtChainLeft,
      XtNright, (XtEdgeType) XtChainLeft,
      XtNtop, (XtEdgeType) XtChainTop,
      XtNbottom, (XtEdgeType) XtChainTop,
      NULL);

    /* Stage 0: domain adjustment */
    box_stage[0] = XtVaCreateManagedWidget("Close",
      boxWidgetClass, box_tforms,
      XtNorientation, (XtOrientation) XtorientHorizontal,
      NULL);

    domain_ind = DOMAIN_OK;
    build_labelled_menu(&domain_menu_box, &domain_menu_lab, domain_menu_str,
      &domain_menu_cmd, &domain_menu, domain_menu_btn,
      domain_menu_btn_label, domain_menu_btn_label,  /* no nicknames */
      NDOMAINBTNS, domain_ind, box_stage[0], NULL,
      XtorientHorizontal, appdata.font, "Transformations", xg);
    for (k=0; k<NDOMAINBTNS; k++)
      XtAddCallback(domain_menu_btn[k],  XtNcallback,
        (XtCallbackProc) set_domain_incr_cback, (XtPointer) xg);

    /* Stage 1: most transformations */
    box_stage[1] = XtVaCreateManagedWidget("Close",
      boxWidgetClass, box_tforms,
      NULL);

    XtVaCreateManagedWidget("Transformation Stage 1:",
      labelWidgetClass, box_stage[1],
      NULL);

    /* Power family */
    power_form = XtVaCreateManagedWidget("Form",
      boxWidgetClass, box_stage[1],
      XtNresizable, False,
      XtNorientation, XtorientHorizontal,
      XtNhSpace, 1,
      XtNvSpace, 1,
      NULL);
    power_lbl = XtVaCreateManagedWidget("Label",
      labelWidgetClass, power_form,
      XtNlabel, (String) "Power (Box-Cox)",
      NULL);

    reduce_exp_arr = XtVaCreateManagedWidget("Icon",
      commandWidgetClass, power_form,
      XtNinternalHeight, (Dimension) 0,
      XtNinternalWidth, (Dimension) 0,
      XtNborderColor, (Pixel) appdata.fg,
      XtNbitmap, (Pixmap) leftarr,
      NULL);
    if (mono) set_mono(reduce_exp_arr);
    XtAddCallback(reduce_exp_arr, XtNcallback,
     (XtCallbackProc) reduce_exp_cback, (XtPointer) xg);

    sprintf(str, "Power Tformations");
    width = XTextWidth(appdata.font, str, strlen(str));

    exp_sbar = XtVaCreateManagedWidget("Scrollbar",
      scrollbarWidgetClass, power_form,
      XtNhorizDistance, (Dimension) 0,
      XtNwidth, (Dimension) width,
      XtNheight, (Dimension) 20,
      XtNorientation, (XtOrientation) XtorientHorizontal,
      NULL);
    if (mono) set_mono(exp_sbar);
    XtAddCallback(exp_sbar, XtNjumpProc,
     (XtCallbackProc) exp_sbar_cback, (XtPointer) xg);

    /*
     * -4,4 -> 0,1
    */

    XawScrollbarSetThumb(exp_sbar, (1. + 4)/8., -1.);

    increase_exp_arr = XtVaCreateManagedWidget("Icon",
      commandWidgetClass, power_form,
      XtNinternalHeight, (Dimension) 0,
      XtNinternalWidth, (Dimension) 0,
      XtNborderColor, (Pixel) appdata.fg,
      XtNhorizDistance, (Dimension) 0,
      XtNbitmap, (Pixmap) rightarr,
      NULL);
    if (mono) set_mono(increase_exp_arr);
    XtAddCallback(increase_exp_arr, XtNcallback,
     (XtCallbackProc) increase_exp_cback, (XtPointer) xg);

    (void) sprintf(str, "%1.2f", -9.99);
    width = XTextWidth(appdata.font, str,
        strlen(str) + 2*ASCII_TEXT_BORDER_WIDTH);
    (void) sprintf(str, "%4.2f", exponent);
    exponent_lbl = XtVaCreateManagedWidget("StdizeLabel",
      labelWidgetClass,  power_form,
      XtNlabel, (String) str,
      XtNwidth, width,
      NULL);
    if (mono) set_mono(exponent_lbl);

    /* From Absolute Value to Normal Score */
    for (j=0; j<NTFORMS; j++) {
      if (j != 2) {
        tform_cmd[j] = CreateCommand(xg, tform_names[j], True, 
          NULL, NULL, box_stage[1], "Transformations");
        XtManageChild(tform_cmd[j]);
        XtAddCallback(tform_cmd[j], XtNcallback,
          (XtCallbackProc) tform_cback, (XtPointer) xg);
      }
    }

    /* Stage 2: permutation, sorting, sphering */
/*
    box_stage[2] = XtVaCreateManagedWidget("Close",
      boxWidgetClass, box_tforms,
      NULL);
*/


    vport = XtVaCreateManagedWidget("ViewPort",
      viewportWidgetClass, form0,
      XtNallowHoriz, False,
      XtNallowVert,  True,
      XtNfromHoriz, box_tforms,
      XtNleft, (XtEdgeType) XtChainLeft,
      XtNmappedWhenManaged, False,
      NULL);
    vport_child = XtVaCreateManagedWidget("Box",
      boxWidgetClass, vport,
      XtNorientation, (XtOrientation) XtorientHorizontal,
      NULL);

    box_variables = XtVaCreateManagedWidget("Box",
      formWidgetClass, vport_child,
      XtNleft, (XtEdgeType) XtChainLeft,
      XtNright, (XtEdgeType) XtChainLeft,
      XtNtop, (XtEdgeType) XtChainTop,
      XtNbottom, (XtEdgeType) XtChainTop,
      NULL);
    for (j=0; j<xg->ncols-1; j++) {
      var_cbox[j] = CreateToggle(xg, xg->collab[j], True, 
        NULL, NULL, NULL, False, ANY_OF_MANY, box_variables, "Transformations");
      if (j>0)
        XtVaSetValues(var_cbox[j], XtNfromVert, var_cbox[j-1], NULL);
      XtAddCallback(var_cbox[j], XtNcallback,
       (XtCallbackProc) var_cback, (XtPointer) xg);
    }
    XtManageChildren(var_cbox, xg->ncols-1);

    reset_vgroups_cmd = CreateCommand(xg, "Reset vgroups", True, 
      box_tforms, vport, form0, "Transformations");
    XtManageChild(reset_vgroups_cmd);
    XtAddCallback(reset_vgroups_cmd, XtNcallback,
     (XtCallbackProc) reset_vgroups_cback, (XtPointer) xg);

    box_varlabels = XtVaCreateManagedWidget("Form",
      formWidgetClass, vport_child,
      XtNorientation, (XtOrientation) XtorientVertical,
      XtNfromHoriz, box_variables,
      NULL);

    for (j=0; j<xg->ncols-1; j++) {
      width = XTextWidth(appdata.font, xg->collab_tform[j],
        strlen(xg->collab_tform[j])) + 2*ASCII_TEXT_BORDER_WIDTH;
      maxwidth = (width > maxwidth) ? width : maxwidth;
    }
    width = XTextWidth(appdata.font, "normsc()", strlen("normsc()")) +
      2*ASCII_TEXT_BORDER_WIDTH;
    maxwidth += width;

    varlabel = (Widget *) XtMalloc((Cardinal) (xg->ncols-1) * sizeof(Widget));
    for (j=0; j<xg->ncols-1; j++) {
      varlabel[j] = XtVaCreateWidget("Label",
        labelWidgetClass, box_varlabels,
        XtNlabel, (String) xg->collab_tform[j],
        XtNfromVert, (j>0) ? varlabel[j-1] : NULL,
        XtNwidth, maxwidth,
        NULL);
    }
    XtManageChildren(varlabel, xg->ncols-1);

    close = XtVaCreateManagedWidget("Close",
      commandWidgetClass, tpane,
      XtNshowGrip, (Boolean) False,
      XtNskipAdjust, (Boolean) True,
      XtNlabel, (String) "Click here to dismiss",
      NULL);
    if (mono) set_mono(close);
    XtAddCallback(close, XtNcallback,
      (XtCallbackProc) close_cback, (XtPointer) xg);

    set_initial_variable(xg);

    doit = true;
  }

  XtPopup(tpopup, (XtGrabKind) XtGrabNone);
  set_wm_protocols(tpopup);
  XRaiseWindow(display, XtWindow(tpopup));

  if (doit) {
    Dimension hgt;
    XtVaGetValues(box_tforms, XtNheight, &hgt, NULL);
    XtVaSetValues(vport, XtNheight, hgt, NULL);
    XtMapWidget(vport);
  }
}

/**** Called from elsewhere ****/

/********** do the reverse transform; needed by move_points *******/

#define signum(x) (((x)<0.0)?(-1.0):(((x)>0.0)?(1.0):(0.0)))
#define INV_DOMAIN_ERROR sprintf(message, \
  "Data outside the domain of function; can\'t do the transformation.\n")

float
inv_transform(int icase, int jvar, xgobidata *xg) {
  double tx = xg->tform_data[icase][jvar];
  double rx = xg->raw_data[icase][jvar];
  double new_rx = tx;

  float min, max, diff;
  int *cols;
  int ncols = 0, i, j, n;
  float mean, stddev;

  cols = (int *) XtMalloc((Cardinal) (xg->ncols-1) * sizeof(int));
  ncols = which_cols(cols, jvar, xg);

  switch (tform_tp[jvar].tform) {
    case RESTORE:
      /* new_rx = tx; */
      break;

    case POWER:

      if (fabs(tform_tp[jvar].param) < .001)
        new_rx = exp(tx);
      else
        new_rx = pow(
          (double) ((tx * tform_tp[jvar].param) + 1),
          (double) (1.0/tform_tp[jvar].param));

      if (!finite(new_rx)) {
        INV_DOMAIN_ERROR;
        show_message(message, xg);
        new_rx = tx;
      }
      break;

    case ABSVALUE:
      new_rx = tx * signum(rx);
      break;

    case INVERSE:
      if (tx == 0) {
        INV_DOMAIN_ERROR;
        show_message(message, xg);
      } else
        new_rx = 1.0/tx;
      break;

    case LOG10:    /* Base 10 log */
      new_rx = pow(10.0, tx);
      break;

    case SCALE:    /* Map onto [0,1] */

      min_max(xg, xg->raw_data, cols, ncols, &min, &max);
      adjust_limits(&min, &max);
      diff = max - min;

      new_rx = (tx * diff) + min;
      break;

    case STANDARDIZE:    /* (x-mean)/sigma */
      mean_stddev(xg, cols, ncols, tform_tp[jvar].inv_domain_adj,
        &mean, &stddev);
      new_rx = (tx * stddev) + mean;
      break;

    case DISCRETE2:    /* x>median */
      show_message(
        "Sorry, I\'m unable to perform the transformation for this point\n", 
        xg);
      break;

/*
 * for the rest of these, I've apparently just decided to restore
 * the data to its untransformed state.  Let's worry about it later ...
*/

    case PERMUTE:    
      break;

    case SORT:    
      break;

    case NORMSCORE:
      cols = (int *) XtMalloc((Cardinal) (xg->ncols-1) * sizeof(int));
      ncols = which_cols(cols, jvar, xg);

      for (n=0; n<ncols; n++) {
        j = cols[n];
        for (i=0; i<xg->nrows; i++)
          xg->tform_data[i][j] = xg->raw_data[i][j];
      }
      break;

    default:
      /* new_rx = tx; */
      break;
  }

  inv_domain_adj = tform_tp[jvar].inv_domain_adj;
  return((float) ((*inv_domain_adj)(new_rx)));
}

void
permute_again(xgobidata *xg) {
  if (xg->is_touring && (xg->is_princ_comp || xg->is_pp))
    ;
  else {
    if (ntform_cols > 0 && tform_cols != NULL) {
      if ( transform(xg, tform_cols, ntform_cols, domain_incr, domain_adj,
        PERMUTE, 0.0) )
      {
        tform_response(xg, tform_cols, ntform_cols, domain_incr,
          domain_adj, inv_domain_adj, PERMUTE, 0.0);
      }
    }
  }
}

void
restore_variables(xgobidata *xg) {
  if (ntform_cols > 0 && tform_cols != NULL) {
    if ( transform(xg, tform_cols, ntform_cols, domain_incr, domain_adj,
      RESTORE, 0.0) )
    {
      tform_response(xg, tform_cols, ntform_cols, domain_incr,
        domain_adj, inv_domain_adj, RESTORE, 0.0);
    }
  }
}

void
transform_all(xgobidata *xg) {
  int j, cols[1];

  for (j=0; j<xg->ncols-1; j++) {
    cols[0] = j;
    if (tform_tp != (TFormType *) NULL && &tform_tp[j] != (TFormType *) NULL)
      transform(xg, cols, 1,
        tform_tp[j].domain_incr,
        tform_tp[j].domain_adj,
        tform_tp[j].tform, tform_tp[j].param);
  }
}
