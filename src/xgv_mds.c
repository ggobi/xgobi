/* AB: changed loss function, fixed old mds_respow =2
   May 18, 1998: reused mds_respow as power of ||x_i-x_j||^mds_respow
   May 25, 1998: reused mds_respow for w_ij = D_ij^(mds_respow-2)
*/
/* dfs:
   Oct 15, 1998: using glyph instead of label for anchoring
*/

/* mds.c: multidimensional scaling (a la andreas) ML 2/92. */
/* molecule.c: run as part of xgvis ML 2/92. */
/*
 * molecule.c: find layout with attraction/repulsion.
 * dists.c: writes out distance matrix.
 * ML: 11/91.
*/

/* Includes. */
#include <sys/types.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include <X11/keysym.h>
#include "xgvis.h"

/* Functions. */
double distance(double *, double *, int, double);    /* In molecule.c. */
double normalize(double *, int, double);

extern void clear_array(struct array *);
extern void copy_array(struct array *, struct array *);
extern void add_stress_value(double);
extern void make_empty_array(struct array *, int, int);
extern void draw_stress(void);
extern void zero_array(struct array *);
extern void scale_array_mean(struct array *, int, int);
extern void update_shepard_labels(int);

/* Macros. */
#define DOUBLE(x) ((double)(x))
#define FLOAT(x) ((float)(x))
#define SAMEGLYPH(i,j) \
( xgobi.color_now[(i)] == xgobi.color_now[(j)] && \
  xgobi.glyph_now[(i)].type == xgobi.glyph_now[(j)].type && \
  xgobi.glyph_now[(i)].size == xgobi.glyph_now[(j)].size ) \

#define CURRENTGLYPH(i) \
( xgobi.color_now[(i)] == xgobi.color_id && \
  xgobi.glyph_now[(i)].type == xgobi.glyph_id.type && \
  xgobi.glyph_now[(i)].size == xgobi.glyph_id.size ) \

/*
static Boolean
sticky_case(int jcase)
{
  int j;
  Boolean sticky = false;
  for (j=0; j<xgobi.nsticky_ids; j++) {
    if (xgobi.sticky_ids[j] == jcase) {
      sticky = True;
      break;
    }
  }
  return (sticky);
}
*/

/* MDS vars:
   Boolean mds_group  INIT(= FALSE);
   double mds_lnorm  INIT(= 2.0);
   double mds_power  INIT(= 1.0);
   double mds_weightpow INIT(= 0.0);
   int mds_n_iters   INIT(= 1);
   double mds_stepsize  INIT(= 0.005?);
   int mds_dims INIT(= 3);
*/

/*
                     p                q
  Stress = sum  |  D   - ||x - x ||  |
           ij       ij      i   j  l

D = d
p = mds_power
l = mds_lnorm
q = mds_weightpow

*/

/*
 * n_iters iterations of mds.  Implement this as if it
 * were a work proc.
*/

double delta = 1.0/100000.0;
#define signum(x) (((x)<0.0)?(-1.0):(((x)>0.0)?(1.0):(0.0)))

/*
 * Perform one loop of the iterative mds function.
 *
 * If doit is False, then we really want to determine the
 * stress function without doing anything to the gradient
*/
void
mds_once(Boolean doit, Boolean shepard, FILE* fpdat, FILE* fprow)
{
  int i, j, k;
  static struct array pos_grad;
  static Boolean cleared_p = False;
  static int nactive_distances = -1, prev_nactive_distances = -1;

  double d;
  double step_mag, step_dir;
  double gsum, psum, gfactor;
  int nchanged;
  double dist_goal;
  double resid;
  double weight;
  double stress, stress_dx = 1.0, stress_dd = 1.0, stress_xx = 1.0;
  double accum_dx, accum_dd, accum_xx;
  double f1, f2;
  xgobidata *xg = (xgobidata *) &xgobi;

  extern int moving_point;
  Boolean keep_going;

  /* tried out linking D_ij^mds_power to ||x_i-x_j||^mds_weightpow:
     worked on morsecodes for powers between 1 and 2, 
     but higher than 2 bombed, even 2 was highly artifactual
  */

  accum_dx = accum_dd = accum_xx = 0.0;
  nactive_distances = 0;

  if (!cleared_p) {    /* Clear out pos_grad. */
    clear_array(&pos_grad);
    cleared_p = True;
    make_empty_array(&pos_grad, pos.nrows, pos.ncols);
  }

/*
 * I can't tell how to handle this; I think because we're
 * rescaling in here somewhere, moving is getting messed up.
*/
  if (xg->is_point_moving && xg->nearest_point != -1) {
    k = xg->nearest_point;
    for (j=0; j<xg->ncols_used; j++)
      pos.data[k][j] = xg->raw_data[k][j] ;
  }

  nchanged = 0;

  /* Zero out the gradient matrix. */
  if (doit)
    zero_array(&pos_grad);

  /* Find new locations based on saved ones. */
  for (i = 0; i < pos.nrows; i++) {    /* For each node. */
    keep_going = True;

/* testing dfs September */
    if (xg->erased[i] == 1) {
      keep_going = False;
      for (k=0; k<mds_dims; k++)
        pos.data[i][k] = 0;
    }
    if (!keep_going)
      continue;  /* next i */

/* testing dfs September -- this is how you find out if a point is excluded */
    if (keep_going && xg->ncols == xg->ncols_used) {
      if (xg->clusv[(int)GROUPID(i)].excluded == 1) {
        keep_going = False;
        for (k=0; k<mds_dims; k++)
          pos.data[i][k] = 0;
      }
    }
    if (!keep_going)
      continue;  /* next i */

    for (j = 0; j < pos.nrows; j++) {  /* Check every pair ... */

      /* skip diagonal elements */
      if (i == j)
        ;

      /* skip erased points */
      else if(xg->erased[j] == 1)
        ;

      /* if the target distances are both missing, skip */
      else if (dist.data[i][j] == DBL_MAX)
        ;

      /*
       * if point i is selected for motion, skip over it:
       * that is, leave pos_grad[i][k] set to 0.
      */
      else if (xg->is_point_moving && xg->nearest_point != -1 &&
          moving_point == i)
        ;

      /*
       * if we're using groups, and these two points don't meet
       * our criteria, skip the pair
       */
      else if (mds_group_ind == within && !SAMEGLYPH(i,j))
        ;
      else if (mds_group_ind == between && SAMEGLYPH(i,j))
        ;
      else if (mds_group_ind == anchored && !CURRENTGLYPH(j))
        ;


      /*
       * If the target distance is within the thresholds
       * set using the barplot of distances, keep going.
       */
      else if (dist.data[i][j] < mds_threshold_low || 
           dist.data[i][j] > mds_threshold_high)
        ;
      else {
        nchanged++;  /* for threshold diagnostics, mostly */

        /* current pairwise distance */
        d = distance(pos.data[i], pos.data[j], mds_dims, mds_lnorm);

    /*
     * dist_goal could be precomputed whenever p changes and stored;
     * there's a lot of extra calculation going on here
     */

    /*
    printf("i=%d: %d  j=%d: %d \n", i, xg->erased[i], j, xg->erased[j]);
    */

        resid = 0.0;
        weight = 0.0;
        dist_goal = 0.0;

        if (dist.data[i][j] != DBL_MAX) {
          dist_goal = pow(dist.data[i][j], mds_power);
          if(dist.data[i][j] > 1E-5) 
            weight = pow(dist.data[i][j], mds_weightpow);
          else weight = pow(1E-5, mds_weightpow);
          resid = weight * (dist_goal - d * stress_dx / stress_xx);

          /* write to file if mds_once is called to generate shepard diagram */
          if(shepard) {
            fprintf(fpdat, "%5.5g %5.5g %5.5g %d %d\n",
            d, dist_goal, weight, i, j);
            fprintf(fprow, "%s,%s\n", xg->rowlab[ i ], xg->rowlab[ j ]);
          }
          nactive_distances++;
        }

        /*
         * If the residuals are both effectively equal to zero,
         * there's no point in working further with this pair of
         * points; move on.
        */
        if (fabs(resid) > delta) {

          /*
           * terms for stress function, actually squared correlation ...
           */
          accum_dx += dist_goal * d * weight;
          accum_xx += d * d * weight;
          accum_dd += dist_goal * dist_goal * weight;
          /*
           * Gradient ...
           */
          if (doit && fabs(d) > delta) {  /* can't have d = 0 */
            /*
             * All the variables in the numerator of this expression
             * were initialized to zero, so if there's a missing
             * distance, it should have no impact.
             */
            step_mag = resid * pow(d, 1.0 - mds_lnorm);
            /*
              printf("%d (%f %f) (%f %f) %f\n", i, dtmp, resid, 
              d, pow(d, mds_lnorm-1.0), step_mag);
            */

            for (k = 0; k < mds_dims; k++) {
              step_dir = step_mag *
                pow(fabs(pos.data[i][k]-pos.data[j][k]), mds_lnorm-1.0) *
                signum(pos.data[i][k]-pos.data[j][k]);
              pos_grad.data[i][k] += step_dir;
            } /* accumulate gradient */

          } /* distance non-zero */
        } /* residual non-zero */
      } /* target is within thresholds. */
    } /* j */
  } /* i */

  if (doit && nchanged > 0) {

    /* Find the gradient normalizing factor */
    gsum = psum = 0.0 ;
    for (i=0; i<pos_grad.nrows; i++) {
      f1 = normalize(pos_grad.data[i], mds_dims, mds_lnorm);
      f2 = normalize(pos.data[i], mds_dims, mds_lnorm);
      gsum += f1;
      psum += f2;
    }
    if (gsum < delta) gfactor = 0.0;
    else gfactor = mds_stepsize * psum/gsum;

    /* Add the gradient matrix to the position matrix */
    for (i=0; i<pos.nrows; i++) {
      for (k=0; k<mds_dims; k++) {
        pos.data[i][k] += (gfactor * pos_grad.data[i][k]);
      }
    }

/* testing July */
/*
  if (xg->is_point_moving && xg->nearest_point != -1 && moving_point != -1)
  {
    k = moving_point;
    for (j=0; j<xg->ncols_used; j++)
      pos.data[k][j] = xg->raw_data[k][j] ;
  }
*/

    /* center at variable means and scale globally to -1,+1 */
    scale_array_mean(&pos, pos.nrows, mds_dims);
  }

  stress_dx = accum_dx;  stress_dd = accum_dd;  stress_xx = accum_xx;
  if (stress_dd > delta && stress_xx > delta) {
    stress = pow( 1.0 - stress_dx * stress_dx / stress_xx / stress_dd, 0.5);
    add_stress_value(stress);
    draw_stress();

    /*
      printf("dims=%d lnorm=%f respow=%f gsum=%f psum=%f dx=%f dd=%f xx=%f\n", 
      mds_dims, mds_lnorm, mds_weightpow, 
      gsum, psum, 
      stress_dx, stress_dd, stress_xx
      );
    */

  }

  /* update Shepard labels */
  if(nactive_distances != prev_nactive_distances) {
    update_shepard_labels(nactive_distances);
    prev_nactive_distances = nactive_distances;
  }

}
