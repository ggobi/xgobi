/* call:
  ./xgvis /usr/andreas/XGVIS/MORSE/morsecodes &
  ./xgvis /usr/andreas/XGVIS/GRAPHS/perm5 &
 */

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
extern void scale_array_mean(struct array *, int, int, xgobidata);
extern void update_shepard_labels(int);

extern Widget metric_cmd[2];

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


void
scale_array_meanx(struct array *arrp, int nr, int nc, int mp)
{
  extern double delta;  /* in mds.c */
  double mean, dsum = 0.0;
  double max, min, scl;
  int i, j, ni;

  if (arrp->nrows < nr || arrp->ncols < nc)
    fprintf(stderr, "This array is smaller than nr or nc\n");
  else {

    for (j=0; j<nc; j++) {
      dsum = 0.0; ni = 0;
      for (i=0; i<nr; i++) { dsum += arrp->data[i][j]; ni++; }
      mean = dsum / ni;
      for (i=0; i<nr; i++) arrp->data[i][j] -= mean;
    }

    max = -1000000; min = 1000000;
    for (j=0; j<nc; j++)
      for (i=0; i<nr; i++){
        if (arrp->data[i][j] < min) min = arrp->data[i][j];
        if (arrp->data[i][j] > max) max = arrp->data[i][j];
      }
    if((max-min)<1E-5) 
      printf("scale_array_mean: max-min too small = %e",max-min);
    if(max > -min) scl = max; else scl = -min;
    for (j=0; j<nc; j++)
      for (i=0; i<nr; i++) 
        arrp->data[i][j] = arrp->data[i][j]/scl;

  }
}


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

static double **tmpVector;
static double *tmpMisVector;

int realCompare(const void* aPtr, const void* bPtr)
{
  int aIndex = *(int*)aPtr;
  int bIndex = *(int*)bPtr;

  double aReal = tmpVector[aIndex/dist.ncols][aIndex%dist.ncols];
  double bReal = tmpVector[bIndex/dist.ncols][bIndex%dist.ncols];
  if ((tmpMisVector[aIndex] == DBL_MAX) && (tmpMisVector[bIndex] == DBL_MAX))  return 0;
  if (tmpMisVector[aIndex] == DBL_MAX) return 1;
  if (tmpMisVector[bIndex] == DBL_MAX) return -1;

  if (aReal < bReal) return -1;
  else if (aReal == bReal) return 0;
  else return 1; 
}

/*
 * Perform one loop of the iterative mds function.
 *
 * If doit is False, then we really want to determine the
 * stress function without doing anything to the gradient
*/
void
mds_once(Boolean doit, Boolean shepard, FILE* fpdat, FILE* fprow)
{
  int i, j, k, ii;
  static struct array pos_grad;
  static Boolean cleared_p = False;
  static int nactive_distances = -1, prev_nactive_distances = -1;

  static double *trans_dist = NULL;
  static int *trans_dist_index = NULL;
  static int *b = NULL; /* blocklength for isotonic regression */

  static int sortnecessary = 0;
  short notfinished;
  short stop;

  double d, sum;
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
  int mp;

  /* May 16 2001: preparation for isotonic regression */
  if (trans_dist == NULL){
    trans_dist = (double *) XtMalloc(dist.nrows * dist.ncols * sizeof(double));
    trans_dist_index = (int *) XtMalloc(dist.nrows * dist.ncols * sizeof(int));
    b = (int *) XtMalloc(dist.nrows * dist.ncols * sizeof(int));
  }
  for (i = 0 ; i < dist.nrows; i++) {
    for (j = 0; j < dist.ncols; j++) {
      trans_dist[i*dist.ncols+j] = DBL_MAX;
      /*      trans_dist_index[i*dist.ncols+j] = i*dist.ncols+j;*/
    } 
  }
  

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
  mp = -1;
  if (xg->is_point_moving && xg->nearest_point != -1) {
    mp = xg->nearest_point;
    for (j=0; j<xg->ncols_used; j++)
      pos.data[mp][j] = xg->raw_data[mp][j] ;
  }

  nchanged = 0;

  /* Zero out the gradient matrix. */
  if (doit)
    zero_array(&pos_grad);

  /* Find new locations based on saved ones. */
  for (i = 0; i < pos.nrows; i++) {    /* For each node. */

/* testing dfs September */
    if (xg->erased[i] == 1) continue;

/* testing dfs September -- this is how you find out if a point is excluded */
    if (xg->ncols == xg->ncols_used) 
      if (xg->clusv[(int)GROUPID(i)].excluded == 1) 
        continue;

    for (j = 0; j < pos.nrows; j++) {  /* Check every pair ... */

      /* skip diagonal elements */
      if (i == j) continue; 

      /* skip erased points */
      if(xg->erased[j] == 1) continue;

      if (xg->ncols == xg->ncols_used)
        if (xg->clusv[(int)GROUPID(j)].excluded == 1) 
          continue;

      /* if the target distances are both missing, skip */
      if (dist.data[i][j] == DBL_MAX) continue;

      /*
       * if point i is selected for motion, skip over it:
       * that is, leave pos_grad[i][k] set to 0.
      */
      /* AB try applying the gradient also to the moving point...
      else if (xg->is_point_moving && xg->nearest_point != -1 &&
          moving_point == i)
        ;
      */

      /*
       * if we're using groups, and these two points don't meet
       * our criteria, skip the pair
       */
      if (mds_group_ind == within && !SAMEGLYPH(i,j)) continue;
      if (mds_group_ind == between && SAMEGLYPH(i,j)) continue;
      if (mds_group_ind == anchored && !CURRENTGLYPH(j)) continue;

      /*
       * If the target distance is within the thresholds
       * set using the barplot of distances, keep going.
       */
      if (dist.data[i][j] < mds_threshold_low || 
          dist.data[i][j] > mds_threshold_high) continue;


        nchanged++;  /* for threshold diagnostics, mostly */

        /* current pairwise distance */
        trans_dist[i*dist.ncols+j] = distance(pos.data[i], pos.data[j], mds_dims, mds_lnorm);

    } /* j */
  } /* i */

  if (nchanged > 0) {
    if (scaling_method == NONMETRIC) {
      /*      printf("before sorting:\n");
       for(i = 0; i <dist.nrows; i++){
        for (j = 0; j < dist.ncols; j++) {
            printf ("%.3f\t",dist.data[trans_dist_index[i*dist.ncols+j]/dist.ncols][trans_dist_index[i*dist.ncols+j]%dist.ncols]);
        }
        printf("\n");
        }*/
      

      tmpVector = dist.data;
      tmpMisVector = trans_dist;
      if (sortnecessary != nchanged) {
        /*printf ("sorting data ...");*/
        for (i = 0 ; i < dist.nrows; i++) {
          for (j = 0; j < dist.ncols; j++) {
            trans_dist_index[i*dist.ncols+j] = i*dist.ncols+j;
          } 
        }
        Myqsort(trans_dist_index, dist.ncols*dist.nrows, sizeof(int),
              realCompare);
        sortnecessary = nchanged;
        /*printf ("done\n");*/
      }

       /* start isotonic regression */
      /* initialize block length */
      i = 0;
      for ( ; i < dist.nrows*dist.ncols; i++) {
        if (trans_dist[trans_dist_index[i]] == DBL_MAX) continue;
        ii = i+1;
        while ((ii <  dist.nrows*dist.ncols)&&((trans_dist[trans_dist_index[ii]] == DBL_MAX) || 
                (dist.data[trans_dist_index[ii]/dist.ncols][trans_dist_index[ii]%dist.ncols]
                 ==dist.data[trans_dist_index[i]/dist.ncols][trans_dist_index[i]%dist.ncols]))){
          ii++;
        }
        /* ii points to start of the next block */
        b[i] = ii-i;
        sum = 0;

        for (k = i; k < ii; k++) {
          if (trans_dist[trans_dist_index[k]] != DBL_MAX) {
            sum += trans_dist[trans_dist_index[k]];
          }
        }

        if (b[i] > 1) {
          sum /= b[i];
          for (k = i; k < ii; k++) {
            if (trans_dist[trans_dist_index[k]] != DBL_MAX)
              trans_dist[trans_dist_index[k]] = sum;
          }
          i += b[i]-1;
        }
      }

      notfinished = 1;
      while (notfinished > 0) {
        notfinished = 0;
        i = 0;

        stop = 0;
        while (!stop) {
          ii = i+b[i];
          if (ii < dist.nrows*dist.ncols) {
            if (trans_dist[trans_dist_index[i]] > 
                trans_dist[trans_dist_index[ii]])
            {
              notfinished = 1;
              trans_dist[trans_dist_index[i]] = 
                (trans_dist[trans_dist_index[i]]*b[i] + 
                 trans_dist[trans_dist_index[ii]]*b[ii])/(b[i] + b[ii]);

              b[i] += b[ii];
            } else {
              i += b[i];
            }
          } else stop = 1;
        }
      }
 
      for (i=0; i < dist.nrows*dist.ncols; i++) {
        for (j = i+1; j < i+b[i]; j++) {
          trans_dist[trans_dist_index[j]] = trans_dist[trans_dist_index[i]];
        }
        i += b[i]-1;
      }
      
      for (i = 0; i < dist.nrows; i++) {
        for (j = 0; j < dist.ncols; j++) {
          if (trans_dist[i*dist.ncols+j] != trans_dist[j*dist.ncols+i])
            printf ("assymetric in f(d): %d %d \n",i,j);
        }
      }
    }

    for (i = 0; i < pos.nrows; i++) {
      for (j = 0; j < pos.nrows; j++) {
        if (trans_dist[i*dist.ncols+j]  ==  DBL_MAX) continue;
        d = distance(pos.data[i], pos.data[j], mds_dims, mds_lnorm);

        resid = 0.0;
        weight = 0.0;
        dist_goal = 0.0;

        if (dist.data[i][j] != DBL_MAX) {
          if (scaling_method == NONMETRIC) dist_goal = trans_dist[i*dist.ncols+j];
          else dist_goal = pow(dist.data[i][j], mds_power);
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
      } /* loop over j end */
    } /* loop over i end */
  }


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
    /* AB scaling gets in the way somewhere... */
    scale_array_meanx(&pos, pos.nrows, mds_dims, mp);
    /*    */
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
