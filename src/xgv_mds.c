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
extern void update_shepard_labels(int);

extern Widget metric_cmd[2];

/* Macros. */
#define DOUBLE(x) ((double)(x))
#define FLOAT(x) ((float)(x))

#define SAMEGLYPH(i,j) \
( xg->color_now[(i)] == xg->color_now[(j)] && \
  xg->glyph_now[(i)].type == xg->glyph_now[(j)].type && \
  xg->glyph_now[(i)].size == xg->glyph_now[(j)].size ) \

#define CURRENTGLYPH(i) \
( xg->color_now[(i)] == xg->color_now[point_midbutton] && \
  xg->glyph_now[(i)].type == xg->glyph_now[point_midbutton].type && \
  xg->glyph_now[(i)].size == xg->glyph_now[point_midbutton].size ) \

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
  int i, j, k, ii, n;
  static struct array pos_grad;
  static Boolean cleared_p = False;
  static int prev_active_dist = -1,  prev_nonmetric_active_dist = -1;

  static double *trans_dist = NULL;
  static int *trans_dist_index = NULL;
  static int *bl = NULL; /* blocklength for isotonic regression */

  short notfinished;
  short stop;

  double d, sum;
  double step_mag, step_dir;
  double gsum, psum, gfactor, tsum, tmp;
  double dist_goal;
  double resid;
  double weight;
  double stress, stress_dx = 1.0, stress_dd = 1.0, stress_xx = 1.0;
  double accum_dx, accum_dd, accum_xx;
  double f1, f2;
  xgobidata *xg = (xgobidata *) &xgobi;
  extern int moving_point;
  extern int move_type;
  extern int point_midbutton;

  /* tried out linking D_ij^mds_power to ||x_i-x_j||^mds_weightpow:
     worked on morsecodes for powers between 1 and 2, 
     but higher than 2 bombed, even 2 was highly artifactual  -- AB
  */

  /* components of stress */
  accum_dx = accum_dd = accum_xx = 0.0;

  if (!cleared_p) {    /* Clear out pos_grad. */
    clear_array(&pos_grad);
    cleared_p = True;
    make_empty_array(&pos_grad, pos.nrows, pos.ncols);
  }

  /* count active dissimilarities/distances */
  num_active_dist = 0;

  /* Zero out the gradient matrix. */
  if (doit)
    zero_array(&pos_grad);

  /* preparation for isotonic regression */
  if (trans_dist == NULL){
    trans_dist = (double *) XtMalloc(dist.nrows * dist.ncols * sizeof(double));
    trans_dist_index = (int *) XtMalloc(dist.nrows * dist.ncols * sizeof(int));
    bl = (int *) XtMalloc(dist.nrows * dist.ncols * sizeof(int));  /* block lengths */
  }
  for (i = 0 ; i < dist.nrows; i++) {
    for (j = 0; j < dist.ncols; j++) {
      trans_dist[i*dist.ncols+j] = DBL_MAX;
    } 
  }


  /* -------------------- collect active dissimilarities --------------- */
  for (i = 0; i < pos.nrows; i++) {    /* For each node. */

    /* skip erased points */
    if (xg->erased[i] == 1) continue;

    /* skip excluded points */
    if (xg->ncols == xg->ncols_used) 
      if (xg->clusv[(int)GROUPID(i)].excluded == 1) 
        continue;

    /* AB: do not exclude moving i: 
     * in nonmetric MDS it matters what the set of distances is! 
     */

    for (j = 0; j < pos.nrows; j++) {  /* Check every pair ... */

      /* skip diagonal elements */
      if (i == j) continue; 

      /* skip erased points */
      if(xg->erased[j] == 1) continue;

      /* AB: what is this? */
      if (xg->ncols == xg->ncols_used)
        if (xg->clusv[(int)GROUPID(j)].excluded == 1) 
          continue;

      /* if the target distances are both missing, skip */
      if (dist.data[i][j] == DBL_MAX) continue;

      /*
       * if we're using groups, and these two points don't meet
       * our criteria, skip the pair
       */
      if (mds_group_ind == within && !SAMEGLYPH(i,j)) continue;
      if (mds_group_ind == between && SAMEGLYPH(i,j)) continue;
      if (mds_group_ind == anchorscales && !CURRENTGLYPH(j)) continue;
      if (mds_group_ind == anchorfixed && (CURRENTGLYPH(i) || !CURRENTGLYPH(j))) continue;

      /*
       * If the target distance is within the thresholds
       * set using the barplot of distances, keep going.
       */
      if (dist.data[i][j] < mds_threshold_low || 
          dist.data[i][j] > mds_threshold_high) continue;

      /* another active dissimilarity */
      num_active_dist++;  

      /* current pairwise distance */
      trans_dist[i*dist.ncols+j] = distance(pos.data[i], pos.data[j], mds_dims, mds_lnorm);

    } /* j */
  } /* i */
  /* -------------------- end collecting active dissimilarities ------------------ */

  /* --------------- for active dissimilarities, do the work --------------------- */ 
  if (num_active_dist > 0) {

    /* ----------------begin isotonic regression ----------- */
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
      /* need to sort? */
      if (num_active_dist != prev_nonmetric_active_dist) {
        /*printf ("sorting data ...");*/
        for (i = 0 ; i < dist.nrows; i++) {
          for (j = 0; j < dist.ncols; j++) {
            trans_dist_index[i*dist.ncols+j] = i*dist.ncols+j;
          } 
        }
        Myqsort(trans_dist_index, dist.ncols*dist.nrows, sizeof(int),
              realCompare);
        /*printf ("done\n");*/
      }
      prev_nonmetric_active_dist = num_active_dist;

      /* initialize block length */
      i = 0;
      for ( ; i < dist.nrows*dist.ncols; i++) {
        if (trans_dist[trans_dist_index[i]] == DBL_MAX) continue;
        ii = i+1;
        while ((ii < dist.nrows*dist.ncols) && 
	       ((trans_dist[trans_dist_index[ii]] == DBL_MAX) || 
                (dist.data[trans_dist_index[ii] / dist.ncols][trans_dist_index[ii] % dist.ncols]
                 == dist.data[trans_dist_index[i] / dist.ncols][trans_dist_index[i] % dist.ncols]))){
          ii++;
        }
        /* ii points to start of the next block */
        bl[i] = ii-i;
        sum = 0;

        for (k = i; k < ii; k++) {
          if (trans_dist[trans_dist_index[k]] != DBL_MAX) {
            sum += trans_dist[trans_dist_index[k]];
          }
        }

        if (bl[i] > 1) {
          sum /= bl[i];
          for (k = i; k < ii; k++) {
            if (trans_dist[trans_dist_index[k]] != DBL_MAX)
              trans_dist[trans_dist_index[k]] = sum;
          }
          i += bl[i]-1;
        }
      }

      notfinished = 1;
      while (notfinished > 0) {
        notfinished = 0;
        i = 0;

        stop = 0;
        while (!stop) {
          ii = i+bl[i];
          if (ii < dist.nrows*dist.ncols) {
            if (trans_dist[trans_dist_index[i]] > 
                trans_dist[trans_dist_index[ii]])
            {
              notfinished = 1;
              trans_dist[trans_dist_index[i]] = 
                (trans_dist[trans_dist_index[i]]*bl[i] + 
                 trans_dist[trans_dist_index[ii]]*bl[ii])/(bl[i] + bl[ii]);

              bl[i] += bl[ii];
            } else {
              i += bl[i];
            }
          } else stop = 1;
        }
      }

      for (i=0; i < dist.nrows*dist.ncols; i++) {
        for (j = i+1; j < i+bl[i]; j++) {
          trans_dist[trans_dist_index[j]] = trans_dist[trans_dist_index[i]];
        }
        i += bl[i]-1;
      }

      for (i = 0; i < dist.nrows; i++) {
        for (j = 0; j < dist.ncols; j++) {
          if (trans_dist[i*dist.ncols+j] != trans_dist[j*dist.ncols+i])
            printf ("asymmetric in f(d): %d %d \n",i,j);
        }
      }

      /* rescale trans_dist to same sum of squares as dist, so the config does not shrink */
      tsum = 0.0;
      for(n=0; n < num_active_dist; n++) {tmp = trans_dist[trans_dist_index[n]]; tsum += tmp*tmp; }
      tsum = pow(tsum/num_active_dist, 0.5);
      for(n=0; n < num_active_dist; n++) trans_dist[trans_dist_index[n]] /= tsum;

    } /* ------------------ end isotonic regression ------------- */


    /* ------------------ gradient descent: j's push i's ----------- */
    for (i = 0; i < pos.nrows; i++) {
      for (j = 0; j < pos.nrows; j++) {
        if (trans_dist[i*dist.ncols+j]  ==  DBL_MAX) continue;

        d = distance(pos.data[i], pos.data[j], mds_dims, mds_lnorm);

        resid = 0.0;
        weight = 0.0;
        dist_goal = 0.0;

	if (scaling_method == NONMETRIC) dist_goal = trans_dist[i*dist.ncols+j];
	else dist_goal = pow(dist.data[i][j], mds_power);

	if(dist.data[i][j] > 1E-5) weight = pow(dist.data[i][j], mds_weightpow);
	else weight = pow(1E-5, mds_weightpow);

	resid = weight * (dist_goal - d * stress_dx / stress_xx);

          /* write to file if mds_once is called to generate shepard diagram */
	if(shepard) {
	  fprintf(fpdat, "%5.5g %5.5g %5.5g %d %d\n",
		  d, dist_goal, weight, i, j);
	  fprintf(fprow, "%s,%s\n", xg->rowlab[ i ], xg->rowlab[ j ]);
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

	    /* accumulate gradient */
            for (k = 0; k < mds_dims; k++) {
              step_dir = step_mag *
                pow(fabs(pos.data[i][k]-pos.data[j][k]), mds_lnorm-1.0) *
                signum(pos.data[i][k]-pos.data[j][k]);
              pos_grad.data[i][k] += step_dir;
            }

          } /* distance non-zero */
        } /* residual non-zero */
      } /* loop over j end */
    } /* loop over i end */
  } /* num_active_dist */

  /* ------------------------- end active dissimilarities ------------------- */

  if (doit && num_active_dist > 0) {

    /* Find the gradient normalizing factor */
    gsum = psum = 0.0 ;
    for (i=0; i<pos_grad.nrows; i++) {
      f1 = normalize(pos_grad.data[i], mds_dims, mds_lnorm);
      f2 = normalize(pos.data[i], mds_dims, mds_lnorm);
      gsum += f1;
      psum += f2;
    }
    if (gsum < delta) gfactor = 0.0;
    else gfactor = 0.1 * mds_stepsize * psum/gsum;

    /* Add the gradient matrix to the position matrix */
    for (i=0; i<pos.nrows; i++) {
      for (k=0; k<mds_dims; k++) {
        pos.data[i][k] += (gfactor * pos_grad.data[i][k]);
      }
    }

  } /* end  if (doit && num_active_dist > 0) { */

/*
 * set moving points where the mouse moved them:
*/
  if (xg->is_point_moving && moving_point != -1) {
    if(move_type==0) {
      for (k=0; k < mds_dims; k++) {
	pos.data[moving_point][k] = xg->raw_data[moving_point][k] ;
    }}
    if(move_type==1) {
      for (i=0; i<xg->nrows_in_plot; i++) {
	n = xg->rows_in_plot[i];
	if (!xg->erased[n] && SAMEGLYPH(n,moving_point)) {
	  for (k=0; k < mds_dims; k++) {
	    pos.data[n][k] = xg->raw_data[n][k] ;
    }}}}
    if(move_type==2) {
      for (i=0; i<xg->nrows_in_plot; i++) {
	n = xg->rows_in_plot[i];
	if (!xg->erased[n]) {
	  for (k=0; k < mds_dims; k++) {
	    pos.data[n][k] = xg->raw_data[n][k] ;
    }}}}
  }  /* end if (xg->is_point_moving && moving_point != -1) { */

  /* calculate stress and draw it */
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
  if(num_active_dist != prev_active_dist) {
    update_shepard_labels(num_active_dist);
    prev_active_dist = num_active_dist;
  }

}

/* ---------------------------------------------------------------- */
