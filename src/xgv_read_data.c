
/*
 * Some routines for getting the data into the system.  ML, DFS 2/92.
 * Adapted from xgobi/read_array.c.
 */

#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"
#include "xgvis.h"

#define UNIFORM 0
#define NORMAL  1

/* use DBL_MAX to represent missing values */
#define NO_NAS    0
#define ALLOW_NAS 1

extern connect_lines *xg_lines;    /* Edges from xgvis.c. */
extern int xg_nlines;

/* Define. */
#define USAGE "xgvis basename \n\
Please supply the base name for a file (ending in .pos, .dist, .edges) \n\
or for a larger set of files. \n\
"

/* Functions. */
void clear_array(struct array *);
void scale_array_median(struct array *, int, int);
void scale_array_mean(struct array *, int, int);
void scale_array_max(struct array *, int, int);
void copy_array(struct array *, struct array *);
void make_empty_array(struct array *, int, int);
void do_exit(int);
void set_dist_matrix_from_edges(struct array *, struct array *, int);
void set_dist_matrix_from_pos(struct array *, struct array *, double);
void set_dist_matrix_from_pos_dot(struct array *, struct array *, int);
double vect_distance(double *, double *, int, int);
int read_labels(char *, char ***, int);
void symmetrize_dists(struct array *);
void report_asymmetry(struct array *);
double drandval(int);
void reset_data(void);

extern double distance(double *, double *, int, double);
extern void update_dissim_plot(void);

/* Macro. */
/* #define XtRealloc(x,y) (((x)==NULL)?XtMalloc((y)):XtRealloc((x),(y))) */


void
scale_array_max(struct array *arrp, int nr, int nc)
{
  extern double delta;  /* in mds.c */
  double max;
  int i, j;

  if (arrp->nrows < nr || arrp->ncols < nc)
    fprintf(stderr, "This array is smaller than nr or nc\n");
  else {

    max = 0.0;
    for (j=0; j<nc; j++) {
      for (i=0; i<nr; i++) {
        if (arrp->data[i][j] != DBL_MAX) {
          if (arrp->data[i][j] > max) max = arrp->data[i][j];
          if (arrp->data[i][j] < 0.0) 
            printf("Negative value %e in dist array at i=%d, j=%d\n",
            arrp->data[i][j], i, j);
        }
      }
    }

    if (max < 1E-10) 
      printf("Range of dist array too small: max=%e\n", max);

    for (j=0; j<nc; j++) {
      for (i=0; i<nr; i++) {
        if(arrp->data[i][j] != DBL_MAX)
          arrp->data[i][j] /= max;
      }
    }
  }
}


void
print_array(struct array *arrp, char *fname)
{
  int i, j;
  FILE *fp;

  if (fname == NULL)
    fp = stderr;
  else
    fp = fopen(fname, "w");

  for (i = 0; i < arrp->nrows; i++) {
    for (j = 0; j < arrp->ncols; j++) {
      fprintf(fp, "%f ", arrp->data[i][j]);
    }
    fprintf(fp, "\n");
  }
  fflush(fp);
}

/* General routine for reading an array from a file (given by name). */
/* Expects array structure to never have been allocated before. */
int
read_array_struct(char *data_in, struct array *arrp, Boolean allow_na)
{
  FILE *fp;
  char ch;
  double row1[10*NCOLS];
  char fname[100];
  int i, j, jrows, jcols, fs;
  int nblocks;
  int nitems;
  char word[64];
  int numna;

  numna = 0;

/*
 * Check file exists and open it - for stdin no open needs to be done
 * only assigning fp to be stdin.
*/
  if ((strcmp(data_in,"stdin") == 0) || (strcmp(data_in,"-") == 0))
    fp = stdin;
  else {
    if ((fp = fopen(data_in, "r")) == NULL) {
      strcpy(fname, data_in);
      strcat(fname, ".dat");
      if ((fp = fopen(fname, "r")) == NULL) {
        (void) printf("data file does not exist\n");
        return(0);
      }
    }
  }

/*
 * Read in the first row of the data file and calculate ncols.
*/
  arrp->ncols = 0;
  while ( (ch = getc(fp)) != '\n') {
    if (ch == '\t' || ch == ' ')
      ;
    else {
      if ( ungetc(ch, fp) == EOF || fscanf(fp, "%s", word) < 0 ) {

        /* fscanf(fp, "%lg", &row1[arrp->ncols++]) < 0 ) */

        perror("ungetc or fscanf");
        return(0);
      } else {
        if (allow_na && (strcasecmp(word, "na") == 0 || strcmp(word, ".") == 0))
          row1[arrp->ncols++] = (double) DBL_MAX;  /* math.h */
        else
          row1[arrp->ncols++] = atof(word);  /* returns double */

        if (arrp->ncols == 10*NCOLS) {
          fprintf(stderr,
            "Recompile: use a larger value of NCOLS in xgobi,\n");
          fprintf(stderr,
            "or a larger multiple of NCOLS in xgvis.\n");
          do_exit(-1);
        }

      }
    }
  }

/*
 * Allocate space for 500 rows.
*/
  arrp->data = (double **) 
    XtMalloc((unsigned int) 500 * sizeof(double *));
  for (i=0; i<500; i++)
    arrp->data[i] = (double *)
      XtMalloc((unsigned int) arrp->ncols * sizeof(double));
    
/*
 * Fill in the first row
*/
  for (j=0; j < arrp->ncols; j++)
    arrp->data[0][j] = row1[j];

/*
 * Read data, reallocating whenever 500 more rows are read in.
 * Determine nrows.
*/
  nblocks = 1;
  nitems = arrp->ncols;
  arrp->nrows = jrows = 1;
  jcols = 0;
  while (1) {
    /* fs = fscanf(fp, "%lg", &arrp->data[arrp->nrows][jcols]); */

    fs = fscanf(fp, "%s", word);
      
    if (fs == EOF)
      break;
    else if (fs < 0) {
      perror("fscanf in read_array");
      return(0);
    }
    else {
      if (allow_na && 
        (strcasecmp(word, "na") == 0 || 
         strcasecmp(word,"NA")  == 0 || 
         strcmp(word, ".") == 0))
      {
        arrp->data[arrp->nrows][jcols] = DBL_MAX; /* values.h */
        numna++;
      }
      else
        arrp->data[arrp->nrows][jcols] = atof(word);  /* returns double */

      jcols++;  /* moved from above */
      nitems++;

      if (jcols == arrp->ncols) {
        jcols = 0;
        arrp->nrows++;
        jrows++;
      }
      if (jrows == 500) {
        jrows = 0;
        nblocks++;
      
        arrp->data = (double **) XtRealloc((char *) arrp->data,
          (unsigned) (nblocks * 500) * sizeof(double *));
        for (i = 500*(nblocks-1); i<500*nblocks; i++)
          arrp->data[i] = (double *) XtMalloc((unsigned int)
            arrp->ncols * sizeof(double));
      }
    }
  }
/*
 * Close the data file
*/
  if (fclose(fp) == EOF)
    perror("read_array");

  if ( nitems != arrp->nrows * arrp->ncols ) {
    (void) fprintf(stderr,
    "read_array_struct: nrows*ncols != nitems, nrows=%d, ncols=%d, nitems=%d\n",
    arrp->nrows, arrp->ncols, nitems );
    return(0);
  }
  else {
    printf("read_array_struct: read nrows=%d, ncols=%d, nitems=%d, numna=%d\n",
      arrp->nrows, arrp->ncols, nitems, numna);

    return(1);
  }
}

/* Parses command line arguments, prepares the original data */
/* matrices.  Interprets defaults. */
void
initialize_data(int argc, char *argv[])
{
  int j, biggest_link;
  char fname[128];
  struct stat buf;
  Boolean found;
  extern double delta;  /* in mds.c */

  /* Clear things out to start. */
  clear_array(&dist_orig);
  clear_array(&dist);
  clear_array(&edges_orig);
  clear_array(&edges);
  clear_array(&pos_orig);
  clear_array(&pos);
  clear_array(&lines);
  linesptr = (struct array *) NULL;
  xg_lines = (connect_lines *) NULL;

  /*sprintf(xgv_basename, ""); *//* pointer, no space yet allocated */
  sprintf(glyphname, "");
  sprintf(pcolorname, "");
  sprintf(lcolorname, "");

/* parsing command line updated, May 99 */
  if (argc < 2) {
    fprintf(stderr, "Please supply a base filename.\n");
    exit(0);
  }

  for(j=0; argc>1; argc--,j++) {
    if (argv[j][0]=='-') {

      /*
       * -dims: specify the dimension of the embedding space
      */
      if (strcmp(argv[j], "-dims") == 0) {
        mds_dims = atoi(argv[j+1]);
        argc--; j++;
      }

      /*
       * -stepsize:  specify the stepsize on (0.,1.]
      */
      if (strcmp(argv[j], "-stepsize") == 0) {
        float ssize = (float) atof(argv[j+1]);
        if (ssize > 0. && ssize < 1.)
          mds_stepsize = ssize;
        argc--; j++;
      }

      /*
       * -mono:  force black and white display
      */
      else if (strcmp(argv[j], "-mono") == 0)
        mono = 1;
    }
  }
  xgv_basename = XtMalloc( (strlen(argv[j])+1) * sizeof(char));
  strcpy(xgv_basename, argv[j]);
/* */

  found = False;
  sprintf(fname, "%s.pos", xgv_basename);
  if (stat(fname, &buf) == 0) {
    found = True;
  } else {
    sprintf(fname, "%s.dat", xgv_basename);
    if (stat(fname, &buf) == 0) {
      found = True;
    }
  }
  if (found) {
    /* Reading position array. */
    if (pos_orig.nrows != 0)
      fprintf(stderr, "Redundant data-- ignoring positions from %s.\n",
        fname);
    else
      fprintf(stderr, "reading positions: %s\n", fname);
      read_array_struct(fname, &pos_orig, NO_NAS);
      scale_array_mean(&pos_orig, pos_orig.nrows, pos_orig.ncols);
  }

  sprintf(fname, "%s.edges", xgv_basename);
  if (stat(fname, &buf) == 0) {
    /* Reading edges array. */
    fprintf(stderr, "reading edges array: %s\n", fname);
    if (edges_orig.nrows != 0) 
      fprintf(stderr,
        "Redundant data-- ignoring edges from %s.\n", fname);
    else
      printf("reading .edges file:\n");
    read_array_struct(fname, &edges_orig, NO_NAS);
  }

  sprintf(fname, "%s.dist", xgv_basename);
  if (stat(fname, &buf) == 0) {
    /*
     * Reading distance array.
    */
    if (dist_orig.nrows != 0) 
      fprintf(stderr,
        "Redundant data-- ignoring dists from %s.\n", fname);
    else {
      fprintf(stderr, "reading distance array: %s\n", fname);
      read_array_struct(fname, &dist_orig, ALLOW_NAS);
      scale_array_max(&dist_orig, dist_orig.nrows, dist_orig.ncols);
    }
  }

/*
 * Sometimes one might want to have an edges file which is
 * used in determining inter-point distances, and sometimes
 * one is only interested in adding lines for decoration.
 * If there's only an edges file, it serves both purposes;
 * if there's both an edges file and a lines file, the edges
 * file is used for determining distances and the lines file
 * is shipped to xgobi for display.
*/
  sprintf(fname, "%s.lines", xgv_basename);
  if (stat(fname, &buf) == 0) {
    /* Reading lines array. */
    fprintf(stderr, "reading lines array: %s\n", fname);
    if (lines.nrows != 0) 
      fprintf(stderr,
        "Redundant data-- ignoring lines from %s.\n", fname);
    else
      read_array_struct(fname, &lines, NO_NAS);
  }

  if (lines.nrows > 0)
    linesptr = &lines;
  else if (edges_orig.nrows > 0)
    linesptr = &edges_orig;
  else
    linesptr = (struct array *) NULL;
  
  /* For later reference, let's compute the biggest link in the edges file. */
  biggest_link = -1;
  if (edges_orig.ncols >= 2)    /* Make sure dimensions are sensible. */
    for (j = 0; j < edges_orig.nrows; j++) {
      biggest_link = MAX(biggest_link, edges_orig.data[j][0]);
      biggest_link = MAX(biggest_link, edges_orig.data[j][1]);
    }

  /*
   * If there is a position matrix, scale the values overall
   * onto [-1, 1]
  */
  if (pos_orig.nrows > 0) {
    scale_array_mean(&pos_orig, pos_orig.nrows, pos_orig.ncols);
  }

/* If a distance matrix was passed in, let's symmetrize it. */
/*
  if (dist_orig.nrows != 0)
    symmetrize_dists(&dist_orig);
*/

  if (dist_orig.nrows != 0)
    report_asymmetry(&dist_orig);

  /* Default out any arguments the user didn't include. */
  if ((pos_orig.nrows == 0)&&(edges_orig.nrows == 0)&&(dist_orig.nrows == 0)) {
    fprintf(stderr, "No data given.\n");
    fprintf(stderr, USAGE);
    do_exit(-1);
  }
  else if ((pos_orig.nrows == 0)&&(edges_orig.nrows == 0)) {
    /* Use size of distance matrix to create the correct size pos. */
    make_empty_array(&pos_orig, dist_orig.nrows, 0);
  }
  else if ((pos_orig.nrows == 0)&&(dist_orig.nrows == 0)) {
    /* Use edges to create a distance matrix. */
    set_dist_matrix_from_edges(&dist_orig, &edges_orig, biggest_link);
    dist_type = LINK;

    /* Use size of distance matrix to create the correct size pos. */
    make_empty_array(&pos_orig, dist_orig.nrows, 0);
  }
  else if ((edges_orig.nrows == 0)&&(dist_orig.nrows == 0)) {
    /* Use full dimension positions to set distance matrix. */
    set_dist_matrix_from_pos(&dist_orig, &pos_orig, 2.0);
    dist_type = EUCLIDIAN;
  }
  else if (pos_orig.nrows == 0) {
    /* Create correct size position matrix from bounds of edges. */
    make_empty_array(&pos_orig, dist_orig.nrows, 0);
  }
  else if (edges_orig.nrows == 0) {
    /* Now problem. */
  }
  else if (dist_orig.nrows == 0) {
    /* Use edges to create a distance matrix. */
    set_dist_matrix_from_edges(&dist_orig, &edges_orig, pos_orig.nrows);
    dist_type = EUCLIDIAN;
  }

  /*print_array(&dist_orig, "foo");*/

  /* Other random sanity checking. */
  if (dist_orig.nrows != dist_orig.ncols) {
    fprintf(stderr, "Non square distance matrix!\n");
    do_exit(-1);
  }
  if (dist_orig.nrows != pos_orig.nrows) {
    fprintf(stderr,
      "Size of distance matrix different from number of positions!\n");
    do_exit(-1);
  }
  if ((edges_orig.nrows != 0) && (biggest_link > pos_orig.nrows)) {
    fprintf(stderr,
      "At least one link points off the end of the position matrix!\n");
    do_exit(-1);
  }
  if ((edges_orig.nrows != 0) &&
      (edges_orig.ncols != 2 && edges_orig.ncols != 3))
  {
    fprintf(stderr, "Edges file must have either 2 or 3 columns!\n");
    do_exit(-1);
  }
}

void
do_exit(int code)
{
  fprintf(stderr, "Exiting with code %d\n", code);
  exit(code);
}

/* Clear out an uninitialized array. */
void
clear_array(struct array *arrp)
{
  arrp->nrows = 0;
  arrp->ncols = 0;
  arrp->data = (double **) NULL;
}

/* Zero an array. */
void
zero_array(struct array *arrp)
{
  int i, j;
  for (i=0; i<arrp->nrows; i++) {
    for (j=0; j<arrp->ncols; j++) {
      arrp->data[i][j] = 0.0;
    }
  }
}

void
configure_pos_data(void)
{
  static double config_factor = 1.0;
  static double config_factor_prev = 1.0;
  double d, dsum = 0;
  int i, j, n;
  static int mds_dims_prev = 3;  /* initial value for mds_dims */
  
/*
 * Using pos_orig.nrows and mds_dims, find the average
 * distance between points.  Divide data in mds_dims columns
 * by that value -- having first multiplied by the previous
 * configuration factor?  Run whenever data or dimension is
 * reset: scrambled, reset, mds_ndims, mds_lnorm
*/

/* 
 * First remove the effect of the previous configuration factor,
 * using the previous value of mds_dims.
*/
  for (i = 0; i < pos.nrows; i++)
    for (j = 0; j < MIN(mds_dims_prev, pos.ncols); j++)
        pos.data[i][j] *= config_factor_prev;

/* Now find the new configuration factor, using current dims and lnorm */
  n = 0;
  for (i = 0; i < pos.nrows-1; i++) {  /* For each node. */
    for (j = i+1; j < pos.nrows; j++) {
      d = distance(pos.data[i], pos.data[j], mds_dims, mds_lnorm);
      dsum += d;
      n++;
    }
  }
  config_factor = dsum / (double) n;

/* Divide by the new configuration factor */
  for (i = 0; i < pos.nrows; i++) {
    for (j = 0; j < MIN(mds_dims, pos.ncols); j++) {
      pos.data[i][j] /= config_factor;
    }
  }

  config_factor_prev = config_factor;
  mds_dims_prev = mds_dims;
}

void
set_vgroups()
{
  int i, j;
  for (i=0; i<mds_dims; i++)
    xgobi.vgroup_ids[i] = 0;
  for (i=mds_dims, j=1; i<xgobi.ncols_used; i++, j++)
    xgobi.vgroup_ids[i] = j;
}

int
dcompare(const void *x1, const void *x2)
{
  int val = 0;
  double *d1 = (double *) x1;
  double *d2 = (double *) x2;

  if (*d1 < *d2)
    val = -1;
  else if (*d1 > *d2)
    val = 1;

  return(val);
}


void
set_distance_factor(void) {
  double d, dsum = 0.0;
  int i, j, count;
  static int initd = 0;
  extern Boolean dissim_initd;

  ndistances = dist.nrows * dist.ncols;
  if (!initd) {
    distance_vector = (double *) XtMalloc(ndistances * sizeof(double));
    distance_vector_sort = (double *) XtMalloc(ndistances * sizeof(double));
    initd = 1;
  }

  count = 0;
  for (i=0; i<dist.nrows; i++) {
    for (j=0; j<dist.ncols; j++) {
      if (i != j) {
        if (dist.data[i][j] == DBL_MAX) 
          ;
        else {
          d = pow(dist.data[i][j], mds_power);
          if (count == ndistances) {
            fprintf(stderr, "Trouble: distance_vector overfilled\n");
            break;
          }
          distance_vector[count] = 
            distance_vector_sort[count] = (double) d;
          count++;
          dsum += d;
        }
      }
    }
  }

  /* reset the (global) value of ndistances */
  ndistances = count;

  /* sort distance_vector */
  qsort((char *) distance_vector_sort, ndistances, sizeof(double), dcompare);

  distance_factor = (double) dsum / (double) ndistances;

/*
 * This value is now always 1.
 *
 *  max_distance = distance_vector_sort[ndistances - 1];
*/

/*
 * Update the dissimilarities/distances plot to reflect the
 * changes in the distance vector.
*/
  update_dissim_plot();
}

/* Copy original data into working buffers.  Make sure to leave space */
/* for 1:n. */
void
reset_data(void)
{
  int i, j;

  copy_array(&edges_orig, &edges);
  copy_array(&dist_orig, &dist);

  /* Copy positions truncating (or padding) to MAXDIMS. */
  if (pos.nrows == 0) {
    pos.nrows = pos_orig.nrows;
    pos.data = (double **)
      XtMalloc((unsigned int) pos_orig.nrows * sizeof(double *));
    for (i = 0; i < pos.nrows; i++)
      pos.data[i] = NULL;
  }

  if (pos_orig.nrows != pos.nrows) {
    /* Free the main spine. */
    fprintf(stderr,
      "uh oh, trouble.  copying from nrows %d to %d.  not implemented yet.\n",
      pos_orig.nrows, pos.nrows);
    do_exit(-1);
  }

  /* Copy each row. */
  for (i = 0; i < pos_orig.nrows; i++) {
    /* Make sure right number of columns. */
    if (pos.ncols != MAXDIMS+1)
      pos.data[i] = (double *)
        XtRealloc((char *) pos.data[i], 
        (unsigned int) (MAXDIMS+1) * sizeof(double));
    /* Copy the data. */
    for (j = 0; j < MIN(pos_orig.ncols, MAXDIMS); j++)
      pos.data[i][j] = pos_orig.data[i][j];

    /* Fill in the rest with random. */
    for (j = MIN(pos_orig.ncols, MAXDIMS); j < MAXDIMS; j++)
      pos.data[i][j] = drandval(UNIFORM);
  }

  pos.ncols = MAXDIMS;

  read_labels(xgv_basename, &rowlab, pos.nrows);

/*
 * Find the distance factor.  This sets (and then resets, if
 * necessary) ndistances, as well as building the distance_vector.
*/
  set_distance_factor();

/*
 * Adjust position data
*/
  configure_pos_data();

/*
 * Initialize the threshold value so that it's greater than
 * the largest distance -- which is now always 1.
*/
  mds_threshold_high = 2.;
  mds_threshold_low = 0.;

  /* Copy edges into xgobi lines structure. */
  if (linesptr != (struct array *) NULL && linesptr->nrows > 0) {
    xg_lines = (connect_lines *)
      XtRealloc((char *) xg_lines,
      (unsigned) linesptr->nrows*sizeof(connect_lines));
    xg_nlines = linesptr->nrows;
    for (i = 0; i < linesptr->nrows; i++) {
      xg_lines[i].a = linesptr->data[i][0];
      xg_lines[i].b = linesptr->data[i][1];
    }
  }

}

/* Copy an array into an existing structure. */
void
copy_array(struct array *from_arrp, struct array *to_arrp)
{
  int i, j;

  if (to_arrp->nrows == 0) {
    to_arrp->nrows = from_arrp->nrows;
    to_arrp->data = (double **)
      XtMalloc((unsigned int) to_arrp->nrows * sizeof(double *));
    for (i = 0; i < to_arrp->nrows; i++)    /* Zero things out. */
      to_arrp->data[i] = NULL;
  }

  if (to_arrp->nrows != from_arrp->nrows) {
    /* Free the main spine. */
    fprintf(stderr,
      "uh oh, trouble.  copying from nrows %d to %d.  not implemented yet.\n",
      to_arrp->nrows, from_arrp->nrows);
    do_exit(-1);
  }

  /* Copy each row. */
  for (i = 0; i < from_arrp->nrows; i++) {
    /* Make sure right number of columns. */
    if (from_arrp->ncols != to_arrp->ncols)
      to_arrp->data[i] = (double *)
        XtRealloc((char *) to_arrp->data[i], 
        (unsigned int) from_arrp->ncols * sizeof(double));
    /* Copy the data. */
    for (j = 0; j < from_arrp->ncols; j++)
      to_arrp->data[i][j] = from_arrp->data[i][j];
  }

  to_arrp->ncols = from_arrp->ncols;
}

void
scale_array_median(struct array *arrp, int nr, int nc)
{
  extern double delta;  /* in mds.c */
  double min, max, fac, median;
  int i, j;

  if (arrp->nrows < nr || arrp->ncols < nc)
    fprintf(stderr, "This array is smaller than nr or nc\n");
  else {

    /* First find overall max and min ... and sum */
    max = min = arrp->data[0][0] ;
    for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
        if (arrp->data[i][j] < min) min = arrp->data[i][j];
        if (arrp->data[i][j] > max) max = arrp->data[i][j];
      }
    }

    median = (max+min)/2.0 ;
    fac = max - median;
    /* printf("min %f max %f median %f fac %f\n", min, max, median, fac); */
    if (fabs(fac) < delta) {
      printf("potential trouble in scaling: median = %f, fac = %f\n",
        median, fac);
      printf("... resetting fac to 1.0\n");
      fac = 1;
    }

    for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
        /* center values at 0 */
        arrp->data[i][j] -= median;
        /* scale values onto [-1, 1] */  
        arrp->data[i][j] /= fac;  /* extreme danger of overflow! */
        if (arrp->data[i][j] >= DBL_MAX)
          arrp->data[i][j] = DBL_MAX/1000.;
      }
    }
  }
}

void
scale_array_mean(struct array *arrp, int nr, int nc)
{
  extern double delta;  /* in mds.c */
  double mean, dsum = 0.0;
  double max, min;
  int i, j;

  if (arrp->nrows < nr || arrp->ncols < nc)
    fprintf(stderr, "This array is smaller than nr or nc\n");
  else {

    for (j=0; j<nc; j++) {
      dsum = 0.0;
      for (i=0; i<nr; i++) dsum += arrp->data[i][j];
      mean = dsum / nr;
      for (i=0; i<nr; i++) arrp->data[i][j] -= mean;
    }

    max = min = arrp->data[0][0] ;
    for (j=0; j<nc; j++)
      for (i=0; i<nr; i++){
        if (arrp->data[i][j] < min) min = arrp->data[i][j];
        if (arrp->data[i][j] > max) max = arrp->data[i][j];
      }
    if((max-min)<1E-5) 
      printf("scale_array_mean: max-min too small = %e",max-min);
    for (j=0; j<nc; j++)
      for (i=0; i<nr; i++) 
        arrp->data[i][j] = (arrp->data[i][j] - min)/(max - min)*2 - 1.0;

  }
}


/* Fill in a possibly empty distance matrix by computing shortest */
/* paths (assume a single link is length 1).  Matrix will be nxn. */
void
set_dist_matrix_from_edges(struct array *dist_arrp, struct array *edges_arrp,
int n)
{
  int i, j;
  int infinity = 2*n;
  int changing;
  int end1, end2, end3;
  double d12;  /* weight */

  if (dist_arrp->nrows == 0)
    make_empty_array(dist_arrp, n, n);
  if ((dist_arrp->nrows != n)||(dist_arrp->ncols != n)) {
    fprintf(stderr, "Don't know how to change the size of distance matrix.\n");
    do_exit(-1);
  }

  /* Ok, we have a nice distance matrix, let's fill it in with infinity. */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      dist_arrp->data[i][j] = DOUBLE(infinity);
    dist_arrp->data[i][i] = DOUBLE(0);
  }

  /* As long as we find a shorter path using the edges, keep going. */
  changing = 1;
  while (changing) {
    changing = 0;
    for (i = 0; i < edges_arrp->nrows; i++) {
      end1 = edges_arrp->data[i][0]-1;
      end2 = edges_arrp->data[i][1]-1;
      d12 = (edges_arrp->ncols > 2) ? edges_arrp->data[i][2] : 1.0;
      for (end3 = 0; end3 < n; end3++) {
        /* So we have a direct link from end1 to end2.  Can this be */
        /* used to shortcut a path from end1 to end3 or end2 to end3? */
        if (dist_arrp->data[end1][end3] > d12 + dist_arrp->data[end2][end3]) {
          dist_arrp->data[end3][end1] =
            dist_arrp->data[end1][end3] =
            d12 + dist_arrp->data[end2][end3];
          changing = 1;
        }
        if (dist_arrp->data[end2][end3] > d12 + dist_arrp->data[end1][end3]) {
          dist_arrp->data[end3][end2] =
            dist_arrp->data[end2][end3] =
            d12 + dist_arrp->data[end1][end3];
          changing = 1;
        }
      }        /* end3 */
    }    /* end1 and end2 */
  }    /* while changing. */
  
  scale_array_max(dist_arrp, n, n);
}
  
void
set_dist_matrix_from_pos(struct array *dist_arrp, struct array *pos_arrp,
double l_norm)
{
  int n = pos_arrp->nrows;
  int i, j;
  double d;

  fprintf(stderr, "Creating %d^2 distance matrix.\n", n);

  if (dist_arrp->nrows == 0)
    make_empty_array(dist_arrp, n, n);
  if ((dist_arrp->nrows != n)||(dist_arrp->ncols != n)) {
    fprintf(stderr, "Don't know how to change the size of distance matrix.\n");
    do_exit(-1);
  }

  /* Ok, we have a nice distance matrix, let's fill it in with distances. */
  for (i = 0; i < n-1; i++) {
    dist_arrp->data[i][i] = 0.0;
    for (j = i+1; j < n; j++) {
      d = distance(pos_arrp->data[i], pos_arrp->data[j],
         pos_arrp->ncols, l_norm);
      dist_arrp->data[i][j] = d;
      dist_arrp->data[j][i] = d;
    }
    dist_arrp->data[n-1][n-1] = 0.0;
  }

  scale_array_max(dist_arrp, n, n);
}

void
set_dist_matrix_from_pos_dot(struct array *dist_arrp, struct array *pos_arrp,
int type)
{
  int n = pos_arrp->nrows;
  int i, j;
  double d;

  fprintf(stderr, "Creating %d^2 distance matrix (cos/dot).\n", n);

  if (dist_arrp->nrows == 0)
    make_empty_array(dist_arrp, n, n);
  if ((dist_arrp->nrows != n)||(dist_arrp->ncols != n)) {
    fprintf(stderr, "Don't know how to change the size of distance matrix.\n");
    do_exit(-1);
  }

  /* Ok, we have a nice distance matrix, let's fill it in with distances. */
  for (i = 0; i < n-1; i++) {
    dist_arrp->data[i][i] = 0.0;
    for (j = i+1; j < n; j++) {
      d = vect_distance(pos_arrp->data[i], pos_arrp->data[j],
          pos_arrp->ncols, type);
      dist_arrp->data[i][j] = d;
      dist_arrp->data[j][i] = d;
    }
    dist_arrp->data[n-1][n-1] = 0.0;
  }

  scale_array_max(dist_arrp, n, n);
}

double
vect_distance (double *v1, double *v2, int n, int type)
  /* n = Number of dimensions. */
{
  double self1, self2, dot, denom;
  register int i;

  for (i = dot = self1 = self2 = 0; i < n; i++) {
    self1 += v1 [i] * v1 [i];
    self2 += v2 [i] * v2 [i];
    dot   += v1 [i] * v2 [i];
  }
  
  if (type == DOTPROD)    /* Dot product. */
    return (dot);
  /* Else, cosine. */

  if ((denom = sqrt (self1 * self2)) == 0)
    return (0.0);
  
  return (dot / denom);
}

/* Make an empty array. */
void
make_empty_array(struct array *arrp, int nr, int nc)
{
  int i;

  if ((arrp->nrows != 0)||(arrp->ncols != 0)) {
    fprintf(stderr, "Asked to allocate over a filled array, not yet.\n");
    do_exit(-1);
  }

  arrp->data = (double **) 
    XtMalloc((unsigned int) nr * sizeof(double *));
  for (i = 0; i < nr; i++)
    arrp->data[i] = (double *)
      XtMalloc((unsigned int) nc * sizeof(double));
  arrp->nrows = nr;
  arrp->ncols = nc;
}

/* Stolen from xgobi read_rowlabels. */
int
read_labels(char *rootname, char ***rowlabp, int nrows)
{
  int i, ncase;
  FILE *fp;
  char fname[128];
  Boolean found = False;
  static char *suffix[] = {
    ".row", ".rowlab", ".case", ".labels"
  };


  (*rowlabp) = (char **) XtMalloc((unsigned) nrows * sizeof (char *));
  for (i = 0; i < nrows; i++)
    (*rowlabp)[i] = (char *) XtMalloc((unsigned) ROWLABLEN * sizeof(char));

  i = 0;
  while (found == False && i<4 ) { 
    (void) sprintf(fname, "%s%s", rootname, suffix[i++]);
    if ( (fp = fopen(fname, "r")) != NULL) {
      found = True;
      break;
    }
  }

  if (found) {
    /*
     * Read in case labels or initiate them to generic if no label
     * file exists
     */
    int len;
    ncase = 0;

    while (fgets((*rowlabp)[ncase], ROWLABLEN, fp) != NULL) {

      /*
       * Get rid of terminating newline character, if one exists.
       */
      len = strlen((*rowlabp)[ncase])-1;
      if ((*rowlabp)[ncase][len] == '\n')
        (*rowlabp)[ncase][len] = '\0';

      ncase++;
      if (ncase >= nrows)
        break;
    }
    if (ncase != nrows) {
      (void) printf("number of labels = %d, number of rows = %d\n",
         ncase, nrows);
      for (i = ncase; i < nrows; i++)
        (void) sprintf((*rowlabp)[ncase], "%d", i+1);
    }
  }
  /* Generic labels. */
  else {
    for (i = 0; i < nrows; i++)
      (void) sprintf((*rowlabp)[i], "%d", i+1);
  }
  return(1);
}

void
symmetrize_dists(struct array *dist_arrp)
{
  int i, j;
  int n = dist_arrp->nrows;
  double mean;

  if (n != dist_arrp->ncols) {
    fprintf(stderr, "Non square distance matrix!\n");
    do_exit(-1);
  }

  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      mean = (dist_arrp->data[i][j]+dist_arrp->data[j][i])/2;
      dist_arrp->data[i][j] = dist_arrp->data[j][i] = mean;
    }
  }
}

void
report_asymmetry(struct array *dist_arrp)
{
  int i, j;
  int n = dist_arrp->nrows;

  if (n != dist_arrp->ncols) {
    fprintf(stderr, "Non square distance matrix!\n");
    do_exit(-1);
  }

  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      if (dist_arrp->data[i][j] != dist_arrp->data[j][i]) {
        fprintf(stderr, "The input distance matrix is asymmetric.\n");
        return;
      }
    }
  }
}

double
drandval(int type)
{
/*
 * generate a random value on approximately -1, 1
*/
  double drand;
  static double dsave;
  static Boolean isave = false;

  if (type == UNIFORM) {
      drand = randvalue();

      /*
       * Center and scale to [-1, 1]
      */
      drand = (drand - .5) * 2;

  } else if (type == NORMAL) {

    Boolean check = true;
    double d, dfac;

    if (isave) {
      isave = false;
      /* prepare to return the previously saved value */
      drand = dsave;
    } else {
      isave = true;
      while (check) {

        rnorm2(&drand, &dsave);
        d = drand*drand + dsave*dsave;

        if (d < 1.0) {
          check = false;
          dfac = sqrt(-2. * log(d)/d);
          drand = drand * dfac;
          dsave = dsave * dfac;
        }
      } /* end while */
    } /* end else */

    /*
     * Already centered; scale to approximately [-1, 1]
    */
    drand = (drand / 3.0);
  }
  return(drand);
}

void
scramble_data(void) {
  int i, j, k;
/*
 * Fill in all data with normal random values ...
*/

  for (i = 0; i < pos_orig.nrows; i++) {
    /* for (j = 0; j < MAXDIMS; j++, k++) { */
    for (j = 0; j < mds_dims; j++, k++) {
      pos.data[i][j] = drandval(UNIFORM);  /* on [-1,1] */
    }
  }
}


#ifdef COMMENT

   spring_once() using edges instead of g.

   on reset_data, create list of connected vertices (for spring)

#endif