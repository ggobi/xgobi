/* read_data.c */
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

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"

static Widget epopup = (Widget) NULL;
static Widget newfname_txt;

#define NCHOICES 12
Widget read_choice[NCHOICES];
#define PTCOLOR_TGL   read_choice[0]
#define PTGLYPH_TGL   read_choice[1]
#define PTERASE_TGL   read_choice[2]
#define ROWLABEL_TGL  read_choice[3]
#define COLLABEL_TGL  read_choice[4]
#define LINE_TGL      read_choice[5]
#define LINECOLOR_TGL read_choice[6]
#define VGROUPS_TGL   read_choice[7]
#define RGROUPS_TGL   read_choice[8]
#define JITTER_TGL    read_choice[9]
#define NLINK_TGL     read_choice[10]
#define NDATA_TGL     read_choice[11]  /* nrows x ncols */

#define INITSTRSIZE 512

/* ARGSUSED */
static XtCallbackProc
read_toggle_cback(Widget w, xgobidata *xg, XtPointer cback_data)
{
  Boolean state;
  XtVaGetValues(w, XtNstate, &state, NULL);
  setToggleBitmap(w, state);
}

/***********************************************************/
/*********************** for scatmat ***********************/
/***********************************************************/

Boolean
read_collab_to_row(char *data_in, Boolean init, xgobidata *xg)
{
  char lab_file[128];
  static char *suffix[] = {
    ".col", ".column", ".collab", ".var"
  };
  char initstr[INITSTRSIZE];
  int i, nvar = 0;
  Boolean found = False;
  FILE *fp;

  /* Now read in the variable labels to put in the last
       ncols rows, to label the diagonal */
  /*
   * Check if variable label file exists, and open if so.
  */
  if (data_in != NULL && data_in != "" && strcmp(data_in,"stdin") != 0) {
    i = 0;
    while (found == False && i<4 ) {
      (void) strcpy(lab_file, data_in);
      (void) strcat(lab_file, suffix[i++]);
      if ( (fp = fopen(lab_file,"r")) != NULL)
        found = True;
    }
  }
 
  /*
   * Read in variable labels or initiate them to generic if no label
   * file exists
  */
  if (found)
  {
    int len;
    nvar = 0;
  
    while (fgets(initstr, INITSTRSIZE-1, fp) != NULL)
    {
      len = MIN(INT(strlen(initstr)), COLLABLEN-1) ;
      while (initstr[len-1] == '\n' || initstr[len-1] == ' ')
        len-- ;
      strncpy(xg->rowlab[xg->nrows-xg->sm_ncols+nvar], initstr, len) ;
      xg->rowlab[xg->nrows-xg->sm_ncols+nvar][len] = '\0' ;
 
      nvar++;
      if (nvar >= xg->sm_ncols)
        break;
    }
  
    if (init && nvar != xg->sm_ncols) {
      (void) fprintf(stderr,
        "number of labels = %d, number of cols = %d\n",
        nvar, xg->sm_ncols);
      for (i=nvar; i<xg->sm_ncols; i++)
        (void) sprintf(xg->rowlab[xg->nrows-xg->sm_ncols+i], "Var %d", i+1);
    }
    else
    {
      /*
       * If this is being run within S, remove the labels file that was
       * created by the S function.
      */
      if (xg->data_mode == Sprocess)
        if (unlink(lab_file) != 0)
          fprintf(stderr, "read_collabels: error in unlink");
    }
  }
  else
  {
    if (init)
      for (i=0; i<xg->sm_ncols; i++)
        (void) sprintf(xg->rowlab[xg->nrows-xg->sm_ncols+i], "Var %d", i+1);
  }

  return(found);
}

Boolean
read_rowlabels(char *data_in, Boolean init, xgobidata *xg)
{
  int i, j, k;
  char lab_file[128];
  static char *suffix[] = {
    ".row", ".rowlab", ".case"
  };
  char initstr[INITSTRSIZE];
  int ncase;
  Boolean found = False;
  FILE *fp;

  if (init) {
    xg->rowlab = (char **) XtMalloc((Cardinal) xg->nrows * sizeof (char *));
    for (i=0; i<xg->nrows; i++)
      xg->rowlab[i] = (char *) XtMalloc((Cardinal) ROWLABLEN * sizeof(char));
  }

  if (!xg->is_scatmat)
  {
  /* Check if case label file exists, and open if so. */
  
    if (data_in != NULL && data_in != "" && strcmp(data_in,"stdin") != 0)
    {
      i = 0;
      while (found == False && i<3 )
      {
        (void) strcpy(lab_file, data_in);
        (void) strcat(lab_file, suffix[i++]);
        if ( (fp = fopen(lab_file,"r")) != NULL)
          found = True;
      }
    }

  /*
   * Read in case labels or initiate them to generic if no label
   * file exists
  */
    if (found)
    {
      int k, len;
      ncase = 0;

      k = 0;  /* k is the file row */
      while (fgets(initstr, INITSTRSIZE-1, fp) != NULL)
      {
        if (xg->file_read_type == read_all ||
            k == xg->file_rows_sampled[ncase])
        {
          len = MIN(INT(strlen(initstr)), ROWLABLEN-1) ;

          /* trim trailing blanks, and eliminate the carriage return */
          while (initstr[len-1] == ' ' || initstr[len-1] == '\n')
            len-- ;
          strncpy(xg->rowlab[ncase], initstr, len);
          xg->rowlab[ncase][len] = '\0' ;
  
          ncase++;
          if (ncase >= xg->nrows)
            break;
        }
        k++;  /* read the next row ... */
      }
  
      /*
       * If there aren't enough labels, use blank labels for
       * the remainder.
      */
      if (init && ncase != xg->nrows) {
        (void) fprintf(stderr, "number of labels = %d, number of rows = %d\n",
          ncase, xg->nrows);
        for (i=ncase; i<xg->nrows; i++)
          (void) sprintf(xg->rowlab[i], " ");
      }
      /*
       * If this is being run within S, remove the labels file that was
       * created by the S function.
      */
      else
      {
        if (xg->data_mode == Sprocess)
          if (unlink(lab_file) != 0)
            fprintf(stderr, "trouble in read_rowlabels");
      }
    }
    else
    {
      if (init) {  /* apply defaults if initialization; else, do nothing */

        for (i=0; i<xg->nrows; i++) {
          if (xg->file_read_type == read_all)
            (void) sprintf(xg->rowlab[i], "%d", i+1);
          else
            (void) sprintf(xg->rowlab[i], "%d", xg->file_rows_sampled[i]+1);
        }
      }
    }
  }
  else  /* scatterplot matrix */
  {
    
    if (data_in != NULL && data_in != "" && strcmp(data_in,"stdin") != 0)
    {
      i = 0;
      while (found == False && i<3 )
      {
        (void) strcpy(lab_file, data_in);
        (void) strcat(lab_file, suffix[i++]);
        if ( (fp = fopen(lab_file,"r")) != NULL)
          found = True;
      }
    }
    /*
     * Read in case labels or initiate them to generic if no label
     * file exists
    */
    if (found)
    {
      int len;
      ncase = 0;
  
      while (fgets(initstr, INITSTRSIZE-1, fp) != NULL)
      {
        len = MIN(INT(strlen(initstr)), ROWLABLEN-1) ;
        while (initstr[len-1] == ' ')
          len-- ;
        strncpy(xg->rowlab[ncase], initstr, len);
        xg->rowlab[ncase][len] = '\0' ;
  
        ncase++;
        if (ncase >= xg->sm_nrows)
          break;
      }
  
      /*
       * If there aren't enough labels, use blank labels for
       * the remainder.
      */
      if (init && ncase != xg->sm_nrows) {
        (void) fprintf(stderr, "number of labels = %d, number of rows = %d\n",
          ncase, xg->sm_nrows);
        for (i=ncase; i<xg->sm_nrows; i++)
          (void) sprintf(xg->rowlab[i], " ");
      }
      /*
       * If this is being run within S, remove the labels file that was
       * created by the S function.
      */
      else {
        if (xg->data_mode == Sprocess)
          if (unlink(lab_file) != 0)
            fprintf(stderr, "trouble in read_rowlabels");
      }
    } else {

      if (init)  /* apply defaults if initialization; else, do nothing */
        for (i=0; i<xg->sm_nrows; i++)
          (void) sprintf(xg->rowlab[i], "%d", i+1);
    }

    /* Fill out the rest of the labels */
    k=xg->sm_nrows;
    for (j=0; j<(xg->sm_ncols*(xg->sm_ncols-1)/2-1); j++)
      for (i=0; i<xg->sm_nrows; i++)
      {
        (void) sprintf(xg->rowlab[k], "%s", xg->rowlab[i]);
        k++;
      }

    (void) read_collab_to_row(data_in, init, xg);
  
  }
  return(found);
}

Boolean
read_collabels(char *data_in, Boolean init, xgobidata *xg)
{
  char lab_file[128];
  static char *suffix[] = {
    ".col", ".column", ".collab", ".var"
  };
  int i, j, nvar = 0;
  Boolean found = False;
  FILE *fp;
  char initstr[INITSTRSIZE];
  char *lbl;
  char *lbl_short;

  /*
   * Use ncols here; allocate space for the brushing
   * groups variable label.
  */
  xg->collab = (char **) XtMalloc(
    (Cardinal) xg->ncols * sizeof (char *));
  xg->collab_short = (char **) XtMalloc(
    (Cardinal) xg->ncols * sizeof (char *));
  xg->collab_tform1 = (char **) XtMalloc(
    (Cardinal) xg->ncols * sizeof (char *));
  xg->collab_tform2 = (char **) XtMalloc(
    (Cardinal) xg->ncols * sizeof (char *));

  for (j=0; j<xg->ncols; j++) {
    xg->collab[j] = (char *) XtMalloc(
      (Cardinal) COLLABLEN * sizeof(char));
    xg->collab_short[j] = (char *) XtMalloc(
      (Cardinal) COLLABLEN * sizeof(char));
    xg->collab_tform1[j] = (char *) XtMalloc(
      (Cardinal) (COLLABLEN+16) * sizeof(char));
    xg->collab_tform2[j] = (char *) XtMalloc(
      (Cardinal) (COLLABLEN+2*16) * sizeof(char));
  }


  if (!xg->is_scatmat)
  {

    /*
     * Check if variable label file exists, and open if so.
    */
    if (data_in != NULL && data_in != "" && strcmp(data_in,"stdin") != 0) {
      i = 0;
      while (found == False && i<4 ) {
        (void) strcpy(lab_file, data_in);
        (void) strcat(lab_file, suffix[i++]);
        if ( (fp = fopen(lab_file,"r")) != NULL)
          found = True;
      }
    }
  
    /*
     * Read in variable labels or initiate them to generic if no label
     * file exists
    */
    if (found)
    {
      int len;
      nvar = 0;

      lbl = XtMalloc(INITSTRSIZE * sizeof(char));
      lbl_short = XtMalloc(INITSTRSIZE * sizeof(char));
  
      while (fgets(initstr, INITSTRSIZE-1, fp) != NULL) {
        lbl = strtok((char *)initstr, "|");
        lbl_short = strtok((char *) NULL, "|");

        len = MIN(INT(strlen(lbl)), COLLABLEN-1) ;
        while (lbl[len-1] == '\n' || lbl[len-1] == ' ')
          len-- ;
        strncpy(xg->collab[nvar], lbl, len) ;
        xg->collab[nvar][len] = '\0' ;

        if (lbl_short == NULL || strlen(lbl_short) == 0)
          strcpy(xg->collab_short[nvar], xg->collab[nvar]) ;
        else {
          len = MIN(INT(strlen(lbl_short)), COLLABLEN-1) ;
          while (lbl_short[len-1] == '\n' || lbl_short[len-1] == ' ')
            len-- ;
          strncpy(xg->collab_short[nvar], lbl_short, len) ;
          xg->collab_short[nvar][len] = '\0' ;
        }
  
        nvar++;
        if (nvar >= xg->ncols_used)
          break;
      }

      if (init && nvar != xg->ncols_used)
      {
        (void) fprintf(stderr,
          "number of labels = %d, number of cols = %d\n",
          nvar, xg->ncols_used);

        if (xg->single_column) {
          (void) sprintf(xg->collab[1], "%s", xg->collab[0]);
          (void) sprintf(xg->collab[0], "Index");

        } else {
          for (i=nvar; i<xg->ncols_used; i++) {
            (void) sprintf(xg->collab[i], "Var %d", i+1);
            (void) sprintf(xg->collab_short[i], "V%d", i+1);
          }
        }
      }
      else
      {
      /*
       * If this is being run within S, remove the labels file that was
       * created by the S function.
      */
        if (xg->data_mode == Sprocess)
          if (unlink(lab_file) != 0)
            fprintf(stderr, "read_collabels: error in unlink");
      }

/*
 * You would think I needed these lines, but they cause a
 * core dump.
 *
      XtFree((XtPointer) lbl);
      XtFree((XtPointer) lbl_short);
*/
    }
    else
    {
      if (init) {
        for (i=0; i<xg->ncols_used; i++) {
          (void) sprintf(xg->collab[i], "Var %d", i+1);
          (void) sprintf(xg->collab_short[i], "V%d", i+1);
        }
      }
    }
  
    /*
     * Put "group" in the extra (unmapped) column label.
    */
    if (init) {
      strcpy(xg->collab[xg->ncols-1], "group");
      for (j=0; j<xg->ncols; j++) {
        (void) strcpy(xg->collab_tform1[j], xg->collab[j]);
        (void) strcpy(xg->collab_tform2[j], xg->collab[j]);
      }
    }
  }
  else
  {
    for (j=0; j<xg->ncols_used; j++)
      (void) sprintf(xg->collab[j], "Var %d", j+1);
    strcpy(xg->collab[xg->ncols-1], "group");
    for (j=0; j<xg->ncols; j++) {
      (void) strcpy(xg->collab_tform1[j], xg->collab[j]);
      (void) strcpy(xg->collab_tform2[j], xg->collab[j]);
    }
  }

  return(found);
}

Boolean reread_dat(char *rootname, xgobidata *xg)
{
  int i, j;
  char *suffix = ".dat";
  char fname[128];
  Boolean found = False;
  FILE *fp;
  Boolean caught_error = False;
  char word[64];
  int fs;
  Boolean show_missings_warning = False;

  sprintf(fname, rootname);

  if (fname != NULL) {
    if ( (fp = fopen(fname,"r")) != NULL)
      found = True;
  }
  if (!found) {
    strcat(fname, suffix);
    if ( (fp = fopen(fname,"r")) != NULL)
      found = True;
  }
  if (found) {
    i = 0;
    xg->nmissing = 0;
    while (i < xg->nrows) {

      /* skip over comment lines */
      if (find_data_start(fp) == False) {
        caught_error = True;
        break;
      }

      if (caught_error) break;
      for (j=0; j<xg->ncols-1; j++) {
        /*
         * This doesn't handle the special case in which the input is
         * a single column.
        */

        fs = fscanf(fp, "%s", word);
        /*
         * Check for missings, update is_missing if it exists.
         * Don't try to support adding missings at this late date ...
        */
        if (fs == EOF) {
          show_message("Insufficient data\n", xg);
          caught_error = True;
          break;
        } else if (fs < 0) {
          show_message("Problem with input data", xg);
          caught_error = True;
          break;
        } else {
          if ( strcasecmp(word, "na") == 0 || strcmp(word, ".") == 0 ) {
            if (xg->missing_values_present) {
              xg->nmissing++;
              xg->is_missing[i][j] = 1;
              xg->raw_data[i][j] = 0.0;
            } else {
              show_missings_warning = True;
              xg->raw_data[i][j] = 0.0;
            }
          } else {
            xg->raw_data[i][j] = (float) atof(word);
          }
        }
      }
      i++ ;
    }

    /*
     * If the previous xgobi had missings but this one doesn't,
     * there could be some serious resetting to do.  For the
     * moment, just reset one flag.
    */
    if (xg->missing_values_present && xg->nmissing == 0)
      xg->missing_values_present = False;

    /*
     * If the previous xgobi had no missing values, but this data 
     * contains missings flags, show this warning message.
    */
    if (show_missings_warning) {
      char message[MSGLENGTH];
      sprintf(message,
        "We cannot initiate missing values structures here;\n");
      strcat(message,
        "simply setting missing values to zero.");
      show_message(message, xg);
    }

    fclose(fp);  /* close or fclose? */

    copy_raw_to_tform(xg);
    /*
     * Reset transformations to the default
    */
    reset_tform(xg);

    update_lims(xg);
    update_world(xg);
    world_to_plane(xg);
    plane_to_screen(xg);
  
    init_tickdelta(xg);
    if (xg->is_xyplotting)
      init_ticks(&xg->xy_vars, xg);
    else if (xg->is_plotting1d) 
      init_ticks(&xg->plot1d_vars, xg); 
  }
  return(found);
}

Boolean
read_vgroups(char *data_in, Boolean init, xgobidata *xg)
/*
 * Read in the grouping numbers for joint scaling of variables
*/
{
  char lab_file[115], suffix[15];
  int itmp, i, j;
  Boolean found = False;
  FILE *fp;

  if (init) {
    xg->vgroup_ids_ori = (int *) XtMalloc((Cardinal) xg->ncols * sizeof(int));
    xg->vgroup_ids = (int *) XtMalloc((Cardinal) xg->ncols * sizeof(int));
  }

  (void) strcpy(suffix, ".vgroups");

  if (data_in != NULL && data_in != "" && strcmp(data_in, "stdin") != 0) {
    i = 0;
    (void) strcpy(lab_file, data_in);
    (void) strcat(lab_file, suffix);
    if ( (fp = fopen(lab_file,"r")) != NULL)
      found = True;
  }

  if (found) {
    i = 0;
    while ((fscanf(fp, "%d", &itmp) != EOF) && (i < xg->ncols_used))
      xg->vgroup_ids_ori[i++] = itmp;

    /*
     * Add the vgroup_id value for the extra column.
    */
    if (i == xg->ncols_used)
      xg->vgroup_ids_ori[i] = i;

    if (init && i < xg->ncols_used) {
      (void) fprintf(stderr,
        "Number of variables and number of group types do not match.\n");
      (void) fprintf(stderr,
        "Creating extra generic groups.\n");
      for (j=i; j<xg->ncols; j++)
        xg->vgroup_ids_ori[j] = j;
    }

/*
 * Find maximum vgroup id.
    maxid = 0;
    for (i=1; i<xg->ncols; i++)
      if (xg->vgroup_ids_ori[i] > maxid)
        maxid = xg->vgroup_ids_ori[i];
*/

/*
     * Find minimum vgroup id, set it to 0.  Find next, set it to 1; etc.
    id = 0;
    newid = -1;
    while (id <= maxid) {
      found = 0;
      for (j=0; j<xg->ncols; j++) {
        if (xg->vgroup_ids_ori[j] == id) {
          newid++;
          found = 1;
          break;
        }
      }
      if (found)
        for (j=0; j<xg->ncols; j++)
          if (xg->vgroup_ids_ori[j] == id)
            xg->vgroup_ids_ori[j] = newid;
      id++;
    }
*/
    resort_vgroup_ids(xg, xg->vgroup_ids_ori);
  }

  else {
    if (init)
      for (i=0; i<xg->ncols; i++)
        xg->vgroup_ids_ori[i] = i;
  }

  for (i=0; i<xg->ncols; i++)
    xg->vgroup_ids[i] = xg->vgroup_ids_ori[i];

  return(found);
}

void
free_rgroups(xgobidata *xg) {
  int i, j;
  for (i=0; i<xg->nrgroups; i++)
    for (j=0; j<xg->rgroups[i].nels; j++)
      XtFree((XtPointer) xg->rgroups[i].els);

  XtFree((XtPointer) xg->rgroups);
  XtFree((XtPointer) xg->rgroup_ids);
}

Boolean
read_rgroups(char *data_in, Boolean init, xgobidata *xg)
/*
 * Read in the grouping numbers for joint scaling of variables
*/
{
  char lab_file[115], suffix[15];
  int itmp, i, j, k;
  Boolean found = False;
  Boolean found_rg;
  FILE *fp;
  long *nels;
  int nr;

  if (!xg->is_scatmat)
  {
    (void) strcpy(suffix, ".rgroups");
  
    if (data_in != NULL && data_in != "" && strcmp(data_in, "stdin") != 0)
    {
      i = 0;
      sprintf(lab_file, "%s%s", data_in, suffix);
      if ( (fp = fopen(lab_file,"r")) != NULL)
        found = True;
    }
  
    if (!found) {
      xg->nrgroups = 0;
    } else {
  
      /*
       * If this isn't the first time we've read files, then
       * see if the rgroups structures should be freed.
      */
      if (!init)
        if (xg->nrgroups > 0)
          free_rgroups(xg);
  
      /* rgroup_ids starts by containing the values in the file */
      xg->rgroup_ids = (long *) XtMalloc(xg->nrows * sizeof(long));
      nels = (long *) XtMalloc(xg->nrows * sizeof(long));
     
      i = 0;
      while ((fscanf(fp, "%d", &itmp) != EOF) && (i < xg->nrows))
        xg->rgroup_ids[i++] = itmp;
  
      /* check the number of group ids read -- should be nrows */
      if (init && i < xg->nrows)
      {
        (void) fprintf(stderr,
          "Number of rows and number of row group types do not match.\n");
        (void) fprintf(stderr,
          "Creating extra generic groups.\n");
        for (k=i; k<xg->nrows; k++)
          xg->rgroup_ids[k] = k;
      }
  
      /*
       * Initialize the global variables: nrows row groups,
       * nrows/10 elements in each group
      */
      xg->rgroups = (rg_struct *) XtMalloc(xg->nrows * sizeof(rg_struct));
      for (i=0; i<xg->nrows; i++) {
        nels[i] = xg->nrows/10;
        xg->rgroups[i].els = (long *)
          XtMalloc((unsigned int) nels[i] * sizeof(long));
        xg->rgroups[i].nels = 0;
        xg->rgroups[i].excluded = False;
      }
      xg->nrgroups = 0;
  
/*
 * For now, only assign linkable points to rgroups.
*/
      /*
       * On this sweep, find out how many groups there are and how
       * many elements are in each group
      */
      nr = (xg->nlinkable < xg->nrows) ? xg->nlinkable : xg->nrows;

      for (i=0; i<nr; i++) {
        found_rg = False;
        for (k=0; k<xg->nrgroups; k++) {
  
          /* if we've found this id before ... */
          if (xg->rgroup_ids[i] == xg->rgroups[k].id) {
  
            /* Reallocate els[k] if necessary */
            if (xg->rgroups[k].nels == nels[k]) {
              nels[k] *= 2;
              xg->rgroups[k].els = (long *)
                XtRealloc((XtPointer) xg->rgroups[k].els,
                  (unsigned) (nels[k] * sizeof(long)));
            }
  
            /* Add the element, increment the element counter */
            xg->rgroups[k].els[ xg->rgroups[k].nels ] = i;
            xg->rgroups[k].nels++;
  
            /*
             * Now the value in rgroup_ids has to change so that
             * it can point to the correct member in the array of
             * rgroups structures
            */
            xg->rgroup_ids[i] = k;
  
            found_rg = True;
            break;
          }
        }
  
        /* If it's a new group id, add it */
        if (!found_rg) {
          xg->rgroups[xg->nrgroups].id = xg->rgroup_ids[i]; /* from file */
          xg->rgroups[xg->nrgroups].nels = 1;
          xg->rgroups[xg->nrgroups].els[0] = i;
          xg->rgroup_ids[i] = xg->nrgroups;  /* rgroup_ids reset to index */
          xg->nrgroups++;
        }
      }
      xg->nrgroups_in_plot = xg->nrgroups;
  
      /* Reallocate everything now that we know how many there are */
      xg->rgroups = (rg_struct *) XtRealloc((XtPointer) xg->rgroups,
        xg->nrgroups * sizeof(rg_struct));
  
      /* Now reallocate the arrays within each rgroups structure */
      for (k=0; k<xg->nrgroups; k++)
        xg->rgroups[k].els = (long *)
          XtRealloc((XtPointer) xg->rgroups[k].els,
            (unsigned) (xg->rgroups[k].nels * sizeof(long)));
  
      XtFree((XtPointer) nels);
    }

  }
  else  /* is_scatmat = True */
  {
    /* rgroup_ids starts by containing the values in the file */
    xg->rgroup_ids = (long *) XtMalloc(xg->nrows * sizeof(long));
    xg->rgroups = (rg_struct *) XtMalloc(xg->nrows * sizeof(rg_struct));
    k=0;
    xg->nrgroups_in_plot = xg->nrgroups = xg->sm_nrows; 

    for (j=0; j<xg->sm_ncols*(xg->sm_ncols-1)/2; j++)
      for (i=0; i<xg->sm_nrows; i++)
        xg->rgroup_ids[k++] = i;

    for (j=0; j<xg->nrgroups; j++)
    {
      xg->rgroups[j].nels = xg->sm_ncols*(xg->sm_ncols-1)/2;
      xg->rgroups[j].els = (long *)
        XtMalloc((unsigned int) xg->rgroups[j].nels * sizeof(long));
      for (k=0; k<xg->rgroups[j].nels; k++)
        xg->rgroups[j].els[k] = k*xg->sm_nrows+j;
    }

    /* fill in the label rows */
/*
    xg->nrgroups++;
    xg->nrgroups_in_plot++;
    for (j=0; j<xg->sm_ncols; j++) {
      xg->rgroup_ids[xg->nrows-xg->sm_ncols+j] = xg->sm_nrows;
    }
    xg->rgroups[xg->nrgroups-1].nels = xg->sm_ncols;
    xg->rgroups[xg->nrgroups-1].els = (long *)
      XtMalloc((unsigned int) xg->rgroups[xg->nrgroups-1].nels * sizeof(long));
    for (k=0; k<xg->rgroups[xg->nrgroups-1].nels; k++)
      xg->rgroups[xg->nrgroups-1].els[k] = xg->nrows-xg->sm_ncols+k;
*/

    /*
     * The value of nrgroups here is one lower than the true
     * number of groups, but that's because the dummy points to
     * anchor the variable labels constitute an extra rgroup
     * that we don't want to be noticed during linked brushing
     * to a non-scatmat xgobi.
    */

    for (j=0; j<xg->sm_ncols; j++) {
      xg->rgroup_ids[xg->nrows-xg->sm_ncols+j] = xg->sm_nrows;
    }
    xg->rgroups[xg->nrgroups].nels = xg->sm_ncols;
    xg->rgroups[xg->nrgroups].els = (long *)
      XtMalloc((unsigned int) xg->rgroups[xg->nrgroups].nels * sizeof(long));
    for (k=0; k<xg->rgroups[xg->nrgroups].nels; k++)
      xg->rgroups[xg->nrgroups].els[k] = xg->nrows-xg->sm_ncols+k;

  }

  if (xg->nlinkable != xg->nrows)
    fprintf(stderr, "xg->nlinkable=%d xg->nrows=%d\n",xg->nlinkable,xg->nrows);

  if (xg->nrgroups != 0)
    fprintf(stderr, "xg->nrgroups=%d\n",xg->nrgroups);

  return(found);
}

void
readGlyphErr(void) {
  fprintf(stderr,
    "The .glyphs file must contain either one number per line,\n");
  fprintf(stderr,
     "with the number between 1 and %d; using defaults,\n", NGLYPHS);
  fprintf(stderr,
     "or a string and a number, with the string being one of\n");
  fprintf(stderr,
     "+, x, or, ft, oc, fc, .  and the number between 1 and 5.\n");
}

Boolean
read_point_glyphs(char *data_in, Boolean addsuffix, Boolean reinit,
xgobidata *xg)
{
  Boolean ok = True;
  char lab_file[128];
  int i, j, k;
  Boolean found;
  FILE *fp;
  int gid;
  glyphv glyph;
  Boolean use_defaults = False;

  if (strcmp(data_in, "stdin") != 0) {

    if (strcmp(data_in, "") == 0)
      found = False;
    else {
      /*
       * Check if glyphs file exists.
      */
      (void) strcpy(lab_file, data_in);
      if (addsuffix)
        (void) strcat(lab_file, ".glyphs");
      if ( (fp = fopen(lab_file,"r")) != NULL)
        found = True;
      else
        found = False;
    }

    if (!found && reinit == True)
      init_glyph_ids(xg);
    else {
      if (!xg->is_scatmat) {
        enum { typeAndSize, glyphNumber } glyph_format;
        int c, retval, gsize;
        char *gtype;
        gtype = XtMalloc(16 * sizeof(char));

        /*
         * For the first row, find out if we're going to be reading
         * %s %d (typeAndSize) or %d (glyphNumber)
        */
        c = getc(fp);
        glyph_format = isdigit(c) ? glyphNumber : typeAndSize;

        ungetc(c, fp);
        i = 0; k = 0;
        while (i < xg->nrows) {  /* should there be a test on k as well? */

          if (glyph_format == glyphNumber) {
            retval = fscanf(fp, "%d", &gid);
          } else {
            /*fscanf(fp, "%s", &gtype);*/
            fscanf(fp, "%s", gtype);
            gsize = 1;
            if (strcmp(gtype, ".") != 0)
              fscanf(fp, "%d", &gsize);
          }

          if (retval <= 0) {
            /* not using show_message() here; reading before xgobi startup */
            (void) fprintf(stderr, "!!Error in reading %s; using defaults.\n",
              lab_file);
            use_defaults = True;
            break;
          }

          if (xg->file_read_type == read_all ||
             (xg->file_rows_sampled != NULL && k == xg->file_rows_sampled[i]))
          {
            /* not using show_message() here; reading before xgobi startup */

            /*
             * If the input is a single number on a line
            */
            if (glyph_format == glyphNumber) {

              if (gid < 1 || gid > NGLYPHS) {
                use_defaults = True;
                break;
              }
              (void) find_glyph_type_and_size(gid, &glyph);

            /*
             * Else if the input is a string and a number
            */
            } else {
              if (strcmp(gtype, "+") == 0)
                glyph.type = 1;
              else if (strcasecmp(gtype, "x") == 0)
                glyph.type = 2;
              else if (strcasecmp(gtype, "or") == 0)
                glyph.type = 3;
              else if (strcasecmp(gtype, "fr") == 0)
                glyph.type = 4;
              else if (strcasecmp(gtype, "oc") == 0)
                glyph.type = 5;
              else if (strcasecmp(gtype, "fc") == 0)
                glyph.type = 6;
              else if (strcasecmp(gtype, ".") == 0)
                glyph.type = 7;
              else {
                readGlyphErr();
                use_defaults = True;
                break;
              }

              glyph.size = gsize;
              if (gsize < 1 || gsize > 5) {
                use_defaults = True;
                readGlyphErr();
              }
            }

            if (use_defaults) {
              break;
            }

            xg->glyph_ids[i].type = xg->glyph_now[i].type =
              xg->glyph_prev[i].type = glyph.type;
            xg->glyph_ids[i].size = xg->glyph_now[i].size =
              xg->glyph_prev[i].size = glyph.size;

            i++;  /* increment the array index */
          }
          k++;  /* increment the file's row counter */
        }
        if (use_defaults)
          init_glyph_ids(xg);
        fclose(fp);
      }

      else  /* scatterplot matrix */
      {
        Boolean use_defaults = False;
        for (i=0; i<xg->sm_nrows; i++)
        {
          if (fscanf(fp, "%d", &gid) > 0)
          {
            /* not using show_message() here; reading before xgobi startup */
            if (gid < 1 || gid > NGLYPHS)
            {
              (void) fprintf(stderr,
                "Sorry, %d is not a legal glyph number;\n", gid);
              (void) fprintf(stderr,
                "glyph numbers must be between 1 and %d; using defaults.\n",
                NGLYPHS);
              use_defaults = True;
              break;
            }
            (void) find_glyph_type_and_size(gid, &glyph);
            xg->glyph_ids[i].type = xg->glyph_now[i].type =
              xg->glyph_prev[i].type = glyph.type;
            xg->glyph_ids[i].size = xg->glyph_now[i].size =
              xg->glyph_prev[i].size = glyph.size;
          }
          else
          {
            /* not using show_message() here; reading before xgobi startup */
            (void) fprintf(stderr, "!!Error in reading %s; using defaults.\n",
              lab_file);
            use_defaults = True;
            break;
          }
        }
        if (use_defaults)
          init_glyph_ids(xg);
        fclose(fp);

        /* Fill out the rest of the labels */
        k=xg->sm_nrows;
        for (j=0; j<(xg->sm_ncols*(xg->sm_ncols-1)/2-1); j++)
          for (i=0; i<xg->sm_nrows; i++)
          {
            xg->glyph_ids[k].type = xg->glyph_now[k].type =
              xg->glyph_prev[k].type = xg->glyph_ids[i].type;
            xg->glyph_ids[k].size = xg->glyph_now[k].size =
              xg->glyph_prev[k].size = xg->glyph_ids[i].size;
            k++;
          }
        for (j=0; j<xg->sm_ncols; j++)
        {
         gid = 31;
         (void) find_glyph_type_and_size(gid, &glyph);
          xg->glyph_ids[k+j].type = xg->glyph_now[k+j].type = 
            xg->glyph_prev[k+j].type = glyph.type;
            xg->glyph_ids[k+j].size = xg->glyph_now[k+j].size =
              xg->glyph_prev[k+j].size = glyph.size;
        }/*scatmat*/
      }
    }
  }
  else
    init_glyph_ids(xg);

  return(ok);
}

Boolean
read_point_colors(char *data_in, Boolean addsuffix, Boolean reinit,
xgobidata *xg)
{
  Boolean ok = True;
  char lab_file[128];
  int i, j, k, ncases;
  Boolean found;
  Colormap cmap = DefaultColormap(display, DefaultScreen(display));
  XColor exact;
  char color_name[32];
  FILE *fp;

  if (!strcmp(data_in, ""))
    return(False);

  if (!mono) {
    if (strcmp(data_in, "stdin") != 0) {

      /*
       * If color, check if colors file exists.
      */

      (void) strcpy(lab_file, data_in);
      if (addsuffix)
        (void) strcat(lab_file, ".colors");
      if ( (fp = fopen(lab_file,"r")) != NULL)
        found = 1;
      else
        found = 0;

      if (!found && reinit == True)
        init_color_ids(xg);
      else {
        /*
         * Keep track of the color names read in;
         * assume that there aren't going to be more than 64 colors
        */
        int nc = 0, n;
        struct {char cname[32]; long pix;} cnp[64]; 
        Boolean color_read;

        if (!xg->is_scatmat)
        {
          int k = 0;
          int retval;

          ncases = 0;
          xg->got_new_paint = True;

          while (ncases < xg->nrows) {  /* should there be a test on k? */
            retval = fscanf(fp, "%s", color_name);
            if (retval <= 0)
              break;

            if (xg->file_read_type == read_all ||
               xg->file_rows_sampled != NULL && k == xg->file_rows_sampled[ncases])
            {

              if (strcmp(color_name, "Default") == 0) {
                xg->color_ids[ncases] = xg->color_now[ncases] =
                  xg->color_prev[ncases] = plotcolors.fg;
              }
              /*
               * If the color name is one
               * that I've read before, then I know the pixel value
               * and there's no need to call XParseColor.
              */
              else {
                color_read = False;
                for (n=0; n<nc; n++) {
                  if (strcmp(color_name, cnp[n].cname) == 0) {
                    color_read = True;
                    break;
                  }
                }
                if (color_read)
                  xg->color_ids[ncases] = xg->color_now[ncases] =
                    xg->color_prev[ncases] = cnp[n].pix;
                else {
                  if (XParseColor(display, cmap, color_name, &exact) &&
                      XAllocColor(display, cmap, &exact) )
                  {
                    xg->color_ids[ncases] = xg->color_now[ncases] =
                      xg->color_prev[ncases] = exact.pixel;
                  }
                  else
                  {
                    fprintf(stderr, "Error in reading %s; using defaults.\n",
                      lab_file);
                    init_color_ids(xg);
                    break;
                  }
                  /*
                   * if we had not encountered that color before,
                   * add the color name and exact.pixel to the list.
                  */
                  strcpy(cnp[nc].cname, color_name);
                  cnp[nc].pix = (long) exact.pixel;
                  nc++;
  
                }
              }
              ncases++;
            }
            k++;  /* increment the file's row counter */
          }
  
          /*
           * If there aren't enough colors supplied, let the rest
           * of the points be assigned the default color.
          */
          if (ncases != xg->nrows) {
            fprintf(stderr,
              "Your .colors file does not contain enough colors; assigning\n");
            fprintf(stderr,
              "the default color to the remainder.\n");
            for (i=ncases; i<xg->nrows; i++)
              xg->color_ids[i] = xg->color_now[i] =
                  xg->color_prev[i] = plotcolors.fg;
          }
  
          /*
           * The colors that were supplied in the file (xg->color_now)
           * should match the brushing colors (color_nums)
           * If they don't, instruct the user to set the color
           * resources.
          */
          for (i=0; i<xg->nrows; i++) {
            found = False;
  
            for (j=0; j<ncolors; j++) {
              if (xg->color_now[i] == color_nums[j]) {
                found = True;
                break;
              }
            }
            if (!found)
              break;
          }
          if (!found)
          {
            fprintf(stderr,
              "Warning:  Your .colors file contains colors that are not\n");
            fprintf(stderr,
              "available as brushing colors, which will cause weird\n");
            fprintf(stderr,
              "behavior during brushing.  To read how to set the brushing\n");
            fprintf(stderr,
              "colors, click right on the Color Menu button.\n");
          }
  
          fclose(fp);
        }
        else
        {
          ncases = 0;
          xg->got_new_paint = True;
          while (fscanf(fp, "%s", color_name) != EOF) {
            if (strcmp(color_name, "Default") == 0) {
              xg->color_ids[ncases] = xg->color_now[ncases] =
                xg->color_prev[ncases] = plotcolors.fg;
            }
            /*
             * If the color name is one
             * that I've read before, then I know the pixel value
             * and there's no need to call XParseColor.
            */
            else {
              color_read = False;
              for (n=0; n<nc; n++) {
                if (strcmp(color_name, cnp[n].cname) == 0) {
                  color_read = True;
                  break;
                }
              }
              if (color_read)
                xg->color_ids[ncases] = xg->color_now[ncases] =
                  xg->color_prev[ncases] = cnp[n].pix;
              else {
                if (XParseColor(display, cmap, color_name, &exact) &&
                    XAllocColor(display, cmap, &exact) )
                {
                  xg->color_ids[ncases] = xg->color_now[ncases] =
                    xg->color_prev[ncases] = exact.pixel;
                }
                else
                {
                  fprintf(stderr, "Error in reading %s; using defaults.\n",
                    lab_file);
                  init_color_ids(xg);
                  break;
                }
                /*
                 * if we had not encountered that color before,
                 * add the color name and exact.pixel to the list.
                */
                strcpy(cnp[nc].cname, color_name);
                cnp[nc].pix = (long) exact.pixel;
                nc++;
  
              }
            }
            ncases++;
          }
  
          /*
           * If there aren't enough colors supplied, let the rest
           * of the points be assigned the default color.
          */
          if (ncases != xg->sm_nrows) {
            fprintf(stderr,
              "Your .colors file does not contain enough colors; assigning\n");
            fprintf(stderr,
              "the default color to the remainder.\n");
            for (i=ncases; i<xg->sm_nrows; i++)
              xg->color_ids[i] = xg->color_now[i] =
                  xg->color_prev[i] = plotcolors.fg;
          }
  
          /*
           * The colors that were supplied in the file (xg->color_now)
           * should match the brushing colors (color_nums)
           * If they don't, instruct the user to set the color
           * resources.
          */
          for (i=0; i<xg->sm_nrows; i++) {
            found = False;
  
            for (j=0; j<ncolors; j++) {
              if (xg->color_now[i] == color_nums[j]) {
                found = True;
                break;
              }
            }
            if (!found)
              break;
          }
          if (!found)
          {
            fprintf(stderr,
              "Warning:  Your .colors file contains colors that are not\n");
            fprintf(stderr,
              "available as brushing colors, which will cause weird\n");
            fprintf(stderr,
              "behavior during brushing.  To read how to set the brushing\n");
            fprintf(stderr,
              "colors, click right on the Color Menu button.\n");
          }
  
          fclose(fp);

          /* Fill out the rest of the labels */
          k=xg->sm_nrows;
          for (j=0; j<(xg->sm_ncols*(xg->sm_ncols-1)/2-1); j++)
            for (i=0; i<xg->sm_nrows; i++)
            {
              xg->color_ids[k] = xg->color_now[k] =
                xg->color_prev[k] = xg->color_ids[i];
              k++;
            }
          for (j=0; j<xg->sm_ncols; j++)
          {
            xg->color_ids[k+j] = xg->color_now[k+j] = 
            xg->color_prev[k+j] = plotcolors.fg;
          } /*scatmat*/

        }
      }

    }
    else
      init_color_ids(xg);

  }
  return(ok);
}

Boolean
read_erase(char *data_in, Boolean reinit, xgobidata *xg)
/*
 * Read in the erase vector
*/
{
  char lab_file[115];
  int itmp, i, j, k, found = False;
  FILE *fp;

  if (data_in != NULL && strcmp(data_in, "stdin") != 0) {

    if (strcmp(data_in, "") == 0)
      found = False;
    else {
      /*
       * Check whether erase file exists.
      */

      (void) strcpy(lab_file, data_in);
      (void) strcat(lab_file, ".erase");
      if ( (fp = fopen(lab_file,"r")) != NULL)
        found = True;
    }
  }

  xg->erased = (unsigned short *) XtRealloc((char *) xg->erased,
    (Cardinal) xg->nrows * sizeof(unsigned short));

  if (found)
  {
    if (!xg->is_scatmat)
    {
      int k = 0;  /* k is the file row, used if file_read_type != read_all */
      i = 0;
      while ((fscanf(fp, "%d", &itmp) != EOF) && (i < xg->nrows)) {
        if (xg->file_read_type == read_all || k == xg->file_rows_sampled[i])
        {
          xg->erased[i++] = (unsigned short) itmp;
        }
        k++;
      }
  
      if (i < xg->nrows) {
        (void) fprintf(stderr, "Problem in reading file %s; \n", lab_file);
        (void) fprintf(stderr, "not enough rows\n");
      }
    }
    else  /* scatterplot matrix */
    {
      i = 0;
      while ((fscanf(fp, "%d", &itmp) != EOF) && (i < xg->sm_nrows))
        xg->erased[i++] = (unsigned short) itmp;
    
      if (i < xg->sm_nrows)
      {
        (void) fprintf(stderr, "Problem in reading file %s; \n", lab_file);
        (void) fprintf(stderr, "not enough rows\n");
      }

      k=xg->sm_nrows;
      for (j=0; j<-(xg->sm_ncols*(xg->sm_ncols-1)/2); j++)
        for (i=0; i<-xg->sm_nrows; i++)
      {
          xg->erased[k] = xg->erased[i];
          k++;
      }
    }
  }
  else
  {
    if (reinit)
      for (i=0; i<xg->nrows; i++)
        xg->erased[i] = 0;
  }

  return(found);
}

/*ARGSUSED*/
Boolean
read_connecting_lines(char *rootname, Boolean startup, xgobidata *xg)
  /* startup - Initializing xgobi? */
{
  int fs, nblocks, bsize = 500;
  Boolean ok = True;
  int jlinks = 0;
  FILE *fp;
  char *fname;

  if ((rootname == NULL) || (strcmp(rootname, "") == 0) || 
      strcmp(rootname, "stdin") == 0) {
    create_default_lines(xg);
    return(ok);
  } else {
    fname = XtMalloc(128 * sizeof(char));
    /* This is for the in-process case */
    if (rootname == (char *) NULL)
      strcpy(fname, xg->datafilename);
    /* This is for the startup case */
    else
      strcpy(fname, rootname);
    strcat(fname, ".lines");
  }

  if ((fp = fopen(fname, "r")) != NULL)
  {
    int a, b;

    xg->nlines = 0;
    /*
     * Allocate space for <bsize> connecting lines.
    */
    xg->connecting_lines = (connect_lines *) XtMalloc(
      (Cardinal) bsize * sizeof(connect_lines));
    nblocks = 1;
    while (1)
    {
      fs = fscanf(fp, "%d %d", &a, &b);
      if (fs == EOF)
        break;
      else if (fs < 0) {
        ok = False;
/* won't change this one; called at startup */
        fprintf(stderr, "Error in reading .lines file\n");
        exit(1);
      }

      if (a < 1 || b > xg->nrows) {
        ok = False;
        fprintf(stderr, "Entry in .lines file > number of rows or < 1\n");
        exit(1);
      }
      else {
        /*
         * Sort lines data such that a <= b
        */
        if (a <= b) {
          xg->connecting_lines[xg->nlines].a = a;
          xg->connecting_lines[xg->nlines].b = b;
        } else {
          xg->connecting_lines[xg->nlines].a = b;
          xg->connecting_lines[xg->nlines].b = a;
        }

        (xg->nlines)++;
        jlinks++;
        if (jlinks == bsize) {
        /*
         * Allocate space for <bsize> more connecting links.
        */
          nblocks++;

          xg->connecting_lines = (connect_lines *)
            XtRealloc((XtPointer) xg->connecting_lines,
            (unsigned) (nblocks*bsize) *
            sizeof(connect_lines));
          jlinks = 0;
        }
      }
    } /* end while */
    /*
     * Close the data file
    */
    if (fclose(fp) == EOF)
      fprintf(stderr, "Error in closing .lines file");
  }
  else /* Create defaults */
    create_default_lines(xg);

  if (fname != (char *) NULL)
    XtFree((char *) fname);
  return(ok);
}

Boolean
read_line_colors(char *rootname, Boolean addsuffix, Boolean startup,
xgobidata *xg)
  /* startup --  Initializing xgobi? */
{
  int i;
  Boolean ok = True;
  char *fname;
  char color_name[32];
  FILE *fp;
  Colormap cmap = DefaultColormap(display, DefaultScreen(display));
  XColor exact;
  char message[512];

  fname = XtMalloc(128 * sizeof(char));
  /* This is for the in-process case */
  if (rootname == (char *) NULL || rootname == "")
    strcpy(fname, xg->datafilename);
  /* This is for the startup case */
  else if (strcmp(rootname, "stdin") != 0)
    strcpy(fname, rootname);
  if (addsuffix)
    strcat(fname, ".linecolors");

  if (!mono) {
    /*
     * Check if line colors file exists.
    */
    if ( (fp = fopen(fname, "r")) == NULL)
    {
      if (!startup) {
        sprintf(message,
          "The file '%s' can't be opened for reading\n", fname);
        show_message(message, xg);
      }
      ok = False;
    }
    else   /*  if (fp != NULL) */
    {
      /*
       * Keep track of the color names read in;
       * assume that there aren't going to be more than 64 colors
      */
      int nc = 0, n;
      struct {char cname[32]; long pix;} cnp[64]; 
      Boolean color_read;

      xg->got_new_paint = True;
      for (i=0; i<xg->nlines; i++) {
        if (fscanf(fp, "%s", color_name) > 0) {
          if (strcmp(color_name, "Default") == 0) {
            xg->line_color_ids[i] = xg->line_color_now[i] =
              xg->line_color_prev[i] = plotcolors.fg;
          }
          else {
            color_read = False;
            for (n=0; n<nc; n++) {
              if (strcmp(color_name, cnp[n].cname) == 0) {
                color_read = True;
                break;
              }
            }
            if (color_read) {
              xg->line_color_ids[i] = xg->line_color_now[i] =
              xg->line_color_prev[i] = cnp[n].pix;
            }
            else {
              if (XParseColor(display, cmap, color_name, &exact) &&
                  XAllocColor(display, cmap, &exact) )
              {
                xg->line_color_ids[i] = xg->line_color_now[i] =
                xg->line_color_prev[i] = exact.pixel;
              }
              else {
                fprintf(stderr,
                  "Error in reading %s; using defaults.\n", fname);
                ok = False;
                break;
              }

              /*
               * if we had not encountered that color before,
               * add the color name and exact.pixel to the list.
              */
              strcpy(cnp[nc].cname, color_name);
              cnp[nc].pix = (long) exact.pixel;
              nc++;
            }
          }
        }
      }
      fclose(fp);
    }
  }
  XtFree(fname);
  return(ok);
}

/* ARGSUSED */
XtCallbackProc
read_selected_xgobi_files(Widget w, xgobidata *xg, XtPointer callback_data)
{
  Boolean set;
  char *rootname = (char *) NULL;
  char message[512];
  Boolean read_files = False;

/* Step 1: get the rootname */
  rootname = (char *) XtMalloc((Cardinal) 132 * sizeof(char));
  XtVaGetValues(newfname_txt, XtNstring, (String) &rootname, NULL);
  /* Having trouble with blanks ... */
  strip_blanks(rootname);

  if (rootname == (char *) NULL) {
    sprintf(message, "Is the file name blank?");
    show_message(message, xg);
    return((XtCallbackProc) 0);
  }

  XtVaGetValues(PTCOLOR_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_point_colors(rootname, True, True, xg) == 0)
      return((XtCallbackProc) 0);
    else
      xg->got_new_paint = True;

  XtVaGetValues(PTGLYPH_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_point_glyphs(rootname, True, True, xg) == 0)
      return((XtCallbackProc) 0);
    else
      xg->got_new_paint = True;

  XtVaGetValues(PTERASE_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_erase(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);
    else
      xg->got_new_paint = True;

  XtVaGetValues(ROWLABEL_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_rowlabels(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);
    else {  /* Reinitialize case list */
      extern Widget caselist_popup;
      if (caselist_popup != NULL)
        XtDestroyWidget(caselist_popup);
      caselist_popup = NULL;
      free_caselist(xg);
      build_caselist(xg);
    }

  XtVaGetValues(COLLABEL_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set) {
    if (read_collabels(rootname, False, xg) == 0) {
      return((XtCallbackProc) 0);
    }
    else {
      int j;
      extern Widget varlist_popup;
      for (j=0; j<xg->ncols_used; j++) {
        (void) strcpy(xg->collab_tform1[j], xg->collab[j]);
        (void) strcpy(xg->collab_tform2[j], xg->collab[j]);
        XtVaSetValues(xg->varlabw[j],
          XtNlabel, xg->collab_tform2[j],
          NULL);
      }
      /* Reinitialize variable list */
      if (varlist_popup != NULL)
        XtDestroyWidget(varlist_popup);
      varlist_popup = NULL;
      free_varlist(xg);
      build_varlist(xg);
    }
  }

  XtVaGetValues(LINE_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_connecting_lines(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);
    else {
      realloc_lines(xg);
      init_line_colors(xg);
    }

  XtVaGetValues(LINECOLOR_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_line_colors(rootname, True, False, xg) == 0)
      return((XtCallbackProc) 0);

  XtVaGetValues(VGROUPS_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set) {
    if (read_vgroups(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);
    else {
      /* reset the axes; this may require some fine tuning */
      update_lims(xg);
      update_world(xg);
      world_to_plane(xg);
      plane_to_screen(xg);

      init_tickdelta(xg);
      if (xg->is_xyplotting)
        init_ticks(&xg->xy_vars, xg);
      else if (xg->is_plotting1d)
        init_ticks(&xg->plot1d_vars, xg);
    }
  }

  XtVaGetValues(RGROUPS_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set) {
    if (read_rgroups(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);
    else {
      set_lgroups(True, xg);
    }
  }

  XtVaGetValues(JITTER_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_jitter_values(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);

  XtVaGetValues(NLINK_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (read_nlinkable(rootname, False, xg) == 0)
      return((XtCallbackProc) 0);

  XtVaGetValues(NDATA_TGL, XtNstate, &set, NULL);
  read_files = read_files || set;
  if (set)
    if (reread_dat(rootname, xg) == 0)
      return((XtCallbackProc) 0);

/* Testing for Juergen */
  if (!read_files)
    read_new_data(newfname_txt, xg);
/* end of test */

  /* A replot is likely to be needed */
  plot_once(xg);
}


/* ARGSUSED */
static XtCallbackProc
close_epopup_cback(Widget w, xgobidata *xg, XtPointer cback_data)
{
  XtPopdown(epopup);
}

/* ARGSUSED */
XtCallbackProc
open_import_xgobi_popup_cback(Widget w, xgobidata *xg, XtPointer cback_data)
{
  static Boolean initd = False;
  Widget fname_lab;
  Dimension width, height;
  Position x, y;
  Widget oframe, panel, choice_panel, box;
  Widget close, doit;
  char str[64];

  if (epopup == (Widget) NULL) {

    XtVaGetValues(w,
      XtNwidth, &width,
      XtNheight, &height, NULL);
    XtTranslateCoords(w,
      (Position) (width/2), (Position) (height/2), &x, &y);

    epopup = XtVaCreatePopupShell("ReadData",
      /*transientShellWidgetClass, XtParent(w),*/
      topLevelShellWidgetClass, XtParent(w),
      XtNinput,            (Boolean) True,
      XtNallowShellResize, (Boolean) True,
      XtNtitle,            (String) "Read XGobi files",
      XtNiconName,         (String) "ExtendData",
      XtNx,                x,
      XtNy,                y,
      NULL);
    if (mono) set_mono(epopup);

    oframe = XtVaCreateManagedWidget("Form",
      panedWidgetClass, epopup,
      XtNorientation, (XtOrientation) XtorientVertical,
      NULL);
    panel = XtVaCreateManagedWidget("Panel",
      formWidgetClass, oframe,
      NULL);
    if (mono) set_mono(panel);

    /* Label and text widget to capture the new file name */
    box = XtVaCreateManagedWidget("Panel",
      boxWidgetClass, panel,
      XtNorientation, (XtOrientation) XtorientHorizontal,
      XtNhSpace, 1,
      XtNvSpace, 1,
      NULL);
    if (mono) set_mono(box);
    fname_lab = (Widget) XtVaCreateManagedWidget("ReadData",
      labelWidgetClass, box,
      XtNlabel, "Root file name: ",
      XtNresize, False,
      NULL);
    if (mono) set_mono(fname_lab);
    sprintf(str, "MMMMMMMMMMMMMMMMMM");
    width = XTextWidth(appdata.font, str, strlen(str)) +
      2*ASCII_TEXT_BORDER_WIDTH;
    newfname_txt = XtVaCreateManagedWidget("MissingText",
      asciiTextWidgetClass, box,
      XtNeditType, (int) XawtextEdit,
      XtNstring, (String) xg->datafilename,
      XtNwidth, width,
      XtNdisplayCaret, (Boolean) True,
      XtNresize, False,
      NULL);
    if (mono) set_mono(newfname_txt);

    choice_panel = XtVaCreateManagedWidget("Panel",
      formWidgetClass, panel,
      XtNfromVert, (Widget) box,
      NULL);
    if (mono) set_mono(choice_panel);

    PTCOLOR_TGL = CreateToggle(xg, ".colors file", True,
      (Widget) NULL, (Widget) NULL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(PTCOLOR_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    PTGLYPH_TGL = CreateToggle(xg, ".glyphs file", True,
      (Widget) NULL, PTCOLOR_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(PTGLYPH_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    PTERASE_TGL = CreateToggle(xg, ".erase file", True,
      (Widget) NULL, PTGLYPH_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(PTERASE_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    ROWLABEL_TGL = CreateToggle(xg, ".row file", True,
      (Widget) NULL, PTERASE_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(ROWLABEL_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    COLLABEL_TGL = CreateToggle(xg, ".col file", True,
      (Widget) NULL, ROWLABEL_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(COLLABEL_TGL ,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    LINE_TGL = CreateToggle(xg, ".lines file", True,
      (Widget) NULL, COLLABEL_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(LINE_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    LINECOLOR_TGL = CreateToggle(xg, ".linecolors file", True,
      (Widget) NULL, LINE_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(LINECOLOR_TGL ,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    VGROUPS_TGL = CreateToggle(xg, ".vgroups file", True,
      (Widget) NULL, LINECOLOR_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(VGROUPS_TGL ,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    RGROUPS_TGL = CreateToggle(xg, ".rgroups file", True,
      (Widget) NULL, VGROUPS_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(RGROUPS_TGL ,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    JITTER_TGL = CreateToggle(xg, ".jit file", True,
      (Widget) NULL, RGROUPS_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(JITTER_TGL ,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    NLINK_TGL = CreateToggle(xg, ".nlinkable file", True,
      (Widget) NULL, JITTER_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(NLINK_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    sprintf(str, ".dat file (%d x %d)", xg->nrows, xg->ncols-1);
    NDATA_TGL = CreateToggle(xg, str, True,
      (Widget) NULL, NLINK_TGL, (Widget) NULL, False,
      ANY_OF_MANY, choice_panel, "ReadData");
    XtAddCallback(NDATA_TGL,  XtNcallback,
      (XtCallbackProc) read_toggle_cback, (XtPointer) xg);

    XtManageChildren(read_choice, NCHOICES);

    doit = XtVaCreateManagedWidget("ReadData",
      commandWidgetClass, panel,
      XtNlabel, (String) "Read files",
      XtNfromVert, (Widget) choice_panel,
      NULL);
    if (mono) set_mono(doit);
    XtAddCallback(doit, XtNcallback,
      (XtCallbackProc) read_selected_xgobi_files, (XtPointer) xg);

    close = XtVaCreateManagedWidget("ReadData",
      commandWidgetClass, oframe,
      XtNshowGrip, (Boolean) False,
      XtNskipAdjust, (Boolean) True,
      XtNlabel, (String) "Click here to dismiss",
      NULL);
    if (mono) set_mono(close);
    XtAddCallback(close, XtNcallback,
      (XtCallbackProc) close_epopup_cback, (XtPointer) NULL);
  }

  XtPopup(epopup, (XtGrabKind) XtGrabNone);
  XRaiseWindow(display, XtWindow(epopup));

  if (!initd)
  {
    set_wm_protocols(epopup);
    initd = True;
  }

}

Boolean
read_nlinkable(char *data_in, Boolean init, xgobidata *xg)
/*
 * Read in the number of rows to be linked.
*/
{
  char lab_file[115], suffix[15];
  int itmp;
  Boolean found = False;
  FILE *fp;

  (void) strcpy(suffix,".nlinkable");

  if (!xg->is_scatmat) {

    if ((data_in != NULL) &&
        (strcmp(data_in, "") != 0) &&
        (strcmp(data_in, "stdin") != 0))
    {
      (void) strcpy(lab_file, data_in);
      (void) strcat(lab_file, suffix);
      if ( (fp = fopen(lab_file,"r")) != NULL)
        found = True;
    }

    /*
     * Initialize nlinkable to be all the rows; if not
     * initialization, leave its value alone.
    */
    if (init)
      xg->nlinkable = xg->nrows;
  
    
    if (found)
    {
      fscanf(fp, "%d", &itmp);
      if (itmp > 0 && itmp <= xg->nrows)
        xg->nlinkable = itmp;
      fclose(fp);
    }
  }
  else
  {
    xg->nlinkable = xg->nrows - xg->sm_ncols;
  }

  if (xg->nrows_in_plot == xg->nrows) xg->nlinkable_in_plot = xg->nlinkable;
  else {
    int i;
    xg->nlinkable_in_plot = 0;
    for (i=0; i<xg->nlinkable; i++)
      if (!xg->excluded[i])
        xg->nlinkable_in_plot++;
  }

  if (xg->nlinkable != xg->nrows)
    fprintf(stdout, "nlinkable = %d\n", xg->nlinkable);

  return(found);
}

Boolean
read_jitter_values(char *data_in, Boolean reinit, xgobidata *xg)
/*
 * Read in a .jit file of jittered values, nrows by ncols
*/
{
  char lab_file[115], suffix[15];
  long ltmp;
  int i, j;
  Boolean found = False;
  FILE *fp;

  /* XtCalloc initializes to zero */
  xg->jitter_data = (long **) XtRealloc((XtPointer)
    xg->jitter_data, xg->nrows * sizeof(long *));
  for (i=0; i<xg->nrows; i++)
    xg->jitter_data[i] = (long *) XtCalloc(xg->ncols, sizeof(long));

  (void) strcpy(suffix, ".jit");

  if (data_in != NULL && data_in != "" && strcmp(data_in, "stdin") != 0)
  {
    i = 0;
    (void) strcpy(lab_file, data_in);
    (void) strcat(lab_file, suffix);
    if ( (fp = fopen(lab_file,"r")) != NULL)
      found = True;
  }

  if (found)
  {
    int icount = 0;
    int nr = xg->nrows;
    int nc = xg->ncols_used;
    i = j = 0;
    while ((fscanf(fp, "%ld", &ltmp) != EOF) && (icount < nr*nc)) {
      icount++;
      xg->jitter_data[i][j] = ltmp;
      j++;
      if (j == nc) {
        i++;
        j = 0;
      }
    }

    if (i < xg->nrows)
    {
      (void) fprintf(stderr, "Problem in reading file %s; \n", lab_file);
      (void) fprintf(stderr, "not enough rows\n");
    }
  }
  else
  {
    if (reinit)
      for (i=0; i<xg->nrows; i++)
        for (j=0; j<xg->ncols_used; j++)
          xg->jitter_data[i][j] = 0;
  }

  update_world(xg);
  world_to_plane(xg);
  plane_to_screen(xg);

  return(found);
}
