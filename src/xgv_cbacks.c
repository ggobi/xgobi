#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <X11/keysym.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"
#include "xgvis.h"

extern void update_plot(xgobidata *);
extern void configure_pos_data(void);
extern void mds_once(Boolean, Boolean, FILE *, FILE *);
extern void set_vgroups(void);
extern void set_dist_matrix_from_edges(struct array *, struct array *, int);
extern void set_dist_matrix_from_pos(struct array *, struct array *, double);
extern void set_dist_matrix_from_pos_dot(struct array *, struct array *, int);
extern void reset_data(void);
extern void reinit_stress(void);
extern void scramble_data(void);
extern void set_distance_factor(void);

static void
reset_dims_label(void) {
  char str[32];
  extern Widget mds_dims_label;

  sprintf(str, "Dim (k): %d", mds_dims);
  XtVaSetValues(mds_dims_label,
    XtNstring, (String) str,
    NULL);
}

/* ARGSUSED */
XtCallbackProc
mds_dimsleft_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  if (mds_dims > 1) {
    mds_dims--;
    reset_dims_label();

    configure_pos_data();
    set_vgroups();

    update_plot(&xgobi);
    plot_once(&xgobi);

    mds_once(False, False, NULL, NULL);
  }
}

/* ARGSUSED */
XtCallbackProc
mds_dimsright_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  if (mds_dims < xgobi.ncols_used) {
    mds_dims++;
    reset_dims_label();

    configure_pos_data();
    set_vgroups();

    update_plot(&xgobi);
    plot_once(&xgobi);

    mds_once(False, False, NULL, NULL);
  }
}

/* ARGSUSED */
XtCallbackProc
PopUpDistMenu(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * Pop up the distance matrix menu.
*/
{
  Dimension width, height;
  Position x, y;
  static int initd = 0;

  if (!initd)
  {
    XtVaGetValues(w,
      XtNwidth, &width,
      XtNheight, &height, NULL);
    XtTranslateCoords(w,
      (Position) (width/2), (Position) (height/2), &x, &y);

    XtVaSetValues(dist_popup,
      XtNx, x, XtNy, y, NULL); 

    initd = 1;
  }

  XtPopup(dist_popup, XtGrabNone);
}

/* ARGSUSED */
XtCallbackProc
PopDownDistMenu(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * Close the distance matrix menu.
*/
{
    XtPopdown(dist_popup);
}

/* ARGSUSED */
XtCallbackProc
choose_dist_cback(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * Close the distance matrix menu.
*/
{
  int i;
  Arg args[1];
  Boolean selected = False;

  for (i=0; i<NDISTTYPES-1; i++)
  {
    XtSetArg(args[0], XtNstate, &selected);
    XtGetValues(dist_types[i], args, 1);
    if (selected)
    {
      dist_type = i ;
      break;
    }
  }

  if (dist_type < 0 || dist_type > NDISTTYPES-1)
  {
    fprintf(stderr, "Sorry, that's not a valid dist_type.\n");
  }
  else
  {
    /* dist_type has been set, let's call the right routine. */
    switch (dist_type) {
    case LINK:
      set_dist_matrix_from_edges(&dist, &edges, pos.nrows);
      break;
    case EUCLIDIAN:
      set_dist_matrix_from_pos(&dist, &pos, 2.0);
      break;
    case MANHATTAN:
      set_dist_matrix_from_pos(&dist, &pos, 1.0);
      break;
    case DOTPROD:
    case COSDIST:
      set_dist_matrix_from_pos_dot(&dist, &pos, dist_type);
      break;
    case USER_SUPPLIED:
    case ADJACENCY:
    case MAHALANOBIS:
      fprintf(stderr, "Not implemented yet.\n");
    }
  }
}

/* ARGSUSED */
XtCallbackProc
reset_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  reset_data();

  configure_pos_data();

  update_plot(&xgobi);
  plot_once(&xgobi);

  reinit_stress();
}

/* ARGSUSED */
XtCallbackProc
scramble_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  scramble_data();

  configure_pos_data();

  update_plot(&xgobi);
  plot_once(&xgobi);

  reinit_stress();
}

/* ARGSUSED */
XtCallbackProc
run_cback(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * This callback defines the actions associated with the run button.
*/
{
  xgv_is_running = !xgv_is_running;

  if (xgv_is_running) {

    /*
     * In case points have been moved in xgobi, we'd better
     * copy in the current values in raw_data.
    */
    int i, j;
    for (i=0; i<xgobi.nrows; i++)
      for (j=0; j<xgobi.ncols_used; j++)
        pos.data[i][j] = xgobi.raw_data[i][j] ;

    (void) XtAppAddWorkProc(app_con, RunWorkProcs, NULL);
  }

  setToggleBitmap(w, xgv_is_running);
}

/* ARGSUSED */
static XtTimerCallbackProc
turn_off_scaling(xgobidata *xg, XtIntervalId id)
{
    is_rescale = 0;
}

/* ARGSUSED */
XtCallbackProc
xgv_rescale_cback(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * This callback defines the actions associated with the rescale button.
*/
{
  XtIntervalId rescale_timeout_id;
  is_rescale = 1;

  rescale_timeout_id = XtAppAddTimeOut(app_con,
    500, (XtTimerCallbackProc) turn_off_scaling, NULL );
}

/* ARGSUSED */
XtCallbackProc
Quit(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * This callback defines the actions associated with the 'exit' button.
*/
{
    extern void do_exit(int);

    do_exit(0);
}

/* ARGSUSED */
/*
XtCallbackProc
mds_group_cback(Widget w, XtPointer client_data, XtPointer callback_data)
 * This callback turns grouping in MDS on and off.
{
  mds_group = !mds_group;
}
*/

/* ARGSUSED */
/*
XtCallbackProc
mds_casewise_cback(Widget w, XtPointer client_data, XtPointer callback_data)
 * This callback turns the use of persistent (sticky) labels
 * for MDS case diagnostics on and off.
{
  mds_casewise = !mds_casewise;
}
*/


/* ARGSUSED */
XtCallbackProc
mds_lnorm_cback (Widget w, XtPointer client_data, XtPointer slideposp)
{
  Arg args[1];
  char str[30];
  extern Widget mds_lnorm_label;

  float slidepos = * (float *) slideposp;

  mds_lnorm = (double) (5.0 * slidepos) + 1.0 ;
  sprintf(str, "%s: %3.1f ", "Minkowski norm (m)", mds_lnorm);
  XtSetArg(args[0], XtNstring, str);
  XtSetValues(mds_lnorm_label, args, 1);

  configure_pos_data();
  set_vgroups();

  update_plot(&xgobi);
  plot_once(&xgobi);
}

/* ARGSUSED */
XtCallbackProc
mds_power_cback (Widget w, XtPointer client_data, XtPointer slideposp)
/*
 * Adjust the power to which we raise the distance matrix.
*/
{
  Arg args[1];
  char str[30];
  extern Widget mds_power_label;

  float slidepos = * (float *) slideposp;

  mds_power = 6. * slidepos ;
  sprintf(str, "%s: %3.1f ", "Power (p) of D", mds_power);
  XtSetArg(args[0], XtNstring, str);
  XtSetValues(mds_power_label, args, 1);

/*
 * When power changes, reset the distance factor.
*/
  set_distance_factor();

/*
 * The third column of the diagnostics matrix has to be
 * reset as well, and this may be the easiest way to do it
 * because the indices are not simple.
*/
  mds_once(False, False, NULL, NULL);
}

/* ARGSUSED */
XtCallbackProc
mds_weightpow_cback (Widget w, XtPointer client_data, XtPointer slideposp)
/*
 * Adjust the power to which we raise the distance matrix.
*/
{
  Arg args[1];
  char str[30];
  extern Widget mds_weightpow_label;

  float slidepos = * (float *) slideposp;

  mds_weightpow = (slidepos - 0.5)*8.0 ;
  sprintf(str, "%s: %4.1f ", "Weight power", mds_weightpow);
  XtSetArg(args[0], XtNstring, str);
  XtSetValues(mds_weightpow_label, args, 1);
}


/* ARGSUSED */
XtCallbackProc
mds_stepsize_cback (Widget w, XtPointer client_data, XtPointer slideposp)
/*
 * Adjust the stepsize (range: currently, 0.000:1.000).
*/
{
  Arg args[1];
  char str[30];
  extern Widget mds_stepsize_label;

  float slidepos = * (float *) slideposp;  /* 0:1 */

  mds_stepsize = slidepos*slidepos;
  sprintf(str, "%s: %3.3f ", "Stepsize", mds_stepsize);
  XtSetArg(args[0], XtNstring, str);
  XtSetValues(mds_stepsize_label, args, 1);
}

/* ARGSUSED */
XtCallbackProc
mds_iterate_cback (Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * Step: one step through the mds loop.
*/
{
  mds_once(True, False, NULL, NULL);
  update_plot(&xgobi);
  RunWorkProc((xgobidata *) &xgobi);
  plot_once(&xgobi);
}

/* ARGSUSED */
static XtCallbackProc
fcancel_cback(Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * If the plot window is fully or partially exposed, clear and redraw.
*/
{
  XtDestroyWidget(XtParent(XtParent(w)));
}

/* ARGSUSED */
XtEventHandler
save_distance_matrix_go (Widget w, XtPointer cldata, XEvent *event)
{
  XKeyPressedEvent *evnt = (XKeyPressedEvent *) event;
  KeySym key;
  char *fname;
  FILE *fp;

  key = XLookupKeysym(evnt, 0);
  if (key == XK_Return)
  {
    /*
     * w is ftext; XtParent(w) = fform
     * XtParent(XtParent(w)) = fpopup
     * XtParent(XtParent(XtParent(w))) = the parent we're interested in
    */
    XtVaSetValues(XtParent(XtParent(XtParent(w))),
      XtNstate, (Boolean) False,
      NULL);

    XtVaGetValues(w, XtNstring, &fname, NULL);
    if ( (fp = fopen(fname, "w")) == NULL)
    {
      char message[MSGLENGTH];
      sprintf(message, "Failed to open the file '%s' for writing.\n", fname);
      show_message(message, &xgobi);
    }
    else
    {
      int i, j;
      for (i = 0; i < dist.nrows; i++) {
        for (j = 0; j < dist.ncols; j++) {
          fprintf(fp, "%2.3f ", dist.data[i][j]);
        }
        fprintf(fp, "\n");
      }
      fflush(fp);
    }

    XtDestroyWidget(XtParent(XtParent(w)));
  }
}

/* ARGSUSED */
XtCallbackProc
save_distance_matrix (Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * Write out the distance matrix.
*/
{
/*
 * Create a popup window to get the name; then call
 * the save routine.
*/
  Widget fpopup, fform, flabel, ftext, fcancel;
  Dimension width, height;
  Position x, y;
  Cursor text_cursor = XCreateFontCursor(display, XC_xterm);

  XtVaGetValues(shell,
    XtNwidth, &width,
    XtNheight, &height,
    NULL);
  XtTranslateCoords(w,
    (Position) (width/2), (Position) (height/2), &x, &y);

/*
 * Create the popup itself.
*/
  fpopup = XtVaCreatePopupShell("FSavePopup",
    /*
     * If this is a topLevelShell, the user is asked to
     * place it; if it's transient, it pops up where we
     * tell it to.
    */
    /*topLevelShellWidgetClass, w,*/
    transientShellWidgetClass, shell,
    XtNx, (Position) x,
    XtNy, (Position) y,
    XtNinput, (Boolean) True,
    XtNallowShellResize, (Boolean) True,
    XtNtitle, (String) "Solicit File Name",
    NULL);

/*
 * Create the form widget.
*/
  fform = XtVaCreateManagedWidget("FSaveForm",
    formWidgetClass, fpopup,
    NULL);

  flabel = XtVaCreateManagedWidget("FSaveText",
    labelWidgetClass, fform,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNright, (XtEdgeType) XtChainLeft,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNbottom, (XtEdgeType) XtChainTop,
    XtNlabel, (String) "Enter file name for distance matrix:",
    NULL);

/*
 * Create the text widget to solicit the filename.
*/
  ftext = XtVaCreateManagedWidget("FSaveName",
    asciiTextWidgetClass, fform,
    XtNfromVert, (Widget) flabel,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNright, (XtEdgeType) XtChainRight,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNbottom, (XtEdgeType) XtChainTop,
    XtNresizable, (Boolean) True,
    XtNeditType, (int) XawtextEdit,
    XtNresize, (XawTextResizeMode) XawtextResizeWidth,
    NULL);

  XtAddEventHandler(ftext, KeyPressMask, FALSE,
     (XtEventHandler) save_distance_matrix_go, (XtPointer) NULL);
/*
 * Add a cancel button
*/
  fcancel = XtVaCreateManagedWidget("Command",
    commandWidgetClass, fform,
    XtNfromVert, ftext,
    XtNlabel, "Cancel",
    NULL);
  XtAddCallback(fcancel, XtNcallback,
    (XtCallbackProc) fcancel_cback, (XtPointer) NULL);

  XtPopup(fpopup, XtGrabExclusive);
  XRaiseWindow(display, XtWindow(fpopup));

  XDefineCursor(display, XtWindow(ftext), text_cursor);
/*
 * Should do something more clever here -- get the size
 * of the window and place the cursor that way.
*/
  XWarpPointer(display, None, XtWindow(ftext), 0,0,0,0, 10,40);
}


/* ARGSUSED */
XtCallbackProc
mds_launch_cback (Widget w, XtPointer client_data, XtPointer callback_data)
/*
 * Launch the xgobi child containing diagnostic data.
*/
{
  extern xgobidata xgobi, xgobi_diag;
  extern Widget mds_launch_ntxt;
  int i, j;
  char **col_name;
  int nc = 5;
  static char *clab[] = {"d_config", "D^p", "Weight", "i", "j"};
  char fname[512], config_basename[512];
  FILE *fp, *fpdat, *fprow;
  static int iter = 0;
  char message[MSGLENGTH];
  char command[512];
  char xgobi_exec[512];
  char *xgobidir;
  struct stat buf;
  char *subset_str;
  int subset_size;

  XtVaGetValues(mds_launch_ntxt, XtNstring, &subset_str, NULL);
  if (strlen(subset_str) == 0)
    subset_size = ndistances;
  else
    subset_size = atoi(subset_str);
  if (subset_size == 0)
    subset_size = ndistances;
  else subset_size = MIN(subset_size, ndistances);

  col_name = (char **) XtMalloc(
    (Cardinal) nc * sizeof (char *));
  for (j=0; j<nc; j++)
    col_name[j] = clab[j];

  sprintf(config_basename, "Shepard_diagram_%d", iter);
  
  /* Write out the data */
  sprintf(fname, "%s.dat", config_basename);
  if ( (fpdat = fopen(fname, "w")) == NULL) {
    sprintf(message,
      "The file '%s' can not be created\n", fname);
    show_message(message, &xgobi);
    return(0);
  } else {
    sprintf(fname, "%s.row", config_basename);
    if ( (fprow = fopen(fname, "w")) == NULL) {
      sprintf(message,
        "The file '%s' can not be created\n", fname);
      show_message(message, &xgobi);
      fclose(fpdat);
      return(0);
    }
  }


  /*
   * This takes care of writing out the data and the row labels
  */
  mds_once(False, True, fpdat, fprow);
  fclose(fpdat);
  fclose(fprow);

  /* Write out the column labels */
  sprintf(fname, "%s.col", config_basename);
  if ( (fp = fopen(fname, "w")) == NULL) {
    sprintf(message,
      "The file '%s' can not be created\n", fname);
    show_message(message, &xgobi);
    return(0);
  } else {

    for (i=0; i<nc; i++) {
      fprintf(fp, "%s\n", col_name[i]);
    }
    fclose(fp);
  }
  
  XtFree((char *) col_name);

  xgobidir = getenv("XGOBID");
  if(xgobidir && strlen(xgobidir) > 0) { 
    sprintf(xgobi_exec, "%s/bin/xgobi", xgobidir);

    /* If no luck there, then just try 'xgobi' without a path name */
    if (stat(xgobi_exec, &buf) != 0)
      sprintf(xgobi_exec, "xgobi");
  } else
    sprintf(xgobi_exec, "xgobi");

  if (mono)
    strcat(xgobi_exec, " -mono");
  sprintf(command,
    "%s -subset %d %s &", xgobi_exec, subset_size, config_basename);
  fprintf(stderr, "%s\n", command);

  system (command);

  iter++;
}
