/***************************************************************************

    This file defines all the widgets .

*************************************************************************** */

#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"
#include <X11/keysym.h>
#include "xgvis.h"
#include "../bitmaps/leftarrow.xbm"
#include "../bitmaps/rightarrow.xbm"

extern XtCallbackProc PopUpDistMenu();
extern XtCallbackProc PopDownDistMenu();
extern XtCallbackProc choose_dist_cback();
extern XtCallbackProc PopUpThresholdPanel();

extern XtCallbackProc reset_cback();
extern XtCallbackProc scramble_cback();
extern XtCallbackProc xgvis_help_cback();
extern XtCallbackProc run_cback();
extern XtCallbackProc Quit();

extern XtCallbackProc mds_group_cback();
extern XtCallbackProc mds_lnorm_cback();
extern XtCallbackProc mds_power_cback();
extern XtCallbackProc mds_weightpow_cback();
extern XtCallbackProc mds_stepsize_cback();
extern XtCallbackProc mds_iterate_cback();
extern XtCallbackProc mds_dimsleft_cback();
extern XtCallbackProc mds_dimsright_cback();
extern XtCallbackProc save_distance_matrix();
extern XtCallbackProc mds_casewise_cback();
extern XtCallbackProc mds_launch_cback();

extern void build_stress_plotwin(Widget, Widget, Widget);
extern void build_dissim_plotwin(Widget);

Widget mds_dims_label, mds_lnorm_label, mds_power_label, mds_power_sbar;
Widget mds_weightpow_label, mds_stepsize_label;
Widget mds_launch_ntxt, mds_launch_nlbl;
Widget run_cmd[5];
#define RUN   run_cmd[0]
#define STEP  run_cmd[1]
#define RESET run_cmd[2]
#define SCRAM run_cmd[3]
#define HELP  run_cmd[4]

/* MDS with groups */
#define NGROUPBTNS 4
Widget group_menu_cmd, group_menu_btn[NGROUPBTNS] ;
Widget group_menu_box, group_menu_lab, group_menu;
static char *group_menu_btn_label[] = {
  "Within groups",
  "Between groups",
  "Anchored",
  "None"
};

/* ARGSUSED */
static XtCallbackProc
set_mds_group_cback(w, cldata, cdata)
  Widget w;
  XtPointer cldata, cdata;
{
  int k, btn;
  for (k=0; k<NGROUPBTNS; k++)
    if (group_menu_btn[k] == w) {
      btn = k;
      break;
    }
  XtVaSetValues(group_menu_cmd,
    XtNlabel, group_menu_btn_label[btn],
    NULL);

  switch (btn) {
    case 0:
      mds_group_ind = within;
      break;
    case 1:
      mds_group_ind = between;
      break;
    case 2:
      mds_group_ind = anchored;
      break;
    case 3:
      mds_group_ind = none;
      break;
    default:
      mds_group_ind = none;
  }

  fprintf(stderr, "%d\n", mds_group_ind);
}

static Widget metric_cmd[2], metricPanel;
/* ARGSUSED */
XtCallbackProc
setMetricCback(Widget w, XtPointer cd, int *which) {

  scaling_method = METRIC;
  
  setToggleBitmap(w, True);
  setToggleBitmap(metric_cmd[1], False);

  XtVaSetValues(mds_power_sbar, XtNsensitive, True, NULL);
}
/* ARGSUSED */
XtCallbackProc
setNonmetricCback(Widget w, XtPointer cd, int *which) {

  scaling_method = NONMETRIC;

  setToggleBitmap(w, True);
  setToggleBitmap(metric_cmd[0], False);

  XtVaSetValues(mds_power_sbar, XtNsensitive, False, NULL);
}

void
make_xgvis_widgets()
{
  char str[30];
  Dimension width, dwidth;
  Widget runPanel;
  Widget mdsPanel;
  Widget mds_lnorm_sbar;
  Widget mds_weightpow_sbar;
  Widget mds_stepsize_sbar;
  static Widget file_menu_cmd, file_menu, file_menu_btn[2];
  int j;
  static char *file_menu_names[] = {
    "Save distance matrix ...",
    "Exit",
  };
  static char *group_menu_str = "MDS with groups:";
  Widget mds_launch_panel;
  Widget mds_launch_cmd, mds_launch_label;

/*
  Pixmap xgvleftarr, xgvrightarr;
  Drawable root_window = RootWindowOfScreen(XtScreen(shell));
  xgvleftarr = XCreatePixmapFromBitmapData(display, root_window,
    leftarrow_bits, leftarrow_width, leftarrow_height,
    appdata.fg, appdata.bg, depth);
  xgvrightarr = XCreatePixmapFromBitmapData(display, root_window,
    rightarrow_bits, rightarrow_width, rightarrow_height,
    appdata.fg, appdata.bg, depth);
*/

  form0 = XtVaCreateManagedWidget("Form0",
    formWidgetClass, shell,
    /*XtNfont, panel_data.Font,*/
    XtNfont, appdata.font,
    NULL);

/* Define Run Panel, to contain the reset and run buttons. */
  runPanel = XtVaCreateManagedWidget("Panel",
    formWidgetClass, form0,
    XtNhorizDistance, 5,
    XtNvertDistance, 5,
    XtNorientation, (XtOrientation) XtorientVertical,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNright, (XtEdgeType) XtChainLeft,
    XtNbottom, (XtEdgeType) XtChainTop,
    NULL);

/* File menu */

  file_menu_cmd = XtVaCreateManagedWidget("Command",
    menuButtonWidgetClass, runPanel,
    XtNlabel, (String) "File",
    XtNmenuName, (String) "Menu",
    NULL);
  if (mono) set_mono(file_menu_cmd);

  file_menu = XtVaCreatePopupShell("Menu",
    simpleMenuWidgetClass, file_menu_cmd,
    NULL);
  if (mono) set_mono(file_menu);

  file_menu_btn[0] = XtVaCreateManagedWidget("Command",
    smeBSBObjectClass, file_menu,
    XtNlabel, (String) file_menu_names[0],
    NULL);
  if (mono) set_mono(file_menu_btn[0]);
  file_menu_btn[1] = XtVaCreateManagedWidget("Command",
    smeBSBObjectClass, file_menu,
    XtNlabel, (String) file_menu_names[1],
    NULL);
  if (mono) set_mono(file_menu_btn[1]);

  XtAddCallback(file_menu_btn[0], XtNcallback,
    (XtCallbackProc) save_distance_matrix, (XtPointer) NULL);
  XtAddCallback(file_menu_btn[1], XtNcallback,
    (XtCallbackProc) Quit, (XtPointer) NULL);


/* Run */
  RUN = XtVaCreateManagedWidget("Command",
    toggleWidgetClass, runPanel,
    XtNlabel,       "Run MDS",
    XtNstate,        False,
    XtNfromVert,     file_menu_cmd,
    XtNvertDistance, 10,
    NULL);
  if (mono) set_mono(RUN);
  setToggleBitmap(RUN, False);
  XtAddCallback(RUN, XtNcallback,
    (XtCallbackProc) run_cback, (XtPointer) NULL);

/* Choice of metric or non-metric MDS */
  metricPanel = XtVaCreateManagedWidget("Panel",
    boxWidgetClass, runPanel,
    XtNhorizDistance, 5,
    XtNvertDistance, 5,
    XtNorientation, (XtOrientation) XtorientVertical,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNright, (XtEdgeType) XtChainLeft,
    XtNbottom, (XtEdgeType) XtChainTop,
    XtNfromVert,   RUN,
    NULL);
  metric_cmd[0] = XtVaCreateWidget("XGVToggle",
    toggleWidgetClass, metricPanel,
    XtNstate, (Boolean) True,
    XtNlabel, (String) "Metric MDS",
    NULL);
  metric_cmd[1] = XtVaCreateWidget("XGVToggle",
    toggleWidgetClass, metricPanel,
    XtNlabel, (String) "Non-metric MDS",
    XtNradioGroup, metric_cmd[0],
    NULL);
  XtManageChildren(metric_cmd, 2);

  setToggleBitmap(metric_cmd[0], True);
  setToggleBitmap(metric_cmd[1], False);
  XtAddCallback(metric_cmd[0], XtNcallback,
    (XtCallbackProc) setMetricCback, (XtPointer) NULL);
  XtAddCallback(metric_cmd[1], XtNcallback,
    (XtCallbackProc) setNonmetricCback, (XtPointer) NULL);


/* Step */
  STEP = XtVaCreateManagedWidget("Command",
    commandWidgetClass, runPanel,
    XtNlabel, (String) "Step",
    XtNfromVert,     metricPanel,
    XtNvertDistance, 10,
    NULL);
  if (mono) set_mono(STEP);
  XtAddCallback(STEP, XtNcallback,
    (XtCallbackProc) mds_iterate_cback, (XtPointer) NULL );

/* Reset */
  RESET = XtVaCreateManagedWidget("Command",
    commandWidgetClass, runPanel,
    XtNlabel,     "Reset",
    XtNsensitive,  pos_orig.nrows != 0,
    XtNfromVert,   STEP,
    XtNvertDistance, 5,
    NULL);
  if (mono) set_mono(RESET);
  XtAddCallback(RESET, XtNcallback,
    (XtCallbackProc) reset_cback, (XtPointer) NULL);

  /*
   * If there was no position matrix passed in, the Reset button
   * should be insensitive.
  if (pos_orig.nrows == 0)
    XtVaSetValues(RESET, XtNsensitive, False, NULL);
  */

/* Scramble */
  SCRAM = XtVaCreateManagedWidget("Command",
    commandWidgetClass, runPanel,
    XtNlabel, (String) "Scramble",
    XtNfromVert,     RESET,
    XtNvertDistance, 5,
    NULL);
  if (mono) set_mono(SCRAM);
  XtAddCallback(SCRAM, XtNcallback,
    (XtCallbackProc) scramble_cback, (XtPointer) NULL);

/* Help */
  HELP = XtVaCreateManagedWidget("Command",
    commandWidgetClass, runPanel,
    XtNlabel, (String) "Help ...",
    XtNfromVert,     SCRAM,
    XtNvertDistance, 5,
    NULL);
  if (mono) set_mono(HELP);
  XtAddCallback(HELP, XtNcallback,
    (XtCallbackProc) xgvis_help_cback, (XtPointer) NULL);

/*
 * Define MDS Panel, to contain the controls for MDS.
*/

  mdsPanel = XtVaCreateManagedWidget("Panel",
    formWidgetClass, form0,
    XtNfromHoriz, (Widget) runPanel,
    XtNhorizDistance, 5,
    XtNvertDistance, 5,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNright, (XtEdgeType) XtChainLeft,
    XtNbottom, (XtEdgeType) XtChainTop,
    NULL);
  if (mono) set_mono(mdsPanel);

/* Number of dimensions */

    dims_left = XtVaCreateManagedWidget("Icon",
      labelWidgetClass, mdsPanel,
      XtNhorizDistance, 5,
      XtNvertDistance, 5,
      XtNinternalHeight, (Dimension) 0,
      XtNinternalWidth, (Dimension) 0,
      XtNborderColor, (Pixel) appdata.fg,
      XtNbitmap, (Pixmap) leftarr,
      NULL);
    if (mono) set_mono(dims_left);
    XtAddEventHandler(dims_left, ButtonPressMask, FALSE,
      (XtEventHandler) mds_dimsleft_cback, (XtPointer) NULL);

    sprintf(str, "Dim (k): 999");
    width = XTextWidth(panel_data.Font, str, strlen(str));
    sprintf(str, "Dim (k): %d", mds_dims);
    mds_dims_label = XtVaCreateManagedWidget("Label",
        asciiTextWidgetClass, mdsPanel,
        XtNhorizDistance, 0,
        XtNvertDistance, 5,
        XtNstring, (String) str,
        XtNdisplayCaret, (Boolean) False,
        XtNwidth, (Dimension) width,
        /*XtNheight, (Dimension) 16,*/
        XtNfromHoriz, dims_left,
        NULL);
    if (mono) set_mono(mds_dims_label);

    dims_right = XtVaCreateManagedWidget("Icon",
      labelWidgetClass, mdsPanel,
      XtNfromHoriz, mds_dims_label,
      XtNhorizDistance, 0,
      XtNvertDistance, 5,
      XtNinternalHeight, (Dimension) 0,
      XtNinternalWidth, (Dimension) 0,
      XtNborderColor, (Pixel) appdata.fg,
      XtNbitmap, (Pixmap) rightarr,
      NULL);
    if (mono) set_mono(dims_right);
    XtAddEventHandler(dims_right, ButtonPressMask, FALSE,
      (XtEventHandler) mds_dimsright_cback, (XtPointer) NULL);

    /* first define width of label fields */
    sprintf(str, "Minkowski norm (m): %3.2f", 2.99);
    width = XTextWidth(panel_data.Font, str, strlen(str));

    /* Stepsize label and scrollbar for mds method */

    sprintf(str, "Stepsize: %3.3f", 999.999);
    sprintf(str, "Stepsize: %3.3f", mds_stepsize);

    mds_stepsize_label = XtVaCreateManagedWidget("Label",
        asciiTextWidgetClass, mdsPanel,
        XtNstring, (String) str,
        XtNdisplayCaret, (Boolean) False,
        XtNfromVert, (Widget) mds_dims_label,
        XtNwidth, (Dimension) width,
        NULL);
    if (mono) set_mono(mds_stepsize_label);

    mds_stepsize_sbar = XtVaCreateManagedWidget("Scrollbar",
        scrollbarWidgetClass, mdsPanel,
        XtNorientation, (XtOrientation) XtorientHorizontal,
        XtNfromVert, (Widget) mds_stepsize_label,
        XtNvertDistance, 0,
        XtNwidth, (Dimension) width,
        NULL);
    if (mono) set_mono(mds_stepsize_sbar);
    XawScrollbarSetThumb(mds_stepsize_sbar,
     (float) pow((double)(mds_stepsize), .5), -1.);
    XtAddCallback(mds_stepsize_sbar, XtNjumpProc,
        (XtCallbackProc) mds_stepsize_cback, (XtPointer) NULL);

    /* Exponent of distance matrix for mds method */

    sprintf(str, "Power (p) of D: %3.3f", 2.9);
    sprintf(str, "Power (p) of D: %3.1f", mds_power);

    mds_power_label = XtVaCreateManagedWidget("Label",
        asciiTextWidgetClass, mdsPanel,
        XtNstring, (String) str,
        XtNdisplayCaret, (Boolean) False,
        XtNfromVert, (Widget) mds_stepsize_sbar,
        XtNwidth, width,
        NULL);
    if (mono) set_mono(mds_power_label);

    mds_power_sbar = XtVaCreateManagedWidget("Scrollbar",
        scrollbarWidgetClass, mdsPanel,
        XtNorientation, (XtOrientation) XtorientHorizontal,
        XtNfromVert, (Widget) mds_power_label,
        XtNvertDistance, 0,
        XtNwidth, (Dimension) width,
        NULL);
    if (mono) set_mono(mds_power_sbar);
    XawScrollbarSetThumb(mds_power_sbar, mds_power/6.0, -1.);
    XtAddCallback(mds_power_sbar, XtNjumpProc,
        (XtCallbackProc) mds_power_cback, (XtPointer) NULL);

    /* Exponent of weights wij=(Dij^p)^r for mds method */
    sprintf(str, "Weight power (q): %3.3f", 2.9);
    sprintf(str, "Weight power (q): %3.1f", mds_weightpow);

    mds_weightpow_label = XtVaCreateManagedWidget("Label",
        asciiTextWidgetClass, mdsPanel,
        XtNstring, (String) str,
        XtNdisplayCaret, (Boolean) False,
        XtNfromVert, (Widget) mds_power_sbar,
        XtNwidth, width,
        NULL);
    if (mono) set_mono(mds_weightpow_label);

    mds_weightpow_sbar = XtVaCreateManagedWidget("Scrollbar",
      scrollbarWidgetClass, mdsPanel,
      XtNorientation, (XtOrientation) XtorientHorizontal,
      XtNfromVert, (Widget) mds_weightpow_label,
      XtNvertDistance, 0,
      XtNwidth, (Dimension) width,
      NULL);
    if (mono) set_mono(mds_weightpow_sbar);
    /* range should be -4 to +4 */
    XawScrollbarSetThumb(mds_weightpow_sbar, mds_weightpow/8.0+0.5, -1.);
    XtAddCallback(mds_weightpow_sbar, XtNjumpProc,
      (XtCallbackProc) mds_weightpow_cback, (XtPointer) NULL);

    /* Minkowski norm label and scrollbar for mds method */
    sprintf(str, "Minkowski norm (m): %3.3f", 2.99);
    sprintf(str, "%s: %3.1f", "Minkowski norm (m)", mds_lnorm);

    mds_lnorm_label = XtVaCreateManagedWidget("Label",
        asciiTextWidgetClass, mdsPanel,
        XtNstring, (String) str,
        XtNdisplayCaret, (Boolean) False,
        XtNfromVert, (Widget) mds_weightpow_sbar,
        XtNwidth, width,
        NULL);
    if (mono) set_mono(mds_lnorm_label);

    mds_lnorm_sbar = XtVaCreateManagedWidget("Scrollbar",
        scrollbarWidgetClass, mdsPanel,
        XtNorientation, (XtOrientation) XtorientHorizontal,
        XtNfromVert, (Widget) mds_lnorm_label,
        XtNvertDistance, 0,
        XtNwidth, (Dimension) width,
        NULL);
    if (mono) set_mono(mds_lnorm_sbar);

    /* Range: 1:6 */
    XawScrollbarSetThumb(mds_lnorm_sbar, (mds_lnorm-1.0)/5.0, -1.);
    XtAddCallback(mds_lnorm_sbar, XtNjumpProc,
        (XtCallbackProc) mds_lnorm_cback, (XtPointer) NULL);

/* Turn grouping on and off */
/*
    mds_group_label = XtVaCreateManagedWidget("Label",
      labelWidgetClass, mdsPanel,
      XtNlabel, (String) "MDS within subsets:",
      XtNdisplayCaret, (Boolean) False,
      XtNfromVert, (Widget) mds_lnorm_sbar,
      NULL);
    if (mono) set_mono(mds_group_label);
    mds_group_cmd = XtVaCreateManagedWidget("Command",
      toggleWidgetClass, mdsPanel,
      XtNfromVert, mds_group_label,
      XtNvertDistance, 0,
      XtNstate, (Boolean) mds_group,
      XtNlabel, (String) "use brush groups",
      NULL);
    if (mono) set_mono(mds_group_cmd);
    XtAddCallback(mds_group_cmd, XtNcallback,
      (XtCallbackProc) mds_group_cback,
      (XtPointer) NULL );
*/

/* Turn on and off the use of labels */
/*
    mds_casewise_label = XtVaCreateManagedWidget("Label",
      labelWidgetClass, mdsPanel,
      XtNlabel, (String) "MDS casewise:",
      XtNfromVert, (Widget) mds_group_cmd,
      NULL);
    if (mono) set_mono(mds_casewise_label);
    mds_casewise_cmd = XtVaCreateManagedWidget("Command",
      toggleWidgetClass, mdsPanel,
      XtNfromVert, mds_casewise_label,
      XtNvertDistance, 0,
      XtNstate, (Boolean) mds_casewise,
      XtNlabel, (String) "use current glyph&color",
      NULL);
    if (mono) set_mono(mds_casewise_cmd);
    XtAddCallback(mds_casewise_cmd, XtNcallback,
      (XtCallbackProc) mds_casewise_cback,
      (XtPointer) NULL );
*/

    mds_group_ind = none;
    build_labelled_menu(&group_menu_box, &group_menu_lab, group_menu_str,
      &group_menu_cmd,
      &group_menu, group_menu_btn,
      group_menu_btn_label, group_menu_btn_label, /* no nicknames */
      NGROUPBTNS, mds_group_ind, mdsPanel, mds_lnorm_sbar,
      XtorientVertical, panel_data.Font, "MDSGroups", NULL);
      /* Pack the box more tightly */
      XtVaSetValues(group_menu_box, 
        XtNborderWidth,  0,
        XtNhSpace,       0,
        XtNvSpace,       0,
        NULL);
    for (j=0; j<NGROUPBTNS; j++)
      XtAddCallback(group_menu_btn[j],  XtNcallback,
        (XtCallbackProc) set_mds_group_cback, (XtPointer) NULL);

/* Button to launch xgobi to contain diagnostic data */
    mds_launch_label = XtVaCreateManagedWidget("Label",
      labelWidgetClass, mdsPanel,
      XtNlabel,    "Shepard diagram:",
      /*XtNfromVert, mds_casewise_cmd,*/
      XtNfromVert, group_menu_box,
      NULL);
    if (mono) set_mono(mds_launch_label);
    mds_launch_panel = XtVaCreateManagedWidget("Panel",
      formWidgetClass, mdsPanel,
      XtNfromVert,     mds_launch_label,
      XtNvertDistance, 0,
      NULL);
    sprintf(str, "%d", ndistances);
    dwidth = XTextWidth(panel_data.Font, str, strlen(str)) + 6;
    mds_launch_ntxt = XtVaCreateManagedWidget("Text",
      asciiTextWidgetClass, mds_launch_panel,
      XtNleft,   (XtEdgeType) XtChainLeft,
      XtNright,  (XtEdgeType) XtChainLeft,
      XtNtop,    (XtEdgeType) XtChainTop,
      XtNbottom, (XtEdgeType) XtChainTop,
      XtNresizable, False,
      XtNeditType,  (int) XawtextEdit,
      XtNresize,    (XawTextResizeMode) XawtextResizeWidth,
      XtNstring,    str,
      XtNwidth,     dwidth,
      NULL);
    if (mono) set_mono(mds_launch_ntxt);
    sprintf(str, "of %d dists", ndistances);
    mds_launch_nlbl = XtVaCreateManagedWidget("Label",
      labelWidgetClass, mds_launch_panel,
      XtNlabel,     str,
      XtNfromHoriz, mds_launch_ntxt,
      NULL);
    if (mono) set_mono(mds_launch_nlbl);
    mds_launch_cmd = XtVaCreateManagedWidget("Command",
      commandWidgetClass, mdsPanel,
      XtNlabel,        "launch xgobi ...",
      XtNfromVert,     mds_launch_panel,
      XtNvertDistance, 0,
      NULL);
    XtAddCallback(mds_launch_cmd, XtNcallback,
      (XtCallbackProc) mds_launch_cback,
      (XtPointer) NULL);
    if (mono) set_mono(mds_launch_cmd);

/* Button to control menu of distance matrix types */
/*
    dist_cmd = XtVaCreateManagedWidget("Command",
      commandWidgetClass, mdsPanel,
      XtNlabel, (String) "distance metric ...",
      XtNfromVert, (Widget) mds_group_cmd,
      XtNhorizDistance, 5,
      XtNvertDistance, 10,
      NULL);
    if (mono) set_mono(dist_cmd);
    XtAddCallback(dist_cmd, XtNcallback,
     (XtCallbackProc) PopUpDistMenu, (XtPointer) NULL);

    dist_popup = XtVaCreatePopupShell("DistMenu",
      transientShellWidgetClass, dist_cmd,
      XtNinput, True,
      XtNtitle, "Dist metrics",
      NULL);
    if (mono) set_mono(dist_popup);

    dist_mgr = XtVaCreateManagedWidget("Panel",
      formWidgetClass, dist_popup,
      NULL);
    if (mono) set_mono(dist_mgr);

    dist_types[USER_SUPPLIED] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNstate, (Boolean) True,
        XtNlabel, (String) "D is user-supplied",
        NULL);
    if (mono) set_mono(dist_types[USER_SUPPLIED]);

    dist_types[LINK] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Link distances",
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNfromVert, (Widget) dist_types[USER_SUPPLIED],
        NULL);
    if (mono) set_mono(dist_types[LINK]);

    dist_types[ADJACENCY] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Adjacency matrix",
        XtNfromVert, (Widget) dist_types[LINK],
        NULL);
    if (mono) set_mono(dist_types[ADJACENCY]);

    dist_types[EUCLIDIAN] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Euclidian",
        XtNfromVert, (Widget) dist_types[ADJACENCY],
        NULL);
    if (mono) set_mono(dist_types[EUCLIDIAN]);

    dist_types[MANHATTAN] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Manhattan",
        XtNfromVert, (Widget) dist_types[EUCLIDIAN],
        NULL);
    if (mono) set_mono(dist_types[MANHATTAN]);

    dist_types[MAHALANOBIS] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Mahalanobis",
        XtNfromVert, (Widget) dist_types[MANHATTAN],
        NULL);
    if (mono) set_mono( dist_types[MAHALANOBIS]);

    dist_types[DOTPROD] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Dot Product",
        XtNfromVert, (Widget) dist_types[MAHALANOBIS],
        NULL);
    if (mono) set_mono(dist_types[DOTPROD]);

    dist_types[COSDIST] = XtVaCreateWidget("Command",
        toggleWidgetClass, dist_mgr,
        XtNradioGroup, (Widget) dist_types[USER_SUPPLIED],
        XtNstate, (Boolean) False,
        XtNlabel, (String) "Cosine",
        XtNfromVert, (Widget) dist_types[DOTPROD],
        NULL);
    if (mono) set_mono( dist_types[COSDIST]);

    XtManageChildren(dist_types, NDISTTYPES-1);

    apply_dist_cmd = XtVaCreateManagedWidget("Command",
        commandWidgetClass, dist_mgr,
        XtNlabel, (String) "Apply",
        XtNfromVert, (Widget) dist_types[COSDIST],
        XtNvertDistance, 15,
        NULL);
    if (mono) set_mono(apply_dist_cmd);
    XtAddCallback(apply_dist_cmd, XtNcallback,
        (XtCallbackProc) choose_dist_cback, (XtPointer) NULL );

    dist_types[CLOSE] = XtVaCreateManagedWidget("Command",
        commandWidgetClass, dist_mgr,
        XtNlabel, (String) "Close",
        XtNfromVert, (Widget) dist_types[COSDIST],
        XtNvertDistance, 15,
        XtNfromHoriz, (Widget) apply_dist_cmd,
        NULL);
    if (mono) set_mono(dist_types[CLOSE]);
    XtAddCallback(dist_types[CLOSE], XtNcallback,
        (XtCallbackProc) PopDownDistMenu, (XtPointer) NULL );
*/

  build_stress_plotwin(form0, NULL, mdsPanel);  /* parent, href, vref */

  build_dissim_plotwin(form0);  /* parent */
}

void
update_shepard_labels(int maxn) {
  String nstr;
  int n;
  char strtmp[64];

  /*
   * Set the value of the user-settable text to the maximum
   * of the current value and the maximum value.
  */
  XtVaGetValues(mds_launch_ntxt, XtNstring, &nstr, NULL);
  if (strlen(nstr) == 0)
    n = maxn;
  else {
    n = atoi(nstr);
    n = MIN(maxn, n);
  }
  sprintf(strtmp, "%d", n);
  XtVaSetValues(mds_launch_ntxt, XtNstring, strtmp, NULL);

  /*
   * Set the value of the label to the maximum value
  */
  sprintf(strtmp, "of %d dists", maxn);
  XtVaSetValues(mds_launch_nlbl, XtNlabel, strtmp, NULL);
}

