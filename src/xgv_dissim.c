#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"
#include "xgvis.h"

#define THRWIN_VMARGIN 5
#define THRWIN_HMARGIN 15
/* This reserved space for a bar with grips */
#define THRWIN_GRIP_SPACE 20
#define THRWIN_GRIP_SIZE 10  /* should be even number */

static int min_grip_pos;
static int max_grip_pos;
static int grip_pos[2];
static Boolean grip_pressed[2];
Boolean *bars_included;

Boolean initd = false;

Widget thr_wksp;
Dimension thr_width = 250, thr_height = 100;
int thr_pwidth, thr_xmin, thr_xmax;
int thr_pheight, thr_ymin, thr_ymax;
int thr_bwidth, thr_nbins, *thr_bins;
XRectangle *thr_bars;
Drawable thr_window;
Pixmap thr_pixmap = (Pixmap) NULL;

/*
 * These have the same meaning as the values used in
 * mds.c, that is, the mds_threshold raised to the power
 * mds_power.
*/

static double dthresh_high = 1.0, dthresh_low = 0.0;

void
make_thr_barchart(void)
{
/*
 * Scale the counts into (xmin, ymin), (xmax, ymax);
 * construct XRectangle structures for drawing.
*/
  int maxcount = 0;
  int barheight;
  int i;

  thr_bars = (XRectangle *) XtRealloc((char *) thr_bars,
    thr_nbins * sizeof(XRectangle));

  /* find maximum count */
  for (i=0; i<thr_nbins; i++)
    if (thr_bins[i] > maxcount) maxcount = thr_bins[i];

  for (i=0; i<thr_nbins; i++) {
    barheight = (int) (thr_pheight * (double)thr_bins[i]/(double)maxcount);
    /* x,y: upper lefthand corner */
    thr_bars[i].x = thr_xmin + (i * thr_bwidth);
    thr_bars[i].y = thr_ymax - barheight;
    thr_bars[i].width = thr_bwidth;
    thr_bars[i].height = barheight;
  }

}

void
make_thr_pixmap(void)
{
  thr_pixmap = XCreatePixmap(display, thr_window,
    thr_width, thr_height, depth);
}

void
clear_thr_pixmap(void)
{
  XFillRectangle(display, thr_pixmap, clear_GC,
    0, 0, thr_width, thr_height);
}

void
copy_thr_pixmap(void)
{
  /* copy the pixmap to the screen */
  XCopyArea(display, thr_pixmap, thr_window, copy_GC,
    0, 0, thr_width, thr_height, 0, 0 );
}

void
draw_grip_control(void) {
  XGCValues *gcv, gcv_inst;
  XRectangle grips[2];
  int ypos;
  static Boolean initd = false;

  /* The grip positions are at the center of the grips */

  /* Suppose I allow the grips to go a bit into the margins ... */
  min_grip_pos = thr_xmin - THRWIN_HMARGIN/2;
  max_grip_pos = thr_xmax + THRWIN_HMARGIN/2;

  if (!initd) {
    grip_pos[0] = min_grip_pos;
    grip_pos[1] = max_grip_pos;
    initd = true;
  }

  XGetGCValues(display, varpanel_xor_GC, GCCapStyle|GCJoinStyle, &gcv_inst);
  gcv = &gcv_inst;
  XSetLineAttributes(display, varpanel_xor_GC, 2, LineSolid,
    gcv->cap_style, gcv->join_style);

  ypos = thr_height - THRWIN_VMARGIN - THRWIN_GRIP_SIZE/2;

/*
  XDrawLine(display, thr_pixmap, copy_GC,
    thr_xmin, ypos, thr_xmax, ypos);
*/
  XDrawLine(display, thr_pixmap, copy_GC,
    min_grip_pos, ypos, max_grip_pos, ypos);

  XSetLineAttributes(display, varpanel_xor_GC, 1, LineSolid,
    gcv->cap_style, gcv->join_style);

  grips[0].x = grip_pos[0] - THRWIN_GRIP_SIZE/2;
  grips[0].y = ypos - THRWIN_GRIP_SIZE/2;
  grips[0].width = grips[0].height = THRWIN_GRIP_SIZE;

  grips[1].x = grip_pos[1] - THRWIN_GRIP_SIZE/2;
  grips[1].y = ypos - THRWIN_GRIP_SIZE/2;
  grips[1].width = grips[1].height = THRWIN_GRIP_SIZE;
  
  XFillRectangles(display, thr_pixmap, copy_GC,
    grips, 2);

}

void
draw_thr(void)
{
  int i;

  if (thr_pixmap == (Pixmap) NULL)
    return;

  clear_thr_pixmap();

  XSetForeground(display, copy_GC, plotcolors.fg);

  for (i=0; i<thr_nbins; i++) {

    if (bars_included[i]) {
      XDrawRectangle(display, thr_pixmap, copy_GC,
        thr_bars[i].x, thr_bars[i].y, thr_bars[i].width, thr_bars[i].height);
      XFillRectangle(display, thr_pixmap, copy_GC,
        thr_bars[i].x, thr_bars[i].y, thr_bars[i].width, thr_bars[i].height);
    } else {
      if (i>0)  /* leading vertical edge */
        XDrawLine(display, thr_pixmap, copy_GC,
          thr_bars[i].x, thr_bars[i-1].y,
          thr_bars[i].x, thr_bars[i].y);

      /* horizontal line */
      XDrawLine(display, thr_pixmap, copy_GC,
        thr_bars[i].x, thr_bars[i].y,
        thr_bars[i].x+thr_bars[i].width, thr_bars[i].y);

      if (i<thr_nbins-1)  /* trailing vertical edge */
        XDrawLine(display, thr_pixmap, copy_GC,
          thr_bars[i].x+thr_bars[i].width, thr_bars[i].y,
          thr_bars[i].x+thr_bars[i].width, thr_bars[i+1].y);
    }
  }

  XDrawLine(display, thr_pixmap, copy_GC,
    thr_bars[thr_nbins-1].x+thr_bars[thr_nbins-1].width,
      thr_bars[thr_nbins-1].y,
    thr_bars[thr_nbins-1].x+thr_bars[thr_nbins-1].width,
      thr_ymax);

  draw_grip_control();

  copy_thr_pixmap();
}

/* ARGSUSED */
XtCallbackProc
thr_expose_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  copy_thr_pixmap();
}

void
set_threshold(void)
{
  int k;

  for (k=0; k<thr_nbins; k++) {

    /* bars are included if their left edge is greater than
       the left grip's center and their right edge is less than
       the right grip's center.
    */
    if (thr_bars[k].x >= grip_pos[0] &&
        thr_bars[k].x + thr_bars[k].width <= grip_pos[1])
    {
      bars_included[k] = true;
    }
    else
      bars_included[k] = false;
  }

  /* Now specify the thresholds in data terms */

  /*
   * Scale the grip position onto [0,1]
  */

  dthresh_low = MAX(0, (double)(grip_pos[0] - thr_xmin) /
    (double) (thr_xmax - thr_xmin) );
  dthresh_high = MIN(1, (double)(grip_pos[1] - thr_xmin) /
    (double) (thr_xmax - thr_xmin) );

  mds_threshold_low = pow(dthresh_low, 1.0/mds_power);
  mds_threshold_high = pow(dthresh_high, 1.0/mds_power);

}

/* ARGSUSED */
XtEventHandler
thr_motion(Widget w, XtPointer client_data, XEvent *evnt, Boolean *cont)
{
  XMotionEvent *xmotion = (XMotionEvent *) evnt;
  if (grip_pressed[0] &&
      xmotion->x + THRWIN_GRIP_SIZE < grip_pos[1] &&
      xmotion->x >= min_grip_pos)
  {
    grip_pos[0] = xmotion->x;
  }
  else if (grip_pressed[1] &&
           xmotion->x > grip_pos[0] + THRWIN_GRIP_SIZE &&
           xmotion->x <= max_grip_pos)
  {
    grip_pos[1] = xmotion->x;
  }

  set_threshold();
  draw_thr();
}

/* ARGSUSED */
XtEventHandler
thr_button(Widget w, XtPointer client_data, XEvent *evnt)
{
/*
 * Using only the x coordinate, see if we're on top of one
 * of the two grips.
*/
  XButtonEvent *xbutton = (XButtonEvent *) evnt;
  if (xbutton->type == ButtonRelease)
    grip_pressed[0] = grip_pressed[1] = false;
  else {
    if (xbutton->x >= grip_pos[0] - THRWIN_GRIP_SIZE/2 &&
        xbutton->x <= grip_pos[0] + THRWIN_GRIP_SIZE/2)
    {
      grip_pressed[0] = true;
    }
    else if (xbutton->x >= grip_pos[1] - THRWIN_GRIP_SIZE/2 &&
             xbutton->x <= grip_pos[1] + THRWIN_GRIP_SIZE/2)
    {
      grip_pressed[1] = true;
    }
  }

}

void
reset_thr_bins(void)
{
  int i, k;
  double edge;

  thr_bwidth = 5;  /* Try a fixed binwidth of 5 for a minute */

  /* width of plot divided by binwidth */
  thr_nbins = (int) ( (double)thr_pwidth / (double)thr_bwidth );

  /* fix any rounding errors */
  thr_pwidth = thr_nbins * thr_bwidth; 
  thr_xmax = thr_xmin + thr_pwidth;

  thr_bins = (int *) XtRealloc((char *) thr_bins, thr_nbins * sizeof(int));

  for (i=0; i<thr_nbins; i++)
    thr_bins[i] = 0;

  i = k = 0;
  edge = (double)(k+1) / (double) thr_nbins ;
  while (i < ndistances) {  /* length of distance_vector */
    while (distance_vector_sort[i] >= edge && k<(thr_nbins-1)) {
      k++;
      edge = (double)k / (double) thr_nbins ;
    }
    if (k>thr_nbins-1)
      fprintf(stderr, "warning: k %d thr_nbins %d\n", k, thr_nbins);
    (thr_bins[k])++;
    i++;
  }
}

void
set_pheight(void)
{
  thr_ymin = THRWIN_VMARGIN;
  /*
   * subtract two lower margins, above and below the grip control,
   * and the height of the grip control.
  */
  thr_ymax = thr_height - 2*THRWIN_VMARGIN - THRWIN_GRIP_SPACE;
  if (thr_ymax < thr_ymin) {
    thr_ymin = 0;
    thr_ymax = thr_height;
  }
  thr_pheight = thr_ymax - thr_ymin;
}
void
set_pwidth(void)
{
  thr_xmin = THRWIN_HMARGIN;
  thr_xmax = thr_width - THRWIN_HMARGIN;
  if (thr_xmax < thr_xmin) {
    thr_xmin = 0;
    thr_xmax = thr_width;
  }
  thr_pwidth = thr_xmax - thr_xmin;

}

/* ARGSUSED */
XtCallbackProc
thr_resize_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  XtVaGetValues(thr_wksp,
    XtNwidth, &thr_width,
    XtNheight, &thr_height,
    NULL);

  set_pwidth();
  set_pheight();
  reset_thr_bins();

  if (thr_nbins > 0) {
    /* Free the pixmap; create a new one */
    XFreePixmap(display, thr_pixmap);
    make_thr_pixmap();

    make_thr_barchart();
    draw_thr();
  }
}

void
init_thresholding(void)
{
  int i=0;

  thr_bars = (XRectangle *) XtMalloc(50 * sizeof(XRectangle));
  bars_included = (Boolean *) XtMalloc(50 * sizeof(Boolean));
  thr_bins = (int *) XtMalloc(50 * sizeof(int));

  for (i=0; i<50; i++) 
    bars_included[i] = true;
}

void
build_dissim_plotwin(Widget parent)
{
  Widget thr_label, thr_form;
  extern Widget stress_form;
  static String lbl = "Distribution of D^p";

  thr_width = MAX(thr_width,
      XTextWidth(panel_data.Font, lbl,
        strlen(lbl) + ASCII_TEXT_BORDER_WIDTH));

  init_thresholding();

  thr_form = XtVaCreateManagedWidget("Form",
    formWidgetClass, parent,
    XtNfromVert, stress_form,
    NULL);
  if (mono) set_mono(thr_form);

  thr_label = XtVaCreateManagedWidget("Label",
    labelWidgetClass, thr_form,
    XtNlabel, lbl,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNright, (XtEdgeType) XtChainRight,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNbottom, (XtEdgeType) XtChainTop,
    XtNwidth, thr_width,
    XtNborderWidth, 0,
    NULL);
  
  thr_wksp = XtVaCreateManagedWidget("Thresholding",
    labelWidgetClass, thr_form,
    XtNlabel, (String) "",
    XtNforeground, (Pixel) plotcolors.fg,
    XtNbackground, (Pixel) plotcolors.bg,
    XtNwidth, (Dimension) thr_width,
    XtNheight, (Dimension) thr_height,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNright, (XtEdgeType) XtChainRight,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNbottom, (XtEdgeType) XtChainBottom,
    XtNfromVert, thr_label,
    NULL);
  if (mono) set_mono(thr_wksp);

  XtAddEventHandler(thr_wksp, ExposureMask,
    FALSE, (XtEventHandler) thr_expose_cback, (XtPointer) NULL);
  XtAddEventHandler(thr_wksp, StructureNotifyMask,
    FALSE, (XtEventHandler) thr_resize_cback, (XtPointer) NULL);
  XtAddEventHandler(thr_wksp, ButtonPressMask | ButtonReleaseMask,
    FALSE, (XtEventHandler) thr_button, (XtPointer) NULL);
  XtAddEventHandler(thr_wksp, Button1MotionMask,
    FALSE, (XtEventHandler) thr_motion, (XtPointer) NULL);
}

void
init_dissim(void)
{
  set_pwidth();
  set_pheight();

  initd = true;

  mds_threshold_low = dthresh_low = 0;
  mds_threshold_high = dthresh_high = 1;

  reset_thr_bins();
  make_thr_barchart();

  thr_window = XtWindow( thr_wksp );
  make_thr_pixmap();

  draw_thr();
}

void
update_dissim_plot() {
  int k;
  /*
   * I can't explain this, but somehow this is being called
   * once before init_dissim has been called.
  */
  if (!initd)
    return;

  reset_thr_bins();

  /* Given that mds_power has changed, recalculate the dthresh values */

  dthresh_low = MAX(pow(mds_threshold_low, mds_power), 0.0);
  dthresh_high = MIN(pow(mds_threshold_high, mds_power), 1.0);

  /* And given that the dthresh values have changed, reset grip_pos[] */

  grip_pos[0] = (int)
    (dthresh_low * (double) (thr_xmax - thr_xmin) + (double) thr_xmin);
  grip_pos[1] = (int)
    (dthresh_high * (double) (thr_xmax - thr_xmin) + (double) thr_xmin);

  /* Given that the grip has changed, reset the included bars */

  for (k=0; k<thr_nbins; k++) {
    if (thr_bars[k].x >= grip_pos[0] &&
        thr_bars[k].x + thr_bars[k].width <= grip_pos[1])
      bars_included[k] = true;
    else
      bars_included[k] = false;
  }

  make_thr_barchart();
  draw_thr();
}
