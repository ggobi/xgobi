/* inference.c */
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

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"

static Widget infpopup = (Widget) NULL;
static Boolean infer_vgroup = True;

/* ARGSUSED */
static XtCallbackProc
close_infer_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  XtPopdown(infpopup);
}

/* ARGSUSED */
static XtCallbackProc
infer_vgroups_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  infer_vgroup = !infer_vgroup;
  setToggleBitmap(w, infer_vgroup);
}

/* ARGSUSED */
XtCallbackProc
open_infer_popup_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  static Boolean initd = False;
  Widget mframe, mform;
  Widget infer_vgroups_cmd, close_cmd;

  if (!initd)
  {
    Dimension width, height;
    Position x, y;

    XtVaGetValues(w,
      XtNwidth, &width,
      XtNheight, &height, NULL);
    XtTranslateCoords(w,
      (Position) (width/2), (Position) (height/2), &x, &y);

    /*
     * Create the missing data popup
    */
    infpopup = XtVaCreatePopupShell("Inference",
      /*transientShellWidgetClass, XtParent(w),*/
      topLevelShellWidgetClass, XtParent(w),
      XtNinput,            (Boolean) True,
      XtNallowShellResize, (Boolean) True,
      XtNtitle,            (String) "Inference",
      XtNiconName,         (String) "Inference",
      XtNx,                x,
      XtNy,                y,
      NULL);
    if (mono) set_mono(infpopup);

    /*
     * Create a paned widget so the 'Click here ...'
     * can be all across the bottom.
    */
    mframe = XtVaCreateManagedWidget("Form",
      panedWidgetClass, infpopup,
      XtNorientation, (XtOrientation) XtorientVertical,
      NULL);
    /*
     * Create the form widget.
    */
    mform = XtVaCreateManagedWidget("Inference",
      formWidgetClass, mframe,
      NULL);
    if (mono) set_mono(mform);

    infer_vgroups_cmd = (Widget) CreateToggle(xg, "Use vgroups",
      True, (Widget) NULL, (Widget) NULL, (Widget) NULL,
      infer_vgroup, ANY_OF_MANY, mform, "Inference");
    XtManageChild(infer_vgroups_cmd);
    XtAddCallback(infer_vgroups_cmd, XtNcallback,
      (XtCallbackProc) infer_vgroups_cback, (XtPointer) xg);

    close_cmd = XtVaCreateManagedWidget("Close",
      commandWidgetClass, mframe,
      XtNshowGrip, (Boolean) False,
      XtNskipAdjust, (Boolean) True,
      XtNlabel, (String) "Click here to dismiss",
      NULL);
    if (mono) set_mono(close_cmd);
    XtAddCallback(close_cmd, XtNcallback,
      (XtCallbackProc) close_infer_cback, (XtPointer) xg);
  }

  XtPopup(infpopup, (XtGrabKind) XtGrabNone);
  XRaiseWindow(display, XtWindow(infpopup));

  if (!initd)
  {
    set_wm_protocols(infpopup);
    initd = True;
  }
}
