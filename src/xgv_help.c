#include <stdlib.h>
#include <stdio.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"
#include "xgvis.h"

#include "../bitmaps/stress.xbm"

/*
 * This is all very well, but it doesn't allow the inclusion
 * of the bitmap.
*/

/*
XtCallbackProc
xgvis_help_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  help(XtParent(XtParent(w)), "xgvis", &xgobi);
}
*/


/* ARGSUSED */
XtCallbackProc
xgvis_help_done_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  XtDestroyWidget(XtParent(XtParent(w)));
}

/* ARGSUSED */
XtCallbackProc
xgvis_help_cback(Widget w, XtPointer client_data, XtPointer callback_data)
{
  char fname[100];
  char message[MSGLENGTH];
  char *xgobidir;
  FILE *fp;
  extern xgobidata xgobi;  /* defined in xgvis.c */
  Dimension /*width,*/ height;

  xgobidir = getenv("XGOBID");
  if (xgobidir == NULL || strlen(xgobidir) == 0)
  {
    xgobidir = (char *) XtMalloc((Cardinal) 150 * sizeof(char));
    (void) strcpy(xgobidir, DEFAULTDIR);
    if (xgobidir == NULL || strlen(xgobidir) == 0)
    {
      sprintf(message,
       "XGOBID is not defined in your environment, and\n");
      strcat(message,
       "DEFAULTDIR is not defined in the XGobi Makefile;\n");
      strcat(message,
        "see the person who installed XGobi for help.\n");
      show_message(message, &xgobi);
      return;
    }
    else
    {
      (void) strcpy(fname, xgobidir);
      XtFree((XtPointer) xgobidir);
    }
  }
  else
  {
    (void) strcpy(fname, xgobidir);
  }

/*
 * Check that the file is good.
*/
  (void) strcat(fname, "/help/xgvis");
  if ((fp = fopen(fname, "r")) == NULL)
  {
    sprintf(message,
      "Unable to open %s.\n", fname);
    strcat(message,
      "Is the shell variable XGOBID the name of the directory\n");
    strcat(message,
      "which contains the help subdirectory? \n");
    show_message(message, &xgobi);
    return((XtCallbackProc) 0);
  }
  else {
    Widget hpopup, hframe, hform, hfunc, htext, hdone;
    Screen *scrn;
    Pixmap func_pix;
    Pixel funcfg, funcbg;

    /*
     * Read in the bitmap of the function.
    */
    scrn = XtScreen(shell);
    funcbg = WhitePixelOfScreen(scrn);
    funcfg = BlackPixelOfScreen(scrn);

    func_pix = XCreatePixmapFromBitmapData(display,
      RootWindowOfScreen(scrn),
      stress_bits, stress_width, stress_height,
      funcfg, funcbg,
      depth);

    /* 80 columns x 20 rows, I hope */
    /* width = 85 * XTextWidth(appdata.helpFont, "M", 1) ; */
    height = 20 * FONTHEIGHT(appdata.helpFont);

  /*
   * Create the popup itself.
  */
    hpopup = XtVaCreatePopupShell("Help",
      topLevelShellWidgetClass, shell,
      XtNtitle, (String) "XGvis Help Window",
      XtNiconName, (String) "XGvis Help Window",
      NULL);
    if (mono) set_mono(hpopup);
  /*
   * Create the paned widget.
  */
    hframe = XtVaCreateManagedWidget("Form",
      panedWidgetClass, hpopup,
      XtNorientation, (XtOrientation) XtorientVertical,
      NULL);
    if (mono) set_mono(hframe);

  /*
   * Create a form widget to hold the label and text.
  */
    hform = XtVaCreateManagedWidget("Help",
      formWidgetClass, hframe,
      XtNbackground, funcbg,
      NULL);

    hfunc = XtVaCreateManagedWidget("Help",
      labelWidgetClass, hform,
      XtNinternalHeight, (Dimension) 0,
      XtNinternalWidth, (Dimension) 0,
      XtNbitmap, (Pixmap) func_pix,
      XtNborderWidth, 0,
      NULL);

  /*
   * Create the text widget.
  */
    htext = XtVaCreateManagedWidget("Text",
      asciiTextWidgetClass, hframe,
      XtNallowResize, (Boolean) True,
      XtNshowGrip, (Boolean) False,
      XtNtype, (XawAsciiType) XawAsciiFile,
      XtNstring, (String) fname,
      XtNscrollVertical, (XawTextScrollMode) XawtextScrollWhenNeeded,
      XtNdisplayCaret, (Boolean) False,
      XtNfont, (XFontStruct *) appdata.helpFont,
/*
      XtNwidth, (Dimension) width,
*/
      XtNheight, (Dimension) height,
      XtNfromVert, hfunc,
      NULL);
    if (mono) set_mono(htext);

  /*
   * Create the Done button.
  */
    hdone = XtVaCreateManagedWidget("Done",
      commandWidgetClass, hframe,
      XtNshowGrip, (Boolean) False,
      XtNskipAdjust, (Boolean) True,
      XtNlabel, (String) "Click here to dismiss",
      NULL);
    if (mono) set_mono(hdone);

    XtAddCallback(hdone, XtNcallback,
      (XtCallbackProc) xgvis_help_done_cback, (XtPointer) NULL);

    XtPopup(hpopup, XtGrabNone);
    XRaiseWindow(display, XtWindow(hpopup));

    set_wm_protocols(hpopup);

    fclose(fp);
  }
}
