/* move_points.c */
/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  the xgobi code and associated documentation,  provided  *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *     Copyright (c)  AT&T Labs (1998).                     *
 *  All Rights Reserved.                                    *
 *                                                          *
 *  We welcome your questions and comments, and request     *
 *  that you share any modifications with us.               *
 *                                                          *
 *      Deborah Swayne             Andreas Buja             *
 *    dfs@research.att.com     andreas@research.att.com     *
 *                                                          *
 ************************************************************/

#include <stdio.h>
#include "xincludes.h"
#include "xgobitypes.h"
#include "xgobivars.h"
#include "xgobiexterns.h"

int moving_point = -1;

/* choose a cursor at some point */

#define UP 0
#define DOWN 1
static Boolean buttonpos = UP;
static int last_moved = -1;
static int last_x = -1, last_y = -1;
lcoords eps;

static Widget mp_panel;
static Widget reset_all_cmd, reset_one_cmd;

static Boolean use_brush_groups = False;

static Widget mpdir_menu_label, mpdir_menu_cmd, mpdir_menu, mpdir_menu_btn[3];
static char *mpdir_menu_btn_label[] = {"Both", "Vert", "Horiz"};
static enum {both, vert, horiz} mpdir_type = both;

/********** Reverse pipeline in pipeline.c ********************/

void
init_point_moving(xgobidata *xg)
{
  xg->is_point_moving = False;
}

/* ARGSUSED */
XtEventHandler
mp_button(Widget w, xgobidata *xg, XEvent *evnt, Boolean *cont)
{
  XButtonEvent *xbutton = (XButtonEvent *) evnt;

  if (xbutton->button == 1) {
    if (xbutton->type == ButtonPress) {
      buttonpos = DOWN;
      last_moved = xg->nearest_point;
      last_x = xg->screen[xg->nearest_point].x;
      last_y = xg->screen[xg->nearest_point].y;
    }
    else if (xbutton->type == ButtonRelease) {
      buttonpos = UP;

      moving_point = -1;
    }
  }
}

void
move_pt_pipeline(int id, lcoords *eps, xgobidata *xg)
{
  extern void plane_to_world(xgobidata *, int, lcoords *);
  extern void world_to_raw(xgobidata *, int);

  /* run the earliest parts of the pipeline backwards */
  plane_to_world(xg, id, eps);
  world_to_raw(xg, id);
}

void
move_pt(int id, int x, int y, xgobidata *xg) {
  int i, k;

  if (mpdir_type == horiz || mpdir_type == both)
    xg->screen[id].x = x;
  if (mpdir_type == vert || mpdir_type == both)
    xg->screen[id].y = y;

  /* Move the selected point */
  screen_to_plane(xg, id, &eps,
    (mpdir_type == horiz || mpdir_type == both),
    (mpdir_type == vert || mpdir_type == both));
  move_pt_pipeline(id, &eps, xg);  /*-- eps won't be changed here --*/

  if (use_brush_groups) {
    if (xg->nclust > 1) {
      int cur_clust = (int) xg->raw_data[id][xg->ncols_used-1];

      /*
       * Move all points which belong to the same cluster
       * as the selected point.
      */
      for (i=0; i<xg->nrows_in_plot; i++) {
        k = xg->rows_in_plot[i];
        if (k == id)
          ;
        else {
          if ((int)xg->raw_data[k][xg->ncols_used-1] == cur_clust) {
            if (!xg->erased[k]) {   /* ignore erased values altogether */

              if (mpdir_type == horiz || mpdir_type == both) {
                xg->planar[k].x += eps.x;
              }
              if (mpdir_type == vert || mpdir_type == both) {
                xg->planar[k].y += eps.y;
              }
              move_pt_pipeline(k, &eps, xg);
            }
          }
        }
      }
    }
  }

  /* and now forward again, all the way ... */
  update_world(xg);
  world_to_plane(xg);
  plane_to_screen(xg);
  if (xg->is_brushing) {
    assign_points_to_bins(xg);
    if (xg->brush_mode == transient)
      reinit_transient_brushing(xg);
  }
  plot_once(xg);
  if (xg->is_cprof_plotting && xg->link_cprof_plotting)
    update_cprof_plot(xg);
}

void
move_points_proc(xgobidata *xg)
{
  int root_x, root_y;
  unsigned int kb;
  Window root, child;
  static int ocpos_x = 0, ocpos_y = 0;
  icoords cpos;
  static int inwindow = 1;
  int wasinwindow;
  Boolean pointer_moved = False;

  wasinwindow = inwindow;

/*
 * Find the nearest point, just as in identification, and highlight it
*/

/*
 * Get the current pointer position.
*/
  if (XQueryPointer(display, xg->plot_window, &root, &child,
            &root_x, &root_y, &cpos.x, &cpos.y, &kb))
  {
    inwindow = (0 < cpos.x && cpos.x < xg->max.x &&
                0 < cpos.y && cpos.y < xg->max.y) ;
    /*
     * If the pointer is inside the plotting region ...
    */
    if (inwindow)
    {
      /*
       * If the pointer has moved ...
      */
      if ((cpos.x != ocpos_x) || (cpos.y != ocpos_y)) {
        pointer_moved = True;
        ocpos_x = cpos.x;
        ocpos_y = cpos.y;
      }
    }
  }

  if (pointer_moved) {
    
    if (buttonpos == UP) {
      moving_point = -1;
      if ( (xg->nearest_point = find_nearest_point(&cpos, xg)) != -1) {
        quickplot_once(xg);
      }
    }
    else {
      /*
       * If the left button is down, move the point: compute the
       * data pipeline in reverse, (then run it forward again?) and
       * draw the plot.
      */
      if (buttonpos == DOWN && xg->nearest_point != -1) {
        moving_point = xg->nearest_point;
        move_pt(xg->nearest_point, cpos.x, cpos.y, xg);
      }
    }
  }

  if (!inwindow && wasinwindow) {
    /*
     * Don't draw the diamond if the pointer leaves the plot window.
    */
    xg->nearest_point = -1;
    quickplot_once(xg);
  }
}

static void
map_move_points(Boolean on) {
  if (on) 
    XtMapWidget(mp_panel);
  else
    XtUnmapWidget(mp_panel);
}

/* ARGSUSED */
static XtCallbackProc
reset_all_cback(Widget w, xgobidata *xg, XtPointer cb_data) {
/*
 * read the raw data back in here instead of having to go
 * to the file menu to accomplish it.
*/
  Boolean reset = False;

  if (xg->datafilename != "") {
    if (reread_dat(xg->datafilename, xg)) {
      plot_once(xg);
      reset = True;
    }
  }

  if (!reset)
    show_message(
     "Sorry, I can\'t reset because I can\'t reread the data file\n", xg);
}

/* ARGSUSED */
static XtCallbackProc
reset_one_cback(Widget w, xgobidata *xg, XtPointer cb_data) {
/*
 * move the most recently moved point back to its position
 * when the mouse went down
*/
  if (last_moved >=0 && last_moved < xg->nrows_in_plot) {
    if ((last_x > -1 && last_x < xg->plotsize.width) &&
        (last_y > -1 && last_y < xg->plotsize.height))
    {
      move_pt(last_moved, last_x, last_y, xg);
    }
  }
}

/* ARGSUSED */
static XtCallbackProc
use_groups_cback(Widget w, xgobidata *xg, XtPointer callback_data)
{
  use_brush_groups = !use_brush_groups;

  if (use_brush_groups) {
    if (xg->ncols_used < xg->ncols) {
      save_brush_groups(xg);
    }
  }
  setToggleBitmap(w, use_brush_groups);

}

/* ARGSUSED */
XtCallbackProc
set_mpdir_cback(Widget w, xgobidata *xg, XtPointer cb)
{
  int k;
  int id;

  for (k=0; k<3; k++) {
    if (w == mpdir_menu_btn[k]) {
      id = k;
      break;
    }
  }

  XtVaSetValues(mpdir_menu_cmd,
    XtNlabel, mpdir_menu_btn_label[id],
    NULL);

  if (id == 0)
    mpdir_type = both;
  else if (id == 1)
    mpdir_type = vert;
  else if (id == 2)
    mpdir_type = horiz;
}

void
make_move_points(xgobidata *xg) {
  Widget use_groups_cmd;
  Widget mpdir_box;
  int k;

  mp_panel = XtVaCreateManagedWidget("MvPtsPanel",
    boxWidgetClass, xg->box0,
    XtNleft, (XtEdgeType) XtChainLeft,
    XtNright, (XtEdgeType) XtChainLeft,
    XtNtop, (XtEdgeType) XtChainTop,
    XtNbottom, (XtEdgeType) XtChainTop,
    XtNmappedWhenManaged, (Boolean) False,
    NULL);
  if (mono) set_mono(mp_panel);

  build_labelled_menu(&mpdir_box, &mpdir_menu_label, "MotionDir:",
    &mpdir_menu_cmd, &mpdir_menu, mpdir_menu_btn,
    mpdir_menu_btn_label, mpdir_menu_btn_label,  /* no nicknames */
    3, 0, mp_panel, mpdir_box,
    XtorientHorizontal, appdata.font, "MovePts", xg);
  for (k=0; k<3; k++)
    XtAddCallback(mpdir_menu_btn[k],  XtNcallback,
      (XtCallbackProc) set_mpdir_cback, (XtPointer) xg);

  use_groups_cmd = CreateToggle(xg, "Use 'group' var",
    True, (Widget) NULL, (Widget) NULL, (Widget) NULL, use_brush_groups,
    ANY_OF_MANY, mp_panel, "MP_Groups");
  XtManageChild(use_groups_cmd);
  XtAddCallback(use_groups_cmd, XtNcallback,
    (XtCallbackProc) use_groups_cback, (XtPointer) xg);

  reset_all_cmd = CreateCommand(xg, "Reset all",
    True, (Widget) NULL, (Widget) NULL,
    mp_panel, "MP_Reset");
  XtManageChild(reset_all_cmd);
  if (xg->progname == NULL)
    XtSetSensitive(reset_all_cmd, False);
  XtAddCallback(reset_all_cmd, XtNcallback,
    (XtCallbackProc) reset_all_cback, (XtPointer) xg);

  reset_one_cmd = CreateCommand(xg, "Undo last",
    True, (Widget) NULL, (Widget) NULL,
    mp_panel, "MP_Reset");
  XtManageChild(reset_one_cmd);
  XtAddCallback(reset_one_cmd, XtNcallback,
    (XtCallbackProc) reset_one_cback, (XtPointer) xg);
}

void move_points_on (xgobidata *xg)
{
/*
 *  If this mode is currently selected, turn it off.
*/
  if (xg->prev_plot_mode == MOVEPTS_MODE && xg->plot_mode != MOVEPTS_MODE)
  {
    xg->is_point_moving = False;
    XtRemoveEventHandler(xg->workspace,
      XtAllEvents, TRUE,
      (XtEventHandler) mp_button, (XtPointer) xg);
    map_move_points(False);
    xg->nearest_point = -1;
    plot_once(xg);
  }
  /* Else -- if xyplotting -- turn it on */
  else if (xg->plot_mode == MOVEPTS_MODE &&
           xg->prev_plot_mode != MOVEPTS_MODE)
  {
    if (xg->is_xyplotting || xg->is_spinning ||
      (xg->is_touring && !xg->is_pp) ||
      (xg->is_corr_touring && !xg->is_corr_pursuit))
    {
      xg->is_point_moving = True;
      map_move_points(True);
      XtAddEventHandler(xg->workspace,
        ButtonPressMask | ButtonReleaseMask, FALSE,
        (XtEventHandler) mp_button, (XtPointer) xg);
      (void) XtAppAddWorkProc(app_con, RunWorkProcs, (XtPointer) NULL);
      /* Should use a different cursor ... */

    } else {
      show_message("Moving points won\'t work in this mode; sorry\n", xg);

      /*
       * It may be necessary to disable point moving if xg->std_type != 0
      */
    }
  }
}
