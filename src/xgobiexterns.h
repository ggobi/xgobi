/* xgobiexterns.h */

#include <stdio.h>

extern Boolean set_plot1dvar(xgobidata *, int);
extern Boolean set_xyplotxvar(xgobidata *, int);
extern Boolean set_xyplotyvar(xgobidata *, int);

extern void init_missing_array(int nr, int nc, xgobidata *xg);
extern void alloc_transform_tp(xgobidata *xg);
extern void update_nrgroups_in_plot(xgobidata *xg);
extern void resort_vgroup_ids(xgobidata *, int *);
extern double randvalue(void);
extern void rnorm2(double *, double *);
extern void AllocBox(xgobidata *);
extern Widget CreateCommand(xgobidata *, char *, Boolean, Widget, Widget, Widget, char *);
extern Widget CreateToggle(xgobidata *, char *, Boolean, Widget, Widget, Widget, Boolean, int, Widget, char *);
extern float NiceValue(float);
extern XtCallbackProc PopupCaselist(Widget, xgobidata *, XtPointer);
extern XtCallbackProc PopupVarlist(Widget, xgobidata *, XtPointer);
extern Boolean RunWorkProc(xgobidata *);
extern Boolean RunWorkProcs(void *);
extern void SetNiceRange(int, xgobidata *);
extern void setToggleBitmap(Widget, Boolean);
extern int Sread_array(xgobidata *);
extern XtCallbackProc about_xgobi_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc reset_clusters_cback(Widget, xgobidata *, XtPointer);
extern void active_paint_lines(xgobidata *);
extern Boolean active_paint_points(xgobidata *);
extern void add_menupb_help(int *n, Widget w, char *fname);
extern void add_pb_help(int *n, Widget w, char *fname);
extern void add_sbar_help(int *n, Widget w, char *fname);
extern XtEventHandler add_tform_menu(Widget, xgobidata *, XEvent *);
extern void add_under_brush_label(xgobidata *xg);
extern void add_x_axis(xgobidata *, icoords *, int, tickinfo *);
extern void add_x_gridlines(xgobidata *, int, int, tickinfo *);
extern void add_y_axis(xgobidata *, icoords *, int, tickinfo *);
extern void add_y_gridlines(xgobidata *, int, int, tickinfo *);
extern void adjust_limits(float *, float *);
extern void all_tour_reproject(xgobidata *);
extern void alloc_axis_arrays(xgobidata *);
extern void alloc_bin_entropy(int);
extern void alloc_bin_fts(int);
extern void alloc_brush_arrays(xgobidata *);
extern void alloc_central_mass(int);
extern void alloc_corr(xgobidata *);
extern void alloc_hermite(int, int);
extern void alloc_holes(int);
extern void alloc_legendre(int, int);
extern void alloc_line_edit_arrays(xgobidata *);
extern void alloc_natural_hermite(int, int);
extern void alloc_pipeline_arrays(xgobidata *);
extern void alloc_plot_arrays(xgobidata *);
extern void alloc_rotate_arrays(xgobidata *);
extern void alloc_skewness(int);
extern void alloc_smooth_arrays(xgobidata *);
extern void alloc_std_vars(xgobidata *);
extern void alloc_tour(xgobidata *);
extern void alloc_transform_types(xgobidata *);
extern void announce_brush_data(xgobidata *);
extern void announce_ids(xgobidata *);
extern void announce_line_brush_data(xgobidata *);
extern void announce_rows_in_plot(xgobidata *);
extern void announce_tour_coefs(xgobidata *);
extern void assign_points_to_bins(xgobidata *);
extern void ax_rot_reproject(xgobidata *);
extern void basis_dir_ang(xgobidata *);
extern float bin_entropy_index(float **, int, int *, float);
extern void bin_fts_deriv(float **, float **, float *, float *, float **, int, int *, int, int, int *, float);
extern float bin_fts_index(float **, int, int *, float);
extern void bm_cancel(void);
extern XtCallbackProc br_complement_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc br_linkopt_menu_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc br_opt_menu_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc br_perst_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc br_save_to_group_var(Widget, xgobidata *, XtPointer);
extern void save_brush_groups(xgobidata *xg);
extern XtCallbackProc br_trans_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc br_undo_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc br_update_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc brush_active_cback(Widget, xgobidata *, XtPointer);
extern XtEventHandler brush_button(Widget, xgobidata *, XEvent *);
extern XtCallbackProc brush_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc brush_lines_cback(Widget, xgobidata *, XtPointer);
extern XtEventHandler brush_motion(Widget, xgobidata *, XEvent *, Boolean *);
extern void brush_on(xgobidata *);
extern Boolean brush_once(xgobidata *, Boolean force);
extern XtCallbackProc brush_points_cback(Widget, xgobidata *, XtPointer);
extern void build_caselist(xgobidata *);
extern void build_circle(icoords *, int, XArc *, int, short);
extern void build_glyph(xgobidata *, int, icoords *, int, XPoint *, int *, XSegment *, int *, XRectangle *, int *, XRectangle *, int *, XArc *, int *, XArc *, int *);
extern void build_labelled_menu(Widget *, Widget *, char *, Widget *, Widget *, Widget *, char **, char **, int, int, Widget, Widget, XtOrientation, XFontStruct *, char *, xgobidata *);
extern void build_plus(icoords *, int, XSegment *, int, short);
extern void build_rect(icoords *, int, XRectangle *, int, short);
extern void build_varlist(xgobidata *);
extern void build_x(icoords *, int, XSegment *, int, short);
extern float calc_norm(float *, int);
extern float calc_norm_sq(float *, int);
extern void calc_var_xy(xgobidata *, int, int, int *, int *);
extern void carry_corr_vars(xgobidata *);
extern void carry_plot1d_vars(xgobidata *);
extern void carry_spin_vars(xgobidata *);
extern void carry_tour_vars(xgobidata *);
extern void carry_xyplot_vars(xgobidata *);
extern void central_mass_deriv(float **, float **, float *, float *, float **, int, int *, int, int, int *);
extern float central_mass_index(float **, int, int *);
extern void change_first_projection(int, int, xgobidata *);
extern void change_spin_direction(void);
extern void check_cprof_fac_and_offsets(float *, float *, float *, float *, float *, float *, int);
extern void check_deselection(int, xgobidata *);
extern void check_pp_fac_and_offsets(xgobidata *, float *, float *, float *, float *, float *, float *, int);
extern int check_proximity(float **, float **, int);
extern Boolean check_singular_vc(void);
extern int check_x_axis(xgobidata *, int, tickinfo *);
extern int check_y_axis(xgobidata *, int, tickinfo *);
extern XtCallbackProc choose_tour_interp_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc choose_tour_io_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc clone_xgobi_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc scatmat_xgobi_cback(Widget, xgobidata *, XtPointer);
extern void compute_vc_matrix(xgobidata *);
extern void convert_axes(tickinfo *, xgobidata *);
extern void convert_ticks(int, int, tickinfo *, xgobidata *);
extern void copy_basis(float **, float **, int);
extern void copy_brushinfo_to_senddata(xgobidata *);
extern Boolean copy_impnames(char *, xgobidata *);
extern Boolean copy_imputations(char *, xgobidata *);
extern void copy_matrix(float **, float **, int, int);
extern void copy_raw_to_tform(xgobidata *);
extern void copy_resources(char *, xgobidata *);
extern void copy_tform_to_sphered(xgobidata *);
extern void copy_u0_to_pd0(xgobidata *);
extern void copy_u1_to_pd1(xgobidata *);
extern Boolean corr_backtracking(xgobidata *);
extern void corr_event_handlers(xgobidata *, Boolean);
extern void corr_proc(xgobidata *);
extern void corr_reproject(xgobidata *);
extern void corr_span_planes(xgobidata *);
extern void corr_tour_on(xgobidata *);
extern int corr_varselect(int, int, int, xgobidata *);
extern void cp_index(xgobidata *, int, int);
extern XtCallbackProc cprof_plot_cback(Widget, xgobidata *, XtPointer);
extern void cprof_plot_once(xgobidata *);
extern void create_default_lines(xgobidata *);
extern void destroy_varsel_widgets(int, xgobidata *);
extern void plot1d_cycle_proc(xgobidata *);
extern void plot1d_reproject(xgobidata *);
extern void plot1d_texture_var(xgobidata *);
extern Boolean plot1d_varselect(int, int, int, xgobidata *);
extern void draw_ax_spin_axes(xgobidata *);
extern void draw_brush(xgobidata *);
extern void draw_cmanip_var(xgobidata *);/* interactive gt*/
extern void draw_connecting_lines(xgobidata *);
extern void draw_corr_axes(xgobidata *);
extern void draw_current_glyph(xgobidata *);
extern void draw_diamond_around_point(int, Drawable, xgobidata *);
extern void draw_editing_lines(xgobidata *);
extern void draw_last_touched(xgobidata *);
extern void draw_manip_var(xgobidata *);/* interactive gt*/
extern void draw_ob_spin_axes(xgobidata *);
extern void draw_ob_var_lines(xgobidata *);
extern void draw_tour_axes(xgobidata *);
extern int dsvd(float **, int, int, float *, float **);
extern void entropy_deriv(float **, float **, float **, int, int *, int, float *, float *, float, int, int *);
extern float entropy_index(float **, int, int *, float);
extern XtCallbackProc exit_panel_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc exit_solo_cback(Widget, xgobidata *, XtPointer);
extern XtEventHandler expose_cback(Widget, xgobidata *, XEvent *, Boolean *);
extern void extend_axes(int, int, tickinfo *, xgobidata *);
extern void fill_extra_column(xgobidata *);
extern void find_Rmat(xgobidata *);
extern void find_extended_brush_corners(icoords *, icoords *, xgobidata *);
extern int find_gid(glyphv *);
extern void find_glyph_type_and_size(int, glyphv *);
extern int find_mono(void);
extern int find_nearest_point(icoords *, xgobidata *);
extern void find_plot_center(xgobidata *);
extern void find_quadrant(xgobidata *);
extern void find_root_name_of_data(char *, char *);
extern void find_tick_label(xgobidata *, int, int, float *, char *);
extern void fname_popup(Widget, xgobidata *);
extern void free_axis_arrays(xgobidata *);
extern void free_bin_entropy(void);
extern void free_bin_fts(void);
extern void free_brush_arrays(xgobidata *);
extern void free_caselist(xgobidata *);
extern void free_central_mass(void);
extern void free_corr(xgobidata *);
extern void free_hermite(int);
extern void free_holes(void);
extern void free_legendre(int);
extern void free_lgroups(xgobidata *);
extern void free_line_edit_arrays(xgobidata *);
extern void free_natural_hermite(int);
extern void free_pipeline_arrays(xgobidata *);
extern void free_plot_arrays(void);
extern void free_rgroups(xgobidata *);
extern void free_rotate_arrays(void);
extern void free_skewness(void);
extern void free_smooth_arrays(void);
extern void free_std_vars(xgobidata *);
extern void free_tour(xgobidata *);
extern void free_txtr_var(void);
extern void free_varlist(xgobidata *);
extern void fts_deriv(float **, float **, float **, int, int *, int, float *, float *, float, int, int *);
extern float fts_index(float **, int, int *, float);
extern void gen_norm_variates(int, int, float *);
extern void get_cprof_win_dims(float *, float *, xgobidata *);
extern void get_pp_win_dims(xgobidata *, float *, float *);
extern int glyph_color_pointtype(xgobidata *, int);
extern void gram_schmidt(float *, float *, int);
extern void grand_tour_on(xgobidata *);
extern void help(Widget, char *, xgobidata *);
extern XtCallbackProc help_cback(Widget, xgobidata *, XtPointer);
extern void hermite_deriv1(float **, float **, float *, float *, float **, int, int *, int, int *, int, int);
extern float hermite_index1(float **, int, int *, int);
extern void highlight_backtrack_cmd(void);
extern void highlight_pause_cmd(void);
extern void holes_deriv(float **, float **, float *, float *, float **, int, int *, int, int, int *);
extern float holes_index(float **, int, int *);
extern void id_proc(xgobidata *xg);
extern void identify_on(xgobidata *);
extern void init_GCs(xgobidata *);
extern void init_V(xgobidata *);
extern void init_axes(xgobidata *, Boolean);
extern void init_basis(xgobidata *);
extern void init_brush_colors(AppData *);
extern void init_brush_menus(void);
extern void init_brush_size(xgobidata *);
extern void init_brush_vars(xgobidata *);
extern void init_color_ids(xgobidata *);
extern void init_corr(xgobidata *);
extern void init_plot1d_vars(xgobidata *);
extern void init_erase(xgobidata *);
extern void init_glyph_ids(xgobidata *);
extern void init_help(xgobidata *);
extern void init_identify_vars(xgobidata *);
extern void init_line_colors(xgobidata *);
extern void init_line_edit_vars(xgobidata *);
extern void init_missing_groupvar(xgobidata *);
extern void init_ob_rotate(xgobidata *);
extern void init_options(xgobidata *);
extern void init_parcoords(xgobidata *);
extern void init_plotwindow_vars(xgobidata *, int);
extern void init_point_moving(xgobidata *);
extern void init_rotate_vars(xgobidata *);
extern void init_scale_vars(xgobidata *);
extern void init_stdview_menu(xgobidata *);
extern void init_tickdelta(xgobidata *);
extern void init_ticks(icoords *, xgobidata *);
extern void init_tour(xgobidata *, int);
extern void init_tour_interp_menu(void);
extern void init_tour_pp_GCs(xgobidata *);
extern void init_tour_pp_menu(void);
extern void init_transform_types(xgobidata *);
extern float inv_transform(int, int, xgobidata *);
extern void init_trig(xgobidata *);
extern void init_xyplot_vars(xgobidata *);
extern float inner_prod(float *, float *, int);
extern void interp_proc(xgobidata *);
extern void invert_proj_coords(xgobidata *);
extern XtCallbackProc launch_missing_cback(Widget, xgobidata *, XtPointer);
extern void legendre_deriv(float **, float **, float *, float *, float **, int, int *, int, int *, int, int);
extern float legendre_index(float **, int, int *, int);
extern void line_edit_proc(xgobidata *);
extern void line_editor_on(xgobidata *);
extern void make_arrows(xgobidata *);
extern void make_brush(xgobidata *);
extern void make_corr(xgobidata *);
extern void make_plot1d(xgobidata *);
extern void make_erase_menu(xgobidata *, Widget);
extern void make_identify(xgobidata *);
extern void make_line_editor(xgobidata *);
extern void make_move_points(xgobidata *);
extern void make_plot_window(xgobidata *);
extern void make_plotwindow_mouse_labels(xgobidata *);
extern void make_pp_panel(xgobidata *, Widget);
extern void make_rotate(xgobidata *);
extern void make_scaling(xgobidata *);
extern void make_section_panel(xgobidata *, Widget, Dimension);
extern void make_stdview(xgobidata *, Widget);
extern void make_tour(xgobidata *);
extern void make_varpanel(xgobidata *);
extern void make_widgets(xgobidata *);
extern void make_xyplot(xgobidata *);
extern void map_brush(xgobidata *, Boolean);
extern void map_group_var(xgobidata *);
extern void map_section_panel(Boolean);
extern void map_tour_panel(xgobidata *, Boolean);
extern void map_tour_pp_panel(Boolean);
extern void map_xyplot(Boolean);
extern void matrix_mult(float **, float **, float **, int, int, int);
extern float mean_fn(float *, int, int *);
extern float mean_fn2(float *, float *, int, int *);
extern void min_max(xgobidata *, float **, int *, int, float *, float *);
extern float median_largest_dist(xgobidata *, float **, int *, int, float *, float *);
extern float mean_largest_dist(xgobidata *, float **, int *, int, float *, float *);
extern void move_points_on(xgobidata *);
extern void move_points_proc(xgobidata *);
extern void natural_hermite_deriv(float **, float **, float *, float *, float **, int, int *, int, int *, int, int);
extern float natural_hermite_index(float **, int, int *, int);
extern void nback_update_label(xgobidata *);
extern void new_basis(xgobidata *);
extern void norm(float *, int);
extern void number_length(float, int *);
extern XtEventHandler ob_button(Widget, xgobidata *, XEvent *);
extern void ob_rot_reproject(xgobidata *);
extern void ob_rotate_proc(xgobidata *);
extern XtCallbackProc open_exclusion_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_tform_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_extend_xgobi_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_import_xgobi_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_imputation_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_jitter_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_infer_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_new_xgobi_popup_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc open_smooth_popup_cback(Widget, xgobidata *, XtPointer);
extern XtConvertSelectionProc pack_brush_data(Widget, Atom *, Atom *, Atom *, XtPointer *, unsigned long *, int *);
extern XtConvertSelectionProc pack_erase_data(Widget, Atom *, Atom *, Atom *, XtPointer *, unsigned long *, int *);
extern XtSelectionCallbackProc unpack_erase_data(Widget, XtPointer,
  Atom *, Atom *, XtPointer, unsigned long *, int *);
extern XtSelectionDoneProc pack_brush_done(Widget, Atom *, Atom *);
extern XtSelectionDoneProc pack_erase_done(Widget, Atom *, Atom *);
extern XtLoseSelectionProc pack_brush_lose(Widget, Atom *, XtPointer);
extern XtLoseSelectionProc pack_erase_lose(Widget, Atom *, XtPointer);
extern XtConvertSelectionProc pack_ids(Widget, Atom *, Atom *, Atom *, XtPointer *, unsigned long *, int *);
extern XtSelectionDoneProc pack_ids_done(Widget, Atom *, Atom *);
extern XtLoseSelectionProc pack_ids_lose(Widget, Atom *, XtPointer);
extern XtConvertSelectionProc pack_line_brush_data(Widget, Atom *, Atom *, Atom *, XtPointer *, unsigned long *, int *);
extern XtSelectionDoneProc pack_line_brush_done(Widget, Atom *, Atom *);
extern XtLoseSelectionProc pack_line_brush_lose(Widget, Atom *, XtPointer);
extern XtConvertSelectionProc pack_rowsinplot_data(Widget, Atom *, Atom *, Atom *, XtPointer *, unsigned long *, int *);
extern void passive_update_cprof_plot(xgobidata *);
extern XtCallbackProc pc_axes_cback(Widget, xgobidata *, XtPointer);/* interactive gt*/
extern void plane_to_screen(xgobidata *);
extern void plot_bins(xgobidata *);
extern void plot_nearest_id(xgobidata *, Drawable);
extern void plot_once(xgobidata *);
extern void possibly_draw_bandwidth(xgobidata *);
extern void pp_dir(xgobidata *);
extern void pp_index(xgobidata *, int, int);
extern float pp_index_retval(xgobidata *);
extern void princ_angs(xgobidata *);
extern XtCallbackProc princ_comp_cback(Widget, xgobidata *, XtPointer);
extern void princ_dirs(xgobidata *);
extern XtCallbackProc print_panel_cback(Widget, xgobidata *, XtPointer);
extern void quickplot_once(xgobidata *);
extern void read_array(xgobidata *);
extern void read_bm_filename(Widget, xgobidata *);
extern Boolean read_collabels(char *, Boolean, xgobidata *);
extern Boolean read_connecting_lines(char *, Boolean, xgobidata *);
extern Boolean read_erase(char *, Boolean, xgobidata *);
extern int read_extra_resources(char *);
extern void read_ids(xgobidata *);
extern void read_imputation_data(xgobidata *);
extern void read_imputation_names(xgobidata *);
extern Boolean read_jitter_values(char *, Boolean, xgobidata *);
extern Boolean read_line_colors(char *, Boolean, Boolean, xgobidata *);
extern void read_line_paint(xgobidata *);
extern Boolean read_nlinkable(char *, Boolean, xgobidata *);
extern void read_paint(xgobidata *);
extern void read_erase_sent(xgobidata *);
extern Boolean read_point_colors(char *, Boolean, Boolean, xgobidata *);
extern Boolean read_point_glyphs(char *, Boolean, Boolean, xgobidata *);
extern Boolean read_rowlabels(char *, Boolean, xgobidata *);
extern void read_rows_in_plot(xgobidata *);
extern void read_tour_coefs(xgobidata *);
extern Boolean read_vgroups(char *, Boolean, xgobidata *);
extern void realloc_lines(xgobidata *);
extern void realloc_tform(xgobidata *);
extern void recalc_vc(int, xgobidata *);
extern void refresh_all_var_labels(xgobidata *);
extern void refresh_vbox(xgobidata *, int, int);
extern void refresh_vboxes(xgobidata *);
extern void refresh_vlab(int, xgobidata *);
extern void reinit_brush_colors(xgobidata *);
extern void reinit_corr(xgobidata *);
extern void reinit_nsteps(void);
extern void reinit_spin(xgobidata *);
extern void reinit_tour(xgobidata *);
extern void reinit_tour_hist(xgobidata *);
extern void reinit_transient_brushing(xgobidata *);
extern void reset_3d_cmds(xgobidata *);
extern void reset_backtrack_cmd(Boolean, Boolean, Boolean, Boolean);
extern void reset_br_types(xgobidata *);
extern XtCallbackProc reset_brush_cback(Widget, xgobidata *, XtPointer);
extern void reset_corr_pursuit_cmd(xgobidata *, Boolean);
extern void reset_corr_pause_cmd(xgobidata *);
extern void reset_cp_plot(void);
extern void reset_cycleback_cmd(Boolean, Boolean, char *);
extern XtCallbackProc reset_glyphs_cback(Widget, xgobidata *, XtPointer);
extern void reset_interp_cmd(int);
extern void reset_interp_proc(void);
extern void reset_last_touched(xgobidata *);
extern XtCallbackProc reset_line_colors_cback(Widget, xgobidata *, XtPointer);
extern void reset_nvars_cprof_plot(xgobidata *);
extern void reset_one_label(xgobidata *, int, int);
extern void reset_spin_pause_cmd(xgobidata *);
extern void reset_tour_pause_cmd(xgobidata *);
extern XtCallbackProc reset_point_colors_cback(Widget, xgobidata *, XtPointer);
extern void reset_pp_plot(void);
extern void reset_princ_comp(Boolean, xgobidata *);
extern void reset_rows_in_plot(xgobidata *, Boolean);
extern void reset_tour_cont_fact_menu(xgobidata *);/* interactive gt*/
extern void reset_tour_link_menu(xgobidata *);
extern void reset_tourhist_cmds(xgobidata *, int);
extern void reset_var_labels(xgobidata *, int);
extern void reset_var_panel(xgobidata *);
extern XtEventHandler resize_cback(Widget, xgobidata *, XEvent *, Boolean *);
extern void retrieve_basis(xgobidata *);
extern void rock_proc(xgobidata *);
extern void rotate_on(xgobidata *);
extern int save_collabels(char *, int *, int, int, xgobidata *);
extern Boolean save_missing(char *, int *, int, int *, int, xgobidata *);
extern int save_resources(char *, xgobidata *);
extern int save_rowlabels(char *, int *, int, xgobidata *);
extern XtEventHandler scale_button_event(Widget, xgobidata *, XEvent *);
extern void scale_proc(xgobidata *);
extern void scaling_on(xgobidata *);
extern void scaling_proc(xgobidata *);
extern void set_br_linkopt_menu_marks(xgobidata *);
extern void set_br_opt_menu_marks(xgobidata *);
extern void set_bt_firsttime(void);
extern void set_counting_to_stop(int);
extern int set_deci(float);
/*extern void set_deletion(Boolean);*/
extern void set_display_menu_marks(xgobidata *);
extern void set_fac_and_offsets(float, float, float, float *, float *, float *);
extern void set_fading_var(int);
extern void set_id_linkopt_menu_marks(xgobidata *);
extern void set_local_scan_dir_in(xgobidata *);
extern void set_manip_var(xgobidata *, int);/* interactive gt*/
extern void set_tourvar(xgobidata *, int);/* interactive gt*/
extern void set_xcorrvar(xgobidata *, int);/* interactive gt*/
extern void set_ycorrvar(xgobidata *, int);/* interactive gt*/
extern void set_cmanip_var(xgobidata *, unsigned int, int);/* interactive gt*/
extern void set_cxfrozen_var(xgobidata *, int);/* interactive gt*/
extern void set_cyfrozen_var(xgobidata *, int);/* interactive gt*/
extern void set_mono(Widget w);
extern void set_ready_to_stop_now(Boolean);
extern void set_sens_go(Boolean);
extern void set_sens_interp(Boolean);
extern void set_sens_io(xgobidata *, int, int, int);
extern void set_sens_link_menu(int);
extern void set_sens_localscan(Boolean);
extern void set_sens_optimz(int);
extern void set_sens_pp_btn(xgobidata *xg, int);
extern void set_sens_princ_comp(xgobidata *xg, int);
extern void set_sens_reinit(Boolean);
extern void set_sens_speed(Boolean);
extern void set_sens_step(Boolean);
extern void set_sens_tour_update(int);
extern void set_shift_wrld0(xgobidata *);
extern void set_showarrows_option(Boolean, xgobidata *);
extern void set_showlines(Boolean);
extern void set_showlines_option(Boolean, xgobidata *);
extern XtEventHandler set_sticky(Widget, xgobidata *, XEvent *);
extern void set_title_and_icon(char *, xgobidata *);
extern XtCallbackProc set_tour_cont_fact_cback(Widget, xgobidata *, XtPointer);/* interactive gt*/
extern XtCallbackProc set_tour_link_state_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc set_tour_manip_type_cback(Widget, xgobidata *, XtPointer);/* interactive gt*/
extern void set_tour_section_eps(float, xgobidata *, int);
extern void set_varpanel_for_receive_tour(xgobidata *);
extern void set_varsel_label(xgobidata *);
extern void set_wm_protocols(Widget);
extern void shift_proc(xgobidata *);
extern void show_message(String, xgobidata *);
extern XtCallbackProc showlines_cback(Widget, XtPointer, XtPointer);
extern XtCallbackProc showpoints_cback(Widget, XtPointer, XtPointer);
extern void skewness_deriv(float **, float **, float *, float *, float **, int, int *, int, int, int *);
extern float skewness_index(float **, int, int *);
extern XtCallbackProc smooth_cback(Widget, XtPointer, XtPointer);
extern void span_planes(xgobidata *);
extern void spherize_data(xgobidata *);
extern XtCallbackProc spin_place_read_rmat_popup(Widget, xgobidata *, XtPointer);
extern XtCallbackProc spin_place_save_coefs_popup(Widget, xgobidata *, XtPointer);
extern XtCallbackProc spin_place_save_rmat_popup(Widget, xgobidata *, XtPointer);
extern void spin_proc(xgobidata *);
extern void spin_var_lines(xgobidata *);
extern int spin_varselect(int, int, int, xgobidata *) ;
extern float sq_inner_prod(float *, float *, int);
extern void start_corr_proc(xgobidata *);
extern void start_spin_proc(xgobidata *);
extern void start_tour_proc(xgobidata *);
extern void stop_spin_proc(xgobidata *);
extern void stop_tour_proc(xgobidata *);
extern void store_Rmat(xgobidata *);
extern void store_basis(xgobidata *, float **);
extern void strip_blanks(char *);
extern XtCallbackProc subset_panel_cback(Widget, xgobidata *, XtPointer);
extern void textur(float *, float *, int, int, float, int);
extern XtCallbackProc tour_backtrack_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_cycleback_cback(Widget, xgobidata *, XtPointer);
extern void tour_event_handlers(xgobidata *, Boolean);
extern XtCallbackProc tour_local_cback(Widget, xgobidata *, XtPointer);
extern Boolean tour_on(xgobidata *);
extern XtCallbackProc tour_pause_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_pp_cback(Widget, xgobidata *, XtPointer);
extern void tour_proc(xgobidata *);
extern void tour_read_hist(Widget, xgobidata *);
extern XtCallbackProc tour_reinit_cback(Widget, xgobidata *, XtPointer);
extern void tour_save_coefs(Widget, xgobidata *);
extern void tour_save_hist(Widget, xgobidata *);
extern void tour_section_calcs(xgobidata *, int);
extern XtCallbackProc tour_section_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_speed_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_step_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_step_go_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_storbas_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc tour_update_cback(Widget, xgobidata *, XtPointer);
extern void tour_var_lines(xgobidata *);
extern int tour_varselect(int, int, xgobidata *);
extern XtCallbackProc tourhist_on_cback(Widget, xgobidata *, XtPointer);
extern void turn_off_cprof_plotting(xgobidata *);
extern void turn_off_local_scan(xgobidata *);
extern void turn_off_optimz(xgobidata *);
extern void turn_off_pc(xgobidata *);
extern void turn_off_pp(xgobidata *);
extern void turn_off_section_tour(xgobidata *);
extern void turn_off_stepping(void);
extern void turn_on_showlines_option(xgobidata *);
extern void turn_on_xyplotting(xgobidata *);
extern int under_brush(xgobidata *, int);
extern void update_cprof_plot(xgobidata *xg);
extern void update_cprof_selectedvars(xgobidata *);
extern void update_lims(xgobidata *);
extern void update_list_selection(xgobidata *, int, Boolean);
extern void update_sphered(xgobidata *xg, int *, int);
extern void update_sticky_ids(xgobidata *);
extern int update_vc_active_and_do_svd(xgobidata *);
extern void update_world(xgobidata *);
extern XtEventHandler varexpose(Widget, xgobidata *, XEvent *, Boolean *);
extern void varlist_add_group_var(xgobidata *);
extern XtEventHandler varselect(Widget, xgobidata *, XEvent *, Boolean *);
extern void world_to_plane(xgobidata *);
extern int write_ascii_data(char *, int *, int, int *, int, int, int, xgobidata *);
extern int write_binary_data(char *, int *, int, int *, int, int, xgobidata *);
extern int xed_by_brush(xgobidata *, int, int);
extern void xfer_brushinfo(xgobidata *);
extern void xy_cycle_proc(xgobidata *);
extern void xy_reproject(xgobidata *);
extern Boolean xy_varselect(int, int, int, xgobidata *);
extern void zero_ab(xgobidata *xg);
extern void zero_corr_taus(void);
extern void zero_princ_angles(xgobidata *xg);
extern void zero_tau(xgobidata *xg);
extern void zero_tinc(xgobidata *);
extern void plot1d_on(xgobidata *);
extern void xyplot_on(xgobidata *);
extern void get_planar_range(xgobidata *, long *, long *, long *, long *);
extern void alloc_corr_index(xgobidata *);
extern void free_corr_index(xgobidata *);
extern void get_corr_index(xgobidata *, float *, float *, float *);
extern XtCallbackProc corr_pursuit_cback(Widget, xgobidata *, XtPointer);
extern void make_cp_plot(xgobidata *, Widget);
extern void make_cp_panel(xgobidata *, Widget);
extern void reset_corr_sphere_cmd(xgobidata *, Boolean);
extern void set_sens_corr_reinit_cmd(xgobidata *, Boolean);
extern void cp_dir(xgobidata *);
extern void stop_corr_proc(xgobidata *);
extern void init_corr_basis(xgobidata *);
extern void read_missing_values(xgobidata *);
extern void strip_suffixes(xgobidata *);
extern Boolean read_rgroups(char *, Boolean, xgobidata *);
extern void set_lgroups(Boolean, xgobidata *);
extern void jitter_data(xgobidata *);
extern void init_jitfac(xgobidata *);
extern int sample_xgobi(int, xgobidata *);
extern int plot_too_big(xgobidata *);
extern int brush_save_colors(char *, int *, int, xgobidata *);
extern int brush_save_glyphs(char *, int *, int, xgobidata *);
extern int brush_save_erase(char *, int *, int, xgobidata *);
extern int save_lines(char *, int *, int, xgobidata *);
extern int save_line_colors(char *, int *, int, xgobidata *);
extern void Clone_XGobi(void);
extern void Clone_XGobi_Scatmat(void);
extern void spin_save_coefs(Widget, xgobidata *);
extern void spin_save_rmat(Widget, xgobidata *);
extern void generate_ticks(int, int, lims *, tickinfo *, xgobidata *);
extern void spin_read_rmat(Widget, xgobidata *);
extern int derivs_equal_zero(xgobidata *);
extern void zero_indx_prev(void);
extern float mean_fn4(float *, float *, float *, float *, int, int *);
extern int find_selected_cols(xgobidata *xg, int *cols);
extern int numvargroups(xgobidata *);
extern void add_vgroups(xgobidata *, int *, int *);
extern void jitter_one_value(int, int, xgobidata *);
extern void set_sens_missing_menu_btns(Boolean);
extern void reset_tform(xgobidata *);
extern Boolean reread_dat(char *, xgobidata *);
extern void read_new_data(Widget w, xgobidata *xg);
extern int find_pgsize(int, int);
extern int make_xgobi(Boolean, char *, float **, char *, Boolean, short **, Boolean, int, int, char **, int, char **, int, connect_lines *, xgobidata *, Widget);
extern void set_sens_direction(Boolean);

extern void read_binary(FILE *, xgobidata *);
extern Boolean find_data_start(FILE *);
extern void print_cprof_win(xgobidata *, FILE *, float, float, float, float, float, XColor *, XColor *, unsigned int);
extern void print_pp_win(xgobidata *, FILE *, float, float, float, float, float, XColor *);

/* juergen */
extern void set_brush_menu_cmd(Boolean);
extern void reset_io_line_read_menu_cmd(xgobidata *);
extern void set_Edit_Lines_cmd(xgobidata *, Boolean);

#ifdef XPLORE
extern XtCallbackProc start_xplore_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc stop_xplore_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc pass_data_xplore_cback(Widget, xgobidata *, XtPointer);
extern XtCallbackProc pass_projection_xplore_cback(Widget, xgobidata *, XtPointer);
#endif

extern void gen_unif_variates(int, int, float *, float); /* jitter */

/* inna & sigbert */
extern void median_smoother(long *, long *, int, long *, int, unsigned long *,
  int *, int);
extern void mean_smoother(long *, long *, int, long *, int, unsigned long *,
  int *, int);
extern void nadaraya_watson_smoother(long *, long *, int, long *, int,
  unsigned long *, int *);
extern void spline_smoother(long *, long *, int, long *, int, unsigned long *,
  int *);
extern void smooth_data(xgobidata *);
extern void init_smooth_vars(xgobidata *);

/* frozen variables */
extern int add_variable(xgobidata *, int);
extern void remove_variable(xgobidata *, int);
extern Boolean var_frozen(xgobidata *, int);
extern void draw_frozen_var(xgobidata *);
extern void set_frozen_var(xgobidata *, int);
extern void draw_cfrozen_var(xgobidata *);
extern void set_cxfrozen_var(xgobidata *, int);
extern void set_cyfrozen_var(xgobidata *, int);
extern void add_xvariable(xgobidata *, int);
extern void remove_xvariable(xgobidata *, int, int);
extern void copy_vector(float *, float *, int);

