#ifndef XGVISINTERN
#define XGVIS_ extern 
#define INIT(x)
#else
#define XGVIS_
#define INIT(x) x
#endif

XGVIS_ xgobidata xgobi;

typedef struct {
    XFontStruct *Font;
    int mdsDimension;
} PanelData, *PanelDataPtr;
XGVIS_ PanelData panel_data;

XGVIS_ Widget shell;
XGVIS_ Widget form0;
XGVIS_ Widget dims_left, dims_right;

/* Defines. */
#define METRIC    0
#define NONMETRIC 1

#define USER_SUPPLIED 0
#define LINK 1
#define ADJACENCY 2
#define EUCLIDIAN 3
#define MANHATTAN 4
#define MAHALANOBIS 5
#define DOTPROD 6
#define COSDIST 7
#define CLOSE 8
#define NDISTTYPES (CLOSE+1)
XGVIS_ Widget apply_dist_cmd;

#define MAXDIMS 13

#define EPSILON .00001

/* Global variables. */ 
XGVIS_ int xgv_is_running INIT(= 0);
XGVIS_ int max_dims INIT(= 10);
XGVIS_ int is_rescale  INIT(= 0);
XGVIS_ Widget dist_cmd, dist_popup, dist_mgr, dist_types[NDISTTYPES] ;
XGVIS_ int dist_type INIT(= 0) ;

XGVIS_ enum {within, between, anchored, none} mds_group_ind;
XGVIS_ double mds_power  INIT(= 1.0);
XGVIS_ double mds_weightpow INIT(= 0.0);
XGVIS_ double mds_lnorm  INIT(= 2.0);
XGVIS_ int mds_n_iters  INIT(= 1);
XGVIS_ double mds_stepsize  INIT(= 0.1);
XGVIS_ double mds_threshold_high  INIT(= 0.0);
XGVIS_ double mds_threshold_low  INIT(= 0.0);  /* what initial value? */
XGVIS_ int mds_dims INIT(= 3);

/* Used in scaling during each mds loop; set in reset_data */
XGVIS_ double distance_factor;
XGVIS_ double *distance_vector, *distance_vector_sort;
XGVIS_ double *config_distances;
XGVIS_ int ndistances;

XGVIS_ double configuration_factor;

/* Used to hold matrix structures. */
struct array {
  double **data;
  int nrows;
  int ncols;
};

/* Global data structures. */
XGVIS_ struct array dist_orig;
XGVIS_ struct array dist;
XGVIS_ struct array edges_orig;
XGVIS_ struct array edges;
XGVIS_ struct array lines;
XGVIS_ struct array *linesptr; /* points to lines, if present; else to edges */
XGVIS_ struct array pos_orig;
XGVIS_ struct array pos;
XGVIS_ char **rowlab INIT(= NULL);	/* Row labels. */
XGVIS_ char *xgv_basename;

/* for diagnostic xgobi, for distance residuals */
XGVIS_ struct array diagnostics;
XGVIS_ char **drowlab INIT(= NULL);

XGVIS_ char pcolorname[128];
XGVIS_ char lcolorname[128];
XGVIS_ char glyphname[128];

XGVIS_ int scaling_method INIT(= METRIC);
