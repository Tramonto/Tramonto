/* This file was automatically generated.  Do not edit! */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#include <string.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define LOCAL_N 1
int loc_find(int iunk,int inode,int flag);
int position_to_node(double *NodePos);
extern int Max_IJK[3];
extern int Nnodes;
extern int Min_IJK[3];
void safe_free(void **ptr);
void safe_free(void **ptr);
#define LB_LINEAR    0
extern int Load_Bal_Flag;
void node_to_position(int inode,double *NodePos);
extern int Nunk_per_node;
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
extern int Nnodes_per_proc;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Num_Proc;
double gsum_double(double c);
int gmin_int(int c);
int gmax_int(int c);
void load_balance(int flag,double *fill_time,int *N_update,int **update);
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
#define FLAG_BULK   -888
#define FLAG_1DBC   -999
extern int *Mesh_coarsen_flag;
extern int L1D_bc;
extern int Mesh_coarsening;
extern int **Zero_density_TF;
#define IN_WALL             -1
#define PERIODIC             1
#define NDIM_MAX  3
extern int Nodes_x[NDIM_MAX];
extern double Esize_x[NDIM_MAX];
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
extern int Ndim;
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_attr;
extern int Ncomp;
int node_to_node_box(int inode);
void node_to_ijk(int node,int *ijk);
double set_weight_for_node(int inode);
void sort_int_array(int n,int ra[]);
#define SCREEN_VERBOSE     3 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
void rcb_median_merge(void *in,void *inout,int *len,MPI_Datatype *dptr);
void rcb_box_merge(void *in,void *inout,int *len,MPI_Datatype *dptr);
typedef struct rcb_dot rcb_dot;
typedef struct rcb_box rcb_box;
void rcb_stats(double timetotal,struct rcb_dot *dotpt,int dotnum,double *timers,int *counters,struct rcb_box *rcbbox,int reuse);
void rcb_check(struct rcb_dot *dotpt,int dotnum,int dotorig,struct rcb_box *rcbbox);
void rcb_error(int size);
typedef struct rcb_median rcb_median;
struct rcb_median {          /* RCB cut info */
  double    lototal, hitotal;   /* weight in each half of active partition */
  double    valuelo, valuehi;	/* position of dot(s) nearest to cut */
  double    wtlo, wthi;         /* total weight of dot(s) at that position */
  int       countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;	/* unique proc who owns a nearest dot */
};
typedef struct rcb_tree rcb_tree;
struct rcb_tree {	     /* tree of RCB cuts */
  double    cut;        	/* position of cut */
  int       dim;	        /* dimension (012) of cut */
};
struct rcb_box {       	     /* bounding box */
  double    lo[3], hi[3];	/* xyz lo/hi bounds */
};
struct rcb_dot {	        /* dot = point in 3-space */
  double    x[3];		/* location of dot */
  double    weight;             /* weight of dot - if used must be > 0 */
};
void rcb(struct rcb_dot **dots,int *pdotnum,int *pdotmax,int wtflag,double overalloc,int reuse,int *pdottop,struct rcb_box *rcbbox,struct rcb_tree *treept);
extern int RCB_STATS;
extern int RCB_CHECK;
