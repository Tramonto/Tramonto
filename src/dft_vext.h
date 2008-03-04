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
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern char Vext_file2[40];
#define READ_VEXT_SUMTWO     2
int ijk_to_node(int *ijk);
int ijk_to_node(int *ijk);
int round_to_int(double x);
int round_to_int(double x);
extern char Vext_file[40];
extern double **Vext_static;
#define READ_VEXT_STATIC     3
extern int Restart_Vext;
void read_external_field_n();
extern int *Comm_offset_node;
extern int *Comm_node_proc;
void comm_loc_to_glob_vec(int *n_loc,int *in_loc_vec,int *out_glob_vec);
int el_box_to_el(int iel_box);
int el_box_to_el(int iel_box);
double integrate_potential(double param1,double param2,double param3,double param4,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
int element_to_node(int ielement);
int element_to_node(int ielement);
void set_gauss_quad(int ngp,double *gp,double *gw);
extern double Vol_el;
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define NWALL_MAX_TYPE 50 
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define point_surface                   4
#define atomic_centers                  3
extern int Surface_type[NWALL_MAX_TYPE];
void setup_vext_HS_atomic(int iwall);
double pairPot_deriv_switch(double r,double x,double param1,double param2,double param3,double param4,int typePairPot);
double pairPot_switch(double r,double param1,double param2,double param3,double param4,int typePairPot);
#define WALL_FLUID  1
extern int Type_vext3D;
double pairPotparams_switch(int typePairPot,int context,int i,int j,double *param1,double *param2,double *param3,double *param4);
extern int *B2G_node;
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern double **X_wall2;
extern double **X_wall;
double Vext_1D_dash(double x,int icomp,int iwall_type);
void find_images(int idim,double cut,int *image,double **image_pos,double *node_image,double *node_ref);
double Vext_1D(double x,int icomp,int iwall_type);
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double WallParam[NWALL_MAX_TYPE];
void find_images_1D(int idim,double cut,int *image,double **image_pos,double *node_image,double *node_ref);
#define REFLECT              2
#define PERIODIC             1
extern int Type_bc[NDIM_MAX][2];
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Orientation[NWALL_MAX_TYPE];
extern double Size_x[NDIM_MAX];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double **Vext_membrane;
int ijk_box_to_node_box(int *ijk_box);
int ijk_box_to_node_box(int *ijk_box);
extern int Nodes_x[NDIM_MAX];
extern int Min_IJK_box[3];
extern int Nnodes_per_el_V;
void node_box_to_ijk_box(int node_box,int *ijk_box);
void node_box_to_ijk_box(int node_box,int *ijk_box);
int element_box_to_node_box(int iel_box);
int element_box_to_node_box(int iel_box);
extern int **Lsemiperm;
void node_to_position(int inode,double *NodePos);
void node_to_position(int inode,double *NodePos);
extern int *L2G_node;
extern double VEXT_MAX;
extern int *L2B_node;
extern int Nnodes;
int gsum_int(int c);
extern int *B2L_node;
extern int Nnodes_box;
double gmin_double(double c);
double gmax_double(double c);
void correct_zeroTF_array();
extern int Num_Proc;
extern int **Zero_density_TF;
int node_to_node_box_no_bound(int inode);
int node_to_node_box_no_bound(int inode);
void node_to_ijk(int node,int *ijk);
void node_to_ijk(int node,int *ijk);
void comm_vext_max(int *nnodes_vext_max,int **nodes_vext_max);
void setup_vext_LJ_atomic(int iwall);
#define VEXT_ATOMIC        6  /* 3D potential for 3D problems */
void safe_free(void **ptr);
void safe_free(void **ptr);
void setup_integrated_LJ_walls(int iwall,int *nelems_w_per_w,int **elems_w_per_w);
void comm_wall_els(int iwall,int **nelems_w_per_w,int ***elems_w_per_w,int *nelems_w_per_w_global,int **elems_w_per_w_global);
extern int Nlists_HW;
#define VEXT_3D_INTEGRATED      5  /* more proper 3D integration potential for funny geometries */
#define VEXT_1D_ORIENTATION     4  /* crude 1D-like treatment of funny geometries */
void setup_1Dvext_xmin(int iwall);
#define VEXT_1D_XMIN     3  /* crude 1D-like treatment of funny geometries */
void setup_1Dvext(int iwall);
#define VEXT_1D          2  /* 1D potential for 1D systems */
void setup_zero();
void setup_vext_max();
void setup_semiperm(int **nelems_w_per_w,int ***elems_w_per_w);
#define UNINIT_VEC -200.0
extern int Ndim;
extern double ***Vext_dash;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define VEXT_HARD        1
#define VEXT_NONE          0
extern int WallType[NWALL_MAX];
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern int Nwall;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Lvext_dash;
extern double **Vext_set;
extern int Ncomp;
extern int Nnodes_per_proc;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern double **Vext;
#define VERBOSE      3 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void setup_external_field_n(int **nelems_w_per_w,int ***elems_w_per_w);
