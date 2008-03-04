/* This file was automatically generated.  Do not edit! */
int node_to_node_box(int inode);
int node_to_node_box(int inode);
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
extern int *B2L_node;
#define NWALL_MAX_TYPE 50 
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define NCOMP_MAX 5
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double WallParam_4[NWALL_MAX_TYPE];
extern double Dielec_pore;
extern double Dielec_X;
#define DIELEC_WF_PORE     2
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int Type_poly;
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
extern int Lhard_surf;
#define MAX_ROUGH_BLOCK 100
extern double Rough_precalc[NWALL_MAX_TYPE][MAX_ROUGH_BLOCK][MAX_ROUGH_BLOCK];
extern double WallParam[NWALL_MAX_TYPE];
extern double Rough_length[NWALL_MAX_TYPE];
extern int Lrough_surf[NWALL_MAX_TYPE];
#define PI    M_PI
extern double WallParam_3[NWALL_MAX_TYPE];
extern double WallParam_2[NWALL_MAX_TYPE];
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
void node_to_position(int inode,double *NodePos);
void node_to_position(int inode,double *NodePos);
int element_to_node(int ielement);
int element_to_node(int ielement);
int el_box_to_el(int iel_box);
int el_box_to_el(int iel_box);
extern int **List_wall_node;
extern int **Wall_touch_node;
extern int *Nwall_touch_node;
extern int *Nodes_wall_box;
extern int *Index_wall_nodes;
int ijk_box_to_node_box(int *ijk_box);
int ijk_box_to_node_box(int *ijk_box);
extern int Pflag[3];
void node_box_to_ijk_box(int node_box,int *ijk_box);
void node_box_to_ijk_box(int node_box,int *ijk_box);
int element_box_to_node_box(int iel_box);
int element_box_to_node_box(int iel_box);
extern int Nodes_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
void node_to_ijk(int node,int *ijk);
extern double *Dielec_wall;
extern double *Dielec;
#define DIELEC_CONST       0
extern int Type_dielec;
#define COULOMB      1
extern int Ipot_ff_c;
#define WALL_EL        0 
extern int ***Wall_owners;
extern int **Nwall_owners;
void flag_wall_el(int inode,int ilist,int iwall,int iel_box,int **L_wall,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type);
extern double Size_x[NDIM_MAX];
#define REFLECT              2
#define PERIODIC             1
extern int Type_bc[NDIM_MAX][2];
void safe_free(void **ptr);
void safe_free(void **ptr);
extern double Rmax_zone[5];
extern int Nnodes_wall_box;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define VERBOSE      3 
extern int Iwrite;
void els_cone_pore_3D(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
void els_cone_pore_2D(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
#define tapered_pore                    9
void els_cyl_pore_3D(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
void els_slit_pore_2D(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
#define cyl3D_slit2D_pore               8
void els_cyl_pores(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
#define cyl2D_sphere3D_pore             7
void els_atomic_centers(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
#define atomic_centers                  3
void els_cyls_cos_3D(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
void els_cyls_3D(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
void els_spheres(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
#define point_surface                   4
#define colloids_cyl_sphere             2
void els_finite_planar(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
#define finite_planar_wall              1
void els_planar(int iwall,int real_wall,int itype,int **L_wall,double *x_min,int **nelems_w_per_w,int ***elems_w_per_w,int ***el_type,double **image_pos);
void find_wall_images(int idim,int *image,double **image_pos,double *pos);
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
#define cyl_periodic_3D                 6
#define finite_cyl_3D                   5
extern int Link[NWALL_MAX];
extern int **Xtest_reflect_TF;
extern int Orientation[NWALL_MAX_TYPE];
#define smooth_planar_wall              0
extern int Surface_type[NWALL_MAX_TYPE];
extern int WallType[NWALL_MAX];
#define FLUID_EL       1 
extern int **Wall_elems;
extern int Nzone;
extern int Coarser_jac;
extern double **X_wall2;
extern int Nnodes_per_proc;
extern double **X_wall;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define VEXT_1D_XMIN     3  /* crude 1D-like treatment of funny geometries */
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern int Nwall_type;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Ndim;
extern int Nwall;
extern int ****Touch_domain_boundary;
extern int Nelements_box;
extern int Nlists_HW;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
void setup_surface(FILE *fp2,int *nelems_f,int **nelems_w_per_w,int **elems_f,int ***elems_w_per_w,int *elem_zones,int ***el_type);
