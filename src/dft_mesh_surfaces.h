/* This file was automatically generated.  Do not edit! */
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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern int *B2L_node;
#define NWALL_MAX_TYPE 10 
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define NCOMP_MAX 5
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Dielec_pore;
extern double Dielec_X;
#define DIELEC_WF_PORE     2
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int Lhard_surf;
extern int Type_poly;
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
#define MAX_ROUGH_BLOCK 100
extern double Rough_precalc[NWALL_MAX_TYPE][MAX_ROUGH_BLOCK][MAX_ROUGH_BLOCK];
extern double Rough_length[NWALL_MAX_TYPE];
#define PI    3.141592653589793238462643383279502884197169399375
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
void node_to_position(int inode,double *NodePos);
int element_to_node(int ielement);
int el_box_to_el(int iel_box);
extern int **List_wall_node;
extern int **Wall_touch_node;
extern int *Nwall_touch_node;
extern int *Nodes_wall_box;
extern int *Index_wall_nodes;
int ijk_box_to_node_box(int *ijk_box);
extern int Pflag[3];
void node_box_to_ijk_box(int node_box,int *ijk_box);
int element_box_to_node_box(int iel_box);
extern int Nodes_x[NDIM_MAX];
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
void find_wall_images(int idim,int *image,double **image_pos,double *pos);
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
#define finite_cyl_3D                   5
extern int Link[NWALL_MAX];
extern int **Xtest_reflect_TF;
#define smooth_planar_wall              0
typedef struct SurfaceGeom_Struct SurfaceGeom_Struct;
extern struct SurfaceGeom_Struct *SGeom;
extern int WallType[NWALL_MAX];
#define FLUID_EL       1 
extern int **Wall_elems;
extern int Nzone;
extern int Coarser_jac;
extern int Nnodes_per_proc;
extern double **X_wall;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
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
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int Lwedge_cutout[NWALL_MAX];
extern int Lperiodic_overlay[NWALL_MAX_TYPE];
extern int Nperiodic_overlay[NWALL_MAX_TYPE];
#define NPERIODIC_MAX 4
extern int Llinear_overlay[NWALL_MAX_TYPE];
extern int Nlinear_overlay[NWALL_MAX_TYPE];
struct SurfaceGeom_Struct {
  int       surfaceTypeID;    /* ID of the type of surface */
  int       orientation;  /* orientation of the surface */
  double    *halfwidth;   /* planar surface params given in halfwidth */
  double    radius;       /* radius of spherical or cylindrical surface */
  double    halflength;   /* length of finite cylinders and pores */
  double    radius2;      /* a second radius for tapered pores or cylinders */
  int       Lwedge_cutout;    /* TRUE or FALSE for wedge cutout from basic surfac */
  double    angle_wedge_start;    /* angle as measured from x0 axis */
  double    angle_wedge_end;    /* angle as measured from x0 axis */
  int       Lrough_surface;    /* TRUE or FALSE for surface roughness */
  double    roughness;          /* maximum roughness amplitude */
  double    roughness_length;    /* lengthscale for the rougness */
  int       Lperiodic_overlay;    /* TRUE or FALSE for periodic function added to surface */
  int       Nperiodic_overlay;     /* The number of periodic functions to apply */
  double    orientation_periodic[NPERIODIC_MAX];    /* maximum amplitude for a cosine wave superimposed on a surface */
  double    amplitude[NPERIODIC_MAX];    /* maximum amplitude for a cosine wave superimposed on a surface */
  double    wavelength[NPERIODIC_MAX];    /* desired wavelength of cosine wave superimposed on a surface */
  double    origin_PeriodicFunc[NPERIODIC_MAX];     /* The origin of periodic functions to apply */
  int       Llinear_overlay;    /* TRUE or FALSE for linear function added to surface */
  int       Nlinear_overlay;     /* The number of linear functions to apply */
  double    orientation_linear[NPERIODIC_MAX];    /* maximum amplitude for a linear function superimposed on a surface */
  double    slope[NPERIODIC_MAX];    /* maximum amplitude for a linear function superimposed on a surface */
  double    origin_LinearFunc[NPERIODIC_MAX];     /* The origin of linear functions to apply */
  double    endpoint_LinearFunc[NPERIODIC_MAX];     /* The end point of linear functions to apply */
  int    *ReflectionsAreIndependent;  /* TRUE or FALSE for treating special boundary conditions */
};
void setup_surface(FILE *fp2,int *nelems_f,int **nelems_w_per_w,int **elems_f,int ***elems_w_per_w,int *elem_zones,int ***el_type);
void setup_surface(FILE *fp2,int *nelems_f,int **nelems_w_per_w,int **elems_f,int ***elems_w_per_w,int *elem_zones,int ***el_type);
