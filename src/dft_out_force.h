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
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern double Vol_el;
extern double ***Vext_dash;
#define LAST_NODE_RESTART    4
#define LAST_NODE            3
#define IN_BULK              0
#define REFLECT              2
int el_to_el_box(int iel);
extern int **Wall_elems;
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int Nnodes_per_el_V;
double calc_local_pressure(double **x,int iden_first,int inode_box);
#define IDEAL_GAS    0
extern int Ipot_ff_n;
#define NDIM_MAX  3
extern int Nodes_x[NDIM_MAX];
void node_to_ijk(int node,int *ijk);
double sum_rho_midplane(double **x);
int node_to_node_box(int inode);
void find_pot_derivs(double **x,double *psi_deriv);
#define DOWN_FRONT   7
#define DOWN_BACK    3
#define UP_FRONT     4
#define UP_BACK      0
#define LEFT_DOWN   11 
#define LEFT_UP      9
#define LEFT_FRONT   6
#define LEFT_BACK    2
#define RIGHT_DOWN  10 
#define RIGHT_UP     8
#define RIGHT_FRONT  5
#define RIGHT_BACK   1
extern double Esize_x[NDIM_MAX];
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
void find_offset(int el_type,int jdim,int *offset);
double calc_deriv(int idim,int inode0,int flag,int *blocked,double **x,int ilist);
extern int **Surf_elem_type;
#define POISSON        1
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern int *L2G_node;
void force_elec(double **x,double **Sum_dphi_dx);
extern int Nnodes_per_el_S;
extern double Area_surf_el[3];
extern int ***Surf_normal;
extern int ***Surf_elem_to_wall;
extern int **Nelems_S;
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
void node_to_position(int inode,double *NodePos);
extern int Graft_wall[NCOMP_MAX];
extern int Grafted_TypeID[NCOMP_MAX];
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
extern int **Nodes_2_boundary_wall;
extern int *L2B_node;
extern int Nnodes_per_proc;
#define NMER_MAX     200
extern int Unk2Comp[NMER_MAX];
extern int Ncomp;
extern int Nseg_tot;
extern int Lseg_densities;
void safe_free(void **ptr);
void safe_free(void **ptr);
extern double Betap;
void print_to_file(FILE *fp,double val,char *var_label,int first);
#define SCREEN_BASIC       1
#define SCREEN_VERBOSE     3 
extern int Nlists_HW;
extern double **S_area_tot;
extern int Lper_area;
#define PI    3.141592653589793238462643383279502884197169399375
extern double Temp_elec;
double gsum_double(double c);
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern double Size_x[NDIM_MAX];
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
#define IN_WALL             -1
extern int Type_bc[NDIM_MAX][2];
#define smooth_planar_wall              0
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Link[NWALL_MAX];
typedef struct SurfaceGeom_Struct SurfaceGeom_Struct;
extern struct SurfaceGeom_Struct *SGeom;
extern int WallType[NWALL_MAX];
extern int Nlink;
void sum_rho_wall(double **x,double **Sum_rho);
extern int Grafted_Logical;
#define WJDC3        5 
extern int Type_poly;
extern int Lhard_surf;
void integrate_rho_vdash(double **x,double **rho_vdash);
extern int Lvext_dash;
extern int Ipot_wf_c;
extern int Ndim;
extern int Nwall;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int Lwedge_cutout[NWALL_MAX];
#define NWALL_MAX_TYPE 20 
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
#define FILENAME_LENGTH 4096
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void calc_force(FILE *fp,double **x,double fac_area);
