/* This file was automatically generated.  Do not edit! */
void safe_free(void **ptr);
void safe_free(void **ptr);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define ATOMIC_CHARGE     3
#define NWALL_MAX_TYPE 50 
extern int Type_bc_elec[NWALL_MAX_TYPE];
double integrate_potential(int flag,double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
double integrate_potential(int flag,double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
void find_images(int idim,double cut,int *image,double **image_pos,double *node_image,double *node_ref);
extern double VEXT_MAX;
void find_images2(int idim,double cut,int *image,double **image_pos,double *node_image,int iwall,int iside);
#define NDIM_MAX  3
extern int Nodes_x[NDIM_MAX];
#define IN_WALL     -1
void node_to_ijk(int node,int *ijk);
int node_to_node_box(int inode);
extern double Esize_x[NDIM_MAX];
void node_to_position(int inode,double *NodePos);
int element_to_node(int ielement);
#define VERBOSE      3 
extern int Iwrite;
void set_gauss_quad(int ngp,double *gp,double *gw);
extern int Nlists_HW;
extern double Vol_el;
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int Nwall_type;
void setup_ww_integrated(int iwall,int *nelems_w_per_w,int **elems_w_per_w,int L_LJ,int L_COULOMB);
#define NWALL_MAX 600 
extern double Elec_param_w[NWALL_MAX];
#define REFLECT      2
extern int Type_bc[NDIM_MAX][2];
extern double Size_x[NDIM_MAX];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
double pairPot_switch(double r,double param1,double param2,double param3);
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Ndim;
void setup_coulomb_atomic(int iwall,int jwall);
extern int Type_coul;
void setup_lj_atomic(int iwall,int jwall);
#define ATOM_CENTERS_WW    1 
#define NO_WW              0
extern int Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int WallType[NWALL_MAX];
extern int Link[NWALL_MAX];
extern int Nlink;
extern double **Uww_link;
extern int Nwall;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern double **Uww;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void setup_wall_wall_potentials(int **nelems_w_per_w,int ***elems_w_per_w);
