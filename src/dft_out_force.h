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
#define NDIM_MAX  3
extern int Type_bc[NDIM_MAX][2];
int el_to_el_box(int iel);
extern int **Wall_elems;
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int Nnodes_per_el_V;
double calc_local_pressure(double **x,int iden_first,int inode_box);
#define IDEAL_GAS    0
extern int Ipot_ff_n;
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
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern double Betap;
void print_to_file(FILE *fp,double val,char *var_label,int first);
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
#define NWALL_MAX_TYPE 50 
extern int Orientation[NWALL_MAX_TYPE];
extern int Nlists_HW;
extern double **S_area_tot;
extern int Lper_area;
#define PI    3.141592653589793238462643383279502884197169399375
extern double Temp_elec;
double gsum_double(double c);
#define NO_SCREEN    4 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Link[NWALL_MAX];
extern int Nlink;
void sum_rho_wall(double **x,double **Sum_rho);
extern int Lhard_surf;
void integrate_rho_vdash(double **x,double **rho_vdash);
extern int Lvext_dash;
extern int Ipot_wf_c;
extern int Ndim;
extern int Nwall;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void calc_force(FILE *fp,double **x,double fac_area);
