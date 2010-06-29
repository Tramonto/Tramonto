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
#define LAST_NODE_RESTART    4
#define LAST_NODE            3
#define IN_BULK              0
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
extern int **Lsemiperm;
extern int **Wall_elems;
int el_to_el_box(int iel);
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
extern int Imax;
extern int List[2];
extern int Lhard_surf;
#define NMER_MAX     200
extern int Unk2Comp[NMER_MAX];
extern int **Nel_hit2;
extern int Nnodes_box;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int **Nel_hit;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void setup_integrals();
extern int Link[NWALL_MAX];
extern double **S_area_tot;
extern int Nlink;
extern int Nwall;
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Lcount_reflect;
#define REFLECT              2
void setup_domain_multipliers();
#define PERIODIC             1
extern int Type_bc[NDIM_MAX][2];
extern int Nodes_x[NDIM_MAX];
extern int Ndim;
extern int *L2G_node;
void node_to_ijk(int node,int *ijk);
extern int Nlists_HW;
extern int **Nodes_2_boundary_wall;
double integrateOverSurface(double(*fp_integrand)(int,int,int,double **),int iunk,double **x,double *profile);
extern int Nseg_tot;
extern int Lseg_densities;
extern int Ncomp;
double integrateInSpace_SumInComp(double(*fp_integrand)(int,int,double **),int **nelhit,double **x,double *profile);
extern double Fac_area;
extern double Fac_vol;
extern double Area;
extern int Lper_area;
double gsum_double(double c);
extern int Nnodes_per_el_V;
extern double Vol_el;
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int *L2B_node;
extern int Nnodes_per_proc;
double integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,int **nelhit,double **x,double *profile);
