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
#define NCOMP_MAX 6
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
void print_to_screen(double val,char *var_label);
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define SCREEN_VERBOSE     3 
#define SCREEN_ERRORS_ONLY  0 
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define NDIM_MAX  3
extern double Size_x[NDIM_MAX];
extern int Ndim;
extern int Nwall;
extern int Grafted_Logical;
extern int Grafted_SegID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
extern double *Gsum_graft;
extern int Ncomp;
extern int Nmer[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
#define DENSITY        0
#define NEQ_TYPE       12 
#define WJDC_FIELD     8
#define G_CHAIN       11 
extern int Phys2Unk_first[NEQ_TYPE];
extern int **Nel_hit;
extern int **Nel_hit2;
double graft_freen(int npol,int iunk,int icomp,double **x);
double WJDCgraft_freen_bulk(double **x);
extern int Npol_comp;
extern int Lseg_densities;
#define NMER_MAX     200
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int **Bonds_SegAll;
extern int *Nbonds_SegAll;
extern int **Poly_to_Unk_SegAll;
extern int SegAll_to_Poly[NMER_MAX];
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern double Rho_g[NCOMP_MAX];
double integrateInSpace(double(*fp_integrand)(int,int,double **),int iunk,int **nelhit,double **x,double *profile);
double integrand_WJDCcomp_freen_bulk(int iunk,int inode_box,double **x);
double integrand_WJDCcomp_freen(int iunk,int inode_box,double **x);
double integrand_WJDC_freen_bulk(int iunk,int inode_box,double **x);
double integrand_WJDC_freen(int iunk,int inode_box,double **x);
extern int     Nnodes_per_proc;
extern int     *L2B_node;
extern double  Vol_el;
extern int     Nnodes_per_el_V;
double gsum_double(double c);
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Nseg_type_pol;
extern double Betamu_chain[NMER_MAX];
extern int Grafted_TypeID[NCOMP_MAX];
#define NWALL_MAX 600 
extern double *Total_area_graft;
#define REFLECT              2
extern int Type_bc[NDIM_MAX][2];
extern int **NodesS_GID_global;
void node_to_position(int inode,double *NodePos);
extern double Esize_x[NDIM_MAX];
extern double Area_surf_el[3];
extern int WallType[NWALL_MAX];
extern int ***Surf_elem_to_wall_global;
extern int ***Surf_normal_global;
extern int **NelemsS_global;
extern int **S2B_node;
extern int *NodesS_global;
extern int Nlists_HW;
extern double **S_area_tot;
#define GRAFT_DENSITY 1




