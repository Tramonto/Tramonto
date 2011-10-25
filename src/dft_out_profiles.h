/* This file was automatically generated.  Do not edit! */
void print_vext(double **vext,char *output_file);
int element_to_node(int ielement);
int node_box_to_elem_box_reflect(int inode_box,int local_node,int *reflect_flag);
int node_to_elem(int inode_all,int local_node,int *reflect_flag);
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
#define PERIODIC             1
#define NDIM_MAX  3
extern int Type_bc[NDIM_MAX][2];
extern int Nelements;
void print_charge_vol(double *charge_els,char *output_file);
void print_freen_profile_1D(double *freen_profile_1D,char *output_file);
void print_charge_surf(double **charge_w_sum,char *output_file);
void print_Nodes_to_zone(int *node_to_zone,char *output_file);
void print_zeroTF(int **zero_TF,char *output_file);
extern double VEXT_MAX;
extern int Lprint_gofr;
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern double Size_x[NDIM_MAX];
extern double **Charge_x;
extern int Nlocal_charge;
extern int Nwall;
extern int L_HSperturbation;
void print_gofr(char *output_file6,double *xold);
extern int Nodes_x[NDIM_MAX];
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define NCOMP_MAX 5
#define NMER_MAX     200
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern double Betamu_chain[NMER_MAX];
#define DIFFUSIVE_INTERFACE 1
extern int Type_interface;
#define NBLOCK_MAX   20 
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
extern int ***Poly_to_Unk;
extern int Geqn_start[NCOMP_MAX];
extern int **Nbond;
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int Nmer[NCOMP_MAX];
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
#define G_CHAIN       11 
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Nseg_type_pol;
extern int Npol_comp;
extern double Temp;
extern double Mass[NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int LDeBroglie;
extern int Unk2Comp[NMER_MAX];
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Unk2Phys[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
extern double Esize_x[NDIM_MAX];
extern int Ndim;
void node_to_ijk(int node,int *ijk);
extern double Charge_f[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Temp_elec;
#define PI    3.141592653589793238462643383279502884197169399375
#define COULOMB      1
extern int Ipot_ff_c;
extern int Npoisson;
#define SCF_FIELD	  10
#define SCF_CONSTR	   9
#define BONDWTC        5
#define CAVWTC         4
#define HSRHOBAR       2
#define WJDC_FIELD     8
#define CMS_FIELD      7
#define DIFFUSION      6
#define POISSON        1
#define MF_EQ          3
extern int Phys2Nunk[NEQ_TYPE];
extern int Lseg_densities;
#define DENSITY        0
#define EXTENDED     2
#define WTC          2
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
#define CMS_SCFT     1
#define CMS          0
extern int Type_poly;
#define VERBOSE      3 
extern int Iwrite;
void print_profile(char *output_file4,double *xold);
extern double *X_old;
void print_profile_box(double **x,char *outfile);
extern double *Vext_old;
extern int Num_Proc;
extern double **Vext;
extern int Ncomp;
void collect_vext_old();
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int *Comm_offset_unk;
extern int *Comm_unk_proc;
extern int *Comm_offset_node;
extern int *Comm_node_proc;
extern int *L2G_node;
extern int Nnodes;
extern int Nunknowns;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int *L2B_node;
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int Nunk_per_node;
extern int Nnodes_per_proc;
void collect_x_old(double **x,double *xold);
