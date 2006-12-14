/* This file was automatically generated.  Do not edit! */
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
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
double gammln(double xx);
void gser(double *gamser,double a,double x,double *gln);
void gcf(double *gammcf,double a,double x,double *gln);
double gauss(double r,int i,int j);
extern double Gauss_k;
extern double Gauss_a;
double int_cr(double r_low,double r_upp,double slope_dr,int icomp,int jcomp,int irmin,double zsq,double *rx_low);
extern double ***Rism_cr;
extern int Last_nz_cr;
extern double Deltar_cr;
double deltaC_MSA(double r,int i,int j);
double uLJatt_n(double r,int i,int j);
double deltaC_MSA_int(double r,int i,int j);
double uLJatt_n_int(double r,int i,int j);
double uLJatt_n_noshift(double r,int i,int j);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
#define PI    M_PI
#define NCOMP_MAX 5
extern double Cr_rad[NCOMP_MAX][NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
typedef struct Stencil_Struct Stencil_Struct;
void sort_stencil(struct Stencil_Struct *sten);
#define NDIM_MAX  3
void sort2_int_double_array(int n,int ra[],double rb[],int *rc[],int ndim);
#define PERIODIC     1
extern int Type_bc[NDIM_MAX][2];
extern double Size_x[NDIM_MAX];
int round_to_int(double x);
int ijk_to_isten_index(int *ijk,int *el_in_radius);
extern double Jac_threshold;
void print_out_stencil(int isten,int izone,int icomp,int jcomp,FILE *ifp);
void renormalize_stencil(struct Stencil_Struct *sten,double vol_sten);
void shorten_stencil(struct Stencil_Struct *sten);
extern int Lcut_jac;
void safe_free(void **ptr);
void safe_free(void **ptr);
void merge_load_stencil(struct Stencil_Struct *sten,int **el_offsets,double *el_weights,int *el_in_radius,int *index_sten);
double get_weight_from_stencil(int isten,int icomp,int jcomp,double rsq,double R,int ngpu,double *gpu,double *gwu);
int calc_in_out_on(double *x_l,double *x_r,double sten_rad);
extern int Lhard_surf;
double calc_sten_rad(int isten,int icomp,int jcomp);
double calc_sten_vol(int isten,int i,int j);
#define DELTA_FN_BOND  7
#define POLYMER_CR     4
#define THETA_CHARGE   3
#define U_ATTRACT      2
#define THETA_FN_SIG   6
#define THETA_FN       1
#define WTC          3
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_poly;
#define DELTA_FN       0
#define POLYMER_GAUSS  5
#define NSTEN        8
extern int Sten_Type[NSTEN];
#if defined(DEC_ALPHA)
#define POW_INT powii
#endif
#if !(defined(DEC_ALPHA))
#define POW_INT (int)pow
#endif
extern double Esize_x[NDIM_MAX];
extern double Jac_grid;
extern int Coarser_jac;
extern int Sten_length_hs[3];
extern int Max_sten_length[3];
void set_gauss_quad(int ngp,double *gp,double *gw);
extern int Ndim;
extern int Nnodes_per_el_V;
extern int Ncomp;
extern int Nzone;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern struct Stencil_Struct ***Stencil;
#define VERBOSE      3 
#define NO_SCREEN    2 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
struct Stencil_Struct {
  int        Length;      /* Number of nodes that interact with current 
                             node through this stencil                    */
  int      **Offset;      /* 2D array to be set to size [Length][Ndim] that
                             gives relative position of interacting node  */
  double    *Weight;      /* 1D array of size [Length] that gives weight
                             of interacting node to total stencil         */
  double   **HW_Weight;   /* 2D array of size [Length][Nnodes_per_el] that
                             holds the weights based on what element they
                             are being contributed from. Only used for Hard
                             Walls when stencil point is a boundary node  */
};
void calc_stencils(void);
