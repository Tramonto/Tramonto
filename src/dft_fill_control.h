/* This file was automatically generated.  Do not edit! */
void safe_free(void **ptr);
void safe_free(void **ptr);
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
extern int Nunk_per_node;
#define NODAL_FLAG -999
typedef struct RB_Struct RB_Struct;
void fill_resid_and_matrix(double **x,struct RB_Struct *dphi_drb,int iter,int resid_only_flag,int unk_flag);
#define MATRIX_FILL_NODAL 1   /* set to zero for physics based ordering */
void FMT1stDeriv_switch(double **x,struct RB_Struct *dphi_drb);
extern int Nnodes_box;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
void fill_resid_and_matrix_control(double **x,int iter,int resid_only_flag);
