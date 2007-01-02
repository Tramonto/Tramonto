/* This file was automatically generated.  Do not edit! */
typedef struct RB_Struct RB_Struct;
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
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
struct RB_Struct d2phi_drb2_theta_rb_FMT3(double *n);
struct RB_Struct d2phi_drb2_theta_rb_FMT2(double *n);
struct RB_Struct d2phi_drb2_theta_rb_FMT1(double *n);
struct RB_Struct FMT2ndDerivTheta_switch(double *n);
struct RB_Struct d2phi_drb2_delta_rb_FMT3(double *n,int *offset,double *sign,int icomp);
struct RB_Struct d2phi_drb2_delta_rb_FMT2(double *n,int *offset,double *sign,int icomp);
struct RB_Struct d2phi_drb2_delta_rb_FMT1(double *n,int *offset,double *sign,int icomp);
struct RB_Struct FMT2ndDerivDelta_switch(double *n,int *offset,double *sign,int icomp);
void FMT3_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
void FMT2_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
void FMT1_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
void calc_FMT_derivatives(void(*fp_FMTderiv)(double *,double,double,double *,double *),int inode_box,double **x,struct RB_Struct *dphi_drb);
extern int Nnodes_box;
void FMT1stDeriv_switch(int inode_box,double **x,struct RB_Struct *dphi_drb);
double FMT3_energy_density(double *n);
#define FMT3       2
double FMT2_energy_density(double *n);
#define FMT2       1
double FMT1_energy_density(double *n);
#define FMT1       0
extern int Type_func;
double phispt_switch(double *n);
