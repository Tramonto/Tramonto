/* This file was automatically generated.  Do not edit! */
typedef struct RB_Struct RB_Struct;
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
#define NDIM_MAX  3
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
struct RB_Struct d2phi_drb2_theta_rb_FMT1(double *n);
#define NCOMP_MAX 5
extern double Inv_4pirsq[NCOMP_MAX];
extern double Inv_4pi;
extern double Inv_4pir[NCOMP_MAX];
extern double Inv_rad[NCOMP_MAX];
extern double Esize_x[NDIM_MAX];
struct RB_Struct d2phi_drb2_delta_rb_FMT1(double *n,int *offset,double *sign,int icomp);
void FMT1_1stderiv(double *n,double DOT_12,double DOT_22,double *inv_n3,double *dphi_drb_loc);
#define PI    3.141592653589793238462643383279502884197169399375
extern int Nrho_bar_s;
extern int Ndim;
double FMT1_energy_density(double *n);
