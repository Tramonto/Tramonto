/* This file was automatically generated.  Do not edit! */
int StenTheta_GENmsa_NquadPtsGaussIntegrand(double r);
int StenTheta_RPMmsa_NquadPtsGaussIntegrand(double r);
int StenTheta_uattr_NquadPtsGaussIntegrand(double r);
int stencil_quadGaussIntegrand_switch(int sten,double r);
int StenDelta_Bond_NquadPtsGauss(double r);
int StenTheta_CrCMS_NquadPtsGauss(double r);
int StenTheta_GENmsa_NquadPtsGauss(double r);
int StenTheta_RPMmsa_NquadPtsGauss(double r);
int StenTheta_uattr_NquadPtsGauss(double r);
int StenTheta_Sigma_NquadPtsGauss(double r);
int StenTheta_R_NquadPtsGauss(double r);
int StenDelta_R_NquadPtsGauss(double r);
int stencil_quadGauss_switch(int sten,double r);
int StenDelta_Bond_NquadPtsBoundary();
int StenTheta_CrCMS_NquadPtsBoundary();
int StenTheta_GENmsa_NquadPtsBoundary();
int StenTheta_RPMmsa_NquadPtsBoundary();
int StenTheta_uattr_NquadPtsBoundary();
int StenTheta_Sigma_NquadPtsBoundary();
int StenTheta_R_NquadPtsBoundary();
int StenDelta_R_NquadPtsBoundary();
int stencil_quadBoundaryEl_switch(int sten);
double StenDelta_Bond_GetWeightFromSten(double rsq,double R);
double StenTheta_CrCMS_GetWeightFromSten(int icomp,int jcomp,double rsq,double R);
double StenTheta_GENmsa_GetWeightFromSten(int icomp,int jcomp,double rsq,int ngpu,double *gpu,double *gwu);
double StenTheta_RPMmsa_GetWeightFromSten(int icomp,int jcomp,double rsq,int ngpu,double *gpu,double *gwu);
double StenTheta_uattr_GetWeightFromSten(int icomp,int jcomp,double rsq,int ngpu,double *gpu,double *gwu);
double StenTheta_Sigma_GetWeightFromSten(double rsq,double R);
double StenTheta_R_GetWeightFromSten(double rsq,double R);
double StenDelta_R_GetWeightFromSten(double rsq,double R);
double stencil_GetWeight_switch(int sten,int icomp,int jcomp,double rsq,double sten_rad,int ngpu,double *gpu,double *gwu);
double StenDelta_Bond_sten_vol(int icomp,int jcomp);
double StenTheta_CrCMS_sten_vol(int icomp,int jcomp);
double StenTheta_GENmsa_sten_vol(int i,int j);
double StenTheta_RPMmsa_sten_vol(int i,int j);
double StenTheta_uattr_sten_vol(int i,int j);
double StenTheta_Sigma_sten_vol(int icomp);
double StenTheta_R_sten_vol(int i);
double StenDelta_R_sten_vol(int icomp);
double stencil_volume_switch(int sten,int icomp,int jcomp);
double StenDelta_Bond_sten_rad(int icomp,int jcomp);
double StenTheta_CrCMS_sten_rad(int icomp,int jcomp);
double StenTheta_GENmsa_sten_rad(int icomp,int jcomp);
double StenTheta_RPMmsa_sten_rad(int icomp);
double StenTheta_uattr_sten_rad(int icomp,int jcomp);
double StenTheta_Sigma_sten_rad(int icomp);
double StenTheta_R_sten_rad(int icomp);
double StenDelta_R_sten_rad(int icomp);
double stencil_radius_switch(int sten,int icomp,int jcomp);
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
#define SCREEN_NONE       -1 
extern int Iwrite_screen;
int StenDelta_Bond_Njcomp();
#define DELTA_FN_BOND         6
int StenTheta_CrCMS_Njcomp();
#define THETA_CR_DATA         4
int StenTheta_GENmsa_Njcomp();
#define THETA_CR_GENERAL_MSA  7
int StenTheta_RPMmsa_Njcomp();
#define THETA_CR_RPM_MSA      3
int StenTheta_uattr_Njcomp();
#define THETA_PAIRPOT_RCUT    2
int StenTheta_Sigma_Njcomp();
#define THETA_FN_SIG          5
int StenTheta_R_Njcomp();
#define THETA_FN_R            1
int StenDelta_R_Njcomp();
#define DELTA_FN_R            0
int stencil_Njcomp_switch(int sten);
