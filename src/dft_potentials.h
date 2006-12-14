/* This file was automatically generated.  Do not edit! */
double uCOULOMB_att_int(double r,int i,int j);
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
#define NCOMP_MAX 5
extern double Charge_f[NCOMP_MAX];
double uCOULOMB_att(double r,int i,int j);
double uLJatt_n_noshift(double r,int i,int j);
double uLJatt_n_int(double r,int i,int j);
double uLJatt_n(double r,int i,int j);
extern double Temp_elec;
double uLJ_wp(double r,int icomp,int iwall_type);
double uderiv_LJ12_6(double r,double x,double sigma,double eps,double rcut);
#if defined(DEC_ALPHA)
#define POW_DOUBLE_INT powi
#endif
#if !(defined(DEC_ALPHA))
#define POW_DOUBLE_INT pow
#endif
double uCOULOMB(double r,double z1,double z2);
#define PI    M_PI
double uLJ12_6_cut(double r,double sigma,double eps,double rcut);
double get_wt_from_sten_coul3D(double r,double z1,double z2);
extern double Vol_el;
double get_wt_from_sten(double r,double sigma,double eps,double rcut,int ngpu,double *gpu,double *gwu);
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern int Ndim;
double integrate_potential(int flag,double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
#define NWALL_MAX_TYPE 50 
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_w[NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_w[NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int Nwall_type;
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_poly;
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
#define IDEAL_GAS    0
extern int Ipot_ff_n;
extern int Ncomp;
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
void pot_parameters(char *output_file1);
