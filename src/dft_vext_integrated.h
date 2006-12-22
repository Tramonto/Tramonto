/* This file was automatically generated.  Do not edit! */
double uCOULOMB(double r,double z1,double z2);
double uCOULOMB(double r,double z1,double z2);
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
#define PI    M_PI
double pairPot_switch(double r,double param1,double param2,double param3);
double get_wt_from_sten_coul3D(double r,double z1,double z2);
double get_wt_from_sten_coul3D(double r,double z1,double z2);
extern double Vol_el;
double get_wt_from_sten(double r,double sigma,double eps,double rcut,int ngpu,double *gpu,double *gwu);
double get_wt_from_sten(double r,double sigma,double eps,double rcut,int ngpu,double *gpu,double *gwu);
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern int Ndim;
double integrate_potential(int flag,double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
double integrate_potential(int flag,double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
