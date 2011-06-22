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
#define PI    3.141592653589793238462643383279502884197169399375
double pairPot_switch(double r,double param1,double param2,double param3,double param4,double param5,double param6,int typePairPot);
extern double Vol_el;
double get_wt_from_sten(int typePot,double r,double param1,double param2,double param3,double param4,double param5,double param6,int ngpu,double *gpu,double *gwu);
#define PAIR_COULOMB          2
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern int Ndim;
double integrate_potential(int typePot,double param1,double param2,double param3,double param4,double param5,double param6,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
