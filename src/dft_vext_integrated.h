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
#define PI    M_PI
double pairPot_switch(double r,double param1,double param2,double param3,int typePairPot);
extern double Vol_el;
#define PAIR_COULOMB       2
double get_wt_from_sten(double r,double param1,double param2,double rcut,int ngpu,double *gpu,double *gwu);
#define PAIR_COULOMB_CS    1
#define PAIR_LJ12_6_CS     0
extern int Type_vext3D;
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern int Ndim;
double integrate_potential(double param1,double param2,double param3,int ngp,int ngpu,double *gp,double *gpu,double *gw,double *gwu,double *node_pos,double *node_pos_f);
