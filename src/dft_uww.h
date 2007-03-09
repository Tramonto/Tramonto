/* This file was automatically generated.  Do not edit! */
double pairPot_switch(double r,double param1,double param2,double param3,double param4,int typePairPot);
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
#define NWALL_MAX 600 
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Ndim;
#define WALL_WALL   2
double pairPotparams_switch(int typePairPot,int context,int i,int j,double *param1,double *param2,double *param3,double *param4);
#define PAIR_COULOMB       2
extern int Type_coul;
void setup_atomic_ww(int iwall,int jwall,int type_uwwpot);
#define PAIR_HARD         -1 
extern int Type_uwwPot;
#define ATOM_CENTERS_WW    1 
#define NWALL_MAX_TYPE 50 
extern int Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int WallType[NWALL_MAX];
extern int Link[NWALL_MAX];
extern int Nlink;
extern double **Uww_link;
extern int Nwall;
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
extern double **Uww;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void setup_wall_wall_potentials();
