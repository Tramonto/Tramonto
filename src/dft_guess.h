/* This file was automatically generated.  Do not edit! */
void check_zero_densities(double **xOwned);
void chop_profile(double **xOwned,int iguess);
void setup_polymer_G(double **xOwned);
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
#define CMS_G          2 
void setup_polymer_field(double **xOwned,int iguess);
#define CMS_FIELD      1
void setup_BondWTC(double **xOwned);
#define BONDWTC       7
void setup_Xi_cavWTC(double **xOwned);
#define CAVWTC     6
void setup_chem_pot(double **xOwned);
#define DIFFUSION      5
void setup_elec_pot(double **xOwned,int iguess);
#define POISSON        3
void setup_rho_bar(double **xOwned);
#define HSRHOBAR       4
void setup_density(double **xOwned,int iguess);
void setup_polymer_rho(double **xOwned,int iguess);
#define WTC          3
#define NONE       -1
#define NONE      -1
#define NONE -1
extern int Type_poly;
#define NEQ_TYPE       8
extern int Phys2Nunk[NEQ_TYPE];
#define DENSITY        0
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
void guess_restart_from_files(int start_no_info,int iguess,double **xOwned);
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Restart_field[NEQ_TYPE];
extern int Imain_loop;
extern int Restart;
#define VERBOSE      3 
extern int Iwrite;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void set_initial_guess(int iguess,double **xOwned);
