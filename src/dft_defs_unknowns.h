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
#define VERBOSE      3 
extern int Iwrite;
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int Npol_comp;
#define NO_UNK        -888
extern int Nunk_per_node;
#define NMER_MAX     100
extern int Unk2Phys[3 *NCOMP_MAX+NMER_MAX+NMER_MAX *NMER_MAX+13];
#define NEQ_TYPE       8
extern int Phys2Unk_last[NEQ_TYPE];
extern int Phys2Unk_first[NEQ_TYPE];
extern int Ngeqn_tot;
#define CMS_G          2 
#define CMS_FIELD      1
extern int Nbonds;
extern int Nrho_bar_bond;
#define BONDWTC       7
extern int Nrho_bar_cavity;
#define CAVWTC     6
extern int Lsteady_state;
extern int Ndiffusion;
#define DIFFUSION      5
#define NONE       -1
#define NONE      -1
#define NONE -1
#define NONE        -1
extern int Type_coul;
extern int Npoisson;
#define POISSON        3
extern int Nrho_bar_s;
extern int Ndim;
#define IDEAL_GAS    0
extern int Ipot_ff_n;
extern int Nrho_bar;
#define HSRHOBAR       4
extern int Ncomp;
extern int Nseg_tot;
extern int Phys2Nunk[NEQ_TYPE];
#define DENSITY        0
extern int L_HSperturbation;
#define CMS_SCFT     1
#define CMS          0
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Lseg_densities;
#define WTC          2
extern int Type_poly;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
void setup_nunk_per_node(char *output_file1);
