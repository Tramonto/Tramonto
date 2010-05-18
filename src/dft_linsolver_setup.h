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
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
extern int **Bonds_SegAll;
extern int **Poly_to_Unk_SegAll;
extern int *Nbonds_SegAll;
extern int Nseg_tot;
#define NEQ_TYPE       12 
extern int Phys2Nunk[NEQ_TYPE];
void linsolver_setup_WJDCTYPE_LINEARONLY();
#define WJDC_FIELD     8
extern int Type_coul;
extern int Mesh_coarsening;
#define WTC          2
extern int Type_attr;
#define CAVWTC         4
#define BONDWTC        5
extern int Ndim;
extern int Nrho_bar_s;
extern int Phys2Unk_first[NEQ_TYPE];
#define HSRHOBAR       2
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Nlists_HW;
extern int Lhard_surf;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
int discover_G_ordering_LT(int *geq);
void safe_free(void **ptr);
void safe_free(void **ptr);
extern int Ncomp;
#define MF_EQ          3
#define POISSON        1
#define NCOMP_MAX 5
extern int Geqn_start[NCOMP_MAX];
extern int *Pol_Sym;
#define G_CHAIN       11 
#define CMS_FIELD      7
#define DENSITY        0
#define NMER_MAX     200
extern int Unk2Phys[3 *NCOMP_MAX+2 *NMER_MAX+NMER_MAX *NMER_MAX+13];
#if defined(__STDC__)
void *array_alloc(int numdim,...);
#endif
void *array_alloc(int numdim,...);
#if !(defined(__STDC__))
void *array_alloc(...);
#endif
void linsolver_setup_CMSTYPE_LINEARONLY();
extern void *ParameterList_list;
extern int Nunk_per_node;
extern void *LinProbMgr_manager;
void linsolver_setup_HSTYPE();
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Type_func;
void linsolver_setup_WJDCTYPE();
#define WJDC3        5 
#define WJDC         3
void linsolver_setup_CMSTYPE();
#define CMS          0
extern int Type_poly;
extern int L_Schur;
void linsolver_setup_control();
