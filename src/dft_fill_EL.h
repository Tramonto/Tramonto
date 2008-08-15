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
#define NCOMP_MAX 5
extern double Pol[NCOMP_MAX];
#define NDIM_MAX  3
extern double Esize_x[NDIM_MAX];
extern int ***Surf_normal;
extern int **Nodes_2_boundary_wall;
#define KAPPA_H2O 78.5
#define PI    M_PI
extern double Temp_elec;
extern int Ndim;
extern int Lpolarize[NCOMP_MAX];
extern double Charge_f[NCOMP_MAX];
#define POISSON        1
extern int Nwall;
#define DIFFUSION      6
extern double Betamu[NCOMP_MAX];
#define NMER_MAX     100
extern double Betamu_seg[NMER_MAX];
extern double *Deltac_b;
extern double Betamu_wtc[NMER_MAX];
extern double Betamu_att[NCOMP_MAX];
#define LJ12_6       2
extern double Betamu_hs_ex[NCOMP_MAX];
#define IDEAL_GAS    0
extern int Ipot_ff_n;
extern int Iliq_vap;
extern int Lsteady_state;
extern int SegAll_to_Poly[NMER_MAX];
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern int **Nseg_type_pol;
extern int Npol_comp;
extern double Field_WJDC_b[NMER_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Rho_seg_b[NMER_MAX];
extern int Lseg_densities;
extern double VEXT_MAX;
extern double **Vext;
extern int **Zero_density_TF;
#define BONDWTC        5
extern int *Pol_Sym;
extern int **Poly_to_Unk_SegAll;
extern int *Nbonds_SegAll;
double load_polyWJDC_cavityEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
double load_polyTC_cavityEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
double load_polyTC_bondEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
double load_polyTC_diagEL(int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
#define THETA_CR_RPM_MSA      3
#define DELTAC     1 
#define THETA_PAIRPOT_RCUT    2
double load_mean_field(int sten_type,int iunk,int loc_inode,int icomp,int izone,int *ijk_box,double **x,int resid_only_flag);
extern void *LinProbMgr_manager;
#define CALC_RESID_ONLY  3
#define MF_EQ          3
#define MF_VARIABLE  2
extern int Type_attr;
#define THETA_FN_R            1
#define DELTA_FN_R            0
typedef struct RB_Struct RB_Struct;
double load_nonlocal_hs_rosen_rb(int sten_type,int iunk,int loc_inode,int inode_box,int icomp,int izone,int *ijk_box,double **x,struct RB_Struct *dphi_drb,int resid_only_flag);
extern int Type_func;
#define FLAG_PBELEC -777
double fill_EL_elec_field(int iunk,int icomp,int loc_inode,int inode_box,double **x,int resid_only_flag);
#define NONE       -1
#define NONE      -1
#define NONE        -1
#define NONE        -1
extern int Type_coul;
double fill_EL_ext_field(int iunk,int icomp,int loc_inode,int resid_only_flag);
double fill_EL_chem_pot(int iunk,int icomp,int iseg,int loc_inode,int inode_box,int mesh_coarsen_flag_i,double **x,int resid_only_flag);
double fill_EL_ideal_gas(int iunk,int icomp,int loc_inode,int inode_box,double **x,int resid_only_flag);
double fill_sym_WTC(int iunk,int iseg,int loc_inode,int inode_box,double **x,int resid_only_flag);
extern int *Pol_Sym_Seg;
#define INIT_GUESS_FLAG  2
double fill_bulk_density(int iunk,int icomp,int iseg,int loc_inode,int inode_box,double **x,int resid_only_flag);
double fill_bulk_field(int iunk,int icomp,int iseg,int loc_inode,int inode_box,double **x,int resid_only_flag);
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define FLAG_BULK   -888
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
double fill_zero_value(int iunk,int loc_inode,int inode_box,double **x,int resid_only_flag);
int check_zero_density_EL(int iunk,int icomp,int iseg,int loc_inode,int inode_box,double **x);
extern int Unk2Comp[NMER_MAX];
#define WTC          2
#define DENSITY        0
#define WJDC_FIELD     8
#define NEQ_TYPE       13 
extern int Phys2Unk_first[NEQ_TYPE];
#define WJDC3        5 
#define WJDC2        4 
#define WJDC         3
extern int Type_poly;
struct RB_Struct {
  double    S0;      /*   1/(4*pi*Ri*Ri) * Delta_fn   */
  double    S1;      /*   1/(4*pi*Ri)    * Delta_fn   */
  double    S2;      /*                    Delta_fn   */
  double    S3;      /*                    Theta_fn   */
  double    V1[NDIM_MAX];      /*  1/(4*pi*Ri) * unit_vec * Delta_Fn   */
  double    V2[NDIM_MAX];      /*                unit_vec * Delta_Fn   */
};
double load_euler_lagrange(int iunk,int loc_inode,int inode_box,int *ijk_box,int izone,double **x,struct RB_Struct *dphi_drb,int mesh_coarsen_flag_i,int resid_only_flag);
