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
extern int LBulk;
#define NCONT_MAX          2 /* the maximum number of solutions possible for use with Loca */
extern int Cont_ID[NCONT_MAX][2];
extern int NID_Cont;
extern int Lbinodal;
typedef struct Loca_Struct Loca_Struct;
struct Loca_Struct {
  int    method;      /* Continuation method                          */
  int    cont_type1;  /* flag specifying the continuation parameter   */
  int    cont_type2;  /* flag specifying the second (free) parameter  */
  int    num_steps;   /* maximum number of continuation steps to take */
  double aggr;        /* step size control parameter                  */
  double step_size;   /* initial continuation step size               */
};
extern struct Loca_Struct Loca;
#define USE_LOCA
extern int Pos_new_nodes;
extern int Plane_new_nodes;
#define NWALL_MAX_TYPE 20 
extern double Del_1[NWALL_MAX_TYPE];
extern int Nruns;
extern double Az_tolerance;
extern int Max_gmres_iter;
extern double Az_ilut_fill_param;
extern int Az_preconditioner;
extern int Az_scaling;
extern int Az_kspace;
extern int Az_solver;
extern int L_Schur;
extern int Load_Bal_Flag;
extern double NL_abs_tol_picard,NL_rel_tol_picard;
extern double NL_update_scalingParam;
extern double NL_abs_tol,NL_rel_tol;
extern int Analyt_WJDC_Jac;
extern int ATTInA22Block;
extern int Physics_scaling;
extern int NL_Solver;
extern int Max_NL_iter;
#define NCOMP_MAX 5
extern double Scale_fac_WJDC[NCOMP_MAX][NCOMP_MAX];
extern double X_1D_bc;
extern int Dim_1Dbc;
extern int L1D_bc;
extern double Jac_threshold;
extern int Lcut_jac;
extern double Jac_grid;
extern int Coarser_jac;
extern int Mesh_coarsening;
extern double Rmax_zone[5];
extern int Nzone;
#define FILENAME_LENGTH 300
extern char EchoInputFile_array[FILENAME_LENGTH];
extern char Outpath_array[FILENAME_LENGTH];
extern char Runpath_array[FILENAME_LENGTH];
extern int Iwrite_files;
extern int Iwrite_screen;
extern int Iwrite;
extern int Print_rho_switch;
extern int Print_rho_type;
extern int Lper_area;
extern double Rho_max;
extern int Restart_Vext;
extern int Nmissing_densities;
extern int Restart;
#define NSTEPS_MAX 10
extern double Rho_step[NCOMP_MAX][NSTEPS_MAX];
extern double Xend_step[NSTEPS_MAX];
extern double Xstart_step[NSTEPS_MAX];
extern int Orientation_step[NSTEPS_MAX];
extern int Nsteps;
#define CHOP_RHO_STEP    4
#define CHOP_RHO         3
#define STEP_PROFILE     2
extern int Iguess_fields;
extern int Iguess;
extern double *Lseg_IC;
extern double *Pore_rad_R_IC;
extern double *Pore_rad_L_IC;
extern int Nseg_IC;
extern int Geom_flag;
extern double Velocity;
extern double D_coef[NCOMP_MAX];
#define DIFFUSIVE_INTERFACE 1
extern int Linear_transport;
extern double *Dielec_wall;
extern double Dielec_X;
extern double Dielec_pore;
extern double Dielec_bulk;
#define DIELEC_WF_PORE     2
extern double DielecConst_plasma;
extern double Temp_K_plasma;
extern double Sigma_Angstroms_plasma;
extern int Type_dielec;
extern int Ipot_wf_c;
extern int Charge_type_local;
extern int Charge_type_atoms;
#define ATOMIC_CHARGE    3
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern double **Charge_x;
extern double *Charge_Diam;
extern double *Charge;
extern int Nlocal_charge;
extern int Type_bc_elec[NWALL_MAX_TYPE];
#define PAIR_COULOMB          2
#define PAIR_COULOMB_CS       1
extern double X_const_mu;
extern double Elec_pot_RTF;
extern double Elec_pot_LBB;
#define COULOMB      1
extern double Rho_b_RTF[NCOMP_MAX];
extern double Rho_b_LBB[NCOMP_MAX];
#define UNIFORM_INTERFACE  0
extern double Rho_b[NCOMP_MAX];
extern int Lconstrain_interface;
extern int Grad_dim;
extern int Type_interface;
extern double **Vext_membrane;
extern int **Lsemiperm;
extern double Crfac;
extern int Ncr_files;
extern double Cr_break[2];
extern double Cr_rad_hs[NCOMP_MAX][NCOMP_MAX];
#define CMS          0
extern int ***pol_sym_tmp;
#define NBOND_MAX 4
extern int ***Bonds;
extern int **Nbond;
extern int Type_poly_arch;
extern double Rho_g[NCOMP_MAX];
extern int Graft_wall[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define NBLOCK_MAX   20 
extern int Type_mer_to_Pol[NBLOCK_MAX];
extern int Poly_to_Type[NCOMP_MAX][NBLOCK_MAX];
extern int Poly_to_Ntype[NCOMP_MAX];
#define NMER_MAX     200
extern int Type_mer[NCOMP_MAX][NMER_MAX];
extern int SegChain2SegAll[NCOMP_MAX][NMER_MAX];
extern int Nmer_t_total[NBLOCK_MAX];
extern int Nmer_t[NCOMP_MAX][NBLOCK_MAX];
extern int Ntype_mer;
extern int Nmer[NCOMP_MAX];
extern int Nblock[NCOMP_MAX];
extern int Npol_comp;
#define NONE       -1
#define NONE          -1
#define NONE        -1
#define NONE        -1
extern int Unk2Comp[NMER_MAX];
extern double YukawaK_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double EpsYukawa_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Eps_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Sigma_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double YukawaK_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double EpsYukawa_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Eps_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Sigma_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern double Rho_w[NWALL_MAX_TYPE];
extern double Npow_ff[NCOMP_MAX][NCOMP_MAX];
extern double YukawaK_ff[NCOMP_MAX][NCOMP_MAX];
extern double EpsYukawa_ff[NCOMP_MAX][NCOMP_MAX];
#define PAIR_rNandYUKAWA_CS   9
#define PAIR_r18andYUKAWA_CS  8
#define PAIR_r12andYUKAWA_CS  7
#define PAIR_LJandYUKAWA_CS   6
#define PAIR_EXP_CS	      4
#define PAIR_YUKAWA_CS        3
extern double Bond_ff[NCOMP_MAX][NCOMP_MAX];
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double Eps_ff[NCOMP_MAX][NCOMP_MAX];
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int Lpolarize[NCOMP_MAX];
extern double Pol[NCOMP_MAX];
extern double Charge_f[NCOMP_MAX];
extern double Mass[NCOMP_MAX];
extern double HS_diam[NCOMP_MAX];
#define MANUAL_HS_DIAM         2
extern int Mix_type;
extern int Ncomp;
extern int Type_uwwPot;
extern int Ipot_ww_n[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
extern int Vext_PotentialID[NWALL_MAX_TYPE];
extern int Type_vext[NWALL_MAX_TYPE];
extern int Ipot_wf_n[NWALL_MAX_TYPE];
extern int Lhard_surf;
#define NWALL_MAX 600 
extern double Angle_wedge_end[NWALL_MAX];
extern double Angle_wedge_start[NWALL_MAX];
extern int Lwedge_cutout[NWALL_MAX];
extern int read_wedge;
#define NPERIODIC_MAX 4
extern double EndpointLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double OriginLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double SlopeLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int OrientationLinearFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int Llinear_overlay[NWALL_MAX_TYPE];
extern int Nlinear_overlay[NWALL_MAX_TYPE];
extern int read_linear;
extern double OriginPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double WavelengthPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern double AmplitudePeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int OrientationPeriodicFunc[NWALL_MAX_TYPE][NPERIODIC_MAX];
extern int Lperiodic_overlay[NWALL_MAX_TYPE];
extern int Nperiodic_overlay[NWALL_MAX_TYPE];
extern int read_periodic;
extern double Rough_length[NWALL_MAX_TYPE];
extern double Rough_param_max[NWALL_MAX_TYPE];
#define TRUE  1
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
extern int Lrough_surf[NWALL_MAX_TYPE];
extern int read_rough;
extern int Lapply_offset[3];
extern double WallParam_3[NWALL_MAX_TYPE];
extern double WallParam_2[NWALL_MAX_TYPE];
extern double WallParam[NWALL_MAX_TYPE];
extern double Elec_param_w[NWALL_MAX];
#define NDIM_MAX  3
extern double WallPos[NDIM_MAX][NWALL_MAX];
extern int Link[NWALL_MAX];
extern int WallType[NWALL_MAX];
extern int Orientation[NWALL_MAX_TYPE];
extern int Surface_type[NWALL_MAX_TYPE];
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
void *array_alloc(int numdim,...);
extern int **Xtest_reflect_TF;
extern int Proc;
#if defined(DEBUG)
extern int Proc;
#endif
extern int Nlink;
extern int Nwall;
extern int Nwall_type;
extern int Ipot_ff_c;
extern int Ipot_ff_n;
extern int Type_poly;
extern int Type_coul;
extern int Type_pairPot;
extern int Type_attr;
extern int Type_hsdiam;
extern int Type_func;
extern int Type_bc[NDIM_MAX][2];
extern double Esize_x[NDIM_MAX];
extern double Size_x[NDIM_MAX];
extern int Ndim;
extern double VEXT_MAX;
extern double Dielec_ref;
extern double Temp;
extern double Density_ref;
extern double Length_ref;
void broadcast_input();
