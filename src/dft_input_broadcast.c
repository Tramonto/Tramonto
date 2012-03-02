/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_input_broadcast.c
 *
 *  This file contains broadcast directives to start up multiprocessor runs
 *  for all input parameters.
 */

#include "dft_input_broadcast.h"

/******************************************************************************/

void broadcast_input()

/*
 *    Author:  Laura Frink
 */

{
   /* Local variable declarations */
   
   int icomp, iwall,iwall_type, idim, pol_number,nseg,nmer_max,iseg,ibond,ncharge,logical;
  
  /********************** BEGIN EXECUTION ************************************/

  /* Broadcast all parameters that have been set up either in the input file
     or the GUI - in either case only one proc obtains the input so it needs to 
     be passed to all processors */


  /**********************************/
  /* Broadcast Dimension Parameters */
  /**********************************/
  MPI_Bcast(&Length_ref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Density_ref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Temp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Dielec_ref,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&VEXT_MAX,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  /*****************************/
  /* Broadcast Mesh Parameters */
  /*****************************/

  MPI_Bcast(&Ndim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(Size_x,NDIM_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Esize_x,NDIM_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Type_bc,Ndim*2,MPI_INT,0,MPI_COMM_WORLD);

  /***********************************/
  /* Broadcast Functional Selections */
  /***********************************/

                               /* hard sphere functionals */
  MPI_Bcast(&Type_func,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Type_hsdiam,1,MPI_INT,0,MPI_COMM_WORLD);

                               /* attractive functionals */
  MPI_Bcast(&Type_attr,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Type_pairPot,1,MPI_INT,0,MPI_COMM_WORLD);

                               /* coulombic functionals */
  MPI_Bcast(&Type_coul,1,MPI_INT,0,MPI_COMM_WORLD);

                               /* polymer functionals */
  MPI_Bcast(&Type_poly,1,MPI_INT,0,MPI_COMM_WORLD);

                                /* other potential types related to functional choices */
  MPI_Bcast(&Ipot_ff_n,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Ipot_ff_c,1,MPI_INT,0,MPI_COMM_WORLD);

  /********************************/
  /* Broadcast Surface Parameters */
  /********************************/
  MPI_Bcast(&Nwall_type,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nwall,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nlink,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Nwall>0){
     if (Proc>0) Xtest_reflect_TF = (int **) array_alloc (2, Nlink,Ndim, sizeof(int));
     MPI_Bcast(*Xtest_reflect_TF,Nlink*Ndim,MPI_INT,0,MPI_COMM_WORLD);
  }

  if(Nwall_type > 0){
      MPI_Bcast(Surface_type,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Orientation,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(WallType,NWALL_MAX,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Link,NWALL_MAX,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(WallPos,NDIM_MAX*NWALL_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(Elec_param_w,NWALL_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(WallParam,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(WallParam_2,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(WallParam_3,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

                 /* Rough Surface Definition Params */
      MPI_Bcast(Lapply_offset,3,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&read_rough,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Lrough_surf,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      if (read_rough==TRUE){ 
          MPI_Bcast(Rough_param_max,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Bcast(Rough_length,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
                 /* Periodic Overlay Function Params */
      MPI_Bcast(&read_periodic,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Nperiodic_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Lperiodic_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      if (read_periodic==TRUE) {
           MPI_Bcast(OrientationPeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_INT,0,MPI_COMM_WORLD);
           MPI_Bcast(AmplitudePeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Bcast(WavelengthPeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Bcast(OriginPeriodicFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
                 /* Linear Overlay Function Params */
      MPI_Bcast(&read_linear,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Nlinear_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Llinear_overlay,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      if (read_linear==TRUE) {
           MPI_Bcast(OrientationLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_INT,0,MPI_COMM_WORLD);
           MPI_Bcast(SlopeLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Bcast(OriginLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Bcast(EndpointLinearFunc,NWALL_MAX_TYPE*NPERIODIC_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
                 /* Angle Cutout Params */
      MPI_Bcast(&read_wedge,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(Lwedge_cutout,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
      if (read_wedge==TRUE) {
          MPI_Bcast(Angle_wedge_start,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Bcast(Angle_wedge_end,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
  }

  /*******************************************************************************/
  /* Broadcast Wall-Fluid and Wall-Wall Interaction Type Parameters    */   
  /*******************************************************************************/

  MPI_Bcast(&Lhard_surf,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Nwall_type > 0){
       MPI_Bcast(Ipot_wf_n,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(Type_vext,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(Vext_PotentialID,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(Ipot_ww_n,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
       if (Ndim==3) MPI_Bcast(&Type_uwwPot,1,MPI_INT,0,MPI_COMM_WORLD);
  }

  /************************************************/
  /* Broadcast Fluid-Fluid Interaction Parameters */ 
  /************************************************/

  MPI_Bcast(&Ncomp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Mix_type,1,MPI_INT,0,MPI_COMM_WORLD);
  if (Type_hsdiam==MANUAL_HS_DIAM) MPI_Bcast(HS_diam,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Mass,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Charge_f,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Pol,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Lpolarize,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(Sigma_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Eps_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Cut_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(Bond_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  if (Type_pairPot==PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS || 
      Type_pairPot==PAIR_LJandYUKAWA_CS || Type_pairPot==PAIR_r12andYUKAWA_CS ||
      Type_pairPot==PAIR_r18andYUKAWA_CS || Type_pairPot==PAIR_rNandYUKAWA_CS) {
          MPI_Bcast(EpsYukawa_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	  MPI_Bcast(YukawaK_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	  MPI_Bcast(Npow_ff,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }

  /**********************************************/
  /* Broadcast Wall-Wall Interaction Parameters */ 
  /**********************************************/

  if (Nwall_type>0){
     MPI_Bcast(Rho_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
     if (Mix_type==1){
            MPI_Bcast(Sigma_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(Eps_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     else {      
            MPI_Bcast(Sigma_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Bcast(Eps_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     MPI_Bcast(Cut_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);

     if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
         Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
         || Type_uwwPot==PAIR_r18andYUKAWA_CS){
        if (Mix_type==1){
           MPI_Bcast(EpsYukawa_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
           MPI_Bcast(YukawaK_ww,NWALL_MAX_TYPE*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        else{
          MPI_Bcast(EpsYukawa_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Bcast(YukawaK_w,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
     }
  }


  /**********************************************/
  /* Broadcast Wall-Fluid Interaction Parameters */ 
  /**********************************************/

                       /* WALL-FLUID PARAMS -- ONLY IF MIX_TYPE == 1 */

  if (Mix_type == 1){
    MPI_Bcast(Sigma_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Eps_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Cut_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(EpsYukawa_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(YukawaK_wf,NCOMP_MAX*NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else{ /* nothing to do - all wall-fluid parameters computing using LB mixing rules */ }

  /**************************************/
  /* Broadcast Polymer Input Parameters */ 
  /**************************************/

  MPI_Bcast(Unk2Comp,NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);

  if (Type_poly != NONE){
    MPI_Bcast(&Npol_comp,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nblock,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Ntype_mer,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer_t,NCOMP_MAX*NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Nmer_t_total,NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(SegChain2SegAll,NCOMP_MAX*NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Type_mer,NCOMP_MAX*NMER_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Poly_to_Ntype,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Poly_to_Type,NCOMP_MAX*NBLOCK_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Type_mer_to_Pol,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Grafted,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Graft_wall,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Rho_g,NCOMP_MAX,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Bcast(&Type_poly_arch,1,MPI_INT,0,MPI_COMM_WORLD);

    nseg=nmer_max=0;
    for (pol_number=0; pol_number<Npol_comp; pol_number++){
         nseg += Nmer[pol_number];
         if (Nmer[pol_number] > nmer_max) nmer_max=Nmer[pol_number];
    }
    if (Proc>0){
       Nbond = (int **) array_alloc (2, Npol_comp,nmer_max,sizeof(int));
       Bonds = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
       pol_sym_tmp = (int ***) array_alloc (3, Npol_comp,nmer_max,NBOND_MAX,sizeof(int));
    }

    for (pol_number=0; pol_number<Npol_comp;pol_number++){
      for (iseg=0;iseg<Nmer[pol_number];iseg++){
          MPI_Bcast(&Nbond[pol_number][iseg],1,MPI_INT,0,MPI_COMM_WORLD);
          for (ibond=0;ibond<Nbond[pol_number][iseg];ibond++){
             MPI_Bcast(&Bonds[pol_number][iseg][ibond],1,MPI_INT,0,MPI_COMM_WORLD);
             MPI_Bcast(&pol_sym_tmp[pol_number][iseg][ibond],1,MPI_INT,0,MPI_COMM_WORLD);
          }
       }
     }

    if (Type_poly == CMS){  /* this bit only applies to the CMS functional */
       MPI_Bcast(Cr_rad_hs,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
       MPI_Bcast(Cr_break,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
       MPI_Bcast(&Ncr_files,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&Crfac,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
  }     
  /* end of polymer input */

  /**********************************************/
  /* Broadcast Input for Semi-Permeable Surfaces */ 
  /**********************************************/

  if (Nwall_type>0){
     if (Proc>0){
         Lsemiperm = (int **) array_alloc (2, Nwall_type,Ncomp,sizeof(int));
         Vext_membrane = (double **) array_alloc (2, Nwall_type,Ncomp,sizeof(double));
     }

     MPI_Bcast(*Lsemiperm,Nwall_type*Ncomp,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(*Vext_membrane,Nwall_type*Ncomp,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  /**********************************************/
  /* Broadcast Input for State Point Parameters */ 
  /**********************************************/

  MPI_Bcast(&Type_interface,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Grad_dim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lconstrain_interface,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(Rho_b,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (Type_interface != UNIFORM_INTERFACE){
     MPI_Bcast(Rho_b_LBB,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(Rho_b_RTF,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  if (Type_interface != UNIFORM_INTERFACE) {
     if (Ipot_ff_c == COULOMB){
         MPI_Bcast(&Elec_pot_LBB,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
         MPI_Bcast(&Elec_pot_RTF,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     MPI_Bcast(&X_const_mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  /**********************************************/
  /* Broadcast Input for Charged Systems        */ 
  /**********************************************/
  if (Type_coul != NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
     if (Nwall_type>0) MPI_Bcast(Type_bc_elec,NWALL_MAX_TYPE,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&Nlocal_charge,1,MPI_INT,0,MPI_COMM_WORLD);
     ncharge=0;
     if      (Nlocal_charge > 0)  ncharge = Nlocal_charge;
     else if (Nlocal_charge < 0)  ncharge = 2;
     if (Nlocal_charge>0){

        if (Proc>0){
           Charge      = (double *) array_alloc (1, ncharge, sizeof(double));
           Charge_Diam = (double *) array_alloc (1, ncharge, sizeof(double));
           Charge_x = (double **) array_alloc (2, Ndim,ncharge,sizeof(double));
        }
        MPI_Bcast(Charge,ncharge,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(Charge_Diam,ncharge,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(*Charge_x,Ndim*ncharge,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
    logical=FALSE;
    for (iwall_type=0; iwall_type<Nwall_type; iwall_type++)
      if (Type_bc_elec[iwall_type] == ATOMIC_CHARGE) logical=TRUE;
      if (ncharge !=0 || logical){
         MPI_Bcast(&Charge_type_atoms,1,MPI_INT,0,MPI_COMM_WORLD);
         MPI_Bcast(&Charge_type_local,1,MPI_INT,0,MPI_COMM_WORLD);
      }
  }
  MPI_Bcast(&Ipot_wf_c,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Type_coul!=NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
    MPI_Bcast(&Type_dielec,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Sigma_Angstroms_plasma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Temp_K_plasma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&DielecConst_plasma,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(Type_dielec != DIELEC_WF_PORE){ 
      MPI_Bcast(&Dielec_bulk,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    }
    else{
      MPI_Bcast(&Dielec_bulk,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&Dielec_pore,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&Dielec_X,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    }
    if (Nwall_type>0){
      if (Proc>0) Dielec_wall = (double *) array_alloc (1, Nwall_type,sizeof(double));
      MPI_Bcast(Dielec_wall,Nwall_type,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
  }
  
  /************************************************/
  /* Broadcast Input for Diffusing Systems        */ 
  /************************************************/
  MPI_Bcast(&Linear_transport,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Type_interface==DIFFUSIVE_INTERFACE) {
    MPI_Bcast(D_coef,NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Velocity,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Geom_flag,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Nseg_IC,1,MPI_INT,0,MPI_COMM_WORLD);
 
    if (Proc>0){
       Pore_rad_L_IC = (double *) array_alloc (1, Nseg_IC,sizeof(double));
       Pore_rad_R_IC = (double *) array_alloc (1, Nseg_IC,sizeof(double));
       Lseg_IC       = (double *) array_alloc (1, Nseg_IC,sizeof(double));
    }
    MPI_Bcast(Pore_rad_L_IC,Nseg_IC,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Pore_rad_R_IC,Nseg_IC,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Lseg_IC,Nseg_IC,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  /***********************************************************/
  /* Broadcast Run Control and Initial Guess Selections      */ 
  /***********************************************************/

  MPI_Bcast(&Iguess,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Iguess_fields,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Iguess==STEP_PROFILE || (Iguess>=CHOP_RHO && Iguess<= CHOP_RHO_STEP)){
    MPI_Bcast(&Nsteps,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Orientation_step,NSTEPS_MAX,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(Xstart_step,NSTEPS_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Xend_step,NSTEPS_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Rho_step,NSTEPS_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  MPI_Bcast(&Restart,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Nmissing_densities,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Restart_Vext,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Rho_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /************************************/
  /* Broadcast Output Format Settings */
  /************************************/

  MPI_Bcast(&Lper_area,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Print_rho_type,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Print_rho_switch,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Iwrite,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Iwrite_screen,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Iwrite_files,1,MPI_INT,0,MPI_COMM_WORLD);

  /********************************************************/
  /* Broadcast Settings for Numerical Methods: Coarsening */
  /********************************************************/

  MPI_Bcast(&Nzone,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(Rmax_zone,5,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Mesh_coarsening,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Coarser_jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Jac_grid,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lcut_jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Jac_threshold,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&L1D_bc,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&X_1D_bc,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /*************************************************************************/
  /* Broadcast Settings for Numerical Methods: Nonlinear Solver Parameters */
  /*************************************************************************/

  MPI_Bcast(Scale_fac_WJDC,NCOMP_MAX*NCOMP_MAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Max_NL_iter,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_Solver,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Physics_scaling,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&ATTInA22Block,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Analyt_WJDC_Jac,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_rel_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_abs_tol,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_update_scalingParam,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_rel_tol_picard,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NL_abs_tol_picard,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Load_Bal_Flag,1,MPI_INT,0,MPI_COMM_WORLD);

  /*************************************************************************/
  /* Broadcast Settings for Linear Solver Parameters                       */
  /*************************************************************************/

  MPI_Bcast(&L_Schur,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_solver,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_kspace,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_scaling,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_preconditioner,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_ilut_fill_param,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Max_gmres_iter,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Az_tolerance,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /******************************************************************/
  /* Broadcast Settings for Mesh Continuation                       */
  /******************************************************************/
  /* READ IN MESH CONTINUATION PARAMETERS */

  MPI_Bcast(&Nruns,1,MPI_INT,0,MPI_COMM_WORLD);

  if (Nruns > 0){
     MPI_Bcast(Del_1,NWALL_MAX_TYPE,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(&Plane_new_nodes,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&Pos_new_nodes,1,MPI_INT,0,MPI_COMM_WORLD);
  }

  /******************************************************************/
  /* Broadcast Settings for LOCA Continuation                       */
  /******************************************************************/
    /* FINALLY READ IN LOCA PARAMETERS */
#ifdef USE_LOCA
  MPI_Bcast(&Loca.method,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Lbinodal,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Loca.cont_type1,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&NID_Cont,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(Cont_ID,NCONT_MAX*2,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Loca.step_size,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Loca.num_steps,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&Loca.aggr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Loca.cont_type2,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&LBulk,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

return;
}
/****************************************************************************/
