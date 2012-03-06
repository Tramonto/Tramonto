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
 *  FILE: dft_setup_params.c
 *
 *  This file contains routines that are required to set up a DFT
 *  run.  The run may start from an input file, or from a GUI.  
 *  Setup includes setting basic variables and computing other
 *  parameters needed by Tramonto.
 */

#include "dft_setup_params.h"

#ifdef BUILD_GUI
#include "GUI/dft_GUI.h"
#endif

/******************************************************************************/
void setup_params_for_dft(char *input_file, char *file_echoinput)
{
  int input_file_exists=TRUE,i,j,nseg,nmer_max,pol_number;
  FILE *fpinput, *fpecho;

  /****************************************************/
  /* first do all the things that only Proc 0 will do */
  /****************************************************/
  if (Proc==0) {
                /*********************************/
                /* attempt to open an input file */
                /*********************************/
    if( (fpinput  = fopen(input_file,"r")) == NULL) {
      printf("Can't open file %s\n", input_file);
#ifdef BUILD_GUI
      printf("Must enter all params from the GUI\n");
      Open_GUI=TRUE;
#else
      exit(-1);
#endif
      input_file_exists=FALSE;
    }
    if( (fpecho = fopen(file_echoinput,"w+")) == NULL) {
      printf("Can't open file %s\n", file_echoinput);
      exit(-1);
    }
    fprintf(fpecho,"test out the printing %s\n",file_echoinput);

    if (input_file_exists==TRUE) read_input_file(fpinput,fpecho);


        /******************************************************************/
        /* read a parameter from file to indicate if GUI should be opened */ 
        /* temporarily set always to TRUE */
        /******************************************************************/

        /*****************************/
        /* open the GUI if requested */ 
        /*****************************/
#ifdef BUILD_GUI
    if (Open_GUI==TRUE) dft_OptikaGUI();
    else{
        if (Temp>0. && Type_coul!=NONE) Flag_mV_elecpot=TRUE;
        else         Flag_mV_elecpot=FALSE;
    }
#else
    if (Open_GUI==TRUE){
       printf("GUI request is ignored because the code was not built with a GUI\n");
    }
    if (Temp>0. && Type_coul!=NONE) Flag_mV_elecpot=TRUE;
    else         Flag_mV_elecpot=FALSE;
#endif 
  }

   /*******************************************************************/
   /* now broadcast all of the basic input - regardless of whether it */
   /*                came from an input file or from the GUI          */ 
   /*******************************************************************/

   if (Num_Proc>1) broadcast_input();
       
   /*********************************************************************/
   /* finally finish setting up some parameters, arrays, and structures */
   /* that are based on in the input parameters and are needed for dft. */ 
   /*********************************************************************/

   if (Length_ref>0.0) make_length_params_dimensionless();
   if (Temp>0.0 && Flag_mV_elecpot){
      make_energy_params_dimensionless();
      if (Type_coul != NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
          Potential_ref=1000.*KBOLTZ*Temp/E_CONST;  /* reference potential in mV */
          make_elecPot_params_dimensionless(); 
      }
   }

   if (Type_coul != NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
     if (Dielec_ref>0.0) make_dielecConst_params_dimensionless();
   }

   if (Density_ref>0.0) make_density_params_dimensionless();

   if (Type_poly != NONE){
        nseg=nmer_max=0;
        for (pol_number=0; pol_number<Npol_comp; pol_number++){
           nseg += Nmer[pol_number];
           if (Nmer[pol_number] > nmer_max) nmer_max=Nmer[pol_number];
        }
        setup_chain_indexing_arrays(nseg,nmer_max,fpecho);
        safe_free((void *) &pol_sym_tmp);
   }

   /*
    * Allocate the SurfaceGeometry variable to be an array of structures
    * So it can be accessed as SGeom[iwall_type].param
    *
    * Then load the structure with surface geometry related parameters.
    */
   if (Nwall_type>0){
     SGeom = (struct SurfaceGeom_Struct *) array_alloc(1,Nwall_type, sizeof(struct SurfaceGeom_Struct));
     fill_surfGeom_struct();
   }

   setup_other_run_constants();
   if (Proc==0 && Iwrite_screen==SCREEN_VERBOSE){
     printf("\n--------------------------------------------------------------------");
     printf("\n    Successfully Completed Problem Setup \n");
     printf("\n--------------------------------------------------------------------\n");
   }
   return;
}
/******************************************************************************/
void setup_other_run_constants()
{
  int i,iwall,iwall_type,iblock,jblock,irand,irand_range,icomp;

  Lmesh_refine=FALSE;    /* used for automated mesh refinement - not currently used */

                         /****************************************************************/
                         /* set up Linked List for surfaces that are bonded or connected */
                         /****************************************************************/
  Link_list = (int **) array_alloc (2, Nlink,Nwall,sizeof(int));
  Nwall_this_link = (int *) array_alloc (1, Nlink,sizeof(int));
  for (i=0; i<Nlink; i++) Nwall_this_link[i]=0;
  for (iwall=0; iwall<Nwall; iwall++){Link_list[Link[iwall]][Nwall_this_link[Link[iwall]]++]=iwall; } 

                /*************************************************************************/
                /* For Rough Surfaces.....populate a roughness array with random numbers */
                /*   that can be used to generate roughness profiles */
                /*************************************************************************/
  
  if (read_rough){
     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
       for (iblock=0;iblock<MAX_ROUGH_BLOCK;iblock++){
         for (jblock=0;jblock<MAX_ROUGH_BLOCK;jblock++){
            if (Proc==0){
            #ifndef _MSC_VER
              irand = random();
            #else
              irand = rand();
            #endif
            }
            MPI_Bcast(&irand,1,MPI_INT,0,MPI_COMM_WORLD);  /* need this for parallel jobs to be identical to serial jobs */
            irand_range = POW_INT(2,31)-1;
            Rough_precalc[iwall_type][iblock][jblock]= Rough_param_max[iwall_type]*(-0.5+( ((double)irand)/((double)irand_range)));
         }
       }
     }
  }

                /***************************************************/
                /* Calculate a Total Density for SCFT calculations */
                /***************************************************/
  if(Type_poly==CMS_SCFT){
    Rho_t = 0.0;
    for(icomp=0; icomp<Ncomp; icomp++) Rho_t += Rho_b[icomp];
  }
  
  if (Type_coul!=NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
       Temp_elec = 4.0*PI*KBOLTZ*Temp_K_plasma*DielecConst_plasma*EPSILON_0*Sigma_Angstroms_plasma*1.e-10/(E_CONST*E_CONST);
      /*Temp_elec = 4*PI*KBOLTZ*298.0*KAPPA_H2O*EPSILON_0*4.25e-10/(E_CONST*E_CONST);  Tang-Davis Paper Parameters*/
       if (Proc==0 && (Iwrite_screen!=SCREEN_NONE || Iwrite_screen!=SCREEN_ERRORS_ONLY)) printf("\t plasma parameter=%9.6f\n",1./Temp_elec);
  }

  Nnodes_per_el_V = POW_INT(2, Ndim);
  Nnodes_per_el_S = POW_INT(2, Ndim-1);

  /* with new units, we want to retain the ability to do continuation
     in temperature independent from the energy parameters Eps_ff or
     Eps_wf.  Therefore, we now set the reduced temperature for
     the purposes of continuation.  It is not used elsewhere in the
     calculations */
   
  if (Ipot_ff_n <= HARD_SPHERE) Temp = 1.0;
  else{ Temp = 1.0/Eps_ff[0][0];}

  if (Coarser_jac == JAC_ZONES_SETFIXED_ESIZE) Nzone += 1;


  return;
}
/******************************************************************************/
void make_length_params_dimensionless()
{
  int idim,iwall,iwall_type,jwall_type,icomp,jcomp,izone;

  for (idim=0; idim < Ndim; ++idim){
       Size_x[idim] /= Length_ref;
       Esize_x[idim] /= Length_ref;
  }

  for (iwall=0;iwall<Nwall;iwall++)
    for (idim=0; idim < Ndim; ++idim){ 
       WallPos[idim][iwall] /=Length_ref;
  }

  for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
     WallParam[iwall_type]/=Length_ref;
     WallParam_2[iwall_type]/=Length_ref;
     WallParam_3[iwall_type]/=Length_ref;
  }

  for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
     if (Mix_type=1){
         for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) Sigma_ww[iwall_type][jwall_type]/=Length_ref;
     }
     else  Sigma_w[iwall_type]/=Length_ref;
     for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) Cut_ww[iwall_type][jwall_type]/=Length_ref;
  }

  if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
      Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
      || Type_uwwPot==PAIR_r18andYUKAWA_CS){

     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
        if (Mix_type=1){
            for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) YukawaK_ww[iwall_type][jwall_type]*=Length_ref;
        }
        else  YukawaK_w[iwall_type]*=Length_ref;
     }
  }

  for (icomp=0;icomp<Ncomp;icomp++)
     for (jcomp=0;jcomp<Ncomp;jcomp++){
        Sigma_ff[icomp][jcomp]/=Length_ref;
        Bond_ff[icomp][jcomp]/=Length_ref;
        Cut_ff[icomp][jcomp]/=Length_ref;

        if (Type_pairPot == PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS ||
            Type_pairPot==PAIR_LJandYUKAWA_CS || Type_pairPot==PAIR_r12andYUKAWA_CS ||
            Type_pairPot==PAIR_r18andYUKAWA_CS || Type_pairPot==PAIR_rNandYUKAWA_CS) YukawaK_ff[icomp][jcomp]*=Length_ref;
     }

  if (Mix_type==1) 
  for (icomp=0;icomp<Ncomp;icomp++){
     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
       Sigma_wf[icomp][iwall_type]/=Length_ref; 
       Cut_wf[icomp][iwall_type]/=Length_ref; 
       YukawaK_wf[icomp][iwall_type]*=Length_ref; 
     }
  }

  if (Type_interface != UNIFORM_INTERFACE){
     X_const_mu /= Length_ref;
  }

  if (Type_coul!=NONE || Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
     for (iwall=0; iwall <  Nwall; ++iwall){
        if(Type_bc_elec[WallType[iwall]]==1) Elec_param_w[iwall]*=(Length_ref*Length_ref); 
     }    
     if (Type_dielec == DIELEC_WF_PORE) Dielec_X /= Length_ref;
  }
 
  if (Nzone !=0) Rmax_zone[izone] /= Length_ref;

  X_1D_bc/=Length_ref;

  return;
}
/******************************************************************************/
void make_energy_params_dimensionless()
{
  int icomp,jcomp,iwall_type,jwall_type,i,j;

  for (icomp=0;icomp<Ncomp;icomp++){
     for (jcomp=0;jcomp<Ncomp;jcomp++){
         Eps_ff[icomp][jcomp]/=Temp;

         if (Type_pairPot == PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS ||
             Type_pairPot==PAIR_LJandYUKAWA_CS || Type_pairPot==PAIR_r12andYUKAWA_CS ||
             Type_pairPot==PAIR_r18andYUKAWA_CS || Type_pairPot==PAIR_rNandYUKAWA_CS) EpsYukawa_ff[icomp][jcomp]/=Temp;
     }
  }

  for (iwall_type=0; iwall_type<Nwall_type; iwall_type++){
     if (Mix_type==1){
         for (jwall_type=0;iwall_type<Nwall_type;iwall_type++) Eps_ww[iwall_type][jwall_type]/=Temp;
     }
     else Eps_w[iwall_type]/=Temp;
  }
   
  if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
      Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
      || Type_uwwPot==PAIR_r18andYUKAWA_CS){

     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
        if (Mix_type=1){
            for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) EpsYukawa_ww[iwall_type][jwall_type]/=Temp;
        }
        else  EpsYukawa_w[iwall_type]/=Temp;
     }
  }

  if (Mix_type==1) 
  for (icomp=0;icomp<Ncomp;icomp++){
     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){ 
         Eps_wf[icomp][iwall_type]/=Temp;
         EpsYukawa_wf[icomp][iwall_type]/=Temp;
         Vext_membrane[iwall_type][icomp]/=Temp;
     }
  }

#ifdef USE_LOCA
  if (Loca.method !=-1){
    if (Loca.cont_type1 == CONT_TEMP) Loca.step_size/=(Eps_ff[0][0]*Temp);
    else if (Loca.cont_type1 == CONT_EPSW_I ||  Loca.cont_type1 == CONT_EPSW_ALL ||
             Loca.cont_type1 == CONT_EPSWF_IJ ||  Loca.cont_type1 == CONT_EPSWF_ALL ||
              Loca.cont_type1 == CONT_EPSFF_IJ || Loca.cont_type1 == CONT_EPSFF_ALL )Loca.step_size /= Temp;
  }
#endif
  
}
/******************************************************************************/
void make_elecPot_params_dimensionless()
{
  int iwall;

  if (Type_interface != UNIFORM_INTERFACE){
    Elec_pot_LBB/=Potential_ref;
    Elec_pot_RTF/=Potential_ref;
  }

   if (Type_pairPot==PAIR_COULOMB_CS || Type_pairPot==PAIR_COULOMB) {
     for (iwall=0; iwall <  Nwall; iwall++){
        if(Type_bc_elec[WallType[iwall]]==1) Elec_param_w[iwall]/=Potential_ref;
     }    
   }
   return;
}
/******************************************************************************/
void make_dielecConst_params_dimensionless()
{
   int iwall_type;

   if (Type_dielec != DIELEC_WF_PORE){
     Dielec_bulk /= Dielec_ref;
   }
   else{
     Dielec_bulk /= Dielec_ref;
     Dielec_pore /= Dielec_ref;
   }
   if (Nwall_type>0) {
      for (iwall_type=0;iwall_type<Nwall_type;iwall_type++) Dielec_wall[iwall_type]/=Dielec_ref;
   }

   return;
}
/******************************************************************************/
void make_density_params_dimensionless()
{
  int iwall_type,icomp,ipol_comp;

  for (iwall_type=0;iwall_type<Nwall_type;iwall_type++) Rho_w[iwall_type] /=Density_ref;

  if (Type_poly==NONE || Ntype_mer==1){
    for (icomp=0;icomp<Ncomp;icomp++) {
        if (Type_interface == UNIFORM_INTERFACE) Rho_b[icomp]/=Density_ref;
        else {
           Rho_b_LBB[icomp]/=Density_ref;
           Rho_b_RTF[icomp]/=Density_ref;
        }
    }
  }

#ifdef USE_LOCA
  if (Loca.method !=-1){
    if (Loca.cont_type1 == CONT_RHO_I || Loca.cont_type1 == CONT_LOG_RHO_I) Loca.step_size/Density_ref;
  }
#endif
  return;
}
/******************************************************************************/
/* fill_SurfGeom_struct : this function pushes all the data from the generic input
formats to the surface geometry structures that will be used to access data */
void fill_surfGeom_struct()
{
  int iw,idim,i;
  struct SurfaceGeom_Struct *sgeom_iw;
  double r,rsq;

  if (Nwall_type>0) Poly_graft_dist = (double *) array_alloc (1, Nwall_type, sizeof(double));

  for (iw=0;iw<Nwall_type;iw++){
    Poly_graft_dist[iw]=0.;
    sgeom_iw = &(SGeom[iw]); 
    sgeom_iw->surfaceTypeID=Surface_type[iw];
    sgeom_iw->Lperiodic_overlay=FALSE;
    sgeom_iw->Llinear_overlay=FALSE;
    sgeom_iw->Lwedge_cutout=FALSE;
    switch(sgeom_iw->surfaceTypeID)
    {
       case smooth_planar_wall:
            sgeom_iw->halfwidth = (double *) array_alloc(1, Ndim, sizeof(double));
            sgeom_iw->orientation=Orientation[iw];
            sgeom_iw->halfwidth[Orientation[iw]]=WallParam[iw];
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Lperiodic_overlay=Lperiodic_overlay[iw];
            sgeom_iw->Nperiodic_overlay=Nperiodic_overlay[iw];
            for (i=0;i<Nperiodic_overlay[iw]; i++){ 
               sgeom_iw->orientation_periodic[i]=OrientationPeriodicFunc[iw][i];
               sgeom_iw->amplitude[i]=AmplitudePeriodicFunc[iw][i];
               sgeom_iw->wavelength[i]=WavelengthPeriodicFunc[iw][i];
               sgeom_iw->origin_PeriodicFunc[i]=OriginPeriodicFunc[iw][i];
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
            }
            Poly_graft_dist[iw]=sgeom_iw->halfwidth[sgeom_iw->orientation];
            break;

       case finite_planar_wall:
            sgeom_iw->halfwidth = (double *) array_alloc(1, Ndim, sizeof(double));
            sgeom_iw->halfwidth[0]=WallParam[iw];
            if (Ndim>1) sgeom_iw->halfwidth[1]=WallParam_2[iw];
            if (Ndim>2) sgeom_iw->halfwidth[2]=WallParam_3[iw];
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
            }
            break;

       case point_surface:
            rsq=0.0;
            for (idim=0;idim<Ndim;idim++) rsq+=Esize_x[idim]*Esize_x[idim];
            r=sqrt(rsq)/2.0;
            sgeom_iw->radius=r;
            sgeom_iw->halfwidth = (double *) array_alloc(1, Ndim, sizeof(double));
            for(idim=0;idim<Ndim;idim++) sgeom_iw->halfwidth[idim]=0.5*Esize_x[idim];
            break;

       case colloids_cyl_sphere:
            sgeom_iw->radius=WallParam[iw];
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            Poly_graft_dist[iw]=sgeom_iw->radius;
            break;

       case finite_cyl_3D:
            sgeom_iw->orientation=Orientation[iw];
            sgeom_iw->radius=WallParam[iw];
            sgeom_iw->halflength=WallParam_2[iw];
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lperiodic_overlay=Lperiodic_overlay[iw];
            sgeom_iw->Nperiodic_overlay=Nperiodic_overlay[iw];
            for (i=0;i<Nperiodic_overlay[iw]; i++){ 
               sgeom_iw->orientation_periodic[i]=OrientationPeriodicFunc[iw][i];
               sgeom_iw->amplitude[i]=AmplitudePeriodicFunc[iw][i];
               sgeom_iw->wavelength[i]=WavelengthPeriodicFunc[iw][i];
               sgeom_iw->origin_PeriodicFunc[i]=OriginPeriodicFunc[iw][i];
               if (OrientationPeriodicFunc[iw][i] != Orientation[iw]){
                   if (Iwrite_screen != SCREEN_NONE) {
                   printf("Orientation of periodic function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylindrical surface.  Adjustments can only be made along length of cylinder\n");
                   printf("Resetting periodic orientation to match the surface orientation\n");
                   }
                   sgeom_iw->orientation_periodic[i]=Orientation[iw];
               }
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
               if (OrientationLinearFunc[iw][i] != Orientation[iw]){
                   if (Iwrite_screen != SCREEN_NONE) {
                   printf("Orientation of linear function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylinder.  Adjustments can only be made along length of the surface.\n");
                   printf("Resetting linear orientation to match the surface orientation\n");
                   }
                   sgeom_iw->orientation_linear[i]=Orientation[iw];
               }
            }
            Poly_graft_dist[iw]=sgeom_iw->radius;
            break;

       case atomic_centers:
            break;

       case cyl2D_sphere3D_pore:
            sgeom_iw->radius=WallParam[iw];
            Poly_graft_dist[iw]=sgeom_iw->radius;
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            break;

       case cyl3D_slit2D_pore:
            sgeom_iw->orientation=Orientation[iw];
            sgeom_iw->radius=WallParam[iw];
            sgeom_iw->halflength=WallParam_2[iw];
            sgeom_iw->Lrough_surface=Lrough_surf[iw];
            if (Lrough_surf[iw]==TRUE){
               sgeom_iw->roughness=Rough_param_max[iw];
               sgeom_iw->roughness_length=Rough_length[iw];
            }
            sgeom_iw->Lperiodic_overlay=Lperiodic_overlay[iw];
            sgeom_iw->Nperiodic_overlay=Nperiodic_overlay[iw];
            for (i=0;i<Nperiodic_overlay[iw]; i++){ 
               sgeom_iw->orientation_periodic[i]=OrientationPeriodicFunc[iw][i];
               sgeom_iw->amplitude[i]=AmplitudePeriodicFunc[iw][i];
               sgeom_iw->wavelength[i]=WavelengthPeriodicFunc[iw][i];
               sgeom_iw->origin_PeriodicFunc[i]=OriginPeriodicFunc[iw][i];
               if (OrientationPeriodicFunc[iw][i] != Orientation[iw]){
                   if (Iwrite_screen != SCREEN_NONE) {
                   printf("Orientation of periodic function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylindrical or 2D slit pore.  adjustments can only be made along length of pore\n");
                   printf("Resetting periodic orientation to match the surface orientation\n");
                   }
                   sgeom_iw->orientation_periodic[i]=Orientation[iw];
               }
            }
            sgeom_iw->Lwedge_cutout=Lwedge_cutout[iw];
            if(Lwedge_cutout[iw]==TRUE){
               sgeom_iw->angle_wedge_start=Angle_wedge_start[iw];
               sgeom_iw->angle_wedge_end=Angle_wedge_end[iw];
            }
            sgeom_iw->Llinear_overlay=Llinear_overlay[iw];
            sgeom_iw->Nlinear_overlay=Nlinear_overlay[iw];
            for (i=0;i<Nlinear_overlay[iw]; i++){ 
               sgeom_iw->orientation_linear[i]=OrientationLinearFunc[iw][i];
               sgeom_iw->slope[i]=SlopeLinearFunc[iw][i];
               sgeom_iw->origin_LinearFunc[i]=OriginLinearFunc[iw][i];
               sgeom_iw->endpoint_LinearFunc[i]=EndpointLinearFunc[iw][i];
               if (OrientationLinearFunc[iw][i] != Orientation[iw]){
                   if (Iwrite_screen != SCREEN_NONE) {
                   printf("Orientation of linear function must be the same as the orientation of the surface\n");
                   printf("for the 3D cylindrical or 2D slit pore.  adjustments can only be made along length of pore\n");
                   printf("Resetting linear orientation to match the surface orientation\n");
                   }
                   sgeom_iw->orientation_linear[i]=Orientation[iw];
               }
            }
            Poly_graft_dist[iw]=sgeom_iw->radius;
            break;

       default: /* No surface found */
            if (Iwrite_screen != SCREEN_NONE)  printf("error with surface type iwall_type=%d not identified\n",iw);
            exit(-1);
            break;
    }
  }
  return;
}
/****************************************************************************************************************/
