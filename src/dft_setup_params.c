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
  char *tmp_string;
  char tmp_string_array[FILENAME_LENGTH];

  Open_GUI=FALSE;
  Read_OLDInput_File=TRUE;
  Read_XMLInput_File=FALSE;


  /****************************************************/
  /* first do all the things that only Proc 0 will do */
  /****************************************************/
  if (Proc==0) {

            /* first set basic selections for starting a problem */
#ifdef BUILD_GUI
        dft_OptikaGUI_control();

        strcpy(Runpath_array, Runpath);
        Runpath=Runpath_array;

        strcpy(Outpath_array, Outpath);
        Outpath=Outpath_array;

         /* add path to filename for this case because it is opened and closed several times */
         /* note that for all other cases, the pathname is added to the filename at the point where fopen is called. */
        strcpy(tmp_string_array,Outpath_array);
        strcpy(EchoInputFile_array,strcat(tmp_string_array,file_echoinput));
        file_echoinput=EchoInputFile_array;

        fpecho  = fopen(file_echoinput,"w+");

        if (Read_OLDInput_File==TRUE){
            strcpy(InputOLDFile_array,InputOLD_File);
            InputOLD_File=InputOLDFile_array;
            input_file=InputOLD_File;
            if( (fpinput  = fopen(input_file,"r")) == NULL) {
                printf("Can't open file %s\n", input_file);
                printf("Must enter all params from the GUI\n");
                Open_GUI=TRUE;
                input_file_exists=FALSE;
                Set_GUIDefaults_to_OLD_File=FALSE;
            }
            else{
               input_file_exists=TRUE;
               Set_GUIDefaults_to_OLD_File=TRUE;
            }
            if(input_file_exists==TRUE) read_input_file(fpinput,fpecho);
        }

        if (Open_GUI==TRUE){
             dft_OptikaGUI();

             /* set up the Density File to be the defined GUE value */
             if (Restart != NORESTART) {
               strcpy(DensityFile_array,DensityFile);
               DensityFile=DensityFile_array;

               if (Lbinodal==TRUE){
                  strcpy(DensityFile2_array,DensityFile2);
                  DensityFile2=DensityFile2_array; 
               }
             }

             if (Type_poly != NONE && Type_poly_arch==POLY_ARCH_FILE){
                 strcpy(poly_file_array,Poly_file_name);
                 Poly_file_name=poly_file_array; 
             }
             if (Type_poly == CMS || Type_poly==CMS_SCFT){
                 strcpy(cr_file_array,Cr_file);
                 Cr_file=cr_file_array; 
                 if (Ncr_files==2){
                    strcpy(cr_file2_array,Cr_file2);
                    Cr_file2=cr_file2_array; 
                 }
             }

             /* set up the chain architecture based on selections in GUI */
             if (Type_poly != NONE) setup_chain_architecture(Poly_file_name,fpecho);
             
        }
        else{
           if (Temp>0. && Type_coul!=NONE) Flag_mV_elecpot=TRUE;
           else         Flag_mV_elecpot=FALSE;

           if (Restart != NORESTART) {
              strcpy(tmp_string_array,Runpath_array);
              strcpy(DensityFile_array,strcat(tmp_string_array,"dft_dens.dat"));
              DensityFile=DensityFile_array;

              if (Lbinodal==TRUE){
                 strcpy(tmp_string_array,Runpath_array);
                 strcpy(DensityFile2_array,strcat(tmp_string_array,"dft_dens2.dat"));
                 DensityFile2=DensityFile2_array;
              }
           }

           if (Type_poly != NONE && Type_poly_arch==POLY_ARCH_FILE){
                 strcpy(tmp_string_array,Runpath_array);
                 strcpy(poly_file_array,strcat(tmp_string_array,Poly_file_name));
                 Poly_file_name=poly_file_array; 
           }
           if (Type_poly == CMS || Type_poly==CMS_SCFT){
                 strcpy(tmp_string_array,Runpath_array);
                 strcpy(cr_file_array,strcat(tmp_string_array,Cr_file));
                 Cr_file=cr_file_array; 
                 if (Ncr_files==2){
                    strcpy(tmp_string_array,Runpath_array);
                    strcpy(cr_file2_array,strcat(tmp_string_array,Cr_file2));
                    Cr_file2=cr_file2_array; 
                 }
           }

        }


#else
        Open_GUI=FALSE;

        strcpy(Runpath_array,"./");
        Runpath=Runpath_array;

        strcpy(tmp_string_array,Runpath_array);
        strcpy(EchoInputFile_array,strcat(tmp_string_array,file_echoinput));
        file_echoinput=EchoInputFile_array;

        fpecho  = fopen(file_echoinput,"w+");

        strcpy(tmp_string_array,Runpath_array);
        if( (fpinput  = fopen(strcat(tmp_string_array,input_file),"r")) == NULL) {
           printf("Can't open file %s\n", strcat(Runpath_array,input_file));
           exit(-1);
        }

        read_input_file(fpinput,fpecho);

        if (Temp>0. && Type_coul!=NONE) Flag_mV_elecpot=TRUE;
        else         Flag_mV_elecpot=FALSE;
        
        strcpy(tmp_string_array,Runpath_array);
        sprintf(wallPos_file_array, strcat(tmp_string_array,"dft_surfaces.dat"));
        WallPos_file_name=wallPos_file_array;

        if (Restart != NORESTART) {
           strcpy(tmp_string_array,Runpath_array);
           strcpy(DensityFile_array,strcat(tmp_string_array,"dft_dens.dat"));
           DensityFile=DensityFile_array;

           if(Lbinodal==TRUE){
              strcpy(tmp_string_array,Runpath_array);
              strcpy(DensityFile2_array,strcat(tmp_string_array,"dft_dens2.dat"));
              DensityFile2=DensityFile2_array;
           }
        }
        if (Type_poly != NONE && Type_poly_arch==POLY_ARCH_FILE){
              strcpy(tmp_string_array,Runpath_array);
              strcpy(poly_file_array,strcat(tmp_string_array,Poly_file_name));
              Poly_file_name=poly_file_array; 
        }
        if (Type_poly == CMS || Type_poly==CMS_SCFT){
                 strcpy(tmp_string_array,Runpath_array);
                 strcpy(cr_file_array,strcat(tmp_string_array,Cr_file));
                 Cr_file=cr_file_array; 
                 if (Ncr_files==2){
                    strcpy(tmp_string_array,Runpath_array);
                    strcpy(cr_file2_array,strcat(tmp_string_array,Cr_file2));
                    Cr_file2=cr_file2_array; 
                 }
         }

        /* Outpath=Runpath */
        strcpy(Outpath_array,Runpath_array);
        Outpath=Outpath_array;
#endif


    if (Nwall >0 && WallPos_file_name!=NULL) readIn_wall_positions_and_charges(WallPos_file_name,fpecho);
    if (Nwall >0 && Lrandom_walls==TRUE) setup_random_wall_positions(fpecho);
    if (Lauto_center == TRUE ) center_the_surfaces(fpecho);
    if (Lauto_size == TRUE  && Mesh_coarsening==2) autosize_the_domain(fpecho);
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
  double kappa,kappa_sq,lambda_B,lambda_D,lambda_GC,charge_sum;
  int test_CC;

  Lmesh_refine=FALSE;    /* used for automated mesh refinement - not currently used */

                         /****************************************************************/
                         /* set up Linked List for surfaces that are bonded or connected */
                         /****************************************************************/
  Link_list = (int **) array_alloc (2, Nlink,Nwall,sizeof(int));
  Nwall_this_link = (int *) array_alloc (1, Nlink,sizeof(int));
  for (i=0; i<Nlink; i++) Nwall_this_link[i]=0;
  for (iwall=0; iwall<Nwall; iwall++){Link_list[Link[iwall]][Nwall_this_link[Link[iwall]]++]=iwall; } 
  irand_range = 2147483647;

                /*************************************************************************/
                /* For Rough Surfaces.....populate a roughness array with random numbers */
                /*   that can be used to generate roughness profiles */
                /*************************************************************************/
  
  if (read_rough){
/*     #ifndef _MSC_VER
       srandom(135649);
     #else
       srand(135649);
     #endif*/
       srand(135649);

     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
       for (iblock=0;iblock<MAX_ROUGH_BLOCK;iblock++){
         for (jblock=0;jblock<MAX_ROUGH_BLOCK;jblock++){
            if (Proc==0){
/*            #ifndef _MSC_VER
              irand = random();
            #else
              irand = rand();
            #endif*/
            irand = rand() % irand_range;
            }
            MPI_Bcast(&irand,1,MPI_INT,0,MPI_COMM_WORLD);  /* need this for parallel jobs to be identical to serial jobs */
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
      if (Proc==0 && (Iwrite_screen!=SCREEN_NONE || Iwrite_screen!=SCREEN_ERRORS_ONLY))  {
          printf("------------------------------------------------------------------------\n");
          printf("Electrostatics parameters .......\n\n");
          printf("\t Electrostatics are set up in standard units rather than Gaussian units\n");
          printf("\t Plasma parameter settings: \n \t\t d=%g Angstroms \n \t\t T=%g K \n \t\t epsilon_s=%g\n",
                    Sigma_Angstroms_plasma,Temp_K_plasma,DielecConst_plasma);

          printf("\t The computed reduced temperature is Telec = 4pi(kT)(epsilon_s epsilon_0)d/e^2 = %9.6f\n", Temp_elec);
          printf("\t (Note: A plasma parameter Gamma=(1/Telec) is defined in some literature)\n");
          printf("\n \t The Bjerrum length:\n");
          printf("\t\t lambda_B = d/Telec = %9.6f Angstroms\n", Sigma_Angstroms_plasma/Temp_elec);
          printf("\t\t lambda_B/d = 1/Telec = %9.6f\n", 1/Temp_elec);
          lambda_B=1./Temp_elec;
          kappa_sq = 0.0;
             for(icomp = 0; icomp<Ncomp; icomp++)
                  kappa_sq += (4.0*PI/Temp_elec)*Rho_b[icomp]*Charge_f[icomp]*Charge_f[icomp];
              kappa = sqrt(kappa_sq);
              lambda_D=1.0/kappa;
          printf("\n \t Debye Length:\n");
          printf("\t\t lambda_D=%9.6f Angstroms\n",Sigma_Angstroms_plasma/kappa);
          printf("\t\t lambda_D/d=%9.6f\n",lambda_D);
          test_CC=TRUE;
          for (iwall_type=0; iwall_type<Nwall_type; ++iwall_type) if (Type_bc_elec[iwall_type]==CONST_POTENTIAL) test_CC=FALSE;
          if (test_CC==TRUE){
               charge_sum=0.0;
               for (iwall=0; iwall<Nwall; ++iwall) charge_sum+=Elec_param_w[iwall];
               printf("\n \t Constant charge surfaces detected....\n");
               printf("\t Total charge expeced in fluid = %9.6f\n",charge_sum);

               printf("\n\t Gouy Chapman length for each surface interacting with an ion of valence 1 is: \n");
               for (iwall=0; iwall<Nwall; ++iwall){ 
                    lambda_GC=1./(2*PI*lambda_B*fabs(Elec_param_w[WallType[iwall]]));
                    printf("\t\t iwall=%d: \n \t\t\t lambda_GC=%9.6f Angstroms \n\t\t\t lambda_GC/d=%9.6f \n\t\t\t lambda_GC/lambda_D=%9.6f\n",
                       iwall,Sigma_Angstroms_plasma*lambda_GC,lambda_GC,lambda_GC/lambda_D);
               }
               printf("\n \t lambda_B/lambda_GC=\%9.6f\n",lambda_B/lambda_GC);
               printf("\t (Recall that weak coupling has lambda_B/lambda_GC<1)\n");

          }
          else printf("\n \t Constant potential surfaces detected....surface charge is unknown\n");
          printf("------------------------------------------------------------------------\n");
      }
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

  if (Type_poly!=NONE && Physics_scaling != FALSE) Lprint_scaleFacWJDC=TRUE;
  else Lprint_scaleFacWJDC=FALSE;

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
     for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) Sigma_ww[iwall_type][jwall_type]/=Length_ref;
     for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) Cut_ww[iwall_type][jwall_type]/=Length_ref;
  }

  if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
      Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
      || Type_uwwPot==PAIR_r18andYUKAWA_CS){

     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++){
        for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) YukawaK_ww[iwall_type][jwall_type]*=Length_ref;
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
     for (jwall_type=0;iwall_type<Nwall_type;iwall_type++) Eps_ww[iwall_type][jwall_type]/=Temp;
  }
   
  if (Type_uwwPot==PAIR_YUKAWA_CS || Type_uwwPot == PAIR_EXP_CS || 
      Type_uwwPot==PAIR_LJandYUKAWA_CS || Type_uwwPot==PAIR_r12andYUKAWA_CS
      || Type_uwwPot==PAIR_r18andYUKAWA_CS){

     for (iwall_type=0;iwall_type<Nwall_type;iwall_type++)
            for (jwall_type=0;jwall_type<Nwall_type;jwall_type++) EpsYukawa_ww[iwall_type][jwall_type]/=Temp;
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
  int iwall_type,icomp,ipol_comp,ncomp_tot;

  if (Type_poly != NONE) ncomp_tot=Npol_comp;
  else ncomp_tot=Ncomp;
  

  for (iwall_type=0;iwall_type<Nwall_type;iwall_type++) Rho_w[iwall_type] /=Density_ref;

  for (icomp=0;icomp<ncomp_tot;icomp++) {
        if (Type_interface == UNIFORM_INTERFACE) Rho_b[icomp]/=Density_ref;
        else {
           Rho_b_LBB[icomp]/=Density_ref;
           Rho_b_RTF[icomp]/=Density_ref;
        }
  }

#ifdef USE_LOCA
  if (Loca.method !=-1){
    if (Loca.cont_type1 == CONT_RHO_I || Loca.cont_type1 == CONT_LOG_RHO_I) Loca.step_size/=Density_ref;
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
printf("iw=%d roughness=%9.4f roughness_length=%9.4f",iw,sgeom_iw->roughness,sgeom_iw->roughness_length);
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
void readIn_wall_positions_and_charges(char *WallPos_file_name,FILE *fpecho)
{
   int iwall,idim,dim_tmp,nlink_chk,new_wall,jwall;
   FILE *fpsurfaces;

    if( (fpsurfaces  = fopen(WallPos_file_name,"r")) == NULL) {
        printf("Can't open file %s and no other file was specified in GUI.\n",WallPos_file_name); exit(-1);
    }
    else{ printf("Read in surface positions and charges from %s\n",WallPos_file_name); }

    for (iwall=0; iwall<Nwall; iwall++){
       fscanf(fpsurfaces,"%d  %d",&WallType[iwall], &Link[iwall]);
       for (idim=0; idim<Ndim; idim++) {
          dim_tmp=idim;
                      /* temporary rotation of coordinates
                                if (idim==0) dim_tmp=1;
                                else if (idim==1) dim_tmp=2;
                                else if (idim==2) dim_tmp=0;*/
                      /* end of temporary code */
          fscanf(fpsurfaces,"%lf",&WallPos[dim_tmp][iwall]);
       }
       fscanf(fpsurfaces,"%lf",&Elec_param_w[iwall]);
     }
     fclose(fpsurfaces);

     nlink_chk = 1;
     for (iwall=1; iwall<Nwall; iwall++){
       new_wall = TRUE;
       for (jwall=0; jwall<iwall; jwall++)
        if (Link[iwall] == Link[jwall]) new_wall=FALSE;
       if (new_wall) nlink_chk++;
     }
     if (nlink_chk != Nlink){
       printf("Check Nlink in dft_input.dat: %d and assignments in %s: %d\n",
             Nlink,WallPos_file_name,nlink_chk);
      exit(-1);
     }


    for (iwall=0; iwall<Nwall; iwall++){
          fprintf(fpecho,"\niwall=%d  WallType=%d  Link=%d Elec_param_w=%g",iwall,WallType[iwall],Link[iwall],Elec_param_w[iwall]);
       for (idim=0; idim<Ndim; idim++)  fprintf(fpecho,"   WallPos[idim=%d]=%g\n",idim,WallPos[idim][iwall]);
    }
        
  return;
}
/****************************************************************************************************/
void setup_random_wall_positions(FILE *fpecho){
                                         /* below is simple code for random placement of surface */
                                         /* coordinates. This could be modified to disallow */
                                         /* overlaps or to apply it only to selected surfaces */
                                         /* to be changed to a proper selection in the input file*/
int irand,irand_range,iwall,idim,dim_tmp;

  irand_range = 2147483647;

/*  #ifndef _MSC_VER
    srandom(135649);
  #else
    srand(135649);
  #endif*/
    srand(135649);

 for (iwall=0; iwall<Nwall; iwall++){
    for (idim=0; idim<Ndim; idim++) {
       dim_tmp=idim;
                      /* temporary rotation of coordinates
                                if (idim==0) dim_tmp=1;
                                else if (idim==1) dim_tmp=2;
                                else if (idim==2) dim_tmp=0;*/
                      /* end of temporary code */
       if (Lrandom_walls==TRUE || fabs(WallPos[dim_tmp][iwall]+9999.0)<1.e-6) {
/*            #ifndef _MSC_VER
              irand = random();
            #else
              irand = rand();
            #endif*/
            irand = rand() % irand_range;
            WallPos[dim_tmp][iwall] = Size_x[idim]*(-0.5+( ((double)irand)/((double)irand_range)));
            fprintf(fpecho,"\n Wall %d dim %d gets WallPos:%g \n",iwall,idim,WallPos[idim][iwall]);
          }  
     }
  }
  return;
}
/****************************************************************************************************/
void center_the_surfaces(FILE *fpecho)
{
  int iwall,idim,dim_tmp;
  double maxpos[3],minpos[3];

  for (idim=0; idim<Ndim; idim++){ minpos[idim] = 1000.; maxpos[idim]=-1000.;}
  for (iwall=0;iwall<Nwall;iwall++){
     for (idim=0;idim<Ndim;idim++){
           dim_tmp=idim;
           if (WallPos[dim_tmp][iwall] < minpos[dim_tmp]) minpos[dim_tmp]=WallPos[dim_tmp][iwall];
           if (WallPos[dim_tmp][iwall] > maxpos[dim_tmp]) maxpos[dim_tmp]=WallPos[dim_tmp][iwall];
     }
  }

  for (iwall=0; iwall<Nwall; iwall++) for (idim=0; idim<Ndim; idim++) {
       WallPos[idim][iwall] -= 0.5*(maxpos[idim] + minpos[idim]);
       fprintf(fpecho,"iwall=%d  idim=%d  (centered)WallPos=%g\n",iwall,idim,WallPos[idim][iwall]);
  }

return;
}
/****************************************************************************************************/
void autosize_the_domain(FILE *fpecho)
{
  int iwall,idim,dim_tmp;
  double maxpos[3],minpos[3];

  for (idim=0; idim<Ndim; idim++){ minpos[idim] = 1000.; maxpos[idim]=-1000.;}
  for (iwall=0;iwall<Nwall;iwall++){
     for (idim=0;idim<Ndim;idim++){
           dim_tmp=idim;
           if (WallPos[dim_tmp][iwall] < minpos[dim_tmp]) minpos[dim_tmp]=WallPos[dim_tmp][iwall];
           if (WallPos[dim_tmp][iwall] > maxpos[dim_tmp]) maxpos[dim_tmp]=WallPos[dim_tmp][iwall];
     }
  }

  for (idim=0; idim<Ndim; idim++){
      Size_x[idim] = maxpos[idim]-minpos[idim] + 2.0*(Rmax_zone[0]+Sigma_ff[0][0]); 
      fprintf(fpecho,"idim=%d  (auto)Size_x=%g\n",idim,Size_x[idim]);
  }

  return;
}
/****************************************************************************************************/
