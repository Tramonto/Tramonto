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
 *  FILE: dft_guess.c
 *
 *  This file contains high level logic for setting up an initial 
 *  guess for the solution vector.
 *
 */

#include "dft_guess.h"
 
void set_initial_guess (int iguess, double** xOwned)
{
  double t1=0.0;
  int i;
  int start_no_info;
  double **xInBox;

  if (Proc==0) { 
    /*  printf("\n%s: Setting an initial guess ... ",yo);*/
      t1 = MPI_Wtime();
  }

  /* temporary array set to size of box node to allow us to compute some variables
     that require integrations over densities */ 
  xInBox = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));

  if (Proc==0 && Iwrite==VERBOSE){
      if (Restart==0 && Imain_loop==0) printf("generating guess from scratch\n");  
      else printf("setting up guess from existing files\n");
  } 

  /* Get some or all parts of unknown vector from previously generated restart files.
     Note that these files do not need to have all of the solution vector and can 
     be patched up in some cases. */

  for (i=0;i<NEQ_TYPE;i++) Restart_field[i]=FALSE;
  if (Restart > 0 || Imain_loop > 0){
       start_no_info = FALSE;
       guess_restart_from_files(start_no_info,iguess,xInBox);  
  } 
  else start_no_info = TRUE;

                          /* this loop sets up portions of the solution vector that were not treated using 
                             restart files.  This can be part or all of the solution vector. */
  for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY:
           if (Phys2Nunk[DENSITY]>0 && start_no_info){
               if (Type_poly == CMS || Type_poly ==CMS_SCFT){
                   setup_polymer_rho(xInBox,iguess);
               }
               else{
                   setup_density(xInBox,iguess);
               }
           }
           else if (Phys2Nunk[DENSITY]>0 && Restart_field[DENSITY]==FALSE){
                printf("we don't have the ability to restart without a density field at this time\n");
                exit(-1);
           }
           break;
         case MF_EQ:
           if (Phys2Nunk[MF_EQ]>0 && (start_no_info || Restart_field[MF_EQ]==FALSE ||Restart==3)){
               calc_init_mf_attract(xInBox); 
               /*setup_mf_attract(xInBox); */
           } break;
         case HSRHOBAR:
           if (Phys2Nunk[HSRHOBAR]>0 && (start_no_info || Restart_field[HSRHOBAR]==FALSE ||Restart==3)) {
                 calc_init_rho_bar(xInBox);
                 /*setup_rho_bar(xInBox);*/ /* bulk setup disabled for now */
           } break;
         case POISSON:
           if (Phys2Nunk[POISSON]>0 && (start_no_info || Restart_field[POISSON]==FALSE ||Restart==3)) setup_elec_pot(xInBox,iguess); break;
         case DIFFUSION: 
           if (Phys2Nunk[DIFFUSION]>0 && (start_no_info || Restart_field[DIFFUSION]==FALSE ||Restart==3)) setup_chem_pot(xInBox); break;

         case CAVWTC:
           if (Phys2Nunk[CAVWTC]>0 && (start_no_info || Restart_field[CAVWTC]==FALSE)){
              calc_init_Xi_cavWTC(xInBox);
              /*setup_Xi_cavWTC(xInBox);*/
           }  break;

         case BONDWTC:
           if (Phys2Nunk[BONDWTC]>0 && (start_no_info || Restart_field[BONDWTC]==FALSE)){
              calc_init_BondWTC(xInBox);
              /*setup_BondWTC(xInBox);*/
           } break;

         case WJDC_FIELD:
           if (Phys2Nunk[WJDC_FIELD]>0 && (start_no_info || Restart_field[WJDC_FIELD]==FALSE)){
                 calc_init_WJDC_field(xInBox);
                 /*setup_polymer_field_wjdc(xInBox); */
           } break;
         case CMS_FIELD:
           if (Phys2Nunk[CMS_FIELD]>0 && (start_no_info || Restart_field[CMS_FIELD]==FALSE)) {
              calc_init_CMSfield(xInBox);
              communicate_to_fill_in_box_values(xInBox);
              /*setup_polymer_field(xInBox,iguess); */
           } break;
         case G_CHAIN:
           if (Phys2Nunk[G_CHAIN]>0 && (start_no_info || Restart_field[G_CHAIN]==FALSE)){
                if (Type_poly==CMS) setup_polymer_G(xInBox); 
                else if (Type_poly==WJDC) setup_polymer_G_wjdc(xInBox);
           }
           break;
     case YW_DENS:
       break;
         default:
           printf("problem with switch in initial guess\n");
           exit(-1);
           break;
     }
  }
  if (Restart == 2) chop_profile(xInBox,iguess);  /* special case for treating wetting problems */

  check_zero_densities(xInBox);              
  translate_xInBox_to_xOwned(xInBox,xOwned);
  safe_free((void **) &xInBox);

/* note need to put the setup_poymer_simple functionality inside the various setup_polymer_CMS routines....*/

  if (Proc==0 && Iwrite==VERBOSE) printf("\n initial guess took %g secs\n",MPI_Wtime()-t1);
  return;
}
/********************************************************************************************************/
void translate_xInBox_to_xOwned(double **xInBox,double **xOwned)
{
   int loc_inode,iunk;
   for (loc_inode=0; loc_inode<Nnodes_per_proc;loc_inode++)
      for (iunk=0;iunk<Nunk_per_node;iunk++){ 
         if (L2B_node[loc_inode]>=0) xOwned[iunk][loc_inode]=xInBox[iunk][L2B_node[loc_inode]];
      }
   return;
}
/********************************************************************************************************/
/*void communicate_to_fill_in_box_values: this routine is needed for CMS and WJDC functionals where the
  field variables can only be computed on local nodes, but need to be known on all of the box coordinates
  in order to enable computation of the chain variables on the local nodes*/
void communicate_to_fill_in_box_values(double **xInBox)
{

  if (Proc == 0) {
    X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
  }
  collect_x_old(xInBox);
  communicate_profile(X_old,xInBox);
  safe_free((void *) &X_old);
  return;
}
/********************************************************************************************************/
