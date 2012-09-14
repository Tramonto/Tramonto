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
 
void set_initial_guess (int guess_type, double** xOwned)
{
  double t1=0.0;
  int i,iunk,inode_box;
  int start_no_info;
  double **xInBox;

/*  if (Proc==0 && Iwrite_screen !=SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY) { 
      printf("Set up Initial Guess ... ");
      t1 = MPI_Wtime();
  }*/

  /* temporary array set to size of box node to allow us to compute some variables
     that require integrations over densities */ 
  if (Type_poly==WJDC3 && Grafted_Logical) xInBox = (double **) array_alloc(2, Nunk_per_node, Nnodes_box_extra, sizeof(double));
  else                                     xInBox = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));
 
  for (iunk=0;iunk< Nunk_per_node; iunk++)
   for (inode_box=0;inode_box<Nnodes_box; inode_box++) xInBox[iunk][inode_box]=0.0;

/*  if (Proc==0 && Iwrite_screen !=SCREEN_NONE && Iwrite_screen != SCREEN_ERRORS_ONLY){
      if (Restart==NORESTART && Imain_loop==0 ) printf("start from scratch\n");  
      else                                      printf("from existing files\n");
  } */

  /* Get some or all parts of unknown vector from previously generated restart files.
     Note that these files do not need to have all of the solution vector and can 
     be patched up in some cases. */

  for (i=0;i<NEQ_TYPE;i++) Restart_field[i]=FALSE;
/*  if (Restart != NORESTART || (Imain_loop > 0 && (Nwall!=0 || LBulk==TRUE)) ){*/
  if (Restart != NORESTART || Imain_loop > 0  ){
       start_no_info = FALSE;
       guess_restart_from_files(start_no_info,guess_type,xInBox);  
       translate_xInBox_to_xOwned(xInBox,xOwned);
  } 
  else start_no_info = TRUE;

                          /* this loop sets up portions of the solution vector that were not treated using 
                             restart files.  This can be part or all of the solution vector. */
  for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY:
           if (Phys2Nunk[DENSITY]>0 && (start_no_info || Restart==RESTART_FEWERCOMP)){
               if (Type_poly == CMS || Type_poly ==CMS_SCFT){
                   setup_polymer_rho(xInBox,xOwned,guess_type);
               }
               else{
                   setup_density(xInBox,xOwned,guess_type);
               }
           }
           else if (Phys2Nunk[DENSITY]>0 && Restart_field[DENSITY]==FALSE){
                if (Iwrite_screen != SCREEN_NONE) printf("we don't have the ability to restart without a density field at this time\n");
                exit(-1);
           }
           break;
         case MF_EQ:
           if (Phys2Nunk[MF_EQ]>0 && (start_no_info || Restart_field[MF_EQ]==FALSE ||Restart==RESTART_DENSONLY || Restart==RESTART_FEWERCOMP)){
               if (Iguess_fields!=BULK)  calc_init_mf_attract(xInBox,xOwned); 
               else                      setup_mf_attract(xOwned); 
           } break;
         case HSRHOBAR:
           if (Phys2Nunk[HSRHOBAR]>0 && (start_no_info || Restart_field[HSRHOBAR]==FALSE ||Restart==RESTART_DENSONLY)) {
                 if (Iguess_fields!=BULK) calc_init_rho_bar(xInBox,xOwned);
                 else                     setup_rho_bar(xOwned);
           } break;
         case POISSON:
           if (Phys2Nunk[POISSON]>0 && (start_no_info || Restart_field[POISSON]==FALSE ||Restart==RESTART_DENSONLY)) {
                 if (Iguess_fields==CALC_ALL_FIELDS){  
                      if (start_no_info==TRUE) setup_elec_pot(xOwned,guess_type);
                      calc_init_elec_pot(xInBox,xOwned);
                 }
                 else setup_elec_pot(xOwned,guess_type); 
           }
           break;
         case DIFFUSION: 
           if (Phys2Nunk[DIFFUSION]>0 && (start_no_info || Restart_field[DIFFUSION]==FALSE ||Restart==RESTART_DENSONLY || Restart==RESTART_FEWERCOMP)) {
                 if (Iguess_fields==CALC_ALL_FIELDS){  
                    if (start_no_info==TRUE) setup_chem_pot(xOwned); 
                    calc_init_chem_pot(xInBox,xOwned);
                 }
                 else  setup_chem_pot(xOwned); 
           }
           break;

         case CAVWTC:
           if (Phys2Nunk[CAVWTC]>0 && (start_no_info || Restart_field[CAVWTC]==FALSE)){
              if (Iguess_fields!=BULK) calc_init_Xi_cavWTC(xInBox,xOwned);
              else                     setup_Xi_cavWTC(xOwned);
           }  break;

         case BONDWTC:
           if (Phys2Nunk[BONDWTC]>0 && (start_no_info || Restart_field[BONDWTC]==FALSE)){
              if (Iguess_fields!=BULK) calc_init_BondWTC(xInBox,xOwned);
              else                     setup_BondWTC(xOwned);
           } break;

         case WJDC_FIELD:
           if (Phys2Nunk[WJDC_FIELD]>0 && (start_no_info || Restart_field[WJDC_FIELD]==FALSE || Restart==RESTART_FEWERCOMP)){
                 if (Iguess_fields==CALC_ALL_FIELDS)  calc_init_WJDC_field(xInBox,xOwned);
                 else                                 setup_polymer_field_wjdc(xOwned); 
           } 
           break;

         case CMS_FIELD:
           if (Phys2Nunk[CMS_FIELD]>0 && (start_no_info || Restart_field[CMS_FIELD]==FALSE || Restart==RESTART_FEWERCOMP)) {
              if (Iguess_fields==CALC_ALL_FIELDS) calc_init_CMSfield(xInBox,xOwned);
              else setup_polymer_field(xInBox,xOwned,guess_type); 
           } 
           break;
         case SCF_FIELD:
           if (Phys2Nunk[SCF_FIELD]>0 && (start_no_info || Restart_field[SCF_FIELD]==FALSE || Restart==RESTART_FEWERCOMP)) {
               if (Iguess_fields==CALC_ALL_FIELDS) calc_init_SCFfield(xInBox,xOwned);
               else setup_polymer_SCF_field(xInBox,xOwned,guess_type); 
           } 
           break;
           case SCF_CONSTR:
           if (Phys2Nunk[SCF_CONSTR]>0 && (start_no_info || Restart_field[SCF_CONSTR]==FALSE)) calc_init_lambda(xInBox,xOwned);
           break; 
         case G_CHAIN:
           if (Phys2Nunk[G_CHAIN]>0 && (start_no_info || Restart_field[G_CHAIN]==FALSE)){
                if (Type_poly==CMS){
                    if (Iguess_fields == CALC_ALL_FIELDS || Iguess_fields == CALC_RHOBAR_AND_G) {
                             calc_init_polymer_G_CMS(xInBox,xOwned);
                    }
                    else     setup_polymer_G(xInBox,xOwned); 
                }
                else if (Type_poly==CMS_SCFT) calc_init_polymer_G_SCF(xInBox,xOwned);
                else if (Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3){
                   if (Iguess_fields == CALC_ALL_FIELDS || Iguess_fields == CALC_RHOBAR_AND_G) {
                          calc_init_polymer_G_wjdc(xInBox,xOwned);
                   }
                   else   setup_polymer_G_wjdc(xOwned);
                }
           }
           break;

	 default:
           if (Iwrite_screen != SCREEN_NONE) printf("problem with switch in initial guess\n");
           exit(-1);
           break;
     }
  }
  if (Restart == RESTART_STEP) chop_profile(xInBox,guess_type);  /* special case for treating wetting problems */
  if (Type_poly==WJDC3 && Grafted_Logical==TRUE) calc_Gsum_new(xInBox);

  check_zero_densities_owned(xOwned);              
  safe_free((void **) &xInBox);

  if (Proc==0 && Iwrite_screen ==SCREEN_VERBOSE){
        printf("\nInitial guess took %g secs\n",MPI_Wtime()-t1);
        printf("-----------------------------------------------------------------\n");
  }

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
