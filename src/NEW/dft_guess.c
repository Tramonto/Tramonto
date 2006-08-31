/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/*
 *  FILE: dft_guess.c
 *
 *  This file contains routines that allocate and pick an initial
 *  guess for the solution vector.
 *
 */

#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
 
double interpolate(double *, int, int, double **, double *);
void chop_profile(double *x, int);
void check_zero_densities(double *);
void communicate_profile(double *,double *);
void shift_the_profile(double *,double);
int  find_length_of_file(char *);
void read_in_a_file(int,char *);
void setup_step_2consts(double *);
void setup_linear_profile(double *);
void setup_exp_density(double *, double *,int,int);
void setup_const_density(double *, double *,int,int);
void setup_stepped_profile(double *);
void setup_elec_pot(double *,int);
void setup_chem_pot(double *);
void setup_polymer(double *);
void setup_polymer_rho(double *,int);
void setup_polymer_field(double *,int);
void setup_polymer_field_step(double *,int);
void setup_polymer_simple(double *, int);
void setup_polymer_G(double *);
void setup_polymer_G_2(double *);
void setup_rho_bar(double *);
int locate_inode_old(int *);

void set_initial_guess (int iguess,double *x,int *niters)
{
 /*
  * Local variable declarations
  */

  char *output_file4 = "rho_init.dat";
  char filename[20];
  double t1=0.0;
  double *x_new;

  int iunk;
  int i1,start_no_info;

  int idim,iwall,iwall_type;
  double sep,fac=0.0;

  
  /********************** BEGIN EXECUTION ************************************/

  if (Proc==0) { 
    /*  printf("\n%s: Setting an initial guess ... ",yo);*/
      t1 = MPI_Wtime();
  }

 /*
  * Allocate array, and point local x to the start of it
  */

  if (Proc==0 && Iwrite==VERBOSE){
    if (Restart==0) printf("generating guess from scratch\n");  
    else printf("setting up guess from existing files\n");
  } 

  if (Restart > 0 || Imain_loop > 0 || (*niters == 1 && Load_Bal_Flag >= LB_TIMINGS) ){
        start_no_info = FALSE;

        x_new = (double *) array_alloc(1, Nnodes*Nunk_per_node, sizeof(double));

        if (Proc == 0) {  /* Proc 0 reads in the data file */

           if ( Imain_loop == 0 && *niters == 0){

              /* START FROM AN OLD FILE - otherwise all of the _old variables
                                    were set in collect_xold (dft_output.c) */

              if (Print_rho_type != PRINT_RHO_0 && Imain_loop>0 ) {
                    if (Imain_loop>0)   i1=Imain_loop-1;
                    else                i1=0;
                    
                    if (Lbinodal && iguess==BINODAL_FLAG)
                          sprintf(filename, "dft_dens2.%0d",i1);
                    else  sprintf(filename, "dft_dens.%0d",i1);
              }
              else{
                    if (Lbinodal && iguess==BINODAL_FLAG)
                         sprintf(filename,"dft_dens2.dat"); 
                    else sprintf(filename,"dft_dens.dat");
              }

              Nodes_old = find_length_of_file(filename);
              if (Restart==5) Nodes_old=Nnodes;  /* going to read in 1D file, need to dimension arrays for 2D or 3D system */
              if (Lbinodal && iguess==BINODAL_FLAG){
                 X2_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
              }
              else
                 X_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
              read_in_a_file(iguess,filename);
           }

           if (Nodes_old != Nnodes) {

if (Proc==0 && Iwrite != NO_SCREEN) printf("Nodes_old=%d  Nnodes=%d\n",Nodes_old,Nnodes);

             idim = Plane_new_nodes;
             iwall = 0;
             iwall_type = WallType[iwall];

             if (Nwall > 0 && Nwall < 3){
             if (Nwall == 1) {
               if (Type_bc[idim][0] == REFLECT)
                   sep = 2.0*(WallPos[idim][iwall] + 0.5*Size_x[0]);
               else if (Type_bc[idim][1] == REFLECT) 
                   sep = 2.0*(0.5*Size_x[0] - WallPos[idim][iwall]);
               else sep = Guess_range[1]; /* use old guess */
             }
             else if (Nwall == 2){
                   sep = (fabs(WallPos[idim][1] - WallPos[idim][0]));
             }
             else sep = Guess_range[1]; /* use old guess */

             if (sep < Guess_range[0])      start_no_info = TRUE;
             else if (sep > Guess_range[1]) fac = 1.0;
             else if (Guess_range[0]==Guess_range[1] || Guess_range[0]<0.0 || Guess_range[1]<0.0) fac=1.0;
             else fac = (sep - Guess_range[0])/(Guess_range[1]-Guess_range[0]);

             if (start_no_info) {
                 if(Iwrite != NO_SCREEN)printf("Starting from bulk intial guess everywhere\n");
             }
             else{
                 if(Iwrite != NO_SCREEN)printf("Mixing old soln and bulk with %9.6f percent of old profile\n",
                          fac*100.);
             }
             }
             else fac=1.0;

             if (!start_no_info){
                shift_the_profile(x_new,fac); 
             }
           }
           else{ 
             if (iguess==Iguess1){
               for (iunk=0; iunk<Nunknowns; iunk++){
                   x_new[iunk] = X_old[iunk];
               }
             }
             else if (iguess==BINODAL_FLAG && Lbinodal)
               for (iunk=0; iunk<Nunknowns; iunk++) x_new[iunk] = X2_old[iunk];
           }
           if (iguess==Iguess1)
              safe_free((void *) &X_old);
           else if (iguess==BINODAL_FLAG && Lbinodal)
              safe_free((void *) &X2_old);
        }

        AZ_broadcast((char *) &start_no_info, sizeof(int), 
                                     Aztec.proc_config, AZ_PACK);
        AZ_broadcast(NULL, 0, Aztec.proc_config, AZ_SEND);

        if (!start_no_info){
            communicate_profile(x_new,x);
            check_zero_densities(x);
            if (Restart == 2) chop_profile(x,iguess);
            if (Restart == 3) {
               if (Lsteady_state)        setup_chem_pot(x);
               if (Matrix_fill_flag >= 3) setup_rho_bar(x);
               if (Ipot_ff_c == COULOMB) setup_elec_pot(x,iguess);
            }
        }
        safe_free((void *) &x_new);


  } /* end of logic for restarting from a file or previous run */

  else start_no_info = TRUE;

  if (start_no_info) {  /* START FROM A SPECIFIED INITIAL GUESS */

 if (!Sten_Type[POLYMER_CR]) {
    switch(iguess){
      case CONST_RHO:    
            setup_const_density(x,Rho_b,Ncomp,0);
            break;

      case CONST_RHO_V:  
            setup_const_density(x,Rho_coex,1,1);
            break;
      case CONST_RHO_L:  
            setup_const_density(x,Rho_coex,1,0);
            break;

      case EXP_RHO:
            setup_exp_density(x,Rho_b,Ncomp,0);
            break;
      case EXP_RHO_V:
            setup_exp_density(x,Rho_coex,1,1);
            break;
      case EXP_RHO_L:
            setup_exp_density(x,Rho_coex,1,0);
            break;

      case STEP_PROFILE:
            setup_stepped_profile(x);

      case CHOP_RHO_L:
            setup_exp_density(x,Rho_coex,1,1);
            chop_profile(x,iguess);
            break;
      case CHOP_RHO_V:
            setup_exp_density(x,Rho_coex,1,0);
            chop_profile(x,iguess);
            break;
      case LINEAR:
            setup_linear_profile(x);
    }  /* end of iguess switch */
 }
 else{
    /*setup_polymer(x);*/

   if (Type_poly == 3)
     setup_polymer_simple(x,iguess);
   else
     setup_polymer_rho(x,iguess);
     /*if(iguess==STEP_PROFILE) setup_polymer_field_step(x,iguess);
     else  setup_polymer_field(x,iguess);*/
     setup_polymer_field_step(x,iguess);
     setup_polymer_G(x);
    /* setup_polymer_G_2(x); */

    check_zero_densities(x);
 }

 if (Matrix_fill_flag >= 3) setup_rho_bar(x);
 if (Ipot_ff_c == COULOMB)  setup_elec_pot(x,iguess);
 if (Lsteady_state)         setup_chem_pot(x);

   } /* end of cases for setting initial guess from scratch */


 /*
  * If requested, write out initial guess
  */
     if (Iwrite == VERBOSE){
         if (Proc == 0){
             X_old = (double *) array_alloc (1, Nnodes*Nunk_per_node, sizeof(double));
             Vext_old = (double *) array_alloc (1, Nnodes*Ncomp, sizeof(double));
         }
             

         collect_x_old(x,0);
         collect_vext_old();

         if (Proc ==0) {
             if (iguess==Iguess1) print_profile(output_file4);
             else if (Lbinodal && iguess==BINODAL_FLAG) print_profile("rho_init2.dat");
             safe_free((void *) &X_old);
             safe_free((void *) &Vext_old);
         }
     }

  if (Proc==0 && Iwrite==VERBOSE) printf("\n initial guess took %g secs\n",MPI_Wtime()-t1);
  return;
}
/*********************************************************/
/*setup_polymer_field_step: in this routine set up the field guesses
                       for the polymers variables    */
void setup_polymer_field_step(double *x, int iguess)
{
  int loc_inode,itype_mer,loc_RHO,loc_i;
  double field;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

     for (itype_mer=0; itype_mer<Ncomp; itype_mer++){
         loc_RHO = Aztec.update_index[Ncomp+itype_mer+Nunk_per_node*loc_inode];
         loc_i = Aztec.update_index[itype_mer+Nunk_per_node*loc_inode];
         if (x[loc_RHO]<1.e-6) field=VEXT_MAX-1.;
         else field=-log(x[loc_RHO]/Rho_b[itype_mer]);
         x[loc_i]=exp(-field);
     }
   }
   return;
}
/*********************************************************/
/*setup_polymer_field: in this routine set up the field guesses
                       for the polymers variables    */
void setup_polymer_field(double *x, int iguess)
{
  int loc_inode,inode_box,ijk_box[3],loc_i,i,inode;
  int reflect_flag[NDIM_MAX];
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight, weight_bulk;
  struct Stencil_Struct *sten;
  double sign,temp;
  int sten_type,izone,jlist,itype_mer,jnode_box,jtype_mer;
  int jnode,loc_j;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);

     sten_type = POLYMER_CR;
     izone = 0;

     for (itype_mer=0; itype_mer<Ncomp; itype_mer++){
        loc_i = Aztec.update_index[itype_mer+Nunk_per_node*loc_inode];
        if (!Zero_density_TF[inode_box][itype_mer]){
           temp = 0.;
           for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
             if (Nlists_HW <= 2) jlist = 0;
             else                jlist = jtype_mer;
           
             sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
             sten_offset = sten->Offset;
             sten_weight = sten->Weight;

             for (isten = 0; isten < sten->Length; isten++) {
                offset = sten_offset[isten];
                weight_bulk = sten_weight[isten];

               /* Find the Stencil point */
               jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

	      if (jnode_box >= 0 ) {
	         if (Lhard_surf) {
		     if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
                        weight = HW_boundary_weight (jtype_mer,jlist,
                                 sten->HW_Weight[isten], jnode_box, reflect_flag);
		     else weight = weight_bulk;
		 }
		 else weight = weight_bulk;

	         if (iguess == CONST_RHO) {
		   if (!Zero_density_TF[jnode_box][jtype_mer])
		     temp -=  (weight - weight_bulk)*Rho_b[jtype_mer];
		   else
		     temp +=  weight_bulk*Rho_b[jtype_mer];
		 }
	         else if (iguess == STEP_PROFILE) {
		   jnode = B2G_node[jnode_box];
                   node_to_position(jnode,nodepos);
                   for (i=0;i<Nsteps;i++){
                        if (nodepos[Orientation_step[i]]>=Xstart_step[i] &&
                            nodepos[Orientation_step[i]]<=Xend_step[i]){
		            if (!Zero_density_TF[jnode_box][jtype_mer]){
		                temp -= (weight*Rho_step[jtype_mer][i] - weight_bulk*Rho_b[jtype_mer]);
                            }
		            else temp +=  weight_bulk*Rho_b[jtype_mer];
                        }
                   }
                 }
              }
	      else if ( jnode_box == -2 ){  /*in wall*/
		     temp +=  weight_bulk*Rho_b[jtype_mer]; 
                   /* jnode_box == -1 = in bulk contributes nothing  */
	      }
	     } /* end loop over stencil */
           }
	   temp += Vext[loc_inode][itype_mer];

           if (temp < -1.0) temp=-1.0;

           if ((Sten_Type[POLYMER_CR] == 1) || (temp > 0.5)) x[loc_i] = exp(-temp);
           else                             x[loc_i] = exp(1. - sqrt(1.-2.*temp));

        } /* end if i not zero density  */
        else   x[loc_i] = DENSITY_MIN; /* zero density - Boltzmann probability = 0 */

     } /* end loop over itype_mer */
  } /* end loop over loc_inode */

  return;
}

/*********************************************************/
/*setup_polymer_simple: in this routine set up the field guesses
                       for the polymers variables for SCF case    */
void setup_polymer_simple(double *x, int iguess)
{
  int loc_inode,inode_box,ijk_box[3],loc_i,loc_j,i;
  double temp;
  int itype_mer,jtype_mer,inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);

     for (itype_mer=0; itype_mer<Ncomp; itype_mer++){
        loc_i = Aztec.update_index[itype_mer+Nunk_per_node*loc_inode];
        if (!Zero_density_TF[inode_box][itype_mer]){
           temp = 0.;
           for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
             loc_j = Aztec.update_index[jtype_mer+Nunk_per_node*loc_inode];
           
	     if (iguess == CONST_RHO) {
	       temp = 0.;
	     }
	     else if (iguess == STEP_PROFILE) {
	       inode = B2G_node[inode_box];
               node_to_position(inode,nodepos);
               for (i=0;i<Nsteps;i++){
                 if ( nodepos[Orientation_step[i]]>=Xstart_step[i] &&
                      nodepos[Orientation_step[i]]<=Xend_step[i] ) {
	            if (!Zero_density_TF[inode_box][jtype_mer])
		      temp += Rism_cr[itype_mer][jtype_mer][0]*(Rho_step[jtype_mer][i]
							   - Rho_b[jtype_mer]);
	            else temp -=  Rism_cr[itype_mer][jtype_mer][0]*Rho_b[jtype_mer];
                 }
               }
	     }

           } /* end loop over jtype_mer */
	   /* ALF: add Vext to initial guess for field */
	   temp += Vext[loc_inode][itype_mer];
           /* if (temp > 1.0) temp=1.0;*/
           if ((Sten_Type[POLYMER_CR]) || (temp > 0.5)) x[loc_i] = exp(temp);
           else                             x[loc_i] = exp(1. - sqrt(1.-2.*temp));
        } /* end if i not zero density  */
        else   x[loc_i] = DENSITY_MIN; /* zero density - Boltzmann probability = 0 */
     } /* end loop over itype_mer */
  } /* end loop over loc_inode */
  return;
}
/*********************************************************/

/*********************************************************/
/*setup_polymer_rho: in this routine set up polymer density profiles    */
void setup_polymer_rho(double *x, int iguess)
{
  int loc_inode,i,inode_box,ijk_box[3],loc_i,icomp;
  int inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);
     if (iguess == CONST_RHO) {
       for (i=0; i<Ncomp; i++){
         loc_i = Aztec.update_index[i+  Ncomp  +Nunk_per_node*loc_inode];
         if (!Zero_density_TF[inode_box][i])
            x[loc_i] = Rho_b[i];
         else x[loc_i] = DENSITY_MIN;
       }
     }
     else if (iguess == STEP_PROFILE) {
       inode = B2G_node[inode_box]; 
       node_to_position(inode,nodepos);
       for (i=0;i<Nsteps;i++){
           if (nodepos[Orientation_step[i]]>=Xstart_step[i] &&
                nodepos[Orientation_step[i]]<=Xend_step[i]){
                for (icomp=0;icomp<Ncomp;icomp++){
if (loc_inode==20) printf("i=%d Xstart_step=%g  Xend_step=%g nodepos[%d]=%g Rho_step[%d]=%g\n",
                            i,Xstart_step[i],Xend_step[i],Orientation_step[i],nodepos[0],icomp,Rho_step[icomp][i] );
	           loc_i = Aztec.update_index[icomp+Ncomp+Nunk_per_node*loc_inode];
       		   if (!Zero_density_TF[inode_box][icomp]) x[loc_i]= Rho_step[icomp][i];
		   else x[loc_i]=0.0;
                }
            }
       }
     }
     else {
       if (Proc==0) printf("invalid initial guess for polymers\n");
       exit(1);
     }
  }
  return;
}
/*********************************************************/
/*setup_polymer_G: in this routine set up guess for the G's   */
/* in this version, guess is simply the Boltzmann factors, with some account taken of hard walls*/
void setup_polymer_G(double *x)
{
  int loc_inode,inode_box,ijk_box[3],loc_i;
  int reflect_flag[NDIM_MAX];
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  int sten_type,izone,jlist,jnode_box,jtype_mer,itype_mer;
  int loc_j,iunk,poln,iseg,ibond,not_done,junk,cycle,loc_B;

     if (Sten_Type[POLYMER_GAUSS]) sten_type = POLYMER_GAUSS;
     else                          sten_type = DELTA_FN;
     izone = 0;

     loc_inode=0;
     for (poln=0; poln < Npol_comp; poln++){
        iunk = Geqn_start[poln];
        for (iseg=0; iseg<Nmer[poln]; iseg++){
             for (ibond=0; ibond<Nbond[poln][iseg]; ibond++){
                 loc_i = Aztec.update_index[iunk + Nunk_per_node*loc_inode];
                 x[loc_i] =999.0;
                 iunk++;
             }
        }
     }
     not_done=TRUE;

     cycle=0;
     while (not_done) {
     for (poln=0; poln < Npol_comp; poln++){
        iunk = Geqn_start[poln];
        for (iseg=0; iseg<Nmer[poln]; iseg++){
           itype_mer =Type_mer[poln][iseg];
           for (ibond=0; ibond<Nbond[poln][iseg]; ibond++){

                             /* TREAT THE END SEGMENTS */
             /* only try to generate the iunk guess if not already filled in */
             if (fabs(x[Aztec.update_index[iunk]]-999.0)<1.e-6){

             if(Bonds[poln][iseg][ibond]== -1){
                 if (Iwrite==VERBOSE) printf("computing end segment poln=%d iseg=%d ibond=%d this cycle=%d \n",poln,iseg,ibond,cycle);
                 for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
                     inode_box = L2B_node[loc_inode];
                     node_box_to_ijk_box(inode_box, ijk_box);
                     loc_i = Aztec.update_index[iunk + Nunk_per_node*loc_inode];
                     loc_j = Aztec.update_index[itype_mer+Nunk_per_node*loc_inode];
                     if (Type_poly == 2) x[loc_i] = 1.0;  
                     else x[loc_i] = x[loc_j];
                  }
             }
             else{
               jtype_mer = Type_mer[poln][Bonds[poln][iseg][ibond]];
               junk = Geqn_start[poln]+2*Bonds[poln][iseg][ibond];
               if (iseg<Bonds[poln][iseg][ibond]) junk += 1;
               /* test if this G equation has been generated yet ... if not, go on */
               if (fabs(x[Aztec.update_index[junk]]-999.0)>1.e-6){
                 if (Iwrite==VERBOSE) printf("computing poln=%d iseg=%d ibond=%d this cycle=%d \n",poln,iseg,ibond,cycle);
                   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
                       inode_box = L2B_node[loc_inode];
                       node_box_to_ijk_box(inode_box, ijk_box);
                       loc_i = Aztec.update_index[iunk + Nunk_per_node*loc_inode];
                       x[loc_i] = 0.;

                   if (!Zero_density_TF[inode_box][jtype_mer]){ 
                       if (Nlists_HW <= 2) jlist = 0;
                       else                jlist = jtype_mer; 
           
		       /* generalize to != bond lengths */
                       sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
                       sten_offset = sten->Offset;
                       sten_weight = sten->Weight;
                       for (isten = 0; isten < sten->Length; isten++) {
                          offset = sten_offset[isten];
                          weight = sten_weight[isten];

                           /* Find the Stencil point */
                           jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

                           if (jnode_box >= -1 ) {  /* (-1 in bulk) */
                              if (Lhard_surf) {
                              if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
                                  weight = HW_boundary_weight 
                                   (jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
                           }
                           if (B2L_node[jnode_box] >-1){ /* node in domain */
                              loc_j = Aztec.update_index[junk+Nunk_per_node*B2L_node[jnode_box]];
                              /*x[loc_i] +=  weight*x[loc_j]; */
                              x[loc_i] +=  weight; 
                           } /* check that node is in domain */
                           else{  /* use the value at loc_inode ... an approximation */
                              loc_j = Aztec.update_index[junk+Nunk_per_node*loc_inode];
                              x[loc_i] +=  weight*x[loc_j]; 
                           }
                         }
                      }

                  loc_B = Aztec.update_index[itype_mer +Nunk_per_node*loc_inode];/*Boltz */
                  x[loc_i] *= x[loc_B];

                   /* simpler guess for debugging purposes LJDF*/
                  if (Type_poly == 2)  x[loc_i] = 1.0;
                    
	           } /*end of Zero_dens_TF test */
                   } /* end of loop over loc_inode */
               }  /* end of test on whether the jtype_mer guess exists */
             } /* end of if Bond test */
             }  /* end of test of whether the itype_mer guess already has been generated */
             iunk++;
           } /* end of loop over bonds on iseg */
        } /* end of loop over each segment on poln */
     } /* end of loop over polymer chains */

     /* test to see if all of the G equations have been generated yet */
     not_done=FALSE;
     for (poln=0; poln < Npol_comp; poln++){
        iunk = Geqn_start[poln];
        for (iseg=0; iseg<Nmer[poln]; iseg++){
             for (ibond=0; ibond<Nbond[poln][iseg]; ibond++){
                 if (fabs(x[Aztec.update_index[iunk]]-999.0)<1.e-6) not_done=TRUE;
                /* else{
                    printf("poln=%d iseg=%d ibond=%d is done \n",poln,iseg,ibond);
                 }*/
                 iunk++;
             }
        }
     }
     cycle++;
    

     } /* end of while test */

  return;
}

/*********************************************************/
/*setup_polymer_G_2: in this routine set up guess for the G's   */
/* version 2: calculate G's from initial field guess
   linear polymers only for now!!! */
void setup_polymer_G_2(double *x)
{
  int loc_inode,inode_box,ijk_box[3],loc_i;
  int reflect_flag[NDIM_MAX];
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  int sten_type,izone,jlist,jnode_box,itype_mer,jtype_mer;
  int loc_j,iunk,poln,iseg,ibond,j_box;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);

     if (Sten_Type[POLYMER_GAUSS]) sten_type = POLYMER_GAUSS;
     else                          sten_type = DELTA_FN;
     izone = 0;

     for (poln=0; poln < Npol_comp; poln++){
       /* do G's first in forward direction */
        iunk = Geqn_start[poln];
	ibond = 0;  /* only doing one set of bonds at a time */
        for (iseg=0; iseg<Nmer[poln]; iseg++){
	  itype_mer = Type_mer[poln][iseg];
	  loc_i = Aztec.update_index[iunk + Nunk_per_node*loc_inode];
	  x[loc_i] = 0.;

                             /* TREAT THE END SEGMENTS */
	  if(Bonds[poln][iseg][ibond]== -1){
	    loc_j = Aztec.update_index[Type_mer[poln][iseg]+Nunk_per_node*loc_inode];
	    x[loc_i] = x[loc_j];
	    if (Type_poly == 2) x[loc_i] = 1.0; 
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
	  else{
	    jtype_mer = Type_mer[poln][Bonds[poln][iseg][ibond]];
	    /* if (!Zero_density_TF[inode_box][jtype_mer]){ */
	      if (Nlists_HW <= 2) jlist = 0;
	      else                jlist = itype_mer;
           
	      /* generalized for != bond lengths */
	      sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
	      sten_offset = sten->Offset;
	      sten_weight = sten->Weight;

	      for (isten = 0; isten < sten->Length; isten++) {
		offset = sten_offset[isten];
		weight = sten_weight[isten];

		/* Find the Stencil point */
		jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

		if (jnode_box >= 0 ) {  
		  if (Lhard_surf) {
		    if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
		      weight = HW_boundary_weight 
			(jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
		  }
		  j_box = jnode_box*Nunk_per_node + iunk;
		  x[loc_i] +=  weight*x[B2L_unknowns[j_box]]; 
		}
		else if (jnode_box == -1 )
		  x[loc_i] += weight;
	      } /* end loop over stencil */

	      /* Boltzmann factor */
	      loc_j = Aztec.update_index[itype_mer +Nunk_per_node*loc_inode];
	      x[loc_i] *= x[loc_j];

	      /* simpler guess for debugging purposes LJDF*/
	      if (Type_poly == 2)  x[loc_i] = 1.0;
                    
	      /* } end of Zero_dens_TF test */
	  } /* end of if Bond test */
	  iunk++;
	  MPI_Barrier(MPI_COMM_WORLD);
        } /* end of loop over each segment on poln */

	/* now do G's in backward direction */
        iunk = Ngeqn_tot;
	ibond = 1;  /* only doing one set of bonds at a time */
        for (iseg=Nmer[poln]-1; iseg>-1; iseg--){
	  itype_mer = Type_mer[poln][iseg];
	  loc_i = Aztec.update_index[iunk + Nunk_per_node*loc_inode];
	  x[loc_i] = 0.;

                             /* TREAT THE END SEGMENTS */
	  if(Bonds[poln][iseg][ibond]== -1){
	    loc_j = Aztec.update_index[Type_mer[poln][iseg]+Nunk_per_node*loc_inode];
	    x[loc_i] = x[loc_j];
	    if (Type_poly == 2) x[loc_i] = 1.0; 
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
	  else{
	    jtype_mer = Type_mer[poln][Bonds[poln][iseg][ibond]];
	    if (!Zero_density_TF[inode_box][jtype_mer]){
	      if (Nlists_HW <= 2) jlist = 0;
	      else                jlist = jtype_mer;
           
	      /* generalized for != bond lengths */
	      sten = &(Stencil[sten_type][izone][itype_mer+Ncomp*jtype_mer]);
	      sten_offset = sten->Offset;
	      sten_weight = sten->Weight;

	      for (isten = 0; isten < sten->Length; isten++) {
		offset = sten_offset[isten];
		weight = sten_weight[isten];

		/* Find the Stencil point */
		jnode_box = offset_to_node_box(ijk_box, offset, reflect_flag);

		if (jnode_box >= 0 ) {  
		  if (Lhard_surf) {
		    if (Nodes_2_boundary_wall[jlist][jnode_box]!=-1) 
		      weight = HW_boundary_weight 
			(jtype_mer,jlist,sten->HW_Weight[isten], jnode_box, reflect_flag);
		  }
		  j_box = jnode_box*Nunk_per_node + iunk;
		  x[loc_i] +=  weight*x[B2L_unknowns[j_box]]; 
		}
		else if (jnode_box == -1 )
		  x[loc_i] += weight;
	      } /* end loop over stencil */

	      /* Boltzmann factor */
	      loc_j = Aztec.update_index[itype_mer +Nunk_per_node*loc_inode];
	      x[loc_i] *= x[loc_j];

	      /* simpler guess for debugging purposes LJDF*/
	      if (Type_poly == 2)  x[loc_i] = 1.0;
                    
	    } /*end of Zero_dens_TF test */
	  } /* end of if Bond test */
	  iunk--;
	  MPI_Barrier(MPI_COMM_WORLD);
        } /* end of loop over each segment on poln */

     } /* end of loop over polymer chains */
  }  /* end of loop over nodes */

  return;
}

/*********************************************************/
/*setup_const_density: in this routine set up a constant
        density profile wherever Zero_density_TF = FALSE */
void setup_const_density(double *x, double *rho,int nloop,int index)
{
  int loc_inode,i,inode_box,loc_i;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0; i<nloop; i++){
         loc_i = Aztec.update_index[i+Nunk_per_node*loc_inode];
         if (!Zero_density_TF[inode_box][i]){
            if (nloop > 1) x[loc_i] = rho[i];
            else           x[loc_i] = rho[index];
         }
         else x[loc_i] = 0.0;
         if (i==1 && !Zero_density_TF[inode_box][1] &&
                      Zero_density_TF[inode_box][0]) x[loc_i] = rho[0];

         /* set up initial guess if chemical potential is an unknown */
         if (Lsteady_state) {
             if (Ipot_ff_c == COULOMB) 
                  loc_i = Aztec.update_index[i+Ncomp+1+Nunk_per_node*loc_inode];
             else loc_i = Aztec.update_index[i+Ncomp+Nunk_per_node*loc_inode];

             if (!Zero_density_TF[inode_box][i]){
                if (nloop > 1) x[loc_i] = Betamu[i];
                else           x[loc_i] = Betamu[index];
             }
             else x[loc_i] = -10.0; /* zero density == -infinite chem.pot*/
         }
      }
  }
  return;
}
/*********************************************************/
/*setup_stepped_profile: in this routine set up a stepped
        density profile wherever Zero_density_TF = FALSE */
void setup_stepped_profile(double *x)
{
  int loc_inode,i,inode_box,loc_i,icomp,inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode = L2G_node[loc_inode];
     node_to_position(inode,nodepos);

     for (i=0;i<Nsteps;i++){
       if (nodepos[Orientation_step[i]] >= Xstart_step[i] &&
           nodepos[Orientation_step[i]] <= Xend_step[i]){
           for (icomp=0; icomp<Ncomp; icomp++){
               loc_i = Aztec.update_index[icomp+Nunk_per_node*loc_inode];
               if (!Zero_density_TF[inode_box][icomp]){
                   x[loc_i]=Rho_step[icomp][i];
               }
               else x[loc_i]=0.0;
           }
       }
    }
  }

  if (Lsteady_state){
     printf("stepped profile not set up for chemical potentials at this time\n");
     exit(-1);
  }
  return;
}
/*********************************************************/
/*setup_exp_density: in this routine set up a density 
                     profile as rho_b*exp(-Vext/kT)*/
void setup_exp_density(double *x, double *rho,int nloop,int index)
{

  int loc_inode,i,inode_box,loc_i;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0; i<nloop; i++) {

        loc_i = Aztec.update_index[i+Nunk_per_node*loc_inode];
        if (!Zero_density_TF[inode_box][i]){
            if (nloop > 1) x[loc_i] = rho[i] * exp(-Vext[loc_inode][i]);
            else           x[loc_i] = rho[index] * exp(-Vext[loc_inode][i]);
        }
        else x[loc_i] = 0.0;

          if (x[loc_i]>1.0) x[loc_i] = 0.99;   /* may need this to
            make sure that rb3<1 always */


         /* set up initial guess if chemical potential is an unknown */
         if (Lsteady_state) {
             if (Ipot_ff_c == COULOMB) 
                  loc_i = Aztec.update_index[i+Ncomp+1+Nunk_per_node*loc_inode];
             else loc_i = Aztec.update_index[i+Ncomp+Nunk_per_node*loc_inode];

             if (!Zero_density_TF[inode_box][i]){
                if (nloop > 1) x[loc_i] = Betamu[i];
                else           x[loc_i] = Betamu[index];
             }
             else x[loc_i] = -10.0; /* zero density == -infinite chem.pot*/
         }
     }
  }
  return;
}
/************************************************************/
/*setup_step_2consts: density profile where the left (LBB)
      and right (RTF) sides of the box have different bulk
      densities.  This can either be for chem.pot. gradients
      or for a liquid-vapor profile.  */

void setup_step_2consts(double *x)
{

  int loc_inode,icomp,inode_box,loc_i,ijk[3],inode;
  double x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0] - 0.5*Size_x[0];

     for (icomp=0; icomp<Ncomp; icomp++) {

        loc_i = Aztec.update_index[icomp+Nunk_per_node*loc_inode];
        if (!Zero_density_TF[inode_box][icomp]){

           if (x_dist< 0.0) {
              if (Lsteady_state) x[loc_i] = Rho_b_LBB[icomp];
              else               x[loc_i] = Rho_coex[1];
           }
           else {
              if (Lsteady_state) x[loc_i] = Rho_b_RTF[icomp];
              else               x[loc_i] = Rho_coex[0];
           }
        }
        else x[loc_i] = 0.0;


       /* set up initial guess if chemical potential is an unknown */
      /*  if (Lsteady_state) {
            if (Ipot_ff_c == COULOMB) 
                 loc_i = Aztec.update_index[icomp+Ncomp+1+Nunk_per_node*loc_inode];
            else loc_i = Aztec.update_index[icomp+Ncomp+Nunk_per_node*loc_inode];

            if (!Zero_density_TF[inode_box][icomp]){
               if (x_dist< 0.0) x[loc_i] = Betamu_LBB[icomp];
               else             x[loc_i] = Betamu_RTF[icomp];
            }
            else x[loc_i] = -10.0;  zero density == -infinite chem.pot
        }*/
     }
  }
  return;
}
/************************************************************/
/*setup_linear profile: density profile that interpolates 
      linearly between densities at the left (LBB) and 
      right (RTF) sides of the box. This can either be for 
      chem.pot. gradients
      or for a liquid-vapor profile.  */

void setup_linear_profile(double *x)
{

  int loc_inode,icomp,inode_box,loc_i,ijk[3],inode;
  double x_dist,x_tot;

  x_tot=Size_x[Grad_dim]-2.0*X_const_mu;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0]-X_const_mu;

     for (icomp=0; icomp<Ncomp; icomp++) {

        loc_i = Aztec.update_index[icomp+Nunk_per_node*loc_inode];
        if (!Zero_density_TF[inode_box][icomp]){

           if (Lsteady_state){
             if (x_dist <0.) x[loc_i]=Rho_b_LBB[icomp];
             else if (x_dist > x_tot) x[loc_i]= Rho_b_RTF[icomp];
             else  x[loc_i] = Rho_b_LBB[icomp] + (Rho_b_RTF[icomp]-Rho_b_LBB[icomp])* x_dist/x_tot;
           }
           else 
             x[loc_i] = Rho_coex[1] + (Rho_coex[0]-Rho_coex[1])*
                                        x_dist/Size_x[Grad_dim];
        }
        else x[loc_i] = 0.0;


     }
  }
  return;
}
/************************************************************/
/*setup_rho_bar: set up the rhobar initial guesses.  For now
                use Rho_b to set initial guess.  Later calculate
                based on rho initial guess. */
void setup_rho_bar(double *x)
{
  int loc_inode,loc_i,inode_box,inode,ijk[3],icomp,idim,iunk;
  double vol,area,x_dist;

  for (iunk=0; iunk<Nrho_bar; iunk++){
    Rhobar_b[iunk] = 0.0;
    Rhobar_b_LBB[iunk] = 0.0;
    Rhobar_b_RTF[iunk] = 0.0;
  }

  for (icomp=0; icomp<Ncomp; icomp++){
    vol = PI*Sigma_ff[icomp][icomp]*
             Sigma_ff[icomp][icomp]*Sigma_ff[icomp][icomp]/6.0;
    area = PI*Sigma_ff[icomp][icomp]*Sigma_ff[icomp][icomp];

    Rhobar_b[0] += vol*Rho_b[icomp];
    Rhobar_b[1] += area*Rho_b[icomp];
    Rhobar_b[2] += area*Rho_b[icomp]*Inv_4pir[icomp];
    Rhobar_b[3] += area*Rho_b[icomp]*Inv_4pirsq[icomp];
    for (idim=0; idim<Ndim; idim++){
        Rhobar_b[4+idim] = 0.0;
        Rhobar_b[4+Ndim+idim] = 0.0;
    }
    if (Lsteady_state){
       Rhobar_b_LBB[0] += vol*Rho_b_LBB[icomp];
       Rhobar_b_LBB[1] += area*Rho_b_LBB[icomp];
       Rhobar_b_LBB[2] += area*Rho_b_LBB[icomp]*Inv_4pir[icomp];
       Rhobar_b_LBB[3] += area*Rho_b_LBB[icomp]*Inv_4pirsq[icomp];
       for (idim=0; idim<Ndim; idim++){
           Rhobar_b_LBB[4+idim] = 0.0;
           Rhobar_b_LBB[4+Ndim+idim] = 0.0;
       }
       Rhobar_b_RTF[0] += vol*Rho_b_RTF[icomp];
       Rhobar_b_RTF[1] += area*Rho_b_RTF[icomp];
       Rhobar_b_RTF[2] += area*Rho_b_RTF[icomp]*Inv_4pir[icomp];
       Rhobar_b_RTF[3] += area*Rho_b_RTF[icomp]*Inv_4pirsq[icomp];
       for (idim=0; idim<Ndim; idim++){
           Rhobar_b_RTF[4+idim] = 0.0;
           Rhobar_b_RTF[4+Ndim+idim] = 0.0;
       }
       
    }
    else if (Nwall == 0 && Iliq_vap == 3){
       Rhobar_b_LBB[0] += vol*Rho_coex[1];
       Rhobar_b_LBB[1] += area*Rho_coex[1];
       Rhobar_b_LBB[2] += area*Rho_coex[1]*Inv_4pir[icomp];
       Rhobar_b_LBB[3] += area*Rho_coex[1]*Inv_4pirsq[icomp];
       for (idim=0; idim<Ndim; idim++){
           Rhobar_b_LBB[4+idim] = 0.0;
           Rhobar_b_LBB[4+Ndim+idim] = 0.0;
       }
       Rhobar_b_RTF[0] += vol*Rho_coex[0];
       Rhobar_b_RTF[1] += area*Rho_coex[0];
       Rhobar_b_RTF[2] += area*Rho_coex[0]*Inv_4pir[icomp];
       Rhobar_b_RTF[3] += area*Rho_coex[0]*Inv_4pirsq[icomp];
       for (idim=0; idim<Ndim; idim++){
           Rhobar_b_RTF[4+idim] = 0.0;
           Rhobar_b_RTF[4+Ndim+idim] = 0.0;
       }
    }
  }

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (iunk = 0; iunk < Nrho_bar; iunk++){
       loc_i = Aztec.update_index[Ncomp+iunk+Nunk_per_node*loc_inode];
       if (Lsteady_state || (Nwall == 0 && Iliq_vap == 3)){
           inode_box = L2B_node[loc_inode];
           inode     = B2G_node[inode_box];
           node_to_ijk(inode,ijk); 
           x_dist = Esize_x[Grad_dim]*ijk[Grad_dim];

           x[loc_i] = Rhobar_b_LBB[iunk] + 
                      (Rhobar_b_RTF[iunk]-Rhobar_b_LBB[iunk])*
                           x_dist/Size_x[Grad_dim];
       }
       else {
          x[loc_i] = Rhobar_b[iunk];
       }
     }
  }
  return;
}
/************************************************************/
/*setup_elec_pot: set up the electrostatic potential initial guess*/
void setup_elec_pot(double *x,int iguess)
{
  int loc_inode,loc_i,inode_box,inode,ijk[3];
  double x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       if (Type_poly==-1) loc_i = Aztec.update_index[Ncomp+Nrho_bar+Ngeqn_tot+Nunk_per_node*loc_inode];
       else               loc_i = Aztec.update_index[2*Ncomp+Ngeqn_tot+Nunk_per_node*loc_inode];
       if (Lsteady_state && iguess==LINEAR){
           inode_box = L2B_node[loc_inode];
           inode     = B2G_node[inode_box];
           node_to_ijk(inode,ijk); 
           x_dist = Esize_x[Grad_dim]*ijk[Grad_dim];

           x[loc_i] = Elec_pot_LBB + 
                      (Elec_pot_RTF-Elec_pot_LBB)*
                           x_dist/Size_x[Grad_dim];
       }
       else x[loc_i] = 1.0;
  }
  return;
}
/************************************************************/
/* setup_chem_pot: for cases with steady state profiles,
   set up an initial guess for (electro)chemical potentials */
void setup_chem_pot(double *x)
{
  int loc_inode,inode_box,inode,ijk[3],icomp,iunk,loc_i;
  double x_dist,x_tot;

  x_tot = Size_x[Grad_dim]-2.*X_const_mu;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[Grad_dim]*ijk[Grad_dim]-X_const_mu;

     for (icomp=0; icomp<Ncomp; icomp++){
        if (Ipot_ff_c == 1) iunk = Ncomp+Nrho_bar+1+icomp;
        else                iunk = Ncomp+Nrho_bar+icomp;
        loc_i = Aztec.update_index[iunk+Nunk_per_node*loc_inode];
       if (!Zero_density_TF[inode_box][icomp]){
        if (Ipot_ff_c == 1){
             x[loc_i] = log(x[Aztec.update_index[icomp+Nunk_per_node*loc_inode]])
                           +Charge_f[icomp]*(x[Aztec.update_index[Ncomp+Nunk_per_node*loc_inode]]);

        }
        else{
            if (x_dist<0.) x[loc_i]=Betamu_LBB[icomp];
            else if (x_dist > x_tot) x[loc_i]=Betamu_RTF[icomp];
            else  x[loc_i] = Betamu_LBB[icomp] + (Betamu_RTF[icomp]-Betamu_LBB[icomp])* x_dist/x_tot;
        }
       }
       else x[loc_i] = -VEXT_MAX;
     }
   
  }
  return;
}
/************************************************************/
/*find_length_of_file: here just establish a file length.*/
int find_length_of_file(char *filename)
{
  int c,nodes_old=0;
  FILE *fp;
  if (Proc==0) {
    fp=fopen(filename,"r");

    nodes_old=0;
    while ( (c=getc(fp)) != EOF){
      if (c == '\n') {
	nodes_old++;
	c=getc(fp);/*so that a following blank line won't increment nodes_old*/
      }
    }
    fclose(fp);
  }
  return(nodes_old);
}
/**************************************************************/
/*read_in_a_file: In this routine, we read in a file containing
                  an old solution and use it for our initial guess */
void read_in_a_file(int iguess,char *filename)
{
  int c;
  int i,iunk,idim, inode,nunk_tot,itype_mer,ipol,iseg,index,dim_tmp;
  int ijk_old[3],ijk_old_max[3],open_now,ndim_max;
  double pos_old,pos_old0[3],tmp;
  char filename2[20];
  FILE *fp5=NULL,*fp6=NULL;

  /* read positions from file find Nodes_x_old[idim] */
  /* read the densities and electrostatic potentials from file */

  open_now=TRUE;
  for (idim=0;idim<Ndim;idim++) ijk_old_max[idim] = 0;
  for (index=0; index<Nodes_old; index++) {
    if (open_now){
       fp5=fopen(filename,"r");
       if (Sten_Type[POLYMER_CR]){
         sprintf(filename2,"%sg",filename);
         fp6=fopen(filename2,"r");
       }
      open_now=FALSE;
    }
    if (Restart==5) ndim_max=1;
    else ndim_max=Ndim;
    for (idim=0;idim<ndim_max;idim++)         {
      fscanf(fp5,"%lf",&pos_old);
      if (index==0) pos_old0[idim]=pos_old;
      pos_old -= pos_old0[idim];
  /*NEW CODE .... use this to start with x guess and transform to y etc.*/
      /*if (idim==0) dim_tmp=1;
      else if (idim==1) dim_tmp=0;
      else*/   dim_tmp=idim;
      ijk_old[dim_tmp] = round_to_int(pos_old/Esize_x[dim_tmp]);
 /*END OF NEW CODE */
      if (Sten_Type[POLYMER_CR]) {
	fscanf(fp6,"%lf",&tmp);
      }
      if (ijk_old[dim_tmp] > ijk_old_max[dim_tmp]){
	ijk_old_max[dim_tmp] = ijk_old[dim_tmp];
      } 
    }
     /* fscanf(fp5,"%lf",&tmp);*/
    inode=ijk_to_node(ijk_old);
    if (Restart==5) inode=index;
    if (Restart == 3) {
        if (Sten_Type[POLYMER_CR]) nunk_tot=2*Ncomp+Ngeqn_tot;
        else                       nunk_tot = Ncomp;
    }
    else              nunk_tot = Nunk_per_node;
    for (iunk=0;iunk<nunk_tot;iunk++) {
      if (Sten_Type[POLYMER_CR]){
	if (iunk < Ncomp){
            if (Lbinodal && iguess==BINODAL_FLAG)
   	         fscanf(fp5,"%lf",&X2_old[Ncomp+inode*Nunk_per_node+iunk]);
            else
   	         fscanf(fp5,"%lf",&X_old[Ncomp+inode*Nunk_per_node+iunk]);
      /*  X_old[Ncomp+inode*Nunk_per_node+iunk]=
	  (X_old[Ncomp+inode*Nunk_per_node+iunk]-0.45)*0.025 + 0.45;
       */
	}
	else if (iunk < 2*Ncomp){
	  fscanf(fp5,"%lf",&tmp);
          if (Lbinodal && iguess==BINODAL_FLAG)
	      X2_old[inode*Nunk_per_node+iunk-Ncomp] = exp(-tmp);
          else{
	      X_old[inode*Nunk_per_node+iunk-Ncomp] = exp(-tmp);
          }
	}
	else if (iunk<2*Ncomp+Ngeqn_tot){
          ipol=Unk_to_Poly[iunk-2*Ncomp];
          iseg=Unk_to_Seg[iunk-2*Ncomp];
          itype_mer=Type_mer[ipol][iseg];

	  if (Type_poly == 2){
	      fscanf(fp6,"%lf",&tmp);
              if (Lbinodal && iguess==BINODAL_FLAG)
	          X2_old[inode*Nunk_per_node+iunk]=tmp/X2_old[inode*Nunk_per_node+itype_mer];
              else
	          X_old[inode*Nunk_per_node+iunk]=tmp/X_old[inode*Nunk_per_node+itype_mer];
	  }
	  else {
            if (Lbinodal && iguess==BINODAL_FLAG)
	          fscanf(fp6,"%lf",&X2_old[inode*Nunk_per_node+iunk]);
	    else  fscanf(fp6,"%lf",&X_old[inode*Nunk_per_node+iunk]);
	  }    
	}
        else{    /* electric field */
	  fscanf(fp5,"%lf",&tmp);
          if (Lbinodal && iguess==BINODAL_FLAG)
	      X2_old[inode*Nunk_per_node+iunk] = tmp;
          else
	      X_old[inode*Nunk_per_node+iunk] = tmp;
        } 
      }
      else
	if (Lbinodal && iguess==BINODAL_FLAG) { 
	  fscanf(fp5,"%lf",&X2_old[inode*Nunk_per_node+iunk]);
	}
	else { 
	  fscanf(fp5,"%lf",&X_old[inode*Nunk_per_node+iunk]);
	}
      
      /* if (iunk<Ncomp && Lsteady_state){
	 if (Lbinodal && iguess==BINODAL_FLAG) 
	 X2_old[inode*Nunk_per_node+iunk] *= Mass[iunk];
	 else
	 X_old[inode*Nunk_per_node+iunk] *= Mass[iunk];
	 }
      */

      
      if (Lsteady_state){
	if ( (Ipot_ff_c == COULOMB && iunk > Ncomp) ||
	     (Ipot_ff_c != COULOMB && iunk>Ncomp-1)){
	  if (Ipot_ff_c == COULOMB) i=iunk-Ncomp-1;
	  else i=iunk-Ncomp;
	  
     /*   if (Ipot_ff_n != IDEAL_GAS){
	  if (Lbinodal && iguess==BINODAL_FLAG) 
	  X2_old[inode*Nunk_per_node+iunk] -=
	  ( 3.0*log(Sigma_ff[i][i]) + 1.5*log(Mass[i]*Temp)  );
	  else
	  X_old[inode*Nunk_per_node+iunk] -=
	  ( 3.0*log(Sigma_ff[i][i]) + 1.5*log(Mass[i]*Temp)  );
	  }
      */
	}
      }
      
    }
    while ((c=getc(fp5)) != EOF && c !='\n') ;
    if (Restart==5 && ijk_old[0]==Nodes_x[0]-1){
       fclose(fp5);
       if (Sten_Type[POLYMER_CR]) fclose(fp6);
       open_now=TRUE;
    }
  }
  for (idim=0; idim<Ndim; idim++) {
    Nodes_x_old[idim] = ijk_old_max[idim] + 1;
  }
  if (Restart!=5){
       fclose(fp5);
       if (Sten_Type[POLYMER_CR]) fclose(fp6);
  }
  return;
}
/*******************************************************************/
/*shift_the_profile: do this if the new mesh and the old mesh
                   have identical Esize, but not identical Nnodes_x */
void shift_the_profile(double *x_new,double fac)
{

  int idim,jdim,iunk,inode,inode_old,ijk[3],ijk_tmp[3],Nadd;
  double x_test,unk_old,unk_1,unk_2;

  idim = Plane_new_nodes;
  Nadd = round_to_int(Del_1[idim]/Esize_x[idim]);

  for (inode=0; inode<Nnodes; inode++){
     node_to_ijk(inode,ijk);
     for (iunk=0; iunk<Nunk_per_node; iunk++){

     switch(Pos_new_nodes){
         case  0:          /*ADDING NODES TO CENTER OF BOX */

           if (ijk[idim] < Nodes_x_old[idim]/2){           /*NODE LBB OF NEW PLANE*/
             inode_old = locate_inode_old(ijk);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else if (ijk[idim] > Nodes_x_old[idim]/2+Nadd){ /*NODE RTF OF NEW PLANE*/
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] -= Nadd;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else {                                    /*NODE IN CENTER OF NEW PLANE*/
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = Nodes_x_old[idim]/2;
             inode_old = locate_inode_old(ijk_tmp);
             unk_1 = X_old[inode_old*Nunk_per_node+iunk];

             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = Nodes_x_old[idim]/2 + 1;
             inode_old = locate_inode_old(ijk_tmp);
             unk_2 = X_old[inode_old*Nunk_per_node+iunk];
      
             unk_old = 0.5*(unk_1+unk_2);
           }
           break;

         case -1:          /*ADDING NODES TO LEFT,BACK,BOTTOM*/
           if (ijk[idim] < Nadd){
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = 0;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else{
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] -= Nadd;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           break;

         case  1:          /*ADDING NODES TO RIGHT,TOP,FRONT*/
           if (ijk[idim] < Nodes_x_old[idim]){
             inode_old = locate_inode_old(ijk);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           else {
             for (jdim=0; jdim<Ndim; jdim++) ijk_tmp[jdim] = ijk[jdim];
             ijk_tmp[idim] = Nodes_x_old[idim]-1;
             inode_old = locate_inode_old(ijk_tmp);
             unk_old = X_old[inode_old*Nunk_per_node+iunk];
           }
           break;

         default:
           printf("Check pos new nodes %d [should be -1,0,1]\n",Pos_new_nodes);
           exit(-1);
           break;
 
     }

     /* check a few limiting values ... and finally set initial guess*/
     if (iunk <Ncomp )
        x_test = AZ_MIN(Rho_max,fac*(unk_old-Rho_b[iunk])+ Rho_b[iunk]);

     else if (iunk<Ncomp+Nrho_bar){
        x_test = fac*(unk_old-Rhobar_b[iunk-Ncomp])+ Rhobar_b[iunk-Ncomp];
        if (iunk == Ncomp && x_test >= 1.0) x_test = Rhobar_b[iunk-Ncomp];
     }
     else if (iunk == Ncomp+Nrho_bar && Ipot_ff_c == COULOMB)
        x_test = fac*(unk_old-1.0) + 1.0;
     else 
        x_test = unk_old;

     x_new[inode*Nunk_per_node+iunk] = x_test;

     }
  }
  return;
}
/*********************************************************************/
/* locate_inode_old: find a given node in the old solution array */
int locate_inode_old(int *ijk)
{
   int inode_old=0;
 
   if (Ndim == 1)      inode_old = ijk[0];
   else if (Ndim == 2) inode_old = ijk[0] + ijk[1] * Nodes_x_old[0];
   else if (Ndim == 3) inode_old = ijk[0] + ijk[1] * Nodes_x_old[0] 
                       + ijk[2] * Nodes_x_old[0]*Nodes_x_old[1];
   return(inode_old);
}
/**********************************************************************/
/* communicate profile: broadcast the x_new profile to all processors,
                        and let each of them pick out the needed entries.*/
void communicate_profile(double *x_new,double *x)
{
   int loc_inode,inode,ijk[3],iunk,loc_i;   
   
    AZ_broadcast((char *) x_new, Nnodes*Nunk_per_node*sizeof(double), 
                                 Aztec.proc_config, AZ_PACK);
    AZ_broadcast(NULL, 0, Aztec.proc_config, AZ_SEND);

    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       inode = L2G_node[loc_inode];
       node_to_ijk(inode,ijk);
       for (iunk=0; iunk<Nunk_per_node; iunk++){
           loc_i = Aztec.update_index[iunk + Nunk_per_node * loc_inode];
           x[loc_i] = x_new[inode*Nunk_per_node+iunk];
       }
    }
    return;
}
/*********************************************************************/
/*check_zero_densities: here just remove zero densities where 
         not appropriate, and make sure x=0 where needed */
void check_zero_densities(double *x)
{

  int loc_inode,zero_all,icomp,inode_box,loc_i,loc_j=0;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      zero_all = TRUE;
      for (icomp=0; icomp<Ncomp; icomp++){
         inode_box = L2B_node[loc_inode];
         loc_i = Aztec.update_index[icomp + Nunk_per_node * loc_inode];
         if (Sten_Type[POLYMER_CR]) loc_j = Aztec.update_index[icomp + Ncomp+ Nunk_per_node * loc_inode];
         if (Zero_density_TF[inode_box][icomp] == TRUE){  
               x[loc_i] = 0.0;
               if (Sten_Type[POLYMER_CR]) x[loc_j]=0.0;
         }
         else if (Zero_density_TF[inode_box][icomp] == FALSE && !Sten_Type[POLYMER_CR]){
             zero_all = FALSE;
             if (x[loc_i]</*DENSITY_MIN*/Rho_b[icomp]*exp(-VEXT_MAX)) {
                  x[loc_i] = /*DENSITY_MIN*/Rho_b[icomp]*exp(-VEXT_MAX);
             }
         }

      }

      if (zero_all == TRUE && Nunk_per_node > Ncomp+Nrho_bar 
              && !Sten_Type[POLYMER_CR] && Type_coul==-1){
        loc_i = Aztec.update_index[Ncomp+Nrho_bar + Nunk_per_node * loc_inode];
        x[loc_i] = 1.0; 
      } 
  }
  return;
}
/**********************************************************************/
/*chop_profile: do the profile chop here. */
void chop_profile(double *x, int iguess)
{
  int loc_inode,icomp,inode_box,iwall,check,loc_i;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++)
    for (icomp=0; icomp<Ncomp; icomp++) {
        inode_box = L2B_node[loc_inode];
        check = 0;
        for (iwall=0; iwall<Nwall; iwall++)
              if (X_wall[loc_inode][iwall] > Xstart_step[0]) check++;
        if (check == Nwall){
              loc_i = Aztec.update_index[icomp+Nunk_per_node*loc_inode];
              if (iguess==CHOP_RHO_L) x[loc_i] = Rho_coex[1];
              else if (iguess==CHOP_RHO_V) x[loc_i] = Rho_coex[0];
              else if (iguess==CHOP_RHO) x[loc_i] = Rho_b[icomp];
              else if (iguess==CHOP_RHO_STEP) x[loc_i] = Rho_step[icomp][0];
        }
    }
  return;
}
/****************************************************************************/
/* interpolate: given a set of densities for an arbitrary number of
                   components, interpolate to find a density at a new
                   location.  For the moment, we assume that the old and new
                   domains are the same size.                               */

double interpolate(double *pos,int icomp, int nodes_old,
                           double **pos_old, double *Esize_old)
{
  /* define local variables */
  int idim,i,imin=0,npoints,match;
  double delx[3],delx_min[3],r12_sq,r12_sq_min,pos_test;
  double x_upper[3],x_lower[3],x_corner[8][3],func[8],fac[3];
  double Lmax[3],interp=0.0;

  /* find the grid spacing in each direction in the old data set */
 
  for(idim=0;idim<Ndim;idim++){
      Lmax[idim] = 0.0;
      for (i=0; i<nodes_old; i++)
         if (pos_old[i][idim] > Lmax[idim]) Lmax[idim]=pos_old[i][idim];
  }

  /* find the old node that is closest to the proposed point */
  r12_sq_min = 10000.;
  for (i=0; i<nodes_old; i++) {

     r12_sq = 0.0;
     for (idim=0; idim<Ndim; idim++){
          delx[idim] = pos_old[i][idim] - pos[idim];
          r12_sq = r12_sq + delx[idim]*delx[idim];
     }

     if (r12_sq < r12_sq_min) {
        r12_sq_min=r12_sq;
        for (idim=0; idim<Ndim; idim++) delx_min[idim]=delx[idim];
        imin = i;
     }
  }

  /* find the positions of the old nodes that surround the proposed point */
  for (idim=0; idim<Ndim; idim++){
     if (delx_min[idim] > 0.0) {
         x_upper[idim] = pos_old[imin][idim];
         x_lower[idim] = pos_old[imin][idim] - Esize_old[idim];
     }
     else if (delx_min[idim] < 0.0) {
         x_upper[idim] = pos_old[imin][idim] + Esize_old[idim];
         x_lower[idim] = pos_old[imin][idim];
     }
     else {
         x_upper[idim] = pos_old[imin][idim] + Esize_old[idim];
         x_lower[idim] = pos_old[imin][idim];
         if (x_upper[idim] > Lmax[idim]) {
            x_upper[idim] = pos_old[imin][idim];
            x_lower[idim] = pos_old[imin][idim] - Esize_old[idim];
         }
     } 
  }

  if (Ndim==1){
     x_corner[0][0] = x_lower[0];
     x_corner[1][0] = x_upper[0];
  }
  else if (Ndim == 2){
     x_corner[0][0] = x_lower[0];
     x_corner[0][1] = x_lower[1];

     x_corner[1][0] = x_upper[0];
     x_corner[1][1] = x_lower[1];

     x_corner[2][0] = x_upper[0];
     x_corner[2][1] = x_upper[1];

     x_corner[3][0] = x_lower[0];
     x_corner[3][1] = x_upper[1];
  }
  else {
     x_corner[0][0] = x_lower[0];
     x_corner[0][1] = x_lower[1];
     x_corner[0][2] = x_lower[2];

     x_corner[1][0] = x_upper[0];
     x_corner[1][1] = x_lower[1];
     x_corner[1][2] = x_lower[2];

     x_corner[2][0] = x_upper[0];
     x_corner[2][1] = x_upper[1];
     x_corner[2][2] = x_lower[2];

     x_corner[3][0] = x_lower[0];
     x_corner[3][1] = x_upper[1];
     x_corner[3][2] = x_lower[2];

     x_corner[4][0] = x_lower[0];
     x_corner[4][1] = x_lower[1];
     x_corner[4][2] = x_upper[2];

     x_corner[5][0] = x_upper[0];
     x_corner[5][1] = x_lower[1];
     x_corner[5][2] = x_upper[2];

     x_corner[6][0] = x_upper[0];
     x_corner[6][1] = x_upper[1];
     x_corner[6][2] = x_upper[2];

     x_corner[7][0] = x_lower[0];
     x_corner[7][1] = x_upper[1];
     x_corner[7][2] = x_upper[2];
  }

  /* find the old values of f(x) at the corner nodes */
  for (npoints = 0; npoints<Nnodes_per_el_V; npoints++) 
                                      func[npoints] = -1.0;

  for (npoints = 0; npoints<Nnodes_per_el_V; npoints++) {
     for (i=0; i<nodes_old; i++) {
        match = 0;
        for (idim=0; idim<Ndim; idim++) {
           if (Type_bc[idim][1]==PERIODIC && x_corner[npoints][idim] > Lmax[idim]) 
                    pos_test = pos_old[i][idim] + Size_x[idim];
           else     pos_test = pos_old[i][idim];
           if (fabs(x_corner[npoints][idim] - pos_test) < 1.e-8) match++;
        }
        if (match == Ndim) func[npoints] = X_old[i*Nunk_per_node+icomp];
     }
  }

  for (npoints = 0; npoints<Nnodes_per_el_V; npoints++) {
     if (func[npoints] == -1.0) {
        if      (npoints == 0 && func[1] != -1.0) func[0] = func[1];
        else if (npoints == 1 && func[0] != -1.0) func[1] = func[0];
        else if (npoints == 2 && func[3] != -1.0) func[2] = func[3];
        else if (npoints == 3 && func[2] != -1.0) func[3] = func[2];
        else if (npoints == 4 && func[5] != -1.0) func[4] = func[5];
        else if (npoints == 5 && func[4] != -1.0) func[5] = func[4];
        else if (npoints == 6 && func[7] != -1.0) func[6] = func[7];
        else if (npoints == 7 && func[6] != -1.0) func[7] = func[6];
        else func[npoints] = 1.0e-5;
/*        else{
           printf("ERROR in interpolate: corner # %d has a value of -1.0\n",npoints);
           printf("Its dimensions are:\n");
           for (idim=0;idim<Ndim; idim++) 
              printf("x%d : %f",idim,x_corner[npoints][idim]);
           exit(-1);
        } */
     }
  }


  /* calculate interpolation factors */
  for (idim=0; idim<Ndim; idim++) 
     fac[idim] = (pos[idim]-x_lower[idim])/(x_upper[idim]-x_lower[idim]);

  /* calculate the density at the proposed point */
  if (Ndim == 1) 
     interp =  (1.0-fac[0])*func[0] + fac[0]*func[1];
  else if (Ndim == 2) 
     interp = (1.0-fac[0])*(1.0-fac[1])*func[0] + fac[0]*(1.0-fac[1])*func[1]
                + fac[1]*(1.0-fac[0])*func[3] + fac[1]*fac[0]*func[2];
  else if (Ndim == 3) {
     printf("ERROR:Haven't coded 3D interpolation yet");
     exit(-1);
  }
  return interp;
}
/****************************************************************************/
