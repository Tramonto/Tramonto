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
#include <string.h>
 
static void chop_profile(double **x, int);
static void check_zero_densities(double **);
static void communicate_profile(double *,double **);
static void shift_the_profile(double *,double);
static int  find_length_of_file(char *);
static void read_in_a_file(int,char *);
static void setup_step_2consts(double **);
static void setup_linear_profile(double **);
static void setup_exp_density(double **, double *,int,int);
static void setup_const_density(double **, double *,int,int);
static void setup_stepped_profile(double **);
static void setup_elec_pot(double **,int);
static void setup_chem_pot(double **);
static void setup_polymer(double **);
static void setup_polymer_rho(double **,int);
static void setup_polymer_field(double **,int);
static void setup_polymer_simple(double **, int);
static void setup_polymer_G(double **);
static void setup_polymer_G_2(double **);
static void setup_rho_bar(double **);
static void setup_exp_density_with_profile(double **);
static int locate_inode_old(int *);
static void setup_BondWTC(double **);
static void setup_Xi_cavWTC(double **);

void set_initial_guess (int iguess, double** xOwned)
{
 /*
  * Local variable declarations
  */

  char filename[20];
  double t1=0.0;
  double *x_new;

  int iunk,i;
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
    if (Restart==0 && Imain_loop==0) printf("generating guess from scratch\n");  
    else printf("setting up guess from existing files\n");
  } 

  if (Restart > 0 || Imain_loop > 0){
       start_no_info = FALSE;
       for (i=0;i<NEQ_TYPE;i++) Restart_field[i]=FALSE;

        x_new = (double *) array_alloc(1, Nnodes*Nunk_per_node, sizeof(double));

        if (Proc == 0) {  /* Proc 0 reads in the data file */

           if ( Imain_loop == 0){

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
              if (Restart==5){ /* going to read in 1D file, need to dimension arrays for 2D or 3D system */
                  Nodes_old=Nnodes;  
if (Proc==0 && Iwrite != NO_SCREEN) printf("set Nodes_old to be Nnodes...Nnodes=%d Nodes_old=%d\n",Nnodes,Nodes_old);
              }
              if (Lbinodal && iguess==BINODAL_FLAG){
                 X2_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
              }
              else
                 X_old = (double *) array_alloc(1, Nodes_old*Nunk_per_node, sizeof(double));
              read_in_a_file(iguess,filename); /* Get X_old */
           }
           else{
              for (iunk=0;iunk<Nunk_per_node;iunk++) Restart_field[Unk2Phys[iunk]]=TRUE;
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

        MPI_Bcast (Restart_field,NEQ_TYPE,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast (&start_no_info,1,MPI_INT,0,MPI_COMM_WORLD);

        if (!start_no_info){
            if (Restart_field[DENSITY]==FALSE) {
                printf("can't restart without density fields yet\n");
                exit(-1);
            }
            communicate_profile(x_new,xOwned);
            check_zero_densities(xOwned);
            if (Lsteady_state && (Restart==3 || Restart_field[DIFFUSION]==FALSE))   setup_chem_pot(xOwned);
            if ((Type_poly==NONE||Type_poly==WTC) &&
                  (Restart==3 || Restart_field[RHOBAR_ROSEN]==FALSE)) setup_rho_bar(xOwned);
            if (Ipot_ff_c == COULOMB && (Restart==3 || Restart_field[POISSON]==FALSE)){
                   printf("setting up electrostatic potential guess....\n");
                   setup_elec_pot(xOwned,iguess);
            }
            if ((Type_poly != NONE && Type_poly!=WTC) && Restart_field[CMS_FIELD]==FALSE) setup_polymer_field(xOwned,iguess);   
            if (Restart == 2) chop_profile(xOwned,iguess);
        /*    if (Restart ==6 && Nwall>0) setup_exp_density_with_profile(xOwned);*/
        }
        safe_free((void *) &x_new);


  } /* end of logic for restarting from a file or previous run */

  else start_no_info = TRUE;

  if (start_no_info) {  /* START FROM A SPECIFIED INITIAL GUESS */

 if ((Type_poly==NONE || Type_poly==WTC)) {
    switch(iguess){
      case CONST_RHO:    
            if (Type_poly==WTC){
                 setup_const_density(xOwned,Rho_seg_b,Nseg_tot,0);
            }
            else              setup_const_density(xOwned,Rho_b,Ncomp,0);
            break;

      case CONST_RHO_V:  
            setup_const_density(xOwned,Rho_coex,1,1);
            break;
      case CONST_RHO_L:  
            setup_const_density(xOwned,Rho_coex,1,0);
            break;

      case EXP_RHO:
            setup_exp_density(xOwned,Rho_b,Ncomp,0);
            break;
      case EXP_RHO_V:
            setup_exp_density(xOwned,Rho_coex,1,1);
            break;
      case EXP_RHO_L:
            setup_exp_density(xOwned,Rho_coex,1,0);
            break;

      case STEP_PROFILE:
            setup_stepped_profile(xOwned);

      case CHOP_RHO_L:
            setup_exp_density(xOwned,Rho_coex,1,1);
            chop_profile(xOwned,iguess);
            break;
      case CHOP_RHO_V:
            setup_exp_density(xOwned,Rho_coex,1,0);
            chop_profile(xOwned,iguess);
            break;
      case LINEAR:
            setup_linear_profile(xOwned);
    }  /* end of iguess switch */

    setup_rho_bar(xOwned);
    if (Lsteady_state)   setup_chem_pot(xOwned);
    if (Type_poly==WTC){ 
         setup_BondWTC(xOwned);
         setup_Xi_cavWTC(xOwned);
    }
 }
 else{
    /*setup_polymer(x);*/

   if (Type_poly == CMS_SCFT)
     setup_polymer_simple(xOwned,iguess);
   else
     setup_polymer_rho(xOwned,iguess);
     setup_polymer_field(xOwned,iguess);
     setup_polymer_G(xOwned);
    /* setup_polymer_G_2(xOwned); */

    check_zero_densities(xOwned);
 }

 if (Ipot_ff_c == COULOMB)  setup_elec_pot(xOwned,iguess);

   } /* end of cases for setting initial guess from scratch */

  if (Proc==0 && Iwrite==VERBOSE) printf("\n initial guess took %g secs\n",MPI_Wtime()-t1);
  return;
}
/*********************************************************/
/*setup_polymer_field: in this routine set up the field guesses
                       for the polymers variables    */
static void setup_polymer_field(double **xOwned, int iguess)
{
  int loc_inode,itype_mer,irho, iunk;
  double field;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){

     for (itype_mer=0; itype_mer<Ncomp; itype_mer++){
	 irho = Phys2Unk_first[DENSITY]+itype_mer;
	 iunk = Phys2Unk_first[CMS_FIELD]+itype_mer;
         if (xOwned[irho][loc_inode]<1.e-6) field=VEXT_MAX-1.;
         else field=-log(xOwned[irho][loc_inode]/Rho_b[itype_mer]);
         xOwned[iunk][loc_inode]=exp(-field);
     }
   }
   return;
}
/*********************************************************/
/*setup_polymer_simple: in this routine set up the field guesses
                       for the polymers variables for SCF case    */
static void setup_polymer_simple(double **xOwned, int iguess)
{
  int loc_inode,inode_box,ijk_box[3],i, iunk,junk;
  double temp;
  int itype_mer,jtype_mer,inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);

     for (itype_mer=0; itype_mer<Ncomp; itype_mer++){
	iunk = Phys2Unk_first[CMS_FIELD]+itype_mer;
        if (!Zero_density_TF[inode_box][itype_mer]){
           temp = 0.;
           for (jtype_mer=0; jtype_mer<Ncomp; jtype_mer++){
	     junk = Phys2Unk_first[CMS_FIELD]+jtype_mer;
           
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
           if ((Sten_Type[POLYMER_CR]) || (temp > 0.5)) xOwned[iunk][loc_inode] = exp(temp);
           else                             xOwned[iunk][loc_inode] = exp(1. - sqrt(1.-2.*temp));
        } /* end if i not zero density  */
        else   xOwned[iunk][loc_inode] = DENSITY_MIN; /* zero density - Boltzmann probability = 0 */
     } /* end loop over itype_mer */
  } /* end loop over loc_inode */
  return;
}
/*********************************************************/
/*setup_polymer_rho: in this routine set up polymer density profiles    */
static void setup_polymer_rho(double **xOwned, int iguess)
{
  int loc_inode,i,inode_box,ijk_box[3],iunk,icomp;
  int inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     node_box_to_ijk_box(inode_box, ijk_box);
     if (iguess == CONST_RHO) {
       for (icomp=0; icomp<Ncomp; icomp++){
	 iunk = Phys2Unk_first[DENSITY]+icomp;
         if (!Zero_density_TF[inode_box][icomp]) xOwned[iunk][loc_inode] = Rho_b[icomp];
         else                                    xOwned[iunk][loc_inode] = 0.0;
       }
     }
     else if (iguess == STEP_PROFILE) {
       inode = B2G_node[inode_box]; 
       node_to_position(inode,nodepos);
       for (i=0;i<Nsteps;i++){
           if (nodepos[Orientation_step[i]]>=Xstart_step[i] &&
                nodepos[Orientation_step[i]]<=Xend_step[i]){
                for (icomp=0;icomp<Ncomp;icomp++){
	           iunk = Phys2Unk_first[DENSITY]+icomp;
       		   if (!Zero_density_TF[inode_box][icomp]) xOwned[iunk][loc_inode]= Rho_step[icomp][i];
		   else xOwned[iunk][loc_inode]=0.0;
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
static void setup_polymer_G(double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],loc_i;
  int reflect_flag[NDIM_MAX];
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  int sten_type,izone,jlist,jnode_box,jtype_mer,itype_mer;
  int iunk,poln,iseg,ibond,not_done,junk,cycle,loc_B;

     if (Sten_Type[POLYMER_GAUSS]) sten_type = POLYMER_GAUSS;
     else                          sten_type = DELTA_FN;
     izone = 0;

     loc_inode=0;
     for (poln=0; poln < Npol_comp; poln++){
        iunk = Geqn_start[poln];
        for (iseg=0; iseg<Nmer[poln]; iseg++){
             for (ibond=0; ibond<Nbond[poln][iseg]; ibond++){
                 xOwned[iunk][loc_inode] =999.0;
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
             if (fabs(xOwned[iunk][0]-999.0)<1.e-6){

             if(Bonds[poln][iseg][ibond]== -1){
                 for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
                     inode_box = L2B_node[loc_inode];
                     node_box_to_ijk_box(inode_box, ijk_box);
                     xOwned[iunk][loc_inode] = xOwned[Phys2Unk_first[CMS_FIELD]+itype_mer][loc_inode];
                  }
             }
             else{
               jtype_mer = Type_mer[poln][Bonds[poln][iseg][ibond]];
               junk = Geqn_start[poln]+2*Bonds[poln][iseg][ibond];
               if (iseg<Bonds[poln][iseg][ibond]) junk += 1;
               /* test if this G equation has been generated yet ... if not, go on */
               if (fabs(xOwned[junk][0]-999.0)>1.e-6){
                   for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
                       inode_box = L2B_node[loc_inode];
                       node_box_to_ijk_box(inode_box, ijk_box);
                       xOwned[iunk][loc_inode] = 0.;

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
                              /*xOwned[iunk][loc_inode] +=  weight * xOwned[junk][B2L_node[jnode_box]]; */
                              xOwned[iunk][loc_inode] +=  weight; 
                           } /* check that node is in domain */
                           else{  /* use the value at loc_inode ... an approximation */
                              xOwned[iunk][loc_inode] +=  weight*xOwned[junk][loc_inode]; 
                           }
                         }
                      }

                  xOwned[iunk][loc_inode] *= xOwned[Phys2Unk_first[CMS_FIELD]+itype_mer][loc_inode];

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
                 if (fabs(xOwned[iunk][0]-999.0)<1.e-6) not_done=TRUE;
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
#ifdef OLDCODE
/* THIS ROUTINE DOES NOT WORK AS WRITTEN IN PARALLEL BECAUSE THE
 * INTEGRATION CAN ACCESS ELEMENTS OFF PROCESSOR. THE OFFENDING
 * COMMANDS ARE   xOwned[][B2L_node()] AND WOULD HAVE A -1 ARGUMENT
 * FOR THE SECOND INDEX. NOT FIXING RIGHT NOW BECAUSE THE CALL TO
 * THIS ROUTINE IS COMMENTED OUT ANYWAY.  AGS 4/05 */
static void setup_polymer_G_2(double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],loc_i;
  int reflect_flag[NDIM_MAX];
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  int sten_type,izone,jlist,jnode_box,itype_mer,jtype_mer;
  int iunk,poln,iseg,ibond,j_box,unk_B;

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
	  xOwned[iunk][loc_inode] = 0.;

                             /* TREAT THE END SEGMENTS */
	  if(Bonds[poln][iseg][ibond]== -1){
            unk_B = Phys2Unk_first[CMS_FIELD] + Type_mer[poln][iseg];
	    xOwned[iunk][loc_inode] = xOwned[unk_B][loc_inode];
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
		  xOwned[iunk][loc_inode] +=  weight*xOwned[iunk][B2L_node[jnode_box]]; 
		}
		else if (jnode_box == -1 )
		  xOwned[iunk][loc_inode] += weight;
	      } /* end loop over stencil */

	      /* Boltzmann factor */
	      xOwned[iunk][loc_inode] *= xOwned[Phys2Unk_first[CMS_FIELD] + itype_mer][loc_inode];

	      /* } end of Zero_dens_TF test */
	  } /* end of if Bond test */
	  iunk++;
        } /* end of loop over each segment on poln */

	/* now do G's in backward direction */
        iunk = Ngeqn_tot-1;
	ibond = 1;  /* only doing one set of bonds at a time */
        for (iseg=Nmer[poln]-1; iseg>-1; iseg--){
	  itype_mer = Type_mer[poln][iseg];
	  xOwned[Phys2Unk_first[CMS_G]+iunk][loc_inode] = 0.;

                             /* TREAT THE END SEGMENTS */
	  if(Bonds[poln][iseg][ibond]== -1){
	    xOwned[Phys2Unk_first[CMS_G]+iunk][loc_inode] = xOwned[Phys2Unk_first[CMS_FIELD]+Type_mer[poln][iseg]][loc_inode];
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
		  xOwned[Phys2Unk_first[CMS_G]+iunk][loc_inode] +=  weight*xOwned[iunk][B2L_node[jnode_box]]; 
		}
		else if (jnode_box == -1 )
		  xOwned[Phys2Unk_first[CMS_G]+iunk][loc_inode] += weight;
	      } /* end loop over stencil */

	      /* Boltzmann factor */
	      xOwned[Phys2Unk_first[CMS_G]+iunk][loc_inode] *= xOwned[Phys2Unk_first[CMS_FIELD]+itype_mer][loc_inode];

	    } /*end of Zero_dens_TF test */
	  } /* end of if Bond test */
	  iunk--;
	  MPI_Barrier(MPI_COMM_WORLD);
        } /* end of loop over each segment on poln */

     } /* end of loop over polymer chains */
  }  /* end of loop over nodes */

  return;
}
#endif

/*********************************************************/
/*setup_const_density: in this routine set up a constant
        density profile wherever Zero_density_TF = FALSE */
static void setup_const_density(double **xOwned, double *rho,int nloop,int index)
{
  int loc_inode,i,inode_box,iunk,zeroTF;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0; i<nloop; i++){
	 iunk = i+Phys2Unk_first[DENSITY];
         if (Type_poly==WTC) zeroTF=Zero_density_TF[inode_box][Unk2Comp[i]];
         else                zeroTF=Zero_density_TF[inode_box][i];
         if (!zeroTF){
            if (nloop > 1) xOwned[iunk][loc_inode] = rho[i];
            else           xOwned[iunk][loc_inode] = rho[index];
         }
         else xOwned[iunk][loc_inode] = 0.0;

         /* set up initial guess if chemical potential is an unknown */
         if (Lsteady_state) {
	     iunk = i+Phys2Unk_first[DIFFUSION];
             if (!zeroTF){
                if (nloop > 1) xOwned[iunk][loc_inode] = Betamu[i];
                else           xOwned[iunk][loc_inode] = Betamu[index];
             }
             else xOwned[iunk][loc_inode] = -10.0; /* zero density == -infinite chem.pot*/
         }
      }
  }
  return;
}
/*********************************************************/
/*setup_stepped_profile: in this routine set up a stepped
        density profile wherever Zero_density_TF = FALSE */
static void setup_stepped_profile(double **xOwned)
{
  int loc_inode,i,inode_box,iunk,icomp,inode;
  double nodepos[3];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode = L2G_node[loc_inode];
     node_to_position(inode,nodepos);

     for (i=0;i<Nsteps;i++){
       if (nodepos[Orientation_step[i]] >= Xstart_step[i] &&
           nodepos[Orientation_step[i]] <= Xend_step[i]){
           for (icomp=0; icomp<Ncomp; icomp++){
	       iunk = Phys2Unk_first[DENSITY]+icomp;
               if (!Zero_density_TF[inode_box][icomp]){
                   xOwned[iunk][loc_inode]=Rho_step[icomp][i];
               }
               else xOwned[iunk][loc_inode]=0.0;
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
/*setup_exp_density_with_profile: in this routine set up a density 
                     profile as rho(x)*exp(-Vext/kT)*/
static void setup_exp_density_with_profile(double **xOwned)
{

  int loc_inode,i,inode_box,iunk,icomp;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0;i<Ncomp;i++){
        iunk = icomp+Phys2Unk_first[DENSITY];
        if (Vext[loc_inode][icomp]>0.0) xOwned[iunk][loc_inode] *= exp(-Vext[loc_inode][i]);

/*        if (Type_poly!=NONE && Type_poly!=WTC)
        xOwned[i+Phys2Unk_first[CMS_FIELD]][loc_inode] = exp(-log(xOwned[i+Phys2Unk_first[CMS_FIELD]][loc_inode])+Vext[loc_inode][i]);*/
     }
  }

  return;
}
/*********************************************************/
/*setup_exp_density: in this routine set up a density 
                     profile as rho_b*exp(-Vext/kT)*/
static void setup_exp_density(double **xOwned, double *rho,int nloop,int index)
{

  int loc_inode,i,inode_box,iunk;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     for (i=0; i<nloop; i++) {

	iunk = i+Phys2Unk_first[DENSITY];
        if (!Zero_density_TF[inode_box][i]){
            if (nloop > 1) xOwned[iunk][loc_inode] = rho[i] * exp(-Vext[loc_inode][i]);
            else           xOwned[iunk][loc_inode] = rho[index] * exp(-Vext[loc_inode][i]);
        }
        else xOwned[iunk][loc_inode] = 0.0;

        if (xOwned[iunk][loc_inode]>1.0) xOwned[iunk][loc_inode] = 0.99;   /* may need this to make sure that rb3<1 always */


         /* set up initial guess if chemical potential is an unknown */
         if (Lsteady_state) {
	    iunk = i+Phys2Unk_first[DIFFUSION];
             if (!Zero_density_TF[inode_box][i]){
                if (nloop > 1) xOwned[iunk][loc_inode] = Betamu[i];
                else           xOwned[iunk][loc_inode] = Betamu[index];
             }
             else xOwned[iunk][loc_inode] = -10.0; /* zero density == -infinite chem.pot*/
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

static void setup_step_2consts(double **xOwned)
{

  int loc_inode,icomp,inode_box,ijk[3],inode;
  int iunk;
  double x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0] - 0.5*Size_x[0];

     for (icomp=0; icomp<Ncomp; icomp++) {

        iunk = icomp+Phys2Unk_first[DENSITY];
        if (!Zero_density_TF[inode_box][icomp]){
           if (x_dist< 0.0)  xOwned[iunk][loc_inode] = constant_boundary(iunk,-3);
           else              xOwned[iunk][loc_inode] = constant_boundary(iunk,-4);
        }
        else xOwned[iunk][loc_inode] = 0.0;

       /* set up initial guess if chemical potential is an unknown */
        if (Lsteady_state) {
            iunk = icomp+Phys2Unk_first[DIFFUSION];
            if (!Zero_density_TF[inode_box][icomp]){
               if (x_dist< 0.0) xOwned[iunk][loc_inode] = constant_boundary(iunk,-3);
               else             xOwned[iunk][loc_inode] = constant_boundary(iunk,-4);
            }
            else xOwned[iunk][loc_inode] = -10.0;  /*zero density == -infinite chem.pot*/
        }
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

static void setup_linear_profile(double **xOwned)
{

  int loc_inode,icomp,inode_box,ijk[3],inode;
  int iunk;
  double x_dist,x_tot,rho_LBB, rho_RTF;

  if (Lsteady_state) x_tot=Size_x[Grad_dim]-2.0*X_const_mu;
  else x_tot = Size_x[Grad_dim];

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode_box = L2B_node[loc_inode];
     inode     = B2G_node[inode_box];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[0]*ijk[0]-X_const_mu;

     for (icomp=0; icomp<Ncomp; icomp++) {

        iunk = Phys2Unk_first[DENSITY]+icomp;
        if (!Zero_density_TF[inode_box][icomp]){
           rho_LBB = constant_boundary(iunk,-3);
           rho_RTF = constant_boundary(iunk,-4);

           if (x_dist <0.) xOwned[iunk][loc_inode] = rho_LBB;
           else if (x_dist > x_tot) xOwned[iunk][loc_inode] = rho_RTF;
           else  xOwned[iunk][loc_inode] = rho_LBB + (rho_RTF-rho_LBB)*x_dist/x_tot;
        }
        else xOwned[iunk][loc_inode] = 0.0;

     }
  }
  return;
}
/************************************************************/
/*setup_rho_bar: set up the rhobar initial guesses.  For now
                use Rho_b to set initial guess.  Later calculate
                based on rho initial guess. */
static void setup_rho_bar(double **xOwned)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,irb;
  double vol,area,x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (irb = 0; irb < Nrho_bar; irb++){
       iunk = Phys2Unk_first[RHOBAR_ROSEN] + irb;
       if (Lsteady_state || (Nwall == 0 && Iliq_vap == 3)){
           inode_box = L2B_node[loc_inode];
           inode     = B2G_node[inode_box];
           node_to_ijk(inode,ijk); 
           x_dist = Esize_x[Grad_dim]*ijk[Grad_dim];

           xOwned[iunk][loc_inode] = Rhobar_b_LBB[irb] + 
                      (Rhobar_b_RTF[irb]-Rhobar_b_LBB[irb])*
                           x_dist/Size_x[Grad_dim];
       }
       else {
          xOwned[iunk][loc_inode] = Rhobar_b[irb];
       }
     }
  }
  return;
}
/************************************************************/
/*setup_Xi_cavWTC: set up the cavity function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
static void setup_Xi_cavWTC(double **xOwned)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,icav;
  double vol,area,x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (icav = 0; icav < Phys2Nunk[CAVITY_WTC]; icav++){
       iunk = Phys2Unk_first[CAVITY_WTC] + icav;
       xOwned[iunk][loc_inode] = Xi_cav_b[icav+2];
     }
  }
  return;
}
/************************************************************/
/*setup_BondWTC: set up the bond function initial guesses for 
                 Wertheim-Tripathi-Chapman functionals.  For now
                use bulk segment densities initial guess.  Later calculate
                based on rho initial guess. */
static void setup_BondWTC(double **xOwned)
{
  int loc_inode,inode_box,inode,ijk[3],iunk,ibond;
  double vol,area,x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     for (ibond = 0; ibond < Nbonds; ibond++){
       iunk = Phys2Unk_first[BOND_WTC] + ibond;
       xOwned[iunk][loc_inode] = BondWTC_b[ibond];
     }
  }
  return;
}
/************************************************************/
/*setup_elec_pot: set up the electrostatic potential initial guess*/
static void setup_elec_pot(double **xOwned,int iguess)
{
  int loc_inode,iunk,inode_box,inode,ijk[3];
  double x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       iunk = Phys2Unk_first[POISSON];
       if (Lsteady_state && iguess==LINEAR){
           inode_box = L2B_node[loc_inode];
           inode     = B2G_node[inode_box];
           node_to_ijk(inode,ijk); 
           x_dist = Esize_x[Grad_dim]*ijk[Grad_dim];

           xOwned[iunk][loc_inode]  = Elec_pot_LBB + 
                                      (Elec_pot_RTF-Elec_pot_LBB)*
                                      x_dist/Size_x[Grad_dim];
       }
       else xOwned[iunk][loc_inode] = 0.0;
  }
  return;
}
/************************************************************/
/* setup_chem_pot: for cases with steady state profiles,
   set up an initial guess for (electro)chemical potentials */
static void setup_chem_pot(double **xOwned)
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
        iunk = Phys2Unk_first[DIFFUSION]+icomp;
        if (!Zero_density_TF[inode_box][icomp]){
           if (Ipot_ff_c == 1){
             xOwned[iunk][loc_inode] = log(xOwned[Phys2Unk_first[DENSITY]+icomp][loc_inode])
                           + Charge_f[icomp]*(xOwned[Phys2Unk_first[POISSON]][loc_inode]);

           }
           else{
               if (x_dist<0.) xOwned[iunk][loc_inode]=Betamu_LBB[icomp];
               else if (x_dist > x_tot) xOwned[iunk][loc_inode]=Betamu_RTF[icomp];
               else  xOwned[iunk][loc_inode] = Betamu_LBB[icomp] + (Betamu_RTF[icomp]-Betamu_LBB[icomp])* x_dist/x_tot;
           }
       }
       else xOwned[iunk][loc_inode] = -VEXT_MAX;
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
static void read_in_a_file(int iguess,char *filename)
{
  int c;
  int i,iunk,junk,idim, inode,itype_mer,ipol,iseg,index,dim_tmp,iunk_file;
  int ijk_old[3],ijk_old_max[3],open_now,ndim_max,node_start;
  int unk_in_file, unk_start_in_file[NEQ_TYPE],header,eq_type;
  int unk_to_eq_in_file[3*NCOMP_MAX+NMER_MAX+NMER_MAX*NMER_MAX+13];
  double pos_old,pos_old0[3],tmp,x_tmp;
  char filename2[20];
  char unk_char[20];
  FILE *fp5=NULL,*fp6=NULL;


                    /* open the dft_dens.dat file */
   fp5=fopen(filename,"r");

                   /* identify which unknowns are found in this file and in what order */
    header=0;
    iunk=0;
    unk_in_file=0;
    for (eq_type=0;eq_type<NEQ_TYPE;eq_type++){
       fgets(unk_char,20,fp5);
       if (strncmp(unk_char,"DENSITY",5)==0) {
             Restart_field[DENSITY]=TRUE;
             header++;
             if (Type_poly==WTC) unk_in_file+=Nseg_tot;
             else              unk_in_file+=Ncomp;
             unk_start_in_file[DENSITY]=iunk;
             if (Type_poly==WTC)
                  for (i=0;i<Nseg_tot;i++) unk_to_eq_in_file[iunk++]=DENSITY;
             else for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=DENSITY;
       }
       else if (strncmp(unk_char,"POISSON",5)==0){
             Restart_field[POISSON]=TRUE;
             header++;
             unk_in_file+=1;
             unk_start_in_file[POISSON]=iunk;
             unk_to_eq_in_file[iunk++]=POISSON;
       }
       else if (strncmp(unk_char,"RHOBAR_ROSEN",5)==0){
             Restart_field[RHOBAR_ROSEN]=TRUE;
             header++;
             unk_in_file+=Nrho_bar;
             unk_start_in_file[RHOBAR_ROSEN]=iunk;
             for (i=0;i<Nrho_bar;i++) unk_to_eq_in_file[iunk++]=RHOBAR_ROSEN;
       }
       else if (strncmp(unk_char,"CMSFIELD",5)==0){
             Restart_field[CMS_FIELD]=TRUE;
             header++;
             unk_in_file+=Ncomp;
             unk_start_in_file[CMS_FIELD]=iunk;
             for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=CMS_FIELD;
       }
       else if (strncmp(unk_char,"CAVITY_WTC",5)==0){
             Restart_field[CAVITY_WTC]=TRUE;
             header++;
             unk_in_file+=2;
             unk_start_in_file[CAVITY_WTC]=iunk;
             for (i=0;i<2;i++) unk_to_eq_in_file[iunk++]=CAVITY_WTC;
       }
       else if (strncmp(unk_char,"BOND_WTC",5)==0){
             Restart_field[BOND_WTC]=TRUE;
             header++;
             unk_in_file+=Nbonds;
             unk_start_in_file[BOND_WTC]=iunk;
             for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=BOND_WTC;
       }
       else if (strncmp(unk_char,"CHEMPOT",5)==0){
             Restart_field[DIFFUSION]=TRUE;
             header++;
             if (Type_poly==WTC) unk_in_file+=Nseg_tot;
             else              unk_in_file+=Ncomp;
             unk_start_in_file[DIFFUSION]=iunk;
             if (Type_poly==WTC)
                 for (i=0;i<Nseg_tot;i++) unk_to_eq_in_file[iunk++]=DIFFUSION;
             else
                 for (i=0;i<Ncomp;i++) unk_to_eq_in_file[iunk++]=DIFFUSION;
       }
       else eq_type=NEQ_TYPE;  /* end loop */
    }
    if (Iwrite != NO_SCREEN) printf("Number of unknowns in the file=%d\n",unk_in_file);

    if (Lsteady_state && Restart_field[DIFFUSION]==FALSE) 
         if (Proc==0 && Iwrite != NO_SCREEN)
           printf("there is no chemical potential data in the restart file\n");
    if (Type_coul != NONE && Restart_field[POISSON]==FALSE)
           printf("there is no electrostatic potential data in the restart file\n");
    if (Type_poly != NONE && Type_poly!=WTC && Restart_field[CMS_FIELD]==FALSE)
           printf("there is no CMS field data in the restart file\n");
    if ((Type_poly==NONE || Type_poly==WTC)  && Restart_field[RHOBAR_ROSEN]==FALSE)
           printf("there is no Rosenfeld nonlocal density data in the restart file\n");
    if (Restart_field[DENSITY]==FALSE)
           printf("there is no density data in the restart file\n");

    fclose(fp5);
    if (Restart != 5) Nodes_old-=header;
    printf("skipping %d lines in the dft_dens.dat file\n",header);

  /* read positions from file find Nodes_x_old[idim] */
  /* read the densities and electrostatic potentials from file */

  open_now=TRUE;
  for (idim=0;idim<Ndim;idim++) ijk_old_max[idim] = 0;
  for (index=0; index<Nodes_old; index++) {

    if (open_now){
       fp5=fopen(filename,"r");
       if (Type_poly != NONE && Type_poly!=WTC){
         sprintf(filename2,"%sg",filename);
         fp6=fopen(filename2,"r");
       }
       open_now=FALSE;
                                      /* discard header when ready to read */
       for (i=0;i<header;i++) while ((c=getc(fp5)) != EOF && c !='\n') ; 
    }

    if (Restart==5) ndim_max=1;  /* again for using 1D solution for 2D/3D guess */
    else ndim_max=Ndim;
    /*printf("Proc=%d  ndim_max=%d\n",Proc,ndim_max);*/

                                 /* find number of nodes in each dimension in the file */
    for (idim=0;idim<ndim_max;idim++)         {
      fscanf(fp5,"%lf",&pos_old);
      if (index==0) pos_old0[idim]=pos_old;
      pos_old -= pos_old0[idim];
                                  /* this code transforms x to y etc in the guess */
                                      /*if (idim==0) dim_tmp=1;
                                        else if (idim==1) dim_tmp=0;
                                        else   dim_tmp=idim; */

      dim_tmp=idim;
      ijk_old[dim_tmp] = round_to_int(pos_old/Esize_x[dim_tmp]);
      if (Type_poly!=NONE && Type_poly!=WTC)  fscanf(fp6,"%lf",&tmp); /* ignore positions in densg files. */

      if (ijk_old[dim_tmp] > ijk_old_max[dim_tmp]) ijk_old_max[dim_tmp] = ijk_old[dim_tmp];
    }
                                    /* identify the node and starting unknown number at that node */
    inode=ijk_to_node(ijk_old);
    if (Restart==5) inode=index;
    node_start=inode*Nunk_per_node;

                                   /* loop over unknows assume the order of input is correct */

    for (iunk_file=0;iunk_file<unk_in_file;iunk_file++) {
          eq_type = unk_to_eq_in_file[iunk_file];

          switch (eq_type){
              case CMS_FIELD: 
   	         fscanf(fp5,"%lf",&tmp);
                 tmp = exp(-tmp); 
                 break;

              case DENSITY:
                               /* icomp = iunk_file-unk_start_in_file[DENSITY]; */
                               /* if (Lsteady_state) tmp *= Mass[icomp]; */
   	         fscanf(fp5,"%lf",&tmp); 
                 break;

              case RHOBAR_ROSEN:
              case POISSON:
              case CAVITY_WTC:
              case BOND_WTC:
   	         fscanf(fp5,"%lf",&tmp); 
                 break;

              case DIFFUSION: 
                               /* icomp = iunk_file-unk_start_in_file[DIFFUSION]; */
                               /* if (Ipot_ff_n != IDEAL_GAS) */  /* Debroglie wavelength contribution */
                               /* tmp -= 3.0*log(Sigma_ff[icomp][icomp]+1.5*log(Mass[icomp]*Temp));*/
   	         fscanf(fp5,"%lf",&tmp); break;
           }
       iunk = Phys2Unk_first[eq_type]+iunk_file-unk_start_in_file[eq_type];
       if (Lbinodal && iguess==BINODAL_FLAG) X2_old[iunk+node_start]=tmp;
       else                                  X_old[iunk+node_start]=tmp;
    }
 
    if (Type_poly != NONE && Type_poly!=WTC){
        for (iunk=Phys2Unk_first[CMS_G];iunk<Phys2Unk_last[CMS_G];iunk++) {

         fscanf(fp6,"%lf",&tmp);
         ipol=Unk_to_Poly[iunk-Phys2Unk_first[CMS_G]];
         iseg=Unk_to_Seg[iunk-Phys2Unk_first[CMS_G]];
         itype_mer=Phys2Unk_first[CMS_FIELD]+Type_mer[ipol][iseg];

         if (Lbinodal && iguess==BINODAL_FLAG) X2_old[iunk+node_start]=tmp;
         else                                  X_old[iunk+node_start]=tmp;

       } /* end of loop over unknowns */
    }  /* end of loading for the CMS G equations. */


             /* read extra variables on the line and ignore them -
                for example the Poisson-Boltzman electrolyte output */
    while ((c=getc(fp5)) != EOF && c !='\n') ;
    if (Restart==5 && ijk_old[0]==Nodes_x[0]-1){
       fclose(fp5);
       if (Type_poly != NONE && Type_poly !=WTC) fclose(fp6);
       open_now=TRUE;
/* if (Proc==0) printf("closing files to read again!\n");*/
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
static void shift_the_profile(double *x_new,double fac)
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
     if (Unk2Phys[iunk]==DENSITY){
        x_test = min(Rho_max,fac*(unk_old-Rho_b[iunk-Phys2Unk_first[DENSITY]])+ Rho_b[iunk-Phys2Unk_first[DENSITY]]);
     }
     else if (Unk2Phys[iunk]==RHOBAR_ROSEN){
        x_test = fac*(unk_old-Rhobar_b[iunk-Phys2Unk_first[RHOBAR_ROSEN]])+ Rhobar_b[iunk-Phys2Unk_first[RHOBAR_ROSEN]];
        if (iunk == Ncomp && x_test >= 1.0) x_test = Rhobar_b[iunk-Phys2Unk_first[RHOBAR_ROSEN]];
     }
     else if (Unk2Phys[iunk]==POISSON){
        x_test = fac*(unk_old-1.0) + 1.0;
     }
     else 
        x_test = unk_old;

     x_new[inode*Nunk_per_node+iunk] = x_test;

     }
  }
  return;
}
/*********************************************************************/
/* locate_inode_old: find a given node in the old solution array */
static int locate_inode_old(int *ijk)
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
static void communicate_profile(double *x_new, double** xOwned)
{
   int loc_inode,inode,iunk;   
   
    MPI_Bcast (x_new, Nnodes*Nunk_per_node,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       inode = L2G_node[loc_inode];
       for (iunk=0; iunk<Nunk_per_node; iunk++){
           xOwned[iunk][loc_inode] = x_new[inode*Nunk_per_node+iunk];
       }
    }
    return;
}
/*********************************************************************/
/*check_zero_densities: here just remove zero densities where 
         not appropriate, and make sure x=0 where needed */
void check_zero_densities(double **xOwned)
{

  int loc_inode,icomp,inode_box,iunk,iloop,nloop;

  if (Type_poly==WTC) nloop=Nseg_tot;
  else nloop=Ncomp;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
      inode_box = L2B_node[loc_inode];
      for (iloop=0; iloop<nloop; iloop++){
         icomp=Unk2Comp[iloop];
	 iunk = Phys2Unk_first[DENSITY]+iloop;
         if (Zero_density_TF[inode_box][icomp])
                 xOwned[iunk][loc_inode] = 0.0;
         else{
           if (Type_poly==WTC)
              if (xOwned[iunk][loc_inode] < Rho_seg_b[iunk]*exp(-VEXT_MAX)) {
                  xOwned[iunk][loc_inode] = Rho_seg_b[iunk]*exp(-VEXT_MAX); /*DENSITY_MIN*/
              }
           else
              if (xOwned[iunk][loc_inode] < Rho_b[icomp]*exp(-VEXT_MAX)) {
                  xOwned[iunk][loc_inode] = Rho_b[icomp]*exp(-VEXT_MAX); /*DENSITY_MIN*/
              }
         }
      }
  }
  return;
}
/**********************************************************************/
/*chop_profile: do the profile chop here. */
static void chop_profile(double **xOwned, int iguess)
{
  int loc_inode,icomp,iwall,check,iunk;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++)
    for (icomp=0; icomp<Ncomp; icomp++) {
        check = 0;
        for (iwall=0; iwall<Nwall; iwall++)
              if (X_wall[loc_inode][iwall] > Xstart_step[0]) check++;
        if (check == Nwall){
	      iunk = Phys2Unk_first[DENSITY]+icomp;
              if (iguess==CHOP_RHO_L) xOwned[iunk][loc_inode] = Rho_coex[1];
              else if (iguess==CHOP_RHO_V) xOwned[iunk][loc_inode] = Rho_coex[0];
              else if (iguess==CHOP_RHO) xOwned[iunk][loc_inode] = Rho_b[icomp];
              else if (iguess==CHOP_RHO_STEP) xOwned[iunk][loc_inode] = Rho_step[icomp][0];
        }
    }
  return;
}
/****************************************************************************/
