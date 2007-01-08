/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_guess_noinfo_CMS.c
 *
 *  This file contains routines that set up an initial guess for the
 *  CMS polymer DFT.
 *
 */

/*#include "dft_globals_const.h"
#include "mpi.h"
#include <string.h>*/
#include "dft_guess_CMS.h"
 
/*********************************************************/
/*setup_polymer_field: in this routine sets up the initial guess for the CMS field variable */
void setup_polymer_field(double **xOwned, int iguess)
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
void setup_polymer_simple(double **xOwned, int iguess)
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
           if ((Type_poly != NONE && Type_poly != WTC) || (temp > 0.5)) xOwned[iunk][loc_inode] = exp(temp);
           else                             xOwned[iunk][loc_inode] = exp(1. - sqrt(1.-2.*temp));
        } /* end if i not zero density  */
        else   xOwned[iunk][loc_inode] = DENSITY_MIN; /* zero density - Boltzmann probability = 0 */
     } /* end loop over itype_mer */
  } /* end loop over loc_inode */
  return;
}
/*********************************************************/
/*setup_polymer_rho: in this routine set up polymer density profiles    */
void setup_polymer_rho(double **xOwned, int iguess)
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
void setup_polymer_G(double **xOwned)
{
  int loc_inode,inode_box,ijk_box[3],loc_i;
  int reflect_flag[NDIM_MAX];
  int   **sten_offset, *offset, isten;
  double *sten_weight,  weight;
  struct Stencil_Struct *sten;
  int sten_type,izone,jlist,jnode_box,jtype_mer,itype_mer;
  int iunk,poln,iseg,ibond,not_done,junk,cycle,loc_B;

     sten_type = DELTA_FN;
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
                              /*xOwned[iunk][loc_inode] +=  weight*xOwned[junk][loc_inode]; */
                              xOwned[iunk][loc_inode] +=  weight; 
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
