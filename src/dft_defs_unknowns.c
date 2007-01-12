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
 *   dft_defs_unknowns.c: set up basic unknown types for the problem of interest.
 */
#include "dft_defs_unknowns.h"
/*****************************************************************************/
/*setup_nunk_per_node:  here we just set up the basic parameters for the number
  of unknowns per node for the run of interest, and some arrays to move between
  unknown number and equation type easily. */
void setup_nunk_per_node(char *output_file1)
{
  int i,iunk,icomp,unk_rel;
  int NCMSField_unk;
  FILE *fp2=NULL;

  if (Proc==0) {
    if( (fp2 = fopen(output_file1,"a+"))==NULL){
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
  }	

  for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY:                  /* unknowns of Euler-Lagrange equation */
            if (Type_poly==WTC) Phys2Nunk[DENSITY]=Nseg_tot;
            else                Phys2Nunk[DENSITY]=Ncomp;
            break;

         case HSRHOBAR:     /* unknowns of Nonlocal Density Eqns for Rosenfeld Functionals */
            Nrho_bar = 0;
            if (Ipot_ff_n != IDEAL_GAS &&( Type_poly == NONE || Type_poly ==WTC)){
                 Nrho_bar = 4 + 2*Ndim;
                 Nrho_bar_s = 4;
            }
            Phys2Nunk[HSRHOBAR]=Nrho_bar;
            break; 

         case POISSON:                             /* unknowns of Poisson Equation */
            Npoisson=0;
            if ( Type_coul !=NONE){
                 Npoisson=1;
            }
            Phys2Nunk[POISSON]=Npoisson;
            break;

         case DIFFUSION:                            /* unknowns of Diffusion Equation */
            Ndiffusion=0;
            if ( Lsteady_state && (Type_poly==NONE || Type_poly==WTC)){
              if (Type_poly==WTC) Ndiffusion=Nseg_tot; 
              else                Ndiffusion=Ncomp;
            }
            Phys2Nunk[DIFFUSION]=Ndiffusion;
            break;

         case CAVWTC:
            Nrho_bar_cavity=0;
            if (Type_poly==WTC){
                 Nrho_bar_cavity = 4;
                 Phys2Nunk[CAVWTC]=Nrho_bar_cavity-2;  
                                   /* strange case because y function only uses two of the 
                                      defined xi's.  But, I'm leaving a placeholder for 4 */
            }
            else Phys2Nunk[CAVWTC]=0;
            break;

         case BONDWTC:
            Nrho_bar_bond=0;
            if (Type_poly==WTC) Nrho_bar_bond = Nbonds;
            Phys2Nunk[BONDWTC]=Nrho_bar_bond;
            break;

         case CMS_FIELD:
            NCMSField_unk = 0;
            if (Type_poly == CMS || Type_poly==CMS_SCFT){
                 NCMSField_unk=Ncomp;
            }
            Phys2Nunk[CMS_FIELD]=NCMSField_unk;
            break;

         case CMS_G:
            if (Type_poly == CMS || Type_poly==CMS_SCFT){
                 Phys2Nunk[CMS_G] = Ngeqn_tot;
            }
            break;

         default:
            printf("problems with defining equation type %d\n",i);
            exit(-1);
            break;
     }
     if (i==0) Phys2Unk_first[i]=0;
     else      Phys2Unk_first[i]=Phys2Unk_last[i-1];
     Phys2Unk_last[i]=Phys2Unk_first[i]+(Phys2Nunk[i]);
     for (iunk=Phys2Unk_first[i]; iunk< Phys2Unk_last[i]; iunk++) Unk2Phys[iunk]=i;
     Nunk_per_node += Phys2Nunk[i];
  }
  for (i=0;i<NEQ_TYPE;i++){
     if (Phys2Nunk[i]==0){
        Phys2Unk_first[i]=NO_UNK;
        Phys2Unk_last[i]=NO_UNK;
     }
     if (i==CMS_G && (Type_poly == CMS || Type_poly==CMS_SCFT)){
        for (icomp=0;icomp<Npol_comp;icomp++) Geqn_start[icomp]+=Phys2Unk_first[i];
     }
  }

/*   for (iunk=0;i<Nunk_per_node;iunk++){
     if (Unk2Phys[iunk]==DENSITY || Unk2Phys[iunk]==DIFFUSION){
         if(Type_poly==WTC){ 
            unk_rel=iunk-Phys2Unk_first[Unk2Phys[iunk]];
            icomp = Unk2Comp[unk_rel];
         }
         else icomp=iunk-Phys2Unk_first[Unk2Phys[iunk]];
         Unk2Comp[iunk]=icomp;
     }
   }*/

   if (Proc==0){
   if (Iwrite==VERBOSE){
        printf("\n******************************************************\n");
        printf("TOTAL Nunk_per_node=%d\n",Nunk_per_node);
        for (i=0;i<NEQ_TYPE;i++) printf("Phys2Nunk[%d]=%d  start_unk=%d  end_unk=%d\n",
                                   i,Phys2Nunk[i],Phys2Unk_first[i],Phys2Unk_last[i]);
        for (iunk=0;iunk<Nunk_per_node;iunk++) printf("iunk=%d equation_type=%d\n",iunk,Unk2Phys[iunk]);
        printf("******************************************************\n");
   }
   fprintf(fp2,"\n******************************************************\n");
   fprintf(fp2,"TOTAL Nunk_per_node=%d\n",Nunk_per_node);
   for (i=0;i<NEQ_TYPE;i++) fprintf(fp2,"Phys2Nunk[%d]=%d  start_unk=%d  end_unk=%d\n",
                                   i,Phys2Nunk[i],Phys2Unk_first[i],Phys2Unk_last[i]);
   for (iunk=0;iunk<Nunk_per_node;iunk++) fprintf(fp2,"iunk=%d equation_type=%d\n",iunk,Unk2Phys[iunk]);
   fprintf(fp2,"******************************************************\n");
   }
   if (Proc==0)fclose(fp2);
   return;
}
/*******************************************************************************/
