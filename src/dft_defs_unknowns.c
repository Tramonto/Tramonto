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
  int NCMSField_unk, NYW_Dens, NWJDCField_unk;
  FILE *fp2=NULL;

  if (Proc==0) {
    if( (fp2 = fopen(output_file1,"a+"))==NULL){
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
  }	


  /* set a couple generic logicals that can be used to toggle between 
     component and segment densities (Lseg_densities), and Hard-Sphere
     perturbation DFTs and other types of DFTs (L_HSperturbation) */

  if (Type_poly==WTC || Type_poly==WJDC) Lseg_densities=TRUE;
  else Lseg_densities=FALSE;

  if (Type_poly == CMS || Type_poly==CMS_SCFT) L_HSperturbation=FALSE;
  else L_HSperturbation=TRUE;

  for (i=0;i<NEQ_TYPE;i++){
     switch(i){
         case DENSITY:                  /* unknowns of Euler-Lagrange equation */
            if (Lseg_densities){ Phys2Nunk[DENSITY]=Nseg_tot; }
            else{                Phys2Nunk[DENSITY]=Ncomp; }
            break;

         case HSRHOBAR:     /* unknowns of Nonlocal Density Eqns for Rosenfeld Functionals */
            Nrho_bar = 0;
            if (Ipot_ff_n != IDEAL_GAS && L_HSperturbation){
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
            if ( Lsteady_state && L_HSperturbation){
              if (Lseg_densities) Ndiffusion=Nseg_tot; 
              else                Ndiffusion=Ncomp;
            }
            Phys2Nunk[DIFFUSION]=Ndiffusion;
            break;

         case CAVWTC:
            Nrho_bar_cavity=0;
            if (Type_poly==WTC || Type_poly==WJDC){
                 Nrho_bar_cavity = 4;
                 Phys2Nunk[CAVWTC]=Nrho_bar_cavity-2;  
                                   /* strange case because y function only uses two of the 
                                      defined xi's.  But, I'm leaving a placeholder for 4 */
            }
            else Phys2Nunk[CAVWTC]=0;
            break;

         case BONDWTC:
            if (Type_poly==WTC){
                Nrho_bar_bond = Nbonds;
            }
            else Nrho_bar_bond=0;
            Phys2Nunk[BONDWTC]=Nrho_bar_bond;
            break;

         case CMS_FIELD:
            if (Type_poly == CMS || Type_poly==CMS_SCFT){
                 NCMSField_unk=Ncomp;
            }
            else NCMSField_unk = 0;
            Phys2Nunk[CMS_FIELD]=NCMSField_unk;
            break;

         case WJDC_FIELD:
            if (Type_poly == WJDC){
                 NWJDCField_unk=Ncomp;
            }
            else NWJDCField_unk = 0;
            Phys2Nunk[WJDC_FIELD]=NWJDCField_unk;
            break;

         case G_CHAIN:
            if (Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC){
                 Phys2Nunk[G_CHAIN] = Ngeqn_tot;
            }
            else Phys2Nunk[G_CHAIN]=0;
            break;

          case YW_DENS:
	    NYW_Dens = 0;
	    Phys2Nunk[YW_DENS] = 0;
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
     if ( (i==G_CHAIN && (Type_poly == CMS || Type_poly==CMS_SCFT || Type_poly==WJDC))){
        for (icomp=0;icomp<Npol_comp;icomp++) Geqn_start[icomp]+=Phys2Unk_first[i];
     }
  }

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
