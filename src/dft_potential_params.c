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
 *  FILE: dft_potentail_params.c
 *
 *  This file contains basic interaction potential parameters.
 */

#include "dft_potential_params.h"
 
/****************************************************************************/
/* pot_parameters: calculate the cross terms (sigmaij,epsilonij,cutoffij) 
   for this potential */
void pot_parameters(char *output_file1)
{        
 int i,j,iw,jw, printproc=FALSE;
 double cut=0.0;
 FILE *fp2=NULL;
 if (Proc==0 && output_file1 != NULL) printproc = TRUE;
 else printproc=FALSE;
 if (printproc) {
   if( (fp2 = fopen(output_file1,"a+")) == NULL) {
     printf("Can't open file %s\n", output_file1);
     exit(1);
   }  
 }       
 if (printproc) fprintf(fp2,"****** SUMMARY OF POTENTIAL PARAMETERS ************\n");
               
 for (i=0; i<Ncomp; i++){
     if (Ipot_ff_n==IDEAL_GAS) Sigma_ff[i][i]=0.0;
     if (Type_poly ==NONE) Bond_ff[i][i]=0.0;
               
     for (j=0; j<Ncomp; j++){
      
        if (i != j){ 
           Sigma_ff[i][j] = 0.5*(Sigma_ff[i][i]+Sigma_ff[j][j]);
           Eps_ff[i][j] = sqrt(Eps_ff[i][i]*Eps_ff[j][j]); 
           Cut_ff[i][j] = 0.5*(Cut_ff[i][i]+Cut_ff[j][j]);
           Bond_ff[i][j] = 0.5*(Bond_ff[i][i]+Bond_ff[j][j]);
           if (Type_pairPot==PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS) 
			   YukawaK_ff[i][j]=0.5*(YukawaK_ff[i][i]+YukawaK_ff[j][j]);
        }
     
        if (printproc) {
           fprintf(fp2,"\ti:%d  j:%d",i,j);
           fprintf(fp2,"   Sigma_ff: %9.6f  Cut_ff: %9.6f  Eps_ff: %9.6f  Bond_ff: %9.6f \n",
                   Sigma_ff[i][j],Cut_ff[i][j],Eps_ff[i][j],Bond_ff[i][j]);
           if (Type_pairPot==PAIR_YUKAWA_CS || Type_pairPot == PAIR_EXP_CS){
               fprintf(fp2,"\ti:%d  j:%d",i,j);
               fprintf(fp2,"   YukawaK_ff: %9.6f\n",YukawaK_ff[i][j]);
           }
        }
     }
 } 
       
 for (iw=0; iw<Nwall_type; iw++){
   for (jw=0; jw<Nwall_type; jw++){
        Sigma_ww[iw][jw] = 0.5*(Sigma_w[iw]+Sigma_w[jw]);
        Eps_ww[iw][jw] = sqrt(Eps_w[iw]*Eps_w[jw]);
        if (Type_vext3D==PAIR_YUKAWA_CS || Type_uwwPot==PAIR_YUKAWA_CS) YukawaK_ww[iw][jw]=0.5*(YukawaK_w[iw]+YukawaK_w[jw]);
	    if (Type_vext3D==PAIR_EXP_CS || Type_uwwPot==PAIR_EXP_CS) YukawaK_ww[iw][jw]=0.5*(YukawaK_w[iw]+YukawaK_w[jw]);
        if (iw != jw){
           Cut_ww[iw][jw] = 0.5*(Cut_ww[iw][iw]+Cut_ww[jw][jw]);
        }
        if (printproc){
          fprintf(fp2,"\tiwall_type:%d  jwall_type:%d",iw,jw);
          fprintf(fp2,"   Sigma_ww: %9.6f  Cut_ww: %9.6f  Eps_ww: %9.6f\n",
                           Sigma_ww[iw][jw],Cut_ww[iw][jw],Eps_ww[iw][jw]);
          if (Type_vext3D==PAIR_YUKAWA_CS || Type_uwwPot==PAIR_YUKAWA_CS || Type_vext3D==PAIR_EXP_CS || Type_uwwPot==PAIR_EXP_CS){
              fprintf(fp2,"\tiwall_type:%d  jwall_type:%d",iw,jw);
              fprintf(fp2,"   YukawaK_ww: %9.6f\n",YukawaK_ww[iw][jw]);
          }
        }
   }
 }

 for (i=0; i<Ncomp; i++){

     for (iw=0; iw<Nwall_type; iw++){
        Sigma_wf[i][iw] = 0.5*(Sigma_ff[i][i] + Sigma_w[iw]);
        Eps_wf[i][iw] = sqrt(Eps_ff[i][i]*Eps_w[iw]);
        Cut_wf[i][iw] = 0.5*(Cut_ff[i][i]+Cut_ww[iw][iw]);
        if (Type_vext3D==PAIR_YUKAWA_CS || Type_vext3D==PAIR_YUKAWA_CS) YukawaK_wf[i][iw] = 0.5*(YukawaK_ff[i][i]+YukawaK_ww[iw][iw]);

        if (printproc) {
          fprintf(fp2,"\ti:%d  iwall_type:%d",i,iw);
          fprintf(fp2,"   Sigma_wf: %9.6f  Cut_wf: %9.6f  Eps_wf: %9.6f\n",
           Sigma_wf[i][iw],Cut_wf[i][iw],Eps_wf[i][iw]);
          if (Type_vext3D==PAIR_YUKAWA_CS || Type_vext3D==PAIR_YUKAWA_CS) {
             fprintf(fp2,"\ti:%d  iwall_type:%d",i,iw);
             fprintf(fp2,"   YukawaK_wf: %9.6f \n",YukawaK_wf[i][iw]);
          }
        }
     }
 }

 if (printproc) {
     fprintf(fp2,"***************************************************\n");
     fclose(fp2);
 }
}
/***********************************************************************/
