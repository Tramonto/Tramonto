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

/* New version of dft_out_ads.c ---- testing a function pointer approach for the output */

#include "dft_out_ads.h"

/**************************************************************************************/
void calc_adsorption(FILE *fp,double **x)
{
  int icomp,iunk,nloop,idim,jcomp;
  double ads[NCOMP_MAX],ads_ex[NCOMP_MAX],ads_b[NCOMP_MAX],volume,total_ads;
  static int first=TRUE;


  Integration_profile=NULL;
  if (Lseg_densities) nloop=Nseg_tot;
  nloop=Ncomp;

  if(!first) {

  if (Proc==0 &&Iwrite != NO_SCREEN) printf("-------------------- ADSORPTION -------------------------------\n");
  for (icomp=0;icomp<nloop;icomp++) ads[icomp]=0.0;
  for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
     if (Lseg_densities) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else                icomp = iunk-Phys2Unk_first[DENSITY];
     ads[icomp]+=integrateInSpace(&integrand_adsorption,iunk,Nel_hit2,x,Integration_profile);
  }
  }

 if (Proc==0 && Iwrite != NO_SCREEN){
     total_ads=0.0;
     for (icomp=0;icomp<nloop;icomp++) total_ads+=ads[icomp];
     volume=1.0;

     for (idim=0;idim<Ndim;idim++) volume*=Size_x[idim];
     for (icomp=0;icomp<nloop;icomp++){
          if(!first){
               if (LBulk) print_to_screen_comp(icomp,ads[icomp]/volume,"DENSITY");
               else       print_to_screen_comp(icomp,ads[icomp],"ADSORPTION");
          }
          if (fp !=NULL){
               if (LBulk){ 
                   print_to_file_comp(fp,icomp,ads[icomp]/volume,"rho",first);
                   if (Ncomp >1){
                       print_to_file_comp(fp,icomp,ads[icomp]/total_ads,"num_frac",first);
                   }
               }
               else       print_to_file_comp(fp,icomp,ads[icomp],"ads",first);
          }
      }
  }    

 if(!first) {
  for (icomp=0;icomp<nloop;icomp++) ads_b[icomp]=0.0;
  for (iunk=Phys2Unk_first[DENSITY];iunk<Phys2Unk_last[DENSITY];iunk++) {
     if (Lseg_densities) icomp=Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else                icomp = iunk-Phys2Unk_first[DENSITY];
     ads_b[icomp]+=integrateInSpace(&integrand_adsorption_bulk,iunk,Nel_hit,x,Integration_profile);
     ads_ex[icomp]=ads[icomp]-ads_b[icomp];
  }
 }

  if (Proc==0 && Iwrite != NO_SCREEN){
     for (icomp=0;icomp<nloop;icomp++){
        if(!first && !LBulk && Type_interface==UNIFORM_INTERFACE) print_to_screen_comp(icomp,ads_ex[icomp],"EXCESS ADSORPTION");
        if (fp !=NULL && !LBulk && Type_interface==UNIFORM_INTERFACE) print_to_file_comp(fp,icomp,ads_ex[icomp],"ads_ex",first);
     }    
  }
  if (first) first=FALSE;
  if (Proc==0 &&Iwrite != NO_SCREEN) printf("---------------------------------------------------------------\n");
  return;
}
/**************************************************************************************/
void calc_fluid_charge(FILE *fp,double **x)
{
 static int first=TRUE;
 double charge;

 if(!first) {
   if (Proc==0&&Iwrite != NO_SCREEN) printf("-------------------- CHARGE     -------------------------------\n");
   Integration_profile=NULL;

   charge=integrateInSpace_SumInComp(&integrand_fluid_charge,Nel_hit2,x,Integration_profile);
 }

 if (Proc==0 && Iwrite != NO_SCREEN){
      if(!first) print_to_screen(charge,"CHARGE IN FLUID");
      if (fp !=NULL) print_to_file(fp,charge,"charge",first);
      if(!first) printf("---------------------------------------------------------------\n");
 }
 if (first) first=FALSE;
 return;
}
/**************************************************************************************/
double integrand_adsorption(int iunk,int inode_box, double **x)
{
     double integrand;
     integrand = x[iunk][inode_box];
     return(integrand);
}
/**************************************************************************************/
double integrand_adsorption_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;

     if (Lseg_densities) rho_bulk = Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
     else                rho_bulk = Rho_b[iunk-Phys2Unk_first[DENSITY]];

     integrand = rho_bulk;
     return(integrand);
}
/**************************************************************************************/
double integrand_fluid_charge(int iunk, int inode_box, double **x)
{
     double integrand;
     int icomp;

     if (Lseg_densities) icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     else                icomp = iunk-Phys2Unk_first[DENSITY];

     integrand = Charge_f[icomp]*x[iunk][inode_box];
     return(integrand);
}
/**************************************************************************************/

