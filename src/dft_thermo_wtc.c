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

/* ----------------------------------------------------------------------
dft_thermo_wtc.c:
Calculate the thermodynamic properties of chain contributions for a 
Wertheim-Tripathi-Chapman bonded fluid.
-------------------------------------------------------------------------*/
#include "dft_thermo_wtc.h"

/****************************************************************************/
/* WTC_thermo_precalc: call all routines needed to process bulk properties of 
                          WTC functionals */
void WTC_thermo_precalc(char *output_file1)
{
  int iseg;

  WTC_overlap();

  /* compute bulk segment densities from polymer component densities */
  for (iseg=0;iseg<Nseg_tot;iseg++){
     if (Type_interface!=UNIFORM_INTERFACE){
        Rho_seg_LBB[iseg]=Rho_b_LBB[Unk2Comp[iseg]]/(double)Nmer_t_total[Unk2Comp[iseg]];
        Rho_seg_RTF[iseg]=Rho_b_RTF[Unk2Comp[iseg]]/(double)Nmer_t_total[Unk2Comp[iseg]];
     }
     else Rho_seg_b[iseg]=Rho_b[Unk2Comp[iseg]]/(double)Nmer_t_total[Unk2Comp[iseg]];
  }
  compute_bulk_nonlocal_wtc_properties(output_file1); 

  return;
}
/****************************************************************************/
/* pressure_WTC: this routine calculates the pressure contribution 
   for the WTC functional */
double pressure_WTC(double *rho_seg,double *xi_cav)
{
  int iseg,icomp,count_bond,ibond;
  double betap_wtc,term_iseg;

  betap_wtc = 0.0;

  /* first include the term for the number of bonds per segment */

  for (iseg=0;iseg<Nseg_tot;iseg++){
     count_bond=0;
     for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
        if (Bonds_SegAll[iseg][ibond]!=-1 && Bonds_SegAll[iseg][ibond]!=-2) count_bond++;
     }  
     betap_wtc -= 0.5*rho_seg[iseg]*count_bond;
  }

  /* now calculate the chain term for the pressure */
  for (iseg=0;iseg<Nseg_tot;iseg++){
     icomp = Unk2Comp[iseg];
     term_iseg=chain_term(iseg,icomp,rho_seg,xi_cav);
     betap_wtc -= term_iseg*rho_seg[iseg];
  }

  return(betap_wtc);
}
/****************************************************************************/
/* chempot_WTC: compute the excess chemical potential contributions due to WTC
   functionals */
/* note that: Betamu_seg = total chemical potential for each segment, including ideal, HS, and att contributions
   Betamu_wtc = contribution to chem. potential from the chain part of the functional only */
void chempot_WTC(double *rho_seg,double *betamu, double *xi_cav)
{
   int icomp,jcomp,kcomp,i,iseg,ibond,jseg,kseg,pol_num,count_comp;
   double y,term_kseg;

   for (iseg=0;iseg<Nseg_tot;iseg++) {
        Betamu_wtc[iseg]=0.0;
        Betamu_seg[iseg]=0.0;
   }
    for (icomp=0; icomp<Ncomp;icomp++){
        for (pol_num=0; pol_num<Npol_comp;pol_num++) Scale_fac_WJDC[pol_num][icomp]=0.0;
    }


      /* first do ideal gas correction for segment densities rather than component densities */ 
   for (iseg=0; iseg<Nseg_tot;iseg++){
      icomp=Unk2Comp[iseg];
      Betamu_seg[iseg]=betamu[icomp]-log(Nmer_t_total[icomp]); 
   }

   /* first compute terms that are based on bond pairs starting from segment iseg*/
   /* note that in the bulk the BondWTC is identical to the site density of the jth segment */
   for (iseg=0; iseg<Nseg_tot;iseg++){
      icomp=Unk2Comp[iseg];
      for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
          jseg=Bonds_SegAll[iseg][ibond];
          if (jseg >=0){
          jcomp=Unk2Comp[jseg];
          y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_cav[2],xi_cav[3]);
          Betamu_seg[iseg] += 0.5*(1.0-Fac_overlap[icomp][jcomp]*log(y)-log(rho_seg[jseg])
                              - rho_seg[jseg]/rho_seg[iseg]);
          Betamu_wtc[iseg] += 0.5*(1.0-Fac_overlap[icomp][jcomp]*log(y)-log(rho_seg[jseg])
                              -rho_seg[jseg]/rho_seg[iseg]);
         }
      }
   }
   /* now add in the term that involves _all_ the bond pairs in the system for each segment via that cavity
      correlation function.  Note that this term is nonzero even for bond pairs on different polymer
      components than the one where the iseg segment is found ! */
   for (kseg=0; kseg<Nseg_tot;kseg++){
     kcomp = Unk2Comp[kseg];
     term_kseg=chain_term(kseg,kcomp,rho_seg,xi_cav);
     Betamu_seg[kseg] -= term_kseg;
     Betamu_wtc[kseg] -= term_kseg;
   }

   if (Physics_scaling &&(Type_poly==WJDC || Type_poly==WJDC2 || Type_poly==WJDC3)){
      if (Proc==0 && Iwrite != NO_SCREEN) printf("Physics scaling is turned on, and the value of the parameter is set to:\n");
      for (icomp=0; icomp<Ncomp;icomp++){ 
         count_comp=0;
         for (iseg=0; iseg<Nseg_tot;iseg++){
            if (Unk2Comp[iseg]==icomp){  
               count_comp++;
               pol_num=SegAll_to_Poly[iseg];
               Scale_fac_WJDC[pol_num][icomp]+=Betamu_seg[iseg];
            }
         }
         if (count_comp>0) Scale_fac_WJDC[pol_num][icomp]/=(double)count_comp;
         /*Scale_fac_WJDC[2][2]=-4.5;*/
         /*Scale_fac_WJDC[2][3]=-4.5;*/
         if (Proc==0 && Iwrite != NO_SCREEN) printf("pol_num=%d icomp=%d Scale_fac_WJDC[pol_num][icomp]=%g\n",icomp,pol_num,Scale_fac_WJDC[pol_num][icomp]);
      }
      if (Proc==0 && Iwrite != NO_SCREEN) printf("NOTE THATE THE SCALING FACTOR IS A HURISTIC THAT MAY NOT BE OPTIMAL IN SOME CASES\n");
      if (Proc==0 && Iwrite != NO_SCREEN) printf("THIS HEURISTIC CAN BE MODIFIED IN dft_thermo_wtc.c\n");
   }

   return;
}
/****************************************************************************/
/* chain_term:  This routine compute the "chain" part of the chemical potential
   that involves _all_ the bond pairs in the system for each segment.  This computation
   is done both for chemical potential and pressure so we separate this part here.*/
double chain_term(int kseg,int kcomp,double *rho_seg,double *xi_cav)
{
  double sig2,sig3,y,dydxi2,dydxi3,term_kseg=0.0;
  int iseg,icomp,ibond,jseg,jcomp;

  sig2=Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp];
  sig3=Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp]*Sigma_ff[kcomp][kcomp];
  for (iseg=0; iseg<Nseg_tot;iseg++){
        icomp=Unk2Comp[iseg];
        for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
          if(Bonds_SegAll[iseg][ibond] != -1 && Bonds_SegAll[iseg][ibond] != -2){
             jseg=Bonds_SegAll[iseg][ibond];
             jcomp=Unk2Comp[jseg];
             y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_cav[2],xi_cav[3]);
             dydxi2 = dy_dxi2_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_cav[2],xi_cav[3]);
             dydxi3 = dy_dxi3_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],xi_cav[2],xi_cav[3]);
             term_kseg += Fac_overlap[icomp][jcomp]*(PI/12.0)*(rho_seg[iseg]/y)*(dydxi2*sig2+dydxi3*sig3); 
          }
        }
  }
  return(term_kseg);
}
/****************************************************************************/
/* WTC_overlap: compute some heuristic scaling constants for WTC functionals when
      bond lengths are shorter than Sigmas and the particles are overlapping */
void WTC_overlap()
{
   int icomp,jcomp;

   /* adjust theory for overlapping bonded spheres */
   for (icomp=0;icomp<Ncomp;icomp++){
     for (jcomp=0;jcomp<Ncomp;jcomp++){
       if (Sigma_ff[icomp][icomp]>=Sigma_ff[jcomp][jcomp]){
         if (Bond_ff[icomp][jcomp] >= 0.5*(Sigma_ff[icomp][icomp]+Sigma_ff[jcomp][jcomp])){
            Fac_overlap[icomp][jcomp]=1.0;
         }
         else if (Bond_ff[icomp][jcomp] <= 0.5*fabs(Sigma_ff[icomp][icomp]-Sigma_ff[jcomp][jcomp])){
            Fac_overlap[icomp][jcomp]=0.0;
         }
         else{
            Fac_overlap[icomp][jcomp]=(Bond_ff[icomp][jcomp]/Sigma_ff[jcomp][jcomp]
                                       -0.5*(Sigma_ff[icomp][icomp]-Sigma_ff[jcomp][jcomp])/Sigma_ff[jcomp][jcomp]);
         }
       }
     }
   }
   Fac_overlap[1][1]=Fac_overlap[0][1]*Fac_overlap[0][1];

   for (icomp=0;icomp<Ncomp;icomp++){
     Fac_overlap_hs[icomp]=0.0;
     for (jcomp=0;jcomp<Ncomp;jcomp++){
       if (Sigma_ff[icomp][icomp]<Sigma_ff[jcomp][jcomp]) Fac_overlap[icomp][jcomp]=Fac_overlap[jcomp][icomp];
       if (Fac_overlap[icomp][jcomp]> Fac_overlap_hs[icomp]) Fac_overlap_hs[icomp]=Fac_overlap[icomp][jcomp];
     }
   }

   /* use this to turn off all of the overlap nonsense */
   Fac_overlap_hs[0]=1.;
   Fac_overlap_hs[1]=1.;
   Fac_overlap[0][0]=1.0;
   Fac_overlap[1][0]=1.0;
   Fac_overlap[0][1]=1.0;
   Fac_overlap[1][1]=1.0;

    return;
}
/****************************************************************************/
   double y_cav(double sigma_1,double sigma_2,double xi_2, double xi_3)
{           
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,y;
                   
   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3);
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;
          
   y = 1./(1.-xi_3) +(3.*sig_m/sig_sum)*(xi_2/one_m_xi3_sq)+
       2.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2*xi_2/one_m_xi3_cb);
                 
   return y;     
}
/********************************************************************************************/
   double dy_dxi2_cav(double sigma_1,double sigma_2,double xi_2, double xi_3)
{                
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,dy_dxi2;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3); 
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3); 
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;
 
  dy_dxi2 = (3.*sig_m/sig_sum)*(1.0/one_m_xi3_sq)+
       4.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2/one_m_xi3_cb);
                 
   return dy_dxi2;
}
/********************************************************************************************/
   double dy_dxi3_cav(double sigma_1,double sigma_2,double xi_2, double xi_3)
{
   double one_m_xi3_sq,one_m_xi3_cb,sig_sum,sig_m,dy_dxi3,one_m_xi3_4th;

   one_m_xi3_sq= (1.-xi_3)*(1.-xi_3);
   one_m_xi3_cb= (1.-xi_3)*(1.-xi_3)*(1.-xi_3);
   one_m_xi3_4th= (1.-xi_3)*one_m_xi3_cb;
   sig_sum=sigma_1+sigma_2;
   sig_m=sigma_1*sigma_2;

  dy_dxi3 = 1./one_m_xi3_sq +(6.*sig_m/sig_sum)*(xi_2/one_m_xi3_cb)+
       6.*(sig_m/sig_sum)*(sig_m/sig_sum)*(xi_2*xi_2/one_m_xi3_4th);

   return dy_dxi3;
}
/*******************************************************************************/
/* compute_bulk_nonlocal_wtc_properties: compute some additional bulk properties we
   need to carry around the calculation. These are based on the input densities
   and particle sizes. */
void compute_bulk_nonlocal_wtc_properties(char *output_file1)
{
  int i,loc_inode,loc_i,inode_box,inode,ijk[3],icomp,idim,iunk,printproc;
  int ibond,iseg,jseg,pol_number,type_jseg,nloop,iloop;
  double vol,area,x_dist;
  FILE *fp2=NULL;
  if (Proc==0 && output_file1 !=NULL) printproc = TRUE;
  else printproc=FALSE;
  if (printproc) {
    if( (fp2 = fopen(output_file1,"a+"))==NULL) {
      printf("Can't open file %s\n", output_file1);
      exit(1);
    }
  }

  /* compute bulk nonlocal densities needed for Wertheim-Tripathi-Chapman functionals */
  for (i=0; i<4; i++){
     Xi_cav_b[i]=0.0;
     Xi_cav_LBB[i]=0.0;
     Xi_cav_RTF[i]=0.0;
  }
  for (icomp=0;icomp<Ncomp;icomp++){
     for (i=0;i<4;i++){
        Xi_cav_b[i]+=(PI/6.0)*Rho_b[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
        if (Type_interface!=UNIFORM_INTERFACE){
          Xi_cav_LBB[i]+=(PI/6.0)*Rho_b_LBB[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
          Xi_cav_RTF[i]+=(PI/6.0)*Rho_b_RTF[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
        }
     }
  }
  if (printproc){
     fprintf(fp2,"Xi_cavity_bulk, LBB, and RTF variables for WTC polymer run are:\n");
     fprintf(fp2,"\t i  Xi_cav_b[i]  Xi_cav_LBB[i]  Xi_cav_RTF[i]\n");
     for (i=0;i<4;i++) fprintf(fp2,"\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
              Xi_cav_b[i], Xi_cav_LBB[i], Xi_cav_RTF[i]);
     if (Iwrite==VERBOSE){
        printf("Xi_cav_bulk, LBB, and RTF variables for WTC polymer run are:\n");
        printf("\t i  Xi_cav_b[i]  Xi_cav_LBB[i]  Xi_cav_RTF[i]\n");
         for (i=0;i<4;i++) printf("\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
              Xi_cav_b[i], Xi_cav_LBB[i], Xi_cav_RTF[i]);
     }
  }

  if (Type_poly==WTC){  /* don't need these for WJDC functionals */
     for (ibond=0; ibond<Nbonds; ibond++){
       BondWTC_b[ibond]=NO_BOND_PAIR;
       BondWTC_LBB[ibond]=NO_BOND_PAIR;
       BondWTC_RTF[ibond]=NO_BOND_PAIR;
     }

     for (ibond=0; ibond<Nbonds; ibond++){
        iseg=Unk_to_Seg[ibond];
        pol_number=Unk_to_Poly[ibond];
        jseg=Bonds[pol_number][iseg][Unk_to_Bond[ibond]]/*+SegChain2SegAll[pol_number][0]*/;
        type_jseg=Type_mer[pol_number][jseg];
        jseg=SegChain2SegAll[pol_number][jseg];
        BondWTC_b[ibond]=Rho_seg_b[jseg];
        if (Type_interface!=UNIFORM_INTERFACE){
           BondWTC_LBB[ibond]=Rho_seg_LBB[jseg];
           BondWTC_RTF[ibond]=Rho_seg_RTF[jseg];
        }
     }
  }
  if (Type_poly==WTC && printproc){
     fprintf(fp2,"BondWTC_bulk, LBB, and RTF variables for WTC polymer run are:\n");
     fprintf(fp2,"\t i  BondWTC_b[i]  BondWTC_LBB[i]  BondWTC_RTF[i]\n");
     for (i=0;i<Nbonds;i++) fprintf(fp2,"\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
              BondWTC_b[i], BondWTC_LBB[i], BondWTC_RTF[i]);
     if (Iwrite==VERBOSE){
        printf("BondWTC_bulk, LBB, and RTF variables for WTC polymer run are:\n");
        printf("\t i  BondWTC_b[i]  BondWTC_LBB[i]  BondWTC_RTF[i]\n");
        for (i=0;i<Nbonds;i++) printf("\t %d \t %9.6f \t %9.6f \t %9.6f\n", i,
              BondWTC_b[i], BondWTC_LBB[i], BondWTC_RTF[i]);
     }
  }
  if (printproc) fclose(fp2);

  return;
}
/*********************************************************************************************/
