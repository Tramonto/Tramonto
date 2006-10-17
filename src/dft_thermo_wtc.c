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

/* ---------------------------------------------------------
Calculate the thermodynamic properties of chain contributions for a 
Wertheim-Tripathi-Chapman bonded fluid.
------------------------------------------------------------*/
#include "dft_globals_const.h"
#include "rf_allo.h"
#include "mpi.h"
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
  if (Type_poly==WTC){  
     for (i=0; i<4; i++){
        Xi_cav_b[i]=0.0;
        Xi_cav_LBB[i]=0.0;
        Xi_cav_RTF[i]=0.0;
     }
     for (icomp=0;icomp<Ncomp;icomp++){
        for (i=0;i<4;i++){
           Xi_cav_b[i]+=(PI/6.0)*Rho_b[icomp]*POW_DOUBLE_INT(Sigma_ff[icomp][icomp],i);
           if (Lsteady_state){
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
        if (Lsteady_state){
           BondWTC_LBB[ibond]=Rho_seg_LBB[jseg];
           BondWTC_RTF[ibond]=Rho_seg_RTF[jseg];
        }
    }
    if (printproc){
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
       fclose(fp2);
     }

  } /* end of Type_poly_WTC rhobars (bulk)*/
  return;
}
/*********************************************************************************************/
