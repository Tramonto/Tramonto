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

#include "dft_energy_wjdc.h"

/****************************************************************************/
double integrand_WJDC_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int iseg,unk_field,ibond,count_ends,pol_num,icomp;

     iseg = iunk-Phys2Unk_first[DENSITY];
     icomp = Unk2Comp[iseg];
     rho_i = x[iunk][inode_box];
     if (Type_poly==WJDC2 || Type_poly==WJDC3) unk_field=Phys2Unk_first[WJDC_FIELD]+icomp;
     else if (Type_poly==WJDC)                 unk_field=Phys2Unk_first[WJDC_FIELD]+iseg;   
     pol_num=SegAll_to_Poly[iseg];

     count_ends=0;
     for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
        if(Bonds_SegAll[iseg][ibond]==-1) count_ends++;
     }

     if (rho_i > 1.e-9){
        integrand = rho_i*(log(x[unk_field][inode_box]/exp(Scale_fac_WJDC[pol_num][icomp])) 
                      + 0.5*(Nbonds_SegAll[iseg]-count_ends)-1.0);
     }
     return(integrand);
}
/****************************************************************************/
double integrand_WJDC_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int iseg,icomp,ibond,count_ends,pol_num;

     iseg = iunk-Phys2Unk_first[DENSITY];
     icomp=Unk2Comp[iseg];
     rho_i=Rho_seg_b[iseg]; /*= Rho_seg_b[icomp];*/
     pol_num=SegAll_to_Poly[iseg];

     count_ends=0;
     for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
        if(Bonds_SegAll[iseg][ibond]==-1) count_ends++;
     }

     if (rho_i > 1.e-9){
        integrand = rho_i*(log(Field_WJDC_b[icomp]/exp(Scale_fac_WJDC[pol_num][icomp]) )
                     + 0.5*(Nbonds_SegAll[iseg]-count_ends)-1.0);
     }
     return(integrand);
}
/****************************************************************************/

/****************************************************************************/
double integrand_WJDCcomp_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i,bondproduct,scale_term,mu;
     int iseg,unk_field,ibond,count_ends,npol,icomp,jcomp,loop_start,loop_end,i,unk_GQ;

     icomp=iunk-Phys2Unk_first[DENSITY];
     unk_field=Phys2Unk_first[WJDC_FIELD]+icomp;

                           /* for WJDC3 polymers where we are concerned with component densities,                      
                           we need to perform a sum over segment of a given type */
     npol = 0;
     while (Nmer_t[npol][icomp]==0){ npol++; }                     
     loop_start=0;      
     loop_end=Nmer[npol];
     integrand=0.0;

     scale_term=0,0;
     for (jcomp=0;jcomp<Ncomp;jcomp++){
        scale_term-=Scale_fac_WJDC[npol][jcomp]*Nseg_type_pol[npol][jcomp];
     }
     
     for (i=loop_start; i<loop_end; i++){
     if (Type_mer[npol][i] == icomp) {
        
        iseg=SegChain2SegAll[npol][i];

        /* compute the segment density here */
        bondproduct=1.0;
        for(ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
           unk_GQ  = Phys2Unk_first[G_CHAIN] + Poly_to_Unk_SegAll[iseg][ibond];
           bondproduct *= x[unk_GQ][inode_box];
        }
        if (fabs(x[unk_field][inode_box])>1.e-12)
        rho_i=bondproduct*POW_DOUBLE_INT(x[unk_field][inode_box],-(Nbonds_SegAll[iseg]-1));        
        else rho_i=0.0;

/*      rho_i*=prefactor_rho_wjdc(iseg);*/
        if (Type_interface==DIFFUSIVE_INTERFACE){
           mu=x[Phys2Unk_first[DIFFUSION]+npol][inode_box];
        }
        else mu=Betamu_chain[npol];

        rho_i*=exp(mu+scale_term);

        count_ends=0;
        for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
          if(Bonds_SegAll[iseg][ibond]==-1) count_ends++;
        }

        if (rho_i > 1.e-9){
           integrand += rho_i*(log(x[unk_field][inode_box]/exp(Scale_fac_WJDC[npol][icomp])) 
                      + 0.5*(Nbonds_SegAll[iseg]-count_ends)-1.0);
        }
     }
     }
     return(integrand);
}
/****************************************************************************/
double integrand_WJDCcomp_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int iseg,icomp,ibond,count_ends,pol_num,npol,loop_start,loop_end,i;

     icomp=iunk-Phys2Unk_first[DENSITY];

                           /* for WJDC3 polymers where we are concerned with component densities,                      
                           we need to perform a sum over segment of a given type */

     if (Grafted_Logical==FALSE || (Type_poly==WJDC3 && Grafted[Icomp_to_polID[icomp]]==FALSE)){

     npol = 0;
     while (Nmer_t[npol][icomp]==0){ npol++; }                     
     loop_start=0;      
     loop_end=Nmer[npol];
     integrand=0.0;
 
     for (i=loop_start; i<loop_end; i++){
     if (Type_mer[npol][i] == icomp) {

        iseg=SegChain2SegAll[npol][i];
        rho_i = Rho_seg_b[iseg];
        pol_num=SegAll_to_Poly[iseg];

        count_ends=0;
        for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
           if(Bonds_SegAll[iseg][ibond]==-1) count_ends++;
        }

        if (rho_i > 1.e-9){
           integrand += rho_i*(log(Field_WJDC_b[icomp]/exp(Scale_fac_WJDC[pol_num][icomp]) )
                        + 0.5*(Nbonds_SegAll[iseg]-count_ends)-1.0);
        }
     }
     }
     }
     else integrand=0.0;
     return(integrand);
}
/****************************************************************************/
