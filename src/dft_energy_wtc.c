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

#include "dft_energy_wtc.h"

/****************************************************************************/
double integrand_WTC_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i,wtc_term,y;
     int icomp,iseg,unk_xi2,unk_xi3,ibond,jseg,jcomp,unk_bond;
     integrand=0.0;

     icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     iseg = iunk-Phys2Unk_first[DENSITY];
     rho_i = x[iunk][inode_box];

     if (rho_i > 1.e-9){

         unk_xi2 = Phys2Unk_first[CAVWTC];
         unk_xi3 = Phys2Unk_first[CAVWTC]+1;
         for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
               jseg = Bonds_SegAll[iseg][ibond];
               jcomp = Unk2Comp[jseg];
               y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],
                          x[unk_xi2][inode_box],x[unk_xi3][inode_box]);
               unk_bond = Poly_to_Unk_SegAll[iseg][ibond]+ Phys2Unk_first[BONDWTC];
               integrand += 1.0-Fac_overlap[icomp][jcomp]*log(y)-log(x[unk_bond][inode_box]) ;
         }
         integrand *= 0.5*rho_i;
     }
     return(integrand);
}
/****************************************************************************/
double integrand_WTC_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_i,wtc_term,y;
     int icomp,iseg,ibond,jseg,jcomp;
     integrand=0.0;

     icomp = Unk2Comp[iunk-Phys2Unk_first[DENSITY]];
     iseg = iunk-Phys2Unk_first[DENSITY];
     rho_i = Rho_seg_b[iseg];

     for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
           jseg = Bonds_SegAll[iseg][ibond];
           jcomp = Unk2Comp[jseg];
           y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],
                        Xi_cav_b[2],Xi_cav_b[3]);

           integrand += 1.0-Fac_overlap[icomp][jcomp]*log(y)-log(BondWTC_b[Poly_to_Unk_SegAll[iseg][ibond]]) ;
     }
     integrand *= 0.5*rho_i;

     return(integrand);
}
/****************************************************************************/

