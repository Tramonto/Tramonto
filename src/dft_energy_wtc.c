#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

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

         unk_xi2 = Phys2Unk_first[CAVITY_WTC];
         unk_xi3 = Phys2Unk_first[CAVITY_WTC]+1;
         for (ibond=0;ibond<Nbonds_SegAll[iseg];ibond++){
               jseg = Bonds_SegAll[iseg][ibond];
               jcomp = Unk2Comp[jseg];
               y = y_cav(Sigma_ff[icomp][icomp],Sigma_ff[jcomp][jcomp],
                          x[unk_xi2][inode_box],x[unk_xi3][inode_box]);
               unk_bond = Poly_to_Unk_SegAll[iseg][ibond]+ Phys2Unk_first[BOND_WTC];
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

