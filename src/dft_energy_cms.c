#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_CMS_freen(int iunk,int inode_box, double **x)
{
     double integrand=0.0,rho_i;
     int icomp,pol_number;

     icomp= iunk-Phys2Unk_first[DENSITY];
     rho_i = x[iunk][inode_box];

     int_stencil(x,inode_box,iunk,POLYMER_CR);

                                            /* note second term is old adsorption */ 
     pol_number=0;
     while (Nmer_t[pol_number][icomp]==0) pol_number++;
     if (rho_i > 0.) integrand = 0.5*rho_i*Temporary_sum-rho_i/Nmer[pol_number];
     return(integrand);
}
/****************************************************************************/
double integrand_CMS_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk,sum;
     int icomp,jcomp,pol_number;

     icomp= iunk-Phys2Unk_first[DENSITY];
     rho_bulk = Rho_b[icomp];

     sum=0;
     for (jcomp=0;jcomp<Ncomp;jcomp++){
       int_stencil_bulk(POLYMER_CR,icomp,jcomp);
       sum+=Temporary_sum*Rho_b[jcomp];
     }

     pol_number=0;
     while (Nmer_t[pol_number][icomp]==0) pol_number++;

     integrand = 0.5*rho_bulk*sum-rho_bulk/Nmer[pol_number];
     return(integrand);
}
/****************************************************************************/
