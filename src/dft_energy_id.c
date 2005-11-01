#include "mpi.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

/****************************************************************************/
double integrand_ideal_gas_freen(int iunk,int inode_box, double **x)
{
     double integrand,rho_i;

     rho_i = x[iunk][inode_box];

     if (rho_i > DENSITY_MIN)
          integrand = rho_i*(log(rho_i)-1.0);
     return(integrand);
}
/****************************************************************************/
double integrand_ideal_gas_freen_bulk(int iunk,int inode_box, double **x)
{
     double integrand,rho_bulk;

     if (Type_poly==WTC) rho_bulk = Rho_seg_b[iunk-Phys2Unk_first[DENSITY]];
     else                rho_bulk = Rho_b[iunk-Phys2Unk_first[DENSITY]];

     integrand = rho_bulk*(log(rho_bulk)-1.0);
     return(integrand);
}
/****************************************************************************/
