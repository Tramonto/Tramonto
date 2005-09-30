
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include <dft_loadManager.hpp>


dft_loadManager::dft_loadManager(int* phys2Nunk) : numBlocks_(0)
{
   for (int i=0; i<NEQ_TYPE; i++) {
     if (phys2Nunk[i] != 0) {
       block2Phys_[numBlocks_] = i;   // Physics type
       block2Nunk_[numBlocks_] = phys2Nunk[i]; // Number of unknowns of this type
       allocatePhysLoader(numBlocks_, i); // Assign a specific physics loader to this block
       numBlocks_++;
     }
   }
}

dft_loadManager::~dft_loadManager()
{ 
   for (i=0; i<numBlocks_; i++)   delete physLoader_[i];
}

void dft_loadManager::loadMatrix(double **x, int iter, int resid_only_flag)
{
   for (i=0; i<numBlocks_; i++)
     physLoader_[i]->loadAll(x, iter, resid_only_flag);
}

void dft_loadManager::allocatePhysLoader(int iBlock, int iPhys)
{
   switch (iphys) {
     case DENSITY:
        if (Type_poly_TC==FALSE)
          physLoader_[iBlock] = new dft_eqLoader_FluidEulerLagrange(iBlock, numBlocks_, block2Phys_);
        if (Type_poly_TC==FALSE)
          physLoader_[iBlock] = new dft_eqLoader_PolyEulerLagrange(iBlock, numBlocks_, block2Phys_);
        break;
     case RHOBAR_ROSEN:
        physLoader_[iBlock] = new dft_eqLoader_RhoBarRosen(iBlock, numBlocks_, block2Phys_); break;
     case POISSON:   
        physLoader_[iBlock] = new dft_eqLoader_Poisson(iBlock, numBlocks_, block2Phys_); break;
     case DIFFUSION:
        physLoader_[iBlock] = new dft_eqLoader_Diffusion(iBlock, numBlocks_, block2Phys_); break;
     case CAVITY_WTC:
        physLoader_[iBlock] = new dft_eqLoader_CavityWTC(iBlock, numBlocks_, block2Phys_); break;
     case BOND_WTC:
        physLoader_[iBlock] = new dft_eqLoader_BondWTC(iBlock, numBlocks_, block2Phys_); break;
     case CMS_FIELD:
        physLoader_[iBlock] = new dft_eqLoader_CMSField(iBlock, numBlocks_, block2Phys_); break;
     case CMS_G:
        physLoader_[iBlock] = new dft_eqLoader_CMSG(iBlock, numBlocks_, block2Phys_); break;
     case NEWPHYS:
        physLoader_[iBlock] = new dft_eqLoader_NewPhys(iBlock, numBlocks_, block2Phys_); break;
     default:  cerr << "dft_loadManager::allocatePhysLoader:  Unknown Phys type "
                    << iPhys << " " < iBlock << endl;  
               exit(-1);
   }
}

