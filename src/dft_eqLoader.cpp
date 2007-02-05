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

  // Top Level class for loading equations / physics. Generic implementations
  // in common to all physics types are implemented here.

dft_eqLoader::dft_eqLoader(int iBlock, int numBlocks, int block2Phys[NEQ_TYPE])
   : iBlock_(iBlock), numBlocks_(numBlocks), block2Phys_(block2Phys)
{}

virtual dft_eqLoader::~dft_eqLoader() {}

void dft_eqLoader::loadAll(double **x_, int iter_, int resid_only_flag_)
{
  // Make these accessible by all fill functions
  x = x_;
  iter=iter_;
  resid_only_flag=resid_only_flag_;

  //AGS: Consider passing just the corect column(s) of x to
  //each block fill routine by pointing into x, &x[Phys2Unk_first[i]])

  for (i=0; i<numBlocks_; i++) {
    jBlock_ = i; 
    switch (block2Phys[i]) {
     case DENSITY:      loadDensity();     break;
     case RHOBAR_ROSEN: loadRhoBarRosen(); break;
     case POISSON:      loadPoisson();     break;
     case DIFFUSION:    loadDiffusion();   break;
     case CAVITY_WTC:   loadCavityWTC();   break;
     case BOND_WTC:     loadBondWTC();     break;
     case CMS_FIELD:    loadCMSField();    break;
     case CMS_G:        loadCMSG();        break;
     case NEWPHYS:      loadNewPhys();     break;
    }
  }
}
