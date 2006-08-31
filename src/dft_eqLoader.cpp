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
