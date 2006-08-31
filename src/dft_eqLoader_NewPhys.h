
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

#ifndef DFT_EQLOADERNEWPHYS_H
#define DFT_EQLOADERNEWPHYS_H

#include <dft_eqLoader.h>

//! dft_eqLoader: Parent class for all equation (physics) loading classes

class dft_eqLoader_NewPhys: public class dft_eqLoader  {
      
 public:

  //@{ \name Constructors.
  dft_eqLoader_NewPhys(int iBlock, int numBlocks, int block2Phys[NEQ_TYPE], int block2Nunk[NEQ_TYPE]);

  //@{ \name Destructor.
    //! Destructor
   ~dft_eqLoader_NewPhys();
  //@}

  //@}
  //@{ \name Public Method
    //! Top-level call to load the entire residual/matrix for this phsics
    //! This is implemented in parent class!
  //!void loadAll(double **x, int iter, int resid_only_flag);
  //@}
  
 private:
   void loadDensity();
   void loadRhoBarRosen();
   void loadPoisson();
   void loadDiffusion();
   void loadCavityWTC();
   void loadBondWTC();
   void loadCMSField();
   void loadCMSG();
   void loadNewPhys();

   void printBlockInfo(string phys);
};

#endif
