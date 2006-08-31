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

#ifndef DFT_LOADMANAGER_H
#define DFT_LOADMANAGER_H

#include <dft_eqLoader.h>

//! dft_loadManager: Top level code for organizing matrix/residual fills

class dft_loadManager:  {
      
 public:

  //@{ \name Constructors.
    //! Figures out what equation blocks are being filled and sets block/Eq numbering

  dft_loadManager(int* Phys2Nunk);
  //@}
  //@{ \name Top-level call to load the entire residual/matrix
  int loadMatrix(double **x, int iter, int resid_only_flag);
  //@}
  //@{ \name Destructor.
    //! Destructor
  ~dft_loadManager();
  //@}
  
private:
  //! Function that allocates the specific dft_eqLoader classes for
  // the requested physics for this problem.
  void allocatePhysLoader(int iBlock, int iPhys)

  int numBlocks_;
  int block2Phys_[NEQ_TYPE];
  int block2Nunk_[NEQ_TYPE];
  dft_eqLoader* physLoader_[NEQ_TYPE];
  
};

#endif 
