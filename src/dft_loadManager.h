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
