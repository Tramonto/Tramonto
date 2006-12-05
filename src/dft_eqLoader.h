
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

#ifndef DFT_EQLOADER_H
#define DFT_EQLOADER_H

//! dft_eqLoader: Parent class for all equation (physics) loading classes

class dft_loadManager:  {
      
 public:

  //@{ \name Constructors.
  dft_eqLoader(int iBlock, int numBlocks, int block2Phys[NEQ_TYPE], int block2Nunk[NEQ_TYPE]);

  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_eqLoader();
  //@}

  //@}
  //@{ \name Public Method
    //! Top-level call to load the entire residual/matrix for this phsics
    //! (This is not virtual since there are several data items that must be set.)
  void loadAll(double **x, int iter, int resid_only_flag);
  //@}
  
 protected:
   virtual void loadDensity()=0;
   virtual void loadRhoBarRosen()=0;
   virtual void loadPoisson()=0;
   virtual void loadDiffusion()=0;
   virtual void loadCavityWTC()=0;
   virtual void loadBondWTC()=0;
   virtual void loadCMSField()=0;
   virtual void loadCMSG()=0;
   virtual void loadNewPhys()=0;

 protected:
   const int iBlock_;
   const int numBlocks_;
   const int block2Phys_[NEQ_TYPE];
   const int block2Nunk_[NEQ_TYPE];

   int jBlock_;
   double **x;
   int iter;
   int resid_fill_flag;
};

#endif
