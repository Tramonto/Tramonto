
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
