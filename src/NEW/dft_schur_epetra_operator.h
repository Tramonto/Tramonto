
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

#ifndef DFT_SCHUR_EPETRA_OPERATOR_H
#define DFT_SCHUR_EPETRA_OPERATOR_H

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"

//! dft_schur_epetra_operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2-by-2 block operator for Tramonto polymer and explicit non-local density problems.
*/    

class dft_schur_epetra_operator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2-by-2 block system

  dft_schur_epetra_operator(Epetra_CrsMatrix * A11, Epetra_CrsMatrix * A12, 
			    Epetra_CrsMatrix * A21, Epetra_CrsMatrix * A22);
  //@}
  //@{ \name Destructor.
    //! Destructor
  ~dft_schur_epetra_operator();
  //@}
  
  //@{ \name Atribute set methods.

    //! Unsupported feature, returns -1.
  int SetUseTranspose(bool UseTranspose){return(-1);};
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a dft_schur_epetra_operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Unsupported feature, returns -1.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {return(-1);};
  
  //! Generate RHS.
  int ComputeRHS(const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
		 Epetra_MultiVector& B) const;
  
  //! Generate X1.
  int ComputeX1(const Epetra_MultiVector& B1, const Epetra_MultiVector& X2, Epetra_MultiVector& X1) const;
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
  */ 
  double NormInf() const {return(0.0);};
  //@}
  
  //@{ \name Atribute access functions

  //! Returns a character string describing the operator
  char * Label() const{return(Label_);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(A11_->Comm());};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(A22_->DomainMap());};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(A22_->DomainMap());};
  //@}
  

  Epetra_CrsMatrix * A11_; /*!< The 1,1 block of the 2 by 2 block matrix */
  Epetra_CrsMatrix * A12_; /*!< The 1,2 block of the 2 by 2 block matrix */
  Epetra_CrsMatrix * A21_; /*!< The 2,1 block of the 2 by 2 block matrix */
  Epetra_CrsMatrix * A22_; /*!< The 2,2 block of the 2 by 2 block matrix */
  char * Label_; /*!< Description of object */
};

#endif /* DFT_SCHUR_EPETRA_OPERATOR_H */
