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

#ifndef DFT_2X2_SCHUR_EPETRA_OPERATOR_H
#define DFT_2X2_SCHUR_EPETRA_OPERATOR_H

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
#include "Epetra_Operator.h"

//! dft_2x2_schur_epetra_operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2-by-2 block operator for Tramonto polymer and explicit non-local density problems.
*/    

class dft_2x2_schur_epetra_operator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2-by-2 block system

  dft_2x2_schur_epetra_operator(Epetra_Operator * A11inv, Epetra_Operator * A12, 
			    Epetra_Operator * A21, Epetra_Operator * A22);
  //@}
  //@{ \name Destructor.
    //! Destructor
  ~dft_2x2_schur_epetra_operator();
  //@}
  
  //@{ \name Atribute set methods.

    //! Unsupported feature, returns -1.
  int SetUseTranspose(bool UseTranspose){return(-1);};
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a dft_2x2_schur_epetra_operator applied to a Epetra_MultiVector X in Y.
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
  
  //! Compute inf-norm residual for global system over all Rhs.
  /*! Compute the inf-norm global residual of the original 2x2 system.
    \param A11 (In) (1,1) block (non-inverted) of the system.
    \param B1 (In) (1) segment of RHS.
    \param B2 (In) (2) segment of RHS.
    \param X1 (In) (1) segment of LHS.
    \param X2 (In) (2) segment of LHS.
  */

  double ComputeGlobalResidual(const Epetra_Operator & A11, const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
			       const Epetra_MultiVector& X1, const Epetra_MultiVector& X2) const;
  
  //! Compute the Rhs given the block matrices and Lhs.
  /*! Compute a Rhs for original 2x2 system.
    \param A11 (In) (1,1) block (non-inverted) of the system.
    \param X1 (In) (1) segment of LHS.
    \param X2 (In) (2) segment of LHS.
    \param B1 (Out) (1) segment of RHS.
    \param B2 (Out) (2) segment of RHS.
  */

  int ComputeGlobalRhs(const Epetra_Operator & A11, const Epetra_MultiVector& X1, const Epetra_MultiVector& X2, 
			  Epetra_MultiVector& B1, Epetra_MultiVector& B2) const;
  
  //! Update block matrices.
  int UpdateBlocks(Epetra_Operator * A11inv, Epetra_Operator * A12, 
		   Epetra_Operator * A21, Epetra_Operator * A22) {A11inv_ = A11inv; A12_ = A12; A21_ = A21; A22_ = A22;}
  
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
  const Epetra_Comm & Comm() const{return(A11inv_->Comm());};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(A22_->DomainMap());};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(A22_->DomainMap());};
  //@}
  

  Epetra_Operator * A11inv_; /*!< The inverse of the 1,1 block of the 2 by 2 block matrix */
  Epetra_Operator * A12_; /*!< The 1,2 block of the 2 by 2 block matrix */
  Epetra_Operator * A21_; /*!< The 2,1 block of the 2 by 2 block matrix */
  Epetra_Operator * A22_; /*!< The 2,2 block of the 2 by 2 block matrix */
  char * Label_; /*!< Description of object */
};

#endif /* DFT_2X2_SCHUR_EPETRA_OPERATOR_H */
