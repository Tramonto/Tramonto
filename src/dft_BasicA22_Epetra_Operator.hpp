
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

#ifndef DFT_BASICA22_EPETRA_OPERATOR_HPP
#define DFT_BASICA22_EPETRA_OPERATOR_HPP

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Ifpack.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <map>
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"

//! dft_BasicA22_Epetra_Operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! This version accepts a collection of equations as a 2-by-2 block set.  The first block is typically the attractions (referred to as CMS).
	The second is the primitive densities.  All of these equations are coalesced into a single matrix and an algebraic preconditioner (specified
	by the usual Tramonto input deck parameters) is computed and used for preconditioning the Schur complement equations:
 
    S = A22 - A21*inv(A11)*A12, Prec(S) = Prec(A22)
    where any influence of A21*inv(A11)*A12 is ignored in the preconditioner.
*/    

class dft_BasicA22_Epetra_Operator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors.
    //! Builds a coalesced sparse matrix from the 2-by-2 block system of cmsMap and densityMap variables
  /* dft_BasicA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map, int * options, double * params);*/

  dft_BasicA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map, Teuchos::ParameterList * parameterList);
  //@}
  //@{ \name Assembly methods.


  virtual int initializeProblemValues();
  virtual int insertMatrixValue(int rowGID, int colGID, double value);
  //! Calling this method with a given rowGID will flag the row as a constant value row and values will not be set to zero in subsequent reuses of the object.
	virtual void assertRowIsStatic(int rowGID) { isStaticRow_[rowGID] = true; return;}
  virtual int finalizeProblemValues();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_BasicA22_Epetra_Operator();
  //@}
  
  //@{ \name Atribute get methods.

  //! Returns an Epetra_Operator pointer that is actually the \e this object, since this class implements Epetra_Operator.
  virtual Epetra_Operator * getA22Inv() {return(this);}
  //@}
  
  //@{ \name Atribute set methods.

    //! Unsupported feature, returns -1.
  int SetUseTranspose(bool UseTranspose){return(-1);};
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a dft_BasicA22_Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of an inverse dft_BasicA22_Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param In
    X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - An Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
  */ 
  double NormInf() const {return(0.0);};

  //! Check for inconsistencies in operators.
  /* \param verbose (In) Print the residual of inv(A22)*A22*x_random.
     
     \return Returns 0 if residual is "small", otherwise it returns -1.
  */ 
  int Check(bool verbose) const;

  //@}
  
  //@{ \name Atribute access functions

  //! Returns a character string describing the operator
  const char * Label() const{return(Label_);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const{return(block2Map_.Comm());};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const {return(block2Map_);};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const {return(block2Map_);};
  //@}
  
protected:

  int F_location_;
  int insertRow();
  Epetra_Map cmsMap_;
  Epetra_Map densityMap_;
  Epetra_Map block2Map_;
  Teuchos::ParameterList * parameterList_;
  //int * options_;
  //double * params_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A22Matrix_;
  char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool firstTime_;
  int curRow_;
  std::map<int, double> curRowValues_;
  Epetra_IntSerialDenseVector indices_;
  Epetra_SerialDenseVector values_;
  std::vector<bool> isStaticRow_;
  Ifpack factory_;
};

#endif /* DFT_BASICA22_EPETRA_OPERATOR_HPP */
