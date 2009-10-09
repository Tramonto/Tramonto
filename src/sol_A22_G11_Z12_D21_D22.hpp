
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

#ifndef SOL_A22_G11_Z12_D21_D22_HPP
#define SOL_A22_G11_Z12_D21_D22_HPP

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <map>
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Ifpack.h"

//! sol_A22_G11_Z12_D21_D22: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/    

class sol_A22_G11_Z12_D21_D22: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors.
    //! Builds a 2-by-2 block Epetra_RowMatrix object for the A22 block of the outer 2-by-2 block global Tramonto linear system.
    /*!
     */

  sol_A22_G11_Z12_D21_D22(const Epetra_Map & map1, const Epetra_Map & map2, const Epetra_Map & A22Map, Teuchos::ParameterList * parameterList);
  //@}
  //@{ \name Assembly methods.

  //! Calling this method with a given rowGID will flag the row as a constant value row and values will not be set to zero in subsequent reuses of the object.
  /*! This method can be called to assert that the specified row of G11 does change from one linear solve to the next.
      In this case, we can avoid filling the associated matrix coefficients on subsequent calls.  Calling this method with "true" will cause the problem manager not reset
      the matrix coefficients for this row and to ignore any values that are submitted for entry in this block.
     \param rowGID (In) Set to true if this row should be define only once.
  */
   void setG11RowStatic(int rowGID) { 
    TEST_FOR_EXCEPTION(!map1_.MyGID(rowGID), std::logic_error, 
                       "You are setting a row of G11 static that is not owned by this matrix"); 
    isG11RowStatic_[map1_.LID(rowGID)] = true; 
    return;
  }

  virtual int initializeProblemValues();
  virtual int insertMatrixValue(int rowGID, int colGID, double value);
  virtual int finalizeProblemValues();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~sol_A22_G11_Z12_D21_D22();
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

    //! Returns the result of a sol_A22_G11_Z12_D21_D22 applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of an inverse sol_A22_G11_Z12_D21_D22 applied to a Epetra_MultiVector X in Y.
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
  virtual const Epetra_Comm & Comm() const{return(A22Map_.Comm());};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const {return(A22Map_);};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const {return(A22Map_);};
  //@}
  
protected:

  int insertRow();
  Epetra_Map map1_;
  Epetra_Map map2_;
  Epetra_Map A22Map_;
  Teuchos::ParameterList * parameterList_;
  Teuchos::RCP<Epetra_CrsMatrix> G11Matrix_;
  Teuchos::RCP<Epetra_Vector> D21Matrix_;
  Teuchos::RCP<Epetra_Vector> D22Matrix_;
  char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool firstTime_;
  int curRow_;
  std::map<int, double> curRowValues_;
  Epetra_IntSerialDenseVector indices_;
  Epetra_SerialDenseVector values_;
  Ifpack factory_;
  Teuchos::RCP<Ifpack_Preconditioner> G11Inverse_;
  std::vector<bool> isG11RowStatic_;

};

#endif /* SOL_A22_G11_Z12_D21_D22_HPP */
