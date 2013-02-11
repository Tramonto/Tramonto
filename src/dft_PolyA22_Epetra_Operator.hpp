
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

#ifndef DFT_POLYA22_EPETRA_OPERATOR_H
#define DFT_POLYA22_EPETRA_OPERATOR_H

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

//! dft_PolyA22_Epetra_Operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/    

class dft_PolyA22_Epetra_Operator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads system
  /* dft_PolyA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map, int * options, double * params);*/

  dft_PolyA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map, Teuchos::ParameterList * parameterList);
  //@}
  //@{ \name Assembly methods.

  //! Assert that field dependence on primitive densities is linear; manager will not reset values between nonlinear solves.
  /*! This method can be called to assert that the field variable dependence on primitive densities does change from one linear solve to the next.
      In this case, we can avoid filling the associated matrix coefficients.  Calling this method with "true" will cause the problem manager not reset
      the matrix coefficients for this block and to ignore any values that are submitted for entry in this block.
     \param isLinear (In) Set to true if the field dependence is linear on primitive densities.
     \warning This method can be called at any time, but should be called before the initializeValues() method is called for the second solve; By default the manager assumes that the relationship is non-linear, so values will be reset to zero and must be refilled before each linear solve.
  */
  int setFieldOnDensityIsLinear(bool isLinear) { 
    isFLinear_ = isLinear;
    return(0);
  }

  virtual int initializeProblemValues();
  virtual int insertMatrixValue(int rowGID, int colGID, double value, int blockColFlag);
  virtual int finalizeProblemValues();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA22_Epetra_Operator();
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

    //! Returns the result of a dft_PolyA22_Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of an inverse dft_PolyA22_Epetra_Operator applied to a Epetra_MultiVector X in Y.
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
  virtual int Check(bool verbose) const;

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
  Teuchos::RefCountPtr<Epetra_CrsMatrix> cmsOnDensityMatrix_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> cmsOnCmsMatrix_;
  Teuchos::RefCountPtr<Epetra_Vector> densityOnDensityMatrix_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> densityOnCmsMatrix_;
  const char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool isFLinear_;
  bool firstTime_;
  int curRow_;
  std::map<int, double> curRowValuesCmsOnDensity_, curRowValuesCmsOnCms_, curRowValuesDensityOnCms_;
  Epetra_IntSerialDenseVector indicesCmsOnDensity_, indicesCmsOnCms_, indicesDensityOnCms_;
  Epetra_SerialDenseVector valuesCmsOnDensity_, valuesCmsOnCms_, valuesDensityOnCms_;
  Teuchos::ParameterList IFList_;
  Teuchos::RCP<Ifpack_Preconditioner> IFPrec;
  string IFPrecType; // incomplete LU
  int IFOverlapLevel;
};

#endif /* DFT_POLYA22_EPETRA_OPERATOR_H */
