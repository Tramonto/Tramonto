
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

#ifndef DFT_POLYA22_COULOMB_EPETRA_OPERATOR_H
#define DFT_POLYA22_COULOMB_EPETRA_OPERATOR_H

class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_BlockMap;
class Epetra_Comm;
class Epetra_LinearProblem;
class AztecOO;
#include "dft_PolyA22_Epetra_Operator.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_LinearProblem.h"
#include "Ifpack.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <map>
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Amesos.h" //may not need from here down
#include "Amesos_BaseSolver.h"
#include "EpetraExt_Reindex_LinearProblem.h"
#include "AztecOO.h"
#include "dft_direct_solver_const.h"
#include "EpetraExt_MultiVectorOut.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_include.h"

//! dft_PolyA22_Epetra_Operator: An implementation of the Epetra_Operator class for Tramonto Schur complements with Coulomb effects.
/*! Special 2*numBeads by 2*numBeads (plus Coulomb) for Tramonto polymer problems.
*/    

class dft_PolyA22_Coulomb_Epetra_Operator: public virtual dft_PolyA22_Epetra_Operator {
      
 public:

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads (plus Coulomb) system

  /* dft_PolyA22_Coulomb_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & densityMap, const Epetra_Map & cmsDensMap, const Epetra_Map & block2Map, int * options, double * params);*/

  dft_PolyA22_Coulomb_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & densityMap, const Epetra_Map & cmsDensMap, const Epetra_Map & block2Map, Teuchos::ParameterList * parameterList);
  //@}
  //@{ \name Assembly methods.
  int initializeProblemValues();
  int insertMatrixValue(int rowGID, int colGID, double value);
  int finalizeProblemValues();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA22_Coulomb_Epetra_Operator();
  //@}
  
  //@{ \name Atribute get methods.

  //! Returns an Epetra_Operator pointer that is actually the \e this object, since this class implements Epetra_Operator.
  Epetra_Operator * getA22Inv() {return(this);} //keep or comment??
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a dft_PolyA22_Coulomb_Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of an inverse dft_PolyA22_Coulomb_Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param In
    X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - An Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Check for inconsistencies in operators.
  /* \param verbose (In) Print the residual of inv(A22)*A22*x_random.
     
     \return Returns 0 if residual is "small", otherwise it returns -1.
  */ 
  int Check(bool verbose) const; //leave or comment?

  //@}
  
  //@{ \name Atribute access functions

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(block2Map_.Comm());};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(block2Map_);};
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(block2Map_);};
  //@}
  
protected:

  int insertPoissonRow();
  int insertCPRow();
  int insertPDRow();
  Epetra_Map poissonMap_;
  Epetra_Map cmsDensMap_;
  Epetra_Map block2Map_;
  // Epetra_Map reindexedMap_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> poissonMatrix_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> cmsOnPoissonMatrix_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> poissonOnDensityMatrix_;
  Teuchos::ParameterList MLList_;
  int curPoissonRow_;
  std::map<int, double> curPoissonRowValues_;
  int curCPRow_;
  std::map<int, double> curCPRowValues_;
  int curPDRow_;
  std::map<int, double> curPDRowValues_;
};

#endif /* DFT_POLYA22_COULOMB_EPETRA_OPERATOR_H */
