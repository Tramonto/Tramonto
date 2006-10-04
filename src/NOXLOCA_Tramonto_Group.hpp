// $Id$ 
// $Source$ 

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

#ifndef NOXLOCA_TRAMONTO_GROUP_H
#define NOXLOCA_TRAMONTO_GROUP_H

#include "LOCA_MultiContinuation_AbstractGroup.H"  // base class

#include "NOX_Common.H"             // class data element (string)
#include "NOXLOCA_Tramonto_Vector.hpp"	    // class data element

//extern "C" {void fill_resid_and_matrix_control(double **, int, int);}

// Forward declares
namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace NOXLOCA {
  namespace Tramonto {

    class Group : public virtual NOX::Abstract::Group {
    //class Group : public virtual LOCA::MultiContinuation::AbstractGroup {

    public:

      //! Constructor
      Group(NOXLOCA::Tramonto::Vector& xOwned, double** xBox);

      //! Copy constructor
      Group(const NOXLOCA::Tramonto::Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      ~Group();

      //Start methods for NOX::ABstract::Group

      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      //! See above.
      NOX::Abstract::Group& operator=(const NOXLOCA::Tramonto::Group& source);

      /** @name "Compute" functions. */
      //@{

      void setX(const NOX::Abstract::Vector& y);
      //! See above
      void setX(const NOXLOCA::Tramonto::Vector& y);

      void computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step);
      //! See above.
      void computeX(const NOXLOCA::Tramonto::Group& grp, const NOXLOCA::Tramonto::Vector& d, double step);

      NOX::Abstract::Group::ReturnType computeF();

      NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

      //@}

      /** @name Jacobian operations.
       *
       * Operations using the Jacobian matrix. These may not be defined in
       * matrix-free scenarios. */

      //@{
  
      NOX::Abstract::Group::ReturnType applyJacobian(const NOXLOCA::Tramonto::Vector& input, 
						     NOXLOCA::Tramonto::Vector& result) const;

      //! See above
      NOX::Abstract::Group::ReturnType applyJacobian(const NOX::Abstract::Vector& input, 
						     NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& params, 
							    const NOXLOCA::Tramonto::Vector& input, 
				Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& params, 
							    const NOX::Abstract::Vector& input, 
							    NOX::Abstract::Vector& result) const;

      //@}

      /** @name "Is" functions
       *
       * Checks to see if various objects have been computed. Returns true
       * if the corresponding "compute" function has been called since the
       * last update to the solution vector (via instantiation or
       * computeX). */
      //@{

      bool isF() const;
      bool isJacobian() const;
      bool isNewton() const;

      //@}

      /** @name "Get" functions 
       *
       * Note that these function do not check whether or not the vectors
       * are valid. Must use the "Is" functions for that purpose. */
      //@{

      const NOX::Abstract::Vector& getX() const;

      const NOX::Abstract::Vector& getF() const;
  
      double getNormF() const;

      const NOX::Abstract::Vector& getNewton() const;

      const NOX::Abstract::Vector& getGradient() const;

      //@}

      virtual Teuchos::RefCountPtr<NOX::Abstract::Group> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Print out the group
      void print() const;

      //Start methods for LOCA::MultiContinuation::ABstractzGroup

    protected:

      //! resets the isValid flags to false
      void resetIsValid();

    protected:

      /** @name Vectors */
      //@{
      //! Solution vector.
      NOXLOCA::Tramonto::Vector xVector;
      //! Right-hand-side vector (function evaluation).
      NOXLOCA::Tramonto::Vector fVector;
      //! Newton direction vector.
      NOXLOCA::Tramonto::Vector newtonVector;
      //! Tramonto Overlap Vector
      double **xBox;
      //@}

      /** @name IsValid flags 
       *  
       * True if the current solution is up-to-date with respect to the
       * currect xVector. */
      //@{
      bool isValidF;
      bool isValidJacobian;
      bool isValidNewton;
      //@}
  
      //! Norm of F
      double normF;
    };

  } // namespace Tramonto
} // namespace NOX


#endif
