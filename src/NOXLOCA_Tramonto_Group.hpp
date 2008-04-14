// $Id$ 
// $Source$ 

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

#ifndef NOXLOCA_TRAMONTO_GROUP_H
#define NOXLOCA_TRAMONTO_GROUP_H

#include "LOCA_Abstract_Group.H"  // base class
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Common.H"             // class data element (string)
#include "NOXLOCA_Tramonto_Vector.hpp"	    // class data element

extern "C" {
double calc_free_energy_conwrap(double **xBox);
}

// Forward declares
namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace NOXLOCA {
  namespace Tramonto {

    //class Group : public virtual NOX::Abstract::Group {
    class Group : public virtual LOCA::Abstract::Group {

    public:

      //! Constructor
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& gD, NOXLOCA::Tramonto::Vector& xOwned,
            double** xBox, const LOCA::ParameterVector& pVec,
            const Teuchos::RefCountPtr<Teuchos::ParameterList>& paramList, bool doPicard);

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

      //! Start methods for LOCA::Abstract::Group

      //@{
      void copy(const NOX::Abstract::Group& source) { *this = source; }

      //! Set the parameter vector in the group to p (pVector = p).
      void setParams(const LOCA::ParameterVector& p);

       //! Set parameter indexed by (integer) paramID
      void setParam(int paramID, double val);

      //! Set parameter indexed by (string) paramID
      void setParam(string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group.
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by (integer) paramID
      double getParam(int paramID) const;

      //! Return copy of parameter indexed by (string) paramID
      double getParam(string paramID) const;

      //! Print Solution method for ouytput. printSolution2 is for x2 of phase trans
      void printSolution(const double param) const;
      void printSolution(const NOX::Abstract::Vector& solution, const double param) const;
      void printSolution2(const NOX::Abstract::Vector& solution, const double param) const;

      //! Calculate the free energy at the given X and Param
      double calcFreeEnergy() const;
      //@}

    protected:

      //! resets the isValid flags to false
      void resetIsValid();

      // Box 2 Local conversions, not using linprob_mgr_importr2c()
      void TV2Box(const NOXLOCA::Tramonto::Vector& xTV, double** xB) const;
      void Box2TV(double** xB, NOXLOCA::Tramonto::Vector& xTV) const;

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
      mutable double **xBox;
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

      LOCA::ParameterVector paramVec;
      mutable int contStep; // Printing index
      mutable bool secondSolution;
      const Teuchos::RefCountPtr<LOCA::GlobalData> globalData;
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& paramList;

      bool doPicard;
    };

  } // namespace Tramonto
} // namespace NOX


#endif
