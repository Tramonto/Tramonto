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

#ifndef NOXLOCA_TRAMONTO_PTVECTOR_H
#define NOXLOCA_TRAMONTO_PTVECTOR_H

#include "NOXLOCA_Tramonto_Vector.hpp" // base class

namespace NOXLOCA {

  //! NOX Tramonto support for Phase Transition vectors (length 2n+1)
  namespace Tramonto {

    //! Implementation of NOX::Abstract::Vector
    class PTVector : public NOX::Abstract::Vector {

    private:

      // SHouldnt be called
      PTVector();

    public:	
      //! Construct a vector of length n from given array
      PTVector(const NOXLOCA::Tramonto::Vector &x1, 
               const NOXLOCA::Tramonto::Vector &x2,
               const double& ptp);
    
      //! Copy constructor
      PTVector(const NOXLOCA::Tramonto::PTVector& source, 
	     NOX::CopyType type = NOX::DeepCopy);
      
      //! Destruct PTVector.
      ~PTVector();

      //@{ \name Initialization methods.

      NOX::Abstract::Vector& init(double value);

      //! Initialize every element of this vector with random values
      virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed =  1);

      // derived
      NOX::Abstract::Vector& operator=(const NOXLOCA::Tramonto::PTVector& y);
      NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
  
      // derived
      NOX::Abstract::Vector& abs(const NOXLOCA::Tramonto::PTVector& y);
      NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);
  
      // derived
      NOX::Abstract::Vector& reciprocal(const NOXLOCA::Tramonto::PTVector& y);
      NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);
  
      //@}
  
      //@{ \name Update methods.
  
      // derived
      NOX::Abstract::Vector& scale(double gamma);
  
      // derived
      NOX::Abstract::Vector& scale(const NOXLOCA::Tramonto::PTVector& a);
      NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOXLOCA::Tramonto::PTVector& a, double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOXLOCA::Tramonto::PTVector& a, 
			       double beta, const NOXLOCA::Tramonto::PTVector& b,
			       double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			       double beta, const NOX::Abstract::Vector& b,
			       double gamma = 0.0);
  
      //@}
  
      //@{ \name Creating new Vectors. 
  
      // derived
      Teuchos::RefCountPtr<NOX::Abstract::Vector> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;
  
      //@}
  
      //@{ \name Norms.
  
      // derived
      double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const;
  
      // derived
      double norm(const NOXLOCA::Tramonto::PTVector& weights) const;
      double norm(const NOX::Abstract::Vector& weights) const;
  
      //@}
  
      //@{ \name Inner products
  
      // derived
      double innerProduct(const NOXLOCA::Tramonto::PTVector& y) const;
      double innerProduct(const NOX::Abstract::Vector& y) const;
  
      //@}
  
      // derived
      NOX::size_type length() const;
  
      // derived
      void print(std::ostream& stream) const;

      // Allow access to Tramonto sub-vectors from the PTVector
      NOXLOCA::Tramonto::Vector& X1();
      const NOXLOCA::Tramonto::Vector& X1() const;
      NOXLOCA::Tramonto::Vector& X2();
      const NOXLOCA::Tramonto::Vector& X2() const;
      double&  PTP();
      const double&  PTP() const;

      std::ostream& leftshift(std::ostream& stream) const;

    private:

      NOXLOCA::Tramonto::Vector x1;
      NOXLOCA::Tramonto::Vector x2;
      double ptp;

    };

  } // namespace Tramonto
} // namespace NOX

//! Function for printing
std::ostream& operator<<(std::ostream& stream, const NOXLOCA::Tramonto::PTVector& v);

#endif
