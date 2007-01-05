//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER

#ifndef NOXLOCA_TRAMONTO_VECTOR_H
#define NOXLOCA_TRAMONTO_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class
#include "NOX_Common.H" // for #include<vector>

#include "dft_globals_const.h"

extern "C" {
//work-arounds for C++ C linkage problems
double gmax_double_conwrap(double);
double gsum_double_conwrap(double);
void fill_resid_and_matrix_control_conwrap(double**, int, int);
double** array_alloc_2d_conwrap(unsigned int, unsigned int, unsigned int);
void safe_free_conwrap(void **);
}


namespace NOXLOCA {

  //! NOX Tramonto Vector support
  namespace Tramonto {

    //! Implementation of NOX::Abstract::Vector 
    class Vector : public NOX::Abstract::Vector {

    public:	

      //! Construct an empty vector
      Vector();

      //! Construct a zero vector of length n
      Vector(int n1, int n2);

      //! Construct a vector of length n from given array
      Vector(int n1, int n2, double **v);
    
      //! Copy constructor
      Vector(const NOXLOCA::Tramonto::Vector& source, 
	     NOX::CopyType type = NOX::DeepCopy);
      
      //! Destruct Vector.
      ~Vector();

      //@{ \name Initialization methods.

      NOX::Abstract::Vector& init(double value);

      //! Initialize every element of this vector with random values
      virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed =  1);

      // derived
      NOX::Abstract::Vector& operator=(const NOXLOCA::Tramonto::Vector& y);
      NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
  
      // derived
      NOX::Abstract::Vector& abs(const NOXLOCA::Tramonto::Vector& y);
      NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);
  
      // derived
      NOX::Abstract::Vector& reciprocal(const NOXLOCA::Tramonto::Vector& y);
      NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);
  
      //@}
  
      //@{ \name Update methods.
  
      // derived
      NOX::Abstract::Vector& scale(double gamma);
  
      // derived
      NOX::Abstract::Vector& scale(const NOXLOCA::Tramonto::Vector& a);
      NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOXLOCA::Tramonto::Vector& a, double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0);
  
      // derived
      NOX::Abstract::Vector& update(double alpha, const NOXLOCA::Tramonto::Vector& a, 
			       double beta, const NOXLOCA::Tramonto::Vector& b,
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
      double norm(const NOXLOCA::Tramonto::Vector& weights) const;
      double norm(const NOX::Abstract::Vector& weights) const;
  
      //@}
  
      //@{ \name Inner products
  
      // derived
      double innerProduct(const NOXLOCA::Tramonto::Vector& y) const;
      double innerProduct(const NOX::Abstract::Vector& y) const;
  
      //@}
  
      // derived
      int length() const;
  
      //! Return the i-th element
      double& operator() (int i);

      //! Return the i-th element (const version)
      const double& operator() (int i) const;

      //! Prints out the vector to the specified stream. 
      /*! 
	For example, a vector would appear as
	\f[ \left[ \; 0.1 \; 2.34 \; 5 \; \right] \f] 
	It will be all on one line, with a single space between each entry, bracketed on either side.
      */
      ostream& leftshift(ostream& stream) const;

      // derived
      void print(std::ostream& stream) const;

      // Get pointer to Tramonto double** vector from the NOXLOCA::Tramonto::Vector
      double** get() const;
    private:

      //! Return the i-th element
      double& operator[] (int i);

      //! Return the i-th element
      const double& operator[] (int i) const;

      //! The used length of vector
      int n1; //first dimention
      int n2; //second dimension
      int n;  //total length, n1*n2

      //! The vector owned by this object
      double** x2d; // Allocated  [n1][n2]
      double* x;    //Always  2d[0], the pointer to the vector of length n

      bool needToDelete;
    };

  } // namespace Tramonto
} // namespace NOX

//! Function for printing
ostream& operator<<(ostream& stream, const NOXLOCA::Tramonto::Vector& v);



#endif
