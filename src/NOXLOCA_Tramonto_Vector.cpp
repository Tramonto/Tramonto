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

#include "NOX_Common.H"
#include "NOXLOCA_Tramonto_Vector.hpp"
#include "NOX_Random.H" // for Random class

//never called
NOXLOCA::Tramonto::Vector::Vector() :
  n1(0),
  n2(0),
  n(n1*n2),
  x2d(0),
  x(0),
  needToDelete(false)
{
}

//view of double** vector
NOXLOCA::Tramonto::Vector::Vector(int N1, int N2, double **v) : 
  n1(N1),
  n2(N2),
  n(N1*N2),
  x2d(v),
  x(v[0]),
  needToDelete(false)
{
}

NOXLOCA::Tramonto::Vector::Vector(const NOXLOCA::Tramonto::Vector& source, 
			    NOX::CopyType type) :
  n1(source.n1),
  n2(source.n2),
  n(source.n),
  x2d(0),
  x(0),
  needToDelete(true)
{
   x2d = (double **) array_alloc_2d_conwrap(n1, n2, sizeof(double));
   x = x2d[0];
   for (int i=0; i < n; i++) x[i] = source.x[i];
}

NOXLOCA::Tramonto::Vector::~Vector()
{
   if (needToDelete) safe_free_conwrap((void**)&x2d);
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::operator=(
					   const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOXLOCA::Tramonto::Vector&>(source));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::operator=(
					   const NOXLOCA::Tramonto::Vector& source)
{
  for (int i=0; i<n; i++) x[i] = source.x[i];
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::init(double value)
{
  for (int i = 0; i < n; i ++)
    x[i] = value;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::random(bool useSeed, int seed) 
{
  if (useSeed)
    NOX::Random::setSeed(seed);

  for (int i = 0; i < n; i ++) 
    x[i] = NOX::Random::number();

  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::abs(
					     const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOXLOCA::Tramonto::Vector&>(base));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::abs(
					     const NOXLOCA::Tramonto::Vector& base)
{
  for (int i = 0; i < n; i ++)
    x[i] = fabs(base[i]);
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::reciprocal(
					    const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOXLOCA::Tramonto::Vector&>(base));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::reciprocal(
					    const NOXLOCA::Tramonto::Vector& base)
{
  for (int i = 0; i < n; i ++)
    x[i] = 1.0 / base[i];
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::scale(double alpha)
{
  for (int i = 0; i <n; i ++)
    x[i] *= alpha;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::update(
					       double alpha, 
					       const NOX::Abstract::Vector& a, 
					       double gamma)
{
  return update(alpha, dynamic_cast<const NOXLOCA::Tramonto::Vector&>(a), gamma);
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::update(
						 double alpha, 
						 const NOXLOCA::Tramonto::Vector& a, 
						 double gamma)
{
  for (int i = 0; i < n; i ++)
    x[i] = alpha * a[i] + gamma * x[i];
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::update(
					      double alpha, 
					      const NOX::Abstract::Vector& a, 
					      double beta, 
					      const NOX::Abstract::Vector& b,
					      double gamma)
{
  return update(alpha, dynamic_cast<const NOXLOCA::Tramonto::Vector&>(a), 
		beta, dynamic_cast<const NOXLOCA::Tramonto::Vector&>(b), gamma);
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::update(
					       double alpha, 
					       const NOXLOCA::Tramonto::Vector& a, 
					       double beta, 
					       const NOXLOCA::Tramonto::Vector& b,
					       double gamma)
{
  for (int i = 0; i < n; i ++)
    x[i] = alpha * a[i] + beta * b[i] + gamma * x[i];
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::scale(
					      const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const NOXLOCA::Tramonto::Vector&>(a));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::Vector::scale(const NOXLOCA::Tramonto::Vector& a)
{  
  for (int i = 0; i < n; i ++)
    x[i] = a[i] * x[i];
  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::Vector> NOXLOCA::Tramonto::Vector::
clone(NOX::CopyType type) const
{
  Teuchos::RefCountPtr<NOX::Abstract::Vector> tmp;
  tmp = Teuchos::rcp(new NOXLOCA::Tramonto::Vector(*this, type));
  return tmp;
}

double NOXLOCA::Tramonto::Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  
  if (n == 0)
    return 0.0;

  int i;			// counter
  double value=0.0;			// final answer

  switch (type) {
  case MaxNorm:
    value = fabs(x[0]);
    for (i = 1; i < n; i ++)
      if (value < fabs(x[i]))
	value = fabs(x[i]);
    value = gmax_double_conwrap(value);
    break;
  case OneNorm:
    for (i = 0; i < n; i ++)
	value += fabs(x[i]);
    value = gsum_double_conwrap(value);
    break;
  case TwoNorm:
    for (i = 0; i < n; i ++)
	value += x[i]*x[i];
    value = gsum_double_conwrap(value);
    value = sqrt(value);
  default:
   break;
  }

  return value;
}

double NOXLOCA::Tramonto::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOXLOCA::Tramonto::Vector&>(weights));
}

double NOXLOCA::Tramonto::Vector::norm(const NOXLOCA::Tramonto::Vector& weights) const
{
  if (weights.length() != n) {
    cerr << "NOXLOCA::Tramonto::Vector::norm - size mismatch for weights vector" << endl;
    throw "NOXLOCA::Tramonto Error";
  }

  double value = 0;		// final answer

  for (int i = 0; i < n; i ++)
    value += weights[i] * x[i] * x[i];

  value = gsum_double_conwrap(value);
  value = sqrt(value);

  return value;
}

double NOXLOCA::Tramonto::Vector::innerProduct(const NOX::Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const NOXLOCA::Tramonto::Vector&>(y));
}

double NOXLOCA::Tramonto::Vector::innerProduct(const NOXLOCA::Tramonto::Vector& y) const
{
  if (y.length() != n) {
    cerr << "NOXLOCA::Tramonto::Vector::innerProduct - size mismatch for y vector" 
	 << endl;
    throw "NOX::Tramonto Error";
  }

  double value = 0;		// final answer

  for (int i = 0; i < n; i ++)
    value +=  x[i] * y[i];

  value = gsum_double_conwrap(value);
  return value;
}

int NOXLOCA::Tramonto::Vector::length() const
{
  return n;
}

double& NOXLOCA::Tramonto::Vector::operator[] (int i)
{
  return x[i];
}

const double& NOXLOCA::Tramonto::Vector::operator[] (int i) const
{
  return x[i];
}

double& NOXLOCA::Tramonto::Vector::operator() (int i)
{
  return x[i];
}

const double& NOXLOCA::Tramonto::Vector::operator() (int i) const
{
  return x[i];
}

ostream& NOXLOCA::Tramonto::Vector::leftshift(ostream& stream) const
{
  stream << "[ ";
  for (int i = 0; i < n; i ++)
    stream << x[i] << " ";
  stream << "]";
  return stream;
}

ostream& operator<<(ostream& stream, const NOXLOCA::Tramonto::Vector& v)
{
  return v.leftshift(stream);
}

void NOXLOCA::Tramonto::Vector::print(std::ostream& stream) const
{
  stream << *this << endl;
}

double** NOXLOCA::Tramonto::Vector::get() const
{
  return x2d;
}

