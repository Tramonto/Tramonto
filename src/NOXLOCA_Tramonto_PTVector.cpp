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
#include "NOXLOCA_Tramonto_PTVector.hpp"
#include "NOX_Random.H" // for Random class

//never called
NOXLOCA::Tramonto::PTVector::PTVector() :
x1(), x2(), ptp(0.0)
{
}

//
NOXLOCA::Tramonto::PTVector::PTVector(const NOXLOCA::Tramonto::Vector& x1_,
                   const NOXLOCA::Tramonto::Vector& x2_, const double& ptp_) : 
  x1(x1_),
  x2(x2_),
  ptp(ptp_)
{
}

NOXLOCA::Tramonto::PTVector::PTVector(const NOXLOCA::Tramonto::PTVector& source, 
			    NOX::CopyType type) :
  x1(source.x1),
  x2(source.x2),
  ptp(source.ptp)
{
}

NOXLOCA::Tramonto::PTVector::~PTVector()
{
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::operator=(
					   const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(source));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::operator=(
					   const NOXLOCA::Tramonto::PTVector& source)
{
  x1 = source.x1;
  x2 = source.x2;
  ptp= source.ptp;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::init(double value)
{
  x1.init(value);
  x2.init(value);
  ptp = value;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::random(bool useSeed, int seed) 
{
  x1.random(useSeed, seed);
  x2.random(useSeed, seed);
  ptp = 0.5;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::abs(
					     const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(base));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::abs(
					     const NOXLOCA::Tramonto::PTVector& base)
{
  x1.abs(base.x1);
  x2.abs(base.x2);
  ptp = fabs(base.ptp);
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::reciprocal(
					    const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(base));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::reciprocal(
					    const NOXLOCA::Tramonto::PTVector& base)
{
  x1.reciprocal(base.x1);
  x2.reciprocal(base.x2);
  ptp = 1.0 / base.ptp;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::scale(double alpha)
{
  x1.scale(alpha);
  x2.scale(alpha);
  ptp *= alpha;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::update(
			       double alpha, 
			       const NOX::Abstract::Vector& a, 
			       double gamma)
{
  return update(alpha, dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(a), gamma);
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::update(
				 double alpha, 
				 const NOXLOCA::Tramonto::PTVector& a, 
				 double gamma)
{
  x1.update(alpha, a.x1, gamma);
  x2.update(alpha, a.x2, gamma);
  ptp = alpha * a.ptp + gamma * ptp;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::update(
			      double alpha, 
			      const NOX::Abstract::Vector& a, 
			      double beta, 
			      const NOX::Abstract::Vector& b,
			      double gamma)
{
  return update(alpha, dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(a), 
		beta, dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(b), gamma);
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::update(
			       double alpha, 
			       const NOXLOCA::Tramonto::PTVector& a, 
			       double beta, 
			       const NOXLOCA::Tramonto::PTVector& b,
			       double gamma)
{

  x1.update(alpha, a.x1, beta, b.x1, gamma);
  x2.update(alpha, a.x2, beta, b.x2, gamma);
  ptp = alpha * a.ptp + + beta * b.ptp + gamma * ptp;
  return *this;
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::scale(
					      const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(a));
}

NOX::Abstract::Vector& NOXLOCA::Tramonto::PTVector::scale(const NOXLOCA::Tramonto::PTVector& a)
{  
  x1.scale(a.x1);
  x2.scale(a.x2);
  ptp *= a.ptp;
  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::Vector> NOXLOCA::Tramonto::PTVector::
clone(NOX::CopyType type) const
{
  Teuchos::RefCountPtr<NOX::Abstract::Vector> tmp;
  tmp = Teuchos::rcp(new NOXLOCA::Tramonto::PTVector(*this, type));
  return tmp;
}

double NOXLOCA::Tramonto::PTVector::norm(NOX::Abstract::Vector::NormType type) const
{
  
  double value=0.0;			// final answer
  double value1 = x1.norm(type);
  double value2 = x2.norm(type);
  double valuep = fabs(ptp);

  switch (type) {
  case MaxNorm:
      value = value1;
      if (value2 > value ) value=value2;
      if (valuep > value ) value=valuep;
    break;
  case OneNorm:
    value = value1 + value2 + valuep;
    break;
  case TwoNorm:
    value = value1*value1 + value2*value2 + valuep*valuep;
    value = sqrt(value);
  default:
   break;
  }

  return value;
}

double NOXLOCA::Tramonto::PTVector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(weights));
}

double NOXLOCA::Tramonto::PTVector::norm(const NOXLOCA::Tramonto::PTVector& weights) const
{
  double value=0.0;			// final answer
  double value1 = x1.norm(weights.X1());
  double value2 = x2.norm(weights.X2());
  double valuep_sq = weights.PTP() * ptp * ptp; 

  value = value1*value1 + value2*value2 + valuep_sq;
  value = sqrt(value);

  return value;
}

double NOXLOCA::Tramonto::PTVector::innerProduct(const NOX::Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(y));
}

double NOXLOCA::Tramonto::PTVector::innerProduct(const NOXLOCA::Tramonto::PTVector& y) const
{
  double value=0.0;			// final answer
  double value1 = x1.innerProduct(y.X1());
  double value2 = x2.innerProduct(y.X2());
  double valuep = ptp * y.PTP();

  value = value1 + value2 + valuep;
  return value;
}

int NOXLOCA::Tramonto::PTVector::length() const
{
  return x1.length() + x2.length() + 1;;
}

NOXLOCA::Tramonto::Vector& NOXLOCA::Tramonto::PTVector::X1()
{ return x1; }
const NOXLOCA::Tramonto::Vector& NOXLOCA::Tramonto::PTVector::X1() const
{ return x1; }
NOXLOCA::Tramonto::Vector& NOXLOCA::Tramonto::PTVector::X2()
{ return x2; }
const NOXLOCA::Tramonto::Vector& NOXLOCA::Tramonto::PTVector::X2() const
{ return x2; }
double& NOXLOCA::Tramonto::PTVector::PTP()
{ return ptp; }
const double& NOXLOCA::Tramonto::PTVector::PTP() const
{ return ptp; }

ostream& NOXLOCA::Tramonto::PTVector::leftshift(ostream& stream) const
{
  stream << "{ ";
  stream <<  x1.leftshift(stream) << " " << x2.leftshift(stream)
         <<  " [ " << ptp << "] ";
  stream << "}";
  return stream;
}

ostream& operator<<(ostream& stream, const NOXLOCA::Tramonto::PTVector& v)
{
  return v.leftshift(stream);
}

void NOXLOCA::Tramonto::PTVector::print(std::ostream& stream) const
{
  stream << *this << endl;
}

