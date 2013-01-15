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
#include "NOXLOCA_Tramonto_PTGroup.hpp"	// class definition

extern "C" {
extern void post_process(double**, char*, int*, double*, int, int);
}

NOXLOCA::Tramonto::PTGroup::PTGroup(
      const Teuchos::RefCountPtr<LOCA::GlobalData> gD,
      const Teuchos::RefCountPtr<NOXLOCA::Tramonto::Group>& grp_,
      const NOXLOCA::Tramonto::PTVector& xVector_):
  LOCA::Abstract::Group(gD),
  trGrp(grp_), // Underlying group for regular system of size n
  xVector(xVector_),	// deep copy      
  fVector(xVector_, NOX::ShapeCopy),	// new vector of same size
  newtonVector(xVector_, NOX::ShapeCopy),	// new vector of same size
  globalData(gD),
  normF(0)
{
}

NOXLOCA::Tramonto::PTGroup::PTGroup(const NOXLOCA::Tramonto::PTGroup& source, NOX::CopyType type) :
  LOCA::Abstract::Group(source.globalData),
  trGrp(source.trGrp),
  xVector(source.xVector, type), 
  fVector(source.fVector, type),  
  newtonVector(source.newtonVector, type),
  globalData(source.globalData)
{
 
  switch (type) {
    
  case NOX::DeepCopy:
    
    isValidF = source.isValidF;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    normF = source.normF;
    break;

  case NOX::ShapeCopy:
    normF = 0.0;
    break;

  default:
    std::cerr << "NOXLOCA:Tramonto::PTGroup - invalid CopyType for copy constructor." << std::endl;
    throw "NOXLOCA Tramonto Error";
  }

}

NOXLOCA::Tramonto::PTGroup::~PTGroup() 
{
}

void NOXLOCA::Tramonto::PTGroup::resetIsValid() //private
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

Teuchos::RefCountPtr<NOX::Abstract::Group> NOXLOCA::Tramonto::PTGroup::
clone(NOX::CopyType type) const 
{
  Teuchos::RefCountPtr<NOX::Abstract::Group> newgrp = 
    Teuchos::rcp(new NOXLOCA::Tramonto::PTGroup(*this, type));
  return newgrp;
}

NOX::Abstract::Group& NOXLOCA::Tramonto::PTGroup::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const NOXLOCA::Tramonto::PTGroup&> (source));
}

NOX::Abstract::Group& NOXLOCA::Tramonto::PTGroup::operator=(const PTGroup& source)
{
  if (this != &source) {

    trGrp = source.trGrp;

    // Copy the xVector
    xVector = source.xVector;

    // Update the isValidVectors
    isValidF = source.isValidF;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    
    // Only copy vectors that are valid
    if (isValidF) {
      fVector = source.fVector;
      normF = source.normF;
    }

    if (isValidNewton)
      newtonVector = source.newtonVector;
    
  }

  return *this;
}

void NOXLOCA::Tramonto::PTGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX(dynamic_cast<const NOXLOCA::Tramonto::PTVector&> (y));
}

void NOXLOCA::Tramonto::PTGroup::setX(const NOXLOCA::Tramonto::PTVector& y) 
{
  resetIsValid();
  xVector = y;
}

void NOXLOCA::Tramonto::PTGroup::computeX(const NOX::Abstract::Group& grp, 
		     const NOX::Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const PTGroup& trgrp = dynamic_cast<const PTGroup&> (grp);
  const NOXLOCA::Tramonto::PTVector& trd = dynamic_cast<const NOXLOCA::Tramonto::PTVector&> (d);
  computeX(trgrp, trd, step); 
}

void NOXLOCA::Tramonto::PTGroup::computeX(const PTGroup& grp, const NOXLOCA::Tramonto::PTVector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
}

NOX::Abstract::Group::ReturnType NOXLOCA::Tramonto::PTGroup::computeF() 
{
  if (isValidF) 
    return NOX::Abstract::Group::Ok;

  trGrp->setParam("BifParam", xVector.PTP());

  trGrp->setX(xVector.X1());
  trGrp->computeF();
  fVector.X1() = trGrp->getF();
  double omega1 = trGrp->calcFreeEnergy();

  trGrp->setX(xVector.X2());
  trGrp->computeF();
  fVector.X2() = trGrp->getF();
  double omega2 = trGrp->calcFreeEnergy();

  fVector.PTP() = omega1 - omega2;

  normF = fVector.norm();

  isValidF = true;
  return (NOX::Abstract::Group::Ok);
}

NOX::Abstract::Group::ReturnType NOXLOCA::Tramonto::PTGroup::computeJacobian() 
{
  // To save memory, only compute Jacobians at withing applyJacobianInverse. 
  return (NOX::Abstract::Group::Ok);
}

NOX::Abstract::Group::ReturnType NOXLOCA::Tramonto::PTGroup::computeNewton(Teuchos::ParameterList& p) 
{
  if (isNewton())
    return NOX::Abstract::Group::Ok;

  if (!isF()) {
    std::cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid F" << std::endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    std::cerr << "ERROR: NOX::Example::Group::computeNewton() - invalid Jacobian" << std::endl;
    throw "NOX Error";
  }

  NOX::Abstract::Group::ReturnType status = applyJacobianInverse(p, fVector, newtonVector);
  isValidNewton = (status == NOX::Abstract::Group::Ok);

  // Scale soln by -1
  newtonVector.scale(-1.0);

  // Return solution
  return status;
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::PTGroup::applyJacobian(const NOX::Abstract::Vector& input, 
				  NOX::Abstract::Vector& result) const
{
  const NOXLOCA::Tramonto::PTVector& lapackinput = dynamic_cast<const NOXLOCA::Tramonto::PTVector&> (input);
  NOXLOCA::Tramonto::PTVector& lapackresult = dynamic_cast<PTVector&> (result);
  return applyJacobian(lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::PTGroup::applyJacobian(const NOXLOCA::Tramonto::PTVector& input,
                                          NOXLOCA::Tramonto::PTVector& result) const
{
  std::cout << "ERROR:  Apply Jacobian not implemented for PTGroup !!!!" << std::endl;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::PTGroup::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  const NOXLOCA::Tramonto::PTVector& lapackinput = dynamic_cast<const NOXLOCA::Tramonto::PTVector&> (input);
  NOXLOCA::Tramonto::PTVector& lapackresult = dynamic_cast<PTVector&> (result); 
  return applyJacobianInverse(p, lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::PTGroup::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const NOXLOCA::Tramonto::PTVector& input, 
					 NOXLOCA::Tramonto::PTVector& result) const 
{
  // Algorithm from Equations 18-26 of Salinger&Frink JChemPhys (2003).
  // This implementation does not assume that input=f so it can solve
  // other RHS's for arclength or bifurcations of PhaseTransitions
  
  const double eps=1.0e-7;
  double perturb = eps * (eps + fabs(xVector.PTP()));
  double pertParam = xVector.PTP() + perturb;
  perturb = pertParam - xVector.PTP();  //improves roundoff accuracy

  // temporary space...
  NOXLOCA::Tramonto::PTVector bdvec(result);
  NOXLOCA::Tramonto::Vector fVec(result.X1());
  NOXLOCA::Tramonto::Vector fPertVec(result.X1());
 
  // First matrix block...

  trGrp->setX(xVector.X1());

  //perturb parameter
  trGrp->setParam("BifParam", pertParam);
  trGrp->computeF();
  fPertVec = trGrp->getF();
  double omega1pert = trGrp->calcFreeEnergy();

  //unperturb parameter, compute df/dp and load Jacobian
  trGrp->setParam("BifParam", xVector.PTP());
  trGrp->computeJacobian(); // computes resid as well
  fVec = trGrp->getF();
  fPertVec.update(-1.0, fVec, 1.0);
  fPertVec.scale(1.0/perturb);
  double omega1 = trGrp->calcFreeEnergy();

  // Do two solves with same matrix for "a" and "b" vectors, Equations 18,19
  trGrp->applyJacobianInverse(p, input.X1(), result.X1());
  trGrp->applyJacobianInverse(p, fPertVec, bdvec.X1());

  // Second matrix block... same but with X2

  trGrp->setX(xVector.X2());

  //perturb parameter
  trGrp->setParam("BifParam", pertParam);
  trGrp->computeF();
  fPertVec = trGrp->getF();
  double omega2pert = trGrp->calcFreeEnergy();

  //unperturb parameter, compute df/dp and load Jacobian
  trGrp->setParam("BifParam", xVector.PTP());
  trGrp->computeJacobian(); // computes resid as well
  fVec = trGrp->getF();
  fPertVec.update(-1.0, fVec, 1.0);
  fPertVec.scale(1.0/perturb);
  double omega2 = trGrp->calcFreeEnergy();

  // Do two solves with same matrix for "c" and "d" vectors, Equations 20,21
  trGrp->applyJacobianInverse(p, input.X2(), result.X2());
  trGrp->applyJacobianInverse(p, fPertVec, bdvec.X2());

  // Compute contributions of equal-energy constraint equation
  
  double g = input.PTP(); // may not be omega1 - omega2 !
  double dgdp = (omega1pert - omega1 - omega2pert + omega2) / perturb;
  
  //For each of 4 terms: calc perturbation, calc perturbed energy, diff
  //Equation 22
  perturb = eps * xVector.X1().norm() / (eps + result.X1().norm());
  fPertVec.update(1.0, xVector.X1(), perturb, result.X1(), 0.0);
  trGrp->setX(fPertVec);
  double dOmdx1a = (trGrp->calcFreeEnergy() - omega1) / perturb;
  
  perturb = eps * xVector.X2().norm() / (eps + result.X2().norm());
  fPertVec.update(1.0, xVector.X2(), perturb, result.X2(), 0.0);
  trGrp->setX(fPertVec);
  double dOmdx2c = (trGrp->calcFreeEnergy() - omega2) / perturb;
  
  perturb = eps * xVector.X1().norm() / (eps + bdvec.X1().norm());
  fPertVec.update(1.0, xVector.X1(), perturb, bdvec.X1(), 0.0);
  trGrp->setX(fPertVec);
  double dOmdx1b = (trGrp->calcFreeEnergy() - omega1) / perturb;
  
  perturb = eps * xVector.X2().norm() / (eps + bdvec.X2().norm());
  fPertVec.update(1.0, xVector.X2(), perturb, bdvec.X2(), 0.0);
  trGrp->setX(fPertVec);
  double dOmdx2d = (trGrp->calcFreeEnergy() - omega2) / perturb;

  // Equations 23 and 24
  double delta_param = ( g - dOmdx1a + dOmdx2c) / (dgdp + dOmdx1b - dOmdx2d);
  bdvec.PTP() = 0.0;
  result.update(delta_param, bdvec, 1.0);
  result.PTP() = delta_param;

  return NOX::Abstract::Group::Ok;
}

bool NOXLOCA::Tramonto::PTGroup::isF() const 
{   
  return isValidF;
}

bool NOXLOCA::Tramonto::PTGroup::isJacobian() const 
{  
  return isValidJacobian;
}

bool NOXLOCA::Tramonto::PTGroup::isNewton() const 
{   
  return isValidNewton;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::PTGroup::getX() const 
{
  return xVector;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::PTGroup::getF() const 
{  
  return fVector;
}

double NOXLOCA::Tramonto::PTGroup::getNormF() const
{
  return normF;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::PTGroup::getNewton() const 
{
  return newtonVector;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::PTGroup::getGradient() const 
{
  std::cout << "ERROR: GRADIENT VECTOR NOT CALCULATEED IN TRAMONTO_GROUP!! " << std::endl;
  return newtonVector;
}


void NOXLOCA::Tramonto::PTGroup::print() const
{
  std::cout << "x = " << xVector << "\n";

  if (isValidF) {
    std::cout << "F(x) = " << fVector << "\n";
    std::cout << "|| F(x) || = " << normF << "\n";
  }
  else
    std::cout << "F(x) has not been computed" << "\n";
  
  std::cout << std::endl;
}

void  NOXLOCA::Tramonto::PTGroup::setParams(const LOCA::ParameterVector& p)
{ trGrp->setParams(p);}

void  NOXLOCA::Tramonto::PTGroup::setParam(std::string paramID, double val)
{ 
  resetIsValid();
  trGrp->setParam(paramID, val);
}

void  NOXLOCA::Tramonto::PTGroup::setParam(int paramID, double val)
{ trGrp->setParam(paramID, val); }

const LOCA::ParameterVector&  NOXLOCA::Tramonto::PTGroup::getParams() const
{ return trGrp->getParams(); }
double  NOXLOCA::Tramonto::PTGroup::getParam(std::string paramID) const
{  return trGrp->getParam(paramID); }
double  NOXLOCA::Tramonto::PTGroup::getParam(int paramID) const
{  return trGrp->getParam(paramID); }


void  NOXLOCA::Tramonto::PTGroup::printSolution(const NOX::Abstract::Vector& sol_,
      const double param) const
{ 

  const NOXLOCA::Tramonto::Vector& sol =
    (dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(sol_)).X1();
  trGrp->printSolution(sol, param);

  const NOXLOCA::Tramonto::Vector& sol2 =
    (dynamic_cast<const NOXLOCA::Tramonto::PTVector&>(sol_)).X2();
  trGrp->printSolution2(sol2, param);
}

void  NOXLOCA::Tramonto::PTGroup::printSolution(const double param) const
{
  printSolution(xVector, param);
}
