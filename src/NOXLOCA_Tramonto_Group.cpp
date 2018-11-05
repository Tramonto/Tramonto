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
#include "NOXLOCA_Tramonto_Group.hpp"	// class definition

extern "C" {
#include "loca_const.h"
#define FROM_LOCA 0
#define FROM_MAIN 1
extern void picard_solver(double** xBox, double** xOwned, int iters);
void box2owned(double**, double**);
extern void post_process(double**, int*, double*, int, int,int);
}


NOXLOCA::Tramonto::Group::Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& gD,
                                NOXLOCA::Tramonto::Vector& xVector_, double** xBox_,
                                const LOCA::ParameterVector& paramVec_,
								const Teuchos::RefCountPtr<Teuchos::ParameterList>& paramList_,
                                bool doPicard_):
  LOCA::Abstract::Group(gD),
  xVector(xVector_),	// deep copy      
  fVector(xVector, NOX::ShapeCopy),	// new vector of same size
  newtonVector(xVector, NOX::ShapeCopy),	// new vector of same size
  xBox(xBox_), 		// tmp memory for overlap vector
  paramVec(paramVec_),  // Copy parameter vector
  contStep(0),
  globalData(gD),
  paramList(paramList_),
  secondSolution (false),
  doPicard(doPicard_)
{
  normF = 0;
  this->setParams(paramVec);
}

NOXLOCA::Tramonto::Group::Group(const NOXLOCA::Tramonto::Group& source, NOX::CopyType type) :
  LOCA::Abstract::Group(source.globalData),
  xVector(source.xVector, type), 
  fVector(source.fVector, type),  
  newtonVector(source.newtonVector, type),
  xBox(source.xBox),
  paramVec(source.paramVec),
  contStep(source.contStep),
  globalData(source.globalData),
  paramList(source.paramList),
  secondSolution(false)
{
 
  switch (type) {
    
  case NOX::DeepCopy:
    
    isValidF = source.isValidF;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    normF = source.normF;
    break;

  case NOX::ShapeCopy:
    this->setParams(paramVec);
    normF = 0.0;
    break;

  default:
    std::cerr << "NOX:Tramonto::Group - invalid CopyType for copy constructor." << std::endl;
    throw "NOX Tramonto Error";
  }

}

NOXLOCA::Tramonto::Group::~Group() 
{
}

void NOXLOCA::Tramonto::Group::resetIsValid() //private
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

Teuchos::RefCountPtr<NOX::Abstract::Group> NOXLOCA::Tramonto::Group::
clone(NOX::CopyType type) const 
{
  Teuchos::RefCountPtr<NOX::Abstract::Group> newgrp =
    Teuchos::rcp(new NOXLOCA::Tramonto::Group(*this, type));
  return newgrp;
}

NOX::Abstract::Group& NOXLOCA::Tramonto::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

NOX::Abstract::Group& NOXLOCA::Tramonto::Group::operator=(const Group& source)
{
  if (this != &source) {

    // Copy the xVector
    xVector = source.xVector;

    paramVec = source.paramVec;

    // Update the isValidVectors
    isValidF = source.isValidF;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    
    contStep = source.contStep;

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

void NOXLOCA::Tramonto::Group::setX(const NOX::Abstract::Vector& y) 
{
  setX(dynamic_cast<const NOXLOCA::Tramonto::Vector&> (y));
}

void NOXLOCA::Tramonto::Group::setX(const NOXLOCA::Tramonto::Vector& y) 
{
  resetIsValid();
  xVector = y;
}

void NOXLOCA::Tramonto::Group::computeX(const NOX::Abstract::Group& grp, 
		     const NOX::Abstract::Vector& d, 
		     double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& lapackgrp = dynamic_cast<const Group&> (grp);
  const NOXLOCA::Tramonto::Vector& lapackd = dynamic_cast<const NOXLOCA::Tramonto::Vector&> (d);
  computeX(lapackgrp, lapackd, step); 
}

void NOXLOCA::Tramonto::Group::computeX(const Group& grp, const NOXLOCA::Tramonto::Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
}

NOX::Abstract::Group::ReturnType NOXLOCA::Tramonto::Group::computeF() 
{
  //if (isValidF) 
   // return NOX::Abstract::Group::Ok;
  if (isValidF) {
    return NOX::Abstract::Group::Ok;
  }

  if (!doPicard) {
    fVector.init(0.0); (void) dft_linprobmgr_setrhs(LinProbMgr_manager, fVector.get());

    (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xVector.get(), xBox);
    fill_resid_and_matrix_control_conwrap(xBox, 0, 1);
    (void) dft_linprobmgr_getrhs(LinProbMgr_manager, fVector.get());
    fVector.scale(-1.0);
  }
  else {
    fVector.init(0.0);
    fVector = xVector;
    TV2Box(xVector, xBox);
    picard_solver(xBox, fVector.get(), 1);
    Box2TV(xBox, fVector);
    fVector.update(-1.0, xVector, 1.0); //resid = delta_x
  }

  normF = fVector.norm();

  isValidF = true;
  return (NOX::Abstract::Group::Ok);
}

NOX::Abstract::Group::ReturnType NOXLOCA::Tramonto::Group::computeJacobian() 
{
  // No matrix for Picard solves
  if (doPicard) isValidJacobian=true;

  // Skip if the Jacobian is already valid 
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  (void) dft_linprobmgr_initializeproblemvalues(LinProbMgr_manager);
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xVector.get(), xBox);
  fill_resid_and_matrix_control_conwrap(xBox, 1, 0);
  (void) dft_linprobmgr_finalizeproblemvalues(LinProbMgr_manager);

  //might as well pull out fVector since it was computed anyway, and set as valid
  (void) dft_linprobmgr_getrhs(LinProbMgr_manager, fVector.get());
  fVector.scale(-1.0);
  normF = fVector.norm();
  isValidF = true;

  isValidJacobian = true;
  return (NOX::Abstract::Group::Ok);
}

NOX::Abstract::Group::ReturnType NOXLOCA::Tramonto::Group::computeNewton(Teuchos::ParameterList& p) 
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
NOXLOCA::Tramonto::Group::applyJacobian(const NOX::Abstract::Vector& input, 
				  NOX::Abstract::Vector& result) const
{
  const NOXLOCA::Tramonto::Vector& lapackinput = dynamic_cast<const NOXLOCA::Tramonto::Vector&> (input);
  NOXLOCA::Tramonto::Vector& lapackresult = dynamic_cast<Vector&> (result);
  return applyJacobian(lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::Group::applyJacobian(const NOXLOCA::Tramonto::Vector& input, NOXLOCA::Tramonto::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian() || doPicard) 
    return NOX::Abstract::Group::BadDependency;

  //DO MATVEC HERE
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, input.get(), xBox);
  (void) dft_linprobmgr_applymatrix(LinProbMgr_manager, xBox, result.get());


  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::Group::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const 
{
  const NOXLOCA::Tramonto::Vector& lapackinput = dynamic_cast<const NOXLOCA::Tramonto::Vector&> (input);
  NOXLOCA::Tramonto::Vector& lapackresult = dynamic_cast<Vector&> (result); 
  return applyJacobianInverse(p, lapackinput, lapackresult);
}

NOX::Abstract::Group::ReturnType 
NOXLOCA::Tramonto::Group::applyJacobianInverse(Teuchos::ParameterList& p, 
					 const NOXLOCA::Tramonto::Vector& input, 
					 NOXLOCA::Tramonto::Vector& result) const 
{

  if (!isJacobian()) {
    std::cerr << "ERROR: NOXLOCA::Tramonto::Group::applyJacobianInverse() - invalid Jacobian" << std::endl;
    throw "NOX Error";
  }

  if (!doPicard) {
    (void) dft_linprobmgr_setrhs(LinProbMgr_manager, input.get());
    (void) dft_linprobmgr_setupsolver(LinProbMgr_manager);
    (void) dft_linprobmgr_solve(LinProbMgr_manager);
    (void) dft_linprobmgr_getlhs(LinProbMgr_manager, xBox);
    (void) box2owned(xBox, result.get());
  }
  else {
    result.init(0.0);
    TV2Box(xVector, xBox);
    result = xVector;
    picard_solver(xBox, result.get(), 100); // 100 Picard steps per "Newton" iter
    Box2TV(xBox, result);
    result.update(1.0, xVector, -1.0); //result = -delta_x
  }

  return NOX::Abstract::Group::Ok;
}

bool NOXLOCA::Tramonto::Group::isF() const 
{   
  return isValidF;
}

bool NOXLOCA::Tramonto::Group::isJacobian() const 
{  
  return isValidJacobian;
}

bool NOXLOCA::Tramonto::Group::isNewton() const 
{   
  return isValidNewton;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::Group::getX() const 
{
  return xVector;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::Group::getF() const 
{  
  return fVector;
}

double NOXLOCA::Tramonto::Group::getNormF() const
{
  return normF;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::Group::getNewton() const 
{
  return newtonVector;
}

const NOX::Abstract::Vector& NOXLOCA::Tramonto::Group::getGradient() const 
{
  std::cout << "ERROR: GRADIENT VECTOR NOT CALCULATEED IN TRAMONTO_GROUP!! " << std::endl;
  return newtonVector;
}


Teuchos::RCP< const NOX::Abstract::Vector >
NOXLOCA::Tramonto::Group::getXPtr() const
{
  std::cout << "ERROR: getXPtr not debugged for NOXLOCA::Tramonto::Group!!!!" << std::endl;
  throw "NOX Error";
  return (*solversVecPtr)[0]->getSolutionGroup().getXPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOXLOCA::Tramonto::Group::getFPtr() const
{
  std::cout << "ERROR: getFPtr not debugged for NOXLOCA::Tramonto::Group!!!!" << std::endl;
  throw "NOX Error";
  return (*solversVecPtr)[0]->getSolutionGroup().getFPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOXLOCA::Tramonto::Group::getGradientPtr() const
{
  std::cout << "ERROR: getGradientPtr not debugged for NOXLOCA::Tramonto::Group!!!!" << std::endl;
  throw "NOX Error";
  return (*solversVecPtr)[0]->getSolutionGroup().getGradientPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOXLOCA::Tramonto::Group::getNewtonPtr() const
{
  std::cout << "ERROR: getNewtonPtr not debugged for NOXLOCA::Tramonto::Group!!!!" << std::endl;
  throw "NOX Error";
  return (*solversVecPtr)[0]->getSolutionGroup().getNewtonPtr();
}

void NOXLOCA::Tramonto::Group::print() const
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

void NOXLOCA::Tramonto::Group::copy(const NOX::Abstract::Group& source)
{ *this = source; }

void  NOXLOCA::Tramonto::Group::setParams(const LOCA::ParameterVector& p)
{ for (int i=0; i<<paramVec.length(); i++) setParam(i, paramVec.getValue(i));}

void  NOXLOCA::Tramonto::Group::setParam(std::string paramID, double val)
{ 
  resetIsValid();
  paramVec.setValue(paramID, val);

  if (paramID=="ConParam") assign_parameter_conwrap(val);
  else if (paramID=="BifParam") assign_bif_parameter_conwrap(val);
}

void  NOXLOCA::Tramonto::Group::setParam(int paramID, double val)
{ setParam(paramVec.getLabel(paramID), val); }

const LOCA::ParameterVector&  NOXLOCA::Tramonto::Group::getParams() const
{ return paramVec; }
double  NOXLOCA::Tramonto::Group::getParam(std::string paramID) const
{  return paramVec.getValue(paramID); }
double  NOXLOCA::Tramonto::Group::getParam(int paramID) const
{  return paramVec.getValue(paramID); }


void  NOXLOCA::Tramonto::Group::printSolution(const NOX::Abstract::Vector& sol_,
      const double param) const
{ 
  double time_save=0.0;
  int num_its=paramList->sublist("NOX").sublist("Output").get("Nonlinear Iterations",-1);

  const NOXLOCA::Tramonto::Vector& solVector =
      dynamic_cast<const NOXLOCA::Tramonto::Vector&>(sol_);
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, solVector.get(), xBox);

//cout<<"calling post_process from NOXLOCA 1 contStep="<<contStep<<"secondSolution="<<secondSolution<<"from_flag="<<FROM_LOCA<<endl;
  post_process(xBox, &num_its, &time_save, contStep, secondSolution, FROM_LOCA);
  contStep++;
}

void  NOXLOCA::Tramonto::Group::printSolution2(const NOX::Abstract::Vector& sol_,
      const double param) const
{ 
  secondSolution = true; // this is the only place this bool is set to true
  contStep--; // cancel out incrementing of this counter for second solution
  printSolution(sol_, param);
  secondSolution = false;
}

void  NOXLOCA::Tramonto::Group::printSolution(const double param) const
{
  double time_save=0.0;
  int num_its=paramList->sublist("NOX").sublist("Output").get("Nonlinear Iterations",-1);

  contStep++;
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xVector.get(), xBox);
  
//cout<<"calling post_process from NOXLOCA 1 contStep="<<contStep<<"secondSolution="<<FALSE<<"from_flag="<<FROM_LOCA<<endl;
  post_process(xBox, &num_its, &time_save, contStep, FALSE, FROM_LOCA);

}

double  NOXLOCA::Tramonto::Group::calcFreeEnergy() const
{
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xVector.get(), xBox);
  return calc_free_energy_conwrap(xBox);
}

void  NOXLOCA::Tramonto::Group::TV2Box(const NOXLOCA::Tramonto::Vector& xTV, double** xB) const
{
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xTV.get(), xB);
}

void  NOXLOCA::Tramonto::Group::Box2TV(double** xB, NOXLOCA::Tramonto::Vector& xTV) const
{
  for (int iunk=0;iunk<Nunk_per_node;iunk++){
    for (int loc_inode=0;loc_inode<Nnodes_per_proc;loc_inode++)
      xTV.get()[iunk][loc_inode] = xB[iunk][L2B_node[loc_inode]];
  }
}

