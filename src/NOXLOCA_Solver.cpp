//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA.H"
#include "NOX_Parameter_List.H"
#include "NOXLOCA_Tramonto_Group.hpp"

extern "C" {

void NOXLOCA_Solver(double ** xBox, double **xOwned)
{
  double alpha = 0.0;
  double beta = 0.0;

  try {

    // Create parameter list
    Teuchos::RefCountPtr<NOX::Parameter::List> paramList = 
      Teuchos::rcp(new NOX::Parameter::List);

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.setParameter("Continuation Method", "Natural");
    stepperList.setParameter("Continuation Method", "Arc Length");
    stepperList.setParameter("Continuation Parameter", "alpha");
    stepperList.setParameter("Initial Value", alpha);
    stepperList.setParameter("Max Value", 5.0);
    stepperList.setParameter("Min Value", 0.0);
    stepperList.setParameter("Max Steps", 50);
    stepperList.setParameter("Max Nonlinear Iterations", Max_Newton_iter);
    stepperList.setParameter("Enable Arc Length Scaling", true);
    stepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
    stepperList.setParameter("Initial Scale Factor", 1.0);
    stepperList.setParameter("Min Scale Factor", 1.0e-8);
    stepperList.setParameter("Enable Tangent Factor Step Size Scaling",true);
    stepperList.setParameter("Min Tangent Factor", -1.0);
    stepperList.setParameter("Tangent Factor Exponent",1.0);
    stepperList.setParameter("Compute Eigenvalues",false);
    stepperList.setParameter("Bordered Solver Method", "Bordering");

    // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Type", "None");

    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    predictorList.setParameter("Method", "Tangent");
    //predictorList.setParameter("Method", "Secant");

    NOX::Parameter::List& firstStepPredictorList = 
      predictorList.sublist("First Step Predictor");
    firstStepPredictorList.setParameter("Method", "Constant");

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.setParameter("Method", "Constant");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.1);
    stepSizeList.setParameter("Min Step Size", 1.0e-3);
    stepSizeList.setParameter("Max Step Size", 10.0);
    stepSizeList.setParameter("Aggressiveness", 0.5); // for adaptive
    stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
    stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList->sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");
      NOX::Parameter::List& lineSearchParameters = nlParams.sublist("Line Search");
      lineSearchParameters.setParameter("Method","Backtrack"); //Cuts step for NANs

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			       NOX::Utils::Details +
			       NOX::Utils::OuterIteration + 
			       NOX::Utils::InnerIteration + 
			       NOX::Utils::Warning + 
			       NOX::Utils::StepperIteration +
			       NOX::Utils::StepperDetails);


    // Create global data object
    /*
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);
    */

    // Set up the problem interface
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.

    NOXLOCA::Tramonto::Vector xTV(Nunk_per_node, Nnodes_per_proc, xOwned);

    for (int i=0; i< 111; i++) cout << "xOwned " << xOwned[0][i] << "xTV  " << xTV(i) << endl;

    //Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp = 
    Teuchos::RefCountPtr<NOX::Abstract::Group> grp = 
      Teuchos::rcp(new NOXLOCA::Tramonto::Group(xTV, xBox));
    
    //grp->setParams(p);

    // Set up the status tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> normF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxIters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(Max_Newton_iter));
    Teuchos::RefCountPtr<NOX::StatusTest::Generic> comboOR = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      normF, 
					      maxIters));


    Teuchos::RefCountPtr<NOX::Parameter::List> nlParamsRCP =
	    Teuchos::rcp(&nlParams, false);

    NOX::Solver::Manager solver(grp, comboOR, nlParamsRCP);
    NOX::StatusTest::StatusType status = solver.solve();


    /*
    // Create the stepper  
    LOCA::NewStepper stepper(globalData, grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    // Check for convergence
    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RefCountPtr<const LOCA::LAPACK::Group> finalGroup = 
      Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution = 
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getParameterList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    destroyGlobalData(globalData);
    */
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

}

}
