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

#include "LOCA.H"
#include "NOXLOCA_Tramonto_Group.hpp"

extern "C" {

void NOXLOCA_Solver(double ** xBox, double **xOwned)
{
  double alpha = 0.0;
  double beta = 0.0;

  try {

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Method", "Arc Length");
    stepperList.set("Continuation Parameter", "alpha");
    stepperList.set("Initial Value", alpha);
    stepperList.set("Max Value", 5.0);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 50);
    stepperList.set("Max Nonlinear Iterations", Max_Newton_iter);
    stepperList.set("Enable Arc Length Scaling", true);
    stepperList.set("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.set("Max Arc Length Parameter Contribution", 0.7);
    stepperList.set("Initial Scale Factor", 1.0);
    stepperList.set("Min Scale Factor", 1.0e-8);
    stepperList.set("Enable Tangent Factor Step Size Scaling",true);
    stepperList.set("Min Tangent Factor", -1.0);
    stepperList.set("Tangent Factor Exponent",1.0);
    stepperList.set("Compute Eigenvalues",false);
    stepperList.set("Bordered Solver Method", "Bordering");

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.set("Method", "Constant");
    predictorList.set("Method", "Tangent");
    //predictorList.set("Method", "Secant");

    Teuchos::ParameterList& firstStepPredictorList = 
      predictorList.sublist("First Step Predictor");
    firstStepPredictorList.set("Method", "Constant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.set("Method", "Constant");
    stepSizeList.set("Method", "Adaptive");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 10.0);
    stepSizeList.set("Aggressiveness", 0.5); // for adaptive
    stepSizeList.set("Failed Step Reduction Factor", 0.5);
    stepSizeList.set("Successful Step Increase Factor", 1.26); // for constant

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    nlParams.set("Nonlinear Solver", "Line Search Based");
      Teuchos::ParameterList& lineSearchParameters = nlParams.sublist("Line Search");
      lineSearchParameters.set("Method","Backtrack"); //Cuts step for NANs

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
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


    Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsRCP =
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
