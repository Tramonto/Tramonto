package Tramonto;
import XML.*;
import java.util.*;
import java.io.*;
import javax.swing.*;


// This class takes an XML message containing input parameters
// to the Tramonto simulation package and can do two things
// First, it validate that the parameters passed in the XML 
// message are acceptable input to Tramonto.  
// Second, it can create an input file on the local file system
// which can be read by Tramonto to start a simulation.

public class TramontoXMLParser
{
    // the xml message containing the input parameters
    private XMLObject inputMessage;
    private XMLObject problemTypeSection;
    private XMLObject functionalsAndInteractionsSection;

    // these strings are used to parse the XML message
    // depending on the dimensionality of the problem
    String problemDimensionPrefix;
    String problemDimensionSuffix;
    String polymericSuffix;
    int problemDimension;

    // data items for the simulation.  These get filled as
    // the xml inputMessage is validated by validateXML 
    // and it's sub routines.
    String inputFileName;
    String workingDirectory;
    boolean realSpaceComputation;               // flag for real/fourier space computation
    // Dimension parameters
    double referenceLength;                     // these reference quantities are set 
    double referenceDensity;                    // in the state parameters section, but
    double referenceTemp;                       // are placed here so that they will be
    double referenceDielectric;                 // writen as one of the first sections to the
    double maximumPotential;                    // to re-read the input file once there are read.
    // Mesh parameters
    int nDim;                                   // number of dimensions in problem (x,y,z)
    double size_x[];                            // size of computational domain
    double esize_x[];                           // mesh spacing in domanan
    int type_bc_lbb[];                          // boundary condition type [left, bottom back]
    int type_bc_rtf[];                          // boundary condition type [right, top, front]
    // Functional switch parameters
    int type_func;                              // type of hard sphere functional 
    int type_attr;                              // type of attractive functional 
    int type_coul;                              // type of functional for coulombic systems
    int type_poly;                              // type polymer functional
    int compareToFastram;                       // flag to do calculations in a way consistent with Fastram
    // Potential type parameters
    int ipot_wf_n;                              // wall-fluid neutral interaction potential
    // Surface type parameters
    int nWall_type;                             // number of different surface types
    int nWall;                                  // total number of surfaces
    int nLink;                                  // number of independent compound surfaces
    int autoCenterSurfInBox;                    // automatically center the surfaces in the computational box
    String independentSurfaceNames[];           // the names of each of the independent surfaces
    int xTest_reflect_TF[][];                   // integer used as 1/0=T/F external surface prop.
    String surfaceName[];                       // the names of the surfaces
    int surf_type[];                            // array of surface types
    int orientation[];                          // array of orientation for each surface
    double wallParam[];                         // extra parameters for surfaces
    double wallParam2[];                        // 2nd extra parameters for surfaces
    double wallParam3[];                        // 3nd extra parameters for surfaces
    // Fluid Particle Parameters
    int nComp;                                  // number of components
    String compName[];                          // name of each component
    double compMass[];                          // mass of each component
    double charge_f[];                         // charges
    // Surface Particle Parameters
    double sigma_w[];                           // interaction diameters
    double eps_w[];                             // interaction energies
    double rho_s[];                             // number densities
    // Mixing Rule Parameters
    int lMix_rule;                              // mixing rule 0=LB 1=manual
    double sigma_ff[][];                        // interaction diameters
    double eps_ff[][];                          // interaction energies
    double cut_ff[][];                          // cut off lengths
    double sigma_wf[][];                        // interaction diameters
    double eps_wf[][];                          // interaction energies
    double cut_wf[][];                          // cut off lengths
    double sigma_ww[][];                        // interaction diameters
    double eps_ww[][];                          // interaction energies
    double cut_ww[][];                          // cut off lengths
    // Polyer Input Parameters
    int nPol_comp;                              // number of polymer components
    String cr_file;                             // file holding c(r) data
    private class polymerBlockInfo {
	int nBlock;                             // number of blocks in this polymer
	int block[];                            // number of segments in each block [nBlock]
	int block_type[];                       // segment types in each block
	String poly_file;                       // file holding polymer connectivity
	double cr_radius;                       // c(r) radius units (sigma)
        double gauss_a;                         // aspect ratio
	polymerBlockInfo( int numBlocks ) {
	    nBlock = numBlocks;
	    block = new int[ numBlocks ];
	    block_type = new int[ numBlocks ];
	}
    };
    java.util.List polymerVector;                // array for storing Polymer objects
    // Semi-Permeable Surface Parameters
    int lSemiPerm[][];                           // logical array 0/1=F/T [nWall_type][nComp] indicating if
                                                 // any wall is permeable to any component
    double vExtMembrane[][];                     // external potential [nWall_type][nComp]
    // State Point Parameters
    double temp;                                 // temperature 
    double tempElec;                             // temperature for electrical calc
    double bulkDensity[];                        // Rho_b [nComp]
    // Charged Surface Boundary Condition Parameters
    int typeBCElectric[];                        // type of electric boundary condition [nWall_type]
    double elecParamW[];                         // value of electrostatic boundary condition [nWall_type]
    int nLocalCharge;                            // number of local charges write -1 for two charges acting like
                                                 // a linear charge profile
    boolean chargeProfile;                       // boolean flag used to indicate case where nLocalCharge should
                                                 // be output as -1
    String chargeName[];                         // string to identify a chage
    double chargeLoc[];                          // total charge associated with a local charge; [nLocalCharge]
    double chargeDiam[];                         // diameter local charge should be spread over; [nLocalCharge]
    double chargeX[][];                          // charge position [nLocalCharge][nDim]
    int lPointChargeLocal;                       // how to treat local charges 
    int lPointChargeAtoms;                       // how to treak atomic charges
    // Dielectric Constant Parameters
    int typeDielec;                              // flag for how dielectric constants will be treated
                                                 // 0 = constant everywhere in domain
                                                 // 1 = surfaces and fluids different dielectrics
                                                 // 2 = give pore fluid and bulk fluid different dielectrics
                                                 // 3 = constant in walls but density dependent in fluids
    double dielecBulk;                           // bulk fluid dielectric
    double dielecPore;                           // pore fluid dielectric
    double dielecX;                              // distance from surfaces that will be considered pore fluid
    double dielecWall[];                         // array of surface dielectrics
    // Steady state boundary conditions
    int lSteadyState;                            // logical 0/1=F/T indicating a steady state problem
    int linearTransport;                         // flag 0 = nonlinear transport, 1 = linear transport
    int gradientDirection;                       // direction of chemical potential gradient (0=x, 1=y, 2=z)
    double xConstMu;                             // distance in gradient direction where chemical potential is constant
    int geometryFlag;                            // flag for one dimensional pores.  0=unit area, 1=cylc pore, 2=vary
    int n1DPoreSegments;                         // number of segments to this 1D pore
    double pore1DRadiusLeft[];                   // for a given segment number, left side radius
    double pore1DRadiusRight[];                  // for a given segment number, right side radius
    double pore1DLength[];                       // for a given segment number, segment length
    double bulkDensityLBB[];                     // for a given component, left, bottom, back side bulk density
    double bulkDensityRTF[];                     // for a given component, right, top, front bulk density
    double diffusionCoef[];                      // for a given component, the diffusion coeffiecient
    double elecPotLBB;                           // the electric potential on the left, bottom, back
    double elecPotRTF;                           // the electric potential on the right, top, front
    double centerOfMassVelocity;                 // the system center of mass velocity
    // Startup Control Parameters
    int iLiqVap;                                 // liquid vapor coexistence flag -2=no coex, -1=none, 1=W-V,
                                                 // 2=W-L, 3=L-v profiles
    int iGuess;                                  // initial guesses. -3,-2,-1= rho; constRho, constRhoL, constRhoV
                                                 // 0, 1, 2= rho*exp(-Vext/kt); expRho, expRhoL, expRhoV
                                                 // 3, 4, 5= rho*theta fnc; stepRho, stepRhoL, stepRhoV
                                                 // 6 = rho*theta fnc for liq-vap prof: stepLV
                                                 // 7 = linear interpolation between left and right
    double thickness;                            // distance from surface at which profile is stepped
    int restart;                                 // restart flag; 0=no, 1=yes, 2=yes w/ step function,
                                                 // 3=yes for densities but not elec.pot or chem.pot
    // mesh continuation parameters
    int numRuns;                                 // number of runs in mesh continuation
    double parameterChange[];                    // how to change box for mesh continuation 
    int planeNewNodes;                           // plane in which to insert new nodes (0=yz, 1=xz, 2=xy)
    int positionNewNodes;                        // position of new nodes (-1=lbb, 0=center, 1=rtf)
    double guessRange[];                         // guess range [0,1] for mesh continuation
    double rhoMax;                               // max density for mesh continuation
    // output format parameters
    int lPerArea;                                // logical 0/1=F/T to output values on a per unit area basis
    int lPrintGr;                                // logical 0/1=F/T to output radial distribution function
    int printRhoType;                            // sets how density data is printed
    int printRhoSwitch;                          // switch for printing density
    int printMeshSwitch;                         // switch for printing mesh
    int printHeader;                             // switch for printing a header
    int printIWrite;                             // formatting flag.
    // Coarsening switches
    int nZone;                                   // number of zones for coarsening
    double rMaxZone[];                           // distances from surface where zones are active
    int coarsenResidual;                         // 0=no, 1=yes
    int coarserJac;                              // type of jacobian coarensing
    double eSizeJacobian;                        // parameter used in jacobian coarsening
    int lJacCut;                                 // logical for jacobian intergral cutoff
    double jacThreshold;                         // value for inclusion of a point in an integration stencil
    int matrixFillFlag;                          // matrix filling options
    // nonlinear solver parameters
    int maxNewtonIter;                           // maximum number of newton iterations
    double newtonRelTol;                         // newton relative tolerance
    double newtonAbsTol;                         // newton absolute tolerance
    int loadBalanceSwitch;                       // load balance flag
    // linear solver parameters
    int solverFlag;                              // flag indicating which linear solver to use
    int kspace;                                  // storage flag for gmres routine.
    int scalingFlag;                             // flag for type of scaling to use
    int preconditionerFlag;                      // flag for type of preconditioner to use
    double iLutFillParam;                        // ilut fill parameter
    int maxLinearSolverIter;                      // storage flag for gmres routine.
    double convergenceTol;                       // convergence tolerance

    // loca continuation parameters
    int locaContinuationMethod;                  // flag -1=none, 0,1,2=0th, 1st, arc-length
                                                 // 3=spinodal (turning point), 4=binodal (phase eq)
    int locaContinuationParameter;               // 1=temp, 2=rho_bulk, 3=eps_wall, 4=cutoff, 
                                                 // 5=elec_param_wall, 6=rho_bulk_lbb, 7=rho_bulk_rtf
                                                 // 8=elec_pot_lbb, 9=elec_pot_rtf, 10=log(rho_bulk)
                                                 // 11=eps_wf[0][0], 12=eps_wf[all][0]
    double locaParameterStep;                    // parameter step size
    int locaNumSteps;                            // number of steps
    double locaStepControl;                      // step control parameter

    public TramontoXMLParser( XMLObject theMessage ) {
	inputMessage = theMessage;
    }

    public void validateXML() throws RuntimeException{
	
	// get the input file name from the xml and store it
	inputFileName = inputMessage.getAttribute("simulationInputFileName");

	// any of these routines can throw a RuntimeException.
	// the calling routine can catch that and let the user
	// know

	// basic layout of information in the xml message
	// dimensional units
	// fluid/component names
	// surface names
	// problem description
	//   problem type (subclass)
	//   surface types
	//   functionals and interactions
	//   surface interactions
	//   boundary conditions
	//   computational parameters

	validateDimensionParameters();
	validateFluidNames();
	validateSurfaceNames();
	validateProblemType();
	validateMeshParameters();
	validateFunctionalSwitchParameters();
	validatePotentialTypeParameters();
	validateSurfaceParameters();
	validateFluidInteractionParameters();
        validateSurfaceInteractionParameters();
	validatePolymerInputParameters();
	validateSemiPermeableSurfaceParameters(); 
	validateStatePointParameters();
	validateChargedSurfaceBoundaryConditionParameters();
	validateDielectricConstantParameters();
	validateSteadyStateBoundaryConditionParameters();
	validateStartupControlParameters();
	validateOutputFormatParameters();
	validateCoarseningSwitches();
	validateNonLinearSolverParameters();
	validateLinearSolverParameters();
	validateMeshContinuationParameters();
	validateLocaContinuationParameters();
    }

    public boolean generateInputFile() {
	boolean returnValue = true;
	try {
	    // use a JFileChooser to make sure we don't overwrite the output file
	    // and that we place the output in the correct place. 
	    JFileChooser chooser = new JFileChooser();
	    chooser.setCurrentDirectory(new File( System.getProperty("user.dir") +
						  System.getProperty("file.separator") ) );
	    chooser.setSelectedFile( new File(System.getProperty("user.dir") +
					     System.getProperty("file.separator") +
					     inputFileName ) );
	    chooser.setDialogTitle("Save generated input file as...");
	    int returnVal = chooser.showSaveDialog(null);
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
		inputFileName = chooser.getSelectedFile().getCanonicalPath();
		workingDirectory = chooser.getSelectedFile().getParent();
		PrintWriter theOutput=null;
		theOutput = new PrintWriter( new FileWriter( inputFileName ) );
		writeDimensionParameters( theOutput );
		writeFluidNames( theOutput );
		writeSurfaceNames( theOutput );
		writeProblemType( theOutput );
		writeMeshParameters( theOutput );
		writeFunctionalSwitches( theOutput );
		writePotentialTypeParameters( theOutput );
		writeSurfaceParameters( theOutput );
		writeFluidSurfaceInteractionParameters( theOutput );
		writePolymerInputParameters( theOutput );
		writeSemiPermeableSurfaceParameters( theOutput );
		writeStatePointParameters( theOutput );
		writeChargedSurfaceBoundaryConditionParameters( theOutput );
		writeDielectricConstantParameters( theOutput );
		writeSteadyStateBoundaryConditionParameters( theOutput );
		writeStartupControlParameters( theOutput );
		writeOutputFormatParameters( theOutput );
		writeCoarseningSwitches( theOutput );
		writeNonLinearSolverParameters( theOutput );
		writeLinearSolverParameters( theOutput );
		writeMeshContinuationParameters( theOutput );
		writeLocaContinuationParameters( theOutput );
		// close the output file
		theOutput.close();
	    }
	    else if(returnVal == JFileChooser.CANCEL_OPTION ) {
		returnValue = false;
	    }
 	}
	catch( IOException ioe ) {
	    returnValue = false;
	} 
	return returnValue;
    }


    //
    // below are private methods used to support this class
    //
    private int findIndexInStringArray( String token, String[] target ) throws RuntimeException {
	for(int i=0; i<target.length; i++) {
	    if( token.equals( target[i] ) ) {
		return i;
	    }
	}
	throw new RuntimeException( "Could not find, \"" + token + ",\" in \"" + target +"\"");
    }


    private boolean isStringInStringArray( String token, String[] target ) {
	for(int i=0; i<target.length; i++) {
	    if( token.equals( target[i] ) ) {
		return true;
	    }
	}
	return false;
    }


    private int guessStringToInteger( String theTag ) {
	int theGuess;  
	
	if( theTag.equals("GuessConstBulk") ) {
	    theGuess = -3;
	} 
	else if( theTag.equals("GuessConstLiquid") ) {
	    theGuess = -2;
	}
	else if( theTag.equals("GuessConstVapor") ) {
	    theGuess = -1;
	}
	else if( theTag.equals("GuessIdealBulk") ) {
	    theGuess = 0;
	}
	else if( theTag.equals("GuessIdealLiquid") ) {
	    theGuess = 1;
	}
	else if( theTag.equals("GuessIdealVapor") ) {
	    theGuess = 2;
	}
	else if( theTag.equals("GuessSteppedBulk") ) {
	    theGuess = 3;
	}
	else if( theTag.equals("GuessSteppedLiquid") ) {
	    theGuess = 4;
	}
	else if( theTag.equals("GuessSteppedVapor") ) {
	    theGuess = 5;
	}
	else if( theTag.equals("GuessSteppedLiquidVapor") ) {
	    theGuess = 6;
	}
	else if( theTag.equals("GuessLinearLiquidVapor") ) {
	    theGuess = 7;
	}
	else {
	    throw new RuntimeException("Invalid Iguess startup control option, \"" + theTag + "\"");
	}
	return theGuess;
    }


    private int mapAttributeToInteger( String token, String[] tokenSet, int[] valueSet ) 
	throws RuntimeException {
	if( tokenSet.length != valueSet.length ) {
	    throw new RuntimeException(
	        "mapAttributeToInteger: tokenSet.lenght != valueSet.length");
	}
        for(int i=0; i<tokenSet.length; i++) {
	    if( token.equals(tokenSet[i]) ) {
		return valueSet[i];
	    }
	}
	throw new RuntimeException("mapAttributeToInteger: token not found in tokenSet");
    }

    private int stringToLocaContinuationParameter( String theParameter ) {
	int returnValue;

	if( theParameter.equals("Temperature") ) {
	    returnValue = 1;
	} 
	else if( theParameter.equals("Bulk density 0") ) {
	    returnValue = 2;
	}
	else if( theParameter.equals("Bulk density all") ) {
	    returnValue = 3;
	}
	else if( theParameter.equals("Log density 0") ) {
	    returnValue = 4;
	}
	else if( theParameter.equals("Log density all") ) {
	    returnValue = 5;
	}
	else if( theParameter.equals("Scale density") ) {
	    returnValue = 6;
	}
	else if( theParameter.equals("Surface-surface interaction energy 0") ) {
	    returnValue = 7;
	}
	else if( theParameter.equals("Surface-surface interaction energy all") ) {
	    returnValue = 8;
	}
	else if( theParameter.equals("Scale surface-surface interaction energy") ) {
	    returnValue = 9;
	}
	else if( theParameter.equals("Surface-fluid interaction energy 0") ) {
	    returnValue = 10;
	}
	else if( theParameter.equals("Surface-fluid interaction energy all") ) {
	    returnValue = 11;
	}
	else if( theParameter.equals("Scale surface-fluid interaction energy") ) {
	    returnValue = 12;
	}
	else if( theParameter.equals("Fluid-fluid interaction energy 0") ) {
	    returnValue = 13;
	}
	else if( theParameter.equals("Fluid-fluid interaction energy all") ) {
	    returnValue = 14;
	}
	else if( theParameter.equals("Scale fluid-fluid interaction energy") ) {
	    returnValue = 15;
	}
	else if( theParameter.equals("Constant scale change") ) {
	    returnValue = 16;
	}
	else {
	    throw new RuntimeException("Unknown continuation parameter, \"" + theParameter + "\"");
	}
	return returnValue;
    }



    /*

      Folowing are the public routines validate...() sorted alphabetically.

    */


    public void validateChargedSurfaceBoundaryConditionParameters() throws RuntimeException {
        // creat arrays for nWall_type dependent data	
	typeBCElectric = new int[nWall_type];
	elecParamW = new double[nWall_type];

	XMLObject boundaryConditionsSection = problemTypeSection.getChild("BoundaryConditions" +
									  problemDimensionSuffix);
	if( boundaryConditionsSection != null ) {
	    XMLObject chargedSurfacesSection = boundaryConditionsSection.getChild("ChargedSurfaces" +
										  problemDimensionSuffix);
	    if( chargedSurfacesSection != null ) {
		if( chargedSurfacesSection.getAttribute("chargeProfile").equalsIgnoreCase("true") ) {
		    chargeProfile = true;
		}
		else {
		    chargeProfile = false;
		}

		if( chargedSurfacesSection.getAttribute("chargeTypeLocal").equals("Point charge") ) {
		    lPointChargeLocal = 0;
		}
		else if(chargedSurfacesSection.getAttribute("chargeTypeLocal").equals("Smeared charge over diameter") ) {
		    lPointChargeLocal = 1;
		}
		else if(chargedSurfacesSection.getAttribute("chargeTypeLocal").equals("Background charge") ) {
		    lPointChargeLocal = 2;
		}

		if( chargedSurfacesSection.getAttribute("chargeTypeAtoms").equals("Point charge") ) {
		    lPointChargeAtoms = 0;
		}
		else if(chargedSurfacesSection.getAttribute("chargeTypeAtoms").equals("Smeared charge over diameter") ) {
		    lPointChargeAtoms = 1;
		}
		else if(chargedSurfacesSection.getAttribute("chargeTypeAtoms").equals("Background charge") ) {
		    lPointChargeAtoms = 2;
		}

		// now there are two chlidren, both tables which we'll have to read
		for( int childIndex=0; childIndex < chargedSurfacesSection.numChildren(); childIndex++) {
		    XMLObject aTable=chargedSurfacesSection.getChild(childIndex);
		    if( aTable.getAttribute("name").equals("ChargedSurfaces") ) {
			for( int entryIndex=0; entryIndex<aTable.numChildren(); entryIndex++) {
			    XMLObject anEntry = aTable.getChild(entryIndex);
			    String aSurfaceName;
			    int aSurfaceIndex;
			    if( anEntry != null ) {
				XMLObject ref = anEntry.getChild("Reference");
				if( ref != null ) {
				    XMLObject panel = ref.getChild("PanelOfSurfaces");
				    if( panel != null ) {
					aSurfaceName = panel.getAttribute("surfaceName");
					try {
					    aSurfaceIndex = findIndexInStringArray( aSurfaceName, surfaceName );
					}
					catch (RuntimeException e) {
					    throw new RuntimeException("Couldn't find surface named, \"" + aSurfaceName +
								       "\" in list of surfaces");
					}
				    }
				    else {
					throw new RuntimeException("Error finding surface panel in charged surfaces section.");
				    }
				} 
				else {
				    throw new RuntimeException("Error finding surface reference in charged surfaces section.");
				}
				// ok, have a surface index now, so we can read the other data
				if( anEntry.getAttribute("typeBCElec").equals("Neutral surface") ) {
				    typeBCElectric[aSurfaceIndex] = 0;
				}
				else if( anEntry.getAttribute("typeBCElec").equals("Constant potential") ) {
				    typeBCElectric[aSurfaceIndex] = 1;
				}
				else if( anEntry.getAttribute("typeBCElec").equals("Constant surface charge") ) {
				    typeBCElectric[aSurfaceIndex] = 2;
				}
				else if( anEntry.getAttribute("typeBCElec").equals("Atomic charges") ) {
				    typeBCElectric[aSurfaceIndex] = 3;
				}
				elecParamW[aSurfaceIndex] = anEntry.getDouble("elecParameter");
			    }
			    else {
				throw new RuntimeException("Error finding table entity in charged surfaces section.");
			    }
			}
		    }
		    
		    else if( aTable.getAttribute("name").equals("LocalCharges") ) {
			nLocalCharge = aTable.numChildren();
			if( (chargeProfile == true) && (nLocalCharge != 2) ) {
			    throw new RuntimeException("A linear charge profile requires " +
						       "that only two local charges be defined.");
			}
			chargeName = new String[nLocalCharge];
			chargeLoc = new double[nLocalCharge];
			chargeDiam = new double[nLocalCharge];
			chargeX = new double[nLocalCharge][nDim];
			for(int charge=0; charge<aTable.numChildren(); charge++) {
			    XMLObject theEntry = aTable.getChild(charge);
			    chargeDiam[ charge ] = theEntry.getDouble("chargeDiameter");
			    chargeLoc[ charge ] = theEntry.getDouble("chargeValue");
			    chargeX[ charge ][0] = theEntry.getDouble("chargeX");
			    if( nDim > 1 ) {
				chargeX[ charge ][1] = theEntry.getDouble("chargeY");
			    }
			    if( nDim > 2 ) {
				chargeX[ charge ][2] = theEntry.getDouble("chargeZ");
			    }
			}
		    }
		}
	    }
	}
	else {
	    throw new RuntimeException("Couldn't find Boundary Conditions Section.");
	}

    }


    public void validateCoarseningSwitches() throws RuntimeException {
	XMLObject computationalParameters;

	if( (computationalParameters = problemTypeSection.getChild("ComputationalParameters" + 
								   problemDimensionSuffix )) != null ) {
	    XMLObject coarseningSwitchesSection;
	    if( (coarseningSwitchesSection = computationalParameters.getChild("CoarseningSwitches")) != null ) {

		// coarsenResidual
		if( coarseningSwitchesSection.getAttribute("coarsenResidual").equals("true") ) {
		    coarsenResidual = 1;
		}
		else {
		    coarsenResidual = 0;
		}

		// matrixFillFlag
		if( coarseningSwitchesSection.getAttribute("matrixFillFlag").equals("Use an exact Jacobian with embeded nonlocal density") ) {
		    matrixFillFlag = 0;
		}
	        else if( coarseningSwitchesSection.getAttribute("matrixFillFlag").equals("Use minimal set Jacobian based on w(3) and w(2) contributions") ) {
		    matrixFillFlag = 1;
		}
	        else if( coarseningSwitchesSection.getAttribute("matrixFillFlag").equals("Use minimal set Jacobian with all scalar terms and no vector terms") ) {
		    matrixFillFlag = 2;
		}
	        else if( coarseningSwitchesSection.getAttribute("matrixFillFlag").equals("Use enumerated nonlocal density formulation") ) {
		    matrixFillFlag = 3;
		}
	        else if( coarseningSwitchesSection.getAttribute("matrixFillFlag").equals("Use enumerated nonlocal density formulation including scalar contributions in the Jacobian") ) {
		    matrixFillFlag = 4;
		}
		else {
		    throw new RuntimeException("Unknown matrix fill flag in coarsening switches.");
		}		

		// loop over the children and gather data from them

		for( int cChild=0; cChild < coarseningSwitchesSection.numChildren(); cChild++ ) {
		    XMLObject theCChild = coarseningSwitchesSection.getChild(cChild);
		    String theTag = theCChild.getTag();

		    // nZone & rMaxZone
                    if( (theTag.equals("Table")) &&
			(theCChild.getAttribute("name").equals("zoneTable")) ) {
			nZone = theCChild.numChildren();
			if( nZone == 0 ) {
			    throw new RuntimeException("At least one coarsening zone must be defined under " +
						       "Computational Parameters, Coarsening Switches." );
			}
			rMaxZone = new double[ nZone ];
			for( int zone=0; zone < nZone; zone++ ) {
			    XMLObject entry = theCChild.getChild( zone );
			    rMaxZone[ zone ] = entry.getDouble("distance");
			}
		    }  
		    
		    // coarserJac
		    else if(theTag.equals("NoJacobianCoarsening")) {
			coarserJac = 0;
			eSizeJacobian = 0;
		    }
		    else if(theTag.equals("FinestZoneJacobianCoarsening") ) {
			coarserJac = 1;
			eSizeJacobian = 0;
		    }
		    else if(theTag.equals("AllButMostCoarseJacobianCoarsening") ) {
			coarserJac = 2;
			eSizeJacobian = 0;
		    }
		    else if(theTag.equals("UseCoarsestZoneJacobianCoarsening") ) {
			coarserJac = 3;
			eSizeJacobian = 0;
		    }
		    else if(theTag.equals("UsePenultimateZoneJacobianCoarsening") ) {
			coarserJac = 4;
			eSizeJacobian = 0;
		    }
		    else if(theTag.equals("UseConstantMeshSpacingJacobianCoarsening") ) {
			coarserJac = 5;
			eSizeJacobian = 0;
			eSizeJacobian = theCChild.getDouble("eSizeJacobian");
		    }
		    
		    // lJacCut
		    else if(theTag.equals("NoJacobianCutOff")) {
			lJacCut = 0;
			jacThreshold = 0;
		    }
		    else if(theTag.equals("YesJacobianCutOff")) {
			lJacCut = 1;
			jacThreshold = theCChild.getDouble("threshold");
		    }
		    
		    // catch bad input messages
		    else {
			throw new RuntimeException("Unknown subsection in coarsening zones section.");
		    }
		}
	    }
	    else {
		throw new RuntimeException("No valid coarsening switches section found.");
	    }
	}
	else {
	    throw new RuntimeException("No valid computational parameters section found.");
	}
    }

    public void validateDielectricConstantParameters() throws RuntimeException {
	dielecWall = new double[nWall_type];

	XMLObject boundaryConditionsSection = problemTypeSection.getChild("BoundaryConditions" +
									  problemDimensionSuffix);
	if( boundaryConditionsSection != null ) {
	    XMLObject dielectricObject = boundaryConditionsSection.getChild("DielectricParameters");
	    if( dielectricObject != null) {
		XMLObject dielectricChild = dielectricObject.getChild(0);
		if( dielectricChild.getTag().equals("DielectricConstantEverywhere") ) {
		    // not much to do here
		    typeDielec = 0;
		    dielecBulk = dielectricChild.getDouble("dielectricBulk");
		}
		else {
		    if( dielectricChild.getTag().equals("DielectricFluidWallsDifferent")) {
			typeDielec = 1;
		    }
		    else if( dielectricChild.getTag().equals("DielectricBulkWallPore")) {
			typeDielec = 2;
		    }
		    else if( dielectricChild.getTag().equals("DielectricConstantWallsVariesFluid")) {
			typeDielec = 3;
		    }

		    if( dielectricChild.hasAttribute("dielectricBulk") ) {
			dielecBulk = dielectricChild.getDouble("dielectricBulk");
		    }
		    if( dielectricChild.hasAttribute("dielectricDistance") ) {
			dielecX = dielectricChild.getDouble("dielectricDistance");
		    }
		    if( dielectricChild.hasAttribute("dielectricWall") ) {
			dielecPore = dielectricChild.getDouble("dielectricWall");
		    }

		    // now grab the wall dielectrics
		    XMLObject aTable = dielectricChild.getChild("Table");
		    if( (aTable != null) &&
			(aTable.getAttribute("name").equals("SurfacesDielectrics")) ) {
			String aSurfaceName;
			int aSurfaceIndex;
			for( int entryIndex=0; entryIndex < aTable.numChildren(); entryIndex++ ) {
			    XMLObject anEntry = aTable.getChild( entryIndex );
			    if( anEntry != null ) {
				XMLObject aRef = anEntry.getChild("Reference");
				if( aRef != null ) {
				    XMLObject panel = aRef.getChild("PanelOfSurfaces");
				    if( panel != null ) {
					aSurfaceName = panel.getAttribute("surfaceName");
					try {
					    aSurfaceIndex = findIndexInStringArray( aSurfaceName, surfaceName );
					}
					catch (RuntimeException e) {
					    throw new RuntimeException("Couldn't find surface, \"" + aSurfaceName +
								       "\" in list of surfaces.");
					}
				    }
				    else {
					throw new RuntimeException("Couldn't find surface panel in dielectrics table.");
				    }
				}
				else {
				    throw new RuntimeException("Couldn't find reference in dielectrics table.");
				}
			    }
			    else {
				throw new RuntimeException("Couldn't find a table entry in the dielectrics table.");
			    }
			    // ok have a surface index, read the dielectric attched to it and place it in the array
			    dielecWall[ aSurfaceIndex ] = anEntry.getDouble("dielectricWall");
			}
		    }
		}
	    }
	    else {
		throw new RuntimeException("Could not find dielectric section.");
	    }
	}
	else {
	    throw new RuntimeException("Could not find Boundary Conditions section.");
	}
    }

    public void validateDimensionParameters() throws RuntimeException {
	XMLObject stateParametersObject;

	// grab the statePareamets section
	if ( (stateParametersObject = inputMessage.getChild("TramontoRealUnits")) != null ) {
	    referenceLength = stateParametersObject.getDouble("refLength");
	    referenceTemp = stateParametersObject.getDouble("refTemperature");
	    referenceDensity = stateParametersObject.getDouble("refDensity");
	    referenceDielectric = stateParametersObject.getDouble("refDielectric");
	    // maximumPotential is read in the validatePotentialTypeParameters()
	}
	else if((stateParametersObject = inputMessage.getChild("TramontoReducedUnits")) != null ) { 
	    referenceLength = -1.0;
	    referenceTemp = -1.0;
	    referenceDensity = -1.0;
	    referenceDielectric = -1.0;
	    // maximumPotential is read in the validatePotentialTypeParameters()
	}
	else {
	    throw new RuntimeException( "Could not find the Dimensional Units  section" );
	}
    }


    public void validateFluidNames() throws RuntimeException {
	XMLObject fluidNamesObject;
	int fluidIndex;

	// grab the Fluid parameters section of the inputmessage
	if( (fluidNamesObject = inputMessage.getChild("FluidNames")) != null ) {
	    XMLObject componentsArray = fluidNamesObject.getChild("Array");
	    if( (componentsArray != null) && 
		(componentsArray.getAttribute("name").equals("componentsArray")) ) {
		nComp = componentsArray.numChildren();
		compName = new String[nComp];

		for(int i=0; i<nComp; i++ ) {
		    XMLObject panelOfComponentData = componentsArray.getChild(i);
		    /* need to check that name is not blank and unique */
		    String potentialCompName = panelOfComponentData.getAttribute("name").trim();
		    if( potentialCompName.equals("") ) {
			throw new RuntimeException( "One of the fluid/components has a blank name. " +
						    "All components must have unique names." );
		    }
		    if( isStringInStringArray(  potentialCompName, compName )  == true ) {
			throw new RuntimeException( "Fluid/Component \"" + potentialCompName + 
						    "\" is duplicated in list of components. " +
						    "All components must have unique names.");
		    }

		    /* made it here so component is named and unique */
		    compName[i] = panelOfComponentData.getAttribute("name").trim();
		}
	    }
	    else {
		throw new RuntimeException( "No components found in Fluid/Components Name Section" );
	    }
	}
	else {
	    throw new RuntimeException( "No valid Fluid/Components names section found." );
	}

	if( nComp == 0 ) {
	    throw new RuntimeException( "No components defined in the Fluid/Components name section." );
	}
    }


    public void validateFluidInteractionParameters() throws RuntimeException {
	int fluidIndex;

	// validateFluidNames has already determined nComp.

	// make the arrays to hold data
	sigma_ff = new double[nComp][nComp];
	eps_ff = new double[nComp][nComp];
	cut_ff = new double[nComp][nComp];
	compMass = new double[nComp];
	charge_f = new double[nComp];

	// make a boolean array to track what data has been read
	boolean dataFound[][] = new boolean[nComp][nComp];
	boolean massDataFound[] = new boolean[nComp];

	// the functionals and interactions section
	// will have a child whose name is of the format:
	// MixingType, LB or UD, Sigma and/or Eps and/or Cut and/or Charge, FluidTable.
	// for example:  MixingTypeLBSigmaChargeFluidTable
	// or MixingTypeUDSigmaEpsCutFluidTable
	// the one exception is the section NoHSNoAttNoCoulomb which if we find
	// we can ignore most fluid interactions.
	if( functionalsAndInteractionsSection.getTag().equals("NoHSNoAttNoCoulomb") ) {
	    // this is the one case where we don't need fluid parameters
	    lMix_rule = 0;
	    for(int j=0; j<nComp; j++) {
		dataFound[j][j] = true;
	    }
	    // check for debroglie mass section
	    XMLObject massSection = functionalsAndInteractionsSection.getChild("NoDeBroglieWavelengthInChemPot");
	    if( massSection != null ) {
		for(int j=0; j<nComp; j++) {
		    compMass[j] = 1.0;
		    massDataFound[j] = true;
		}
	    }
	    massSection = functionalsAndInteractionsSection.getChild("IncludeDeBroglieWavelengthInChemPot");
	    if( massSection != null ) {
		// now read some masses
		XMLObject massFluidTable = massSection.getChild("MassFluidTable");
		if( massFluidTable != null ) {
		    XMLObject theTable = massFluidTable.getChild("Table");
		    if( theTable != null ) {
			for( int j = 0; j<theTable.numChildren(); j++) {
			    XMLObject dataNode = theTable.getChild(j);
			    if( dataNode != null ) {
				// enough digging now we can read some data
				XMLObject component1 = dataNode.getChild(0).getChild(0);
				String component1Name = "";
				int component1Index;
				if( component1 != null ) {
				    component1Name = component1.getAttribute("name");
				    try {
					component1Index = findIndexInStringArray( component1Name, compName );
				    }
				    catch (RuntimeException e) {
					throw new RuntimeException("Functionals and Interaction Data table: " +
								   "Couldn't find component \"" +
								   component1Name + 
								   "\"in the list of coomponent names ");
				    }
				} 
				else {
				    throw new RuntimeException("Couldn't find a reference to a component name " +
							       "in row, " + j + " of the data table in the " +
							       "Functionals and Interaction section");
				}
				if( massDataFound[component1Index] == true ) {
				    throw new RuntimeException("In the Functionals and Interactions mass data table, " +
							       "data for component \"" + component1Name +
							       "\" appears more than once.");
				}
				compMass[ component1Index ] = dataNode.getDouble("massF");
				massDataFound[ component1Index ] = true;
			    }
			    else {
				throw new RuntimeException("Couldn't read table row " + 
							   j + " in Functionals and Interactions mass section.");
			    }
			}
		    }
		    else {
			throw new RuntimeException("Couldn't find mass table in Functionals and interactions section.");
		    }
		}
		else {
		    throw new RuntimeException("");
		}
	    }
	}
	
	// so we'll look at all the children of the functionalsAndInteractionsSection
	// looking for one whose tag starts with MixingType
	for( int i=0; i<functionalsAndInteractionsSection.numChildren(); i++ ) {
	    XMLObject mixingSection = functionalsAndInteractionsSection.getChild(i);
	    if( mixingSection.getTag().matches("MixingType.*") ) {
		// we found the mixing section
		// now test if the user requesed LB or UD mixing
		if( mixingSection.getTag().matches("MixingTypeLB.*") ){
		    lMix_rule = 0;
		}
		else if( mixingSection.getTag().matches("MixingTypeUD.*") ) {
		    lMix_rule = 1;
		}

		if( mixingSection.getTag().matches(".*Polymer.*") ) {
		    // we'll skip the mass section 
		    for(int j=0; j<nComp; j++) {
			compMass[j]=1.0;
			massDataFound[j] = true;
		    }
		}

		// figure out what the table will contain and use
		// these flags indicate what data we will read
		boolean readSigma = false;
		boolean readEps = false;
		boolean readCut = false;
		boolean readCharge = false;
		if( mixingSection.getTag().matches(".*Sigma.*") ) {
		    readSigma = true;
		}
		if( mixingSection.getTag().matches(".*Eps.*") ) {
		    readEps = true;
		}
		if( mixingSection.getTag().matches(".*Cut.*") ) {
		    readCut = true;
		}
		if( mixingSection.getTag().matches(".*Charge.*") ) {
		    readCharge = true;
		}
		
		// now lets read some data!
		// here we should have two children.  One holds the table of 
		// interaction parameters and the other the table of masses
		//
		for( int child=0; child<mixingSection.numChildren(); child++) {
		    XMLObject mixingSectionChild = mixingSection.getChild(child);

		    // look for specific cases then the general interaction parameters
		    if( mixingSectionChild.getTag().equals("IncludeDeBroglieWavelengthInChemPot") ) {
			XMLObject massFluidTable = mixingSectionChild.getChild("MassFluidTable");
			if( massFluidTable != null ) {
			    XMLObject theTable = massFluidTable.getChild("Table");
			    if( theTable != null ) {
				for( int j = 0; j<theTable.numChildren(); j++) {
				    XMLObject dataNode = theTable.getChild(j);
				    if( dataNode != null ) {
					// enough digging now we can read some data
					XMLObject component1 = dataNode.getChild(0).getChild(0);
					String component1Name = "";
					int component1Index;
					if( component1 != null ) {
					    component1Name = component1.getAttribute("name");
					    try {
						component1Index = findIndexInStringArray( component1Name, compName );
					    }
					    catch (RuntimeException e) {
						throw new RuntimeException("Functionals and Interaction Data table: " +
									   "Couldn't find component \"" +
									   component1Name + 
									   "\"in the list of coomponent names ");
					    }
					} 
					else {
					    throw new RuntimeException("Couldn't find a reference to a component name " +
								       "in row, " + j + " of the data table in the " +
								       "Functionals and Interaction section");
					}
					if( massDataFound[component1Index] == true ) {
					    throw new RuntimeException("In the Functionals and Interactions mass data table, " +
								       "data for component \"" + component1Name +
								       "\" appears more than once.");
					}
					compMass[ component1Index ] = dataNode.getDouble("massF");
					massDataFound[ component1Index ] = true;
				    }
				    else {
					throw new RuntimeException("Couldn't read table row " + 
								   i + " in Functionals and Interactions mass section.");
				    }
				}
			    }
			    else {
				throw new RuntimeException("Couldn't find mass table in Functionals and interactions section.");
			    }
			    
			}
			else {
			    throw new RuntimeException("Couldn't find mass fluid table in Functionals and interactions section.");
			}
		    }
		    else if( mixingSectionChild.getTag().equals("NoDeBroglieWavelengthInChemPot") ) {
			// load up compMass array with 1.0's as that's the default value
			// when not using the deBroglie wavelength factor
			for(int j=0; j<nComp; j++) {
			    compMass[j]=1.0;
			    massDataFound[j] = true;
			}
		    }
		    else if( mixingSectionChild != null ) {
			XMLObject dataTable = mixingSectionChild.getChild("Table");
			if( dataTable != null ) {
			    for( int j = 0; j<dataTable.numChildren(); j++) {
				XMLObject dataNode = dataTable.getChild(j);
				if( dataNode != null ) {
				    // enough digging now we can read some data
				    XMLObject component1 = dataNode.getChild(0).getChild(0);
				    String component1Name = "";
				    int component1Index;
				    if( component1 != null ) {
					component1Name = component1.getAttribute("name");
					try {
					    component1Index = findIndexInStringArray( component1Name, compName );
					}
					catch (RuntimeException e) {
					    throw new RuntimeException("Functionals and Interaction Data table: " +
								       "Couldn't find component \"" +
								       component1Name + 
								       "\"in the list of coomponent names ");
					}
				    } 
				    else {
					throw new RuntimeException("Couldn't find a reference to a component name " +
								   "in row, " + j + " of the data table in the " +
								   "Functionals and Interaction section");
				    }
				    XMLObject component2;
				    String component2Name = "";
				    int component2Index;
				    if( lMix_rule == 1 ) {
					// look for component 2
					component2 = dataNode.getChild(1).getChild(0);
					if( component2 != null ) {
					    component2Name = component2.getAttribute("name");
					    try {
						component2Index = findIndexInStringArray( component2Name, compName );
					    }
					    catch (RuntimeException e) {
						throw new RuntimeException("Functionals and Interaction Data table: " +
									   "Couldn't find component \"" +
									   component2Name + 
									   "\"in the list of coomponent names ");
					    }
					} 
					else {
					    throw new RuntimeException("Couldn't find a reference to a component name " +
								       "in row, " + j + " of the data table in the " +
								       "Functionals and Interaction section");
					}
					
				    }
				    else {
					// for regular mixing, only diagonal elements are needed 
					// in the matrix of parameters like sigma_ff[][]
					component2Index = component1Index;
				    }
				    // check that this data hasn't alread been read
				    if( (dataFound[component1Index][component2Index] == true) &&
					(problemTypeSection.getTag().matches(".*Atomic.*")) ) {
					if( lMix_rule == 0 ) {
					    throw new RuntimeException("In the Functionals and Interactions data table, " +
								       "data for component \"" + component1Name +
								   "\" appears more than once.");
					}
					else {
					    throw new RuntimeException("In the Functionals and Interactions data table, " +
								       "data for component pair \"" + component1Name +
								       "\", \"" + component2Name +
								       "\" appears more than once.");
					}
				    }
				    else {
					dataFound[component1Index][component2Index] = true;
				    }
				    
				    
				    // read some data
				    if( readSigma == true ) {
					sigma_ff[component1Index][component2Index] = dataNode.getDouble("sigmaFF");
				    }
				    if( readEps == true ) {
					eps_ff[component1Index][component2Index] = dataNode.getDouble("epsFF");
				    }
				    if( readCut == true ) {
					cut_ff[component1Index][component2Index] = dataNode.getDouble("cutFF");
				    }
				    if( readCharge == true ) {
					charge_f[component1Index] = dataNode.getDouble("chargeF");
				    }
				    
				}
				else {
				    throw new RuntimeException("Couldn't read table row " + 
							       i + " in Functionals and Interactions section.");
				}
			    }
			}
		    }
		    else {
			throw new RuntimeException("Coudn't find child " + i + " in Functionals and interactions section.");
		    }
		}
	    }
	}

	// now check that we found all the needed data
	for(int i=0; i< nComp; i++) {
	    if( lMix_rule != 0 ) {
		for(int j=0; j< nComp; j++) {
		    if( dataFound[i][j] != true ) {
			// check for symmetric case
			if( dataFound[j][i] == true ) {
			    dataFound[i][j] = true;
			    sigma_ff[i][j] = sigma_ff[j][i];
			    eps_ff[i][j] = eps_ff[j][i];
			    cut_ff[i][j] = cut_ff[j][i];
			}
			else {
			    throw new RuntimeException("In Functionals and Interactions, interaction data for was not entered for components \"" +
						       compName[i] + "\", \"" + compName[j] + "\"" );
			}
		    }
		    
		}
	    }
	    else {
		if( dataFound[i][i] != true ) {
		throw new RuntimeException("In Functionals and Interactions, interaction data for was not entered for component \"" +
					   compName[i] + "\"");
		}
	    }
	    if (massDataFound[i] != true ) {
		throw new RuntimeException("In Functionals and Interactions, mass data for was not entered for component \"" +
					   compName[i] + "\"");
	    }
	}

    }
	
    public void validateFunctionalSwitchParameters() throws RuntimeException {
	XMLObject switchData;

	// functional switches are selected by picking one of several 
	// classes as follows:
	//     Class:                            Switch Status  
	//                           type_func  type_attr  type_coul  type_poly
	// NoHSNoAttNoCoulomb           -1         -1         -1         -1
        // NoHSNoAttStrictCoulomb       -1         -1          0         -1
	// NoHSNoAttSecondCoulomb       -1         -1          1         -1
	// NoHSStMFNoCoulomb            -1          0         -1         -1
	// NoHSStMFStrictCoulomb        -1          0          0         -1
	// NoHSStMFSecondCoulomb        -1          0          1         -1
	// RosenOneNoAttNoCoulomb        0         -1         -1         -1
        // RosenOneNoAttStrictCoulomb    0         -1          0         -1
	// RosenOneNoAttSecondCoulomb    0         -1          1         -1
	// RosenOneStMFNoCoulomb         0          0         -1         -1
	// RosenOneStMFStrictCoulomb     0          0          0         -1
	// RosenOneStMFSecondCoulomb     0          0          1         -1
	// RosenTwoNoAttNoCoulomb        1         -1         -1         -1
        // RosenTwoNoAttStrictCoulomb    1         -1          0         -1
	// RosenTwoNoAttSecondCoulomb    1         -1          1         -1
	// RosenTwoStMFNoCoulomb         1          0         -1         -1
	// RosenTwoStMFStrictCoulomb     1          0          0         -1
	// RosenTwoStMFSecondCoulomb     1          0          1         -1
	// PolymerModel0                -1         -1         -1          0
	// PolymerModel1                -1         -1         -1          1
	// PolymerModel2                -1         -1         -1          2

	String classNames[] = { "NoHSNoAttNoCoulomb",        "-1", "-1", "-1", "-1",
				"NoHSNoAttStrictCoulomb",    "-1", "-1",  "0", "-1",
				"NoHSNoAttSecondCoulomb",    "-1", "-1",  "1", "-1",
				"NoHSStMFNoCoulomb",         "-1",  "0", "-1", "-1",
				"NoHSStMFStrictCoulomb",     "-1",  "0",  "0", "-1",
				"NoHSStMFSecondCoulomb",     "-1",  "0",  "1", "-1",
				"RosenOneNoAttNoCoulomb",     "0", "-1", "-1", "-1",
				"RosenOneNoAttStrictCoulomb", "0", "-1",  "0", "-1",
				"RosenOneNoAttSecondCoulomb", "0", "-1",  "1", "-1",
				"RosenOneStMFNoCoulomb",      "0",  "0", "-1", "-1",
				"RosenOneStMFStrictCoulomb",  "0",  "0",  "0", "-1",
				"RosenOneStMFSecondCoulomb",  "0",  "0",  "1", "-1",
				"RosenTwoNoAttNoCoulomb",     "1", "-1", "-1", "-1",
				"RosenTwoNoAttStrictCoulomb", "1", "-1",  "0", "-1",
				"RosenTwoNoAttSecondCoulomb", "1", "-1",  "1", "-1",
				"RosenTwoStMFNoCoulomb",      "1",  "0", "-1", "-1",
				"RosenTwoStMFStrictCoulomb",  "1",  "0",  "0", "-1",
				"RosenTwoStMFSecondCoulomb",  "1",  "0",  "1", "-1",
				"PolymerModel0",             "-1", "-1", "-1",  "0",
				"PolymerModel1",             "-1", "-1", "-1",  "1",
				"PolymerModel2",             "-1", "-1", "-1",  "2" };
	
	boolean foundData = false;

	for( int i=0; i < classNames.length; i+=5 )
	    {
		if( (switchData = problemTypeSection.getChild( classNames[i] )) != null ) {
		    // fouond a matching class.  Remember it so that other routines
                    // can pull data from this seciton
                    functionalsAndInteractionsSection = switchData;

		    // so set the type_ variable flags and exit this loop
		    foundData = true;
		    type_func = Integer.parseInt( classNames[++i] );
		    type_attr = Integer.parseInt( classNames[++i] );
		    type_coul = Integer.parseInt( classNames[++i] );
		    type_poly = Integer.parseInt( classNames[++i] );
		    break;
		}
	    }
	
	if( foundData == false ) {
	    throw new RuntimeException( "No valid functional selections found in the Problem Type section." );
	}

	// finally, look in the computational parameters section for the compare to fastram flag.
	XMLObject computationalSection;
	if( (computationalSection = problemTypeSection.getChild( "ComputationalParameters" +
								 problemDimensionSuffix ) ) != null ) {
	    XMLObject outputFormatSection;
	    if( (outputFormatSection = computationalSection.getChild("OutputFormat")) != null ) {
		// now look for the compare with fastram flag
		if( outputFormatSection.getAttribute("compareFastram").equals("true") ) {
		    compareToFastram = 1;
		}
		else {
		    compareToFastram = 0;
		}
	    }
	    else {
		throw new RuntimeException( "Could not find the Output format section.");
	    }
	}
	else {
	    throw new RuntimeException( "Could not find the Computational parameters section.");
	}
    }


    public void validateLinearSolverParameters() throws RuntimeException {
	XMLObject computationalParameters;

	if( (computationalParameters = problemTypeSection.getChild("ComputationalParameters" + 
							     problemDimensionSuffix )) != null ) {
	    XMLObject linearSolverSection;
	    if( (linearSolverSection = computationalParameters.getChild("LinearSolver")) != null ) {
		
		// solverFlag
		if( linearSolverSection.getAttribute("solver").equals("grmes") ) {
		    solverFlag = 0;
		}
		else if( linearSolverSection.getAttribute("solver").equals("cg") ) {
		    solverFlag = 1;
		}
		else if( linearSolverSection.getAttribute("solver").equals("tfqmr") ) {
		    solverFlag = 2;
		}
		else if( linearSolverSection.getAttribute("solver").equals("cg2") ) {
		    solverFlag = 3;
		}
		else if( linearSolverSection.getAttribute("solver").equals("bicgstab") ) {
		    solverFlag = 4;
		}
		else {
		    throw new RuntimeException("Unknown solver in linear solver in Computational Parameters Section.");
		}

		// scalingFlag
		if( linearSolverSection.getAttribute("scaling").equals("none") ) {
		    scalingFlag = -1;
		}
		else if( linearSolverSection.getAttribute("scaling").equals("row-sum") ) {
		    scalingFlag = 0;
		}
		else if( linearSolverSection.getAttribute("scaling").equals("Jacobian") ) {
		    scalingFlag = 1;
		}
		else if( linearSolverSection.getAttribute("scaling").equals("symmetric row-sum") ) {
		    scalingFlag = 2;
		}
		else {
		    throw new RuntimeException("Unknown scaling in linear solver section.");
		}

		// preconditionerFlag
		if( linearSolverSection.getAttribute("preconditioner").equals("none") ) {
		    preconditionerFlag = -1;
		}
		else if( linearSolverSection.getAttribute("preconditioner").equals("ilu") ) {
		    preconditionerFlag = 0;
		}
		else if( linearSolverSection.getAttribute("preconditioner").equals("Jacobian") ) {
		    preconditionerFlag = 1;
		}
		else if( linearSolverSection.getAttribute("preconditioner").equals("symGS") ) {
		    preconditionerFlag = 2;
		}
		else if( linearSolverSection.getAttribute("preconditioner").equals("LSpoly3") ) {
		    preconditionerFlag = 3;
		}
		else if( linearSolverSection.getAttribute("preconditioner").equals("ilut (3 levels)") ) {
		    preconditionerFlag = 4;
		}
		else if( linearSolverSection.getAttribute("preconditioner").equals("ilut (7 levels)") ) {
		    preconditionerFlag = 5;
		}
		else {
		    throw new RuntimeException("Unknown preconditioner in linear solver section.");
		}
		
		iLutFillParam = linearSolverSection.getDouble("ilutfillparam");
		convergenceTol = linearSolverSection.getDouble("convergenceTolerance");
		maxLinearSolverIter = linearSolverSection.getInt("maxIterations");
		kspace = linearSolverSection.getInt("kspace");
	    }
	    else {
		throw new RuntimeException("No linear solver parameters found in Computational Parameters Section");
	    }
	}
	else {
	    throw new RuntimeException("No Computational Parameters Section found.");
	}

    }


    public void validateLocaContinuationParameters() throws RuntimeException {
	XMLObject computationalParameters;

	if( (computationalParameters = problemTypeSection.getChild("ComputationalParameters" +
								   problemDimensionSuffix )) != null ) {
	    XMLObject locaContinuationSection;
	    if( (locaContinuationSection = computationalParameters.getChild("LocaContinuation")) != null ) {
		XMLObject theContinuation;
		if( (theContinuation = locaContinuationSection.getChild(0)) != null ) {
		    // get the simple parameters first
		    if ( theContinuation.getTag().equals("NoContinuationMethod") ) {
			locaContinuationMethod = -1;
		    }
		    else if ( theContinuation.getTag().equals("ZerothContinuationMethod") ) {
			locaContinuationMethod = 0;
		    }
		    else if ( theContinuation.getTag().equals("FirstContinuationMethod") ) {
			locaContinuationMethod = 1;
		    }
		    else if ( theContinuation.getTag().equals("ArcLengthContinuationMethod") ) {
			locaContinuationMethod = 2;
		    }
		    else if ( theContinuation.getTag().equals("SpinodalContinuationMethod") ) {
			locaContinuationMethod = 3;
		    }
		    else if ( theContinuation.getTag().equals("BinodalContinuationMethod") ) {
			locaContinuationMethod = 4;
		    }
		    else {
			throw new RuntimeException("Could not determine cointinuation type in loca continuation. type: " + 
						   theContinuation.getTag());
		    }
		    if( locaContinuationMethod != -1 ) {
			// if there is a continuation desired, then try to read the remaining parameters.
			locaNumSteps = theContinuation.getInt("numberSteps");
			locaParameterStep = theContinuation.getDouble("parameterStepSize");
			locaStepControl = theContinuation.getDouble("stepControl");
			locaContinuationParameter = stringToLocaContinuationParameter(theContinuation.getAttribute("continuationParameter"));
		    }

		}
		else {
		    throw new RuntimeException("Cannot find a valid continuation type in loca continuation section.");
		}
	    }
	    else {
		throw new RuntimeException("Cannot find a valid loca continuation section.");
	    }
	}
	else {
	    throw new RuntimeException("Cannot find a valid computational parameters section");
	}
    }


    public void validateMeshParameters() throws RuntimeException {
	XMLObject meshData;

	if( (meshData = problemTypeSection.getChild( "ComputationalParameters" + problemDimensionSuffix )) != null ) {
	    // if( (meshData = inputMessage.getChild("TramontoComputationalParameters")) != null ) {
	    realSpaceComputation = true;
            if( meshData.getChild("OneDMesh") != null ) {
		nDim = 1;
		// step into child to get mesh details
		meshData = meshData.getChild("OneDMesh");
	    }
	    else if( meshData.getChild("TwoDMesh") != null ) {
		nDim = 2;
		// step into child to get mesh details
		meshData = meshData.getChild("TwoDMesh");
	    }
	    else if( meshData.getChild("ThreeDMesh") != null ) {
		nDim = 3;
		// step into child to get mesh details
		meshData = meshData.getChild("ThreeDMesh");
	    }
	    else {
		throw new RuntimeException( "Invalid mesh dimension" );
	    }

	    size_x = new double[nDim];
	    esize_x = new double[nDim];
            type_bc_lbb = new int[nDim];
            type_bc_rtf = new int[nDim];


	    String[] children = {"Xdim", "Ydim", "Zdim"};
	    String[] boundaryLocations = {"left", "right", "down", "up", "back", "front"};

	    for(int i=0; i < meshData.numChildren(); i++) {
		XMLObject theChild = meshData.getChild(children[i]);
		double sigma = theChild.getDouble("sigma");
		if( sigma <= 0.0 ) {
		    throw new RuntimeException( "Sigma must be greater than zero in mesh parameters, Computational Parameters Section" );
		}
		double domainSize = theChild.getDouble("domainSize");
		double meshDivisions = theChild.getDouble("meshDivisions");
		size_x[i] = domainSize * sigma;
		esize_x[i] = sigma / meshDivisions;
		type_bc_lbb[i] = mapAttributeToInteger( theChild.getAttribute("bclbb"),
		    new String[] {"Wall", "Bulk fluid", "Periodic", "Reflective"}, 
                    new int[] {-1, 0, 1, 2});
                type_bc_rtf[i] = mapAttributeToInteger( theChild.getAttribute("bcrtf"),
		    new String[] {"Wall", "Bulk fluid", "Periodic", "Reflective"}, 
                    new int[] {-1, 0, 1, 2});
	    }

	}
	else {
	    throw new RuntimeException( "No mesh data found in the Computational Parameters Section" );
	}
    }


    public void validateMeshContinuationParameters() throws RuntimeException {
	XMLObject computationalParameters;

	if( (computationalParameters = problemTypeSection.getChild("ComputationalParameters" +
								   problemDimensionSuffix )) != null ) {
	    XMLObject meshContinuationSection;
	    if( (meshContinuationSection = computationalParameters.getChild("MeshContinuation" +
									    problemDimensionSuffix )) != null ) {
		numRuns = meshContinuationSection.getInt("nRunsFirst");
		if( numRuns <= 0 ) {
		    throw new RuntimeException("Number of runs in Computational Parameters, Mesh Continuation must be greater than zero.");
		}
		
		if( meshContinuationSection.getAttribute("planeNewNodes").equals("x") ) {
		    planeNewNodes = 0;
		}
		else if( meshContinuationSection.getAttribute("planeNewNodes").equals("y") ) {
		    planeNewNodes = 1;
		}
		else if( meshContinuationSection.getAttribute("planeNewNodes").equals("z") ) {
		    planeNewNodes = 2;
		}
		else {
		    throw new RuntimeException("Invalid plane for new nodes in mesh continuation.");
		}

		if( meshContinuationSection.getAttribute("positionNewNodes").matches("Left.*") ) {
		    positionNewNodes = -1;
		}
		else if( meshContinuationSection.getAttribute("positionNewNodes").matches("Center.*") ) {
		    positionNewNodes = 0;
		}
		else if( meshContinuationSection.getAttribute("positionNewNodes").matches("Right.*") ) {
		    positionNewNodes = 1;
		}
		else {
		    throw new RuntimeException("Invalid position for new nodes in mesh continuation.");
		}

		parameterChange = new double[3];
		parameterChange[0] = meshContinuationSection.getDouble("changeX");
		if( nDim > 1) {
		    parameterChange[1] = meshContinuationSection.getDouble("changeY");
		}
		if( nDim > 2) {
		    parameterChange[2] = meshContinuationSection.getDouble("changeZ");
		}

		guessRange = new double[2];
		guessRange[0] = meshContinuationSection.getDouble("guessRange1");
		guessRange[1] = meshContinuationSection.getDouble("guessRange2");

	    }
	}
    }


    public void validateMixingParameters() throws RuntimeException {
	XMLObject mixingParametersObject;

	if ( (mixingParametersObject = inputMessage.getChild("TramontoMixingParameters")) != null ) {
	    XMLObject mixingType;

	    if ( (mixingType = mixingParametersObject.getChild("LorentzBerthelot")) != null ) {
		// default Lorentz Berthelot mixing.  Just the standard interaction parameters 
		lMix_rule = 0;
 
		for(int child=0; child<mixingType.numChildren(); child++) {

		    XMLObject componentsArray = mixingType.getChild(child);

		    // grab the fluid section
		    if( (componentsArray != null) && 
			(componentsArray.getAttribute("name").equals("fluidMixingArray")) ) {
			
			int nMixingComp = componentsArray.numChildren();
			if( nMixingComp < nComp ) {
			    throw new RuntimeException( "In the Mixing Parameters section, each fluid " + 
							"defined in the Fluid Parameters section " +
							"must have mixing data supplied." );
			    // if nMixingComp > nComp is caught in an exception that looks for
			    // duplicate names
			}
			
			// allocate space for fluid mixing data
			sigma_ff = new double[nComp][nComp];
			eps_ff = new double[nComp][nComp];
			cut_ff = new double[nComp][nComp];
			
			// use this to check if we're setting a single set of values twice
			boolean alreadySet[] = new boolean[nComp];
			
			for(int i=0; i<nMixingComp; i++ ) {
			    
			    // each child has the mixing data as well as one child, with one
			    // child which is the reference to the fluid to which this data applies
			    
			    XMLObject panelOfComponentData = componentsArray.getChild(i);
			    String potentialCompName = panelOfComponentData.getChild(0).getChild(0).getAttribute("name").trim();
			    
			    // check that named component is known and if so get its index
			    int index;
			    try {
				index = findIndexInStringArray( potentialCompName, compName );
			    }
			    catch (RuntimeException exc) {
				throw new RuntimeException( "Component \"" + potentialCompName + 
							    "\" in Mixing Parameters is not in the " +
							    "list of components in the Fluid " +
							    "Parameters section.");
			    }
			    
			    // check if data has already been collectd on this component
			    if( alreadySet[index] == true ) {
				throw new RuntimeException( "Component \"" + potentialCompName + 
							    "\" in Mixing Parameters is listed " +
							    "more than once. ");
			    }
			    
			    /* made it here so component is named and unique */
			    
			    sigma_ff[index][index] = panelOfComponentData.getDouble("sigma_f");
			    eps_ff[index][index] = panelOfComponentData.getDouble("eps_f");
			    cut_ff[index][index] = panelOfComponentData.getDouble("cut_f");
			    alreadySet[index] = true;
			}
			
			// now thow an exception if any of the components was skipped
			for(int i=0; i<nComp; i++) {
			    if( alreadySet[i] != true ) {
				throw new RuntimeException( "Component \"" + compName[i] + 
							    "\" in not listed in Mixing Parameters, " +
							    "fluid mixing.");
			    }
			}
		    } 
		    else if( (componentsArray != null) && 
			(componentsArray.getAttribute("name").equals("surfaceMixingArray")) ) {
			// found the surface mixing data
			
			int nMixingComp = componentsArray.numChildren();
			if( nMixingComp < nWall_type ) {
			    throw new RuntimeException( "In the Mixing Parameters section, each surface " + 
							"defined in the Surface Parameters section " +
							"must have mixing data supplied." );
			    // if nMixingComp > nComp is caught in an exception that looks for
			    // duplicate names
			}
			
			// allocate space for fluid mixing data
			sigma_ww = new double[nWall_type][nWall_type];
			eps_ww = new double[nWall_type][nWall_type];
			cut_ww = new double[nWall_type][nWall_type];
			
			// use this to check if we're setting a single set of values twice
			boolean alreadySet[] = new boolean[nWall_type];
			
			for(int i=0; i<nMixingComp; i++ ) {
			    
			    // each child has the mixing data as well as one child, with one
			    // child which is the reference to the fluid to which this data applies
			    
			    XMLObject panelOfComponentData = componentsArray.getChild(i);
			    String potentialSurfaceName = panelOfComponentData.getChild(0).getChild(0).getAttribute("surfaceName").trim();
			    
			    // check that named surface is known and if so get its index
			    int index;
			    try {
				index = findIndexInStringArray( potentialSurfaceName, surfaceName );
			    }
			    catch (RuntimeException exc) {
				throw new RuntimeException( "Surface \"" + potentialSurfaceName + 
							    "\" in Mixing Parameters is not in the " +
							    "list of surfaces in the Surface " +
							    "Parameters section.");
			    }
			    
			    // check if data has already been collectd on this component
			    if( alreadySet[index] == true ) {
				throw new RuntimeException( "Surface \"" + potentialSurfaceName + 
							    "\" in Mixing Parameters is listed " +
							    "more than once. ");
			    }
			    
			    /* made it here so component is named and unique */
			    
			    sigma_ww[index][index] = panelOfComponentData.getDouble("sigma_w");
			    eps_ww[index][index] = panelOfComponentData.getDouble("eps_w");
			    cut_ww[index][index] = panelOfComponentData.getDouble("cut_w");
			    alreadySet[index] = true;
			}
			
			// now thow an exception if any of the components was skipped
			for(int i=0; i<nWall_type; i++) {
			    if( alreadySet[i] != true ) {
				throw new RuntimeException( "Surface \"" + surfaceName[i] + 
							    "\" in not listed in Mixing Parameters, " +
							    "surface mixing.");
			    }
			}
		    }
		    else {
			throw new RuntimeException( "Could not find fluid or surface mixing data " +
						    "sections in Mixing Parameters section.");
		    }
		}
	    }
	    else if ( (mixingType = mixingParametersObject.getChild("UserDefinedMixing")) != null ) {
		// user is supplying all the parameters so try and read everything
		lMix_rule = 1;
 
		for(int child=0; child<mixingType.numChildren(); child++) {

		    XMLObject componentsArray = mixingType.getChild(child);

		    // grab fluid-fluid data
		    if( (componentsArray != null) && 
			(componentsArray.getAttribute("name").equals("fluidFluidMixingArray")) ) {
			
			int nMixingComp = componentsArray.numChildren();
			if( nMixingComp < (nComp*nComp) ) {
			    throw new RuntimeException( "In the Mixing Parameters section, each pair of fluids " + 
							"defined in the Fluid Parameters section " +
							"must have mixing data supplied." );
			    // if nMixingComp > (nComp*nComp) is caught in an exception that looks for
			    // duplicate names
			}
			
			// allocate space for fluid mixing data
			sigma_ff = new double[nComp][nComp];
			eps_ff = new double[nComp][nComp];
			cut_ff = new double[nComp][nComp];
			
			// use this to check if we're setting a single set of values twice
			boolean alreadySet[][] = new boolean[nComp][nComp];
			
			for(int i=0; i<nMixingComp; i++ ) {
			    
			    // each child has the mixing data as well as two children. Each of these
			    // children is a reference and it's child has the component name.
			    
			    XMLObject panelOfComponentData = componentsArray.getChild(i);
			    String potentialCompNameOne = panelOfComponentData.getChild(0).getChild(0).getAttribute("name").trim();
			    String potentialCompNameTwo = panelOfComponentData.getChild(1).getChild(0).getAttribute("name").trim();
			    
			    // check that named component is known and if so get its index
			    int index1;
			    try {
				index1 = findIndexInStringArray( potentialCompNameOne, compName );
			    }
			    catch (RuntimeException exc) {
				throw new RuntimeException( "Component \"" + potentialCompNameOne + 
							    "\" in Mixing Parameters is not in the " +
							    "list of components in the Fluid " +
							    "Parameters section.");
			    }

			    int index2;
			    try {
				index2 = findIndexInStringArray( potentialCompNameTwo, compName );
			    }
			    catch (RuntimeException exc) {
				throw new RuntimeException( "Component \"" + potentialCompNameTwo + 
							    "\" in Mixing Parameters is not in the " +
							    "list of components in the Fluid " +
							    "Parameters section.");
			    }
			    
			    // check if data has already been collectd on this component
			    if( alreadySet[index1][index2] == true ) {
				throw new RuntimeException( "Component pair \"" + potentialCompNameOne + ", " +
							    potentialCompNameTwo +
							    "\" in Mixing Parameters, fluid mixing is listed " +
							    "more than once. ");
			    }
			    
			    /* made it here so component is named and unique */
			    
			    sigma_ff[index1][index2] = panelOfComponentData.getDouble("sigma_ff");
			    eps_ff[index1][index2] = panelOfComponentData.getDouble("eps_ff");
			    cut_ff[index1][index2] = panelOfComponentData.getDouble("cut_ff");
			    alreadySet[index1][index2] = true;
			}
			
			// now thow an exception if any of the components was skipped

			for(int i=0; i<nComp; i++) {
			    for(int j=0; j<nComp; j++) {
				if( alreadySet[i][j] != true ) {
				    throw new RuntimeException( "Component pair \"" + compName[i] + ", " +
								compName[j] + 
								"\" in not listed in Mixing Parameters, " +
								"fluid mixing.");
				}
			    }
			}
		    }

		    
		    // grab surface-surface data
		    if( (componentsArray != null) && 
			(componentsArray.getAttribute("name").equals("surfaceSurfaceMixingArray")) ) {
			
			int nMixingComp = componentsArray.numChildren();
			if( nMixingComp < (nWall_type*nWall_type) ) {
			    throw new RuntimeException( "In the Mixing Parameters section, each pair of surfaces " + 
							"defined in the Surface Parameters section " +
							"must have mixing data supplied." );
			    // if nMixingComp > (nWall_type*nWall_type) is caught in an exception that looks for
			    // duplicate names
			}
			
			// allocate space for fluid mixing data
			sigma_ww = new double[nWall_type][nWall_type];
			eps_ww = new double[nWall_type][nWall_type];
			cut_ww = new double[nWall_type][nWall_type];
			
			// use this to check if we're setting a single set of values twice
			boolean alreadySet[][] = new boolean[nWall_type][nWall_type];
			
			for(int i=0; i<nMixingComp; i++ ) {
			    
			    // each child has the mixing data as well as two children. Each of these
			    // children is a reference and it's child has the component name.
			    
			    XMLObject panelOfComponentData = componentsArray.getChild(i);
			    String potentialSurfaceNameOne = panelOfComponentData.getChild(0).getChild(0).getAttribute("surfaceName").trim();
			    String potentialSurfaceNameTwo = panelOfComponentData.getChild(1).getChild(0).getAttribute("surfaceName").trim();
			    
			    // check that named surface is known and if so get its index
			    int index1;
			    try {
				index1 = findIndexInStringArray( potentialSurfaceNameOne, surfaceName );
			    }
			    catch (RuntimeException exc) {
				throw new RuntimeException( "Surface \"" + potentialSurfaceNameOne + 
							    "\" in Mixing Parameters is not in the " +
							    "list of surfaces in the Surface " +
							    "Parameters section.");
			    }

			    int index2;
			    try {
				index2 = findIndexInStringArray( potentialSurfaceNameTwo, surfaceName );
			    }
			    catch (RuntimeException exc) {
				throw new RuntimeException( "Surface \"" + potentialSurfaceNameTwo + 
							    "\" in Mixing Parameters is not in the " +
							    "list of surfaces in the Surface " +
							    "Parameters section.");
			    }
			    
			    // check if data has already been collectd on this surface
			    if( alreadySet[index1][index2] == true ) {
				throw new RuntimeException( "Surface pair \"" + potentialSurfaceNameOne + ", " +
							    potentialSurfaceNameTwo +
							    "\" in Mixing Parameters, surface mixing is listed " +
							    "more than once. ");
			    }
			    
			    /* made it here so component is named and unique */
			    
			    sigma_ww[index1][index2] = panelOfComponentData.getDouble("sigma_ww");
			    eps_ww[index1][index2] = panelOfComponentData.getDouble("eps_ww");
			    cut_ww[index1][index2] = panelOfComponentData.getDouble("cut_ww");
			    alreadySet[index1][index2] = true;
			}
			
			// now thow an exception if any of the components was skipped

			for(int i=0; i<nWall_type; i++) {
			    for(int j=0; j<nWall_type; j++) {
				if( alreadySet[i][j] != true ) {
				    throw new RuntimeException( "Surface pair \"" + surfaceName[i] + ", " +
								surfaceName[j] + 
								"\" in not listed in Mixing Parameters, " +
								"surface mixing.");
				}
			    }
			}
		    }
		    
		    // grab surface-fluid data 


		}

	    }
	    else {
		throw new RuntimeException( "Could not find the TramontoMixingParameters, mixing type section");
	    }
	}
	else {
	    throw new RuntimeException( "Could not find the TramontoMixingParameters section");
	}

    }


    public void validateNonLinearSolverParameters() throws RuntimeException {
	XMLObject computationalParameters;

	if( (computationalParameters = problemTypeSection.getChild("ComputationalParameters" +
							     problemDimensionSuffix )) != null ) {
	    XMLObject nonLinearSolverSection;
	    if( (nonLinearSolverSection = computationalParameters.getChild("NonlinearSolver")) != null ) {
		maxNewtonIter = nonLinearSolverSection.getInt("maxNewtonIterations");
		newtonRelTol =  nonLinearSolverSection.getDouble("relativeConvergence");
		newtonAbsTol =  nonLinearSolverSection.getDouble("absoluteConvergence");
		if( nonLinearSolverSection.getAttribute("loadBalanceSwitch").equals("Linear") ) {
		    loadBalanceSwitch = 0;
		}
		else if( nonLinearSolverSection.getAttribute("loadBalanceSwitch").equals("Box") ) {
		    loadBalanceSwitch = 1;
		}
		else if( nonLinearSolverSection.getAttribute("loadBalanceSwitch").equals("Weights") ) {
		    loadBalanceSwitch = 2;
		}
		else if( nonLinearSolverSection.getAttribute("loadBalanceSwitch").equals("Timings") ) {
		    loadBalanceSwitch = 3;
		}
		else {
		    throw new RuntimeException("Unknown load balancing switch in nonlinear solver, Computational Parameters Section.");
		}
	    }
	    else {
		throw new RuntimeException("No non-linear solver parameters found in Computational Parameters Section.");
	    }
	}
	else {
	    throw new RuntimeException("No Computational Parameters Section found.");
	}
    }


    public void validateOutputFormatParameters() throws RuntimeException {
	XMLObject computationalSection;
	if( (computationalSection = problemTypeSection.getChild( "ComputationalParameters" +
								 problemDimensionSuffix ) ) != null ) {
	    XMLObject outputFormatSection;
	    if( (outputFormatSection = computationalSection.getChild("OutputFormat")) != null ) {
		
		// lPerArea
		if( outputFormatSection.getAttribute("areaScaling").equals("true") ) {
		    lPerArea = 1;
		}
		else {
		    lPerArea = 0;
		}

		// lPrintGr
		if( outputFormatSection.getAttribute("printGr").equals("true") ) {
		    lPrintGr = 1;
		}
		else {
		    lPrintGr = 0;
		}
		
		// printRhoType
		if( outputFormatSection.getAttribute("densityType").equals("All output in default files (dft_dens.dat, dft_dens2.dat or dft_dens.datg)" ) ) {
		    printRhoType = 0;
		}
		else if( outputFormatSection.getAttribute("densityType").equals("Output from each run in a different file (dft_dens0.0, dft_dens1.0...)" ) ) {
		    printRhoType = 1;
		}
		else {
		    throw new RuntimeException("Unknown density type in output format section,  Computational Parameters Section.");
		}

		// printDensitSwitch
		if( outputFormatSection.getAttribute("densitySwitch").equals("Only densities (rho_b sigma^3)" ) ) {
		    printRhoSwitch = 0;
		}
		else if( outputFormatSection.getAttribute("densitySwitch").equals("Pressure (p / p_sat only valid for a 1 component LJ fluid)" ) ) {
		    printRhoSwitch = 1;
		}
		else if( outputFormatSection.getAttribute("densitySwitch").equals("Debye screening length (1/k only valid for electrolytes)" ) ) {
		    printRhoSwitch = 2;
		}
		else if( outputFormatSection.getAttribute("densitySwitch").equals("Chemical potential (mu / kT)" ) ) {
		    printRhoSwitch = 3;
		}
		else {
		    throw new RuntimeException("Unknown density switch in output format section, Computational Parameters Section.");
		}

		// printMeshSwitch
		if( outputFormatSection.getAttribute("meshSwitch").equals("Print surface separations between all pairs of surfaces" ) ) {
		    printMeshSwitch = 0;
		}
		else if(outputFormatSection.getAttribute("meshSwitch").equals("Print surface positions at the center of each surface" ) ) {
		    printMeshSwitch = 1;
		}
		else {
		    throw new RuntimeException("Unknown mesh switch in output format section, Computational Parameters Section");
		}

		// printer header flag
		if( outputFormatSection.getAttribute("printHeader").equals("true" ) ) {
		    printHeader = 1;
		}
		else {
		    printHeader = 0;
		}

		// output quantity printIWrite
		if( outputFormatSection.getAttribute("outputQuantity").equals("Minimal output, no density profiles" ) ) {
		    printIWrite = 0;
		}
		else if( outputFormatSection.getAttribute("outputQuantity").equals("Minimal output and density profiles" ) ) {
		    printIWrite = 1;
		}
		else if( outputFormatSection.getAttribute("outputQuantity").equals("Verbose" ) ) {
		    printIWrite = 3;
		}
		else {
		    throw new RuntimeException("Unknown output quantity choice in output format section, Computational Parameters Section.");
		}

	    }
	    else {
		throw new RuntimeException("No output format parameters section found, Computational Parameters Section.");
	    }
	}
	else {
	    throw new RuntimeException("No computational Parameters Section found.");
	}
    }

    public void validatePolymerInputParameters() throws RuntimeException {
	XMLObject polymerParametersObject;
	nPol_comp = 0;
	polymerVector = new Vector();
	boolean[] polymerDataFound = new boolean[ nComp ];

	// grab the Fluid parameters section of the inputmessage
	if( (polymerParametersObject = problemTypeSection.getChild("PolymericFluidType")) != null ) {
	    cr_file = polymerParametersObject.getAttribute("cRFile");
	    XMLObject componentsArray = polymerParametersObject.getChild("Array");
	    if( (componentsArray != null) && 
		(componentsArray.getAttribute("name").equals("componentsArray")) ) {
		// loop over each component
		for( int i = 0; i< componentsArray.numChildren(); i++) {
		    XMLObject componentPanel = componentsArray.getChild(i);
		    // look for the reference for a polymer name
		    XMLObject polymerReference = componentPanel.getChild("Reference");
		    String polymerName;
		    int polymerIndex;
		    if( polymerReference != null ) {
			XMLObject polymerNamePanel = polymerReference.getChild(0);
			if( polymerNamePanel != null ) {
			    polymerName = polymerNamePanel.getAttribute("name");
			    try {
				polymerIndex = findIndexInStringArray( polymerName, compName );
			    }
			    catch (RuntimeException e) {
				throw new RuntimeException( "Could not find polymer name, \"" + polymerName +
							    "\" in list of component names.");
			    }
			}
			else {
			    throw new RuntimeException("Could not find polymer name panel.");
			}
		    }
		    else {
			throw new RuntimeException("Could not find polymer reference section.");
		    }
		    // now look for a valid polymer section
		    XMLObject polymerSection = componentPanel;
		    if( polymerSection != null ) {
			// found a polymer section so now grab its data
			XMLObject blockArray = polymerSection.getChild("Array");
			int numBlock=0;
			if( blockArray != null ) {
			    numBlock = blockArray.numChildren();
			}
			polymerBlockInfo thePolymer = new polymerBlockInfo(numBlock);
			thePolymer.poly_file = polymerSection.getAttribute("polyFile");
			thePolymer.cr_radius = polymerSection.getDouble("cRRadius");
			thePolymer.gauss_a = polymerSection.getDouble("gaussA");

			if( (blockArray != null ) &&
			    (blockArray.getAttribute("name").equals("polymerSegmentArray")) ) {
			    for( int j=0; j<blockArray.numChildren(); j++) {
				XMLObject blockData = blockArray.getChild(j);
				int blockNumber = blockData.getInt("blockNumber");
				if( (blockNumber < 0) || (blockNumber > (numBlock-1)) ) {
				    throw new RuntimeException( "Polymer block number, \"" + blockNumber + ",\"" +
								" in fluid, \"" +
								componentPanel.getAttribute("name") + "\"" +
								" in polymer parameters section is out of range 0 to " + 
								(numBlock-1) + "(numBlock-1)" + "in Fluid Parameters Section.");
				}
				thePolymer.block[blockNumber] = blockData.getInt("segments");
				thePolymer.block_type[blockNumber] = blockData.getInt("segmentType");
			    }
			}
			polymerVector.add( polymerIndex, thePolymer );
			if( polymerDataFound[ polymerIndex ] == true ) {
			    throw new RuntimeException("Data for polymer, \"" + compName[polymerIndex] +
						       "\" was found more than once in the Polymer types section.");
			}
			polymerDataFound[ polymerIndex ] = true;
			nPol_comp++;
		    }
		}
	    }
	    else {
		throw new RuntimeException( "No fluid components found in the Fluid Parameters Section." );
	    }

	    // check that everything was entered.
	    for( int i=0; i<nComp; i++) {
		if( polymerDataFound[i] == false ) {
		    throw new RuntimeException("No polymer data was found for polymer component, \"" +
					       compName[i] + "\".");
		}
	    }
	}
	else {
	    if( problemTypeSection.getTag().matches(".*Polymeric.*") ) {
		throw new RuntimeException( "No Fluid Parameters Section found." );
	    }
	    // not an error so do nothing
	}

    }

    public void validatePotentialTypeParameters() throws RuntimeException {
	XMLObject potentialData;

	String classNames[] = { "NoWallSurfaceInteractions1D",             "0",
				"NoWallSurfaceInteractions2D3D",           "0",
				"NoWallSurfaceInteractions1DPolymer",             "0",
				"NoWallSurfaceInteractions2D3DPolymer",           "0",
				"HardWallSurfaceInteractions1D",           "1",
				"HardWallSurfaceInteractions2D3D",         "1",
				"HardWallSurfaceInteractions1DPolymer",           "1",
				"HardWallSurfaceInteractions2D3DPolymer",         "1",
				"LJ93SurfaceInteractions1D",               "2",
				"LJ93SurfaceInteractions2D3D",             "2",
				"LJ93SurfaceInteractions1DPolymer",               "2",
				"LJ93SurfaceInteractions2D3DPolymer",             "2",
				"NormLJ93SurfaceInteractions1D",           "3",
				"NormLJ93SurfaceInteractions2D3D",         "3",
				"NormLJ93SurfaceInteractions1DPolymer",           "3",
				"NormLJ93SurfaceInteractions2D3DPolymer",         "3",
				"LJ126SurfaceInteractions1D",              "4",
				"LJ126SurfaceInteractions2D3D",            "4",
				"LJ126SurfaceInteractions1DPolymer",              "4",
				"LJ126SurfaceInteractions2D3DPolymer",            "4",
				"LJAtomicSurfaceInteractions1D",           "5",
				"LJAtomicSurfaceInteractions2D3D",         "5",
				"LJAtomicSurfaceInteractions1DPolymer",           "5",
				"LJAtomicSurfaceInteractions2D3DPolymer",         "5",
                                "LJStepSurfaceInteractions2D3D",           "6",
                                "LJStepSurfaceInteractions2D3DPolymer",           "6",
				"HardExpSurfaceInteractions1D",            "7",
				"HardExpSurfaceInteractions2D3D",          "7",
				"HardExpSurfaceInteractions1DPolymer",            "7",
				"HardExpSurfaceInteractions2D3DPolymer",          "7",
				"HardChargedAtomsSurfaceInteractions1D",   "8",
				"HardChargedAtomsSurfaceInteractions2D3D", "8",
				"HardChargedAtomsSurfaceInteractions1DPolymer",   "8",
				"HardChargedAtomsSurfaceInteractions2D3DPolymer", "8",
				"ChargedLJAtomsSurfaceInteractions1D",     "9",
				"ChargedLJAtomsSurfaceInteractions2D3D",   "9",
				"ChargedLJAtomsSurfaceInteractions1DPolymer",     "9",
				"ChargedLJAtomsSurfaceInteractions2D3DPolymer",   "9" };


	boolean foundData = false;

	for( int i=0; i < classNames.length; i+=2 )
	    {
		if( (potentialData = problemTypeSection.getChild( classNames[i] )) != null ) {
		    // fouond a matching class,
		    // so set the type_ variable flags and exit this loop
		    foundData = true;
		    ipot_wf_n = Integer.parseInt( classNames[++i] );
		    maximumPotential=potentialData.getDouble("maxPotential");
		    break;
		}
	    }
	
	if( foundData == false ) {
	    throw new RuntimeException( "No valid surface interactions selection found in the Problem Type section." );
	}
    }

    public void	validateProblemType() throws RuntimeException {
	if( (problemTypeSection = inputMessage.getChild("OneDAtomicProblem")) != null ) {
	    problemDimensionPrefix="OneD";
	    problemDimensionSuffix="1D";
	    polymericSuffix="";
	    problemDimension=1;
	    return;
	}
	else if( (problemTypeSection = inputMessage.getChild("TwoDAtomicProblem")) != null ) {
	    problemDimensionPrefix="TwoD";
	    problemDimensionSuffix="2D";
	    polymericSuffix="";
	    problemDimension=2;
	    return;
	}
	else if( (problemTypeSection = inputMessage.getChild("ThreeDAtomicProblem")) != null ) {
	    problemDimensionPrefix="ThreeD";
	    problemDimensionSuffix="3D";
	    polymericSuffix="";
	    problemDimension=3;
	    return;
	}
	else if( (problemTypeSection = inputMessage.getChild("OneDPolymericProblem")) != null ) {
	    problemDimensionPrefix="OneD";
	    problemDimensionSuffix="1D";
	    polymericSuffix="Polymeric";
	    problemDimension=1;
	    return;
	}
	else if( (problemTypeSection = inputMessage.getChild("TwoDPolymericProblem")) != null ) {
	    problemDimensionPrefix="TwoD";
	    problemDimensionSuffix="2D";
	    polymericSuffix="Polymeric";
	    problemDimension=2;
	    return;
	}
	else if( (problemTypeSection = inputMessage.getChild("ThreeDPolymericProblem")) != null ) {
	    problemDimensionPrefix="ThreeD";
	    problemDimensionSuffix="3D";
	    polymericSuffix="Polymeric";
	    problemDimension=3;
	    return;
	}
	else {
	    throw new RuntimeException( "Could not determine the problem type (i.e. 1D,2D or 3D & Atomic or Polymeric)");
	}
    }

    public void validateSemiPermeableSurfaceParameters() throws RuntimeException {
        XMLObject boundaryConditionsSection;

	// allocate the permeable data arrays
	lSemiPerm = new int[nWall_type][nComp];
	vExtMembrane = new double[nWall_type][nComp];

	// find the surface parameters sections
	if( (boundaryConditionsSection = problemTypeSection.getChild("BoundaryConditions" + problemDimensionSuffix)) != null ) {
	    XMLObject semiPermeableSurfacesSection;
	    if( problemTypeSection.getTag().matches(".*Atomic.*") ) {
		semiPermeableSurfacesSection = boundaryConditionsSection.getChild("SemiPermeableSurfacesAtomic");
	    }
	    else {
		semiPermeableSurfacesSection = boundaryConditionsSection.getChild("SemiPermeableSurfacesPolymeric");
	    }
	    if( semiPermeableSurfacesSection != null ) {

	    }
	    else {
		throw new RuntimeException("Couldn't find SemiPermeable Surfaces section.");
	    }
		
	    XMLObject surfTable = semiPermeableSurfacesSection.getChild("Table");
	    if( surfTable != null ) {
		for( int entryNumber=0; entryNumber < surfTable.numChildren(); entryNumber++ ) {
		    XMLObject anEntry = surfTable.getChild(entryNumber);
		    // look at it's references to find out the fluid/surface this applies to
		    XMLObject reference1, reference2;
		    String surface1Name, component1Name;
		    int surface1Index, component1Index;
		    reference1 = anEntry.getChild(0).getChild(0);
		    reference2 = anEntry.getChild(1).getChild(0);
		    if( reference1.hasAttribute("surfaceName") ) {
			surface1Name = reference1.getAttribute("surfaceName");
			component1Name = reference2.getAttribute("name");
		    }
		    else {
			surface1Name = reference2.getAttribute("surfaceName");
			component1Name = reference1.getAttribute("name");
		    }
		    try {
			surface1Index = findIndexInStringArray( surface1Name, surfaceName );
			component1Index = findIndexInStringArray( component1Name, compName );
		    }
		    catch (RuntimeException e) {
			throw new RuntimeException("Couldn't find surface, \"" + surface1Name + "\" and " +
						   "fluid name, \"" + component1Name + "\" in list of defined names.");
		    }
		    
		    lSemiPerm[surface1Index][component1Index] = 1;
		    vExtMembrane[surface1Index][component1Index] = anEntry.getDouble("vExtMembrane");
		}
	    }
	    else {
		throw new RuntimeException("Couldn't find SemiPermeable Surfaces table.");
	    }
	}
    }

    public void validateStartupControlParameters() throws RuntimeException {
	// initial guess types and their values
	String[] guessTypes = { "GuessConstBulk", "-3",
				"GuessConstLiquid", "-2",
				"GuessConstVapor", "-1",
				"GuessIdealBulk", "0",
				"GuessIdealLiquid", "1",
				"GuessIdealVapor", "2",
				"GuessSteppedBulk", "3",
				"GuessSteppedLiquid", "4",
				"GuessSteppedVapor", "5",
				"GuessSteppedLiquidVapor", "6",
				"GuessLinearLiquidVapor", "7" };

	XMLObject computationalParameters = problemTypeSection.getChild("ComputationalParameters" +
									problemDimensionSuffix );

	if( computationalParameters != null ) {
	    XMLObject startupControlSection;
	    if( (startupControlSection = computationalParameters.getChild("StartupControl")) != null ) {
		String liquidVaporControl = startupControlSection.getAttribute("liquidVapor");
		if( liquidVaporControl == null ) {
		    throw new RuntimeException("No liquid-vapor coexistence parameter.");
		}
		else if( liquidVaporControl.equals("Off and don't calculate coexistence properties") ) {
		    iLiqVap = -2;
		}
		else if( liquidVaporControl.equals("Off but calculate bulk coexistence properties") ) {
		    iLiqVap = -1;
		}
		else if( liquidVaporControl.equals("Wall-vapor") ) {
		    iLiqVap = 1;
		}
		else if( liquidVaporControl.equals("Wall-liquid") ) {
		    iLiqVap = 2;
		}
		else if( liquidVaporControl.equals("Liquid-vapor") ) {
		    iLiqVap = 3;
		}
		else {
		    throw new RuntimeException("Invalid liquid-vapor coexistence setting, \"" + liquidVaporControl + "\"");
		}
		
		// look at the type of child to figure out the rest of the startup control parameters
		XMLObject restartObject = startupControlSection.getChild(0);
		if( restartObject == null ) {
		    throw new RuntimeException("No valid restart object found in startup control.");
		}
		else if( restartObject.getTag().equals("RestartNoRestartFile") ) {
		    restart = 0;
		    // look for initial guess type
		    XMLObject guessType = restartObject.getChild(0);
		    if( guessType != null ) {
			int guessIndex;
			try {
			    guessIndex = findIndexInStringArray( guessType.getTag(), guessTypes );
			}
			catch (RuntimeException e) {
			    throw new RuntimeException("Unknown guess type in startup control.");
			}
			iGuess = Integer.parseInt( guessTypes[ (guessIndex + 1) ] );

			// some guess types require an extra parameter, so read that.
			if( iGuess >= 3 ) {
			    thickness = guessType.getDouble("thickness");
			}
			else {
			    thickness = 0;
			}

		    }
		    else {
			throw new RuntimeException("Couldn't read guess type in startup control.");
		    }
		}
		else if( restartObject.getTag().equals("RestartFromFile") ) {
		    restart = 1;
		}
		else if( restartObject.getTag().equals("RestartFromFileWithThickness") ) {
		    restart = 2;
		    thickness = restartObject.getDouble("thickness");
		    String guessType = restartObject.getAttribute("startGuess1");
		    if( guessType == null ) {
			throw new RuntimeException("No startup guess when restarting from file with thickness");
		    }
		    else if( guessType.equals("Stepped (bulk)") ) {
			iGuess = 3;
		    }
		    else if( guessType.equals("Stepped (liquid-coexistence)") ) {
			iGuess = 4;
		    }
		    else if( guessType.equals("Stepped (vapor-coexistence)") ) {
			iGuess = 5;
		    }
		    else if( guessType.equals("Stepped (between liquid and vapor)")) {
			iGuess = 6;
		    }
		    else {
			throw new RuntimeException("Invalid guess when restarting from a file with thickness, \"" +
						   guessType + "\"");
		    }
		    return;
		}
		else if( restartObject.getTag().equals("RestartFromFileWithOldDensity") ) {
		    restart = 3;
		}
		else {
		    throw new RuntimeException("Invalid restart type in startup control");
		}
		// look for rho max
		rhoMax = startupControlSection.getDouble("rhoMax");

	    }
	    else {
		throw new RuntimeException("Could not find a valid startup control section");
	    }
	}
	else {
	    throw new RuntimeException("Could not find a valid computational parameters section");
	}
    }


    public void validateStatePointParameters() throws RuntimeException {
	bulkDensity = new double[ nComp ];

	boolean[] dataFound = new boolean[nComp];
	XMLObject boundaryConditionsSection = problemTypeSection.getChild( "BoundaryConditions" + 
									   problemDimensionSuffix );
	if( boundaryConditionsSection != null ) {
	    XMLObject statePointParametersSection = boundaryConditionsSection.getChild("EquilibriumBoundaryConditions" + 
										       problemDimensionSuffix );
	    if( statePointParametersSection != null ) {
		XMLObject aTable = statePointParametersSection.getChild("Table");
		if( aTable != null ) {
		    for( int entryNum = 0 ; entryNum < aTable.numChildren(); entryNum++) {
			XMLObject anEntry = aTable.getChild(entryNum);
			XMLObject aRef = anEntry.getChild("Reference");
			String componentName;
			int componentIndex;
			if( aRef != null ) {
			    XMLObject refPanel = aRef.getChild(0);
			    if( refPanel != null ) {
				componentName = refPanel.getAttribute("name");
				try {
				    componentIndex = findIndexInStringArray( componentName, compName );
				}
				catch (RuntimeException e) {
				    throw new RuntimeException("Couldn't find component, \"" + componentName + 
							       "\" in list of component names");
				}
				if( dataFound[ componentIndex ] == true ) {
				    throw new RuntimeException("Bulk density data for component, \"" +
							       componentName + "\" appears more than once " +
							       "in the Boundary Conditions section");
				}
				dataFound[ componentIndex ] = true;
				bulkDensity[ componentIndex ] = anEntry.getDouble("rhoB");
			    }
			    else {
				throw new RuntimeException("Couldn't read reference panel in equilibrium boundary conditions table.");
			    }   

			}
			else {
			    throw new RuntimeException("Couldn't read reference in equilibrium boundary conditions table.");
			}
		    }
		    // check that everyone had a bulk density read
		    for(int i=0; i<nComp; i++) {
			if( dataFound[i] == false ) {
			    throw new RuntimeException("Data for component, \"" + compName[i] + "\" was not found " +
						       "in the Boundary Conditions list of bulk densities.");
			}
		    }
		}
		else {
		    throw new RuntimeException("Could not find the table in the equilibrium boundary conditions section.");
		}
	    }
	}
	else {
	    throw new RuntimeException("Could not find the boundary conditions section");
	}
    }


    public void validateSteadyStateBoundaryConditionParameters() throws RuntimeException {
	// assume no steady state
	lSteadyState = 0;

	// now look for the steady state section
	XMLObject boundaryConditions = problemTypeSection.getChild("BoundaryConditions" + problemDimensionSuffix );
	if( boundaryConditions != null ) {
	    
	    XMLObject transportSection = boundaryConditions.getChild("TransportBoundaryConditions" + problemDimensionSuffix );
	    if( transportSection != null ) {

		lSteadyState = 1;
		bulkDensityLBB = new double[nComp];
		bulkDensityRTF = new double[nComp];
		diffusionCoef = new double[nComp];

		if( transportSection.getAttribute("linearTransport").equals("Linear") ) {
		    linearTransport = 1;
		}
		else if( transportSection.getAttribute("linearTransport").equals("Non-Linear") ) {
		    linearTransport = 0;
		}
		else {
		    throw new RuntimeException("Bad transport type, \"" + 
					       transportSection.getAttribute("linearTransport") + "\"");
		}

		if( transportSection.getAttribute("gradientDirection").equals("x") ) {
		    gradientDirection = 0;
		}
		else if(transportSection.getAttribute("gradientDirection").equals("y") ) {
		    gradientDirection = 1;
		}
		else if(transportSection.getAttribute("gradientDirection").equals("z") ) {
		    gradientDirection = 2;
		}
		else {
		    throw new RuntimeException("Bad gradient direction specified, \"" +
					       transportSection.getAttribute("gradientDirection") + "\"");
		}

		xConstMu = transportSection.getDouble("xConstantMu");
		elecPotLBB = transportSection.getDouble("electricPotentialLBB");
		elecPotRTF = transportSection.getDouble("electricPotentialRTF");
		centerOfMassVelocity = transportSection.getDouble("velocity");
		
		// now read the table of bulk densities
		boolean[] dataFound = new boolean[ nComp ];
		XMLObject aTable = transportSection.getChild("Table");
		if( (aTable != null) &&
		    (aTable.getAttribute("name").equals("fluidBulkDensities")) ) {
		    for( int entryNum=0; entryNum < aTable.numChildren(); entryNum++) {
			XMLObject anEntry = aTable.getChild( entryNum );
			String aComponentName;
			int aComponentIndex;
			if( anEntry != null ) {
			    XMLObject aRef = anEntry.getChild("Reference");
			    if( aRef != null ) {
				XMLObject aPanel = aRef.getChild("PanelOfComponents");
				if( aPanel != null ) {
				    aComponentName = aPanel.getAttribute("name");
				    try {
					aComponentIndex = findIndexInStringArray( aComponentName, compName );
				    }
				    catch (RuntimeException e ) {
					throw new RuntimeException( "Could not find component, \"" + aComponentName +
								    "\" in list of components in Transport Boundary Conditions.");
				    }
				}
				else {
				    throw new RuntimeException( "Could not find component panel in Transport Boundary Conditions.");
				}
			    }
			    else {
				throw new RuntimeException("Could not find componet reference in Transport Boundary Conditions.");
			    }
			}
			else {
			    throw new RuntimeException("Error in findind a table entry in Transport Boundary Conditions.");
			}
			// Check if this data has already been set
			if( dataFound[ aComponentIndex ] == true ) {
			    throw new RuntimeException("Data for component \"" + aComponentName + "\" appears more than once " +
						       "in the bulk density table for transport boundary conditions.");
			}
			else {
			    dataFound[ aComponentIndex ] = true;
			}
			// ok we have a component index, now attach the data to it
			bulkDensityLBB[ aComponentIndex ] = anEntry.getDouble("rhoB1");
			bulkDensityRTF[ aComponentIndex ] = anEntry.getDouble("rhoB2");
			diffusionCoef[ aComponentIndex ] = anEntry.getDouble("diffCoef");
		    }
		    // check that all bulk densities have been read
		    for(int i=0; i<nComp; i++ ) {
			if( dataFound[i] != true ) {
			    throw new RuntimeException("Data for component \"" + compName[i] + "\" is missing " +
						       "in the bulk density table for transport boundary conditions." +
						       "  Each component must appear once.");
			}
		    }

		    // finally look for data on a 1d pore
		    XMLObject poreData = transportSection.getChild("NotAOneDPoreGeometry");
		    if( poreData != null ) {
			// not much to do here
			geometryFlag = 0;
			n1DPoreSegments = 0;
		    }
		    poreData = transportSection.getChild("UnitAreaOneDPoreGeometry");
		    if( poreData != null ) {
			geometryFlag = 0;
			n1DPoreSegments=1;
			pore1DRadiusLeft = new double[ n1DPoreSegments ];
			pore1DRadiusRight = new double[ n1DPoreSegments ];
			pore1DLength = new double[ n1DPoreSegments ];
			pore1DRadiusLeft[0] = poreData.getDouble("radius");
			pore1DRadiusRight[0] = pore1DRadiusLeft[0];
			pore1DLength[0] = poreData.getDouble("length");
		    }		    
		    poreData = transportSection.getChild("CylindricalPoreOneDPoreGeometry");
		    if( poreData != null ) {
			geometryFlag = 1;
			n1DPoreSegments=1;
			pore1DRadiusLeft = new double[ n1DPoreSegments ];
			pore1DRadiusRight = new double[ n1DPoreSegments ];
			pore1DLength = new double[ n1DPoreSegments ];
			pore1DRadiusLeft[0] = poreData.getDouble("radius");
			pore1DRadiusRight[0] = pore1DRadiusLeft[0];
			pore1DLength[0] = poreData.getDouble("length");
		    }
		    poreData = transportSection.getChild("VariablePoreOneDPoreGeometry");
		    if( poreData != null ) {
			geometryFlag = 2;
			XMLObject poreTable = poreData.getChild("Table");
			if( poreTable != null ) {
			    n1DPoreSegments=poreTable.numChildren();			    
			    pore1DRadiusLeft = new double[ n1DPoreSegments ];
			    pore1DRadiusRight = new double[ n1DPoreSegments ];
			    pore1DLength = new double[ n1DPoreSegments ];
			    boolean entryFound[] = new boolean[ n1DPoreSegments ];
			    for(int entryNum = 0; entryNum < poreTable.numChildren(); entryNum++ ) {
				XMLObject anEntry = poreTable.getChild( entryNum );
				if( anEntry != null ) {
				    int segIndex = anEntry.getInt("segmentNumber");
				    if( (segIndex >= 0) && (segIndex < n1DPoreSegments) ) {
					if( entryFound[ segIndex ] != true ) {
					    pore1DRadiusLeft[segIndex] = anEntry.getDouble("radiusLeft");
					    pore1DRadiusRight[segIndex] = anEntry.getDouble("radiusRight");
					    pore1DLength[segIndex] = anEntry.getDouble("length");
					    entryFound[ segIndex ] = true;
					}
					else {
					    throw new RuntimeException("In the Transport boundary conditions, variable pore size " +
								       "table, the segment number \"" + segIndex + "\" appers " +
								       "more than once.  Segement numbers must run from 0 to " +
								       "the total number of segments - 1.");
					}
				    }
				    else {
					throw new RuntimeException("Transport boundary conditions, variable 1D pore, " +
								   "segment number must greater than or equal to 0, " +
								   "and less than the total number of segments, \"" +
								   n1DPoreSegments + "\"");
				    }
				}
				else {
				    throw new RuntimeException("Error in finding an entry in the variable 1D pore table.");
				}
			    }
			}
			else {
			    throw new RuntimeException("Couldn't find table of pore segments in variable 1D pore description.");
			}
		    }
		}
		else {
		    throw new RuntimeException("Couldn't find table of bulk densities under Transport Boundary Conditions.");
		}
	    }
	    else {
		// not steady state.  dimension arrays to 1 so that we'll print zeros in the input file
		pore1DRadiusLeft = new double[1];
		pore1DRadiusRight = new double[1];
		pore1DLength = new double[1];
		// dimension these arrays properly for safe input file generation.
		bulkDensityLBB = new double[nComp];
		bulkDensityRTF = new double[nComp];
		diffusionCoef = new double[nComp];
	    }
	}
	else {
	    throw new RuntimeException("Could not find the Boundary Conditions section");
	}
    }


    public void validateSurfaceNames() throws RuntimeException {
	XMLObject surfaceNamesSection;

	if( (surfaceNamesSection = inputMessage.getChild("SurfaceNames")) != null ) {
	    XMLObject surfacesArray = surfaceNamesSection.getChild("Array");
	    if( (surfacesArray != null) && 
		(surfacesArray.getAttribute("name").equals("surfaceArray")) ) {
		nWall_type = surfacesArray.numChildren();
		surfaceName = new String[nWall_type];

		for(int i=0; i<nWall_type; i++ ) {
		    XMLObject panelOfSurfaceData = surfacesArray.getChild(i);
		    /* need to check that name is not blank and unique */
		    String potentialSurfaceName = panelOfSurfaceData.getAttribute("surfaceName").trim();
		    if( potentialSurfaceName.equals("") ) {
			throw new RuntimeException( "One of the surfaces has a blank name. " +
						    "All surfaces must have unique names." );
		    }
		    if( isStringInStringArray(  potentialSurfaceName, surfaceName )  == true ) {
			throw new RuntimeException( "Surface \"" + potentialSurfaceName + 
						    "\" is duplicated in list of surfaces. " +
						    "All surfaces must have unique names.");
		    }

		    /* made it here so component is named and unique */
		    surfaceName[i] = panelOfSurfaceData.getAttribute("surfaceName").trim();
		}
	    }
	    else {
		throw new RuntimeException( "No surfaces found in Surface Name Section" );
	    }
	}	
	else {
	    throw new RuntimeException( "No valid Surface names section found." );
	}
    }


    public void validateSurfaceParameters() throws RuntimeException {
	XMLObject surfaceTypesSection;

	if( (surfaceTypesSection = problemTypeSection.getChild("SurfaceEditor" + problemDimensionSuffix )) != null ) {
	    XMLObject surfaceArray;

	    // read in the total number of surfaces
	    nWall = surfaceTypesSection.getInt("nWall");

	    // read the autocenter flag
	    if( surfaceTypesSection.getBoolean("autoCenter") == true ) {
		autoCenterSurfInBox = 1;
	    }
	    else {
		autoCenterSurfInBox = 0;
	    }

	    if( (surfaceArray = surfaceTypesSection.getChild("Array")) != null ) {

		// each memeber of the surfaceArray is a panel which has two children
		// one child is a referece to the surface name and the other child
		// is the class which describes this surface

		// a quick test to see if we at least have the same number of surface names and surface descriptions
		if( surfaceArray.numChildren() != nWall_type ) {
		    throw new RuntimeException("The number of surfaces defined in the Surface Names section, " +
					       nWall_type +
					       ", is different than the number of surfaces described in the Surface Types section, " +
					       surfaceArray.numChildren() +
					       "." );
		}
		if( nWall_type > nWall ) {
		    throw new RuntimeException( "The total number of surfaces, including duplicates, " +
						nWall +
						", is less than the number of surface types, " +
						nWall_type +
						", as described in the Surface types section");
		}
		boolean surfFound[] = new boolean[nWall_type];
		surf_type = new int[nWall_type];
		orientation = new int[nWall_type];
		wallParam = new double[nWall_type];
		wallParam2 = new double[nWall_type];
		wallParam3 = new double[nWall_type];

		// now loop over all the panels to find each surface description
		for(int i=0; i<surfaceArray.numChildren(); i++) {
		    XMLObject surfacePanel = surfaceArray.getChild(i);

		    // first find the reference object in this panel
		    XMLObject referenceChild;
		    int surfaceIndex;
		    if( (referenceChild = surfacePanel.getChild("Reference")) != null ) {
			String referencedSurface = referenceChild.getAttribute("selection");
			try {
			    surfaceIndex = findIndexInStringArray( referencedSurface, surfaceName );
			}
			catch( RuntimeException e ) {
			    throw new RuntimeException( "Couldn't find surface named \"" +
						    referencedSurface +
						    "\" in the set of surfaces defined in the Surface Names section." );
			}
			if( surfFound[ surfaceIndex ] == true ) {
			    // this surface has already been defined 
			    throw new RuntimeException( "The surface \"" +
						    referencedSurface +
						    "\" has been defined more than once in the Surface Types section.");
			}

		    }
		    else {
			throw new RuntimeException( "Could not find the reference for the " +
						    i +
						    "'th panel in the Surface Types section.");
		    }
		    // remember that we found it;
		    surfFound[ surfaceIndex ] = true;

		    // now scan the children of the Panel to find the
		    // type of surface being described.
		    for( int j=0; j<surfacePanel.numChildren(); j++) {
			XMLObject surfaceObject = surfacePanel.getChild(j);

			boolean foundSurfaceType = false;

			// here we begin an exhaustive search over all of the possible surface classes
			if( surfaceObject.getTag().equals("InfinitePlanarWall") ){
			    surf_type[ surfaceIndex ] = 0;
			    orientation[ surfaceIndex ] = 0;
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("halfThickness");
			    wallParam2[ surfaceIndex ] = 0.0;
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			else if( (surfaceObject.getTag().equals("InfinitePlanarWall2D")) || 
				 (surfaceObject.getTag().equals("InfinitePlanarWall3D")) ) {
			    surf_type[ surfaceIndex ] = 0;
			    orientation[ surfaceIndex ] = mapAttributeToInteger(
										surfaceObject.getAttribute("surfaceNormal"), 
										new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("halfThickness");
			    wallParam2[ surfaceIndex ] = 0.0;
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			else if( surfaceObject.getTag().equals("FinitePlanarWall") ){
			    surf_type[ surfaceIndex ] = 1;
			    orientation[ surfaceIndex ] = 0;
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("halfThickness1");
			    wallParam2[ surfaceIndex ] = surfaceObject.getDouble("halfThickness2");
			    wallParam3[ surfaceIndex ] = surfaceObject.getDouble("halfThickness3");
			    foundSurfaceType = true;
			}
			else if( (surfaceObject.getTag().equals("FinitePlanarWall2D")) ||
				 (surfaceObject.getTag().equals("FinitePlanarWall3D")) ) {
			    surf_type[ surfaceIndex ] = 1;
			    orientation[ surfaceIndex ] = mapAttributeToInteger(
										surfaceObject.getAttribute("surfaceNormal"), 
										new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("halfThickness1");
			    wallParam2[ surfaceIndex ] = surfaceObject.getDouble("halfThickness2");
			    wallParam3[ surfaceIndex ] = surfaceObject.getDouble("halfThickness3");
			    foundSurfaceType = true;
			}

                        /* omitted from latest version 
                        else if( surfaceObject.getTag().equals("BumpyWall") ){
                            surf_type[ surfaceIndex ] = 2;
                            orientation[ surfaceIndex ] = mapAttributeToInteger(
                                surfaceObject.getAttribute("surfaceNormal"), 
                                new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
                            wallParam[ surfaceIndex ] = surfaceObject.getDouble("bumpRadius");
                            wallParam2[ surfaceIndex ] = 0.0;
                            wallParam3[ surfaceIndex ] = 0.0;
                        }
			*/

			else if( (surfaceObject.getTag().equals("Colloids2D")) ||
				 (surfaceObject.getTag().equals("Colloids3D")) ) {
			    surf_type[ surfaceIndex ] = 3;
			    orientation[ surfaceIndex ] = 0;
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("radius");
			    wallParam2[ surfaceIndex ] = 0.0;
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			else if( (surfaceObject.getTag().equals("Pore2D")) ||
				 (surfaceObject.getTag().equals("Pore3D")) ) {
			    surf_type[ surfaceIndex ] = 5;
			    orientation[ surfaceIndex ] = 0;
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("radius");
			    wallParam2[ surfaceIndex ] = 0.0;
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			else if( (surfaceObject.getTag().equals("Atoms2D")) ||
				 (surfaceObject.getTag().equals("Atoms3D")) ) {
			    surf_type[ surfaceIndex ] = 6;
			    orientation[ surfaceIndex ] = 0;
			    wallParam[ surfaceIndex ] = 0.0;
			    wallParam2[ surfaceIndex ] = 0.0;
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			else if( (surfaceObject.getTag().equals("FinitePore2D")) ||
				 (surfaceObject.getTag().equals("FinitePore3D")) ) {
			    surf_type[ surfaceIndex ] = 7;
			    orientation[ surfaceIndex ] =  mapAttributeToInteger(
										 surfaceObject.getAttribute("longAxis"), 
										 new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("radius");
			    wallParam2[ surfaceIndex ] = surfaceObject.getDouble("length");
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			else if( (surfaceObject.getTag().equals("FinitePoreTaper2D")) ||
				 (surfaceObject.getTag().equals("FinitePoreTaper3D")) ){
			    surf_type[ surfaceIndex ] = 8;
			    orientation[ surfaceIndex ] =  mapAttributeToInteger(
										 surfaceObject.getAttribute("longAxis"), 
										 new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("radiusLBB");
			    wallParam2[ surfaceIndex ] = surfaceObject.getDouble("radiusRTF");
			    wallParam3[ surfaceIndex ] = surfaceObject.getDouble("length");
			    foundSurfaceType = true;
			}
			else if( surfaceObject.getTag().equals("FiniteCylinder3D") ){
			    surf_type[ surfaceIndex ] = 4;
			    orientation[ surfaceIndex ] =  mapAttributeToInteger(
										 surfaceObject.getAttribute("longAxis"), 
										 new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("radius");
			    wallParam2[ surfaceIndex ] = surfaceObject.getDouble("length");
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
                        else if( surfaceObject.getTag().equals("PeriodicCylinder3D") ){
                            surf_type[ surfaceIndex ] = 9;
			    orientation[ surfaceIndex ] =  mapAttributeToInteger(
										 surfaceObject.getAttribute("longAxis"), 
										 new String[] {"x", "y", "z"}, new int[] {0, 1, 2} );
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("meanRadius");
			    wallParam2[ surfaceIndex ] = surfaceObject.getDouble("amplitude");
			    wallParam3[ surfaceIndex ] = surfaceObject.getDouble("periodLength");
			    foundSurfaceType = true;
			}
                        else if( (surfaceObject.getTag().equals("RandomCylSphere2D")) ||
				 (surfaceObject.getTag().equals("RandomCylSphere3D")) ){
			    surf_type[ surfaceIndex ] = 10;
			    orientation[ surfaceIndex ] = 0;
			    wallParam[ surfaceIndex ] = surfaceObject.getDouble("radius");
			    wallParam2[ surfaceIndex ] = 0.0;
			    wallParam3[ surfaceIndex ] = 0.0;
			    foundSurfaceType = true;
			}
			if( (foundSurfaceType == false) && 
			    ( !(surfaceObject.getTag().equals("Reference"))) ) {
			    throw new RuntimeException( "Unknown surface type \"" +
							surfaceObject.getTag() +
							"\" found in Surface Types section");
			}
		    }
		}
		for( int i=0; i<nWall_type; i++ ) {
		    if( surfFound[ i ] == false ) {
			 // need to check surfFound for any surfaces we didn't find
			throw new RuntimeException( "Could not find a defination for the surface \"" +
						    surfaceName[i] +
						    "\" in the Surface Types seciton.  Every named surface" +
						    " must be defined." );
		    }
		}
	    }
	    else {
		throw new RuntimeException( "Could not find a set of surfaces in the Surface Types section.");
	    }
	}
	else {
	    throw new RuntimeException("Could not find the Surface Types section.");
	}

	// now try and read the treatment of reflected surfaces from the 
	// boundary conditions / reflected surfaces section
	XMLObject boundaryConditionsSection;
	if( ( boundaryConditionsSection = problemTypeSection.getChild("BoundaryConditions" + 
								      problemDimensionSuffix + 
								      polymericSuffix ) ) != null ) {
	    XMLObject surfaceReflectionsSection;
	    if( (surfaceReflectionsSection = boundaryConditionsSection.getChild( "SurfaceReflections" +
										 problemDimensionSuffix )) != null ) {
		XMLObject linkTable;
		if( (linkTable = surfaceReflectionsSection.getChild("Table")) != null ) {
		    nLink = linkTable.numChildren();
		    if( (nLink == 0) && (nWall_type > 0) ) {
			throw new RuntimeException("If any surfaces are defined, then at lest one Reflected Suface Treatment " +
						   "must be defined under Boundary Conditions / Treatment of reflected surfaces");
		    }
		    xTest_reflect_TF = new int[nLink][problemDimension];
		    boolean xTest_reflect_set[] = new boolean[nLink];
		    for( int i=0; i<nLink; i++ ) {
			XMLObject aLink;
			aLink = linkTable.getChild(i);
			if( aLink != null ) {
			    // find out the link number
			    int linkNumber = aLink.getInt("linkNumber");
			    // check if it's in range
			    if( (linkNumber < 0) || (linkNumber >=nLink) ) {
				throw new RuntimeException( "Link numbers in Boundary Conditions," +
							    "Treatment of reflected surfaces must between " +
							    "0 and total number of links minus 1");
			    }
			    // check that the user hasn't already defined this link
			    if( xTest_reflect_set[linkNumber] == true ) {
				throw new RuntimeException("Link number " + linkNumber + " is specified more than once " +
							   "in the Treatment of reflected surfaces section.  Each link number " +
							   "between 0 and the total number of links minus 1 should only be " +
							   "specified once.");
			    }
			    // remember that the user has defined this link number.
			    xTest_reflect_set[linkNumber] = true;

			    boolean xReflection = aLink.getBoolean("xReflections");
			    if( xReflection == true ) {
				xTest_reflect_TF[linkNumber][0] = 1;
			    }
			    else {
				xTest_reflect_TF[linkNumber][0] = 0;
			    }
			    if( problemDimension > 1 ) {
				boolean yReflection = aLink.getBoolean("yReflections");
				if( yReflection == true ) {
				    xTest_reflect_TF[linkNumber][1] = 1;
				}
				else {
				    xTest_reflect_TF[linkNumber][1] = 0;
				}
			    }
			    if( problemDimension > 2 ) {
				boolean zReflection = aLink.getBoolean("zReflections");
				if( zReflection == true ) {
				    xTest_reflect_TF[linkNumber][2] = 1;
				}
				else {
				    xTest_reflect_TF[linkNumber][2] = 0;
				}
			    }
			}
			else {
			    throw new RuntimeException( "One of the links in Treatment of reflected surfaces was null.");
			}
		    }
		    // check that we set all of the links
		    for( int i=0; i < nLink; i++ ) {
			if( xTest_reflect_set[i] == false ) {
			    throw new RuntimeException("Link number " + i + " was not specified in the Treatment of reflections " +
						       "section.  Each link number between 0 and the total number of links - 1 must " +
						       "be described.");
			}
		    }
		}
		else {
		    throw new RuntimeException( "Could not find the surface link table in the Surface Reflections section.");
		}

	    }
	    else {
		throw new RuntimeException( "Cound not find Surface reflections section.");

	    }
	    

	}
	else {
	    throw new RuntimeException( "Could not find the Boundary conditions section.");
	}
	    
    }




    public void validateSurfaceInteractionParameters() throws RuntimeException {
	// allocate space for the surface-surface and surface-fluid arrays
	sigma_ww = new double[nWall_type][nWall_type];
	eps_ww = new double[nWall_type][nWall_type];
	rho_s = new double[nWall_type];
	cut_ww = new double[nWall_type][nWall_type];
	sigma_wf = new double[nComp][nWall_type];
	eps_wf = new double[nComp][nWall_type];
	cut_wf = new double[nComp][nWall_type];

	boolean[][] foundWWData = new boolean[nWall_type][nWall_type];
	boolean[][] foundWFData = new boolean[nComp][nWall_type];
	
	// from the problem type xml object look at its children until we find
	// the surface iteraction section
	XMLObject theSurfacesSection;
	for(int child=0; child<problemTypeSection.numChildren(); child++) {
	    theSurfacesSection = problemTypeSection.getChild(child);
	    if( theSurfacesSection.getTag().matches("NoWallSurfaceInteractions.*")) {
		// mark that we found data for WW and WF since no data is needed
		// and we need to make it back the checks later on.
		for(int i=0; i<nWall_type; i++) {
		    for(int j=0; j<nWall_type; j++ ) {
			foundWWData[i][j] = true;
		    }
		    for(int k=0; k<nComp; k++) {
			foundWFData[k][i] = true;
		    }
		}
	    }
	    else if( theSurfacesSection.getTag().matches(".*SurfaceInteractions.*") ) {
		// now look at this sections children to see if we can find
		// the type of information stored
		for( int surfChild=0; surfChild < theSurfacesSection.numChildren(); surfChild++) {
		    XMLObject surfChildSection = theSurfacesSection.getChild(surfChild);
		    if( surfChildSection.getTag().matches("MixingType.*") ) {
			// check that mixing type here is the same as it was earlier
			if( ((surfChildSection.getTag().matches("MixingTypeLB.*")) && (lMix_rule == 1 )) ||
			    ((surfChildSection.getTag().matches("MixingTypeUD.*")) && (lMix_rule == 0 )) ) {
			    throw new RuntimeException("Mixing type specified under Functionals and Interactions " +
						       "(Lorentz-Berthelot or User Defined) is different than type " +
						       "specified under Surface Interactions.");
			}
			
			// now look at the name of this child to see what type of data we'll need to read
			boolean readSigma = false;
			boolean readEps = false;
			boolean readCut = false;
			if( surfChildSection.getTag().matches(".*Sigma.*") ) {
			    readSigma = true;
			}
			if( surfChildSection.getTag().matches(".*Eps.*") ) {
			    readEps = true;
			}
			if( surfChildSection.getTag().matches(".*Cut.*") ) {
			    readCut = true;
			}
			
			// the MixingType* section will have one child
			// that child will contian one or two tables for
			// LB or UD mixing respectively
			
			// get next child
			XMLObject tableWrapperClass = surfChildSection.getChild(0);
			if( tableWrapperClass != null ) {
			    for( int tableWrapperChild = 0; tableWrapperChild < tableWrapperClass.numChildren(); tableWrapperChild++) {
				XMLObject aTable = tableWrapperClass.getChild( tableWrapperChild );

				// look for the *WallTable which only as wall-wall data
				if( aTable.getAttribute("name").matches(".*WallTable") ) {
				    for( int entry=0; entry < aTable.numChildren(); entry++) {
					XMLObject anEntry = aTable.getChild( entry );
					String surf1Name, surf2Name;
					int surf1Index, surf2Index;
					// look at the children of this entity to find out which surfaces 
					// are beeing described
					XMLObject entryChild = anEntry.getChild(0);
					if( entryChild != null ) {
					    XMLObject entryChildPanel = entryChild.getChild(0);
					    surf1Name = entryChildPanel.getAttribute("surfaceName");
					    try {
						surf1Index = findIndexInStringArray( surf1Name, surfaceName );
					    }
					    catch( RuntimeException e ) {
						throw new RuntimeException( "Couldn't find surface named \"" +
									    surf1Name +
									    "\" in the set of surfaces defined in the Surface Names section." );
					    }
					}
					else {
					    throw new RuntimeException("Couldn't find child panel in Surface Interactions section");
					}
					if( lMix_rule == 1 ) {
					    // look for other surface name
					    entryChild = anEntry.getChild(1);
					    if( entryChild != null ) {
						XMLObject entryChildPanel = entryChild.getChild(0);
						surf2Name = entryChildPanel.getAttribute("surfaceName");
						try {
						    surf2Index = findIndexInStringArray( surf2Name, surfaceName );
						}
						catch( RuntimeException e ) {
						    throw new RuntimeException( "Couldn't find surface named \"" +
										surf2Name +
										"\" in the set of surfaces defined in the Surface Names section." );
						}
					    }
					    else {
						throw new RuntimeException("Couldn't find child panel in Surface Interactions section");
					    }
					}
					else {
					    // just reading symmetric values
					    surf2Name = surf1Name;
					    surf2Index = surf1Index;
					}

					// we have the surface names, now note that we're reading these
					if( foundWWData[surf1Index][surf2Index] == true ) {
					    if( lMix_rule == 0 ) {
						throw new RuntimeException("Surface \"" + surf1Name + 
									   "\" appears more than once in the " +
									   "Surface Interactions Section");
					    }
					    else {
						throw new RuntimeException("Surface combination \"" + surf1Name + "\" and \"" + 
									   surf2Name + "\" appears more than once in the " +
									   "Surface Interactions Section");
					    }
					}
					else {
					    foundWWData[surf1Index][surf2Index] = true;
					}

					// ok we can read some data now
					if( readSigma == true ) {
					    sigma_ww[surf1Index][surf2Index] = anEntry.getDouble("sigmaWW");
					}
					if( readEps == true ) {
					    eps_ww[surf1Index][surf2Index] = anEntry.getDouble("epsWW");
					    rho_s[surf1Index] = anEntry.getDouble("rhoW");

					}
					if( readCut == true ) {
					    cut_ww[surf1Index][surf2Index] = anEntry.getDouble("cutWW");
					}
				    }
				}

				// look for the *WallFluidTable which has wall-fluid data
				if( aTable.getAttribute("name").matches(".*WallFluidTable") ) {
				    for( int entry=0; entry < aTable.numChildren(); entry++) {
					XMLObject anEntry = aTable.getChild( entry );
					String comp1Name, surf1Name;
					int comp1Index, surf1Index;
					// look at the children of this entity to find out which surface and fluid
					// are beeing described

					XMLObject entryChild1 = anEntry.getChild(0);
					XMLObject entryChild2 = anEntry.getChild(1);

					if( (entryChild1 != null) && (entryChild2 != null) ) {
					    XMLObject entryChild1Panel = entryChild1.getChild(0);
					    XMLObject entryChild2Panel = entryChild2.getChild(0);

					    if( (entryChild1Panel.hasAttribute("surfaceName")) &&
						(entryChild2Panel.hasAttribute("name")) ) {
						surf1Name = entryChild1Panel.getAttribute("surfaceName");
						comp1Name = entryChild2Panel.getAttribute("name");

					    } 
					    else if( (entryChild1Panel.hasAttribute("name")) &&
						(entryChild2Panel.hasAttribute("surfaceName")) ) {
						surf1Name = entryChild2Panel.getAttribute("surfaceName");
						comp1Name = entryChild1Panel.getAttribute("name");

					    }
					    else {
						throw new RuntimeException("Couldn't find surface and componet names in " +
									   "Surface Interactions section");
					    }

					    try {
						surf1Index = findIndexInStringArray( surf1Name, surfaceName );
						comp1Index = findIndexInStringArray( comp1Name, compName );
					    }
					    catch( RuntimeException e ) {
						throw new RuntimeException( "Couldn't find surface named \"" +
									    surf1Name +
									    "\" or component named \"" +
									    comp1Name +
									    "\" in the set of surfaces defined in the Surface Names section." );
					    }
					}
					else {
					    throw new RuntimeException("Couldn't find child panel in Surface Interactions section");
					}

					// we have the surfac/component names, now note that we're reading these
					if( foundWFData[comp1Index][surf1Index] == true ) {
					    throw new RuntimeException("Surface, component combination \"" + surf1Name + "\" and \"" + 
								       comp1Name + "\" appears more than once in the " +
								       "Surface Interactions Section");
					}
					else {
					    foundWFData[comp1Index][surf1Index] = true;
					}

					// ok we can read some data now
					if( readSigma == true ) {
					    sigma_wf[comp1Index][surf1Index] = anEntry.getDouble("sigmaWF");
					}
					if( readEps == true ) {
					    eps_wf[comp1Index][surf1Index] = anEntry.getDouble("epsWF");
					}
					if( readCut == true ) {
					    cut_wf[comp1Index][surf1Index] = anEntry.getDouble("cutWF");
					}
				    }
				}
			    }
			}
			else {
			    throw new RuntimeException("Couldn't find Surface Interactions data table.");
			}
		    }
		}
	    }
	}

	// now check that everything was set
	for( int i=0; i<nWall_type; i++ ) {
	    if( lMix_rule == 0 ) {
		if( foundWWData[i][i] != true ) {
		    throw new RuntimeException("No surface interaction data defined for surface \"" + surfaceName[i] +
					       "\" in Surface Interactions");
		}
	    }
	    else {
		for( int j=0; j<nWall_type; j++ ) {
		    if( foundWWData[i][j] != true ) {
			if( foundWWData[j][i] == true ) {
			    sigma_ww[i][j] = sigma_ww[j][i];
			    eps_ww[i][j] = eps_ww[j][i];
			    cut_ww[i][j] = cut_ww[j][i];
			    foundWWData[i][j] = true;
			}
			else {
			    throw new RuntimeException("No surface interaction data defined for the combination \"" + surfaceName[i] +
						       "\", \"" + surfaceName[j] + "\" in Surface Interactions.");
			}
		    }
		}
		for( int j=0; j< nComp; j++ ) {
		    if( foundWFData[i][j] != true ) {
			throw new RuntimeException("No surface/component interaction data defined for the combination \"" +
						   surfaceName[i] + "\", \"" + compName[j] + "\" in Surface Interactions.");
		    }
		}
	    }
	}
    }


    public void validateSurfacesFile() throws RuntimeException {

	String surfacesFileName = "dft_surfaces.dat";
	// use a JFileChooser to have the user help us
	JFileChooser chooser = new JFileChooser();
	if( (workingDirectory != null) && (workingDirectory != "") ) {
	    chooser.setCurrentDirectory( new File( workingDirectory + 
						   System.getProperty("file.separator") +
						   surfacesFileName ) );
	    chooser.setSelectedFile( new File( workingDirectory +
					     System.getProperty("file.separator") +
					     surfacesFileName ) );
	}
	else {
	    chooser.setCurrentDirectory( new File( System.getProperty("user.dir")  + 
						   System.getProperty("file.separator") +
						   surfacesFileName ) );
	    chooser.setSelectedFile( new File(System.getProperty("user.dir") +
					     System.getProperty("file.separator") +
					     surfacesFileName ) );
	}
	chooser.setDialogTitle("Select a dft_surfaces.dat file to check...");
	int returnVal = chooser.showDialog(null, "Check File");
	if(returnVal == JFileChooser.APPROVE_OPTION) {
	    try {
		surfacesFileName = chooser.getSelectedFile().getCanonicalPath();
	    } 
	    catch (IOException e) {
		throw new RuntimeException( e.toString() );
	    }
	}
	else if( returnVal == JFileChooser.CANCEL_OPTION ) {
	    return;
	}
	else {
	    throw new RuntimeException("Couldn't find \"dft_surfaces.dat\" file");
	}
	File surfacesFile = new File( surfacesFileName );
	if( (surfacesFile.exists() == false) || (surfacesFile.canRead() == false) ) {
	    throw new RuntimeException("Couldn't read the surfaces file \"" + surfacesFileName + "\"");
	}
	
	// now we have a file that we can find and read so analyzie it
	// we'll check for the following things.
	// 1) Each line of the file contains the right number of values
	// 2) That the total number of walls equals the value set as nWall 
	// 3) That the total number of links equals nLink and that each 
	// link is used
	
	
	BufferedReader surfaceData = null;
	try {
	    surfaceData = new BufferedReader( new FileReader( surfacesFile ) );
	}
	catch (FileNotFoundException e) {
	    throw new RuntimeException( e.toString() );
	}
	
	if( surfaceData != null ) {
	    // count the number of walls found and if all the wall
	    // typs and link types are used
	    int numWallsInFile = 0;
	    boolean[] foundSurfaceType = new boolean[ nWall_type ];
	    boolean[] foundLinkType = new boolean[ nLink ];
	    
	    int lineNumber = 0;
	    
	    boolean done = false;
	    while( !done ) {
		String aLine;
		boolean status;
		try {
		    aLine = surfaceData.readLine();
		    status = surfaceData.ready();
		    if( status == false ) {
			done = true;
		    }
		    if( numWallsInFile > 10 ) {
			done = true;
		    }
		} 
		catch (IOException e) {
		    throw new RuntimeException("Error reading \"" + surfacesFileName + "\" " + e.toString() );
		}
		
		// we have a string now, read it's tokens
		lineNumber++;
		StringTokenizer stk = new StringTokenizer( aLine );
		
		// check the number of things on each line
		if( stk.countTokens() == 0 ) {
		    throw new RuntimeException( "Line " + lineNumber +
						" of the surfaces file is blank, which is not allowed.");
		}

		if( stk.countTokens() < 3 ) {
		    throw new RuntimeException( "Line " + lineNumber +
						" of the surfaces file has too few values on it. " +
						"Should be 3, 4 or 5 for 1D, 2D or 3D problem respectively.");
		}
		if( stk.countTokens() > 5 ) {
		    throw new RuntimeException( "Line " + lineNumber +
						" of the surfaces file has too many values on it. " +
						"Should be 3, 4 or 5 for 1D, 2D or 3D problem respectively.");
		}
		switch( nDim ) {
		case 1:
		    if( stk.countTokens() != 3 ) {
			throw new RuntimeException("Line " + lineNumber + 
						   " of the surfaces file has the wrong number of values for a 1D problem.");
		    }
		    break;
		case 2:
		    if( stk.countTokens() != 4 ) {
			throw new RuntimeException("Line " + lineNumber + 
						   " of the surfaces file has the wrong number of values for a 2D problem.");
		    }
		    break;
		case 3:
		    if( stk.countTokens() != 5 ) {
			throw new RuntimeException("Line " + lineNumber + 
						   " of the surfaces file has the wrong number of values for a 3D problem.");
		    }
		    break;
		}
		
		// now try and read various parts of the line
		int aSurfaceType;
		try {
		    aSurfaceType = Integer.parseInt( stk.nextToken() );
		}
		catch( NumberFormatException e ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file has a surface type number that cannot be read.");
		}
		catch( NoSuchElementException e ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file is missing a surface type number.");
		}
		if( (aSurfaceType < 0 ) || (aSurfaceType >= nWall_type) ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file has a surface type of " + aSurfaceType +
					       " which is out of the range of surface types defined in the input file: 0, " +
					       (nWall_type  - 1) );
		}
		foundSurfaceType[ aSurfaceType ] = true;
		numWallsInFile++;
		
		
		int aLinkNumber;
		try {
		    aLinkNumber = Integer.parseInt( stk.nextToken() );
		}
		catch( NumberFormatException e ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file has a link number that cannot be read.");
		}
		catch( NoSuchElementException e ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file is missing a link number.");
		}
		if( (aLinkNumber < 0 ) || (aLinkNumber >= nLink) ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file has a link number of " + aLinkNumber +
					       " which is out of the range of link numbers defined in the input file: 0, " +
					       (nLink  - 1) );
		}
		foundLinkType[ aLinkNumber ] = true;
		
		
		double anX;
		try {
		    anX = Double.parseDouble( stk.nextToken() );
		}
		catch( NumberFormatException e ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file has an X position that cannot be read.");
		}
		catch( NoSuchElementException e ) {
		    throw new RuntimeException("Line " + lineNumber +
					       " of the surfaces file is missing an X position.");
		}
		
		if( nDim >= 2 ) {
		    double anY;
		    try {
			anY = Double.parseDouble( stk.nextToken() );
		    }
		    catch( NumberFormatException e ) {
			throw new RuntimeException("Line " + lineNumber +
						   " of the surfaces file has a Y position that cannot be read.");
		    }
		    catch( NoSuchElementException e ) {
			throw new RuntimeException("Line " + lineNumber +
						   " of the surfaces file is missing a Y position.");
		    }
		}
		
		if( nDim == 3 ) {
		    double anZ;
		    try {
			anZ = Double.parseDouble( stk.nextToken() );
		    }
		    catch( NumberFormatException e ) {
			throw new RuntimeException("Line " + lineNumber +
						   " of the surfaces file has a Z position that cannot be read.");
		    }
		    catch( NoSuchElementException e ) {
			throw new RuntimeException("Line " + lineNumber +
						   " of the surfaces file is missing a Z position.");
		    }
		}
		
		if( numWallsInFile > nWall) {
		    throw new RuntimeException("Total number of walls listed in surfaces file is " + numWallsInFile +
					       ", which is greater than the number listed in the input data, " + nWall );
		}

	    }

	    // done reading so check that all surface types defined were used 
	    // and that all links defined were used as well.
	    for( int i=0; i<nWall_type; i++) {
		if( foundSurfaceType[i] == false ) {
		    throw new RuntimeException("Surface type, " + i + ", or \"" + surfaceName[i] +
					       "\" was defined in the input data , but no such surfaces were " +
					       "found in the sufaces file.");
		}
	    }
	    
	    for( int i=0; i<nLink; i++ ) {
		if( foundLinkType[i] == false ) {
		    throw new RuntimeException("Link number, " + i + ", was defined in the input data but no use of that link " +
					       " was found in the surfaces file.");
		}
	    }
	    
	}
	else {
	    throw new RuntimeException("Couldn't open \"" + surfacesFileName + "\"" );
	}
	
    }
    
    


    /*
      Following are the write...() methods inalphabetical order 
    */


    private void writeChargedSurfaceBoundaryConditionParameters(PrintWriter theOutput) {
	theOutput.println("********* CHARGED SURFACE BOUNDARY CONDITIONS *********************************");
	theOutput.print("@ ");
	for(int wall=0; wall<nWall_type; wall++) {
	    theOutput.print(typeBCElectric[wall] + " ");
	}
	theOutput.println("\tType_bc_elec[iwall_type]:type of electrical boundary condition\n" + 
                          "\t             0=neutral surface, \n" +
                          "\t             1=const potential,\n" +
                          "\t             2=const surface charge, \n" +
                          "\t             3=atomic charges\n" +
                          "\t *Note:  The read statement for this parameter has been removed\n" +
                          "\t         for the special case of hard spheres with charge");
	theOutput.print("@ ");
	for(int wall=0; wall<nWall_type; wall++) {
	    theOutput.print(elecParamW[wall] + " ");
	}
	theOutput.println("\tElec_param_w[nWall_type] valaue of electricostatic boundary condition");

	theOutput.print("@ ");
	if(chargeProfile == true) {
	    theOutput.print(-1);
	}
	else {
	    theOutput.print( nLocalCharge );
	}
	theOutput.println("\tNlocal_charge, the number of local charges");
	theOutput.println("\t\t(-1 for linear profile of charge between");
	theOutput.println("\t\ttwo points aligned with the princible axes. !!!)");

	theOutput.print("@ ");
	for(int charge=0; charge<nLocalCharge; charge++) {
	    theOutput.print(chargeLoc[charge] + " ");
	}
	theOutput.println("\tCharge_loc[Nlocal_charge] total charge at each site");

	theOutput.print("@ ");
	for(int charge=0; charge < nLocalCharge; charge++ ) {
	    theOutput.print(chargeDiam[charge] + " ");
	}
	theOutput.println("\tCharge_Diam[Nlocal_charge] charge diameter");

	theOutput.print("@ ");
	for(int charge=0; charge < nLocalCharge; charge++) {
	    for(int dim=0; dim < nDim; dim++) {
		theOutput.print(chargeX[charge][dim] + " ");
	    }
	}
	theOutput.println("\tCharge_x[0][0], [0][Ndim] ... [NLocal_charge][Ndim]");

	theOutput.println("@ " + lPointChargeAtoms + "  " + lPointChargeLocal + "\tCharge_type_atoms Charge_type_local");
    }


    private void writeCoarseningSwitches(PrintWriter theOutput) {
	theOutput.println("********* COARSENING SWITCHES *************************************************");
	
	theOutput.println("@ " + nZone +
			  "\tNzone (Coarsens Mesh/Jacobian by a factor of 2)");
	theOutput.print("@ ");
	for(int zn=0; zn < rMaxZone.length; zn++ ) {
	    theOutput.print(rMaxZone[zn] + " ");
	}
	theOutput.println("\tRmax_zone[Nzone-1] [0.0 for complete coarsening]");

	theOutput.println("@ " + coarsenResidual +
			  "\tCoarsen_Residual ? (0=No, 1=Yes)");
	theOutput.println("@ " + coarserJac +
			  " " + eSizeJacobian +
			  "\tCoarser_jac; Esize_jacobian" );
        theOutput.println("\t\t0 =Jac. zones are the same as resid zones.\n" +
                          "\t\t1 =coarsen finest Jacobian zone by fac of 2\n" +
                          "\t\t2 =coarsen all but coarsest zone by fac of 2\n" +
                          "\t\t3 =use coarsest zone everywhere\n" +
                          "\t\t4 =use 2nd coarsest zone in all but coarsest zone\n" +
                          "\t\t5 =use Esize_jacobian for all Jacobian integrals.");
	theOutput.println("@ " + lJacCut +
			  " " + jacThreshold +
			  "\tLjac_cut Jac_threshold");
	theOutput.println("@ " + matrixFillFlag +
			  "\tMatrix_fill_flag");
        theOutput.println("\t\t0= Exact Jacobian,\n" +
                          "\t\t1= Jac_Save1 (only rhobar2 and rhobar 3)\n" +
                          "\t\t2= Jac_Save2 (rhobar 2,3 exact: rhobar 0,1 estimated)\n" +
                          "\t\t3= Enumerated rhobar equations (all of them)\n" + 
                          "\t\t4= Enumerated rhobar equations (no vector terms)\n");
    }


    private void writeDielectricConstantParameters(PrintWriter theOutput) {
	theOutput.println("********* DIELECTRIC CONSTANT PARAMETERS **************************************");

	theOutput.println("@ " + typeDielec + "\tDielectric is 0=constant everywhere");
	theOutput.println("\t1 = different between surfaces and fluids");
	theOutput.println("\t2 = different between pore fluid and bulk fluid");
	theOutput.println("\t3 = constant in walls, density dependent in fluid.");
	
	theOutput.println("@ " + dielecBulk + 
			  " " + dielecPore +
			  " " + dielecX + "\tDielec_bulk Dielec_pore Dielec_X" );

	theOutput.print("@ ");
	for(int wall=0; wall < nWall_type; wall++) {
	    theOutput.print( dielecWall[wall] + " " );
	}
	theOutput.println("\tDielec_wall[Nwall_type]");
    }


    private void writeDimensionParameters(PrintWriter theOutput) {
	theOutput.println("********* DIMENSION PARAMETERS ************************************************");
	theOutput.println("@ " + referenceLength + "  " +
			  referenceDensity + "  " +
			  referenceTemp + "  " +
			  referenceDielectric + "  " +
			  maximumPotential + 
			  "  Length_ref  Density_ref  Temp  Dielec_ref  VEXT_MAX");
	theOutput.println("");
    }

    private void writeFluidNames(PrintWriter theOutput) {
	theOutput.println("********* FLUID/COMPONENT NAMES ***********************************************");
	theOutput.print("  Fluid/Component names: " );
	if( compName.length >= 1 ) {
	    theOutput.print( compName[0] );
	}
	for(int i=1; i < compName.length; i++ ) {
	    theOutput.print( ", " + compName[i] );
	}
	theOutput.println("");
	theOutput.println("");
    }

    private void writeFluidSurfaceInteractionParameters(PrintWriter theOutput) {
	theOutput.println("************** FLUID INTERACTION PARAMETERS ***********************************");
	theOutput.println("@ " + nComp + "  " + lMix_rule +
			  "\tNcomp (no. of components) OR for polymers Nblock_tot\n" +
			  "\t\tMixing Rule Type (0=L-B, 1=manual input)");

	for(int i=0; i<nComp; i++) {
	    theOutput.print(compName[i] + " ");
	}
	theOutput.println("\tName of each comonent");

	theOutput.print("@ ");
	for(int i=0; i<nComp; i++) {
	    theOutput.print(compMass[i] + " ");
	}
	theOutput.println("\tMass of each component to include deBroglie effects");

	theOutput.print("@ ");
	for(int i=0; i<nComp; i++) {
	    theOutput.print(charge_f[i] + " ");
	}
	theOutput.println("\tCharge_f(icomp): icomp=1,Ncomp; Valence ");

	
	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	    for(int i=0; i<nComp; i++) {
		theOutput.print(sigma_ff[i][i] + " ");
	    }
	}
	else {
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nComp; j++) {
		    theOutput.print(sigma_ff[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tSigma_ff(icomp, icomp): icomp=1,Ncomp (Nblock_tot)");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	    for(int i=0; i<nComp; i++) {
		theOutput.print(eps_ff[i][i] + " ");
	    }
	}
	else {
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nComp; j++) {
		    theOutput.print(eps_ff[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tEps_ff(icomp, icomp): icomp=1,Ncomp (Nblock_tot)");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	    for(int i=0; i<nComp; i++) {
		theOutput.print(cut_ff[i][i] + " ");
	    }
	}
	else {
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nComp; j++) {
		    theOutput.print(cut_ff[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tCut_ff(icomp, icomp): icomp=1,Ncomp (Nblock_tot)");

	theOutput.println("");
	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	    theOutput.print(rho_s[i] + " ");
	}
	theOutput.println("\tRho_w(iwall): iwall=1,Nwall_type");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	    for(int i=0; i<nWall_type; i++) {
		theOutput.print(sigma_ww[i][i] + " ");
	    }
	}
	else {
	    for(int i=0; i<nWall_type; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(sigma_ww[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tSigma_ww(iwall, iwall): iwall=1,Nwall_type");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	    for(int i=0; i<nWall_type; i++) {
		theOutput.print(eps_ww[i][i] + " ");
	    }
	}
	else {
	    for(int i=0; i<nWall_type; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(eps_ww[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tEps_ww(iwall, iwall): iwall=1,Nwall_type");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	    for(int i=0; i<nWall_type; i++) {
		theOutput.print(cut_ww[i][i] + " ");
	    }
	}
	else {
	    for(int i=0; i<nWall_type; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(cut_ww[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tCut_ww(iwall, iwall): iwall=1,Nwall_type");
	  
	theOutput.println("");
	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	}
	else {
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(sigma_wf[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tSigma_wf(icomp, iwall)");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	}
	else {
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(eps_wf[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tEps_wf(icomp, iwall)");

	theOutput.print("@ ");
	if( lMix_rule == 0 ) {
	}
	else {
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(cut_wf[i][j] + " ");
		}
	    }
	}
	theOutput.println("\tCut_wf(icomp, iwall)");


	theOutput.println("    *Note for polymers: --- treat each segment (or block) TYPE as a distinct");
        theOutput.println("     component in this section => 1 for homopolymer, 2 for diblock or");
        theOutput.println("      ABA triblock, 3 for ABC triblock 3 for diblock with solvent, etc.");

    }


    private void writeFunctionalSwitches(PrintWriter theOutput) {
	theOutput.println("********* FUNCTIONAL SWITCHES *************************************************");
	theOutput.println("@ " + type_func + 
			  "\tType_func (-1=No HS, 0=Rosen1, 1=Rosen2, 2=LDA, 3=GHRM, 4=GVDWM)");
	theOutput.println("@ " + type_attr +
			  "\tType_attr (-1=No attr, 0=strict MF, 1=B2, 2=B2 ions & MF solvent)");
	theOutput.println("@ " + type_coul +
			  "\tType_coul (-1=No coul, 0=strict MF, 1=include 2nd order corr)");
	theOutput.println("@ " + type_poly +
			  "\tType_poly (-1=No polymer, 0-2 different polymer formulations)");
	theOutput.println("@ " + compareToFastram + "\tLcompare_fastram");
    }


    private void writeLinearSolverParameters(PrintWriter theOutput) {
	theOutput.println("********* LINEAR SOLVER PARAMETERS ********************************************");
	
	theOutput.println("@ " + solverFlag + "  " + kspace +
			  "\tSolver (0=gmres, 1=cg, 2=tfqmr, 3=cg2, 4=bicgstab) kspace\n" +
                          "\t\t(for Solver=gmres second parameter is kspace which must be > 0)");
	theOutput.println("@ " + scalingFlag +
			  "\tScaling (0=row sum, 1=Jacobi, 2=symrow sum, -1=none)");
	theOutput.println("@ " + preconditionerFlag + "  " + iLutFillParam + 
			  "\tPreconditioner (0=ilu, 1=Jacobi, 2=symGC, 3=LSPoly3,");
	theOutput.println("\t\t4=ilut (3 levels), 5=ilut (7 levels), -1=none)");
	theOutput.println("@ " + maxLinearSolverIter + "  " + convergenceTol +
			  "\t\tMax iterations and Convergence Tolerance for Linear Solver");
    }


    private void writeLocaContinuationParameters(PrintWriter theOutput) {
        theOutput.println("********* LOCA CONTINUATION PARAMETERS ****************************************");
	theOutput.println("@ " + locaContinuationMethod +
			  "\tContinuation Method (-1=None, 0,1,2=0th, 1st, arc-length");
	theOutput.println("\t3=Spinodal (Turning Point), 4=Binodal (Phase Eq))");

	theOutput.println("@ " + locaContinuationParameter +
			  "\tContinuation Parameter :  Scale_fac (for CONT_SCALE cases only)");

        theOutput.println("\t\tCONST_TEMP           1");
        theOutput.println("\t\tCONST_RHO_0          2");
        theOutput.println("\t\tCONST_RHO_ALL        3");
        theOutput.println("\t\tCONST_LOG_RHO_0      4");
        theOutput.println("\t\tCONST_LOG_RHO_ALL    5");
        theOutput.println("\t\tCONST_SCALE_RHO      6");
        theOutput.println("\t\tCONST_EPSW_0         7");
        theOutput.println("\t\tCONST_EPSW_ALL       8");
        theOutput.println("\t\tCONST_SCALE_EPSW     9");
        theOutput.println("\t\tCONST_EPSWF00        10");
        theOutput.println("\t\tCONST_EPSWF_ALL_0    11");
        theOutput.println("\t\tCONST_SCALE_EPSWF    12");
        theOutput.println("\t\tCONST_EPSFF_00       13");
        theOutput.println("\t\tCONST_EPSFF_ALL      14");
        theOutput.println("\t\tCONST_SCALE_EPSFF    15");
        theOutput.println("\t\tCONST_SCALE_CHG      16");


	theOutput.println("@ " + locaParameterStep +
			  "\tParameter initial step size");

	theOutput.println("@ " + locaNumSteps +
			  " " + locaStepControl +
			  "\tN Steps, Step Control Aggressiveness (0.0 = constant step)");

        theOutput.println("********* END OF INPUT FILE ***************************************************");
    }


    private void writeMeshContinuationParameters(PrintWriter theOutput) {
        theOutput.println("********* MESH CONTINUATION PARAMETERS ****************************************");

        theOutput.println("\t\tHere you enter information for mesh continuation.");
        theOutput.println("\t\tLoca does not do mesh continuation because it is not.");
        theOutput.println("\t\ta continuous variable.");
	
	theOutput.println("@ " + numRuns +  "\tN_runs");
	theOutput.print("@ ");
	for(int i=0; i<nDim; i++ ) {
	    theOutput.print(parameterChange[i] + "  ");
	}
	theOutput.println("\tDel_1[] How much to change parameter");

	theOutput.println("@ " + planeNewNodes +
			  " " + positionNewNodes +
			  "\tPlane_new_nodes Pos_new_nodes");
	theOutput.println("\t(0=yz, 1=xz, 2=xy)  (-1=lbb, 0=center, 1=rtf)");

	theOutput.println("@ " + guessRange[0] +
			  " " + guessRange[1] +
			  "\tXGuess_range[0,1] parameters for mesh continuation\n" +
			  "\t\tGuess_range[0] is the surf separation to stop using 100% Rho_b\n" +
			  "\t\tGuess_range[1] is the surf separation to start using 100% X_old\n");
    }

 
    private void writeMeshParameters(PrintWriter theOutput) {

	if( realSpaceComputation == true ) {
	    theOutput.println("********* MESH PARAMETERS *****************************************************");
	    theOutput.println("@ " + nDim + "\tNdim");
	    theOutput.print("@ ");
	    for(int i=0; i<nDim; i++)
		theOutput.print(size_x[i] + " ");
	    theOutput.println("\tSize_x(idim); idim=1,Ndim");
	    theOutput.print("@ ");
	    for(int i=0; i<nDim; i++)
		theOutput.print(esize_x[i] + " ");
	    theOutput.println("\tEsize_x(idim); idim=1,Ndim");
	    for(int i=0; i<3; i++) {
		theOutput.print("@ ");
		if(i < nDim) {
		    theOutput.print(type_bc_lbb[i] + " " + type_bc_rtf[i]);
		}
		else {
		    theOutput.print("-1 -1");
		}
		theOutput.println("\tType_bc -1=wall, 0=bulk, 1=periodic, 2=ref");
	    }
	} 
    }

    private void writeMixingParameters(PrintWriter theOutput) {

	theOutput.println("...Mixing Parameters...");

	// write fluid or fluid-fluid interactions
	if( lMix_rule == 0 ) {
	    theOutput.print("@ ");
	    for(int i=0; i<nComp; i++) {
		theOutput.print(sigma_ff[i][i] + "  ");
	    }
	    theOutput.println("\tSigma_ff(icomp,icomp):icomp=1,Ncomp (Nblock_tot) only diagonal elements");

	    theOutput.print("@ ");
	    for(int i=0; i<nComp; i++) {
		theOutput.print(eps_ff[i][i] + "  ");
	    }
	    theOutput.println("\tEps_ff(icomp,icomp):icomp=1,Ncomp (Nblock_tot) only diagonal elements");

	    theOutput.print("@ ");
	    for(int i=0; i<nComp; i++) {
		theOutput.print(cut_ff[i][i] + "  ");
	    }
	    theOutput.println("\tcut_ff(icomp,icomp):icomp=1,Ncomp (Nblock_tot) only diagonal elements");
	}
	else {
	    theOutput.print("@ ");
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nComp; j++) {
		    theOutput.print(sigma_ff[i][j] + "  ");
		}
	    }
	    theOutput.println("\tSigma_ff(icomp,icomp):icomp=1,Ncomp (Nblock_tot)");

	    theOutput.print("@ ");
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nComp; j++ ) {
		    theOutput.print(eps_ff[i][j] + "  ");
		}
	    }
	    theOutput.println("\tEps_ff(icomp,icomp):icomp=1,Ncomp (Nblock_tot)");

	    theOutput.print("@ ");
	    for(int i=0; i<nComp; i++) {
		for( int j=0; j<nComp; j++ ) {
		    theOutput.print(cut_ff[i][j] + "  ");
		}
	    }
	    theOutput.println("\tcut_ff(icomp,icomp):icomp=1,Ncomp (Nblock_tot)");
	}

	// write surface or surface-surface interactions
	if( lMix_rule == 0 ) {
	    theOutput.print("@ ");
	    for(int i=0; i<nWall_type; i++) {
		theOutput.print(sigma_ww[i][i] + "  ");
	    }
	    theOutput.println("\tSigma_ww(icomp,icomp):iwall_type=1,Nwall_type only diagonal elements");

	    theOutput.print("@ ");
	    for(int i=0; i<nWall_type; i++) {
		theOutput.print(eps_ww[i][i] + "  ");
	    }
	    theOutput.println("\tEps_ww(icomp,icomp):iwall_type=1,Nwall_type only diagonal elements");

	    theOutput.print("@ ");
	    for(int i=0; i<nWall_type; i++) {
		theOutput.print(cut_ww[i][i] + "  ");
	    }
	    theOutput.println("\tcut_ww(icomp,icomp):iwall_type=1,Nwall_type only diagonal elements");
	}
	else {
	    theOutput.print("@ ");
	    for(int i=0; i<nWall_type; i++) {
		for( int j=0; j<nWall_type; j++) {
		    theOutput.print(sigma_ww[i][j] + "  ");
		}
	    }
	    theOutput.println("\tSigma_ww(nwall_type,nwall_type):nwall_type=1,nWall_type");

	    theOutput.print("@ ");
	    for(int i=0; i<nWall_type; i++) {
		for( int j=0; j<nWall_type; j++ ) {
		    theOutput.print(eps_ww[i][j] + "  ");
		}
	    }
	    theOutput.println("\tEps_ww(nwall_type,nwall_type):nwall_type=1,nWall_type");

	    theOutput.print("@ ");
	    for(int i=0; i<nWall_type; i++) {
		for( int j=0; j<nWall_type; j++ ) {
		    theOutput.print(cut_ww[i][j] + "  ");
		}
	    }
	    theOutput.println("\tcut_ww(nwall_type,nwall_type):nwall_type=1,nWall_type");

	}

	// write surface-fluid interaction 
/*
	// Sigma_ff
	theOutput.print("@ ");
	for(int icomp = 0; icomp < nComp; icomp++) {
	    for(int jcomp = 0; jcomp < nComp; jcomp++) {
		theOutput.print( Sigma_ff[icomp][jcomp] + " ");
	    }
	}
	theOutput.println("\t Sigma_ff[nComp][nComp] order [0][0], [0][1]...");

	// Eps_ff
	theOutput.print("@ ");
	for(int icomp = 0; icomp < nComp; icomp++) {
	    for(int jcomp = 0; jcomp < nComp; jcomp++) {
		theOutput.print( Eps_ff[icomp][jcomp] + " ");
	    }
	}
	theOutput.println("\t Eps_ff[nComp][nComp] order [0][0], [0][1]...");

	// Cut_ff
	theOutput.print("@ ");
	for(int icomp = 0; icomp < nComp; icomp++) {
	    for(int jcomp = 0; jcomp < nComp; jcomp++) {
		theOutput.print( Cut_ff[icomp][jcomp] + " ");
	    }
	}
	theOutput.println("\t Cut_ff[nComp][nComp] order [0][0], [0][1]...");

	// Sigma_wf
	theOutput.print("@ ");
	for(int icomp = 0; icomp < nComp; icomp++) {
	    for(int jWall = 0; jWall < nWall_type; jWall++) {
		theOutput.print( Sigma_wf[icomp][jWall] + " ");
	    }
	}
	theOutput.println("\t Sigma_wf[nComp][nWall_type] order [0][0], [0][1]...");

	// Eps_wf
	theOutput.print("@ ");
	for(int icomp = 0; icomp < nComp; icomp++) {
	    for(int jWall = 0; jWall < nWall_type; jWall++) {
		theOutput.print( Eps_wf[icomp][jWall] + " ");
	    }
	}
	theOutput.println("\t Eps_wf[nComp][nWall_type] order [0][0], [0][1]...");

	// Cut_wf
	theOutput.print("@ ");
	for(int icomp = 0; icomp < nComp; icomp++) {
	    for(int jWall = 0; jWall < nWall_type; jWall++) {
		theOutput.print( Cut_wf[icomp][jWall] + " ");
	    }
	}
	theOutput.println("\t Cut_wf[nComp][nWall_type] order [0][0], [0][1]...");


	// Sigma_ww
	theOutput.print("@ ");
	for(int iWall = 0; iWall < nWall_type; iWall++) {
	    for(int jWall = 0; jWall < nWall_type; jWall++) {
		theOutput.print( Sigma_ww[iWall][jWall] + " ");
	    }
	}
	theOutput.println("\t Sigma_ww[nWall_type][nWall_type] order [0][0], [0][1]...");

	// Eps_ww
	theOutput.print("@ ");
	for(int iWall = 0; iWall < nWall_type; iWall++) {
	    for(int jWall = 0; jWall < nWall_type; jWall++) {
		theOutput.print( Eps_ww[iWall][jWall] + " ");
	    }
	}
	theOutput.println("\t Eps_ww[nWall_type][nWall_type] order [0][0], [0][1]...");

	// Cut_ww
	theOutput.print("@ ");
	for(int iWall = 0; iWall < nWall_type; iWall++) {
	    for(int jWall = 0; jWall < nWall_type; jWall++) {
	    theOutput.print( Cut_ww[iWall][jWall] + " ");
	    }
	}
	theOutput.println("\t Cut_ww[nWall_type][nWall_type] order [0][0], [0][1]...");
*/
    }

    private void writeNonLinearSolverParameters(PrintWriter theOutput) {
	theOutput.println("********* NONLINEAR SOLVER PARAMETERS *****************************************");
	
	theOutput.println("@ " +
			  maxNewtonIter +
			  "\tMaximum # of Newton Iterations");
	theOutput.println("@ " + newtonRelTol +
			  " " + newtonAbsTol +
			  "\tRelative and Absolute convergence tolerances");
	theOutput.println("@ " + loadBalanceSwitch +
			  "\tLoad balance switch (0=linear, 1=box, 2=weights, 3=timings)");
    }

    private void writeOutputFormatParameters(PrintWriter theOutput) {
	theOutput.println("********* OUTPUT FORMAT PARAMETERS ********************************************");
	
	theOutput.println("@ " + lPerArea + "  " + lPrintGr +
			  "\tLper_area lPrint_Gr (radial distribution fxn)");
	theOutput.println("@ " + printRhoType +
			  "\tPrint_rho_type");
	theOutput.println("@ " + printRhoSwitch +
			  " " + printMeshSwitch +
			  "\tPrint_rho_switch Print_mesh_switch");
	theOutput.println("@ " + printHeader + 
			  "\tPrint_header: logical for printing header in the dft_force.dat files 0=false, 1=true");
	theOutput.println("@ " + printIWrite +
			  "\tIWRITE (0=Minimal, 1=Density_Prof, 3=Verbose)");
    }
    
    private void writePolymerInputParameters(PrintWriter theOutput) {
	theOutput.println("********* POLYMER INPUT PARAMETERS ********************************************");
	theOutput.println("@ " + nPol_comp +
			  "\tNpol_comp: Number of (co)polymer components");

	theOutput.print("@ ");
	for(int p=0; p<polymerVector.size(); p++) {
	    theOutput.print( ((polymerBlockInfo) polymerVector.get(p)).nBlock + " ");
	}
	theOutput.println("\tNblock[pol_num]: Number of blocks in each copolymer");

	theOutput.print("@ ");
	for(int p=0; p<polymerVector.size(); p++) {
	    int nBlock =  ((polymerBlockInfo) polymerVector.get(p)).nBlock;
	    int block[] = ((polymerBlockInfo) polymerVector.get(p)).block;
	    for(int b=0; b<nBlock; b++) {
		theOutput.print(  block[b] + " ");
	    }
	}
	theOutput.println("\tblock[pol_num][iblock]: Number of segments in cach block");

	theOutput.print("@ ");
	for(int p=0; p<polymerVector.size(); p++) {
	    int nBlock =  ((polymerBlockInfo) polymerVector.get(p)).nBlock;
	    int block_type[] = ((polymerBlockInfo) polymerVector.get(p)).block_type;
	    for(int b=0; b<nBlock; b++) {
		theOutput.print(  block_type[b] + " ");
	    }
	}
	theOutput.println("\tblock_type[pol_num][iblock]: Segment types in each block");

	theOutput.print("@ ");
	for(int p=0; p<polymerVector.size(); p++) {
	    theOutput.print( ((polymerBlockInfo) polymerVector.get(p)).poly_file + " ");
	}
	theOutput.println("\tpoly_file: File of connectivity data");

	theOutput.print("@ ");
	theOutput.print( cr_file );
	theOutput.println("\tCr_file: c(r) filename");

	theOutput.print("@ ");
	for(int p=0; p<polymerVector.size(); p++) {
	    theOutput.print( ((polymerBlockInfo) polymerVector.get(p)).cr_radius + " ");
	}
	theOutput.println("\tCr_rad: c(r) radius (units of sigma)");

	theOutput.print("@ ");
	for(int p=0; p<polymerVector.size(); p++) {
	    theOutput.print( ((polymerBlockInfo) polymerVector.get(p)).gauss_a + " ");
	}
	theOutput.println("\tGauss_a: Aspect ratio (gauss bl/sigma)");

    }

    private void writePotentialTypeParameters(PrintWriter theOutput) {
	theOutput.println("********* POTENTIAL TYPE SELECTIONS *******************************************");
        theOutput.println( "@ " + ipot_wf_n + 
			   "\tIpot_wf_n (0=No_wall-fluid\n" + 
			   "\t           1=Hard_wall\n" +
			   "\t           2=LJ9_3\n" +
			   "\t           3=normalized LJ9_3\n" +
			   "\t           4=LJ12_6\n" +
			   "\t           5=LJ_ATOMIC: atomic surfaces with\n" +
			   "\t             LJ12-6 wf interactions\n" +
			   "\t           6=stepped 9_3 (for 2D or 3D problems)\n" +
			   "\t           7=HARD_EXP : hard wall + exponential\n" +
                           "\t           8=HARD_ATOMS hard surfaces AND some hard fixed atoms)\n" +
                           "\t           9=LJ9_3 wall plus LJ atoms like case 8 for the hard systems)" );
    }

    private void writeProblemType(PrintWriter theOutput) {
	theOutput.println("********* PROBLEM TYPE ********************************************************");
	theOutput.println("  Problem type: " + problemTypeSection.getTag() );
	theOutput.println("");
    }


    private void writeSemiPermeableSurfaceParameters(PrintWriter theOutput) {
	theOutput.println("********* SEMI-PERMEABLE SURFACE PARAMETERS ***********************************");

	theOutput.print("@ ");
	for( int wall=0; wall < nWall_type; wall++) {
	    for( int comp=0; comp < nComp; comp++) {
		theOutput.print( lSemiPerm[wall][comp] + " " );
	    }
	}
	theOutput.println("\tLsemiperm[iwall_type][icomp]: [0][0], [0][1]...");

	theOutput.print("@ ");
	for( int wall=0; wall < nWall_type; wall++) {
	    for( int comp=0; comp < nComp; comp++) {
		theOutput.print( vExtMembrane[wall][comp] + " " );
	    }
	}
	theOutput.println("\tVext_membrane[iwall_type][icomp]: [0][0], [0][1]...");
	
	theOutput.println("\t*Note for polymers: See note above. Again replace Ncomp with Nblock_tot.");

    }


    private void writeStartupControlParameters(PrintWriter theOutput) {
	theOutput.println("********* STARTUP CONTROL PARAMETERS ******************************************");
	
	theOutput.println("@ " +
			  iLiqVap +
			  "\tIliq_vap (-2=no coex -1=none, 1=W-V, 2=W-L, 3=L-V profiles)");
	theOutput.println("@ " +
			  iGuess + " " +
			  "\tIguess");
	theOutput.println("\tIguess choices ...");
	theOutput.println("\t[-3, -2, -1: rho; CONST_RHO, CONST_RHO_L, CONST_RHO_V");
	theOutput.println("\t  0, 1, 2: rho*exp(-Vext/kt); EXP_RHO, EXP_RHO_L, EXP_RHO_V");
	theOutput.println("\t  3, 4, 5: rho*theta fnc: STEP_RHO, STEP_RHO_L, STEP_RHO_V");
	theOutput.println("\t  6      : rho*theta fnc for liq-vap prof; STEP_LV");
	theOutput.println("\t  7      : LINEAR (interpolation between left and right");

	theOutput.println("@ " +
			  thickness +
			  "\tThickness for stepped and linear profiles");
	theOutput.println("@ " + 
			  restart +
			  "\tRestart (0=no, 1=yes, 2=yes, but w/ step function");
	theOutput.println("\t3=yes for densities but not elec. pot. or chem. pot.)");
	theOutput.println("\t(if 1: dft_dens.dat must exist)");
	
	theOutput.println("@ " + rhoMax + " Rho_max");

    }

    private void writeStatePointParameters(PrintWriter theOutput) {
	theOutput.println("********* STATE POINT PARAMETERS **********************************************");
	theOutput.print("@ ");
	for(int i=0; i<nComp; i++) {
	    theOutput.print(bulkDensity[i] + " ");
	}
	theOutput.println("\tRho_b[icomp], icomp=0, Ncomp-1 (or Npol_comp for polymers)");
	theOutput.println("    *Note for polymers:");
	theOutput.println("     Rho_b is indexed Npol_comp rather than Nblock_tot.  The code automatically");
	theOutput.println("     converts to the density of the different polymer segments.");
	theOutput.println("     For example: For an ABC triblock in solvent you enter Rho_b[0],Rho_b[1]");
	theOutput.println("     corresponding to the polymer density and the solvent density.  The");
	theOutput.println("     code converts them to Rho_b's[0-2] based on the first value, and Rho_b'[3]");
	theOutput.println("     based on the second entry.");
    }

    private void writeSteadyStateBoundaryConditionParameters(PrintWriter theOutput) {
	theOutput.println("********* STEADY STATE BOUNDARY CONDITION PARAMETERS **************************");
	
	theOutput.println("@ " + lSteadyState + 
			  " " + linearTransport +
			  "\tLsteady_state Linear_transport");
	theOutput.println("\tlSteadyState 0=equilibrium problem, 1=steady state problem");
	theOutput.println("\tLinear_transport 0=nonlinear, 1=linear");
	
	theOutput.println("@ " + gradientDirection +
			  "\t Grad_dim direction of gradient (0=x, 1=y, 2=z)");
	
	theOutput.println("@ " + xConstMu +
			  "\tx_const_mu distance at which chemical potential is constant");
	
	theOutput.println("@ " + geometryFlag +
			  " " + n1DPoreSegments +
			  "\tGeom_flag (0=unit area, 1=cyl pore, 2=vary pore) Nseg (num pore segments)");

	theOutput.print("@ ");
	for(int seg=0; seg<n1DPoreSegments; seg++ ) {
	    theOutput.print(pore1DRadiusLeft[seg] + " " + pore1DRadiusRight[seg] + " " + pore1DLength[seg] + " " );
	}
	theOutput.println("\tRadius_L, Radius_R, Length");

	theOutput.print("@ ");
	for(int comp=0; comp<nComp; comp++) {
	    theOutput.print(bulkDensityLBB[comp] + " ");
	}
	theOutput.println("\tRho_b_Left[Ncomp] B.C. on left or bottom or back");

	theOutput.print("@ ");
	for(int comp=0; comp<nComp; comp++) {
	    theOutput.print(bulkDensityRTF[comp] + " ");
	}
	theOutput.println("\tRho_b_Right[Ncomp] B.C. on right or top or front");

	theOutput.print("@ ");
	for(int comp=0; comp<nComp; comp++) {
	    theOutput.print(diffusionCoef[comp] + " ");
	}
	theOutput.println("\tD_coef[Ncomp] Diff Coef per component (cm^2/sec)");

	theOutput.println("@ " + elecPotLBB + 
			  " " + elecPotRTF + 
			  "\tElec_pot_L, Elec_Pot_R B.C. on elec. potential lbb or rtf");

	theOutput.println("@ " + centerOfMassVelocity +
			  "\tCenter of mass velocity" );
    }

    private void writeSurfaceNames(PrintWriter theOutput) {
	theOutput.println("********* SURFACE NAMES *******************************************************");
	theOutput.print("  Surface names: " );
	if( surfaceName.length >= 1 ) {
	    theOutput.print( surfaceName[0] );
	}
	for(int i=1; i < surfaceName.length; i++ ) {
	    theOutput.print( ", " + surfaceName[i] );
	}
	theOutput.println("");
	theOutput.println("");
    }


    private void writeSurfaceParameters(PrintWriter theOutput) {
	theOutput.println("********* SURFACE PARAMETERS **************************************************");
	theOutput.println("@ " + nWall_type + " " + nWall + " " + 
			  nLink + " " + autoCenterSurfInBox +
			  "\tNwall_type Nwall Nlink Lauto_center");
	theOutput.print("@ ");
	for(int i=0; i<nLink; i++) {
	    for(int j=0; j<nDim; j++) {
		theOutput.print(xTest_reflect_TF[i][j] + " ");
	    }
	}
	theOutput.println("\t Xtest_reflect_TF");

	for(int i=0; i<nWall_type; i++) {
	  theOutput.print(surfaceName[i] + " ");
	}
	theOutput.println("\t Name of each surface");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	  theOutput.print(surf_type[i] + " ");
	}
	theOutput.println("\t Surf_type[iwall_type]; iwall_type=0,Nwall_type");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	  theOutput.print(orientation[i] + " ");
	}
	theOutput.println("\t Orientation[iwall_type]; iwall_type=0,Nwall_type");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	  theOutput.print(wallParam[i] + " ");
	}
	theOutput.println("\t WallParam[iwall_type]; iwall_type=0,Nwall_type");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	  theOutput.print(wallParam2[i] + " ");
	}
	theOutput.println("\t WallParam2[iwall_type]; iwall_type=0,Nwall_type");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	  theOutput.print(wallParam3[i] + " ");
	}
	theOutput.println("\t WallParam3[iwall_type]; iwall_type=0,Nwall_type");

    }

    private void writeSurfaceParticleParameters(PrintWriter theOutput) {
	theOutput.println("********** SURFACE PARTICLE PARAMETERS **********");
	theOutput.println("\t***** if Ipot_wf_n > 0 *****");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	    theOutput.print(sigma_w[i] + " ");
	}
	theOutput.println("\tSigma_w(iwall_type): iwall_type=1,Nwall_type");

	theOutput.println("\t***** if Ipot_wf_n > 1 *****");
	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	    theOutput.print(eps_w[i] + " ");
	}
	theOutput.println("\tEps_w(iwall_type): iwall_type=1,Nwall_type");

	theOutput.print("@ ");
	for(int i=0; i<nWall_type; i++) {
	    theOutput.print(rho_s[i] + " ");
	}
	theOutput.println("\tRho_s(iwall_type): iwall_type=1,Nwall_type");

    }

}
