/*
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
/ You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

#include <mpi.h>
/*#include "Optika_StandardConditions.hpp"
#include "Optika_StandardDependencies.hpp"
#include "Optika_DependencySheet.hpp"
*/
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"
/*#include "Teuchos_Version.hpp"*/
#include "Optika_GUI.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_StandardDependencies.hpp"



void dft_GUI_mesh_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List); 
void dft_GUI_mesh_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List); 
void dft_GUI_mesh_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List); 

void dft_GUI_functionals_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List);
void dft_GUI_functionals_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List);
void dft_GUI_functionals_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List);

void dft_GUI_potentialsFF_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List);
void dft_GUI_potentialsFF_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List);
void dft_GUI_potentialsFF_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List);

void dft_GUI_Polymer_set_defaults( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,   
                  Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List);
void dft_GUI_Polymer_set_OldFormat( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,   
                  Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List);
void dft_GUI_Polymer_dependencies( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,   
                  Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List);

void dft_GUI_StatePoint_set_defaults( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                      Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                      Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                      Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List);
void dft_GUI_StatePoint_set_OldFormat( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                      Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                      Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                      Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List);
void dft_GUI_StatePoint_dependencies( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                      Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                      Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                      Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List);

void dft_GUI_surfaces_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List);
void dft_GUI_surfaces_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List);
void dft_GUI_surfaces_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List);

void dft_GUI_surface_geometry_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfGeom_List);
void dft_GUI_surface_geometry_set_OldFormat(int iwall_type,Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfGeom_List);
void dft_GUI_surface_geometry_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfGeom_List);

void dft_GUI_potentialsWW_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List); 
void dft_GUI_potentialsWW_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List); 
void dft_GUI_potentialsWW_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List,
                      Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List); 

void dft_GUI_vextType_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List); 
void dft_GUI_vextType_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List); 
void dft_GUI_vextType_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List); 


void dft_GUI_potentialsWF_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List); 
void dft_GUI_potentialsWF_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List); 
void dft_GUI_potentialsWF_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List); 

void dft_GUI_Continuation_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                     Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                     Teuchos::RCP<Teuchos::ParameterList> Continuation_List);
void dft_GUI_Continuation_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                     Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                     Teuchos::RCP<Teuchos::ParameterList> Continuation_List);
void dft_GUI_Continuation_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                     Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                     Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List,
                     Teuchos::RCP<Teuchos::ParameterList> Continuation_List);

void dft_GUI_toTramonto(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Polymer_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PolymerGraft_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PolymerArch_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List, 
                         Teuchos::RCP<Teuchos::ParameterList> StatePoint_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Diffusion_List, 
                         Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Continuation_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List, 
                         Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PhysicsMethod_List, 
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List, 
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Output_List, 
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Surface_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfacePosition_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List, 
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom0_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom1_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom2_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom3_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom4_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom5_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom6_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom7_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom8_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom9_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom10_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfGeom11_List); 

void dft_GUI_OutputParams_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Output_List);
void dft_GUI_OutputParams_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Output_List);
void dft_GUI_OutputParams_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Continuation_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List,
                         Teuchos::RCP<Teuchos::ParameterList> Output_List);

void dft_GUI_DensityStartupParams_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List);
void dft_GUI_DensityStartupParams_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List);
void dft_GUI_DensityStartupParams_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Continuation_List,
                         Teuchos::RCP<Teuchos::ParameterList> DensProfile_List);

void dft_GUI_NumericalMethods_set_defaults(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                         Teuchos::RCP<Teuchos::ParameterList> PhysicsMethods_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List);
void dft_GUI_NumericalMethods_set_OldFormat(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                         Teuchos::RCP<Teuchos::ParameterList> PhysicsMethods_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List);
void dft_GUI_NumericalMethods_dependencies(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Optika::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> LoadBalance_List,
                         Teuchos::RCP<Teuchos::ParameterList> PhysicsMethods_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List);


