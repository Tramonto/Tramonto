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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Optika_GUI.hpp"


void dft_GUI_mesh(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                  Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Mesh_List); 

void dft_GUI_functionals(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List);

void dft_GUI_surfaces(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceGeometry_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWW_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceParamCharge_List); 

void dft_GUI_vextType(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List); 

void dft_GUI_potentialsFF(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> PotentialsFF_List);

void dft_GUI_potentialsWF(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                      Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                      Teuchos::RCP<Teuchos::ParameterList> Surface_List,
                      Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> PotentialsWF_List); 

void dft_GUI_StatePoint( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                      Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                      Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                      Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                      Teuchos::RCP<Teuchos::ParameterList> StatePoint_List,
                      Teuchos::RCP<Teuchos::ParameterList> Diffusion_List,
                      Teuchos::RCP<Teuchos::ParameterList> ChargedFluid_List);

void dft_GUI_toTramonto(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Mesh_List,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List, 
                         Teuchos::RCP<Teuchos::ParameterList> Surface_List, 
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceGeometry_List); 

void dft_GUI_OutputParams(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                         Teuchos::RCP<Teuchos::ParameterList> SurfaceInteraction_List,
                         Teuchos::RCP<Teuchos::ParameterList> Output_List);

void dft_GUI_NumericalMethods(Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,
                         Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                         Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                         Teuchos::RCP<Teuchos::ParameterList> Solver_List,
                         Teuchos::RCP<Teuchos::ParameterList> Coarsening_List,
                         Teuchos::RCP<Teuchos::ParameterList> NonlinearSolver_List,
                         Teuchos::RCP<Teuchos::ParameterList> LinearSolver_List);

void dft_GUI_Polymer( Teuchos::RCP<Teuchos::ParameterList> Tramonto_List,   
                  Teuchos::RCP<Teuchos::DependencySheet> depSheet_Tramonto,
                  Teuchos::RCP<Teuchos::ParameterList> Functional_List,
                  Teuchos::RCP<Teuchos::ParameterList> Fluid_List,
                  Teuchos::RCP<Teuchos::ParameterList> Polymer_List,
                  Teuchos::RCP<Teuchos::ParameterList> PolymerCMS_List);

