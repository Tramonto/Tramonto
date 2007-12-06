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

dft_eqLoader_NewPhys::dft_eqLoader_NewPhys(int iBlock, int numBlocks, int block2Phys[NEQ_TYPE])
   : iBlock_(iBlock), numBlocks_(numBlocks), block2Phys_(block2Phys)
{}

dft_eqLoader_NewPhys::~dft_eqLoader_NewPhys() {}


void dft_eqLoader_NewPhys::loadDensity()
{
   printBlockInfo("Density");
}

void dft_eqLoader_NewPhys::loadRhoBarRosen()
{
   printBlockInfo("RhoBarRosen");
}

void dft_eqLoader_NewPhys::loadPoisson()
{
   printBlockInfo("Poisson");
}

void dft_eqLoader_NewPhys::loadDiffusion()
{
   printBlockInfo("Diffusion");
}

void dft_eqLoader_NewPhys::loadCavityWTC()
{
   printBlockInfo("CavityWTC");
}

void dft_eqLoader_NewPhys::loadBondWTC()
{
   printBlockInfo("BondWTC");
}

void dft_eqLoader_NewPhys::loadCMSField()
{
   printBlockInfo("CMSField");
}

void dft_eqLoader_NewPhys::loadGCHAIN()
{
   printBlockInfo("GCHAIN");
}

void dft_eqLoader_NewPhys::loadNewPhys()
{
   printBlockInfo("NewPhys");
}

void dft_eqLoader_NewPhys::printBlockInfo(string phys)
{
   cout << "Loading Block (" << iBlock_ << ", " << jBlock_ << ") with d(NewPhys)/d(" 
        << phys << ")." << endl;
}

#endif
