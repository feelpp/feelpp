/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-07-05

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file fsi.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-05
 */

#ifndef FEELPP_MODELS_FSI_H
#define FEELPP_MODELS_FSI_H 1

#include <feel/feelmodels2/fluid/fluidmechanics.hpp>
#include <feel/feelmodels2/solid/solidmechanics.hpp>

#include <feel/feelmodels2/fsi/interpolationfsi.hpp>
#include <feel/feelmodels2/fsi/aitkenrelaxationfsi.hpp>
//#include <feel/feelmodels2/modelcore/modelbase.hpp>
//#include <feel/feelmodels2/modelcore/modelalgebraic.hpp>
#include <feel/feelmodels2/modelcore/modelnumerical.hpp>
#include <feel/feelts/tsbase.hpp>

namespace Feel
{
namespace FeelModels
{


template< class FluidType, class SolidType >
class FSI : public ModelNumerical
{
public :
    typedef ModelNumerical super_type;
    typedef FSI<FluidType,SolidType> self_type;

    typedef FluidType fluid_type;
    typedef SolidType solid_type;
    typedef boost::shared_ptr<fluid_type> fluid_ptrtype;
    typedef boost::shared_ptr<solid_type> solid_ptrtype;

    typedef InterpolationFSI<fluid_type,solid_type> interpolationFSI_type;
    typedef boost::shared_ptr<interpolationFSI_type> interpolationFSI_ptrtype;

    typedef AitkenRelaxationFSI<solid_type> aitkenrelaxationFSI_type;
    typedef boost::shared_ptr<aitkenrelaxationFSI_type> aitkenrelaxationFSI_ptrtype;

    //---------------------------------------------------------------------------------------------------------//

    FSI( std::string prefix, WorldComm const& _worldComm = Environment::worldComm() );
    FSI( self_type const & M ) = default;

    //---------------------------------------------------------------------------------------------------------//

    double meshSize() const { return M_meshSize; }

    fluid_ptrtype const& fluidAppli() const { return M_fluid; }
    solid_ptrtype const& solidAppli() const { return M_solid; }
    void setFluidModel( fluid_ptrtype const& fm ) { M_fluid=fm; }
    void setSolidModel( solid_ptrtype const& sm ) { M_solid=sm; }

    //bool verbose() const { return M_verbose; }
    //bool verboseAllProc() const { return M_verboseAllProc; }

    std::string fsiCouplingType() const { return M_fsiCouplingType; }
    std::string fsiCouplingBoundaryCondition() const { return M_fsiCouplingBoundaryCondition; }
    bool interfaceFSIisConforme() const { return M_interfaceFSIisConforme; }
    double tolPtFixe() const { return M_tolPtFixe; }
    double initialTheta() const { return M_initialTheta; }
    double minTheta() const { return M_minTheta; }
    int fixPointMaxIt() const { return M_fixPointMaxIt; }

    interpolationFSI_ptrtype interpolationTool() { return M_interpolationFSI; }
    interpolationFSI_ptrtype const& interpolationTool() const { return M_interpolationFSI; }
    aitkenrelaxationFSI_ptrtype aitkenRelaxTool() { return M_aitkenFSI; }
    aitkenrelaxationFSI_ptrtype const& aitkenRelaxTool() const { return M_aitkenFSI; }

    //---------------------------------------------------------------------------------------------------------//

    boost::shared_ptr<std::ostringstream> getInfo() const;
#if 0
    void
    printInfo() const
    {
        this->fluidAppli()->printInfo();
        this->solidAppli()->printInfo();

        if ( M_fluid->verboseAllProc() || M_solid->verboseAllProc() ) std::cout << this->getInfo()->str();
        else if (M_fluid->worldComm().globalRank()==M_fluid->worldComm().masterRank() )
            std::cout << this->getInfo()->str();
    }

    //---------------------------------------------------------------------------------------------------------//

    void
    saveInfo() const
    {
        this->fluidAppli()->saveInfo();
        this->solidAppli()->saveInfo();

        if ( this->worldComm().isMasterRank() )
        {
            std::string nameFile = "FSI.info";// prefixvm(this->prefix(),"FSI.info");
            std::ofstream file(nameFile.c_str(), std::ios::out);
            file << this->getInfo()->str();
            file.close();
        }
    }

    //---------------------------------------------------------------------------------------------------------//

    void printAndSaveInfo() const { this->printInfo();this->saveInfo(); }
#endif

    //---------------------------------------------------------------------------------------------------------//

    void createMesh();
    void init();
    void solve();

    //---------------------------------------------------------------------------------------------------------//

    void
    updateTime(double time)
    {
        M_fluid->updateTime(time);
        M_solid->updateTime(time);
    }

    double currentTime() const { return M_fluid->time(); }
    double time() const { return this->currentTime(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->fluidTimeStepBase(); }
    boost::shared_ptr<TSBase> fluidTimeStepBase() const { return this->fluidAppli()->timeStepBase(); }
    boost::shared_ptr<TSBase> solidTimeStepBase() const { return this->solidAppli()->timeStepBase(); }
    void updateTimeStep();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time )
    {
        this->fluidAppli()->exportResults(time);
        this->solidAppli()->exportResults(time);
    }

private :
    void updateBackendOptimisation( bool restartFullStepFluid,bool restartFullStepSolid );
    void solveImpl1();
    void solveImpl2();

    //---------------------------------------------------------------------------------------------------------//

private :

    fluid_ptrtype M_fluid;
    solid_ptrtype M_solid;

    double M_meshSize;
    fs::path M_mshfilepathFluidPart1,M_mshfilepathSolidPart1;
    fs::path M_mshfilepathFluidPartN,M_mshfilepathSolidPartN;
    std::set<std::string> M_markersNameFluid,M_markersNameSolid;
    std::string M_tagFileNameMeshGenerated;

    std::string M_fsiCouplingType; // implicit,semi-implicit
    std::string M_fsiCouplingBoundaryCondition; // dirichlet-neumann, robin-neumann
    bool M_interfaceFSIisConforme;
    double M_tolPtFixe,M_initialTheta,M_minTheta;
    int M_fixPointMaxIt;

    interpolationFSI_ptrtype M_interpolationFSI;
    aitkenrelaxationFSI_ptrtype M_aitkenFSI;

    int M_previousTimeOrder,M_currentTimeOrder;
    bool M_reusePrecOptFluid,M_reusePrecRebuildAtFirstFSIStepOptFluid,M_reuseJacOptFluid,M_reuseJacRebuildAtFirstNewtonStepOptFluid,M_reuseJacRebuildAtFirstFSIStepOptFluid;
    bool M_reusePrecOptSolid,M_reusePrecRebuildAtFirstFSIStepOptSolid,M_reuseJacOptSolid,M_reuseJacRebuildAtFirstNewtonStepOptSolid,M_reuseJacRebuildAtFirstFSIStepOptSolid;
};

} // namespace FeelModels
} // namespace Feel



#endif // FEELPP_MODELS_FSI_H
