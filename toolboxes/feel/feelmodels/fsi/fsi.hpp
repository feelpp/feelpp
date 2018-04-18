/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/solid/solidmechanics.hpp>

#include <feel/feelmodels/fsi/interpolationfsi.hpp>
#include <feel/feelmodels/fsi/aitkenrelaxationfsi.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
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
    typedef FixPointConvergenceFSI<solid_type> fixpointconvergenceFSI_type;
    typedef boost::shared_ptr<fixpointconvergenceFSI_type> fixpointconvergenceFSI_ptrtype;

    //---------------------------------------------------------------------------------------------------------//

    FSI( std::string const& prefix, WorldComm const& _worldComm = Environment::worldComm(),
         std::string const& rootRepository = "" );
    FSI( self_type const & M ) = default;

    static std::string expandStringFromSpec( std::string const& expr );

    //---------------------------------------------------------------------------------------------------------//

    double meshSize() const { return M_meshSize; }

    fluid_ptrtype const& fluidModel() const { return M_fluidModel; }
    solid_ptrtype const& solidModel() const { return M_solidModel; }
    void setFluidModel( fluid_ptrtype const& fm ) { M_fluidModel=fm; }
    void setSolidModel( solid_ptrtype const& sm ) { M_solidModel=sm; }

    std::string fsiCouplingType() const { return M_fsiCouplingType; }
    std::string fsiCouplingBoundaryCondition() const { return M_fsiCouplingBoundaryCondition; }
    bool useFSISemiImplicitScheme() const { return ( this->fsiCouplingType() == "Semi-Implicit" ); }
    bool interfaceFSIisConforme() const { return M_interfaceFSIisConforme; }
    double fixPointTolerance() const { return M_fixPointTolerance; }
    double fixPointInitialTheta() const { return M_fixPointInitialTheta; }
    double fixPointMinTheta() const { return M_fixPointMinTheta; }
    int fixPointMaxIt() const { return M_fixPointMaxIt; }
    int fixPointMinItConvergence() const { return M_fixPointMinItConvergence; }

    interpolationFSI_ptrtype interpolationTool() { return M_interpolationFSI; }
    interpolationFSI_ptrtype const& interpolationTool() const { return M_interpolationFSI; }
    aitkenrelaxationFSI_ptrtype aitkenRelaxTool() { return M_aitkenFSI; }
    aitkenrelaxationFSI_ptrtype const& aitkenRelaxTool() const { return M_aitkenFSI; }

    //---------------------------------------------------------------------------------------------------------//

    boost::shared_ptr<std::ostringstream> getInfo() const;

    //---------------------------------------------------------------------------------------------------------//

    void createMesh();
    void init();
    void solve();

    //---------------------------------------------------------------------------------------------------------//

    void updateTime(double time);

    boost::shared_ptr<TSBase> timeStepBase() const { return this->fluidTimeStepBase(); }
    boost::shared_ptr<TSBase> fluidTimeStepBase() const { return this->fluidModel()->timeStepBase(); }
    boost::shared_ptr<TSBase> solidTimeStepBase() const { return this->solidModel()->timeStepBase(); }
    void updateTimeStep();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time )
    {
        this->fluidModel()->exportResults(time);
        this->solidModel()->exportResults(time);
    }

    //---------------------------------------------------------------------------------------------------------//
    void updateLinearPDE_Fluid( DataUpdateLinear & data ) const;
    void updateJacobian_Fluid( DataUpdateJacobian & data ) const;
    void updateResidual_Fluid( DataUpdateResidual & data ) const;
    //void updateLinearPDEStrongDirichletBC_Fluid( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const;
    //void updateJacobianStrongDirichletBC_Fluid( sparse_matrix_ptrtype& J,vector_ptrtype& RBis ) const;


private :
    void updateBackendOptimisation( int iterationFSI, double lastErrorRelative );
    void solveImpl1();
    void solveImpl2();
    void solveImpl3();

    //---------------------------------------------------------------------------------------------------------//

private :

    fluid_ptrtype M_fluidModel;
    solid_ptrtype M_solidModel;

    double M_meshSize;
    fs::path M_mshfilepathFluidPart1,M_mshfilepathSolidPart1;
    fs::path M_mshfilepathFluidPartN,M_mshfilepathSolidPartN;
    std::set<std::string> M_markersNameFluid,M_markersNameSolid;
    std::string M_tagFileNameMeshGenerated;

    std::string M_fsiCouplingType; // implicit,semi-implicit
    std::string M_fsiCouplingBoundaryCondition; // dirichlet-neumann, robin-robin, ...
    bool M_interfaceFSIisConforme;
    double M_fixPointTolerance, M_fixPointInitialTheta, M_fixPointMinTheta;
    int M_fixPointMaxIt, M_fixPointMinItConvergence;

    interpolationFSI_ptrtype M_interpolationFSI;
    aitkenrelaxationFSI_ptrtype M_aitkenFSI;
    fixpointconvergenceFSI_ptrtype M_fixPointConvergenceFSI;

    int M_previousTimeOrder,M_currentTimeOrder;
    bool M_reusePrecOptFluid,M_reusePrecRebuildAtFirstFSIStepOptFluid,M_reuseJacOptFluid,M_reuseJacRebuildAtFirstNewtonStepOptFluid,M_reuseJacRebuildAtFirstFSIStepOptFluid;
    bool M_reusePrecOptSolid,M_reusePrecRebuildAtFirstFSIStepOptSolid,M_reuseJacOptSolid,M_reuseJacRebuildAtFirstNewtonStepOptSolid,M_reuseJacRebuildAtFirstFSIStepOptSolid;
    int M_reusePrecActivatedAfterNbFsiIterationFluid,M_reusePrecActivatedAfterNbFsiIterationSolid;
    double M_reusePrecActivatedToleranceFluid,M_reusePrecActivatedToleranceSolid;

    double M_couplingNitscheFamily_gamma, M_couplingNitscheFamily_gamma0, M_couplingNitscheFamily_alpha;
    double M_couplingRNG_manualScaling;
    bool M_couplingRNG_useInterfaceOperator;

};

} // namespace FeelModels
} // namespace Feel



#endif // FEELPP_MODELS_FSI_H
