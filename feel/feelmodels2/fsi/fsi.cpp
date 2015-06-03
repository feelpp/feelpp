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
   \file fsi.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-05
 */

#include <feel/feelmodels2/fsi/fsi.hpp>

namespace Feel
{
namespace FeelModels
{

template< class FluidType, class SolidType >
FSI<FluidType,SolidType>::FSI(std::string prefix,WorldComm const& worldComm )
    :
    super_type( prefix, worldComm ),
    //M_fluid(fluid),
    //M_solid(solid),
    //M_verbose(fluid->verbose() || solid->verbose()),
    //M_verboseAllProc(fluid->verboseAllProc() || solid->verboseAllProc()),
    M_fsiCouplingType( option(_name="fsi.coupling-type").template as<std::string>() ),
    M_fsiCouplingBoundaryCondition( option(_name="fsi.coupling-bc").template as<std::string>() ),
    M_interfaceFSIisConforme( option(_name="fsi.conforming-interface").template as<bool>() ),
    M_tolPtFixe( option(_name="fsi.fixpoint.tol").template as<double>() ),
    M_initialTheta( option(_name="fsi.fixpoint.initialtheta").template as<double>() ),
    M_minTheta( option(_name="fsi.fixpoint.min_theta").template as<double>() ),
    M_fixPointMaxIt( option(_name="fsi.fixpoint.maxit").template as<int>() ),
    M_previousTimeOrder(0),M_currentTimeOrder(1),
    M_reusePrecOptFluid(false),
    M_reusePrecRebuildAtFirstFSIStepOptFluid( option(_name="fsi.fluid.reuse-prec.rebuild-at-first-fsi-step").template as<bool>() ),
    M_reuseJacOptFluid(false),
    M_reuseJacRebuildAtFirstNewtonStepOptFluid(false),
    M_reuseJacRebuildAtFirstFSIStepOptFluid( option(_name="fsi.fluid.reuse-jac.rebuild-at-first-fsi-step").template as<bool>() ),
    M_reusePrecOptSolid(false),
    M_reusePrecRebuildAtFirstFSIStepOptSolid( option(_name="fsi.solid.reuse-prec.rebuild-at-first-fsi-step").template as<bool>() ),
    M_reuseJacOptSolid(false),
    M_reuseJacRebuildAtFirstNewtonStepOptSolid(false),
    M_reuseJacRebuildAtFirstFSIStepOptSolid( option(_name="fsi.solid.reuse-jac.rebuild-at-first-fsi-step").template as<bool>() )
{
    this->log("fsi","constructor","start-finish");
}
//---------------------------------------------------------------------------------------------------------//

namespace detail
{

template <typename FluidType,typename SolidType>
typename SolidType::mesh_1dreduced_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, mpl::bool_<false> /**/ )
{
    auto submeshStruct = createSubmesh( FM->getMeshALE()->referenceMesh(), markedfaces( FM->getMeshALE()->referenceMesh(), FM->markersNameMovingBoundary()/*"Paroi"*/) );
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
      submeshStruct->faces().modify( submeshStruct->faceIterator( itp->id() ) , Feel::detail::UpdateMarker( submeshStruct->markerName("Fixe") ) );

    typedef SubMeshData smd_type;
    typedef boost::shared_ptr<smd_type> smd_ptrtype;
    smd_ptrtype smd( new smd_type(FM->mesh()) );
    for ( auto const& e : elements(submeshStruct) )
      {
        auto const& theface = FM->getMeshALE()->referenceMesh()->face( submeshStruct->subMeshToMesh(e.id()) );
        size_type idElt2 = FM->getMeshALE()->dofRelationShipMap()->geoElementMap()[ theface.element0().id() ].first;
        //std::cout << " e.G() " << e.G() << " other.G() " <<  theface.G() << std::endl;
        auto const& theface2 = FM->mesh()->element(idElt2,e.processId()).face(theface.pos_first());
        smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), theface2.id() ) );
      }
    submeshStruct->setSubMeshData( smd );

    return submeshStruct;
}

template <typename FluidType,typename SolidType>
typename FluidType::mesh_type::trace_mesh_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, mpl::bool_<true> /**/ )
{
    auto submeshStruct = createSubmesh( FM->mesh(), markedfaces( FM->mesh(),FM->markersNameMovingBoundary()/*"Paroi"*/) );
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
      submeshStruct->faces().modify( submeshStruct->faceIterator( itp->id() ) , Feel::detail::UpdateMarker( submeshStruct->markerName("Fixe") ) );

    return submeshStruct;
}

template <typename FluidType,typename SolidType>
typename SolidType::mesh_1dreduced_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM )
{
  static const bool hasSameOrderGeo = FluidType::mesh_type::nOrder == SolidType::mesh_1dreduced_type::nOrder;
  return createMeshStruct1dFromFluidMesh<FluidType,SolidType>(FM, mpl::bool_<hasSameOrderGeo>() );
}

} // namespace detail

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::init()
{
    this->log("ToolBoxFSI","init","start");

    if ( !M_fluid )
        M_fluid = fluid_ptrtype( new fluid_type("fluid",true,this->worldComm() ) );
    if ( !M_solid )
    {
        bool doExtractSubmesh = boption(_name="solid-mesh.extract-1d-from-fluid-mesh",_prefix=this->prefix() );
        M_solid = solid_ptrtype( new solid_type("solid",!doExtractSubmesh,this->worldComm() ) );
        if ( doExtractSubmesh )
        {
            CHECK( !M_fluid->markersNameMovingBoundary().empty() ) << "no marker moving boundary in fluid model";

            //std::string optionMarkerSubmesh = prefixvm(this->prefix(),"solid-mesh.extract-1d-from-fluid-mesh.boundary-markers" );
            //CHECK( Environment::vm().count( optionMarkerSubmesh ) ) << "option " << optionMarkerSubmesh << "must also given";

            //std::vector<std::string> configFiles = Environment::vm()[optionMarkerSubmesh.c_str()].template as<std::vector<std::string> >();
            //CHECK(false ) << "TODO";

            auto submeshStruct = detail::createMeshStruct1dFromFluidMesh<fluid_type,solid_type>(M_fluid);
            M_solid->loadMesh(submeshStruct);

        }

    }

    //-------------------------------------------------------------------------//
    // set coupling parameters
    if ( this->fsiCouplingBoundaryCondition()=="dirichlet-neumann" )
    {
        M_fluid->couplingFSIcondition("dirichlet");
        M_solid->couplingFSIcondition("neumann");
    }
    else if (this->fsiCouplingBoundaryCondition()=="robin-neumann")
    {
        M_fluid->couplingFSIcondition("robin");
        M_solid->couplingFSIcondition("neumann");
    }
    else if (this->fsiCouplingBoundaryCondition()=="robin-robin")
    {
        M_fluid->couplingFSIcondition("robin");
        M_solid->couplingFSIcondition("robin");
    }

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        M_fluid->useFSISemiImplicitScheme(true);
        M_solid->useFSISemiImplicitScheme(true);
    }

    M_solid->muFluidFSI( M_fluid->densityViscosityModel()->cstMu() );
    double gammaNitsche = option(_name="fsi.coupling-robin-robin.gamma").template as<double>();
    M_fluid->gammaNitschFSI( gammaNitsche );
    M_solid->gammaNitschFSI( gammaNitsche );
    double gamma0Nitsche = option(_name="fsi.coupling-robin-robin.gamma0").template as<double>();
    M_fluid->gamma0NitschFSI( gamma0Nitsche );

    //-------------------------------------------------------------------------//
    // save if reuse prec option at the begining
    M_reusePrecOptFluid = M_fluid->backend()->reusePrec();
    M_reuseJacOptFluid = M_fluid->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptFluid = M_fluid->backend()->reuseJacRebuildAtFirstNewtonStep();
    M_reusePrecOptSolid = M_solid->backend()->reusePrec();
    M_reuseJacOptSolid = M_solid->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptSolid = M_solid->backend()->reuseJacRebuildAtFirstNewtonStep();
    //-------------------------------------------------------------------------//

    // call init
    M_fluid->init();
    M_solid->init();
    //-------------------------------------------------------------------------//
    //M_solid->createAdditionalFunctionSpacesFSI();
    //-------------------------------------------------------------------------//
    // init interpolation tool
    M_interpolationFSI.reset(new interpolationFSI_type(M_fluid,M_solid));
    //-------------------------------------------------------------------------//
    // init aitken relaxation tool
    //AitkenType aitkenType = AitkenType::AITKEN_METHOD_1;
    std::string aitkenType = "method1";
    if (this->fsiCouplingBoundaryCondition()=="robin-robin" || this->fsiCouplingBoundaryCondition()=="robin-neumann")
    {
        //aitkenType = AitkenType::FIXED_RELAXATION_METHOD;
        aitkenType = "fixed-relaxation";
        M_initialTheta=0;
    }
    M_aitkenFSI.reset(new aitkenrelaxationFSI_type(M_solid,
                                                   aitkenType,//AITKEN_METHOD_1,
                                                   M_initialTheta,
                                                   M_tolPtFixe,
                                                   M_minTheta));
    //-------------------------------------------------------------------------//


    this->log("ToolBoxFSI","init","finish");
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solve()
{
    if ( this->fsiCouplingBoundaryCondition()=="dirichlet-neumann" )
        this->solveImpl1();
    else if (this->fsiCouplingBoundaryCondition()=="robin-robin" ||
             this->fsiCouplingBoundaryCondition()=="robin-neumann" )
        this->solveImpl2();
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateBackendOptimisation( bool restartFullStepFluid,bool restartFullStepSolid )
{
    // ici on impose la reconstruction du preconditionneur au debut du solveur FSI
    // maybe add an option to do that
    if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid) )
    {
        M_fluid->backend()->setReuseJacRebuildAtFirstNewtonStep(restartFullStepFluid);
    }

    if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
    {
        M_fluid->backend()->setReusePrec(!restartFullStepFluid);
    }

    if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
    {
        M_solid->backend()->setReuseJacRebuildAtFirstNewtonStep(restartFullStepSolid);
    }

    if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
    {
        M_solid->backend()->setReusePrec(!restartFullStepSolid);
    }
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl1()
{
    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluid->updateALEmesh();

        M_fluid->setRebuildLinearPartInJacobian(true);M_fluid->setRebuildCstPartInLinearSystem(true);
        M_solid->setRebuildLinearPartInJacobian(true);M_solid->setRebuildCstPartInLinearSystem(true);
        M_fluid->setRebuildCstPartInResidual(true);M_solid->setRebuildCstPartInResidual(true);
    }
#if 0
    // ici on impose la reconstruction du preconditionneur au debut du solveur FSI
    // maybe add an option to do that
    if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid) )
    {
        M_fluid->backend()->setReuseJacRebuildAtFirstNewtonStep(true);
    }

    if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
    {
        M_fluid->backend()->setReusePrec(false);
    }

    if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
    {
        M_solid->backend()->setReuseJacRebuildAtFirstNewtonStep(true);
    }

    if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
    {
        M_solid->backend()->setReusePrec(false);
    }
#else
    this->updateBackendOptimisation(true,true);
#endif
    // predictor disp
    M_solid->predictorDispl();
    // coupling fluid structure
    this->aitkenRelaxTool()->restart();
    //aitkenFSI.setTheta(initialTheta);

    boost::mpi::timer timerCur,timerIter;

    while (!this->aitkenRelaxTool()->isFinished() && this->aitkenRelaxTool()->nIterations()<this->fixPointMaxIt())
    {
        timerIter.restart();
        timerCur.restart();
        //--------------------------------------------------------------//

        this->aitkenRelaxTool()->saveOldSolution();

        //--------------------------------------------------------------//
        // ALE solver
        timerCur.restart();
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluid->updateALEmesh();
            double tALE = timerCur.elapsed();
            this->log("main","update ale","finish in "+(boost::format("%1% s") % tALE).str() );
        }
        else  if (this->fsiCouplingType()=="Semi-Implicit")
        {
            this->interpolationTool()->transfertVelocity();
            double tALE = timerCur.elapsed();
            this->log("main","transfert velocity","finish in "+(boost::format("%1% s") % tALE).str());
        }
        // robin-neumann as fsi condition
        if (this->fsiCouplingBoundaryCondition()=="robin-neumann")
        {
            this->interpolationTool()->transfertStressS2F();
        }

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        M_fluid->solve();
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//

        //--------------------------------------------------------------//
        timerCur.restart();
        // revert ref mesh
        M_fluid->getMeshALE()->revertReferenceMesh();
        // transfert stress
        this->interpolationTool()->transfertStress();
        // revert moving mesh
        M_fluid->getMeshALE()->revertMovingMesh();


        //if (this->fsiCouplingBoundaryCondition()=="robin-robin")
        //    this->interpolationTool()->transfertVelocityF2S(true);

        double t3 = timerCur.elapsed();
        this->log("main","transfert stress","finish in "+(boost::format("%1% s") % t3).str() );
        //--------------------------------------------------------------//

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        M_solid->solve();
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
#if 0
        M_fluid->updateTime(M_fluid->time()+1);
        M_fluid->exportResults();
        M_solid->updateTime(M_solid->time()+1);
        M_solid->exportResults();
#endif

        //--------------------------------------------------------------//
        timerCur.restart();
        //compute and apply aitken relaxation
        //if ( this->fsiCouplingBoundaryCondition()=="dirichlet-neumann" )
        this->aitkenRelaxTool()->applyRelaxation();
        // update velocity and acceleration
        M_solid->updateVelocity();
        // aitken relaxtion
        if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
            this->aitkenRelaxTool()->printInfo();
        this->aitkenRelaxTool()->shiftRight();

        double t4 = timerCur.elapsed();
        this->log("main","apply relax and up vel/acc struct","finish in "+(boost::format("%1% s") % t4).str() );
        //--------------------------------------------------------------//

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluid->setRebuildLinearPartInJacobian(false);M_fluid->setRebuildCstPartInLinearSystem(false);
            M_solid->setRebuildLinearPartInJacobian(false);M_solid->setRebuildCstPartInLinearSystem(false);
            M_fluid->setRebuildCstPartInResidual(false);M_solid->setRebuildCstPartInResidual(false);
        }
        /*else if (this->fsiCouplingType()=="Implicit" &&
         M_reuseJacOptFluid && !M_reuseJacRebuildAtFirstNewtonStepOptFluid &&
         !M_fluid->useLinearJacobianInResidual() )
         {
         M_fluid->setRebuildLinearPartInJacobian(false);
         }*/
#if 0
        if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid ) )
        {
            M_fluid->backend()->setReuseJacRebuildAtFirstNewtonStep(false);
        }
        if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
        {
            M_fluid->backend()->setReusePrec(true);
        }

        if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
        {
            M_solid->backend()->setReuseJacRebuildAtFirstNewtonStep(false);
        }
        if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
        {
            M_solid->backend()->setReusePrec(true);
        }
#else
        this->updateBackendOptimisation(false,false);
#endif
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("main","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
    }

}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl2()
{
    boost::mpi::timer timerIter;

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluid->updateALEmesh();

        M_fluid->setRebuildLinearPartInJacobian(true);M_fluid->setRebuildCstPartInLinearSystem(true);
        M_solid->setRebuildLinearPartInJacobian(true);M_solid->setRebuildCstPartInLinearSystem(true);
        M_fluid->setRebuildCstPartInResidual(true);M_solid->setRebuildCstPartInResidual(true);
    }


    this->updateBackendOptimisation(true,true);


#if 1
    auto oldSolDisp = M_solid->functionSpaceDisplacement()->elementPtr();
    auto residualDisp = M_solid->functionSpaceDisplacement()->elementPtr();
#endif
#if 0
    this->aitkenRelaxTool()->restart();
#endif
    int cptFSI=0;
    double residualConvergence=1;

    bool solveStruct = this->fluidAppli()->timeStepBDF()->iteration() > 1 || this->fixPointMaxIt()==1;

    while ( ( residualConvergence > this->tolPtFixe() || cptFSI <3 ) && cptFSI < this->fixPointMaxIt() )
    {
        timerIter.restart();
        //timerCur.restart();
#if 1
        *oldSolDisp = Feel::vf::project(_space=oldSolDisp->functionSpace(),
                                        _range=elements(oldSolDisp->mesh()),
                                        _expr=vf::idv(M_solid->fieldDisplacement()) );
#endif
#if 0
        this->aitkenRelaxTool()->saveOldSolution();
#endif

        //--------------------------------------------------------------//
        if (solveStruct)
        {
            M_fluid->getMeshALE()->revertReferenceMesh();
            // transfert stress
            this->interpolationTool()->transfertStress();
            // revert moving mesh
            M_fluid->getMeshALE()->revertMovingMesh();
            bool useExtrap = this->fluidAppli()->timeStepBDF()->iteration() > 2;
            if (this->fsiCouplingBoundaryCondition()=="robin-robin")
                this->interpolationTool()->transfertVelocityF2S(useExtrap);
            M_solid->solve();
            M_solid->updateVelocity();
        }

        //--------------------------------------------------------------//
        this->interpolationTool()->transfertVelocity();
        M_fluid->solve();
        //--------------------------------------------------------------//

#if 1
        *residualDisp = vf::project(residualDisp->functionSpace(),
                                    elements(residualDisp->mesh()),
                                    vf::idv(M_solid->fieldDisplacement() ));
        *residualDisp -= *oldSolDisp;
        auto oldEltL2Norm = oldSolDisp->l2Norm();
        if ( oldEltL2Norm > 1e-8 )
            residualConvergence = residualDisp->l2Norm()/oldEltL2Norm;
        else
            residualConvergence = residualDisp->l2Norm();
        if (M_fluid->worldComm().isMasterRank() )
            std::cout << " cptFSI " << cptFSI << " residualConvergence " << residualConvergence << "\n";
#endif
#if 0
        this->aitkenRelaxTool()->applyRelaxation();
        // update velocity and acceleration
        //M_solid->updateVelocity();
        // aitken relaxtion
        if (M_fluid->worldComm().isMasterRank())
            this->aitkenRelaxTool()->printInfo();
        this->aitkenRelaxTool()->shiftRight();
#endif

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluid->setRebuildLinearPartInJacobian(false);M_fluid->setRebuildCstPartInLinearSystem(false);
            M_fluid->setRebuildCstPartInResidual(false);

            if (solveStruct)
            {
                M_solid->setRebuildLinearPartInJacobian(false);M_solid->setRebuildCstPartInLinearSystem(false);
                M_solid->setRebuildCstPartInResidual(false);
            }
        }

        this->updateBackendOptimisation(false,!solveStruct);


        if (!solveStruct) solveStruct=true;

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("main","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");

        ++cptFSI;
    }

}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateTimeStep()
{
    M_previousTimeOrder=M_fluid->timeStepBDF()->timeOrder();

    M_fluid->updateTimeStep();
    M_solid->updateTimeStep();

    M_currentTimeOrder=M_fluid->timeStepBDF()->timeOrder();
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
boost::shared_ptr<std::ostringstream>
FSI<FluidType,SolidType>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||-----------------Info : FSI-------------------||"
           << "\n||==============================================||"
           << "\n   Coupling type : " << this->fsiCouplingType()
           << "\n   Coupling BC : " << this->fsiCouplingBoundaryCondition();
    if (M_interfaceFSIisConforme)
        *_ostr << "\n   Interface property : conformal";
    else
        *_ostr << "\n   Interface property : non conformal";

    *_ostr << "\n   Fix point parameters"
           << "\n     -- method : Aitken"
           << "\n     -- tolerance  : " << M_tolPtFixe
           << "\n     -- initial theta  : " << M_initialTheta
           << "\n     -- min theta  : " << M_minTheta
           << "\n     -- maxit : " << M_fixPointMaxIt;
    *_ostr << "\n||==============================================||"
           << "\n";

    return _ostr;
}


} // namespace FeelModels
} // namespace Feel
