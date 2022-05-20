/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */


#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelvf/fluidmecconvection.hpp>
#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !BuildCstPart;

    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    if (this->verbose()) Feel::FeelModels::Log("--------------------------------------------------\n",
                                               this->prefix()+".FluidMechanics","updateResidual","start"+sc,
                                               this->worldComm(),this->verboseAllProc());

    boost::mpi::timer thetimer, thetimerBis;


    //--------------------------------------------------------------------------------------------------//

    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );
    // if ( !doAssemblyRhs )
    //     std::cout << "hola \n";


    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( M_timeStepping == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - M_timeStepThetaValue;
            else
                timeSteppingScaling = M_timeStepThetaValue;
        }
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");
    size_type startBlockIndexPressure = this->startSubBlockSpaceIndex("pressure");
    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormV_PatternDefault = form1( _test=XhV, _vector=R,
                                            _pattern=size_type(Pattern::DEFAULT),
                                            _rowstart=rowStartInVector+startBlockIndexVelocity );
    auto linearFormV_PatternCoupled = form1( _test=XhV, _vector=R,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=rowStartInVector+startBlockIndexVelocity );
    auto linearFormP = form1( _test=XhP, _vector=R,
                              _pattern=size_type(Pattern::COUPLED),
                              _rowstart=rowStartInVector+startBlockIndexPressure );

    auto u = XhV->element(XVec, rowStartInVector+startBlockIndexVelocity);
    auto const& v = u;
    auto p = XhP->element(XVec, rowStartInVector+startBlockIndexPressure);
    auto const& q = p;

    //--------------------------------------------------------------------------------------------------//

    // identity Matrix
    auto Id = eye<nDim,nDim>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();
    auto const& rho = this->materialProperties()->fieldRho();

    double timeElapsedBis=thetimerBis.elapsed();
    this->log("FluidMechanics","updateResidual","init done in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//

    thetimerBis.restart();
    CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);

        if ( BuildNonCstPart && physicFluidData->equation() == "Navier-Stokes" )
        {
            if ( this->solverName() == "Oseen" ) // call when evaluate residual for time-stepping
            {
                auto const& betaU = *M_fieldConvectionVelocityExtrapolated;
                linearFormV_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= timeSteppingScaling*idv(rho)*trans( gradv(u)*idv(betaU) )*id(v),
                               _geomap=this->geomap() );
                if ( this->doStabConvectionEnergy() )
                    CHECK( false ) << "TODO";
            }
            else
            {
                if ( this->doStabConvectionEnergy() )
                {
                    linearFormV_PatternCoupled +=
                        integrate( _range=M_rangeMeshElements,
                                   //_expr= /*idv(*M_P0Rho)**/inner( Feel::vf::FSI::fluidMecConvection(u,*M_P0Rho) + idv(*M_P0Rho)*0.5*divv(u)*idv(u), id(v) ),
                                   _expr=timeSteppingScaling*inner( Feel::FeelModels::fluidMecConvectionWithEnergyStab(u,rho), id(v) ),
                                   _geomap=this->geomap() );
                }
                else
                {
                    // convection term
                    // auto convecTerm = val( idv(rho)*trans( gradv(u)*idv(u) ))*id(v);
                    auto convecTerm = inner( Feel::FeelModels::fluidMecConvection(u,rho),id(v) );
                    linearFormV_PatternCoupled +=
                        integrate( _range=M_rangeMeshElements,
                                   _expr=timeSteppingScaling*convecTerm,
                                   _geomap=this->geomap() );
                }
            }
        }


        if ( !BuildCstPart && !UseJacobianLinearTerms && physicFluidData->equation() == "Navier-Stokes" && data.hasVectorInfo( "explicit-part-of-solution" ) )
        {
            auto uExplicitPartOfSolution = XhV->element( data.vectorInfo( "explicit-part-of-solution" ), rowStartInVector+startBlockIndexVelocity );
            linearFormV_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= timeSteppingScaling*val( idv(rho)*trans( gradv(u)*idv(uExplicitPartOfSolution) + gradv(uExplicitPartOfSolution )*idv(u)  ))*id(v),
                           _geomap=this->geomap() );
        }


        //! gravity force
        if ( doAssemblyRhs && physicFluidData->gravityForceEnabled() )
        {
            auto const& gravityForce = physicFluidData->gravityForceExpr();
            bool assembleGravityTerm = gravityForce.expression().isNumericExpression()? BuildCstPart : BuildNonCstPart;
            if ( assembleGravityTerm )
            {
                auto const& gravityForceExpr = gravityForce; // TODO apply symbols expr
                linearFormV_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= -timeSteppingScaling*idv(rho)*inner(gravityForceExpr,id(u)),
                               _geomap=this->geomap() );
            }
        }

    } // foreach physic

    timeElapsedBis=thetimerBis.elapsed();thetimerBis.restart();
    this->log("FluidMechanics","updateResidual","build convective--1-- term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( M_isMoveDomain && !BuildCstPart && !UseJacobianLinearTerms )
    {
        // mesh velocity (convection) term
        linearFormV_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -timeSteppingScaling*val(idv(rho)*trans( gradv(u)*( idv( this->meshVelocity() ))))*id(v),
                       _geomap=this->geomap() );
        timeElapsedBis=thetimerBis.elapsed();
        this->log("FluidMechanics","updateResidual","build convective--2-- term in "+(boost::format("%1% s") % timeElapsedBis ).str() );
    }
#endif

    //--------------------------------------------------------------------------------------------------//
    // sigma : grad(v) on Omega
    for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& dynamicViscosity = this->materialProperties()->dynamicViscosity(matName);

        if ( !timeSteppingEvaluateResidualWithoutTimeDerivative && BuildNonCstPart && !UseJacobianLinearTerms )
        {
            linearFormV_PatternCoupled +=
                integrate( _range=range,
                           _expr= -idv(p)*div(v),
                           _geomap=this->geomap() );
        }
        if ( ( dynamicViscosity.isConstant() && BuildNonCstPart && !UseJacobianLinearTerms ) ||
             ( !dynamicViscosity.isConstant() && BuildNonCstPart ) )
        {
            auto const StressTensorExpr = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(u),idv(p),*this->materialProperties(),matName,false/*true*/);
            // sigma : grad(v) on Omega
            linearFormV_PatternCoupled +=
                integrate( _range=range,
                           _expr= timeSteppingScaling*inner( StressTensorExpr,grad(v) ),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // take into account that div u != 0
    if (!this->velocityDivIsEqualToZero() && BuildCstPart && doAssemblyRhs )
    {
        linearFormP +=
            integrate( _range=M_rangeMeshElements,
                       _expr= idv(this->velocityDiv())*id(q),
                       _geomap=this->geomap() );

        auto coeffDiv = (2./3.)*idv(this->materialProperties()->fieldMu()); //(eps-2mu/3)
        linearFormV_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= val(-timeSteppingScaling*coeffDiv*gradv(this->velocityDiv()))*id(v),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    // incompressibility term
    if (!BuildCstPart && !UseJacobianLinearTerms )
    {
        linearFormP +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -divv(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    // body forces
    if (BuildCstPart && doAssemblyRhs)
    {
        if ( this->M_overwritemethod_updateSourceTermResidual != NULL )
        {
            this->M_overwritemethod_updateSourceTermResidual(R);
        }
        else
        {
            for( auto const& d : this->M_volumicForcesProperties )
            {
                auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
                linearFormV_PatternCoupled +=
                    integrate( _range=rangeBodyForceUsed,
                               _expr= -timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                               _geomap=this->geomap() );
            }
        }

        if (M_haveSourceAdded)
        {
            linearFormV_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -timeSteppingScaling*trans(idv(*M_SourceAdded))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //------------------------------------------------------------------------------------//

    //transients terms
    if ( !this->isStationaryModel() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
    {
        bool Build_TransientTerm = !BuildCstPart;
        if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=!BuildCstPart && !UseJacobianLinearTerms;

        if (Build_TransientTerm) //  !BuildCstPart && !UseJacobianLinearTerms )
        {
            linearFormV_PatternDefault +=
                integrate( _range=M_rangeMeshElements,
                           _expr= val(idv(rho)*trans(idv(u))*M_bdfVelocity->polyDerivCoefficient(0))*id(v),
                           _geomap=this->geomap() );
        }

        if (BuildCstPart && doAssemblyRhs)
        {
            auto buzz = M_bdfVelocity->polyDeriv();
            linearFormV_PatternDefault +=
                integrate( _range=M_rangeMeshElements,
                           _expr= val(-idv(rho)*trans(idv(buzz)))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && !BuildCstPart && !UseJacobianLinearTerms )
        {
            double beta = this->definePressureCstPenalisationBeta();
            for ( auto const& rangeElt : M_definePressureCstMeshRanges )
                linearFormP +=
                    integrate( _range=rangeElt,
                               _expr=beta*idv(p)*id(q),
                               _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("define-pressure-cst-lm") ) << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startBlockIndexDefinePressureCstLM = this->startSubBlockSpaceIndex("define-pressure-cst-lm");

            if ( !BuildCstPart && !UseJacobianLinearTerms )
            {
                for ( int k=0;k<M_XhMeanPressureLM.size();++k )
                {
                    auto lambda = M_XhMeanPressureLM[k]->element(XVec,rowStartInVector+startBlockIndexDefinePressureCstLM+k);
                    //M_blockVectorSolution.setSubVector( lambda, *XVec, rowStartInVector+startBlockIndexDefinePressureCstLM+k );
                    //for ( size_type k=0;k<M_XhMeanPressureLM->nLocalDofWithGhost();++k )
                    //    lambda( k ) = XVec->operator()( startDofIndexDefinePressureCstLM + k);

                    form1( _test=M_XhMeanPressureLM[k],_vector=R,
                           _rowstart=rowStartInVector+startBlockIndexDefinePressureCstLM+k ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= id(p)*idv(lambda) + idv(p)*id(lambda),
                                   _geomap=this->geomap() );
                }
            }
#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE)
            if ( BuildCstPart && doAssemblyRhs )
            {
                for ( int k=0;k<M_XhMeanPressureLM.size();++k )
                {
                    auto lambda = M_XhMeanPressureLM[k]->element();
                    form1( _test=M_XhMeanPressureLM[k],_vector=R,
                           _rowstart=rowStartInVector+startDofIndexDefinePressureCstLM+k ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= -(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE(this->shared_from_this()))*id(lambda),
                                   _geomap=this->geomap() );
                }
            }
#endif
        }
    }


    //------------------------------------------------------------------------------------//

    this->updateResidualStabilisation( data, u,p );

    //------------------------------------------------------------------------------------//
#if 0
    if ( UsePeriodicity && !BuildCstPart )
    {
        std::string marker1 = soption(_name="periodicity.marker1",_prefix=this->prefix());
        double pressureJump = doption(_name="periodicity.pressure-jump",_prefix=this->prefix());
        linearForm_PatternCoupled +=
            integrate( _range=markedfaces( this->mesh(),this->mesh()->markerName(marker1) ),
                       _expr=-inner(pressureJump*N(),id(v) ) );
    }
#endif
    //------------------------------------------------------------------------------------//

    this->updateResidualWeakBC( data, u,p );

    //------------------------------------------------------------------------------------//

    double timeElapsed=thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidual",
                                               "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str()+
                                               "\n--------------------------------------------------",
                                               this->worldComm(),this->verboseAllProc());


} // updateResidual

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    this->log("FluidMechanics","updateNewtonInitialGuess","start");

    vector_ptrtype& U = data.initialGuess();
    size_type rowStartInVector = this->rowStartInVector();
    auto XhV = this->functionSpaceVelocity();
    auto u = XhV->element( U, rowStartInVector );
    auto mesh = this->mesh();

    if ( this->hasMarkerDirichletBCelimination() )
    {
        // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
        std::map<std::string, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > mapMarkerBCToEntitiesMeshMarker;
        for( auto const& d : M_bcDirichlet )
        {
            mapMarkerBCToEntitiesMeshMarker[name(d)] =
                detail::distributeMarkerListOnSubEntity( mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) );
        }
        std::map<std::pair<std::string,ComponentType>, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > mapCompMarkerBCToEntitiesMeshMarker;
        for ( auto const& bcDirComp : M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                mapCompMarkerBCToEntitiesMeshMarker[std::make_pair(name(d),comp)] =
                    detail::distributeMarkerListOnSubEntity(mesh, this->markerDirichletBCByNameId( "elimination",name(d), comp )   );
            }
        }

        // strong Dirichlet bc with vectorial velocity
        for( auto const& d : M_bcDirichlet )
        {
            auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
            if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
                continue;
            auto const& listMarkerFaces = std::get<0>( itFindMarker->second );
            if ( !listMarkerFaces.empty() )
                u.on(_range=markedfaces(mesh,listMarkerFaces ),
                     _expr=expression(d,this->symbolsExpr()) );
            auto const& listMarkerEdges = std::get<1>( itFindMarker->second );
            if ( !listMarkerEdges.empty() )
                u.on(_range=markededges(mesh,listMarkerEdges ),
                     _expr=expression(d,this->symbolsExpr()) );
            auto const& listMarkerPoints = std::get<2>( itFindMarker->second );
            if ( !listMarkerPoints.empty() )
                u.on(_range=markedpoints(mesh,listMarkerPoints ),
                     _expr=expression(d,this->symbolsExpr()) );
        }
        // strong Dirichlet bc with velocity components
        for ( auto const& bcDirComp : M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                auto itFindMarker = mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(name(d),comp) );
                if ( itFindMarker == mapCompMarkerBCToEntitiesMeshMarker.end() )
                    continue;
                auto const& listMarkerFaces = std::get<0>( itFindMarker->second );
                if ( !listMarkerFaces.empty() )
                    u[comp].on(_range=markedfaces(mesh,listMarkerFaces ),
                               _expr=expression(d,this->symbolsExpr()) );
                auto const& listMarkerEdges = std::get<1>( itFindMarker->second );
                if ( !listMarkerEdges.empty() )
                    u[comp].on(_range=markededges(mesh,listMarkerEdges ),
                               _expr=expression(d,this->symbolsExpr()) );
                auto const& listMarkerPoints = std::get<2>( itFindMarker->second );
                if ( !listMarkerPoints.empty() )
                    u[comp].on(_range=markedpoints(mesh,listMarkerPoints ),
                               _expr=expression(d,this->symbolsExpr()) );
            }
        }
    }

    for ( auto const& inletbc : M_fluidInletDesc )
    {
        std::string const& marker = std::get<0>( inletbc );
        auto itFindMark = M_fluidInletVelocityInterpolated.find(marker);
        if ( itFindMark == M_fluidInletVelocityInterpolated.end() )
            continue;
        auto const& inletVel = std::get<0>( itFindMark->second );
        u.on(_range=markedfaces(mesh, marker),
             _expr=-idv(inletVel)*N() );
    }

    for( auto const& d : M_bcMovingBoundaryImposed )
    {
        auto listMarkerFaces = M_bcMarkersMovingBoundaryImposed.markerDirichletBCByNameId( "elimination",name(d) );
        u.on( _range=markedfaces(this->mesh(),listMarkerFaces),
              _expr=idv(M_meshALE->velocity()) );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "velocity", data );

    // dofs imposed with pressure Dirichlet bc
    if ( this->hasMarkerPressureBC() )
    {
        auto rangePressureBC = boundaryfaces(M_meshLagrangeMultiplierPressureBC);
        size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");
        auto plm1 = M_spaceLagrangeMultiplierPressureBC->element( U, rowStartInVector+startBlockIndexPressureLM1 );
        plm1.on( _range=rangePressureBC,_expr=cst(0.) );
        this->updateDofEliminationIds( "pressurelm1", this->dofEliminationIds( "pressurebc-lm" ), data );
        if ( nDim == 3 )
        {
            size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");
            auto plm2 = M_spaceLagrangeMultiplierPressureBC->element( U, rowStartInVector+startBlockIndexPressureLM2 );
            plm2.on( _range=rangePressureBC,_expr=cst(0.) );
            this->updateDofEliminationIds( "pressurelm2", this->dofEliminationIds( "pressurebc-lm" ), data );
        }
    }

    // body bc
    for ( auto const& [bpname,bpbc] : M_bodySetBC )
    {
        if ( bpbc.hasTranslationalVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".translational-velocity";
            size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex( spaceName );
            auto uBodyTranslationalVelocity = bpbc.spaceTranslationalVelocity()->element( U, rowStartInVector+startBlockIndexTranslationalVelocity );
            uBodyTranslationalVelocity.on(_range=elements(bpbc.mesh()), _expr=bpbc.translationalVelocityExpr() );
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.translational-velocity" ), data );
        }
        if ( bpbc.hasAngularVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".angular-velocity";
            size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex( spaceName );
            auto uBodyAngularVelocity = bpbc.spaceAngularVelocity()->element( U, rowStartInVector+startBlockIndexAngularVelocity );
            uBodyAngularVelocity.on(_range=elements(bpbc.mesh()), _expr=bpbc.angularVelocityExpr() );
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.angular-velocity" ), data );
        }
    }

    // imposed mean pressure (TODO use updateDofEliminationIds)
    if ( this->definePressureCst() && this->definePressureCstMethod() == "algebraic" )
    {
        auto pSol = this->functionSpacePressure()->element( U, this->rowStartInVector()+1 );
        CHECK( !M_definePressureCstAlgebraicOperatorMeanPressure.empty() ) << "mean pressure operator does not init";

        for ( int k=0;k<M_definePressureCstAlgebraicOperatorMeanPressure.size();++k )
        {
            double meanPressureImposed = 0;
            double meanPressureCurrent = inner_product( *M_definePressureCstAlgebraicOperatorMeanPressure[k].first, pSol );
            for ( size_type dofId : M_definePressureCstAlgebraicOperatorMeanPressure[k].second )
                pSol(dofId) += (meanPressureImposed - meanPressureCurrent);
        }
        sync( pSol, "=" );
    }

    this->log("FluidMechanics","updateNewtonInitialGuess","finish");
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !this->hasStrongDirichletBC() ) return;

    this->updateDofEliminationIds( "velocity", data );

    if ( this->hasMarkerPressureBC() )
    {
        this->updateDofEliminationIds( "pressurelm1", this->dofEliminationIds( "pressurebc-lm" ), data );
        if ( nDim == 3 )
            this->updateDofEliminationIds( "pressurelm2", this->dofEliminationIds( "pressurebc-lm" ), data );
    }

    for ( auto const& [bpname,bpbc] : M_bodySetBC )
    {
        if ( bpbc.hasTranslationalVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".translational-velocity";
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.translational-velocity" ), data );
        }
        if ( bpbc.hasAngularVelocityExpr() )
        {
            std::string spaceName = "body-bc."+bpbc.name()+".angular-velocity";
            this->updateDofEliminationIds( spaceName, this->dofEliminationIds( "body-bc.angular-velocity" ), data );
        }
    }
}




} // namespace FeelModels
} // namespace Feel


