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
    using namespace Feel::vf;
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

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    size_type rowStartInVector = this->rowStartInVector();
    auto linearForm_PatternDefault = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::DEFAULT),
                                            _rowstart=rowStartInVector );
    auto linearForm_PatternCoupled = form1( _test=Xh, _vector=R,
                                            _pattern=size_type(Pattern::COUPLED),
                                            _rowstart=rowStartInVector );

    auto U = Xh->element(XVec, rowStartInVector);
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    //--------------------------------------------------------------------------------------------------//

    // strain tensor (trial)
    auto defv = sym(gradv(u));
    // identity Matrix
    auto Id = eye<nDim,nDim>();
    // dynamic viscosity
    auto const& mu = this->materialProperties()->fieldMu();
    auto const& rho = this->materialProperties()->fieldRho();
    // stress tensor (eval)
    auto Sigmav = -idv(p)*Id + 2*idv(mu)*defv;

    double timeElapsedBis=thetimerBis.elapsed();
    this->log("FluidMechanics","updateResidual","init done in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//

    thetimerBis.restart();
    if ( BuildNonCstPart )
    {
        if ( this->doStabConvectionEnergy() )
        {
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           //_expr= /*idv(*M_P0Rho)**/inner( Feel::vf::FSI::fluidMecConvection(u,*M_P0Rho) + idv(*M_P0Rho)*0.5*divv(u)*idv(u), id(v) ),
                           _expr=inner( Feel::FeelModels::fluidMecConvectionWithEnergyStab(u,rho), id(v) ),
                           _geomap=this->geomap() );

            /*if (this->isMoveDomain()  && !BuildCstPart && !UseJacobianLinearTerms)
             {
             linearForm_PatternCoupled +=
             integrate( _range=elements(mesh),
             _expr= -0.5*idv(M_P0Rho)*divv(this->meshVelocity())*trans(idv(u))*id(v),
             _geomap=this->geomap() );
             }*/
        }
        else
        {
            // convection term
#if 0
            auto convecTerm = val( idv(rho)*trans( gradv(u)*idv(u) ))*id(v);
#else
            auto convecTerm = inner( Feel::FeelModels::fluidMecConvection(u,rho),id(v) );
#endif
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr=convecTerm,
                           _geomap=this->geomap() );
        }

    }

    timeElapsedBis=thetimerBis.elapsed();thetimerBis.restart();
    this->log("FluidMechanics","updateResidual","build convective--1-- term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( M_isMoveDomain && !BuildCstPart && !UseJacobianLinearTerms )
    {
        // mesh velocity (convection) term
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -val(idv(rho)*trans( gradv(u)*( idv( this->meshVelocity() ))))*id(v),
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
        if ( dynamicViscosity.isNewtonianLaw() )
        {
            // sigma : grad(v) on Omega
            if ( BuildNonCstPart && !UseJacobianLinearTerms )
            {
                this->log("FluidMechanics","updateResidualModel","assembly with newtonian viscosity" );
                auto const mu_newtonian = idv(mu);
                auto const Sigmav_newtonian = -idv(p)*Id + 2*mu_newtonian*defv;
#if 1
                linearForm_PatternCoupled +=
                    integrate( _range=range,
                               //_expr= inner( StressTensorExpr,grad(v) ),
                               _expr= inner( val(Sigmav_newtonian),grad(v) ),
                               _geomap=this->geomap() );
#else
                form1( Xh, R ) +=
                    integrate( _range=range,
                               _expr= 2*idv(*M_P0Mu)*trace(trans(defv)*grad(v)),
                               _geomap=this->geomap() );
                form1( Xh, R ) +=
                    integrate( _range=range,
                               _expr= -idv(p)*div(v),
                               _geomap=this->geomap() );
#endif
            }
        }
        else
        {
            if ( BuildNonCstPart && !UseJacobianLinearTerms )
            {
                linearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= -idv(p)*div(v),
                               _geomap=this->geomap() );
            }
            if ( BuildNonCstPart )
            {
                auto const StressTensorExpr = Feel::FeelModels::fluidMecNewtonianStressTensor<2*nOrderVelocity>(u,p,*this->materialProperties(),matName,false/*true*/);
                // sigma : grad(v) on Omega
                linearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= inner( StressTensorExpr,grad(v) ),
                               _geomap=this->geomap() );
            }
        } // non newtonian
    }

    //--------------------------------------------------------------------------------------------------//
    // take into account that div u != 0
    if (!this->velocityDivIsEqualToZero() && BuildCstPart)
    {
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= idv(this->velocityDiv())*id(q),
                       _geomap=this->geomap() );

        auto coeffDiv = (2./3.)*idv(this->materialProperties()->fieldMu()); //(eps-2mu/3)
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= val(-coeffDiv*gradv(this->velocityDiv()))*id(v),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    // incompressibility term
    if (!BuildCstPart && !UseJacobianLinearTerms )
    {
        linearForm_PatternCoupled +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -divv(u)*id(q),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//

    // body forces
    if (BuildCstPart)
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
                linearForm_PatternCoupled +=
                    integrate( _range=rangeBodyForceUsed,
                               _expr= -inner( expression(d,this->symbolsExpr()),id(v) ),
                               _geomap=this->geomap() );
            }
        }

        if (M_haveSourceAdded)
        {
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -trans(idv(*M_SourceAdded))*id(v),
                           _geomap=this->geomap() );
        }
        if ( M_useGravityForce )
        {
            linearForm_PatternCoupled +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -idv(rho)*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );
        }
    }

    //------------------------------------------------------------------------------------//

    //transients terms
    if (!this->isStationaryModel())
    {
        bool Build_TransientTerm = !BuildCstPart;
        if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT ) Build_TransientTerm=!BuildCstPart && !UseJacobianLinearTerms;

        if (Build_TransientTerm) //  !BuildCstPart && !UseJacobianLinearTerms )
        {
            linearForm_PatternDefault +=
                integrate( _range=M_rangeMeshElements,
                           _expr= val(idv(rho)*trans(idv(u))*M_bdf_fluid->polyDerivCoefficient(0))*id(v),
                           _geomap=this->geomap() );
        }

        if (BuildCstPart)
        {
            auto Buzz = M_bdf_fluid->polyDeriv();
            auto buzz = Buzz.template element<0>();
            linearForm_PatternDefault +=
                integrate( _range=M_rangeMeshElements,
                           _expr= val(-idv(rho)*trans(idv(buzz)))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // user-defined additional terms
    this->updateResidualAdditional( R, BuildCstPart );

    //------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && !BuildCstPart && !UseJacobianLinearTerms )
        {
            double beta = this->definePressureCstPenalisationBeta();
            for ( auto const& rangeElt : M_definePressureCstMeshRanges )
                linearForm_PatternCoupled +=
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
            if ( BuildCstPart )
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

    this->updateResidualStabilisation( data, U );

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

    this->updateResidualWeakBC( data, U );

    //------------------------------------------------------------------------------------//

    double timeElapsed=thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateResidual",
                                               "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str()+
                                               "\n--------------------------------------------------",
                                               this->worldComm(),this->verboseAllProc());


} // updateResidual

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess(vector_ptrtype& U) const
{
    this->log("FluidMechanics","updateNewtonInitialGuess","start");

    size_type rowStartInVector = this->rowStartInVector();
    auto Xh = this->functionSpace();
    auto up = Xh->element( U, rowStartInVector );
    auto u = up.template element<0>();
    auto mesh = this->mesh();

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann")
    {
        this->log("FluidMechanics","updateNewtonInitialGuess","update moving boundary with strong Dirichlet");
        u.on(_range=markedfaces(mesh, this->markersNameMovingBoundary()),
             _expr=idv(this->meshVelocity2()) );
    }
#endif

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
    // synchronize velocity dof on interprocess
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("velocity");
    if ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )
        sync( u, "=", itFindDofsWithValueImposed->second );

    if ( this->definePressureCst() && this->definePressureCstMethod() == "algebraic" )
    {
        auto upSol = this->functionSpace()->element( U, this->rowStartInVector() );
        auto pSol = upSol.template element<1>();
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

    vector_ptrtype& R = data.residual();
    auto Xh = this->spaceVelocityPressure();
    size_type rowStartInVector = this->rowStartInVector();

    auto resFeView = Xh->element(R,rowStartInVector);
    auto resFeViewVelocity = resFeView.template element<0>();

    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("velocity");
    auto const& dofsWithValueImposedVelocity = ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposed->second : std::set<size_type>();
    for ( size_type thedof : dofsWithValueImposedVelocity )
        resFeViewVelocity.set( thedof,0. );
    sync( resFeViewVelocity, "=", dofsWithValueImposedVelocity );


    if ( this->hasMarkerPressureBC() )
    {
#if 0
        size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");
        auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM1 );
        lambdaPressure1.on(_range=boundaryfaces(M_meshLagrangeMultiplierPressureBC),
                           _expr=vf::zero<1,1>() );
        if ( nDim == 3 )
        {
            size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");
            auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM2 );
            lambdaPressure2.on(_range=boundaryfaces(M_meshLagrangeMultiplierPressureBC),
                               _expr=vf::zero<1,1>() );
        }
#else
        auto itFindDofsWithValueImposedPressureBC = M_dofsWithValueImposed.find("pressurebc-lm");
        auto const& dofsWithValueImposedPressureBC = ( itFindDofsWithValueImposedPressureBC != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposedPressureBC->second : std::set<size_type>();
        size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");
        auto lambdaPressure1 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM1 );
        for ( size_type thedof : dofsWithValueImposedPressureBC )
            lambdaPressure1.set( thedof,0. );
        sync( lambdaPressure1, "=", dofsWithValueImposedPressureBC );
        if ( nDim == 3 )
        {
            size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");
            auto lambdaPressure2 = M_spaceLagrangeMultiplierPressureBC->element( R/*XVec*/, rowStartInVector+startBlockIndexPressureLM2 );
            for ( size_type thedof : dofsWithValueImposedPressureBC )
                lambdaPressure2.set( thedof,0. );
            sync( lambdaPressure2, "=", dofsWithValueImposedPressureBC );
        }

#endif
    }

}




} // namespace FeelModels
} // namespace Feel


