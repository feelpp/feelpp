/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>


namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updatePicard( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool build_ConvectiveTerm = BuildNonCstPart;
    bool build_Form2TransientTerm = BuildNonCstPart;
    bool build_Form1TransientTerm = BuildNonCstPart;
    bool build_SourceTerm = BuildNonCstPart;
    bool build_StressTensorNonNewtonian = BuildNonCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        build_Form2TransientTerm=BuildCstPart;
    }
    if (this->useFSISemiImplicitScheme())
    {
        build_StressTensorNonNewtonian = BuildCstPart;
        if ( this->solverName() == "Oseen" )
            build_ConvectiveTerm=BuildCstPart;
        build_Form2TransientTerm=BuildCstPart;
        build_Form1TransientTerm=BuildCstPart;
        build_SourceTerm=BuildCstPart;
    }

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("FluidMechanics","updateLinearPDE", "start"+sc );
    this->timerTool("Solve").start();

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearFormVV_PatternDefault = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::DEFAULT),
                                              _rowstart=rowStartInMatrix+0,
                                              _colstart=colStartInMatrix+0 );
    auto bilinearFormVV_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix+0,
                                              _colstart=colStartInMatrix+0 );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix+0,
                                 _colstart=colStartInMatrix+1 );
    auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=rowStartInMatrix+1,
                                 _colstart=colStartInMatrix+0 );

    auto myLinearFormV =form1( _test=XhV, _vector=F,
                               _rowstart=rowStartInVector+0 );


    //std::shared_ptr<element_fluid_external_storage_type> fieldVelocityPressureExtrapolated;
    element_velocity_ptrtype fieldVelocityPressureExtrapolated;
    if ( this->solverName() == "Picard" )
    {
        fieldVelocityPressureExtrapolated = XhV->elementPtr();
        *fieldVelocityPressureExtrapolated = *XhV->elementPtr(*vecCurrentPicardSolution, rowStartInVector);
    }
    else if ( !this->isStationary() )
    {
        fieldVelocityPressureExtrapolated = M_fieldConvectionVelocityExtrapolated;
    }

    auto const& u = this->fieldVelocity();
    auto const& v = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    auto const& q = this->fieldPressure();

    // strain tensor (trial)
    auto deft = sym(gradt(u));
    //auto deft = 0.5*gradt(u);
    // density
    auto const& rho = this->materialProperties()->fieldRho();
    // identity matrix
    auto const Id = eye<nDim,nDim>();

    //--------------------------------------------------------------------------------------------------//
    this->timerTool("Solve").start();

    CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        // stress tensor sigma : grad(v)
        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& dynamicViscosity = this->materialProperties()->dynamicViscosity(matName);

            if ( ( dynamicViscosity.isConstant() && BuildCstPart ) ||
                 ( !dynamicViscosity.isConstant() && build_StressTensorNonNewtonian ) )
            {
                if ( fieldVelocityPressureExtrapolated )
                {
                    auto const& betaU = *fieldVelocityPressureExtrapolated;
                    auto myViscosity = Feel::FeelModels::fluidMecViscosity(gradv(betaU),*this->materialProperties(),matName);
                    bilinearFormVV_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*2*myViscosity*inner(deft,grad(v)),
                                   _geomap=this->geomap() );
                }
                else
                {
                    // case with steady Stokes
                    CHECK( dynamicViscosity.isNewtonianLaw() ) << "not allow with non newtonian law";

                    auto myViscosity = Feel::FeelModels::fluidMecViscosity( vf::zero<nDim,nDim>(),*this->materialProperties(),matName);
                    bilinearFormVV_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*2*myViscosity*inner(deft,grad(v)),
                                   _geomap=this->geomap() );
                }
            }
            if ( BuildCstPart )
            {
                bilinearFormVP +=
                    integrate( _range=range,
                               _expr= -div(v)*idt(p),
                               _geomap=this->geomap() );
            }

            //! gravity force
            if ( physicFluidData->gravityForceEnabled() )
            {
                auto const& gravityForce = physicFluidData->gravityForceExpr();
                bool assembleGravityTerm = gravityForce.expression().isNumericExpression()? BuildCstPart : BuildNonCstPart;
                if ( assembleGravityTerm )
                {
                    auto const& gravityForceExpr = gravityForce; // TODO apply symbols expr
                    myLinearFormV +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*idv(rho)*inner(gravityForceExpr,id(v)),
                                   _geomap=this->geomap() );
                }
            }

        } // foreach material

        // incompressibility term
        if ( BuildCstPart )
        {
            bilinearFormPV +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -divt(u)*id(q),
                           _geomap=this->geomap() );
        }

        double timeElapsedStressTensor = this->timerTool("Solve").stop();
        this->log("FluidMechanics","updateLinearPDE","assembly stress tensor + incompressibility in "+(boost::format("%1% s") %timeElapsedStressTensor).str() );

        //--------------------------------------------------------------------------------------------------//
        // convection
        if ( physicFluidData->equation() == "Navier-Stokes" && build_ConvectiveTerm )
        {
            this->timerTool("Solve").start();

            CHECK( this->solverName() == "Oseen" || this->solverName() == "Picard" ) << "invalid solver name " << this->solverName();
            auto const& betaU = *fieldVelocityPressureExtrapolated;
#if 0
            //velocityExprFromFields
            double myvelX=0;
            for ( auto const& [bpname,bpbc] : M_bodySetBC )
            {
                myvelX = bpbc.fieldTranslationalVelocityPtr()->operator()( 0 );
                break;
            }
            auto myVelXEXPR = vec( cst(myvelX), cst(0.) );
#endif
            if ( this->isMoveDomain() )
            {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                bilinearFormVV_PatternDefault +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= timeSteppingScaling*idv(rho)*trans( gradt(u)*( idv(betaU) -idv( this->meshVelocity() )   /*-  myVelXEXPR*/  ))*id(v),
                               _geomap=this->geomap() );
#endif
            }
            else
            {
                bilinearFormVV_PatternDefault +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= timeSteppingScaling*idv(rho)*trans( gradt(u)*idv(betaU) )*id(v),
                               _geomap=this->geomap() );
            }

            if ( this->doStabConvectionEnergy() )
            {
                bilinearFormVV_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= timeSteppingScaling*0.5*idv(rho)*divt(u)*trans(idv(betaU))*id(v),
                               _geomap=this->geomap() );
            }

            double timeElapsedConvection = this->timerTool("Solve").stop();
            this->log("FluidMechanics","updateLinearPDE","assembly convection in "+(boost::format("%1% s") %timeElapsedConvection).str() );
        }
        else if ( (  physicFluidData->equation() == "Stokes" ||  physicFluidData->equation() == "StokesTransient")
                  && build_ConvectiveTerm && this->isMoveDomain() )
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            bilinearFormVV_PatternDefault +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -timeSteppingScaling*idv(rho)*trans( gradt(u)*(idv( this->meshVelocity() )))*id(v),
                           _geomap=this->geomap() );
#endif
        }


        //--------------------------------------------------------------------------------------------------//
        //transients terms
        if ( !this->isStationary() && physicFluidData->equation() != "Stokes" )  //!this->isStationaryModel())
        {
            if (build_Form2TransientTerm)
            {
                bilinearFormVV_PatternDefault +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= idv(rho)*trans(idt(u))*id(v)*M_bdfVelocity->polyDerivCoefficient(0),
                               _geomap=this->geomap() );
            }

            if (build_Form1TransientTerm)
            {
                auto buzz = M_bdfVelocity->polyDeriv();
                myLinearFormV +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= idv(rho)*trans(idv(buzz))*id(v),
                               _geomap=this->geomap() );
            }
        }
    } //  for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    //--------------------------------------------------------------------------------------------------//
    // body forces
    if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        this->M_overwritemethod_updateSourceTermLinearPDE(F,BuildCstPart);
    }
    else
    {
        if ( build_SourceTerm )
        {
            for( auto const& d : this->M_volumicForcesProperties )
            {
                auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
                myLinearFormV +=
                    integrate( _range=rangeBodyForceUsed,
                               _expr= timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                               _geomap=this->geomap() );
            }
        }
    }

    // source given by user
    if ( M_haveSourceAdded && BuildNonCstPart)
    {
        myLinearFormV +=
            integrate( _range=M_rangeMeshElements,
                       _expr= timeSteppingScaling*trans(idv(*M_SourceAdded))*id(v),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && BuildCstPart  )
        {
            double beta = this->definePressureCstPenalisationBeta();
            auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=A,
                                         _pattern=size_type(Pattern::COUPLED),
                                         _rowstart=rowStartInMatrix+1,
                                         _colstart=colStartInMatrix+1 );
            for ( auto const& rangeElt : M_definePressureCstMeshRanges )
                bilinearFormPP +=
                    integrate( _range=rangeElt,
                               _expr=beta*idt(p)*id(q),
                               _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" )
        {
            CHECK( this->hasStartSubBlockSpaceIndex("define-pressure-cst-lm") ) << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startBlockIndexDefinePressureCstLM = this->startSubBlockSpaceIndex("define-pressure-cst-lm");

            if (BuildCstPart)
            {
                for ( int k=0;k<M_XhMeanPressureLM.size();++k )
                {
                    auto lambda = M_XhMeanPressureLM[k]->element();
                    form2( _test=XhP, _trial=M_XhMeanPressureLM[k], _matrix=A,
                           _rowstart=this->rowStartInMatrix()+1,
                           _colstart=this->colStartInMatrix()+startBlockIndexDefinePressureCstLM+k ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= id(p)*idt(lambda) /*+ idt(p)*id(lambda)*/,
                                   _geomap=this->geomap() );

                    form2( _test=M_XhMeanPressureLM[k], _trial=XhP, _matrix=A,
                           _rowstart=this->rowStartInMatrix()+startBlockIndexDefinePressureCstLM+k,
                           _colstart=this->colStartInMatrix()+1 ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= + idt(p)*id(lambda),
                                   _geomap=this->geomap() );
                }
            }

#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE)
            if (BuildNonCstPart)
            {
                this->log("FluidMechanics","updateLinearPDE", "also add nonzero MEANPRESSURE" );
                for ( int k=0;k<M_XhMeanPressureLM.size();++k )
                {
                    auto lambda = M_XhMeanPressureLM[k]->element();
                    form1( _test=M_XhMeanPressureLM[k], _vector=F,
                           _rowstart=this->rowStartInMatrix()+startBlockIndexDefinePressureCstLM+k ) +=
                        integrate( _range=M_definePressureCstMeshRanges[k],
                                   _expr= FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE(this->shared_from_this())*id(lambda),
                                   _geomap=this->geomap() );
                }
            }
#endif
        } // if ( this->definePressureCstMethod() == "lagrange-multiplier" )
    } // if ( this->definePressureCst() )

    //--------------------------------------------------------------------------------------------------//
    // div u != 0
    if (!this->velocityDivIsEqualToZero() && BuildNonCstPart)
    {
        auto myLinearFormP =form1( _test=XhP, _vector=F,
                                   _rowstart=rowStartInVector+1 );
        myLinearFormP +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -idv(this->velocityDiv())*id(q),
                       _geomap=this->geomap() );

        auto coeffDiv = (2./3.)*idv(this->materialProperties()->fieldMu()); //(eps-2mu/3)
        myLinearFormV +=
            integrate( _range=M_rangeMeshElements,
                       _expr= val(timeSteppingScaling*coeffDiv*gradv(this->velocityDiv()))*id(v),
                       _geomap=this->geomap() );

        if ( this->doStabConvectionEnergy() )
        {
            auto const& betaU = *fieldVelocityPressureExtrapolated;
            myLinearFormV +=
                integrate( _range=M_rangeMeshElements,
                           _expr= timeSteppingScaling*0.5*idv(rho)*idv(this->velocityDiv())*trans(idv(betaU))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->updateLinearPDEStabilisation( data );

    //--------------------------------------------------------------------------------------------------//

    // others bc
    this->updateLinearPDEWeakBC( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("FluidMechanics","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );

} // updateLinearPDE


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    if ( !this->hasStrongDirichletBC() ) return;
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateLinearPDEDofElimination", "start",
                                               this->worldComm(),this->verboseAllProc());
    this->timerTool("Solve").start();

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto XhV = this->functionSpaceVelocity();
    auto mesh = this->mesh();
    size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");
    auto bilinearFormVV = form2( _test=XhV,_trial=XhV,_matrix=A,
                                 _rowstart=this->rowStartInMatrix()+startBlockIndexVelocity,
                                 _colstart=this->colStartInMatrix()+startBlockIndexVelocity );
    auto const& u = this->fieldVelocity();

    // store markers for each entities in order to apply strong bc with priority (points erase edges erace faces)
    std::map<std::string, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > mapMarkerBCToEntitiesMeshMarker;
    for( auto const& d : this->M_bcDirichlet )
    {
        mapMarkerBCToEntitiesMeshMarker[name(d)] =
            detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",name(d) ) );
    }
    std::map<std::pair<std::string,ComponentType>, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > mapCompMarkerBCToEntitiesMeshMarker;
    for ( auto const& bcDirComp : this->M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            mapCompMarkerBCToEntitiesMeshMarker[std::make_pair(name(d),comp)] =
                detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",name(d),comp ) );
        }
    }

    // apply strong Dirichle bc on velocity field
    for( auto const& d : this->M_bcDirichlet )
    {
        auto itFindMarker = mapMarkerBCToEntitiesMeshMarker.find( name(d) );
        if ( itFindMarker == mapMarkerBCToEntitiesMeshMarker.end() )
            continue;
        auto const& listMarkerFaces = std::get<0>( itFindMarker->second );
        if ( !listMarkerFaces.empty() )
            bilinearFormVV +=
                on( _range=markedfaces( mesh, listMarkerFaces ),
                    _element=u, _rhs=F, _expr=expression(d,this->symbolsExpr()) );
        auto const& listMarkerEdges = std::get<1>( itFindMarker->second );
        if ( !listMarkerEdges.empty() )
            bilinearFormVV +=
                on( _range=markededges( mesh, listMarkerEdges ),
                    _element=u, _rhs=F, _expr=expression(d,this->symbolsExpr()) );
        auto const& listMarkerPoints = std::get<2>( itFindMarker->second );
        if ( !listMarkerPoints.empty() )
            bilinearFormVV +=
                on( _range=markedpoints( mesh, listMarkerPoints ),
                    _element=u, _rhs=F, _expr=expression(d,this->symbolsExpr()) );
    }
    // apply strong Dirichle bc on velocity component
    for ( auto const& bcDirComp : this->M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto itFindMarker = mapCompMarkerBCToEntitiesMeshMarker.find( std::make_pair(name(d),comp) );
            if ( itFindMarker == mapCompMarkerBCToEntitiesMeshMarker.end() )
                continue;
            auto const& listMarkerFaces = std::get<0>( itFindMarker->second );
            if ( !listMarkerFaces.empty() )
                bilinearFormVV +=
                    on( _range=markedfaces( mesh, listMarkerFaces ),
                        _element=u[comp],
                        _rhs=F, _expr=expression(d,this->symbolsExpr()) );
            auto const& listMarkerEdges = std::get<1>( itFindMarker->second );
            if ( !listMarkerEdges.empty() )
                bilinearFormVV +=
                    on( _range=markededges( this->mesh(), listMarkerEdges ),
                        _element=u[comp],
                        _rhs=F, _expr=expression(d,this->symbolsExpr()) );
            auto const& listMarkerPoints = std::get<2>( itFindMarker->second );
            if ( !listMarkerPoints.empty() )
                bilinearFormVV +=
                    on( _range=markedpoints( mesh, listMarkerPoints ),
                        _element=u[comp],
                        _rhs=F, _expr=expression(d,this->symbolsExpr()) );
        }
    }

    for ( auto const& inletbc : M_fluidInletDesc )
    {
        std::string const& marker = std::get<0>( inletbc );
        auto const& inletVel = std::get<0>( M_fluidInletVelocityInterpolated.find(marker)->second );

        bilinearFormVV +=
            on( _range=markedfaces(this->mesh(),marker),
                _element=u, _rhs=F,
                _expr=-idv(inletVel)*N() );
    }

    for( auto const& d : M_bcMovingBoundaryImposed )
    {
        auto listMarkerFaces = M_bcMarkersMovingBoundaryImposed.markerDirichletBCByNameId( "elimination",name(d) );
        bilinearFormVV +=
            on( _range=markedfaces(this->mesh(),listMarkerFaces),
                _element=u, _rhs=F,
                _expr=idv(M_meshALE->velocity()) );
    }

    if ( this->hasMarkerPressureBC() )
    {
        auto rangePressureBC = boundaryfaces(M_meshLagrangeMultiplierPressureBC);
        size_type startBlockIndexPressureLM1 = this->startSubBlockSpaceIndex("pressurelm1");
        form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=A,
               _rowstart=this->rowStartInMatrix()+startBlockIndexPressureLM1,
               _colstart=this->colStartInMatrix()+startBlockIndexPressureLM1 ) +=
            on( _range=rangePressureBC, _rhs=F,
                _element=*M_fieldLagrangeMultiplierPressureBC1, _expr=cst(0.));
        if constexpr ( nDim == 3 )
        {
            size_type startBlockIndexPressureLM2 = this->startSubBlockSpaceIndex("pressurelm2");
            form2( _test=M_spaceLagrangeMultiplierPressureBC,_trial=M_spaceLagrangeMultiplierPressureBC,_matrix=A,
                   _rowstart=this->rowStartInMatrix()+startBlockIndexPressureLM2,
                   _colstart=this->colStartInMatrix()+startBlockIndexPressureLM2 ) +=
                on( _range=rangePressureBC, _rhs=F,
                    _element=*M_fieldLagrangeMultiplierPressureBC2, _expr=cst(0.));
        }
    }

    for ( auto const& [bpname,bpbc] : M_bodySetBC )
    {
        if ( bpbc.hasTranslationalVelocityExpr() )
        {
            size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
            form2( _test= bpbc.spaceTranslationalVelocity(),_trial=bpbc.spaceTranslationalVelocity(),_matrix=A,
                   _rowstart=this->rowStartInMatrix()+startBlockIndexTranslationalVelocity,
                   _colstart=this->colStartInMatrix()+startBlockIndexTranslationalVelocity ) +=
                on( _range=elements(bpbc.mesh()), _rhs=F,
                    _element=*bpbc.fieldTranslationalVelocityPtr(), _expr=bpbc.translationalVelocityExpr() );
        }
        if ( bpbc.hasAngularVelocityExpr() )
        {
            size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
            form2( _test=bpbc.spaceAngularVelocity(),_trial=bpbc.spaceAngularVelocity(),_matrix=A,
                   _rowstart=this->rowStartInMatrix()+startBlockIndexAngularVelocity,
                   _colstart=this->colStartInMatrix()+startBlockIndexAngularVelocity ) +=
                on( _range=elements(bpbc.mesh()), _rhs=F,
                    _element=*bpbc.fieldAngularVelocityPtr(), _expr=bpbc.angularVelocityExpr() );
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("FluidMechanics","updateLinearPDEDofElimination","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}



} // end namespace FeelModels
} // end namespace Feel


