#ifndef FEELPP_TOOLBOXES_HEAT_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_HEAT_ASSEMBLY_HPP 1

#include <feel/feelmodels/modelcore/diffsymbolicexpr.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType>
void
Heat<ConvexType,BasisTemperatureType>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );
    bool doAssemblyLhs = !data.hasInfo( "ignore-assembly.lhs" );


    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();

    bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
    bool BuildNonCstPart_Form1TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        BuildNonCstPart_Form2TransientTerm = buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto const& symbolsExpr = mctx.symbolsExpr();

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto const& v = this->fieldTemperature();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector() );

    //--------------------------------------------------------------------------------------------------//

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
        auto const& matProps = this->materialsProperties()->materialProperties( matName );
        auto const& thermalConductivity = this->materialsProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isMatrix() )
        {
            auto const& kappa = expr( thermalConductivity.template expr<nDim,nDim>(), symbolsExpr );
            bool buildDiffusion = kappa.expression().isConstant()? buildCstPart : buildNonCstPart;
            if ( doAssemblyLhs && buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*grad(v)*(kappa*trans(gradt(u))),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto const& kappa = expr( thermalConductivity.expr(), symbolsExpr );
            bool buildDiffusion = kappa.expression().isConstant()? buildCstPart : buildNonCstPart;
            if ( doAssemblyLhs && buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*kappa*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }


        if ( this->hasVelocityConvectionExpr( matName ) || !this->isStationary() )
        {
            auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
            auto const& rhoHeatCapacityExpr = expr( rhoHeatCapacity.expr(), symbolsExpr );
            if ( doAssemblyLhs && buildNonCstPart && this->hasVelocityConvectionExpr( matName ) )
            {
                auto velConvExpr = expr( this->velocityConvectionExpr( matName ), symbolsExpr );
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoHeatCapacityExpr*(gradt(u)*velConvExpr)*id(v),
                               _geomap=this->geomap() );
            }

            if ( !this->isStationary() )
            {
                if ( doAssemblyLhs && BuildNonCstPart_Form2TransientTerm )
                {
                    auto thecoeff = rhoHeatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= thecoeff*idt(u)*id(v),
                                   _geomap=this->geomap() );
                }
                if ( doAssemblyRhs && BuildNonCstPart_Form1TransientTerm )
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= rhoHeatCapacityExpr*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }

            // update stabilization gls
            if ( M_stabilizationGLS && buildNonCstPart && this->hasVelocityConvectionExpr( matName ) )
            {
                this->updateLinearPDEStabilizationGLS( data, mctx, *physicData, matProps, range );
            }
        }
    }

#if 0
    if ( this->fieldVelocityConvectionIsUsedAndOperational() && !buildCstPart )
    {
        //viscous dissipation
        if ( false/*true*/ )
        {
            double mu = 1.;
            auto defv = sym(gradv( this->fieldVelocityConvection() ) );
            auto defv2 = inner(defv,defv);

            if ( !this->fieldVelocityConvectionIsIncompressible() )
            {
#if 0
                bilinearForm_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= thecoeff*(idt(u)*divv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
#endif
                myLinearForm +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= 2*mu*defv2*id(v),
                               _geomap=this->geomap() );
            }
            else
            {
                auto incomp2 = pow( divv( this->fieldVelocityConvection() ),2 );
                myLinearForm +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= 2*mu*(defv2-(1./3.)*incomp2)*id(v),
                               _geomap=this->geomap() );
            }
        }
    }
#endif

    //--------------------------------------------------------------------------------------------------//
    // update source term
    for( auto const& d : this->M_volumicForcesProperties )
    {
        auto theExpr = expression( d,symbolsExpr );
        bool buildSourceTerm = theExpr.expression().isConstant()? buildCstPart : buildNonCstPart;
        if ( doAssemblyRhs && buildSourceTerm )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // update weak bc
    if ( buildNonCstPart )
    {
        for( auto const& d : this->M_bcNeumann )
        {
            auto theExpr = expression( d,symbolsExpr );
            if ( doAssemblyRhs )
            {
                myLinearForm +=
                    integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                               _expr= timeSteppingScaling*theExpr*id(v),
                               _geomap=this->geomap() );
            }
        }

        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expression1( d,symbolsExpr );
            if ( doAssemblyLhs )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                               _geomap=this->geomap() );
            }
            if ( doAssemblyRhs )
            {
                auto theExpr2 = expression2( d,symbolsExpr );
                myLinearForm +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr1*theExpr2*id(v),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Heat","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );

}

template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType>
void
Heat<ConvexType,BasisTemperatureType>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("Heat","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcDirichlet )
    {
        auto theExpr = expression(d,se);
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=theExpr );
    }

    this->log("Heat","updateLinearPDEDofElimination","finish" );
}


template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType>
void
Heat<ConvexType,BasisTemperatureType>::updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Heat","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexTemperature = this->startSubBlockSpaceIndex( "temperature" );
    auto u = this->spaceTemperature()->element( U, this->rowStartInVector()+startBlockIndexTemperature );
    auto const& se = mctx.symbolsExpr();

    for( auto const& d : M_bcDirichlet )
    {
        auto theExpr = expression(d,se);
        u.on(_range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=theExpr );
    }

    // update info for synchronization
    this->updateDofEliminationIds( "temperature", data );

    this->log("Heat","updateNewtonInitialGuess","finish" );
}

template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType>
void
Heat<ConvexType,BasisTemperatureType>::updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool _BuildCstPart = data.buildCstPart();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateJacobian", "start"+sc);

    bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        BuildNonCstPart_Form2TransientTerm = buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = mctx.field( FieldTag::temperature(this), "temperature" );
    auto const& v = this->fieldTemperature();
    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
        auto const& matProps = this->materialsProperties()->materialProperties( matName );
        auto const& thermalConductivity = this->materialsProperties()->thermalConductivity( matName );

        bool thermalConductivityDependOnTrialSymbol = thermalConductivity.hasSymbolDependency( trialSymbolNames,se );
        if ( thermalConductivity.template hasExpr<nDim,nDim>() )
        {
            auto const& kappaExpr = expr( thermalConductivity.template expr<nDim,nDim>(), se );
            bool buildDiffusion = kappaExpr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*grad(v)*(kappaExpr*trans(gradt(u))),
                               _geomap=this->geomap() );
            }

            if ( thermalConductivityDependOnTrialSymbol && buildNonCstPart )
            {
                CHECK( false ) << "TODO";
            }
        }
        else
        {
            auto kappaExpr = expr( thermalConductivity.expr(), se );
            bool buildDiffusion = kappaExpr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*kappaExpr*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
            if ( thermalConductivityDependOnTrialSymbol && buildNonCstPart )
            {
                hana::for_each( tse.map(), [this,&kappaExpr,&u,&v,&J,&range,&Xh,&timeSteppingScaling]( auto const& e )
                    {
                        // NOTE : a strange compilation error related to boost fusion if we use [trialXh,trialBlockIndex] in the loop for
                        for ( auto const& trialSpacePair /*[trialXh,trialBlockIndex]*/ : hana::second(e).blockSpaceIndex() )
                        {
                            auto trialXh = trialSpacePair.first;
                            auto trialBlockIndex = trialSpacePair.second;

                            auto kappaDiffExpr = diffSymbolicExpr( kappaExpr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                            if ( !kappaDiffExpr.expression().hasExpr() )
                                continue;

                            form2( _test=Xh,_trial=trialXh,_matrix=J,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=trialBlockIndex ) +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*kappaDiffExpr*inner(gradv(u),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    });
            }
        }


        if ( this->hasVelocityConvectionExpr( matName ) || !this->isStationary() )
        {
            auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
            auto rhoHeatCapacityExpr = expr(rhoHeatCapacity.expr(),se);

            if ( buildNonCstPart && this->hasVelocityConvectionExpr( matName ) )
            {
                auto velConvExpr = expr( this->velocityConvectionExpr( matName ), se );
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoHeatCapacityExpr*(gradt(u)*velConvExpr)*id(v),
                               _geomap=this->geomap() );
            }

            if ( !this->isStationary() )
            {
                if ( BuildNonCstPart_Form2TransientTerm )
                {
                    auto thecoeff = rhoHeatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= thecoeff*idt(u)*id(v),
                                   _geomap=this->geomap() );
                }
            }

            // update stabilization gls
            if ( M_stabilizationGLS && buildNonCstPart && this->hasVelocityConvectionExpr( matName ) )
            {
                this->updateJacobianStabilizationGLS( data, mctx, *physicData, matProps, range );
            }

        }

    }

    //--------------------------------------------------------------------------------------------------//
    // update weak bc
    if ( buildNonCstPart )
    {
        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expression1( d,se );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
        }

        for( auto const& d : M_bcNeumann )
        {
            auto neumannExprBase = expression( d );
            bool neumannnBcDependOnUnknown = neumannExprBase.hasSymbolDependency( trialSymbolNames, se );
            if ( neumannnBcDependOnUnknown )
            {
                auto neumannExpr = expr( neumannExprBase, se );
                hana::for_each( tse.map(), [this,&d,&neumannExpr,&u,&v,&J,&Xh,&timeSteppingScaling]( auto const& e )
                {
                    for ( auto const& trialSpacePair : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.first;
                        auto trialBlockIndex = trialSpacePair.second;

                        auto neumannDiffExpr = diffSymbolicExpr( neumannExpr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                        if ( !neumannDiffExpr.expression().hasExpr() )
                            continue;

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                                       _expr= -timeSteppingScaling*inner(neumannDiffExpr, id(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }

    }

}

template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType>
void
Heat<ConvexType,BasisTemperatureType>::updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateResidual", "start"+sc);

    bool Build_TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        Build_TransientTerm=buildNonCstPart && !UseJacobianLinearTerms;

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm( this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( M_timeStepping == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - M_timeStepThetaValue;
            else
                timeSteppingScaling = M_timeStepThetaValue;
        }
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }


    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& v = this->fieldTemperature();
    auto const& u = mctx.field( FieldTag::temperature(this), "temperature" );
    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
        auto const& matProps = this->materialsProperties()->materialProperties( matName );
        auto const& thermalConductivity = this->materialsProperties()->thermalConductivity( matName );

        if ( thermalConductivity.template hasExpr<nDim,nDim>() )
        {
            auto const& kappaExpr = expr( thermalConductivity.template expr<nDim,nDim>(), se );
            bool buildDiffusion = kappaExpr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( buildDiffusion )
            {
                myLinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*grad(v)*(kappaExpr*trans(gradv(u))),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto kappaExpr = expr(thermalConductivity.expr(),se);
            bool buildDiffusion = kappaExpr.expression().isNumericExpression()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
            if ( buildDiffusion )
            {
                myLinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*kappaExpr*inner(gradv(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }

        if ( this->hasVelocityConvectionExpr( matName ) || !this->isStationary() )
        {
            auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
            auto rhoHeatCapacityExpr = expr(rhoHeatCapacity.expr(),se);
            if ( buildNonCstPart && this->hasVelocityConvectionExpr( matName ) )
            {
                auto velConvExpr = expr( this->velocityConvectionExpr( matName ), se );
                myLinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoHeatCapacityExpr*(gradv(u)*velConvExpr)*id(v),
                               _geomap=this->geomap() );
            }
            if ( !this->isStationary() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
            {
                if ( Build_TransientTerm )
                {
                    auto thecoeff = rhoHeatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= thecoeff*idv(u)*id(v),
                                   _geomap=this->geomap() );
                }
                if (buildCstPart)
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= -rhoHeatCapacityExpr*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }

            // update stabilization gls
            if ( M_stabilizationGLS && buildNonCstPart && this->hasVelocityConvectionExpr( matName ) )
            {
                this->updateResidualStabilizationGLS( data, mctx, *physicData, matProps, range );
            }

        }

    }

    //--------------------------------------------------------------------------------------------------//
    // update source term
    if ( buildCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto theExpr = expression(d,se);
            auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= -timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // update weak bc
    for( auto const& d : this->M_bcNeumann )
    {
        auto neumannExprBase = expression( d );
        bool neumannnBcDependOnUnknown = neumannExprBase.hasSymbolDependency( trialSymbolNames, se );
        bool assembleNeumannBcTerm = neumannnBcDependOnUnknown? buildNonCstPart : buildCstPart;
        if ( assembleNeumannBcTerm )
        {
            auto theExpr = expr( neumannExprBase, se );
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

    for( auto const& d : this->M_bcRobin )
    {
        auto theExpr1 = expression1( d,se );
        if ( buildNonCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idv(u)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            auto theExpr2 = expression2( d,se );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                           _expr= -timeSteppingScaling*theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->log("Heat","updateResidual", "finish");

}



} // namespace Feel
} // namespace FeelModels

#endif
