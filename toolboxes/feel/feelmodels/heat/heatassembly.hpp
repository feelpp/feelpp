#ifndef FEELPP_TOOLBOXES_HEAT_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_HEAT_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisTemperatureType >
template <typename SymbolsExpr>
void
Heat<ConvexType,BasisTemperatureType>::updateLinearPDE( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

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

    std::string symbolStr = "heat_T";
    //--------------------------------------------------------------------------------------------------//

    for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isMatrix() )
        {
            auto const& kappa = expr( thermalConductivity.template exprMatrix<nDim,nDim>(), symbolsExpr );
            bool buildDiffusion = kappa.expression().isConstant()? buildCstPart : buildNonCstPart;
            if ( buildDiffusion )
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
            if ( buildDiffusion )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*kappa*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }


        if ( this->fieldVelocityConvectionIsUsedAndOperational() || !this->isStationary() )
        {
            auto const& rhoHeatCapacity = this->thermalProperties()->rhoHeatCapacity( matName );
            auto const& rhoHeatCapacityExpr = expr( rhoHeatCapacity.expr(), symbolsExpr );
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoHeatCapacityExpr*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
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
                if ( BuildNonCstPart_Form1TransientTerm )
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= rhoHeatCapacityExpr*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }

            // update stabilization gls
            if ( M_stabilizationGLS && buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                CHECK( !thermalConductivity.isMatrix() ) << "NotImplemented";
                auto const& kappa = expr( thermalConductivity.expr(), symbolsExpr );
                this->updateLinearPDEStabilizationGLS(  rhoHeatCapacityExpr, kappa, idv(this->fieldVelocityConvection()), range, data );
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
        if ( buildSourceTerm )
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
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }

        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expression1( d,symbolsExpr );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
            auto theExpr2 = expression2( d,symbolsExpr );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Heat","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );

}



template< typename ConvexType, typename BasisTemperatureType >
template <typename SymbolsExpr>
void
Heat<ConvexType,BasisTemperatureType>::updateJacobian( DataUpdateJacobian & data, SymbolsExpr const& symbolsExpr ) const
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
    //auto const& u = this->fieldTemperature();
    auto const u = Xh->element(XVec, this->rowStartInVector());
    auto const& v = this->fieldTemperature();


    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isConstant() )
        {
            if ( buildCstPart )
            {
                double kappa = thermalConductivity.value();
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*kappa*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto kappaExpr = expr( thermalConductivity.expr(), symbolsExpr );
            std::string symbolStr = "heat_T";
            if ( kappaExpr.expression().hasSymbol( symbolStr ) )
            {
                if ( buildNonCstPart )
                {
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*kappaExpr*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                    auto kappaDiffExpr = diff( kappaExpr,symbolStr,1,"",this->worldComm(),this->repository().expr());
                    //auto kappaDiffEval = expr( kappaDiff.expression(), symbolsExpr );
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*kappaDiffExpr*idt(u)*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                if ( buildCstPart )
                {
                    //auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*kappaExpr*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }


        if ( this->fieldVelocityConvectionIsUsedAndOperational() || !this->isStationary() )
        {
            auto const& rhoHeatCapacity = this->thermalProperties()->rhoHeatCapacity( matName );
            auto rhoHeatCapacityExpr = expr(rhoHeatCapacity.expr(),symbolsExpr);

            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoHeatCapacityExpr*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
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
            if ( M_stabilizationGLS && buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                CHECK( !thermalConductivity.isMatrix() ) << "NotImplemented";
                auto const& kappa = expr(thermalConductivity.expr(),symbolsExpr);
                this->updateJacobianStabilizationGLS(  rhoHeatCapacityExpr, kappa, idv(this->fieldVelocityConvection()), range, data );
            }
        }

    }

    //--------------------------------------------------------------------------------------------------//
    // update weak bc
    if ( buildNonCstPart )
    {
        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expression1( d,symbolsExpr );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->M_bcRadiation )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRadiationBC( name(d) ) ),
                           _expr= timeSteppingScaling*expression1(d)*4*pow(idv(v),3)*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
    }

}

template< typename ConvexType, typename BasisTemperatureType >
template <typename SymbolsExpr>
void
Heat<ConvexType,BasisTemperatureType>::updateResidual( DataUpdateResidual & data, SymbolsExpr const& symbolsExpr ) const
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
    auto const u = Xh->element(XVec, this->rowStartInVector());

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() );

    for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isConstant() )
        {
            if ( buildNonCstPart && !UseJacobianLinearTerms )
            {
                double kappa = thermalConductivity.value();
                myLinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*kappa*inner(gradv(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto kappaExpr = expr(thermalConductivity.expr(),symbolsExpr);
            std::string symbolStr = "heat_T";
            if ( kappaExpr.expression().hasSymbol( symbolStr ) )
            {
                if ( buildNonCstPart )
                {
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*kappaExpr*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                if ( buildNonCstPart && !UseJacobianLinearTerms )
                {
                    //auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*kappaExpr*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }


        if ( this->fieldVelocityConvectionIsUsedAndOperational() || !this->isStationary() )
        {
            auto const& rhoHeatCapacity = this->thermalProperties()->rhoHeatCapacity( matName );
            auto rhoHeatCapacityExpr = expr(rhoHeatCapacity.expr(),symbolsExpr);
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                myLinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoHeatCapacityExpr*(gradv(u)*idv(this->fieldVelocityConvection()))*id(v),
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
            if ( M_stabilizationGLS && buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                CHECK( !thermalConductivity.isMatrix() ) << "NotImplemented";
                auto const& kappa = expr(thermalConductivity.expr(),symbolsExpr);
                this->updateResidualStabilizationGLS( rhoHeatCapacityExpr, kappa, idv(this->fieldVelocityConvection()), range, data );
            }
        }

    }

    //--------------------------------------------------------------------------------------------------//
    // update source term
    if ( buildCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto theExpr = expression(d,symbolsExpr);
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
        if ( buildCstPart )
        {
            auto theExpr = expression(d,symbolsExpr);
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

    for( auto const& d : this->M_bcRobin )
    {
        auto theExpr1 = expression1( d,symbolsExpr );
        if ( buildNonCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= timeSteppingScaling*theExpr1*idv(u)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            auto theExpr2 = expression2( d,symbolsExpr );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= -timeSteppingScaling*theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }
    }
    for( auto const& d : this->M_bcRadiation )
    {
        if ( !buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRadiationBC( name(d) ) ),
                           _expr= timeSteppingScaling*expression1(d)*pow(idv(u),4)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRadiationBC( name(d) ) ),
                           _expr= -timeSteppingScaling*expression1(d)*pow(expression2(d),4)*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->log("Heat","updateResidual", "finish");

}



} // namespace Feel
} // namespace FeelModels

#endif
