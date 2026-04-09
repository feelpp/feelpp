#ifndef FEELPP_TOOLBOXES_MAGNETIC_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_MAGNETIC_ASSEMBLY_HPP

#include <feel/feelmodels/modelcore/diffsymbolicexpr.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisMagneticVectorPotentialType>
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );
    bool doAssemblyLhs = !data.hasInfo( "ignore-assembly.lhs" );


    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Magnetic","updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();


#if 0
    bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
    bool BuildNonCstPart_Form1TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        BuildNonCstPart_Form2TransientTerm = buildCstPart;
#endif
    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
#if 0
      if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
#endif
    }

    auto const& symbolsExpr = mctx.symbolsExpr();

    auto mesh = this->mesh();
    auto Xh = this->spaceVectorPotential();
    auto const& u = this->fieldVectorPotential();
    auto const& v = this->fieldVectorPotential();

    size_type startBlockIndexVectorPotential = this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() );

    auto bilinearForm_A_A = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix()+startBlockIndexVectorPotential,
                                   _colstart=this->colStartInMatrix()+startBlockIndexVectorPotential );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector()+startBlockIndexVectorPotential );
    //--------------------------------------------------------------------------------------------------//

    double mu_0 = 1.25663706127e-6;

    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicMagneticData = std::static_pointer_cast<ModelPhysicMagnetic<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto const& magneticRelativePermeability = this->materialsProperties()->materialProperty( matName, "magnetic-relative-permeability" );
            if ( magneticRelativePermeability.isMatrix() )
              {
                if constexpr (  nDim == 3 )
                  {
                    auto const& mu_r = expr( magneticRelativePermeability.template expr<nDim,nDim>(), symbolsExpr );
                    bool buildRotRot = mu_r.expression().isConstant()? buildCstPart : buildNonCstPart;
                    if ( doAssemblyLhs && buildRotRot )
                      {
                        bilinearForm_A_A +=
                          integrate( _range=range,
                                     _expr= timeSteppingScaling*(1./mu_0)*inner(inv(mu_r)*curlt(u),curl(v)),
                                     _geomap=this->geomap() );
                      }
                  }
                else
                  CHECK( false ) << "material property magnetic-relative-permeability should be a scalar in 2D";
            }
            else
            {
              auto const& mu_r = expr( magneticRelativePermeability.expr(), symbolsExpr );
              bool buildRotRot = mu_r.expression().isConstant()? buildCstPart : buildNonCstPart;
              if ( doAssemblyLhs && buildRotRot )
                {
                    bilinearForm_A_A +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*(1./(mu_0*mu_r))*inner(curlt(u),curl(v)),
                                   _geomap=this->geomap() );
                }
            }

            // // additional term in regularized formulation
            // if ( buildCstPart )
            //   {
            //     double epsilonPenal = 1.;
            //     bilinearForm_PatternCoupled +=
            //       integrate( _range=range,
            //                  _expr= timeSteppingScaling*epsilonPenal*inner(idt(u),id(v)),
            //                  _geomap=this->geomap() );

            //   }
            // current density sources
            for ( auto const& currentDensitySource : physicMagneticData->currentDensitySources() )
            {
                auto theExpr = currentDensitySource.expr( symbolsExpr );
                bool buildSourceTerm = theExpr.expression().isConstant()? buildCstPart : buildNonCstPart;
                if ( doAssemblyRhs && buildSourceTerm )
                  {
                    myLinearForm +=
                      integrate( _range=range,
                                 _expr= timeSteppingScaling*inner(theExpr,id(v)),
                                 _geomap=this->geomap() );
                  }
              }
        }
    }


    // additional term in regularized formulation
    if ( M_nullSpaceMethod == "regularized-formulation" && buildCstPart )
      {
        double epsilonPenal = 1.;
        bilinearForm_A_A +=
          integrate( _range=this->rangeMeshElements(),
                     _expr= timeSteppingScaling*epsilonPenal*inner(idt(u),id(v)),
                     _geomap=this->geomap() );
      }

    if ( M_nullSpaceMethod == "saddle-point" && buildCstPart )
      {
        auto XhLm = this->spaceLagrangeMultiplierCoulombGauge();
        auto const& p = this->fieldLagrangeMultiplierCoulombGauge();
        size_type startBlockIndexLmCoulombGauge = this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() );
        auto bilinearForm_lm_A = form2( _test=XhLm,_trial=Xh,_matrix=A,
                                        _pattern=size_type(Pattern::COUPLED),
                                        _rowstart=this->rowStartInMatrix()+startBlockIndexLmCoulombGauge,
                                        _colstart=this->colStartInMatrix()+startBlockIndexVectorPotential );
        auto bilinearForm_A_lm = form2( _test=Xh,_trial=XhLm,_matrix=A,
                                        _pattern=size_type(Pattern::COUPLED),
                                        _rowstart=this->rowStartInMatrix()+startBlockIndexVectorPotential,
                                        _colstart=this->colStartInMatrix()+startBlockIndexLmCoulombGauge );
        bilinearForm_A_lm +=
          integrate( _range=this->rangeMeshElements(),
                     _expr= timeSteppingScaling*inner(id(v),trans(gradt(p))),
                     _geomap=this->geomap() );

        bilinearForm_lm_A +=
          integrate( _range=this->rangeMeshElements(),
                     _expr= timeSteppingScaling*inner(idt(u),trans(grad(p))),
                     _geomap=this->geomap() );
      }


    //--------------------------------------------------------------------------------------------------//
    // update weak bc
#if 0
    if ( buildNonCstPart )
    {
        for ( auto const& [bcName,bcData] : M_boundaryConditions->heatFlux() )
        {
            auto theExpr = bcData->expr( symbolsExpr );
            if ( doAssemblyRhs )
            {
                myLinearForm +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= timeSteppingScaling*theExpr*id(v),
                               _geomap=this->geomap() );
            }
        }
        for ( auto const& [bcName,bcData] : M_boundaryConditions->convectiveHeatFlux() )
        {
            auto theExpr_h = bcData->expr_h( symbolsExpr );
            auto theExpr_Text = bcData->expr_Text( symbolsExpr );
            if ( doAssemblyLhs )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=markedfaces(mesh,bcData->markers()),
                               _expr= timeSteppingScaling*theExpr_h*idt(v)*id(v),
                               _geomap=this->geomap() );
            }
            if ( doAssemblyRhs )
            {
                myLinearForm +=
                    integrate( _range=markedfaces(mesh,bcData->markers()),
                               _expr= timeSteppingScaling*theExpr_h*theExpr_Text*id(v),
                               _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
#endif

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Magnetic","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );
}

template< typename ConvexType, typename BasisMagneticVectorPotentialType >
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;
    this->log("Magnetic","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceVectorPotential();
    auto const& u = this->fieldVectorPotential();
    size_type startBlockIndexVectorPotential = this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() );
    auto bilinearForm_A_A = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix()+startBlockIndexVectorPotential,
                                   _colstart=this->colStartInMatrix()+startBlockIndexVectorPotential );
    M_boundaryConditions->applyDofEliminationLinear( bilinearForm_A_A, F, mesh, u, se );

    if ( M_nullSpaceMethod == "saddle-point" )
      {
        auto XhLm = this->spaceLagrangeMultiplierCoulombGauge();
        auto const& p = this->fieldLagrangeMultiplierCoulombGauge();
        size_type startBlockIndexLmCoulombGauge = this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() );

        form2( _test=XhLm,_trial=XhLm,_matrix=A,
               _rowstart=this->rowStartInMatrix()+startBlockIndexLmCoulombGauge,
               _colstart=this->colStartInMatrix()+startBlockIndexLmCoulombGauge ) +=
          on( _range=boundaryfaces( support( this->spaceLagrangeMultiplierCoulombGauge() ) ),
              _element=p,_rhs=F,_expr=cst(0.),
              _vm=this->clovm(),_prefix=this->prefix() );
      }

    this->log("Magnetic","updateLinearPDEDofElimination","finish" );
}


template< typename ConvexType, typename BasisMagneticVectorPotentialType >
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const
{
#if 0
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;
    this->log("Magnetic","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexTemperature = this->startSubBlockSpaceIndex( "temperature" );
    auto u = this->spaceTemperature()->element( U, this->rowStartInVector()+startBlockIndexTemperature );
    auto const& se = mctx.symbolsExpr();

    M_boundaryConditions->applyNewtonInitialGuess( mesh, u, se );

    // update info for synchronization
    this->updateDofEliminationIds( "temperature", data );

    this->log("Magnetic","updateNewtonInitialGuess","finish" );
#endif
}

template< typename ConvexType, typename BasisMagneticVectorPotentialType >
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const
{
#if 0
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool _BuildCstPart = data.buildCstPart();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Magnetic","updateJacobian", "start"+sc);

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

    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicMagneticData = std::static_pointer_cast<ModelPhysicHeat<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
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
                            auto trialXh = trialSpacePair.second;
                            auto trialBlockIndex = trialSpacePair.first;

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

            if ( physicHeatData->hasConvectionEnabled() || !this->isStationary() )
            {
                auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
                auto rhoHeatCapacityExpr = expr(rhoHeatCapacity.expr(),se);

                if ( buildNonCstPart && physicHeatData->hasConvectionEnabled() )
                {
                    auto velConvExpr = physicHeatData->convection().expr( se );
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
                if ( M_stabilizationGLS && buildNonCstPart && physicHeatData->hasConvectionEnabled() )
                {
                    this->updateJacobianStabilizationGLS( data, mctx, *physicHeatData, matProps, range );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // update weak bc
    if ( buildNonCstPart )
    {
        for ( auto const& [bcName,bcData] : M_boundaryConditions->convectiveHeatFlux() )
        {
            auto theExpr_h = bcData->expr_h( se );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,bcData->markers()),
                           _expr= timeSteppingScaling*theExpr_h*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
        for ( auto const& bcDataPair : M_boundaryConditions->heatFlux() )
        {
            auto const& bcData = bcDataPair.second;
            auto neumannExprBase = bcData->expr();
            bool neumannnBcDependOnUnknown = neumannExprBase.hasSymbolDependency( trialSymbolNames, se );
            if ( neumannnBcDependOnUnknown )
            {
                auto neumannExpr = bcData->expr( se );
                //auto neumannExpr = expr( neumannExprBase, se );
                hana::for_each( tse.map(), [this,&bcData,&neumannExpr,&u,&v,&J,&Xh,&timeSteppingScaling]( auto const& e )
                {
                    for ( auto const& trialSpacePair : hana::second(e).blockSpaceIndex() )
                    {
                        auto trialXh = trialSpacePair.second;
                        auto trialBlockIndex = trialSpacePair.first;

                        auto neumannDiffExpr = diffSymbolicExpr( neumannExpr, hana::second(e), trialXh, trialBlockIndex, this->worldComm(), this->repository().expr() );

                        if ( !neumannDiffExpr.expression().hasExpr() )
                            continue;

                        form2( _test=Xh,_trial=trialXh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=trialBlockIndex ) +=
                            integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                       _expr= -timeSteppingScaling*inner(neumannDiffExpr, id(v)),
                                       _geomap=this->geomap() );
                    }
                });
            }
        }
    }
#endif
}

template< typename ConvexType, typename BasisMagneticVectorPotentialType >
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const
{
#if 0
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Magnetic","updateResidual", "start"+sc);

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

    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicHeatData = std::static_pointer_cast<ModelPhysicHeat<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
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

            for ( auto const& heatSource : physicHeatData->heatSources() )
            {
                auto theExpr = heatSource.expr( se );
                bool buildSourceTerm = buildCstPart;
                if ( buildSourceTerm )
                {
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= -timeSteppingScaling*theExpr*id(v),
                                   _geomap=this->geomap() );
                }
            }

            if ( physicHeatData->hasConvectionEnabled() || !this->isStationary() )
            {
                auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
                auto rhoHeatCapacityExpr = expr(rhoHeatCapacity.expr(),se);
                if ( buildNonCstPart && physicHeatData->hasConvectionEnabled() )
                {
                    auto velConvExpr = physicHeatData->convection().expr( se );
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
                if ( M_stabilizationGLS && buildNonCstPart && physicHeatData->hasConvectionEnabled() )
                {
                    this->updateResidualStabilizationGLS( data, mctx, *physicHeatData, matProps, range );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//
    // update weak bc
    for ( auto const& [bcName,bcData] : M_boundaryConditions->heatFlux() )
    {
        //auto theExpr = bcData.expr( se );
        auto neumannExprBase = bcData->expr();
        bool neumannnBcDependOnUnknown = neumannExprBase.hasSymbolDependency( trialSymbolNames, se );
        bool assembleNeumannBcTerm = neumannnBcDependOnUnknown? buildNonCstPart : buildCstPart;
        if ( assembleNeumannBcTerm )
        {
            //auto theExpr = expr( neumannExprBase, se );
            auto theExpr = bcData->expr( se );
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                           _expr= -timeSteppingScaling*theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }
    for ( auto const& [bcName,bcData] : M_boundaryConditions->convectiveHeatFlux() )
    {
        auto theExpr_h = bcData->expr_h( se );
        if ( buildNonCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,bcData->markers()),
                           _expr= timeSteppingScaling*theExpr_h*idv(u)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            auto theExpr_Text = bcData->expr_Text( se );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,bcData->markers()),
                           _expr= -timeSteppingScaling*theExpr_h*theExpr_Text*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    this->log("Magnetic","updateResidual", "finish");
#endif
}



} // namespace Feel
} // namespace FeelModels

#endif
