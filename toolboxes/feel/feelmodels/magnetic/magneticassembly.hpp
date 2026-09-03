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

    double mu_0_cst = ModelPhysicMagnetic<nDim>::vacuumPermeabilityConstant();
    double equationScaling = M_equationScalingUseVacuumPermeability? mu_0_cst : 1.0;

    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicMagneticData = std::static_pointer_cast<ModelPhysicMagnetic<nDim>>(physicData);
        //auto mu_0 = physicMagneticData->vacuumPermeabilityExpr();
        //double mu_0_cst = physicMagneticData->vacuumPermeabilityConstant();
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto const& magneticRelativePermeability = this->materialsProperties()->materialProperty( matName, "magnetic-relative-permeability" );
            if ( magneticRelativePermeability.isMatrix() )
            {
                if constexpr (  nDim == 3 )
                {
                    auto mu_r = expr( magneticRelativePermeability.template expr<nDim,nDim>(), symbolsExpr );
                    bool buildRotRot = mu_r.expression().isConstant()? buildCstPart : buildNonCstPart;
                    if ( doAssemblyLhs && buildRotRot )
                    {
                        bilinearForm_A_A +=
                            integrate( _range=range,
                                       _expr= equationScaling*timeSteppingScaling*(1./mu_0_cst)*inner(inv(mu_r)*curlt(u),curl(v)),
                                       _geomap=this->geomap() );
                    }
                }
                else
                    CHECK( false ) << "material property magnetic-relative-permeability should be a scalar in 2D";
            }
            else
            {
                auto mu_r = expr( magneticRelativePermeability.expr(), symbolsExpr );
                bool buildRotRot = mu_r.expression().isConstant()? buildCstPart : buildNonCstPart;
                if ( doAssemblyLhs && buildRotRot )
                {
                    bilinearForm_A_A +=
                        integrate( _range=range,
                                   _expr= equationScaling*timeSteppingScaling*(1./(mu_0_cst*mu_r))*inner(curlt(u),curl(v)),
                                   _geomap=this->geomap() );
                }
            }

            // current density sources
            for ( auto const& currentDensitySource : physicMagneticData->currentDensitySources() )
            {
                auto theExpr = currentDensitySource.expr( symbolsExpr );
                bool buildSourceTerm = theExpr.expression().isConstant()? buildCstPart : buildNonCstPart;
                if ( doAssemblyRhs && buildSourceTerm )
                {
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= equationScaling*timeSteppingScaling*inner(theExpr,id(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
    }


    // additional term in regularized formulation
    if ( M_nullSpaceMethod == "regularized-formulation" && buildCstPart )
    {
        double epsilonPenal = M_nullSpaceRegularizationEpsilon;
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
    if ( buildNonCstPart )
    {
        double M_penaldir = 1e6;
        auto reluctivityExpr = this->reluctivityExpr( symbolsExpr );
        for ( auto const& [bcName,bcData] : M_boundaryConditions->magneticPotentialImposed() )
        {
            if ( !bcData->isMethodNitsche() )
                continue;
            auto theExpr = bcData->expr( symbolsExpr );
            if ( doAssemblyRhs )
            {
                if constexpr (  nDim == 3 ) // TODO 2D cases
                {
                    myLinearForm +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*equationScaling*inner(curl(v),reluctivityExpr*cross(theExpr,N()))
                                   + timeSteppingScaling*equationScaling*M_penaldir*trans( reluctivityExpr*cross(theExpr,N()) )*cross(id(v),N())/hFace(),
                                   _geomap=this->geomap() );
                }
            }
        }

        if ( doAssemblyRhs )
        {
            std::set<std::string> allmarkers;
            for ( auto const& [bcId,bcData] : M_boundaryConditions->anyBcWithVectorPotentialImposed( boundary_conditions_type::vector_potential_imposed_base_type::Method::nitsche ) )
                allmarkers.insert( bcData->markers().begin(), bcData->markers().end() );

            if constexpr (  nDim == 3 ) // TODO 2D cases
            {
                if ( !allmarkers.empty() )
                    bilinearForm_A_A +=
                        integrate(_range=markedfaces(this->mesh(),allmarkers),
                                  _expr=-equationScaling*timeSteppingScaling*inner(reluctivityExpr*curlt(u),cross(id(v),N()))
                                  - equationScaling*timeSteppingScaling*inner(curl(v),reluctivityExpr*cross(idt(u),N()))
                                  + equationScaling*timeSteppingScaling*M_penaldir*trans( reluctivityExpr*cross(idt(u),N()) )*cross(id(v),N())/hFace(),
                                  _geomap=this->geomap() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------//


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
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;
    this->log("Magnetic","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto mesh = this->mesh();
    size_type startBlockIndexVectorPotential = this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() );
    auto u = this->spaceVectorPotential()->element( U, this->rowStartInVector()+startBlockIndexVectorPotential );
    auto const& se = mctx.symbolsExpr();

    M_boundaryConditions->applyNewtonInitialGuess( mesh, u, se );
    // update info for synchronization
    this->updateDofEliminationIds( FieldTag::vectorPotential(this).identifierString(), data );

    if ( M_nullSpaceMethod == "saddle-point" )
    {
        size_type startBlockIndexLmCoulombGauge = this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() );
        auto p = this->spaceLagrangeMultiplierCoulombGauge()->element( U, this->rowStartInVector()+startBlockIndexLmCoulombGauge );
        p.on( _range=boundaryfaces( support( this->spaceLagrangeMultiplierCoulombGauge() ) ), _expr=cst(0.), _close=false );
        this->updateDofEliminationIds( FieldTag::lagrangeMultiplierCoulombGauge(this).identifierString(), data );
    }

    this->log("Magnetic","updateNewtonInitialGuess","finish" );
}

template< typename ConvexType, typename BasisMagneticVectorPotentialType >
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool _BuildCstPart = data.buildCstPart();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Magnetic","updateJacobian", "start"+sc);
    this->timerTool("Solve").start();
#if 0
    bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
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

    auto mesh = this->mesh();
    auto Xh = this->spaceVectorPotential();
    auto const& u = mctx.field( FieldTag::vectorPotential(this), FieldTag::vectorPotential(this).identifier() );
    auto const& v = this->fieldVectorPotential();
    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();
    size_type startBlockIndexVectorPotential = this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() );

    auto bilinearForm_A_A = form2( _test=Xh,_trial=Xh,_matrix=J,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix()+startBlockIndexVectorPotential,
                                   _colstart=this->colStartInMatrix()+startBlockIndexVectorPotential );



    double mu_0_cst = ModelPhysicMagnetic<nDim>::vacuumPermeabilityConstant();
    double equationScaling = M_equationScalingUseVacuumPermeability? mu_0_cst : 1.0;

    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicMagneticData = std::static_pointer_cast<ModelPhysicMagnetic<nDim>>(physicData);
        auto mu_0 = physicMagneticData->vacuumPermeabilityExpr();
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto const& magneticRelativePermeability = this->materialsProperties()->materialProperty( matName, "magnetic-relative-permeability" );
            bool magneticRelativePermeabilityDependOnTrialSymbol = magneticRelativePermeability.hasSymbolDependency( trialSymbolNames,se );

            if ( magneticRelativePermeability.isMatrix() )
            {
                if constexpr (  nDim == 3 )
                {
                    auto mu_r = expr( magneticRelativePermeability.template expr<nDim,nDim>(), se );
                    bool buildRotRot = mu_r.expression().isEvaluable()? buildCstPart : buildNonCstPart;
                    if ( buildRotRot )
                    {
                        bilinearForm_A_A +=
                            integrate( _range=range,
                                       _expr= equationScaling*timeSteppingScaling*(1./mu_0)*inner(inv(mu_r)*curlt(u),curl(v)),
                                       _geomap=this->geomap() );
                    }
                }
                else
                    CHECK( false ) << "material property magnetic-relative-permeability should be a scalar in 2D";

                if ( magneticRelativePermeabilityDependOnTrialSymbol && buildNonCstPart )
                {
                    throw std::runtime_error( "TODO: magneticRelativePermeabilityDependOnTrialSymbol aniso" );
                }
            }
            else
            {
                auto mu_r = expr( magneticRelativePermeability.expr(), se );
                bool buildRotRot = mu_r.expression().isEvaluable()? buildCstPart : buildNonCstPart;
                if ( buildRotRot )
                {
                    bilinearForm_A_A +=
                        integrate( _range=range,
                                   _expr= equationScaling*timeSteppingScaling*(1./(mu_0*mu_r))*inner(curlt(u),curl(v)),
                                   _geomap=this->geomap() );
                }

                if ( magneticRelativePermeabilityDependOnTrialSymbol && buildNonCstPart )
                {
                    throw std::runtime_error( "TODO: magneticRelativePermeabilityDependOnTrialSymbol" );
                }
            }
        }
    }


        // additional term in regularized formulation
    if ( M_nullSpaceMethod == "regularized-formulation" && buildCstPart )
    {
        double epsilonPenal = M_nullSpaceRegularizationEpsilon;
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
        auto bilinearForm_lm_A = form2( _test=XhLm,_trial=Xh,_matrix=J,
                                        _pattern=size_type(Pattern::COUPLED),
                                        _rowstart=this->rowStartInMatrix()+startBlockIndexLmCoulombGauge,
                                        _colstart=this->colStartInMatrix()+startBlockIndexVectorPotential );
        auto bilinearForm_A_lm = form2( _test=Xh,_trial=XhLm,_matrix=J,
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


    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Magnetic","updateJacobian",
              fmt::format("finish {} in {} s",sc,timeElapsed) );
}

template< typename ConvexType, typename BasisMagneticVectorPotentialType >
template <typename ModelContextType>
void
Magnetic<ConvexType,BasisMagneticVectorPotentialType>::updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Magnetic","updateResidual", "start"+sc);
    this->timerTool("Solve").start();
#if 0
    bool Build_TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        Build_TransientTerm=buildNonCstPart && !UseJacobianLinearTerms;
#endif
    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationary() )
    {
#if 0
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm( this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( M_timeStepping == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - M_timeStepThetaValue;
            else
                timeSteppingScaling = M_timeStepThetaValue;
        }
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
#endif
    }



    auto mesh = this->mesh();
    auto Xh = this->spaceVectorPotential();
    auto const& u = mctx.field( FieldTag::vectorPotential(this), FieldTag::vectorPotential(this).identifier() );
    auto const& v = this->fieldVectorPotential();
    auto const& se = mctx.symbolsExpr();
    auto const& tse = mctx.trialSymbolsExpr();
    auto trialSymbolNames = tse.names();

    size_type startBlockIndexVectorPotential = this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() );

    auto linearForm_A = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() + startBlockIndexVectorPotential );

    double mu_0_cst = ModelPhysicMagnetic<nDim>::vacuumPermeabilityConstant();
    double equationScaling = M_equationScalingUseVacuumPermeability? mu_0_cst : 1.0;

    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicMagneticData = std::static_pointer_cast<ModelPhysicMagnetic<nDim>>(physicData);
        auto mu_0 = physicMagneticData->vacuumPermeabilityExpr();
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto const& magneticRelativePermeability = this->materialsProperties()->materialProperty( matName, "magnetic-relative-permeability" );
            if ( magneticRelativePermeability.isMatrix() )
            {
                if constexpr (  nDim == 3 )
                {
                    auto mu_r = expr( magneticRelativePermeability.template expr<nDim,nDim>(), se );
                    bool buildRotRot = mu_r.expression().isEvaluable()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;

                    if ( buildRotRot )
                    {
                        linearForm_A +=
                            integrate( _range=range,
                                       _expr= equationScaling*timeSteppingScaling*(1./mu_0)*inner(inv(mu_r)*curlv(u),curl(v)),
                                       _geomap=this->geomap() );
                    }
                }
                else
                    CHECK( false ) << "material property magnetic-relative-permeability should be a scalar in 2D";
            }
            else
            {
                auto mu_r = expr( magneticRelativePermeability.expr(), se );
                bool buildRotRot = mu_r.expression().isEvaluable()? buildNonCstPart && !UseJacobianLinearTerms : buildNonCstPart;
                if ( buildRotRot )
                {
                    linearForm_A +=
                        integrate( _range=range,
                                   _expr= equationScaling*timeSteppingScaling*(1./(mu_0*mu_r))*inner(curlv(u),curl(v)),
                                   _geomap=this->geomap() );
                }
            }

            // current density sources
            for ( auto const& currentDensitySource : physicMagneticData->currentDensitySources() )
            {
                auto theExpr = currentDensitySource.expr( se );
                bool buildSourceTerm = theExpr.expression().isEvaluable()? buildCstPart : buildNonCstPart;
                if ( buildSourceTerm )
                {
                    linearForm_A +=
                        integrate( _range=range,
                                   _expr= -equationScaling*timeSteppingScaling*inner(theExpr,id(v)),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    // additional term in regularized formulation
    if ( M_nullSpaceMethod == "regularized-formulation" && buildNonCstPart && !UseJacobianLinearTerms )
    {
        double epsilonPenal = M_nullSpaceRegularizationEpsilon;
        linearForm_A +=
            integrate( _range=this->rangeMeshElements(),
                       _expr= timeSteppingScaling*epsilonPenal*inner(idv(u),id(v)),
                       _geomap=this->geomap() );
    }

    if ( M_nullSpaceMethod == "saddle-point" && buildNonCstPart && !UseJacobianLinearTerms )
    {
        auto XhLm = this->spaceLagrangeMultiplierCoulombGauge();
        auto const& p = mctx.field( FieldTag::lagrangeMultiplierCoulombGauge(this), FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() );
        size_type startBlockIndexLmCoulombGauge = this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() );
        auto linearForm_lm = form1( _test=XhLm, _vector=R,
                                    _rowstart=this->rowStartInVector() + startBlockIndexLmCoulombGauge );

        linearForm_A +=
            integrate( _range=this->rangeMeshElements(),
                       _expr= timeSteppingScaling*inner(id(v),trans(gradv(p))),
                       _geomap=this->geomap() );

        linearForm_lm +=
            integrate( _range=this->rangeMeshElements(),
                       _expr= timeSteppingScaling*inner(idv(u),trans(grad(p))),
                       _geomap=this->geomap() );
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Magnetic","updateResidual",
              fmt::format("finish {} in {} s",sc,timeElapsed) );
}



} // namespace Feel
} // namespace FeelModels

#endif
