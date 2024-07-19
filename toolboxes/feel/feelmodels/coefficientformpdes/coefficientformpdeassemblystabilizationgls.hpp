/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_STABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_STABILIZATIONGLS_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

#include <feel/feelvf/exproptionalconcat.hpp>

#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

#include <feel/feelmodels/modelvf/shockcapturing.hpp>

namespace Feel
{
namespace FeelModels
{

namespace CoefficientFormPDE_detail
{
template<typename CoefficientFormPDEType, typename FieldType, typename SymbolsExprType>
auto
exprResidual( CoefficientFormPDEType const& cfpde, std::shared_ptr<ModelPhysicCoefficientFormPDE<CoefficientFormPDEType::nDim>> physicCFPDEData, std::string const& matName,
              FieldType const& u, SymbolsExprType const& se, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    using Coefficient = typename CoefficientFormPDEType::Coefficient;
    constexpr uint16_type nDim = CoefficientFormPDEType::nDim;
    constexpr bool unknown_is_scalar = CoefficientFormPDEType::unknown_is_scalar;
    constexpr uint16_type nOrderUnknown = CoefficientFormPDEType::nOrderUnknown;

    bool hasConvectionTerm = physicCFPDEData->hasCoefficient( Coefficient::convection );
    bool hasDiffusionTerm = physicCFPDEData->hasCoefficient( Coefficient::diffusion );
    bool hasReactionTerm = physicCFPDEData->hasCoefficient( Coefficient::reaction );
    bool hasSourceTerm = physicCFPDEData->hasCoefficient( Coefficient::source );

    using expr_coeff_convection_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::diffusion, se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,nDim>( Coefficient::diffusion, se ) )>;
    using expr_coeff_first_time_derivative_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se ) )>;
    using expr_coeff_source_type = typename mpl::if_c<unknown_is_scalar,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::source, se ) )>,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::source, se ) )>  >::type;

    using expr_convection_residual_type = std::decay_t<decltype( timeSteppingScaling*(gradv(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idv(u) )>;
    using expr_diffusion_scalar_residual_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplacianv(u) )>;
    using expr_diffusion_part2_scalar_residual_type = typename mpl::if_c<unknown_is_scalar,
                                                                         std::decay_t<decltype( -timeSteppingScaling*inner( grad<nDim>(expr_coeff_diffusion_scalar_type{}), gradv(u) ) )>,
                                                                         std::decay_t<decltype( -timeSteppingScaling*gradv(u)*trans(grad<nDim>(expr_coeff_diffusion_scalar_type{})) )> >::type;
    using expr_diffusion_matrix_residual_type = typename mpl::if_c<unknown_is_scalar,
                                                                   std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hessv(u)) )>,
                                                                   std::decay_t<decltype( cst(0.) )> >::type;
    using expr_first_time_derivative_residual_type = std::decay_t<decltype( expr_coeff_first_time_derivative_type{}*(cfpde.timeStepBdfUnknown()->polyDerivCoefficient(0)*idv(u) - idv(cfpde.timeStepBdfUnknown()->polyDeriv()) ) )>;
    using expr_source_residual_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_source_type{} )>;
    auto residual_full = exprOptionalConcat<expr_convection_residual_type,expr_reaction_residual_type,
                                            expr_diffusion_scalar_residual_type,expr_diffusion_part2_scalar_residual_type,
                                            expr_diffusion_matrix_residual_type,
                                            expr_first_time_derivative_residual_type,expr_source_residual_type>();

    // convection
    if ( hasConvectionTerm )
    {
        auto coeff_beta_expr = physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se );
        residual_full.expression().add( timeSteppingScaling*(gradv(u)*coeff_beta_expr) );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto coeff_a_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se );
        residual_full.expression().add( timeSteppingScaling*coeff_a_expr*idv(u) );
    }

    // diffusion
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = physicCFPDEData->coefficient( Coefficient::diffusion );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_full.expression().add( -timeSteppingScaling*inner(coeff_c_expr,hessv(u)) );
                }
                if ( coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
                {
                    CHECK( false ) << "not implemented"; // diffusion term depend on x,y,z : div( k /nabla v ) != k laplacian v)
                }
            }
            else
                CHECK( false ) << "can not define stabilization with matrix diffusion coefficient and vectorial unknown";
        }
        else
        {
            auto coeff_c_expr = expr( coeff_c.expr(), se );
            if constexpr ( nOrderUnknown > 1 )
            {
                residual_full.expression().add( -timeSteppingScaling*coeff_c_expr*laplacianv(u) );
            }
            if ( cfpde.stabilizationDoAssemblyWithGradDiffusionCoeff() && coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
            {
                auto grad_coeff_c_expr = grad<nDim>(coeff_c_expr,"",cfpde.worldComm(),cfpde.repository().expr());
                if constexpr( unknown_is_scalar )
                    residual_full.expression().add( -timeSteppingScaling*inner( grad_coeff_c_expr, gradv(u) ) );
                else
                    residual_full.expression().add( -timeSteppingScaling*gradv(u)*trans(grad_coeff_c_expr) );
            }
        }
    }

    // first time derivative
    if ( !cfpde.isStationary() && !timeSteppingEvaluateResidualWithoutTimeDerivative && physicCFPDEData->hasCoefficient( Coefficient::firstTimeDerivative ) )
    {
        auto coeff_d_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se );
        residual_full.expression().add( coeff_d_expr*(cfpde.timeStepBdfUnknown()->polyDerivCoefficient(0)*idv(u) - idv(cfpde.timeStepBdfUnknown()->polyDeriv()) ) );
    }

    // source
    if ( hasSourceTerm )
    {
        auto const& coeff_f = physicCFPDEData->coefficient( Coefficient::source );
        auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                           [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                           [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
        residual_full.expression().add( -timeSteppingScaling*coeff_f_expr );
    }

    return residual_full;
}

} // namespace



template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType,typename RangeType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateLinearPDEStabilizationGLS( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx,
                                                                                  std::shared_ptr<ModelPhysicCoefficientFormPDE<nDim>> physicCFPDEData,
                                                                                  std::string const& matName, RangeType const& range ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("CoefficientFormPDE","updateLinearPDEStabilizationGLS", "start"+sc);

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    auto const& se = mctx.symbolsExpr();

    bool hasConvectionTerm = physicCFPDEData->hasCoefficient( Coefficient::convection );
    bool hasDiffusionTerm = physicCFPDEData->hasCoefficient( Coefficient::diffusion );
    bool hasReactionTerm = physicCFPDEData->hasCoefficient( Coefficient::reaction );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_convection_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::diffusion, se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,nDim>( Coefficient::diffusion, se ) )>;
    using expr_coeff_first_time_derivative_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se ) )>;
    using expr_coeff_source_type = typename mpl::if_c<unknown_is_scalar,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::source, se ) )>,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::source, se ) )>  >::type;

    // residual lhs
    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradt(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idt(u) )>;
    using expr_diffusion_scalar_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplaciant(u) )>;
    using expr_diffusion_part2_scalar_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                             std::decay_t<decltype( -timeSteppingScaling*inner( grad<nDim>(expr_coeff_diffusion_scalar_type{}), gradt(u) ) )>,
                                                                             std::decay_t<decltype( -timeSteppingScaling*gradt(u)*trans(grad<nDim>(expr_coeff_diffusion_scalar_type{})) )> >::type;
    using expr_diffusion_matrix_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                       std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hesst(u)) )>,
                                                                       std::decay_t<decltype( cst(0.) )> >::type;
    using expr_first_time_derivative_residual_lhs_type = std::decay_t<decltype( expr_coeff_first_time_derivative_type{}*M_bdfUnknown->polyDerivCoefficient(0)*idt(u) )>;
    auto residual_lhs = exprOptionalConcat<expr_convection_residual_lhs_type,expr_reaction_residual_lhs_type,
                                           expr_diffusion_scalar_residual_lhs_type,expr_diffusion_part2_scalar_residual_lhs_type,
                                           expr_diffusion_matrix_residual_lhs_type,
                                           expr_first_time_derivative_residual_lhs_type>();
    // residual rhs
    using expr_first_time_derivative_residual_rhs_type = std::decay_t<decltype( expr_coeff_first_time_derivative_type{}*idv(M_bdfUnknown->polyDeriv()) )>;
    using expr_source_residual_rhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_source_type{} )>;
    using expr_previous_residual_type = std::decay_t<decltype( -CoefficientFormPDE_detail::exprResidual( *this,physicCFPDEData,matName,u,se,1.0-timeSteppingScaling,true ) )>;
    auto residual_rhs = exprOptionalConcat<expr_first_time_derivative_residual_rhs_type,expr_source_residual_rhs_type,expr_previous_residual_type>();
    // test functions
    using expr_convection_test_type = std::decay_t<decltype( grad(u)*expr_coeff_convection_type{} )>;
    using expr_reaction_test_type = std::decay_t<decltype( coeffNatureStabilization*expr_coeff_reaction_type{}*id(u) )>;
    using expr_diffusion_scalar_test_type = std::decay_t<decltype( -coeffNatureStabilization*expr_coeff_diffusion_scalar_type{}*laplacian(u) )>;
    using expr_diffusion_matrix_test_type = typename mpl::if_c<unknown_is_scalar,
                                                               std::decay_t<decltype( -coeffNatureStabilization*inner(expr_coeff_diffusion_matrix_type{},hess(u)) )>,
                                                               std::decay_t<decltype( Feel::vf::zero<nDim,1>() /*cst(0.)*/ )> // TODO
                                                               >::type;
    using expr_diffusion_variable_scalar_test_type = typename mpl::if_c<unknown_is_scalar,
                                                                        std::decay_t<decltype( -inner( grad<nDim>(expr_coeff_diffusion_scalar_type{}),grad(u)) )>,
                                                                        std::decay_t<decltype( Feel::vf::zero<nDim,1>() )> // TODO
                                                                         >::type;
    auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_reaction_test_type,expr_diffusion_scalar_test_type,expr_diffusion_matrix_test_type,
                                        expr_diffusion_variable_scalar_test_type>();

    using expr_grad_coeff_diffusion_scalar_type = std::decay_t<decltype( grad<nDim>(expr_coeff_diffusion_scalar_type{}) )>;
    using expr_convection_from_diffusion_variable_scalar_type = std::decay_t<decltype( -trans(expr_grad_coeff_diffusion_scalar_type{}) )>;
    auto exprConvectionUseByTauScalarDiffusion = exprOptionalConcat<expr_coeff_convection_type,expr_convection_from_diffusion_variable_scalar_type>();
    using expr_convection_tau_scalar_diffusion_type = std::decay_t<decltype( exprConvectionUseByTauScalarDiffusion )>;


    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_convection_tau_scalar_diffusion_type/*expr_coeff_convection_type*/,expr_coeff_diffusion_scalar_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_matrix_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );

    // convection
    if ( hasConvectionTerm )
    {
        auto coeff_beta_expr = physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se );
        residual_lhs.expression().add( timeSteppingScaling*(gradt(u)*coeff_beta_expr) );
        stab_test.expression().add( grad(u)*coeff_beta_expr );
#if 0
        tauExprScalarDiffusion.expression().setConvection( coeff_beta_expr );
#else
        exprConvectionUseByTauScalarDiffusion.expression().add( coeff_beta_expr );
#endif
        tauExprMatrixDiffusion.expression().setConvection( coeff_beta_expr );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto coeff_a_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se );
        residual_lhs.expression().add( timeSteppingScaling*coeff_a_expr*idt(u) );
        if ( coeffNatureStabilization != 0 )
            stab_test.expression().add( coeffNatureStabilization*coeff_a_expr*id(u) );
        tauExprScalarDiffusion.expression().setReaction( coeff_a_expr );
        tauExprMatrixDiffusion.expression().setReaction( coeff_a_expr );
    }

    // diffusion
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = physicCFPDEData->coefficient( Coefficient::diffusion );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_lhs.expression().add( -timeSteppingScaling*inner(coeff_c_expr,hesst(u)) );
                    if ( coeffNatureStabilization != 0 )
                        stab_test.expression().add( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
                }
                if ( coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
                {
                    CHECK( false ) << "not implemented"; // diffusion term depend on x,y,z : div( k /nabla v ) != k laplacian v)
                }
                tauExprMatrixDiffusion.expression().setDiffusion( coeff_c_expr );
            }
            else
                CHECK( false ) << "can not define stabilization with matrix diffusion coefficient and vectorial unknown";
        }
        else
        {
            auto coeff_c_expr = expr( coeff_c.expr(), se );
            if constexpr ( nOrderUnknown > 1 )
            {
                residual_lhs.expression().add( -timeSteppingScaling*coeff_c_expr*laplaciant(u) );
                if ( coeffNatureStabilization != 0 )
                    stab_test.expression().add( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
            }
            if ( this->stabilizationDoAssemblyWithGradDiffusionCoeff() && coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
            {
                auto grad_coeff_c_expr = grad<nDim>(coeff_c_expr,"",this->worldComm(),this->repository().expr());
                if constexpr( unknown_is_scalar )
                {
                    residual_lhs.expression().add( -timeSteppingScaling*inner( grad_coeff_c_expr, gradt(u) ) );
                    stab_test.expression().add( -inner( grad_coeff_c_expr,grad(u) ) );
                    exprConvectionUseByTauScalarDiffusion.expression().add( -trans(grad_coeff_c_expr) );
                }
                else
                {
                    residual_lhs.expression().add( -timeSteppingScaling*gradt(u)*trans(grad_coeff_c_expr) );
                    // TODO stab_test
                }
            }
            tauExprScalarDiffusion.expression().setDiffusion( coeff_c_expr );
        }
    }

    // Source
    if ( physicCFPDEData->hasCoefficient( Coefficient::source ) )
    {
        auto const& coeff_f = physicCFPDEData->coefficient( Coefficient::source );
        auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                           [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                           [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );

        residual_rhs.expression().add( timeSteppingScaling*coeff_f_expr );
    }

    // first time derivative
    if ( !this->isStationary() && physicCFPDEData->hasCoefficient( Coefficient::firstTimeDerivative ) )
    {
        auto coeff_d_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se );
        residual_lhs.expression().add( coeff_d_expr*M_bdfUnknown->polyDerivCoefficient(0)*idt(u) );
        residual_rhs.expression().add( coeff_d_expr*idv(M_bdfUnknown->polyDeriv()) );
    }

    if ( !this->isStationary() && this->timeStepping() == "Theta" )
    {
        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
        auto const& uPrevious = mctxPrevious.field( FieldTag::unknown(this), this->unknownName() );
        auto const& sePrevious = mctxPrevious.symbolsExpr();
        auto previousResidualExpr = CoefficientFormPDE_detail::exprResidual( *this,physicCFPDEData,matName,uPrevious,sePrevious,1.0-timeSteppingScaling,true );
        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
            previousResidualExpr.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
        residual_rhs.expression().add( -previousResidualExpr );
    }

    if ( exprConvectionUseByTauScalarDiffusion.expression().hasExpr() )
    {
        //if ( thermalConductivityIsMatrix )
        //  tauExprMatrixDiffusion.expression().setConvection( exprCoeffConvection );
        //else
        tauExprScalarDiffusion.expression().setConvection( exprConvectionUseByTauScalarDiffusion );
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();

    if ( hasDiffusionTerm && tauExprMatrixDiffusion.expression().hasDiffusion() )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);

    auto tau = idv(tauFieldPtr);

    if ( residual_lhs.expression().hasExpr() )
    {
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        bilinearForm +=
            integrate( _range=range,
                       _expr=tau*inner(residual_lhs,stab_test),
                       _geomap=this->geomap() );
    }

    if ( residual_rhs.expression().hasExpr() )
    {
        auto linearForm = form1( _test=Xh, _vector=F,
                                 _rowstart=this->rowStartInVector() );
        linearForm +=
            integrate( _range=range,
                       _expr= tau*inner(residual_rhs,stab_test),
                       _geomap=this->geomap() );
    }

}


template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType,typename RangeType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateJacobianStabilizationGLS( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx,
                                                                                 std::shared_ptr<ModelPhysicCoefficientFormPDE<nDim>> physicCFPDEData,
                                                                                 std::string const& matName, RangeType const& range ) const
{
    sparse_matrix_ptrtype& J = data.jacobian();

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    bool hasConvectionTerm = physicCFPDEData->hasCoefficient( Coefficient::convection );
    bool hasDiffusionTerm = physicCFPDEData->hasCoefficient( Coefficient::diffusion );
    bool hasReactionTerm = physicCFPDEData->hasCoefficient( Coefficient::reaction );
    bool hasSourceTerm = physicCFPDEData->hasCoefficient( Coefficient::source );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_convection_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::diffusion, se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,nDim>( Coefficient::diffusion, se ) )>;
    using expr_coeff_first_time_derivative_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se ) )>;
    using expr_coeff_source_type = typename mpl::if_c<unknown_is_scalar,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::source, se ) )>,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::source, se ) )>  >::type;


    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradt(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idt(u) )>;
    using expr_diffusion_scalar_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplaciant(u) )>;
    using expr_diffusion_part2_scalar_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                             std::decay_t<decltype( -timeSteppingScaling*inner( grad<nDim>(expr_coeff_diffusion_scalar_type{}), gradt(u) ) )>,
                                                                             std::decay_t<decltype( -timeSteppingScaling*gradt(u)*trans(grad<nDim>(expr_coeff_diffusion_scalar_type{})) )> >::type;
    using expr_diffusion_matrix_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                       std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hesst(u)) )>,
                                                                       std::decay_t<decltype( cst(0.) )> >::type;
    using expr_first_time_derivative_residual_lhs_type = std::decay_t<decltype( expr_coeff_first_time_derivative_type{}*M_bdfUnknown->polyDerivCoefficient(0)*idt(u) )>;
    auto residual_lhs = exprOptionalConcat<expr_convection_residual_lhs_type,expr_reaction_residual_lhs_type,
                                           expr_diffusion_scalar_residual_lhs_type,expr_diffusion_part2_scalar_residual_lhs_type,
                                           expr_diffusion_matrix_residual_lhs_type,
                                           expr_first_time_derivative_residual_lhs_type>();

    using expr_convection_test_type = std::decay_t<decltype( grad(u)*expr_coeff_convection_type{} )>;
    using expr_reaction_test_type = std::decay_t<decltype( coeffNatureStabilization*expr_coeff_reaction_type{}*id(u) )>;
    using expr_diffusion_scalar_test_type = std::decay_t<decltype( -coeffNatureStabilization*expr_coeff_diffusion_scalar_type{}*laplacian(u) )>;
    using expr_diffusion_matrix_test_type = typename mpl::if_c<unknown_is_scalar,
                                                               std::decay_t<decltype( -coeffNatureStabilization*inner(expr_coeff_diffusion_matrix_type{},hess(u)) )>,
                                                               std::decay_t<decltype( cst(0.) )> >::type;
    auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_reaction_test_type,expr_diffusion_scalar_test_type,expr_diffusion_matrix_test_type>();

    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_scalar_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_matrix_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );

    // residual (evaluation)
    using expr_convection_residual_eval_type = std::decay_t<decltype( timeSteppingScaling*(gradv(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_eval_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idv(u) )>;
    using expr_diffusion_scalar_residual_eval_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplacianv(u) )>;
    using expr_diffusion_part2_scalar_residual_eval_type = typename mpl::if_c<unknown_is_scalar,
                                                                              std::decay_t<decltype( -timeSteppingScaling*inner( grad<nDim>(expr_coeff_diffusion_scalar_type{}), gradv(u) ) )>,
                                                                              std::decay_t<decltype( -timeSteppingScaling*gradv(u)*trans(grad<nDim>(expr_coeff_diffusion_scalar_type{})) )> >::type;
    using expr_diffusion_matrix_residual_eval_type = typename mpl::if_c<unknown_is_scalar,
                                                                            std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hessv(u)) )>,
                                                                            std::decay_t<decltype( cst(0.) )> >::type;
    using expr_first_time_derivative_residual_eval_type = std::decay_t<decltype( expr_coeff_first_time_derivative_type{}*(M_bdfUnknown->polyDerivCoefficient(0)*idv(u) - idv(M_bdfUnknown->polyDeriv()) ) )>;
    using expr_source_residual_eval_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_source_type{} )>;
    auto residual_eval = exprOptionalConcat<expr_convection_residual_eval_type,expr_reaction_residual_eval_type,
                                            expr_diffusion_scalar_residual_eval_type,expr_diffusion_part2_scalar_residual_eval_type,
                                            expr_diffusion_matrix_residual_eval_type,
                                            expr_first_time_derivative_residual_eval_type,expr_source_residual_eval_type>();
    bool useShockCapturing = this->M_stabilizationGLS_applyShockCapturing && hasConvectionTerm;


    // convection
    if ( hasConvectionTerm )
    {
        auto coeff_beta_expr = physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se );
        residual_lhs.expression().add( timeSteppingScaling*(gradt(u)*coeff_beta_expr) );
        if ( useShockCapturing )
            residual_eval.expression().add( timeSteppingScaling*(gradv(u)*coeff_beta_expr) );
        stab_test.expression().add( grad(u)*coeff_beta_expr );
        tauExprScalarDiffusion.expression().setConvection( coeff_beta_expr );
        tauExprMatrixDiffusion.expression().setConvection( coeff_beta_expr );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto coeff_a_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se );
        residual_lhs.expression().add( timeSteppingScaling*coeff_a_expr*idt(u) );
        if ( useShockCapturing )
            residual_eval.expression().add( timeSteppingScaling*coeff_a_expr*idv(u) );
        if ( coeffNatureStabilization != 0 )
            stab_test.expression().add( coeffNatureStabilization*coeff_a_expr*id(u) );
        tauExprScalarDiffusion.expression().setReaction( coeff_a_expr );
        tauExprMatrixDiffusion.expression().setReaction( coeff_a_expr );
    }

    // diffusion
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = physicCFPDEData->coefficient( Coefficient::diffusion );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_lhs.expression().add( -timeSteppingScaling*inner(coeff_c_expr,hesst(u)) );
                    if ( useShockCapturing )
                        residual_eval.expression().add( -timeSteppingScaling*inner(coeff_c_expr,hessv(u)) );
                    if ( coeffNatureStabilization != 0 )
                        stab_test.expression().add( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
                }
                if ( coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
                {
                    CHECK( false ) << "not implemented"; // diffusion term depend on x,y,z : div( k /nabla v ) != k laplacian v)
                }
                tauExprMatrixDiffusion.expression().setDiffusion( coeff_c_expr );
            }
            else
                CHECK( false ) << "can not define stabilization with matrix diffusion coefficient and vectorial unknown";
        }
        else
        {
            auto const& coeff_c_expr = expr( coeff_c.expr(), se );
            if constexpr ( nOrderUnknown > 1 )
            {
                residual_lhs.expression().add( -timeSteppingScaling*coeff_c_expr*laplaciant(u) );
                if ( useShockCapturing )
                    residual_eval.expression().add( -timeSteppingScaling*coeff_c_expr*laplacianv(u) );
                if ( coeffNatureStabilization != 0 )
                    stab_test.expression().add( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
            }
            if ( this->stabilizationDoAssemblyWithGradDiffusionCoeff() && coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
            {
                auto grad_coeff_c_expr = grad<nDim>(coeff_c_expr,"",this->worldComm(),this->repository().expr());
                if constexpr( unknown_is_scalar )
                {
                    residual_lhs.expression().add( -timeSteppingScaling*inner( grad_coeff_c_expr, gradt(u) ) );
                    if ( useShockCapturing )
                        residual_eval.expression().add( -timeSteppingScaling*inner( grad_coeff_c_expr, gradv(u) ) );
                }
                else
                {
                    residual_lhs.expression().add( -timeSteppingScaling*gradt(u)*trans(grad_coeff_c_expr) );
                    if ( useShockCapturing )
                        residual_eval.expression().add( -timeSteppingScaling*gradv(u)*trans(grad_coeff_c_expr) );
                }
            }
            tauExprScalarDiffusion.expression().setDiffusion( coeff_c_expr );
        }
    }

    // first time derivative
    if ( !this->isStationary() && physicCFPDEData->hasCoefficient( Coefficient::firstTimeDerivative ) )
    {
        auto coeff_d_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se );
        residual_lhs.expression().add( coeff_d_expr*M_bdfUnknown->polyDerivCoefficient(0)*idt(u) );
        if ( useShockCapturing )
            residual_eval.expression().add( coeff_d_expr*(M_bdfUnknown->polyDerivCoefficient(0)*idv(u) - idv(M_bdfUnknown->polyDeriv()) ) );
    }


    // source
    if ( hasSourceTerm )
    {
        auto const& coeff_f = physicCFPDEData->coefficient( Coefficient::source );
        auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                           [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                           [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
        residual_eval.expression().add( -timeSteppingScaling*coeff_f_expr );
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();

    if ( hasDiffusionTerm && tauExprMatrixDiffusion.expression().hasDiffusion() )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);

#if 1
    auto tau = idv(tauFieldPtr);
#else
    auto tau = tauExprScalarDiffusion;
#endif

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    bilinearForm +=
        integrate( _range=range,
                   _expr=tau*inner(residual_lhs,stab_test),
                   _geomap=this->geomap() );

    if ( useShockCapturing )
    {
        auto coeff_beta_expr = physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se );
        if constexpr ( unknown_is_scalar )
        {
            // auto tau_sc_FieldPtr = this->stabilizationGLSParameter()->fieldTauPtr()->functionSpace()->elementPtr();
            // tau_sc_FieldPtr->on(_range=range,_expr=shockCapturingParameterExpr(residual_eval,tau,coeff_beta_expr,u));
            // auto tau_sc = idv(tau_sc_FieldPtr);
            uint16_type quadUsed = ioption(_name="stabilization.gls.shock-capturing.quad",_prefix=this->prefix(),_vm=this->clovm());
            auto sc = shockCapturingJacobianExpr(residual_eval,residual_lhs,tau,coeff_beta_expr,u);
            bilinearForm +=
                integrate( _range=range,
                           _expr=inner(sc, trans(grad(v))),
                           _quad= (quadUsed>=0)? quadUsed : quad_order_from_expression,
                           _geomap=this->geomap() );
        }

    }

}

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType,typename RangeType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateResidualStabilizationGLS( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx,
                                                                                 std::shared_ptr<ModelPhysicCoefficientFormPDE<nDim>> physicCFPDEData,
                                                                                 std::string const& matName, RangeType const& range ) const
{
    vector_ptrtype& R = data.residual();

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo("time-stepping.evaluate-residual-without-time-derivative");
    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
        return;

    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    bool hasConvectionTerm = physicCFPDEData->hasCoefficient( Coefficient::convection );
    bool hasDiffusionTerm = physicCFPDEData->hasCoefficient( Coefficient::diffusion );
    bool hasReactionTerm = physicCFPDEData->hasCoefficient( Coefficient::reaction );
    bool hasSourceTerm = physicCFPDEData->hasCoefficient( Coefficient::source );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_convection_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::diffusion, se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,nDim>( Coefficient::diffusion, se ) )>;
    using expr_coeff_first_time_derivative_type = std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se ) )>;
    using expr_coeff_source_type = typename mpl::if_c<unknown_is_scalar,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<1,1>( Coefficient::source, se ) )>,
                                                      std::decay_t<decltype( physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::source, se ) )>  >::type;

    using expr_convection_residual_type = std::decay_t<decltype( timeSteppingScaling*(gradv(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idv(u) )>;
    using expr_diffusion_scalar_residual_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplacianv(u) )>;
    using expr_diffusion_part2_scalar_residual_type = typename mpl::if_c<unknown_is_scalar,
                                                                         std::decay_t<decltype( -timeSteppingScaling*inner( grad<nDim>(expr_coeff_diffusion_scalar_type{}), gradv(u) ) )>,
                                                                         std::decay_t<decltype( -timeSteppingScaling*gradv(u)*trans(grad<nDim>(expr_coeff_diffusion_scalar_type{})) )> >::type;
    using expr_diffusion_matrix_residual_type = typename mpl::if_c<unknown_is_scalar,
                                                                   std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hessv(u)) )>,
                                                                   std::decay_t<decltype( cst(0.) )> >::type;
    using expr_first_time_derivative_residual_type = std::decay_t<decltype( expr_coeff_first_time_derivative_type{}*(M_bdfUnknown->polyDerivCoefficient(0)*idv(u) - idv(M_bdfUnknown->polyDeriv()) ) )>;
    using expr_source_residual_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_source_type{} )>;
    using expr_previous_residual_type = std::decay_t<decltype( CoefficientFormPDE_detail::exprResidual( *this,physicCFPDEData,matName,u,se,1.0-timeSteppingScaling,true ) )>;
    auto residual_full = exprOptionalConcat<expr_convection_residual_type,expr_reaction_residual_type,
                                            expr_diffusion_scalar_residual_type,expr_diffusion_part2_scalar_residual_type,
                                            expr_diffusion_matrix_residual_type,
                                            expr_first_time_derivative_residual_type,expr_source_residual_type,
                                            expr_previous_residual_type>();

    using expr_convection_test_type = std::decay_t<decltype( grad(u)*expr_coeff_convection_type{} )>;
    using expr_reaction_test_type = std::decay_t<decltype( coeffNatureStabilization*expr_coeff_reaction_type{}*id(u) )>;
    using expr_diffusion_scalar_test_type = std::decay_t<decltype( -coeffNatureStabilization*expr_coeff_diffusion_scalar_type{}*laplacian(u) )>;
    using expr_diffusion_matrix_test_type = typename mpl::if_c<unknown_is_scalar,
                                                               std::decay_t<decltype( -coeffNatureStabilization*inner(expr_coeff_diffusion_matrix_type{},hess(u)) )>,
                                                               std::decay_t<decltype( cst(0.) )> >::type;
    auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_reaction_test_type,expr_diffusion_scalar_test_type,expr_diffusion_matrix_test_type>();

    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_scalar_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_matrix_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );

    // convection
    if ( hasConvectionTerm )
    {
        auto coeff_beta_expr = physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se );
        residual_full.expression().add( timeSteppingScaling*(gradv(u)*coeff_beta_expr) );
        stab_test.expression().add( grad(u)*coeff_beta_expr );
        tauExprScalarDiffusion.expression().setConvection( coeff_beta_expr );
        tauExprMatrixDiffusion.expression().setConvection( coeff_beta_expr );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto coeff_a_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::reaction, se );
        residual_full.expression().add( timeSteppingScaling*coeff_a_expr*idv(u) );
        if ( coeffNatureStabilization != 0 )
            stab_test.expression().add( coeffNatureStabilization*coeff_a_expr*id(u) );
        tauExprScalarDiffusion.expression().setReaction( coeff_a_expr );
        tauExprMatrixDiffusion.expression().setReaction( coeff_a_expr );
    }

    // diffusion
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = physicCFPDEData->coefficient( Coefficient::diffusion );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_full.expression().add( -timeSteppingScaling*inner(coeff_c_expr,hessv(u)) );
                    if ( coeffNatureStabilization != 0 )
                        stab_test.expression().add( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
                }
                if ( coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
                {
                    CHECK( false ) << "not implemented"; // diffusion term depend on x,y,z : div( k /nabla v ) != k laplacian v)
                }
                tauExprMatrixDiffusion.expression().setDiffusion( coeff_c_expr );
            }
            else
                CHECK( false ) << "can not define stabilization with matrix diffusion coefficient and vectorial unknown";
        }
        else
        {
            auto const& coeff_c_expr = expr( coeff_c.expr(), se );
            if constexpr ( nOrderUnknown > 1 )
            {
                residual_full.expression().add( -timeSteppingScaling*coeff_c_expr*laplacianv(u) );
                if ( coeffNatureStabilization != 0 )
                    stab_test.expression().add( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
            }
            if ( this->stabilizationDoAssemblyWithGradDiffusionCoeff() && coeff_c_expr.template hasSymbolDependencyOnCoordinatesInSpace<nDim>() )
            {
                auto grad_coeff_c_expr = grad<nDim>(coeff_c_expr,"",this->worldComm(),this->repository().expr());
                if constexpr( unknown_is_scalar )
                    residual_full.expression().add( -timeSteppingScaling*inner( grad_coeff_c_expr, gradv(u) ) );
                else
                    residual_full.expression().add( -timeSteppingScaling*gradv(u)*trans(grad_coeff_c_expr) );
            }

            tauExprScalarDiffusion.expression().setDiffusion( coeff_c_expr );
        }
    }

    // first time derivative
    if ( !this->isStationary() && !timeSteppingEvaluateResidualWithoutTimeDerivative && physicCFPDEData->hasCoefficient( Coefficient::firstTimeDerivative ) )
    {
        auto coeff_d_expr = physicCFPDEData->template coefficientExpr<1,1>( Coefficient::firstTimeDerivative, se );
        residual_full.expression().add( coeff_d_expr*(M_bdfUnknown->polyDerivCoefficient(0)*idv(u) - idv(M_bdfUnknown->polyDeriv()) ) );
    }

    // source
    if ( hasSourceTerm )
    {
        auto const& coeff_f = physicCFPDEData->coefficient( Coefficient::source );
        auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                           [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                           [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
        residual_full.expression().add( -timeSteppingScaling*coeff_f_expr );
    }

    if ( !this->isStationary() && this->timeStepping() == "Theta" )
    {
        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
        auto const& uPrevious = mctxPrevious.field( FieldTag::unknown(this), this->unknownName() );
        auto const& sePrevious = mctxPrevious.symbolsExpr();
        auto previousResidualExpr = CoefficientFormPDE_detail::exprResidual( *this,physicCFPDEData,matName,uPrevious,sePrevious,1.0-timeSteppingScaling,true );
        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
            previousResidualExpr.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
        residual_full.expression().add( previousResidualExpr );
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();

    if ( hasDiffusionTerm && tauExprMatrixDiffusion.expression().hasDiffusion() )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);
#if 1
    auto tau = idv(tauFieldPtr);
#else
    auto tau = tauExprScalarDiffusion;
#endif

    auto linearForm = form1( _test=Xh, _vector=R,
                             _rowstart=this->rowStartInVector() );

    linearForm +=
        integrate( _range=range,
                   _expr=tau*inner(residual_full,stab_test),
                   _geomap=this->geomap() );


    if ( this->M_stabilizationGLS_applyShockCapturing && hasConvectionTerm )
    {
        auto coeff_beta_expr = physicCFPDEData->template coefficientExpr<nDim,1>( Coefficient::convection, se );
        if constexpr ( unknown_is_scalar )
        {
            // auto tau_sc_FieldPtr = this->stabilizationGLSParameter()->fieldTauPtr()->functionSpace()->elementPtr();
            // tau_sc_FieldPtr->on(_range=range,_expr=shockCapturingParameterExpr(residual_full,tau,coeff_beta_expr,u));
            // auto tau_sc = idv(tau_sc_FieldPtr);
            uint16_type quadUsed = ioption(_name="stabilization.gls.shock-capturing.quad",_prefix=this->prefix(),_vm=this->clovm());
            auto sc = shockCapturingExpr(residual_full,tau,coeff_beta_expr,u);
            linearForm +=
                integrate( _range=range,
                           _expr=inner(sc, trans(grad(v))),
                           _quad= (quadUsed>=0)? quadUsed : quad_order_from_expression,
                           _geomap=this->geomap() );
        }

    }
}

} // namespace FeelModels
} // namespace Feel

#endif
