/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_STABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_STABILIZATIONGLS_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

#include <feel/feelmodels/modelvf/exproptionalconcat.hpp>
#include <feel/feelmodels/modelvf/shockcapturing.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType,typename RangeType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateLinearPDEStabilizationGLS( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const
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
    auto const& u = this->fieldUnknown();
    auto const& v = this->fieldUnknown();

    auto const& se = mctx.symbolsExpr();

    bool hasConvectionTerm = this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() );
    bool hasDiffusionTerm = this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() );
    bool hasReactionTerm = this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_convection_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() ).template expr<nDim,1>(), se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).template expr<nDim,nDim>(), se ) )>;

    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradt(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idt(u) )>;
    using expr_diffusion_scalar_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplaciant(u) )>;
    using expr_diffusion_matrix_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                       std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hesst(u)) )>,
                                                                       std::decay_t<decltype( cst(0.) )> >::type;
    auto residual_lhs = exprOptionalConcat<expr_convection_residual_lhs_type,expr_reaction_residual_lhs_type,expr_diffusion_scalar_residual_lhs_type,expr_diffusion_matrix_residual_lhs_type>();

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
        auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
        auto const& coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
        residual_lhs.expression().template set<0>( timeSteppingScaling*(gradt(u)*coeff_beta_expr) );
        stab_test.expression().template set<0>( grad(u)*coeff_beta_expr );
        tauExprScalarDiffusion.expression().setConvection( coeff_beta_expr );
        tauExprMatrixDiffusion.expression().setConvection( coeff_beta_expr );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
        auto const& coeff_a_expr = expr( coeff_a.expr(), se );
        residual_lhs.expression().template set<1>( timeSteppingScaling*coeff_a_expr*idt(u) );
        if ( coeffNatureStabilization != 0 )
            stab_test.expression().template set<1>( coeffNatureStabilization*coeff_a_expr*id(u) );
        tauExprScalarDiffusion.expression().setReaction( coeff_a_expr );
        tauExprMatrixDiffusion.expression().setReaction( coeff_a_expr );
    }

    // diffusion ( TODO : this part is not always correct if diffusion term k depend on x,y,z (i.e. div( k /nabla v ) != k laplacian v) )
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_lhs.expression().template set<3>( -timeSteppingScaling*inner(coeff_c_expr,hesst(u)) );
                    if ( coeffNatureStabilization != 0 )
                        stab_test.expression().template set<3>( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
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
                residual_lhs.expression().template set<2>( -timeSteppingScaling*coeff_c_expr*laplaciant(u) );
                if ( coeffNatureStabilization != 0 )
                    stab_test.expression().template set<2>( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
            }
            tauExprScalarDiffusion.expression().setDiffusion( coeff_c_expr );
        }
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();

    if ( hasDiffusionTerm && tauExprMatrixDiffusion.expression().hasDiffusion() )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);

    auto tau = idv(tauFieldPtr);

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    bilinearForm +=
        integrate( _range=range,
                   _expr=tau*inner(residual_lhs,stab_test),
                   _geomap=this->geomap() );

    // Source
    if ( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) )
    {
        auto const& coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );
        auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                           [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                           [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );

        auto linearForm = form1( _test=Xh, _vector=F,
                                 _rowstart=this->rowStartInVector() );
        linearForm +=
            integrate( _range=range,
                       _expr= tau*inner(timeSteppingScaling*coeff_f_expr,stab_test),
                       _geomap=this->geomap() );

    }

}


template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType,typename RangeType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateJacobianStabilizationGLS( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const
{
    sparse_matrix_ptrtype& J = data.jacobian();
    double timeSteppingScaling = 1.; // TODO

    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    bool hasConvectionTerm = this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() );
    bool hasDiffusionTerm = this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() );
    bool hasReactionTerm = this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_convection_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() ).template expr<nDim,1>(), se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).template expr<nDim,nDim>(), se ) )>;

    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradt(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idt(u) )>;
    using expr_diffusion_scalar_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplaciant(u) )>;
    using expr_diffusion_matrix_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                       std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hesst(u)) )>,
                                                                       std::decay_t<decltype( cst(0.) )> >::type;
    auto residual_lhs = exprOptionalConcat<expr_convection_residual_lhs_type,expr_reaction_residual_lhs_type,expr_diffusion_scalar_residual_lhs_type,expr_diffusion_matrix_residual_lhs_type>();

    using expr_convection_test_type = std::decay_t<decltype( grad(u)*expr_coeff_convection_type{} )>;
    using expr_reaction_test_type = std::decay_t<decltype( coeffNatureStabilization*expr_coeff_reaction_type{}*id(u) )>;
    using expr_diffusion_scalar_test_type = std::decay_t<decltype( -coeffNatureStabilization*expr_coeff_diffusion_scalar_type{}*laplacian(u) )>;
    using expr_diffusion_matrix_test_type = typename mpl::if_c<unknown_is_scalar,
                                                               std::decay_t<decltype( -coeffNatureStabilization*inner(expr_coeff_diffusion_matrix_type{},hess(u)) )>,
                                                               std::decay_t<decltype( cst(0.) )> >::type;
    auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_reaction_test_type,expr_diffusion_scalar_test_type,expr_diffusion_matrix_test_type>();

    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_scalar_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_coeff_diffusion_matrix_type,expr_coeff_reaction_type>( *this->stabilizationGLSParameter() );


    using expr_convection_residual_eval_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradv(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_eval_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idv(u) )>;
    using expr_diffusion_scalar_residual_eval_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplacianv(u) )>;
    using expr_diffusion_matrix_residual_eval_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                            std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hessv(u)) )>,
                                                                            std::decay_t<decltype( cst(0.) )> >::type;
    auto residual_eval_lhs = exprOptionalConcat<expr_convection_residual_eval_lhs_type,expr_reaction_residual_eval_lhs_type,expr_diffusion_scalar_residual_eval_lhs_type,expr_diffusion_matrix_residual_eval_lhs_type>();
    bool useShockCapturing =  this->M_stabilizationGLS_applyShockCapturing && hasConvectionTerm;

    // convection
    if ( hasConvectionTerm )
    {
        auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
        auto const& coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
        residual_lhs.expression().template set<0>( timeSteppingScaling*(gradt(u)*coeff_beta_expr) );
        if ( useShockCapturing )
            residual_eval_lhs.expression().template set<0>( timeSteppingScaling*(gradv(u)*coeff_beta_expr) );
        stab_test.expression().template set<0>( grad(u)*coeff_beta_expr );
        tauExprScalarDiffusion.expression().setConvection( coeff_beta_expr );
        tauExprMatrixDiffusion.expression().setConvection( coeff_beta_expr );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
        auto const& coeff_a_expr = expr( coeff_a.expr(), se );
        residual_lhs.expression().template set<1>( timeSteppingScaling*coeff_a_expr*idt(u) );
        if ( useShockCapturing )
            residual_eval_lhs.expression().template set<1>( timeSteppingScaling*coeff_a_expr*idv(u) );
        if ( coeffNatureStabilization != 0 )
            stab_test.expression().template set<1>( coeffNatureStabilization*coeff_a_expr*id(u) );
        tauExprScalarDiffusion.expression().setReaction( coeff_a_expr );
        tauExprMatrixDiffusion.expression().setReaction( coeff_a_expr );
    }

    // diffusion ( TODO : this part is not always correct if diffusion term k depend on x,y,z (i.e. div( k /nabla v ) != k laplacian v) )
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_lhs.expression().template set<3>( -timeSteppingScaling*inner(coeff_c_expr,hesst(u)) );
                    if ( useShockCapturing )
                        residual_eval_lhs.expression().template set<3>( -timeSteppingScaling*inner(coeff_c_expr,hessv(u)) );
                    if ( coeffNatureStabilization != 0 )
                        stab_test.expression().template set<3>( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
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
                residual_lhs.expression().template set<2>( -timeSteppingScaling*coeff_c_expr*laplaciant(u) );
                if ( useShockCapturing )
                    residual_eval_lhs.expression().template set<2>( -timeSteppingScaling*coeff_c_expr*laplacianv(u) );
                if ( coeffNatureStabilization != 0 )
                    stab_test.expression().template set<2>( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
            }
            tauExprScalarDiffusion.expression().setDiffusion( coeff_c_expr );
        }
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

    if ( this->M_stabilizationGLS_applyShockCapturing && hasConvectionTerm )
    {
        auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
        auto const& coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
        if constexpr ( unknown_is_scalar )
        {
            // auto tau_sc_FieldPtr = this->stabilizationGLSParameter()->fieldTauPtr()->functionSpace()->elementPtr();
            // tau_sc_FieldPtr->on(_range=range,_expr=shockCapturingParameterExpr(residual_eval_lhs,tau,coeff_beta_expr,u));
            // auto tau_sc = idv(tau_sc_FieldPtr);

            auto sc = shockCapturingJacobianExpr(residual_eval_lhs,residual_lhs,tau,coeff_beta_expr,u);
            bilinearForm +=
                integrate( _range=range,
                           _expr=inner(sc, trans(grad(v))),
                           _geomap=this->geomap() );
        }

    }

}

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType,typename RangeType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateResidualStabilizationGLS( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const
{
    vector_ptrtype& R = data.residual();

    double timeSteppingScaling = 1.; // TODO

    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = mctx.field( FieldTag::unknown(this), this->unknownName() );
    auto const& v = this->fieldUnknown();

    bool hasConvectionTerm = this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() );
    bool hasDiffusionTerm = this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() );
    bool hasReactionTerm = this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() );
    bool hasSourceTerm = this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_convection_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() ).template expr<nDim,1>(), se ) )>;
    using expr_coeff_reaction_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).template expr<nDim,nDim>(), se ) )>;
    using expr_coeff_source_type = typename mpl::if_c<unknown_is_scalar,
                                                      std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() ).expr(), se ) )>,
        std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() ).template expr<nDim,1>(), se ) )> >::type;

    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradv(u)*expr_coeff_convection_type{}) )>;
    using expr_reaction_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idv(u) )>;
    using expr_diffusion_scalar_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplacianv(u) )>;
    using expr_diffusion_matrix_residual_lhs_type = typename mpl::if_c<unknown_is_scalar,
                                                                       std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hessv(u)) )>,
                                                                       std::decay_t<decltype( cst(0.) )> >::type;
    using expr_source_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_source_type{} )>;
    // TODO  time discretisation
    auto residual_full = exprOptionalConcat<expr_convection_residual_lhs_type,expr_reaction_residual_lhs_type,expr_diffusion_scalar_residual_lhs_type,expr_diffusion_matrix_residual_lhs_type,expr_source_residual_lhs_type>();

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
        auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
        auto const& coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
        residual_full.expression().template set<0>( timeSteppingScaling*(gradv(u)*coeff_beta_expr) );
        stab_test.expression().template set<0>( grad(u)*coeff_beta_expr );
        tauExprScalarDiffusion.expression().setConvection( coeff_beta_expr );
        tauExprMatrixDiffusion.expression().setConvection( coeff_beta_expr );
    }

    // reaction
    if ( hasReactionTerm )
    {
        auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
        auto const& coeff_a_expr = expr( coeff_a.expr(), se );
        residual_full.expression().template set<1>( timeSteppingScaling*coeff_a_expr*idv(u) );
        if ( coeffNatureStabilization != 0 )
            stab_test.expression().template set<1>( coeffNatureStabilization*coeff_a_expr*id(u) );
        tauExprScalarDiffusion.expression().setReaction( coeff_a_expr );
        tauExprMatrixDiffusion.expression().setReaction( coeff_a_expr );
    }

    // diffusion ( TODO : this part is not always correct if diffusion term k depend on x,y,z (i.e. div( k /nabla v ) != k laplacian v) )
    if ( hasDiffusionTerm )
    {
        auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            if constexpr ( unknown_is_scalar )
            {
                auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                if constexpr ( nOrderUnknown > 1 )
                {
                    residual_full.expression().template set<3>( -timeSteppingScaling*inner(coeff_c_expr,hessv(u)) );
                    if ( coeffNatureStabilization != 0 )
                        stab_test.expression().template set<3>( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
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
                residual_full.expression().template set<2>( -timeSteppingScaling*coeff_c_expr*laplacianv(u) );
                if ( coeffNatureStabilization != 0 )
                    stab_test.expression().template set<2>( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
            }
            tauExprScalarDiffusion.expression().setDiffusion( coeff_c_expr );
        }
    }

    if ( hasSourceTerm )
    {
        auto const& coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );
        auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                           [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                           [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
        residual_full.expression().template set<4>( -timeSteppingScaling*coeff_f_expr );
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
        auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
        auto const& coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
        if constexpr ( unknown_is_scalar )
        {
            // auto tau_sc_FieldPtr = this->stabilizationGLSParameter()->fieldTauPtr()->functionSpace()->elementPtr();
            // tau_sc_FieldPtr->on(_range=range,_expr=shockCapturingParameterExpr(residual_full,tau,coeff_beta_expr,u));
            // auto tau_sc = idv(tau_sc_FieldPtr);

            auto sc = shockCapturingExpr(residual_full,tau,coeff_beta_expr,u);
            linearForm +=
                integrate( _range=range,
                           _expr=inner(sc, trans(grad(v))),
                           _geomap=this->geomap() );
        }

    }
}

} // namespace FeelModels
} // namespace Feel

#endif
