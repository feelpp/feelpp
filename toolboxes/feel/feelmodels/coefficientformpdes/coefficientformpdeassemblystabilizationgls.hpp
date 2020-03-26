/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_STABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_STABILIZATIONGLS_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

#include <feel/feelmodels/modelvf/exproptionalconcat.hpp>

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
    CHECK( hasConvectionTerm ) << "stab onlyilization if has convection term";
    bool hasDiffusionTerm = this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() );
    bool hasReactionTerm = this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() );

    auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
    auto const& coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );

    int coeffNatureStabilization = ( this->stabilizationType() == "supg" )? 0 : (this->stabilizationType() == "gls")? 1 : -1;

    using expr_coeff_reaction_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_scalar_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).expr(), se ) )>;
    using expr_coeff_diffusion_matrix_type = std::decay_t<decltype( expr( this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).template expr<nDim,nDim>(), se ) )>;

    using expr_convection_test_type = std::decay_t<decltype( grad(u)*coeff_beta_expr )>;
    using expr_reaction_test_type = std::decay_t<decltype( coeffNatureStabilization*expr_coeff_reaction_type{}*id(u) )>;
    using expr_diffusion_scalar_test_type = std::decay_t<decltype( -coeffNatureStabilization*expr_coeff_diffusion_scalar_type{}*laplacian(u) )>;
    using expr_diffusion_matrix_test_type = std::decay_t<decltype( -coeffNatureStabilization*inner(expr_coeff_diffusion_matrix_type{},hess(u)) )>;

    auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_reaction_test_type,expr_diffusion_scalar_test_type,expr_diffusion_matrix_test_type>();

    // convection
    stab_test.expression().template set<0>( grad(u)*coeff_beta_expr );

    if ( coeffNatureStabilization !=0 )
    {
        // reaction
        if ( hasReactionTerm )
        {
            auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
            auto const& coeff_a_expr = expr( coeff_a.expr(), se );
            stab_test.expression().template set<1>( coeffNatureStabilization*coeff_a_expr*id(u) );
        }

        // diffusion
        if ( hasDiffusionTerm )
        {
            // TODO : this part is not always correct if diffusion term k depend of x,y,z (i.e. div( k /nabla v ) != k laplacian v)
            if constexpr ( nOrderUnknown > 1 )
            {
                auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
                if ( coeff_c.template hasExpr<nDim,nDim>() )
                {
                    auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                    stab_test.expression().template set<3>( -coeffNatureStabilization*inner(coeff_c_expr,hess(u)) );
                }
                else
                {
                    auto const& coeff_c_expr = expr( coeff_c.expr(), se );
                    stab_test.expression().template set<2>( -coeffNatureStabilization*coeff_c_expr*laplacian(u) );
                }
            }
        }
    }

    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*(gradt(u)*coeff_beta_expr) )>;
    using expr_reaction_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_coeff_reaction_type{}*idt(u) )>;
    using expr_diffusion_scalar_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*expr_coeff_diffusion_scalar_type{}*laplaciant(u) )>;
    using expr_diffusion_matrix_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_coeff_diffusion_matrix_type{},hesst(u)) )>;

    auto residual_lhs = exprOptionalConcat<expr_convection_residual_lhs_type,expr_reaction_residual_lhs_type,expr_diffusion_scalar_residual_lhs_type,expr_diffusion_matrix_residual_lhs_type>();
    residual_lhs.expression().template set<0>( timeSteppingScaling*(gradt(u)*coeff_beta_expr) );

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
    if ( hasDiffusionTerm )
    {
        // convection - diffusion (-reaction)
        auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
        if ( coeff_c.template hasExpr<nDim,nDim>() )
        {
            auto const& coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
            residual_lhs.expression().template set<3>( -timeSteppingScaling*inner(coeff_c_expr,hesst(u)) );
            if ( hasReactionTerm )
            {
                auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
                auto const& coeff_a_expr = expr( coeff_a.expr(), se );
                residual_lhs.expression().template set<1>( -timeSteppingScaling*coeff_a_expr*idt(u) );
                CHECK( false ) << "TODO";
            }
            else
            {
                auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),coeff_beta_expr, coeff_c_expr, true, true );
                tauFieldPtr->on(_range=range,_expr=tauExpr);
            }
        }
        else
        {
            auto const& coeff_c_expr = expr( coeff_c.expr(), se );
            residual_lhs.expression().template set<2>( -timeSteppingScaling*coeff_c_expr*laplaciant(u) );
            if ( hasReactionTerm )
            {
                auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
                auto const& coeff_a_expr = expr( coeff_a.expr(), se );
                residual_lhs.expression().template set<1>( -timeSteppingScaling*coeff_a_expr*idt(u) );
                CHECK( false ) << "TODO";
            }
            else
            {
                auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),coeff_beta_expr, coeff_c_expr, true, true );
                tauFieldPtr->on(_range=range,_expr=tauExpr);
            }
        }
    }
    else if ( hasReactionTerm )
    {
        // convection - reaction
        auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
        auto const& coeff_a_expr = expr( coeff_a.expr(), se );
        residual_lhs.expression().template set<1>( -timeSteppingScaling*coeff_a_expr*idt(u) );
        CHECK( false ) << "TODO";
    }
    else // convection pure
    {
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),coeff_beta_expr, cst(0.), true, false );
        tauFieldPtr->on(_range=range,_expr=tauExpr);
    }

    auto tau = idv(tauFieldPtr);

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    bilinearForm +=
        integrate( _range=range,
                   _expr=tau*inner(residual_lhs,stab_test),
                   _geomap=this->geomap() );

}

} // namespace FeelModels
} // namespace Feel

#endif
