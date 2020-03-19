/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_ASSEMBLY_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateLinearPDE( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("CoefficientFormPDE","updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();

    bool build_ReactionTerm = buildNonCstPart;
    // bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
    // bool BuildNonCstPart_Form1TransientTerm = buildNonCstPart;
    // if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    //     BuildNonCstPart_Form2TransientTerm = buildCstPart;

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
#if 0
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
#endif
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    auto const& se = mctx.symbolsExpr();

    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = this->fieldUnknown();
    auto const& v = this->fieldUnknown();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto linearForm = form1( _test=Xh, _vector=F,
                             _rowstart=this->rowStartInVector() );

    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physic() ) )
    {
        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( matName );
        // Diffusion
        if ( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
        {
            auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
            auto const& coeff_c_expr = expr( coeff_c.expr(), se );
            bool build_DiffusionTerm = coeff_c_expr.expression().isConstant()? buildCstPart : buildNonCstPart;
            if ( /*this->hasDiffusion() &&*/ build_DiffusionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= coeff_c_expr*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }

        if ( this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() ) )
        {
            auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
            auto const& coeff_a_expr = expr( coeff_a.expr(), se );
            if ( build_ReactionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr=coeff_a_expr*inner(idt(u), id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Transient terms
        if ( !this->isStationary() && this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
        {
            auto const& coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
            auto const& coeff_d_expr = expr( coeff_d.expr(), se );

            bool build_Form2TransientTerm = buildNonCstPart;
            bool build_Form1TransientTerm = buildNonCstPart;
            if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
            {
                build_Form2TransientTerm = buildCstPart;
            }

            if (build_Form2TransientTerm)
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr=M_bdfUnknown->polyDerivCoefficient(0)*inner(idt(u),id(v)),
                               _geomap=this->geomap() );
            }
            if (build_Form1TransientTerm)
            {
                linearForm +=
                    integrate( _range=range,
                               _expr=inner(idv(M_bdfUnknown->polyDeriv()),id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Source
        if ( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) )
        {
            auto const& coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );
            auto const& coeff_f_expr = expr( coeff_f.expr(), se );
            bool buildSourceTerm = coeff_f_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( buildSourceTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= coeff_f_expr*id(v),
                               _geomap=this->geomap() );
            }
        }

    } // for each material


}
template< typename ConvexType, typename BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateLinearPDEDofElimination( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( !M_bcDirichletMarkerManagement.hasMarkerDirichletBCelimination() ) return;

    this->log("CoefficientFormPDE","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& se = mctx.symbolsExpr();
    auto mesh = this->mesh();
    auto Xh = this->spaceUnknown();
    auto const& u = this->fieldUnknown();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    for( auto const& d : M_bcDirichlet )
    {
        auto theExpr = expression(d,se);
        bilinearForm +=
            on( _range=markedfaces(mesh, M_bcDirichletMarkerManagement.markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=theExpr );
    }

    this->log("CoefficientFormPDE","updateLinearPDEDofElimination","finish" );
}

} // namespace Feel
} // namespace FeelModels

#endif
