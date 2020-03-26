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

    //bool build_ReactionTerm = buildNonCstPart;
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

        if ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
        {
            auto const& coeff_beta = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
            auto coeff_beta_expr = expr( coeff_beta.template expr<nDim,1>(), se );
            bool build_ConvectionTerm = coeff_beta_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_ConvectionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*inner( gradt(u)*coeff_beta_expr, id(v) ),
                               //_expr= timeSteppingScaling*(gradt(u)*coeff_beta_expr)*id(v),
                               _geomap=this->geomap() );
            }
        }

        // Diffusion
        if ( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
        {
            auto const& coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
            if ( coeff_c.template hasExpr<nDim,nDim>() )
            {
                auto coeff_c_expr = expr( coeff_c.template expr<nDim,nDim>(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    if constexpr ( unknown_is_vectorial )
                        {
                            bilinearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*inner(coeff_c_expr*gradt(u),grad(v)),
                                           _geomap=this->geomap() );
                        }
                    else
                    {
                            bilinearForm +=
                                integrate( _range=range,
                                           _expr= timeSteppingScaling*grad(v)*(coeff_c_expr*trans(gradt(u))),
                                           _geomap=this->geomap() );
                    }
                }
            }
            else
            {
                auto coeff_c_expr = expr( coeff_c.expr(), se );
                bool build_DiffusionTerm = coeff_c_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
                if ( build_DiffusionTerm )
                {
                    bilinearForm +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*coeff_c_expr*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }

        // Reaction
        if ( this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() ) )
        {
            auto const& coeff_a = this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() );
            auto coeff_a_expr = expr( coeff_a.expr(), se );
            bool build_ReactionTerm = coeff_a_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( build_ReactionTerm )
            {
                bilinearForm +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*coeff_a_expr*inner(idt(u), id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Transient terms
        if ( !this->isStationary() && this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
        {
            auto const& coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
            auto coeff_d_expr = expr( coeff_d.expr(), se );

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
            auto coeff_f_expr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                               [&coeff_f,&se] { return expr( coeff_f.expr(), se ); },
                                               [&coeff_f,&se] { return expr( coeff_f.template expr<nDim,1>(), se ); } );
            //auto const& coeff_f_expr = expr( coeff_f.expr(), se );
            bool buildSourceTerm = coeff_f_expr.expression().isNumericExpression()? buildCstPart : buildNonCstPart;
            if ( buildSourceTerm )
            {
                linearForm +=
                    integrate( _range=range,
                               _expr= inner(coeff_f_expr,id(v)),
                               _geomap=this->geomap() );
            }
        }

        // Stab
        if ( this->M_applyStabilization && buildNonCstPart && this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
        {
            this->updateLinearPDEStabilizationGLS( data, mctx, matName, range );
        }

    } // for each material


    // update weak bc
    if ( buildNonCstPart )
    {
        // k \nabla u  n = g
        for( auto const& d : this->M_bcNeumann )
        {
            auto theExpr = hana::eval_if( hana::bool_c<unknown_is_scalar>,
                                          [&d,&se] { return expression( d,se ); },
                                          [&d,&se] { return expression( d,se )*N(); } );
            linearForm +=
                integrate( _range=markedfaces(this->mesh(),M_bcNeumannMarkerManagement.markerNeumannBC(MarkerManagementNeumannBC::NeumannBCShape::SCALAR,name(d)) ),
                           _expr= timeSteppingScaling*inner(theExpr,id(v)),
                           _geomap=this->geomap() );
        }

        // k \nabla u  n + g u = h  -> - k \nabla u  n = gu - h
        for( auto const& d : this->M_bcRobin )
        {
            if constexpr( unknown_is_scalar )
            {
                auto theExpr1 = expression1( d,se );
                bilinearForm +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr1*idt(v)*id(v),
                               _geomap=this->geomap() );
                auto theExpr2 = expression2( d,se );
                linearForm +=
                    integrate( _range=markedfaces(mesh,M_bcRobinMarkerManagement.markerRobinBC( name(d) ) ),
                               _expr= timeSteppingScaling*theExpr2*id(v),
                               _geomap=this->geomap() );
            }
            else
                CHECK( false ) << "robin bc with vectorial unknown can not be implemented for now";
        }
    }


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
