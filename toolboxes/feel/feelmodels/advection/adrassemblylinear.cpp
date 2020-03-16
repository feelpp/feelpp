/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/advection/advection.hpp>

namespace Feel {
namespace FeelModels {

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data, this->modelContext() );
}

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    if ( this->M_bcDirichlet.empty() && this->M_bcInflowMarkers.empty() ) return;

    this->log("AdvDiffReac","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto mesh = this->mesh();
    auto Xh = this->functionSpace();
    auto const& u = this->fieldSolution();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto se = this->symbolsExpr();

    // Dirichlet bc
    for( auto const& d : this->M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=expression(d,se) );
    }

#if 0 // TODO VINCENT
    // Inflow bc
    if( !this->isStationary() )
    {
        for( auto const& bcMarker: this->M_bcInflowMarkers )
        {
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(mesh, bcMarker),
                        _element=u,
                        _rhs=F,
                        _expr=(
                            // Transient part
                            idv(this->timeStepBDF()->polyDeriv())
                            // Advection part
                            - gradv(u)*idv(this->fieldAdvectionVelocity())
                            // Diffusion part
                            + (this->hasDiffusion()) * idv(this->diffusionReactionModel()->fieldDiffusionCoeff())*laplacianv(u)
                            // Reaction part
                            - (this->hasReaction()) * idv(this->diffusionReactionModel()->fieldReactionCoeff())*idv(u)
                            // Source part
                            + (this->hasSourceAdded() || this->hasSourceTerm()) * idv(*this->M_fieldSource)
                            )/(this->timeStepBDF()->polyDerivCoefficient(0))
                  );
        }
    }
#endif
    this->log("AdvDiffReac","updateLinearPDEDofElimination","finish" );
}

} // namespace FeelModels
} // namespace Feel
