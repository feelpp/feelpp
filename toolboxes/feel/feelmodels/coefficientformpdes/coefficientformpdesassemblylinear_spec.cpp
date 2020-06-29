/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include "../coefficientformpdesconfig.h"
#include "coefficientformpdesbasisspecialisation.h"
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>

namespace Feel {
namespace FeelModels {

template <>
template <>
void
COEFFICIENTFORMPDES_CLASS_TYPE::updateLinearPDE_spec<typename COEFFICIENTFORMPDES_CLASS_TYPE::template FilterBasisUnknown< COEFFICIENTFORMPDES_UNKNOWN_BASIS_SPECIALISATION > >( DataUpdateLinear & data, std::any const& mctxAsAny ) const
{
    using the_model_context_type = std::decay_t<decltype(this->modelContext( vector_ptrtype{} ))>;
    try
    {
        auto const& mctx = *std::any_cast<const the_model_context_type*>(mctxAsAny);
        this->updateLinearPDE<typename COEFFICIENTFORMPDES_CLASS_TYPE::template FilterBasisUnknown< COEFFICIENTFORMPDES_UNKNOWN_BASIS_SPECIALISATION> >( data, mctx );
    }
    catch (const std::bad_any_cast& e)
    {
        CHECK( false ) << "ERROR updateLinearPDE_spec : " <<  e.what() << '\n';
    }
}

template <>
template <>
void
COEFFICIENTFORMPDES_CLASS_TYPE::updateLinearPDEDofElimination_spec<typename COEFFICIENTFORMPDES_CLASS_TYPE::template FilterBasisUnknown< COEFFICIENTFORMPDES_UNKNOWN_BASIS_SPECIALISATION > >( DataUpdateLinear & data, std::any const& mctxAsAny ) const
{
    using the_model_context_type = std::decay_t<decltype(this->modelContext( vector_ptrtype{} ))>;
    auto const& mctx = *std::any_cast<const the_model_context_type*>(mctxAsAny);
    this->updateLinearPDEDofElimination<typename COEFFICIENTFORMPDES_CLASS_TYPE::template FilterBasisUnknown< COEFFICIENTFORMPDES_UNKNOWN_BASIS_SPECIALISATION> >( data, mctx );
}

} // namespace FeelModels
} // namespace Feel
