/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>

namespace Feel {
namespace FeelModels {

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    auto mctx = this->modelContext();
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateLinearPDE_spec<CompilerSelectorBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

COEFFICIENTFORMPDES_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDES_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    auto mctx = this->modelContext();
    using the_model_context_type = std::decay_t<decltype(mctx)>;
    auto mctxAsAny = std::make_any<const the_model_context_type*>(&mctx);
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctxAsAny]( auto const& e )
                    {
                        using the_basis_type = typename std::decay_t<decltype(e)>::type;
                        this->updateLinearPDEDofElimination_spec<CompilerSelectorBasisUnknown<the_basis_type>>( data, mctxAsAny );
                    });
}

} // namespace FeelModels
} // namespace Feel
