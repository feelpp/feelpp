/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDES_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDES_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename... BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDEs<ConvexType,BasisUnknownType...>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctx]( auto const& e )
                    {
                        for ( auto const& cfpdeBase : M_coefficientFormPDEs )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                continue;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );

                            cfpde->updateLinearPDE( data, mctx );
                        }
                    });
}

template< typename ConvexType, typename... BasisUnknownType>
template <typename ModelContextType>
void
CoefficientFormPDEs<ConvexType,BasisUnknownType...>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    hana::for_each( tuple_type_unknown_basis, [this,&data,&mctx]( auto const& e )
                    {
                        for ( auto const& cfpdeBase : M_coefficientFormPDEs )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                continue;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );

                            cfpde->updateLinearPDEDofElimination( data, mctx );
                        }
                    });
}

} // namespace Feel
} // namespace FeelModels

#endif
