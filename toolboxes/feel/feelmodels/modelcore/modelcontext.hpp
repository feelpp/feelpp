/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELCONTEXT_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELCONTEXT_H 1

#include <feel/feelmodels/modelcore/modelfields.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename ModelFieldsType, typename SymbolsExprType, typename TrialSymbolsExprType = trials_symbols_expr_empty_t>
struct ModelContext
{
    using model_fields_type = ModelFieldsType;
    using symbol_expr_type = SymbolsExprType;
    using trial_symbol_expr_type = TrialSymbolsExprType;

    ModelContext( model_fields_type const& mfields, symbol_expr_type const& se )
        :
        M_modelFields( mfields ),
        M_symbolsExpr( se )
        {}
    ModelContext( model_fields_type && mfields, symbol_expr_type && se )
        :
        M_modelFields( mfields ),
        M_symbolsExpr( se )
        {}
    ModelContext( model_fields_type && mfields, symbol_expr_type && se, trial_symbol_expr_type && tse )
        :
        M_modelFields( mfields ),
        M_symbolsExpr( se ),
        M_trialSymbolsExpr( tse )
        {}
    ModelContext( ModelContext const& ) = default;
    ModelContext( ModelContext && ) = default;


    template <typename TagType>
    auto const&
    field( TagType const& thetag,std::string const& name ) const
        {
            return M_modelFields.field( thetag, name );
        }

    symbol_expr_type const& symbolsExpr() const { return M_symbolsExpr; }

    trial_symbol_expr_type const& trialSymbolsExpr() const { return M_trialSymbolsExpr; }

private :
    model_fields_type M_modelFields;
    symbol_expr_type M_symbolsExpr;
    trial_symbol_expr_type M_trialSymbolsExpr;
};



template <typename ModelFieldsType, typename SymbolsExprType>
auto modelContext( ModelFieldsType const& mfields, SymbolsExprType const& se )
{
    return ModelContext<ModelFieldsType,SymbolsExprType>( mfields,se );
}

template <typename ModelFieldsType, typename SymbolsExprType>
auto modelContext( ModelFieldsType && mfields, SymbolsExprType && se )
{
    return ModelContext<ModelFieldsType,SymbolsExprType>( std::forward<ModelFieldsType>(mfields), std::forward<SymbolsExprType>(se) );
}

template <typename ModelFieldsType, typename SymbolsExprType, typename TrialSymbolsExprType>
auto modelContext( ModelFieldsType && mfields, SymbolsExprType && se, TrialSymbolsExprType && tse )
{
    return ModelContext<ModelFieldsType,SymbolsExprType,TrialSymbolsExprType>( std::forward<ModelFieldsType>(mfields), std::forward<SymbolsExprType>(se), std::forward<TrialSymbolsExprType>(tse) );
}

} // namespace FeelModels
} // namespace Feel

#endif
