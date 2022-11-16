/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELCONTEXT_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELCONTEXT_H 1

#include <feel/feelmodels/modelcore/modelfields.hpp>

namespace Feel
{
namespace FeelModels
{

/**
 * @brief Context for Models
 * @ingroup ModelCore
 * 
 * @tparam ModelFieldsType 
 * @tparam SymbolsExprType 
 * @tparam TrialSymbolsExprType 
 */
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
        M_modelFields( std::move(mfields) ),
        M_symbolsExpr( std::move(se) )
        {}
    ModelContext( model_fields_type && mfields, symbol_expr_type && se, trial_symbol_expr_type && tse )
        :
        M_modelFields( std::move(mfields) ),
        M_symbolsExpr( std::move(se) ),
        M_trialSymbolsExpr( std::move(tse) )
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

    void setAdditionalContext( std::string const& key, ModelContext<ModelFieldsType,SymbolsExprType> && mctx )
        {
            auto itFindKey = M_mctxAdditional.find( key );
            if ( itFindKey != M_mctxAdditional.end() )
                M_mctxAdditional.erase( itFindKey );
            M_mctxAdditional.emplace( key, std::forward<ModelContext<ModelFieldsType,SymbolsExprType>>( mctx ) );
        }
    bool hasAdditionalContext( std::string const& key ) const { return M_mctxAdditional.find( key ) != M_mctxAdditional.end(); }

    ModelContext<ModelFieldsType,SymbolsExprType> const& additionalContext( std::string const& key ) const
        {
            auto itFindKey = M_mctxAdditional.find( key );
            CHECK( itFindKey != M_mctxAdditional.end() ) << "no additionalContext";
            return itFindKey->second;
        }
private :
    model_fields_type M_modelFields;
    symbol_expr_type M_symbolsExpr;
    trial_symbol_expr_type M_trialSymbolsExpr;

    std::map<std::string,ModelContext<ModelFieldsType,SymbolsExprType>> M_mctxAdditional;

};

/**
 * @brief create a ModelContext
 * @ingroup ModelCore
 * 
 * @tparam ModelFieldsType 
 * @tparam SymbolsExprType 
 * @param mfields 
 * @param se 
 * @return auto 
 */
template <typename ModelFieldsType, typename SymbolsExprType>
auto modelContext( ModelFieldsType && mfields, SymbolsExprType && se )
{
    return ModelContext<std::decay_t<ModelFieldsType>,std::decay_t<SymbolsExprType>>( std::forward<ModelFieldsType>(mfields), std::forward<SymbolsExprType>(se) );
}

/**
 * @brief create a ModelContext
 * @ingroup ModelCore
 *
 * @tparam ModelFieldsType
 * @tparam SymbolsExprType
 * @tparam TrialSymbolsExprType
 * @param mfields
 * @param se
 * @param tse
 * @return auto
 */
template <typename ModelFieldsType, typename SymbolsExprType, typename TrialSymbolsExprType>
auto modelContext( ModelFieldsType && mfields, SymbolsExprType && se, TrialSymbolsExprType && tse )
{
    return ModelContext<std::decay_t<ModelFieldsType>,std::decay_t<SymbolsExprType>,std::decay_t<TrialSymbolsExprType>>( std::forward<ModelFieldsType>(mfields), std::forward<SymbolsExprType>(se), std::forward<TrialSymbolsExprType>(tse) );
}

} // namespace FeelModels
} // namespace Feel

#endif
