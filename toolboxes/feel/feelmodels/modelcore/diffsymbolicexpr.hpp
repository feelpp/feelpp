/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_MODELCORE_DIFFSYMBOLICEXPR_H
#define FEELPP_TOOLBOXES_MODELCORE_DIFFSYMBOLICEXPR_H 1

#include <feel/feelmodels/modelcore/trialssymbolsexpr.hpp>
#include <feel/feelmodels/modelvf/exproptionalconcat.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename SymbolicExprType>
struct TransformDiffSymbolicExpr
{
    TransformDiffSymbolicExpr( SymbolicExprType const& expr ) : M_expr( expr ) {}

    template <typename T>
    constexpr auto operator()(T const& t) const
        {
            auto firstTse = *t.begin();
#if 0
            auto diffExpr = diff( M_expr, firstTse.symbol(),1);
#else
            auto diffExpr = M_expr.template diff<1>( firstTse.symbol()/*,world, dirLibExpr*/ );
#endif
            auto trialExpr = firstTse.expr();
            return trialExpr*diffExpr;
        }

    template <typename... TheType>
    static constexpr auto
    toExpr( hana::tuple<TheType...> const& t )
        {
            return exprOptionalConcat<TheType...>();
        }

private :
    SymbolicExprType const& M_expr;
};


template<typename SymbolicExprType, typename SpaceType, typename TupleTrialSymbolExprType>
auto diffSymbolicExpr( SymbolicExprType const& theExpr, TrialSymbolsExpr<SpaceType, TupleTrialSymbolExprType> const tse,
                       std::shared_ptr<SpaceType> const& trialSpace, size_type trialBlockSpaceIndex,
                       WorldComm const& world = Environment::worldComm(), std::string const& dirLibExpr = Environment::exprRepository() )
{
    using the_expr_type = std::decay_t<decltype( TransformDiffSymbolicExpr<SymbolicExprType>::toExpr( hana::transform( tse.tuple(), TransformDiffSymbolicExpr<SymbolicExprType>(theExpr) ) ) )>;
    the_expr_type resExpr;

    hana::for_each( tse.tuple(), [&theExpr,&trialSpace,&trialBlockSpaceIndex,&world,&dirLibExpr,&resExpr]( auto const& e )
                    {
                        for ( auto const& e2 : e )
                        {
                            if ( !e2.isDefinedOn( trialSpace, trialBlockSpaceIndex ) )
                                continue;

                            if ( !theExpr.hasSymbolDependency( e2.symbol() ) )
                                continue;
#if 0
                            auto diffExpr = diff( theExpr, e2.symbol(),1,"",world, dirLibExpr);
#else
                            auto diffExpr = theExpr.template diff<1>( e2.symbol(),world, dirLibExpr );
#endif
                            // std::cout << "diffExprBIS = " << str( diffExprBIS.expression() ) << std::endl;
                            // std::cout << "diffExpr    = " << str( diffExpr.expression() ) << std::endl;
                            // std::cout << "diffExprBIS se names : " << diffExprBIS.expression().symbolsExpression().names() << std::endl;
                            // std::cout << "diffExpr    se names : " << diffExpr.expression().symbolsExpression().names() << std::endl;

                            auto trialExpr = e2.expr();
                            //std::cout << "diffSymbolicExpr add expr" << std::endl;

                            resExpr.expression().add( trialExpr*diffExpr );
                        }
                    });

    return resExpr;
}

} // namespace FeelModels
} // namespace Feel

#endif
