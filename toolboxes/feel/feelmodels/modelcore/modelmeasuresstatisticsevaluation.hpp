/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_CORE_MEASURE_STATISTICS_EVALUATION_HPP
#define FEELPP_TOOLBOXES_CORE_MEASURE_STATISTICS_EVALUATION_HPP 1

#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/mean.hpp>
#include <feel/feelmodels/modelcore/traits.hpp>
#include <feel/feelcore/tuple_utils.hpp>
#include <feel/feelmodels/modelvf/evalonentities.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename RangeType, typename ExprType >
void
measureStatisticsEvaluationMinMax( RangeType const& range, ExprType const& expr,
                                   ModelPostprocessStatistics const& ppStat,  ModelMeasuresStorage & res,
                                   std::string ppStatType, bool useQuadOrder = true )
{
    //uint16_type quadOrder = (useQuadOrder)? ppStat.quadOrder() : quad_order_from_expression;
    //uint16_type quad1Order = (useQuadOrder)? ppStat.quad1Order() : quad_order_from_expression;
    uint16_type quadOrder = ppStat.quadOrder();
    auto p = minmax(_range=range,_expr=expr,_pset=_Q<>(quadOrder) );
    double pmin = p.min();
    double pmax = p.max();
    auto p_pmin = p.argmin();
    auto p_pmax = p.argmax();

    if ( ppStatType == "min" || ppStatType == "min-max" )
    {
        std::string statNameOutput = (boost::format("Statistics_%1%_min")%ppStat.name() ).str();
        res.setValue(statNameOutput, pmin);
    }
    if ( ppStatType == "max" || ppStatType == "min-max" )
    {
        std::string statNameOutput = (boost::format("Statistics_%1%_max")%ppStat.name() ).str();
        res.setValue(statNameOutput, pmax);
    }
}

template<typename RangeType, typename ExprType >
void
measureStatisticsEvaluationMean( RangeType const& range, ExprType const& expr,
                                 ModelPostprocessStatistics const& ppStat,  ModelMeasuresStorage & res,
                                 bool useQuadOrder = true )
{
    uint16_type quadOrder = (useQuadOrder)? ppStat.quadOrder() : quad_order_from_expression;
    uint16_type quad1Order = (useQuadOrder)? ppStat.quad1Order() : quad_order_from_expression;

    auto statComputed = mean(_range=range,_expr=expr, _quad=quadOrder,_quad1=quad1Order );

    res.setValue( fmt::format("Statistics_{}_mean", ppStat.name()), statComputed);
#if 0
    for ( int i=0;i<statComputed.rows() ;++i )
        for ( int j=0;j<statComputed.cols() ;++j )
        {
            double val = statComputed(i,j);
            std::string statNameOutput;
            if ( statComputed.rows() == 1 && statComputed.cols() == 1 )
                statNameOutput = (boost::format("Statistics_%1%_mean")%ppStat.name() ).str();
            else if ( statComputed.rows() == 1 )
                statNameOutput = (boost::format("Statistics_%1%_mean_%2%")%ppStat.name()%j ).str();
            else if ( statComputed.cols() == 1 )
                statNameOutput = (boost::format("Statistics_%1%_mean_%2%")%ppStat.name()%i ).str();
            else
                statNameOutput = (boost::format("Statistics_%1%_mean_%2%_%3%")%ppStat.name()%i%j ).str();
            res.setValue(statNameOutput, val);
        }
#endif
}

template<typename RangeType, typename ExprType >
void
measureStatisticsEvaluationIntegrate( RangeType const& range, ExprType const& expr,
                                      ModelPostprocessStatistics const& ppStat, ModelMeasuresStorage & res,
                                      bool useQuadOrder = true )
{
    uint16_type quadOrder = (useQuadOrder)? ppStat.quadOrder() : quad_order_from_expression;
    uint16_type quad1Order = (useQuadOrder)? ppStat.quad1Order() : quad_order_from_expression;

    auto exprUsed = evalOnEntities( range,expr,ppStat.requiresMarkersConnection(),ppStat.internalFacesEvalutationType() );
    auto statComputed = integrate(_range=range,_expr=exprUsed/*expr*/, _quad=quadOrder,_quad1=quad1Order ).evaluate();

    res.setValue( fmt::format("Statistics_{}_integrate",ppStat.name()), statComputed);
#if 0
    for ( int i=0;i<statComputed.rows() ;++i )
        for ( int j=0;j<statComputed.cols() ;++j )
        {
            double val = statComputed(i,j);
            std::string statNameOutput;
             if ( statComputed.rows() == 1 && statComputed.cols() == 1 )
                 statNameOutput = (boost::format("Statistics_%1%_integrate")%ppStat.name() ).str();
             else if ( statComputed.rows() == 1 )
                statNameOutput = (boost::format("Statistics_%1%_integrate_%2%")%ppStat.name()%j ).str();
            else if ( statComputed.cols() == 1 )
                statNameOutput = (boost::format("Statistics_%1%_integrate_%2%")%ppStat.name()%i ).str();
            else
                statNameOutput = (boost::format("Statistics_%1%_integrate_%2%_%3%")%ppStat.name()%i%j ).str();
             res.setValue( statNameOutput, val );
        }
#endif
}
template<typename RangeType, typename SymbolsExpr, typename... FieldTupleType >
void
measureStatisticsEvaluation( RangeType const& range,
                             ModelPostprocessStatistics const& ppStat, ModelMeasuresStorage & res,
                             SymbolsExpr const& symbolsExpr, FieldTupleType const& ... fieldTuple )
{
    std::set<std::string> ppStatType;
    bool hasMin = false, hasMax = false;
    for ( std::string const& statType : ppStat.types() )
    {
        if ( statType == "min" )
            hasMin = true;
        else if ( statType == "max" )
            hasMax = true;
        else
            ppStatType.insert( statType );
    }
    if ( hasMin && hasMax )
        ppStatType.insert( "min-max");
    else if ( hasMin )
        ppStatType.insert( "min");
    else if ( hasMax )
        ppStatType.insert( "max");

    if ( ppStat.hasField() )
    {
        ( Feel::for_each( fieldTuple.tuple(), [&]( auto const& e )
                        {
                            if constexpr ( is_iterable_v<decltype(e)> )
                                {
                                    for ( auto const& mfield : e )
                                    {
                                        std::string fieldName = mfield.nameWithPrefix();
                                        auto const& fieldFunc = mfield.field();
                                        if constexpr ( is_shared_ptr<decltype(fieldFunc)>::value )
                                            {
                                                if ( !fieldFunc )
                                                    continue;
                                            }
                                        if ( ppStat.field() == fieldName )
                                        {
                                            mfield.applyUpdateFunction();
                                            for ( std::string const& statType : ppStatType )
                                            {
                                                if ( statType == "min" || statType == "max" || statType == "min-max"  )
                                                    measureStatisticsEvaluationMinMax( range, idv(fieldFunc), ppStat, res, statType, false );
                                                else if ( statType == "mean" )
                                                    measureStatisticsEvaluationMean( range, idv(fieldFunc), ppStat, res, false );
                                                else if ( statType == "integrate" )
                                                    measureStatisticsEvaluationIntegrate( range, idv(fieldFunc), ppStat, res, false );
                                            }
                                        }
                                    }
                                }
                            else
                            {
#if 0
                                if ( ppStat.field() == e.first )
                                {
                                    auto const& fieldFunc = e.second;
                                    if constexpr ( is_shared_ptr<decltype(fieldFunc)>::value )
                                        {
                                            if ( !fieldFunc )
                                                return;
                                        }
                                    for ( std::string const& statType : ppStatType )
                                    {
                                        if ( statType == "min" || statType == "max" || statType == "min-max"  )
                                            measureStatisticsEvaluationMinMax( range, idv(fieldFunc), ppStat, res, statType, false );
                                        else if ( statType == "mean" )
                                            measureStatisticsEvaluationMean( range, idv(fieldFunc), ppStat, res, false );
                                        else if ( statType == "integrate" )
                                            measureStatisticsEvaluationIntegrate( range, idv(fieldFunc), ppStat, res, false );
                                    }
                                }
#endif
                            }
                        }), ... );
    }
    else if ( ppStat.hasExpr() )
    {
        auto exprShape = hana::make_tuple( hana::make_tuple(hana::int_c<1>,hana::int_c<1>),
                                           hana::make_tuple(hana::int_c<2>,hana::int_c<1>),
                                           hana::make_tuple(hana::int_c<3>,hana::int_c<1>),
                                           hana::make_tuple(hana::int_c<2>,hana::int_c<2>),
                                           hana::make_tuple(hana::int_c<3>,hana::int_c<3>)
                                           );
        auto const& exprGeneric = ppStat.expr();
        for ( std::string const& statType : ppStatType )
        {
            hana::for_each( exprShape, [&]( auto const& e )
                            {
                                constexpr int i = std::decay_t<decltype(hana::at_c<0>(e))>::value;
                                constexpr int j = std::decay_t<decltype(hana::at_c<1>(e))>::value;
                                if ( exprGeneric.template hasExpr<i,j>() )
                                {
                                    auto statExpr = expr( exprGeneric.expr<i,j>(), symbolsExpr );
                                    if ( statType == "min" || statType == "max" || statType == "min-max"  )
                                        measureStatisticsEvaluationMinMax( range, statExpr, ppStat, res, statType );
                                    else if ( statType == "mean" )
                                        measureStatisticsEvaluationMean( range, statExpr, ppStat, res );
                                    else if ( statType == "integrate" )
                                        measureStatisticsEvaluationIntegrate( range, statExpr, ppStat, res );
                                }
                            });
        }
    }
}

template<typename MeshType, typename RangeType, typename SymbolsExpr, typename... FieldTupleType >
void
measureStatisticsEvaluation( std::shared_ptr<MeshType> const& mesh, RangeType const& defaultRange,
                             ModelPostprocessStatistics const& ppStat, ModelMeasuresStorage & res,
                             SymbolsExpr const& symbolsExpr, FieldTupleType const& ... fieldTuple )
{
    auto meshMarkers = ppStat.markers();
    if ( meshMarkers.empty() )
        measureStatisticsEvaluation( defaultRange,ppStat,res,symbolsExpr,fieldTuple... );
    else
    {
        std::string firstMarker = *meshMarkers.begin();
        if ( mesh->hasElementMarker( firstMarker ) )
            measureStatisticsEvaluation(  markedelements( mesh,ppStat.markers() ),ppStat,res,symbolsExpr,fieldTuple... );
        else if ( mesh->hasFaceMarker( firstMarker ) )
            measureStatisticsEvaluation(  markedfaces( mesh,ppStat.markers() ),ppStat,res,symbolsExpr,fieldTuple... );
        else if ( mesh->hasEdgeMarker( firstMarker ) || mesh->hasPointMarker( firstMarker ) )
            CHECK( false ) << "not implemented for edges/points";
        else if ( !mesh->hasMarker( firstMarker ) )
            CHECK( false ) << "marker " << firstMarker << " not present in mesh";
    }
}


} // namespace FeelModels
} // namespace Feel

#endif
