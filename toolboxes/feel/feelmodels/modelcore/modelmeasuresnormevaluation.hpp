/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_CORE_MEASURE_NORM_EVALUATION_HPP
#define FEELPP_TOOLBOXES_CORE_MEASURE_NORM_EVALUATION_HPP 1

#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/normh1.hpp>
#include <feel/feelvf/normsemih1.hpp>
#include <feel/feelmodels/modelcore/traits.hpp>
#include <feel/feelcore/tuple_utils.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename RangeType, typename ExprType, typename GradExprType, typename SymbolsExpr>
void
measureNormEvaluationH1( RangeType const& range, ExprType const& idExpr, GradExprType const& gradExpr, std::string const& normType,
                         ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, ModelMeasuresStorage & res, bool useQuadOrder = true,
                         typename std::enable_if< ExprTraits<RangeType,ExprType>::shape::is_scalar>::type* = nullptr )
{
    typedef typename ExprTraits<RangeType,ExprType>::shape shape_type;
    typedef typename ExprTraits<RangeType,ExprType>::element_type element_type;
    uint16_type quadOrder = (useQuadOrder)? ppNorm.quadOrder() : quad_order_from_expression;
    uint16_type quadOrderError = ppNorm.quadOrder();
    uint16_type quad1Order = (useQuadOrder)? ppNorm.quad1Order() : quad_order_from_expression;
    uint16_type quad1OrderError = ppNorm.quad1Order();

    std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
    double normComputed = 0;
    if ( normType == "H1" )
        normComputed = normH1(_range=range,_expr=idExpr, _grad_expr=gradExpr, _quad=quadOrder,_quad1=quad1Order );
    else if ( normType == "SemiH1" )
        normComputed = normSemiH1(_range=range, _grad_expr=gradExpr, _quad=quadOrder,_quad1=quad1Order );
    else if ( normType == "H1-error" )
        normComputed = normH1(_range=range,_expr=idExpr - expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                              _grad_expr= gradExpr -  trans( expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(), symbolsExpr ) ),
                              _quad=quadOrderError,_quad1=quad1OrderError );
    else if ( normType == "SemiH1-error" )
        normComputed = normSemiH1(_range=range,_grad_expr= gradExpr - trans( expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(),symbolsExpr ) ),
                                  _quad=quadOrderError,_quad1=quad1OrderError );
    else
        CHECK( false ) << "not a H1 norm type";
    res.setValue( normNameOutput, normComputed );
}

template<typename RangeType, typename ExprType, typename GradExprType, typename SymbolsExpr>
void
measureNormEvaluationH1( RangeType const& range, ExprType const& idExpr, GradExprType const& gradExpr, std::string const& normType,
                         ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, ModelMeasuresStorage & res, bool useQuadOrder = true,
                         typename std::enable_if< ExprTraits<RangeType,ExprType>::shape::is_vectorial>::type* = nullptr )
{
    typedef typename ExprTraits<RangeType,ExprType>::shape shape_type;
    typedef typename ExprTraits<RangeType,ExprType>::element_type element_type;
    uint16_type quadOrder = (useQuadOrder)? ppNorm.quadOrder() : quad_order_from_expression;
    uint16_type quadOrderError = ppNorm.quadOrder();
    uint16_type quad1Order = (useQuadOrder)? ppNorm.quad1Order() : quad_order_from_expression;
    uint16_type quad1OrderError = ppNorm.quad1Order();

    std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
    double normComputed = 0;
    if ( normType == "H1" )
        normComputed = normH1(_range=range,_expr=idExpr, _grad_expr=gradExpr, _quad=quadOrder,_quad1=quad1Order );
    else if ( normType == "SemiH1" )
        normComputed = normSemiH1(_range=range, _grad_expr=gradExpr, _quad=quadOrder,_quad1=quad1Order );
    else if ( normType == "H1-error" )
        normComputed = normH1(_range=range,_expr=idExpr - expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                              _grad_expr= gradExpr - expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(), symbolsExpr ),
                              _quad=quadOrderError,_quad1=quad1OrderError );
    else if ( normType == "SemiH1-error" )
        normComputed = normSemiH1(_range=range,_grad_expr= gradExpr - expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(), symbolsExpr ),
                                  _quad=quadOrderError,_quad1=quad1OrderError );
    else
        CHECK( false ) << "not a H1 norm type";
    res.setValue( normNameOutput, normComputed );
}

template<typename RangeType, typename ExprType, typename SymbolsExpr>
void
measureNormEvaluationL2( RangeType const& range, ExprType const& idExpr, std::string const& normType,
                         ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, ModelMeasuresStorage & res, bool useQuadOrder = true )
{
    typedef typename ExprTraits<RangeType,ExprType>::shape shape_type;
    typedef typename ExprTraits<RangeType,ExprType>::element_type element_type;
    uint16_type quadOrder = (useQuadOrder)? ppNorm.quadOrder() : quad_order_from_expression;
    uint16_type quadOrderError = ppNorm.quadOrder();
    uint16_type quad1Order = (useQuadOrder)? ppNorm.quad1Order() : quad_order_from_expression;
    uint16_type quad1OrderError = ppNorm.quad1Order();

    std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
    double normComputed = 0;
    if ( normType == "L2" )
        normComputed = normL2(_range=range,_expr=idExpr,_quad=quadOrder,_quad1=quad1Order );
    else if ( normType == "L2-error" || normType == "L2-relative-error" )
    {
        normComputed = normL2(_range=range,_expr=idExpr - expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                              _quad=quadOrderError,_quad1=quad1OrderError );
        if ( normType == "L2-relative-error" )
        {
            double normSolution = normL2(_range=range,_expr=expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                                         _quad=quadOrderError,_quad1=quad1OrderError );
            if( normSolution > 1e-10 )
                normComputed /= normSolution;
        }
    }
    else
        CHECK( false ) << "not a L2 norm type";
    res.setValue( normNameOutput, normComputed );
}

template<typename RangeType, typename FieldType, typename SymbolsExpr>
void
measureNormEvaluationField( RangeType const& range, FieldType const& field, std::string const& normType,
                            ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, ModelMeasuresStorage & res,
                            typename std::enable_if< FieldType::is_scalar || FieldType::is_vectorial >::type* = nullptr )
{
    if ( normType == "L2" || normType == "L2-error" || normType == "L2-relative-error" )
        measureNormEvaluationL2( range, idv(field), normType, ppNorm, symbolsExpr, res );
    else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
        measureNormEvaluationH1( range, idv(field), gradv(field), normType, ppNorm, symbolsExpr, res, false );
    else
        CHECK( false ) << "invalid norm type : " << normType;
}

template<typename RangeType, typename FieldType, typename SymbolsExpr>
void
measureNormEvaluationField( RangeType const& range, FieldType const& field, std::string const& normType,
                            ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, ModelMeasuresStorage & res,
                            typename std::enable_if< FieldType::is_tensor2 || FieldType::is_tensor2symm >::type* = nullptr )
{
    if ( normType == "L2" || normType == "L2-error" || normType == "L2-relative-error" )
        measureNormEvaluationL2( range, idv(field), normType, ppNorm, symbolsExpr, res, false );
    else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
        CHECK( false ) << "normType " << normType << " is not implemented with tensor field";
    else
        CHECK( false ) << "invalid norm type : " << normType;
}


template<typename RangeType, typename SymbolsExpr, typename... FieldTupleType >
void
measureNormEvaluation( RangeType const& range,
                       ModelPostprocessNorm const& ppNorm, ModelMeasuresStorage & res,
                       SymbolsExpr const& symbolsExpr, FieldTupleType const& ... fieldTuple )
{
    typedef typename RangeTraits<RangeType>::element_type element_type;
    constexpr uint16_type nRealDim = element_type::nRealDim;

    if ( ppNorm.hasField() )
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
                                        if ( ppNorm.field() == fieldName )
                                        {
                                            mfield.applyUpdateFunction();
                                            for ( std::string const& normType : ppNorm.types() )
                                                measureNormEvaluationField( range, unwrap_ptr(fieldFunc), normType, ppNorm, symbolsExpr, res );
                                        }
                                    }
                                }
                            else
                            {
#if 0
                                if ( ppNorm.field() == e.first )
                                {
                                    auto const& fieldFunc = e.second;
                                    if constexpr ( is_shared_ptr<decltype(fieldFunc)>::value )
                                        {
                                            if ( !fieldFunc )
                                                return;
                                        }
                                    for ( std::string const& normType : ppNorm.types() )
                                        measureNormEvaluationField( range, unwrap_ptr(fieldFunc), normType, ppNorm, symbolsExpr, res );
                                }
#endif
                            }
                        }), ... );
    }
    else if ( ppNorm.hasExpr() )
    {
        auto const& exprGeneric = ppNorm.expr();
        auto const& gradExprGeneric = ppNorm.gradExpr();
        for ( std::string const& normType : ppNorm.types() )
        {
            if ( exprGeneric.hasExprScalar() )
            {
                auto idExpr = expr( exprGeneric.expr<1,1>(), symbolsExpr );
                if ( normType == "L2" || normType == "L2-error" )
                    measureNormEvaluationL2( range, idExpr, normType, ppNorm, symbolsExpr, res );
                else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
                {
                    auto gradExpr = trans( expr( gradExprGeneric.expr<nRealDim,1>(), symbolsExpr ) );
                    measureNormEvaluationH1( range, idExpr, gradExpr, normType, ppNorm, symbolsExpr, res );
                }
                else
                    CHECK( false ) << "invalid norm type : " << normType;
            }
            else if ( exprGeneric.hasExpr<nRealDim,1>() )
            {
                auto idExpr = expr( exprGeneric.expr<nRealDim,1>(), symbolsExpr );
                if ( normType == "L2" || normType == "L2-error" )
                    measureNormEvaluationL2( range, idExpr, normType, ppNorm, symbolsExpr, res );
                else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
                {
                    auto gradExpr = expr( gradExprGeneric.expr<nRealDim,nRealDim>(), symbolsExpr );
                    measureNormEvaluationH1( range, idExpr, gradExpr, normType, ppNorm, symbolsExpr, res );
                }
                else
                    CHECK( false ) << "invalid norm type : " << normType;
            }
            else if ( exprGeneric.hasExpr<nRealDim,nRealDim>() )
            {
                auto idExpr = expr( exprGeneric.expr<nRealDim,nRealDim>(), symbolsExpr );
                if ( normType == "L2" || normType == "L2-error" )
                    measureNormEvaluationL2( range, idExpr, normType, ppNorm, symbolsExpr, res );
                else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
                    CHECK( false ) << "normType " << normType << " is not implemented with tensor field";
                else
                    CHECK( false ) << "invalid norm type : " << normType;
            }
        }
    }
}


template<typename MeshType, typename RangeType, typename SymbolsExpr, typename... FieldTupleType>
void
measureNormEvaluation( std::shared_ptr<MeshType> const& mesh, RangeType const& defaultRange,
                       ModelPostprocessNorm const& ppNorm, ModelMeasuresStorage & res,
                       SymbolsExpr const& symbolsExpr, FieldTupleType const& ... fieldTuple )
{
    auto meshMarkers = ppNorm.markers();
    if ( meshMarkers.empty() )
        measureNormEvaluation( defaultRange,ppNorm,res,symbolsExpr,fieldTuple... );
    else
    {
        std::string firstMarker = *meshMarkers.begin();
        if ( mesh->hasElementMarker( firstMarker ) )
            measureNormEvaluation(  markedelements( mesh,ppNorm.markers() ),ppNorm,res,symbolsExpr,fieldTuple... );
        else if ( mesh->hasFaceMarker( firstMarker ) )
            measureNormEvaluation(  markedfaces( mesh,ppNorm.markers() ),ppNorm,res,symbolsExpr,fieldTuple... );
        else if ( mesh->hasEdgeMarker( firstMarker ) || mesh->hasPointMarker( firstMarker ) )
            CHECK( false ) << "not implemented for edges/points";
        else if ( !mesh->hasMarker( firstMarker ) )
            CHECK( false ) << "marker " << firstMarker << " not present in mesh";
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
