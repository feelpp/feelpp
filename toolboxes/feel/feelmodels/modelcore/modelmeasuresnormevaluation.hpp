/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_CORE_MEASURE_NORM_EVALUATION_HPP
#define FEELPP_TOOLBOXES_CORE_MEASURE_NORM_EVALUATION_HPP 1

#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/normh1.hpp>
#include <feel/feelvf/normsemih1.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename RangeType>
struct NormEvalRange
{
    typedef typename boost::tuples::template element<1, RangeType>::type element_iterator;
    typedef typename boost::unwrap_reference<typename element_iterator::value_type>::type range_elt_type;
    typedef typename boost::remove_reference<range_elt_type>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_face_element_type;
    typedef typename the_face_element_type::super2::template Element<the_face_element_type>::type element_type;
};
template<typename RangeType, typename ExprType>
struct NormEvalExpr
{
    typedef typename NormEvalRange<RangeType>::element_type element_type;

    typedef typename element_type::gm_type gm_type;
    typedef typename gm_type::template Context<ExprType::context|vm::JACOBIAN, element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename ExprType::template tensor<map_gmc_type> eval_expr_type;
    typedef typename eval_expr_type::shape shape;
};

template<typename RangeType, typename ExprType, typename GradExprType, typename SymbolsExpr>
void
measureNormEvaluationH1( RangeType const& range, ExprType const& idExpr, GradExprType const& gradExpr, std::string const& normType,
                         ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, std::map<std::string,double> & res,
                         typename std::enable_if< NormEvalExpr<RangeType,ExprType>::shape::is_scalar>::type* = nullptr )
{
    typedef typename NormEvalExpr<RangeType,ExprType>::shape shape_type;
    typedef typename NormEvalExpr<RangeType,ExprType>::element_type element_type;

    std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
    double normComputed = 0;
    if ( normType == "H1" )
        normComputed = normH1(_range=range,_expr=idExpr, _grad_expr=gradExpr );
    else if ( normType == "SemiH1" )
        normComputed = normSemiH1(_range=range, _grad_expr=gradExpr );
    else if ( normType == "H1-error" )
        normComputed = normH1(_range=range,_expr=idExpr - expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                              _grad_expr= gradExpr -  trans( expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(), symbolsExpr ) ),
                              _quad=_Q<10>(),_quad1=_Q<10>() );
    else if ( normType == "SemiH1-error" )
        normComputed = normSemiH1(_range=range,_grad_expr= gradExpr - trans( expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(),symbolsExpr ) ),
                                  _quad=_Q<10>(),_quad1=_Q<10>() );
    else
        CHECK( false ) << "not a H1 norm type";
    res[normNameOutput] = normComputed;
}

template<typename RangeType, typename ExprType, typename GradExprType, typename SymbolsExpr>
void
measureNormEvaluationH1( RangeType const& range, ExprType const& idExpr, GradExprType const& gradExpr, std::string const& normType,
                         ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, std::map<std::string,double> & res,
                         typename std::enable_if< NormEvalExpr<RangeType,ExprType>::shape::is_vectorial>::type* = nullptr )
{
    typedef typename NormEvalExpr<RangeType,ExprType>::shape shape_type;
    typedef typename NormEvalExpr<RangeType,ExprType>::element_type element_type;

    std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
    double normComputed = 0;
    if ( normType == "H1" )
        normComputed = normH1(_range=range,_expr=idExpr, _grad_expr=gradExpr );
    else if ( normType == "SemiH1" )
        normComputed = normSemiH1(_range=range, _grad_expr=gradExpr );
    else if ( normType == "H1-error" )
        normComputed = normH1(_range=range,_expr=idExpr - expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                              _grad_expr= gradExpr - expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(), symbolsExpr ),
                              _quad=_Q<10>(),_quad1=_Q<10>() );
    else if ( normType == "SemiH1-error" )
        normComputed = normSemiH1(_range=range,_grad_expr= gradExpr - expr( ppNorm.gradSolution().template expr<element_type::nRealDim,shape_type::M>(), symbolsExpr ),
                                  _quad=_Q<10>(),_quad1=_Q<10>() );
    else
        CHECK( false ) << "not a H1 norm type";
    res[normNameOutput] = normComputed;
}

template<typename RangeType, typename ExprType, typename SymbolsExpr>
void
measureNormEvaluationL2( RangeType const& range, ExprType const& idExpr, std::string const& normType,
                         ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, std::map<std::string,double> & res )
{
    typedef typename NormEvalExpr<RangeType,ExprType>::shape shape_type;
    typedef typename NormEvalExpr<RangeType,ExprType>::element_type element_type;

    std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
    double normComputed = 0;
    if ( normType == "L2" )
        normComputed = normL2(_range=range,_expr=idExpr );
    else if ( normType == "L2-error" )
        normComputed = normL2(_range=range,_expr=idExpr - expr( ppNorm.solution().template expr<shape_type::M,shape_type::N>(), symbolsExpr ),
                              _quad=_Q<10>(),_quad1=_Q<10>() );
    else
        CHECK( false ) << "not a L2 norm type";
    res[normNameOutput] = normComputed;
}

template<typename RangeType, typename FieldType, typename SymbolsExpr>
void
measureNormEvaluationField( RangeType const& range, FieldType const& field, std::string const& normType,
                            ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, std::map<std::string,double> & res,
                            typename std::enable_if< FieldType::is_scalar || FieldType::is_vectorial >::type* = nullptr )
{
    if ( normType == "L2" || normType == "L2-error" )
        measureNormEvaluationL2( range, idv(field), normType, ppNorm, symbolsExpr, res );
    else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
        measureNormEvaluationH1( range, idv(field), gradv(field), normType, ppNorm, symbolsExpr, res );
    else
        CHECK( false ) << "invalid norm type : " << normType;
}

template<typename RangeType, typename FieldType, typename SymbolsExpr>
void
measureNormEvaluationField( RangeType const& range, FieldType const& field, std::string const& normType,
                            ModelPostprocessNorm const& ppNorm, SymbolsExpr const& symbolsExpr, std::map<std::string,double> & res,
                            typename std::enable_if< FieldType::is_tensor2 || FieldType::is_tensor2symm >::type* = nullptr )
{
    if ( normType == "L2" || normType == "L2-error" )
        measureNormEvaluationL2( range, idv(field), normType, ppNorm, symbolsExpr, res );
    else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error" )
        CHECK( false ) << "normType " << normType << " is not implemented with tensor field";
    else
        CHECK( false ) << "invalid norm type : " << normType;
}


template<typename RangeType, typename SymbolsExpr, typename... Args >
void
measureNormEvaluation( RangeType const& range,
                       ModelPostprocessNorm const& ppNorm,  std::map<std::string,double> & res,
                       SymbolsExpr const& symbolsExpr, Args... args )
{
    typedef typename NormEvalRange<RangeType>::element_type element_type;
    static const uint16_type nRealDim = element_type::nRealDim;

    if ( ppNorm.hasField() )
    {
        hana::for_each( hana::make_tuple(args...), [&]( auto const& e )
                        {
                            if ( ppNorm.field() == e.first )
                            {
                                auto const& fieldFunc = e.second;
                                for ( std::string const& normType : ppNorm.types() )
                                    measureNormEvaluationField( range, fieldFunc, normType, ppNorm, symbolsExpr, res );
                            }
                        });
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


template<typename MeshType, typename RangeType, typename SymbolsExpr, typename... Args >
void
measureNormEvaluation( boost::shared_ptr<MeshType> const& mesh, RangeType const& defaultRange,
                       ModelPostprocessNorm const& ppNorm,  std::map<std::string,double> & res,
                       SymbolsExpr const& symbolsExpr, Args... args )
{
    auto meshMarkers = ppNorm.markers();
    if ( meshMarkers.empty() )
        measureNormEvaluation( defaultRange,ppNorm,res,symbolsExpr,args... );
    else
    {
        std::string firstMarker = *meshMarkers.begin();
        if ( mesh->hasElementMarker( firstMarker ) )
            measureNormEvaluation(  markedelements( mesh,ppNorm.markers() ),ppNorm,res,symbolsExpr,args... );
        else if ( mesh->hasFaceMarker( firstMarker ) )
            measureNormEvaluation(  markedfaces( mesh,ppNorm.markers() ),ppNorm,res,symbolsExpr,args... );
        else if ( mesh->hasEdgeMarker( firstMarker ) || mesh->hasPointMarker( firstMarker ) )
            CHECK( false ) << "not implemented for edges/points";
        else if ( !mesh->hasMarker( firstMarker ) )
            CHECK( false ) << "marker " << firstMarker << " not present in mesh";
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
