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

template<typename RangeType, typename FieldType>
void
measureNormEvaluation( RangeType const& range, FieldType const& field, ModelPostprocessNorm const& ppNorm,  std::map<std::string,double> & res,
                       typename std::enable_if< FieldType::is_scalar>::type* = nullptr )
{
    for ( std::string const& normType : ppNorm.types() )
    {
        std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
        double normComputed = 0;
        if ( normType == "L2" )
            normComputed = normL2(_range=range,_expr=idv(field) );
        else if ( normType == "H1" )
            normComputed = normH1(_range=range,_expr=idv(field), _grad_expr=gradv(field) );
        else if ( normType == "SemiH1" )
            normComputed = normSemiH1(_range=range, _grad_expr=gradv(field) );
        else if ( normType == "L2-error" )
            normComputed = normL2(_range=range,_expr=idv(field) - ppNorm.solution().template expr<FieldType::nComponents1,FieldType::nComponents2>(),
                                  _quad=_Q<10>(),_quad1=_Q<10>());
        else if ( normType == "H1-error" )
            normComputed = normH1(_range=range,_expr=idv(field) - ppNorm.solution().template expr<FieldType::nComponents1,FieldType::nComponents2>(),
                                  _grad_expr= gradv(field) -  trans( ppNorm.gradSolution().template expr<FieldType::nRealDim,FieldType::nComponents1>() ),
                                  _quad=_Q<10>(),_quad1=_Q<10>() );
        else if ( normType == "SemiH1-error" )
            normComputed = normSemiH1(_range=range,_grad_expr= gradv(field) - trans(ppNorm.gradSolution().template expr<FieldType::nRealDim,FieldType::nComponents1>()),
                                      _quad=_Q<10>(),_quad1=_Q<10>() );
        else
            CHECK( false ) << "error";
        res[normNameOutput] = normComputed;
    }
}

template<typename RangeType, typename FieldType>
void
measureNormEvaluation( RangeType const& range, FieldType const& field, ModelPostprocessNorm const& ppNorm,  std::map<std::string,double> & res,
                       typename std::enable_if< FieldType::is_vectorial>::type* = nullptr )
{
    for ( std::string const& normType : ppNorm.types() )
    {
        std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
        double normComputed = 0;
        if ( normType == "L2" )
            normComputed = normL2(_range=range,_expr=idv(field) );
        else if ( normType == "H1" )
            normComputed = normH1(_range=range,_expr=idv(field), _grad_expr=gradv(field) );
        else if ( normType == "SemiH1" )
            normComputed = normSemiH1(_range=range, _grad_expr=gradv(field) );
        else if ( normType == "L2-error" )
            normComputed = normL2(_range=range,_expr=idv(field) - ppNorm.solution().template expr<FieldType::nComponents1,FieldType::nComponents2>(),
                                  _quad=_Q<10>(),_quad1=_Q<10>() );
        else if ( normType == "H1-error" )
            normComputed = normH1(_range=range,_expr=idv(field) - ppNorm.solution().template expr<FieldType::nComponents1,FieldType::nComponents2>(),
                                  _grad_expr= gradv(field) -  ppNorm.gradSolution().template expr<FieldType::nComponents1,FieldType::nRealDim>(),
                                  _quad=_Q<10>(),_quad1=_Q<10>() );
        else if ( normType == "SemiH1-error" )
            normComputed = normSemiH1(_range=range,_grad_expr= gradv(field) -  ppNorm.gradSolution().template expr<FieldType::nComponents1,FieldType::nRealDim>(),
                                      _quad=_Q<10>(),_quad1=_Q<10>() );
        else
            CHECK( false ) << "error";
        res[normNameOutput] = normComputed;
    }
}

template<typename RangeType, typename FieldType>
void
measureNormEvaluation( RangeType const& range, FieldType const& field, ModelPostprocessNorm const& ppNorm,  std::map<std::string,double> & res,
                       typename std::enable_if< FieldType::is_tensor2 || FieldType::is_tensor2symm>::type* = nullptr )
{
    for ( std::string const& normType : ppNorm.types() )
    {
        std::string normNameOutput = (boost::format("Norm_%1%_%2%")%ppNorm.name() %normType).str();
        double normComputed = 0;
        if ( normType == "L2" )
            normComputed = normL2(_range=range,_expr=idv(field) );
        else if ( normType == "L2-error" )
            normComputed = normL2(_range=range,_expr=idv(field) - ppNorm.solution().template expr<FieldType::nComponents1,FieldType::nComponents2>(),
                                  _quad=_Q<10>(),_quad1=_Q<10>() );
        else if ( normType == "H1" || normType == "SemiH1" || normType == "H1-error" || normType == "SemiH1-error"   )
            CHECK( false ) << "normType " << normType << " is not implemented with tensor field";
        else
            CHECK( false ) << "error";
        res[normNameOutput] = normComputed;
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
