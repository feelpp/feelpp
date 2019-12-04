/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_HEAT_UPDATESTABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_HEAT_UPDATESTABILIZATIONGLS_HPP 1

#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feeldiscr/pchv.hpp>
//#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelmesh/intersect.hpp>

#include <feel/feelmodels/modelvf/exproperations.hpp>

namespace Feel
{
namespace FeelModels
{

namespace HeatDetail_StabGLS
{

template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientLinearExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                             TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<0> )
{
    return rhocp*(idt(u)*heatToolbox.timeStepBdfTemperature()->polyDerivCoefficient(0) + timeSteppingScaling*gradt(u)*uconv );
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientLinearExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                             TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<1>,
                             typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_scalar >* = nullptr )
{
    return rhocp*(idt(u)*heatToolbox.timeStepBdfTemperature()->polyDerivCoefficient(0) + timeSteppingScaling*gradt(u)*uconv ) - timeSteppingScaling*kappa*laplaciant(u);
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientLinearExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                             TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<1>,
                             typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_tensor2 >* = nullptr )
{
    return rhocp*(idt(u)*heatToolbox.timeStepBdfTemperature()->polyDerivCoefficient(0) + timeSteppingScaling*gradt(u)*uconv ) - timeSteppingScaling*inner(kappa,hesst(u));
}

template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualStationaryLinearExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                              TemperatureFieldType const& u, HeatType const& heatToolbox, mpl::int_<0> )
{
    return rhocp*gradt(u)*uconv;
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualStationaryLinearExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                              TemperatureFieldType const& u, HeatType const& heatToolbox, mpl::int_<1>,
                              typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_scalar >* = nullptr )
{
    return rhocp*gradt(u)*uconv - kappa*laplaciant(u);
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualStationaryLinearExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                              TemperatureFieldType const& u, HeatType const& heatToolbox, mpl::int_<1>,
                              typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_tensor2 >* = nullptr )
{
    return rhocp*gradt(u)*uconv -inner(kappa,hesst(u));
}

template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TimeDerivativeExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientResidualExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv, TimeDerivativeExprType const& dudt,
                               TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<0> )
{
    return rhocp*( dudt + timeSteppingScaling*gradv(u)*uconv );
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TimeDerivativeExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientResidualExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv, TimeDerivativeExprType const& dudt,
                               TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<1>,
                               typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_scalar >* = nullptr )
{
    return rhocp*( dudt + timeSteppingScaling*gradv(u)*uconv ) - timeSteppingScaling*kappa*laplacianv(u);
}

template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TimeDerivativeExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientResidualExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv, TimeDerivativeExprType const& dudt,
                               TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<1>,
                               typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_tensor2 >* = nullptr )
{
    return rhocp*( dudt + timeSteppingScaling*gradv(u)*uconv ) - timeSteppingScaling*inner(kappa,hessv(u));
}

template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualStationaryResidualExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                                TemperatureFieldType const& u, HeatType const& heatToolbox, mpl::int_<0> )
{
    return rhocp*gradv(u)*uconv;
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualStationaryResidualExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                                TemperatureFieldType const& u, HeatType const& heatToolbox, mpl::int_<1>,
                                typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_scalar >* = nullptr )
{
    return rhocp*gradv(u)*uconv - kappa*laplacianv(u);
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualStationaryResidualExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                                TemperatureFieldType const& u, HeatType const& heatToolbox, mpl::int_<1>,
                                typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_tensor2 >* = nullptr )
{
    return rhocp*gradv(u)*uconv - inner(kappa,hessv(u));
}

template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientResidualWithoutTimeDerivativeExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                                                    TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<0> )
{
    return timeSteppingScaling*rhocp*gradv(u)*uconv;
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientResidualWithoutTimeDerivativeExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                                                    TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<1>,
                                                    typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_scalar >* = nullptr )
{
    return timeSteppingScaling*(rhocp*gradv(u)*uconv - kappa*laplacianv(u));
}
template<typename RhoCpExprType, typename KappaExprType, typename VelocityConvectionExprType, typename TemperatureFieldType, typename HeatType >
auto
residualTransientResidualWithoutTimeDerivativeExpr( RhoCpExprType const& rhocp, KappaExprType const& kappa, VelocityConvectionExprType const& uconv,
                                                    TemperatureFieldType const& u, double timeSteppingScaling, HeatType const& heatToolbox, mpl::int_<1>,
                                                    typename std::enable_if_t< KappaExprType::template evaluator_t<typename HeatType::mesh_type::element_type>::shape::is_tensor2 >* = nullptr )
{
    return timeSteppingScaling*(rhocp*gradv(u)*uconv - inner(kappa,hessv(u)) );
}


template <int StabResidualType, typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType, typename TestFunctionExpr, typename HeatToolbox, typename... ExprT>
void
updateResidualStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                Expr<ConductivityExprType> const& kappa,
                                Expr<ConvectionExprType> const& uconv,
                                RangeType const& range,
                                TestFunctionExpr const& stab_test,
                                ModelAlgebraic::DataUpdateResidual & data,
                                HeatToolbox const& heatToolbox,
                                const ExprT&... exprs )
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    //bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !heatToolbox.isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(heatToolbox.prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        timeSteppingScaling = data.doubleInfo( prefixvm(heatToolbox.prefix(),"time-stepping.scaling") );
    }
    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
        return;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    heatToolbox.log("Heat","updateResidualStabilizationGLS", "start"+sc);


    auto mesh = heatToolbox.mesh();
    auto Xh = heatToolbox.spaceTemperature();
    auto u = Xh->element( XVec, heatToolbox.rowStartInVector() );
    auto const& v = heatToolbox.fieldTemperature();

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=heatToolbox.rowStartInVector() );
#if 1
    auto rhocpuconv = rhocp*uconv;
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *heatToolbox.stabilizationGLSParameter(),rhocpuconv/*rhocp*uconv*/, kappa );
    auto tauFieldPtr = heatToolbox.stabilizationGLSParameter()->fieldTauPtr();
    tauFieldPtr->on(_range=range,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);
#else
    auto tau = idv( heatToolbox.stabilizationGLSParameter()->fieldTauPtr() );
#endif


    if (!heatToolbox.isStationary() )
    {
        auto rhsTimeStep = heatToolbox.timeStepBdfTemperature()->polyDeriv();
        auto dudt = idv(u)*heatToolbox.timeStepBdfTemperature()->polyDerivCoefficient(0) - idv( rhsTimeStep );
        auto stab_residual_basic = HeatDetail_StabGLS::residualTransientResidualExpr( rhocp, kappa, uconv, dudt, u, timeSteppingScaling, heatToolbox, mpl::int_<StabResidualType>() );
        auto stab_residual = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( stab_residual_basic, exprs... ) );
        myLinearForm +=
            integrate( _range=range,
                       _expr=val(tau*stab_residual)*stab_test,
                       _geomap=heatToolbox.geomap() );
    }
    else
    {
        auto stab_residual_basic = HeatDetail_StabGLS::residualStationaryResidualExpr( rhocp, kappa, uconv, u, heatToolbox, mpl::int_<StabResidualType>() );
        auto stab_residual = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( stab_residual_basic, exprs... ) );
        myLinearForm +=
            integrate( _range=range,
                       _expr=val(tau*stab_residual)*stab_test,
                       _geomap=heatToolbox.geomap() );
    }
    for( auto const& d : heatToolbox.bodyForces() )
    {
        auto rangeBodyForceUsed = ( markers(d).empty() )? range : intersect( markedelements(mesh,markers(d)), range );
        myLinearForm +=
            integrate( _range=rangeBodyForceUsed,
                       _expr=-timeSteppingScaling*tau*expression(d,heatToolbox.symbolsExpr())*stab_test,
                       _geomap=heatToolbox.geomap() );
    }

    if ( heatToolbox.timeStepping() == "Theta" )
    {
        if ( heatToolbox.fieldVelocityConvectionIsUsedAndOperational() )
        {
            auto previousSol = data.vectorInfo( prefixvm( heatToolbox.prefix(),"time-stepping.previous-solution") );
            auto tOld = Xh->element( previousSol, heatToolbox.rowStartInVector() );
            auto previousConv = data.vectorInfo( prefixvm( heatToolbox.prefix(),"time-stepping.previous-convection-velocity-field") );
            auto uConvOld = heatToolbox.fieldVelocityConvection().functionSpace()->element( previousConv );
            //auto stab_residual_old = (1.0 - timeSteppingScaling)*rhocp*gradv(tOld)*idv(uConvOld);
            auto stab_residual_old = HeatDetail_StabGLS::residualTransientResidualWithoutTimeDerivativeExpr( rhocp, kappa, idv(uConvOld), tOld, 1.0-timeSteppingScaling, heatToolbox, mpl::int_<StabResidualType>() );
            myLinearForm +=
                integrate( _range=range,
                           _expr=tau*stab_residual_old*stab_test,
                           _geomap=heatToolbox.geomap() );
        }
        else if constexpr ( StabResidualType == 1  && Expr<ConvectionExprType>::template evaluator_t<typename HeatToolbox::mesh_type::element_type>::shape::is_scalar ) // TODO kappa is matrix!!!
        {
            auto previousSol = data.vectorInfo( prefixvm( heatToolbox.prefix(),"time-stepping.previous-solution") );
            auto tOld = Xh->element( previousSol, heatToolbox.rowStartInVector() );
            auto stab_residual_old = -(1.0-timeSteppingScaling)*kappa*laplacianv(tOld);
            myLinearForm +=
                integrate( _range=range,
                           _expr=tau*stab_residual_old*stab_test,
                           _geomap=heatToolbox.geomap() );

        }
        // TODO body forces
    }

    heatToolbox.log("Heat","updateResidualStabilizationGLS", "finish"+sc);
}


template <int StabResidualType,typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType,typename TestFunctionExpr,typename HeatToolbox>
void
updateLinearPDEStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                 Expr<ConductivityExprType> const& kappa,
                                 Expr<ConvectionExprType> const& uconv,
                                 RangeType const& range,
                                 TestFunctionExpr const& stab_test,
                                 ModelAlgebraic::DataUpdateLinear & data,
                                 HeatToolbox const& heatToolbox )
{
    //static const int StabResidualType = ( nOrderTemperature>1 )? 1 : 0;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    heatToolbox.log("Heat","updateLinearPDEStabilizationGLS", "start"+sc);

    double timeSteppingScaling = 1.;
    if ( !heatToolbox.isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(heatToolbox.prefix(),"time-stepping.scaling") );

    auto mesh = heatToolbox.mesh();
    auto Xh = heatToolbox.spaceTemperature();
    auto const& u = heatToolbox.fieldTemperature();
    auto const& v = heatToolbox.fieldTemperature();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=heatToolbox.rowStartInMatrix(),
                                              _colstart=heatToolbox.colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=heatToolbox.rowStartInVector() );



    auto rhocpuconv = rhocp*uconv;
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *heatToolbox.stabilizationGLSParameter(),rhocpuconv/*rhocp*uconv*/, kappa );
#if 1
    auto tauFieldPtr = heatToolbox.stabilizationGLSParameter()->fieldTauPtr();
    tauFieldPtr->on(_range=range,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);
#else
    auto tau = tauExpr;
#endif

    if (!heatToolbox.isStationary())
    {
        auto stab_residual_bilinear = HeatDetail_StabGLS::residualTransientLinearExpr( rhocp, kappa, uconv, u, timeSteppingScaling, heatToolbox, mpl::int_<StabResidualType>() );
        bilinearForm_PatternCoupled +=
            integrate( _range=range,
                       _expr=tau*stab_residual_bilinear*stab_test,
                       _geomap=heatToolbox.geomap() );
        auto rhsTimeStep = heatToolbox.timeStepBdfTemperature()->polyDeriv();
        auto stab_residual_linear = rhocp*idv( rhsTimeStep );
        myLinearForm +=
            integrate( _range=range,
                       _expr= val(tau*stab_residual_linear)*stab_test,
                       _geomap=heatToolbox.geomap() );
    }
    else
    {
        auto stab_residual_bilinear = HeatDetail_StabGLS::residualStationaryLinearExpr( rhocp, kappa, uconv, u, heatToolbox, mpl::int_<StabResidualType>() );
        bilinearForm_PatternCoupled +=
            integrate( _range=range,
                       _expr=tau*stab_residual_bilinear*stab_test,
                       _geomap=heatToolbox.geomap() );
    }
    for( auto const& d : heatToolbox.bodyForces() )
    {
        auto rangeBodyForceUsed = ( markers(d).empty() )? range : intersect( markedelements(mesh,markers(d)),range );
        myLinearForm +=
            integrate( _range=rangeBodyForceUsed,
                       _expr=timeSteppingScaling*tau*expression(d,heatToolbox.symbolsExpr())*stab_test,
                       _geomap=heatToolbox.geomap() );
    }

    if ( heatToolbox.timeStepping() == "Theta" )
    {
        auto previousSol = data.vectorInfo( prefixvm( heatToolbox.prefix(),"time-stepping.previous-solution") );
        auto tOld = Xh->element( previousSol, heatToolbox.rowStartInVector() );
        CHECK( heatToolbox.fieldVelocityConvectionIsUsedAndOperational() ) << "something wrong";
        auto previousConv = data.vectorInfo( prefixvm( heatToolbox.prefix(),"time-stepping.previous-convection-velocity-field") );
        auto uConvOld = heatToolbox.fieldVelocityConvection().functionSpace()->element( previousConv );
        //auto stab_residual_old = (1.0 - timeSteppingScaling)*rhocp*gradv(tOld)*idv(uConvOld);
        auto stab_residual_old = HeatDetail_StabGLS::residualTransientResidualWithoutTimeDerivativeExpr( rhocp, kappa, idv(uConvOld), tOld, 1.0-timeSteppingScaling, heatToolbox, mpl::int_<StabResidualType>() );
        myLinearForm +=
            integrate( _range=range,
                       _expr=-tau*stab_residual_old*stab_test,
                       _geomap=heatToolbox.geomap() );
        // TODO body forces
    }

    heatToolbox.log("Heat","updateLinearPDEStabilizationGLS", "finish"+sc);
}


template <int StabResidualType,typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType,typename TestFunctionExpr,typename HeatToolbox>
void
updateJacobianStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                Expr<ConductivityExprType> const& kappa,
                                Expr<ConvectionExprType> const& uconv,
                                RangeType const& range,
                                TestFunctionExpr const& stab_test,
                                ModelAlgebraic::DataUpdateJacobian & data,
                                HeatToolbox const& heatToolbox )
{
    //static const int StabResidualType = ( nOrderTemperature>1 )? 1 : 0;

    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    heatToolbox.log("Heat","updateJacobianStabilizationGLS", "start"+sc);

    double timeSteppingScaling = 1.;
    if ( !heatToolbox.isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(heatToolbox.prefix(),"time-stepping.scaling") );

    auto mesh = heatToolbox.mesh();
    auto Xh = heatToolbox.spaceTemperature();
    auto u = Xh->element( XVec, heatToolbox.rowStartInVector() );
    auto const& v = heatToolbox.fieldTemperature();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=heatToolbox.rowStartInMatrix(),
                               _colstart=heatToolbox.colStartInMatrix() );
#if 0
    auto rhocpuconv = rhocp*uconv;
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *heatToolbox.stabilizationGLSParameter(),rhocpuconv/*rhocp*uconv*/, kappa );
    auto tauFieldPtr = heatToolbox.stabilizationGLSParameter()->fieldTauPtr();
    tauFieldPtr->on(_range=range,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);
#else
    auto tau = idv( heatToolbox.stabilizationGLSParameter()->fieldTauPtr() );
#endif

    if (!heatToolbox.isStationary())
    {
        auto stab_residual_bilinear = HeatDetail_StabGLS::residualTransientLinearExpr( rhocp, kappa, uconv, u, timeSteppingScaling, heatToolbox, mpl::int_<StabResidualType>() );
        bilinearForm +=
            integrate( _range=range,
                       _expr=tau*stab_residual_bilinear*stab_test,
                       _geomap=heatToolbox.geomap() );
    }
    else
    {
        auto stab_residual_bilinear = HeatDetail_StabGLS::residualStationaryLinearExpr( rhocp, kappa, uconv, u, heatToolbox, mpl::int_<StabResidualType>() );
        bilinearForm +=
            integrate( _range=range,
                       _expr=tau*stab_residual_bilinear*stab_test,
                       _geomap=heatToolbox.geomap() );
    }

    heatToolbox.log("Heat","updateJacobianStabilizationGLS", "finish"+sc);

}


} // namespace HeatDetail_StabGLS


template< typename ConvexType, typename BasisTemperatureType >
template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateJacobianStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                                                       Expr<ConductivityExprType> const& kappa,
                                                                       Expr<ConvectionExprType> const& uconv,
                                                                       RangeType const& range,
                                                                       DataUpdateJacobian & data ) const
{
    static const int StabResidualType = ( nOrderTemperature>1 )? 1 : 0;

    auto const& u = this->fieldTemperature();
    if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg" )
    {
        auto stab_test = rhocp*grad(u)*uconv;
        HeatDetail_StabGLS::updateJacobianStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this );
    }
    else
    {
        int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? -1 : 1;
        if constexpr ( Expr<ConductivityExprType>::template evaluator_t<typename mesh_type::element_type>::shape::is_scalar )
            {
                auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
                HeatDetail_StabGLS::updateJacobianStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this );
            }
        else
        {
            auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*inner(kappa,hess(u));
            HeatDetail_StabGLS::updateJacobianStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this );
        }
    }
}

template< typename ConvexType, typename BasisTemperatureType >
template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType, typename... ExprT>
void
Heat<ConvexType,BasisTemperatureType>::updateResidualStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                                                       Expr<ConductivityExprType> const& kappa,
                                                                       Expr<ConvectionExprType> const& uconv,
                                                                       RangeType const& range,
                                                                       DataUpdateResidual & data,
                                                                       const ExprT&... exprs ) const
{
    static const int StabResidualType = ( nOrderTemperature>1 )? 1 : 0;

    auto const& u = this->fieldTemperature();
    if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg" )
    {
        auto stab_test = rhocp*grad(u)*uconv;
        HeatDetail_StabGLS::updateResidualStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this, exprs... );

    }
    else
    {
        int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? -1 : 1;
        if constexpr ( Expr<ConductivityExprType>::template evaluator_t<typename mesh_type::element_type>::shape::is_scalar )
            {
                auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
                HeatDetail_StabGLS::updateResidualStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this, exprs... );
            }
        else
        {
            auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*inner(kappa,hess(u));
            HeatDetail_StabGLS::updateResidualStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this, exprs... );
        }
    }
}

template< typename ConvexType, typename BasisTemperatureType >
template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateLinearPDEStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                                                        Expr<ConductivityExprType> const& kappa,
                                                                        Expr<ConvectionExprType> const& uconv,
                                                                        RangeType const& range,
                                                                        DataUpdateLinear & data ) const
{
    static const int StabResidualType = ( nOrderTemperature>1 )? 1 : 0;
    auto const& u = this->fieldTemperature();
    if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg" )
    {
        auto stab_test = rhocp*grad(u)*uconv;
        HeatDetail_StabGLS::updateLinearPDEStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this );
    }
    else
    {
        int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? -1 : 1;
        if constexpr ( Expr<ConductivityExprType>::template evaluator_t<typename mesh_type::element_type>::shape::is_scalar )
            {
                auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
                HeatDetail_StabGLS::updateLinearPDEStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this );
            }
        else
        {
            auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*inner(kappa,hess(u));
            HeatDetail_StabGLS::updateLinearPDEStabilizationGLS<StabResidualType>( rhocp, kappa, uconv, range, stab_test, data, *this );
        }
    }
}


} // namespace Feel
} // namespace FeelModels

#endif
