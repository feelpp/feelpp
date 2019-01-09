/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_UPDATESTABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_UPDATESTABILIZATIONGLS_HPP 1

#include <feel/feelvf/vf.hpp>
//#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmesh/intersect.hpp>

namespace Feel
{
namespace FeelModels
{
namespace FluidMechanicsDetail
{

enum FModel { Stokes=0, NavierStokes=1 };

template<typename... Dummy>
auto
addExpr( hana::tuple<> const& /**/ )
{
    CHECK( false ) << "not allow";
    return 0*one();
}
template<typename T1>
auto
addExpr( hana::tuple<T1> const& t )
{
    return hana::at_c<0>( t );
}
template<typename T1,typename T2,typename... TOther>
auto
addExpr( hana::tuple<T1,T2,TOther...> const& t )
{
    return hana::at_c<0>( t ) + addExpr( hana::remove_at_c<0>( t ) );
}

template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                               VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + gradt(u)*beta_u );
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                               VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + gradt(u)*beta_u ) - mu*laplaciant(u);
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualStationaryLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                                VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*gradt(u)*beta_u;
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualStationaryLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                                VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*gradt(u)*beta_u - mu*laplaciant(u);
}

template<typename RhoExprType, typename ViscosityExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientJacobianExpr_u( RhoExprType const& rho, ViscosityExprType const& mu,
                                 VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + gradt(u)*idv(u) + gradv(u)*idt(u) );
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientJacobianExpr_u( RhoExprType const& rho, ViscosityExprType const& mu,
                                 VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + gradt(u)*idv(u) + gradv(u)*idt(u) ) - mu*laplaciant(u);
}

template<typename RhoExprType, typename ViscosityExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualStationaryJacobianExpr_u( RhoExprType const& rho, ViscosityExprType const& mu,
                                 VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*(gradt(u)*idv(u) + gradv(u)*idt(u) );
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualStationaryJacobianExpr_u( RhoExprType const& rho, ViscosityExprType const& mu,
                                  VelocityFieldType const& u, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*(gradt(u)*idv(u) + gradv(u)*idt(u) ) - mu*laplaciant(u);
}

template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename TimeDerivativeExprType,
         typename VelocityFieldType, typename PressureFieldType, typename FluidMechanicsType >
auto
residualTransientResidualExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u, TimeDerivativeExprType const& dudt,
                                 VelocityFieldType const& u, PressureFieldType const& p, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*( dudt + gradv(u)*beta_u ) + trans(gradv(p));
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename TimeDerivativeExprType,
         typename VelocityFieldType, typename PressureFieldType, typename FluidMechanicsType >
auto
residualTransientResidualExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u, TimeDerivativeExprType const& dudt,
                                 VelocityFieldType const& u, PressureFieldType const& p, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*( dudt + gradv(u)*beta_u ) - mu*laplacianv(u) + trans(gradv(p));
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename PressureFieldType, typename FluidMechanicsType >
auto
residualStationaryResidualExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                                  VelocityFieldType const& u, PressureFieldType const& p, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*gradv(u)*beta_u + trans(gradv(p));
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename PressureFieldType, typename FluidMechanicsType >
auto
residualStationaryResidualExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                                  VelocityFieldType const& u, PressureFieldType const& p, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*gradv(u)*beta_u - mu*laplacianv(u) + trans(gradv(p));
}



template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType >
auto
stabGLStestLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                         VelocityFieldType const& u, mpl::int_<0> )
{
    return rho*grad(u)*beta_u;
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType >
auto
stabGLStestLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                         VelocityFieldType const& u, mpl::int_<1> )
{
    return rho*grad(u)*beta_u - mu*laplacian(u);
}



template< int FModelType, int StabGLSType,typename FluidMechanicsType, typename DensityExprType, typename ViscosityExprType, typename ConvectionExprType,
          typename AdditionalRhsType, typename AdditionalMatType >
void
updateLinearPDEStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateLinear & data,
                                 Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu, Expr<ConvectionExprType> const& uconv,
                                 std::string const& matName, AdditionalRhsType const& addRhsTuple, AdditionalMatType const& addMatTuple )
{
    static const int StabResidualType = ( FluidMechanicsType::nOrderVelocity>1 )? 1 : 0;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLS", "start" );

    static const bool hasConvection = (FModelType == FModel::NavierStokes);

    auto mesh = fluidmec.mesh();
    auto Xh = fluidmec.functionSpace();
    auto const& U = fluidmec.fieldVelocityPressure();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=fluidmec.rowStartInMatrix(),
                                              _colstart=fluidmec.colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=fluidmec.rowStartInVector() );

    auto rhouconv = rho*uconv;

    bool hasUpdatedTauForConvectionDiffusion = false;
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion( matName );
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderVelocity> stab_gls_parameter_impl_type;
        auto stabGLSParamConvectionDiffusion =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterConvectionDiffusion() );
        stabGLSParamConvectionDiffusion->updateTau(rhouconv, mu, rangeEltConvectionDiffusion);
        hasUpdatedTauForConvectionDiffusion = true;
        auto tau = idv(stabGLSParamConvectionDiffusion->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterConvectionDiffusion(), rhouconv, mu, hasConvection );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif
        auto stab_test = stabGLStestLinearExpr_u( rho,mu,uconv,u, mpl::int_<StabGLSType>() );
        if (!fluidmec.isStationaryModel())
        {
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( rho,mu,uconv,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto stab_residual_linear = rho*idv(rhsTimeDerivative);
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( rho,mu,uconv,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearForm_PatternCoupled +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );

        if ( hana::length( addRhsTuple ).value > 0 )
        {
            auto stab_residual_linear = addExpr( addRhsTuple );
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        hana::for_each( addMatTuple, [&rangeEltConvectionDiffusion,&tau,&stab_test,&fluidmec]( auto e )
                        {
                            e.first +=
                                integrate( _range=rangeEltConvectionDiffusion,
                                           _expr=tau*inner(e.second,stab_test ),
                                           _geomap=fluidmec.geomap() );
                        });

        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? rangeEltConvectionDiffusion : intersect(markedelements(mesh,markers(d)),rangeEltConvectionDiffusion);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }
    }

    if ( fluidmec.stabilizationGLSType() == "pspg" || fluidmec.stabilizationGLSType() == "supg-pspg" || fluidmec.stabilizationGLSType() == "gls" )
    {
        auto rangeEltPressure = fluidmec.stabilizationGLSEltRangePressure( matName );
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderPressure> stab_gls_parameter_impl_type;
        auto stabGLSParamPressure =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterPressure() );
        if ( ( FluidMechanicsType::nOrderPressure != FluidMechanicsType::nOrderVelocity ) || !hasUpdatedTauForConvectionDiffusion )
            stabGLSParamPressure->updateTau(rhouconv, mu, rangeEltPressure);
        auto tau = idv(stabGLSParamPressure->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterPressure(), rhouconv, mu, hasConvection );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterPressure()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltPressure,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif

        auto stab_test = -trans(grad(p));
        if ( !fluidmec.isStationaryModel() )
        {
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( rho,mu,uconv,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto stab_residual_linear = rho*idv(rhsTimeDerivative);
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( rho,mu,uconv,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearForm_PatternCoupled +=
            integrate( _range=rangeEltPressure,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );

        if ( hana::length( addRhsTuple ).value > 0 )
        {
            auto stab_residual_linear = addExpr( addRhsTuple );
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        hana::for_each( addMatTuple, [&rangeEltPressure,&tau,&stab_test,&fluidmec]( auto e )
                        {
                            e.first +=
                                integrate( _range=rangeEltPressure,
                                           _expr=tau*inner(e.second,stab_test ),
                                           _geomap=fluidmec.geomap() );
                        });

        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? rangeEltPressure : intersect(markedelements(mesh,markers(d)),rangeEltPressure);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }
    }
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLS", "finish" );
}


template< int StabGLSType,typename FluidMechanicsType, typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
updateJacobianStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateJacobian & data,
                                typename FluidMechanicsType::element_fluid_external_storage_type const& U,
                                Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                std::string const& matName, const ExprT&... exprs )
{
    static const int StabResidualType = ( FluidMechanicsType::nOrderVelocity>1 )? 1 : 0;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    auto mesh = fluidmec.mesh();
    auto Xh = fluidmec.functionSpace();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=fluidmec.rowStartInMatrix(),
                                              _colstart=fluidmec.colStartInMatrix() );

    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto uconv = rho*idv(u);

    auto additionTermsTuple = hana::make_tuple(exprs...);

    bool hasUpdatedTauForConvectionDiffusion = false;
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion( matName );
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderVelocity> stab_gls_parameter_impl_type;
        auto stabGLSParamConvectionDiffusion =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterConvectionDiffusion() );
        stabGLSParamConvectionDiffusion->updateTau(uconv, mu, rangeEltConvectionDiffusion);
        hasUpdatedTauForConvectionDiffusion = true;
        auto tau = idv(stabGLSParamConvectionDiffusion->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterConvectionDiffusion(), uconv, mu );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif
        auto stab_test = stabGLStestLinearExpr_u( rho,mu,idv(u),u, mpl::int_<StabGLSType>() );
        if (!fluidmec.isStationaryModel())
        {
            auto stab_residual_bilinear_u = residualTransientJacobianExpr_u( rho,mu,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryJacobianExpr_u( rho,mu,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearForm_PatternCoupled +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );
        hana::for_each( additionTermsTuple, [&rangeEltConvectionDiffusion,&tau,&stab_test,&fluidmec]( auto & e )
                        {
                            e.first +=
                                integrate( _range=rangeEltConvectionDiffusion,
                                           _expr=tau*inner(e.second,stab_test ),
                                           _geomap=fluidmec.geomap() );
                        });

    }

    if ( fluidmec.stabilizationGLSType() == "pspg" || fluidmec.stabilizationGLSType() == "supg-pspg" || fluidmec.stabilizationGLSType() == "gls" )
    {
        auto rangeEltPressure = fluidmec.stabilizationGLSEltRangePressure( matName );
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderPressure> stab_gls_parameter_impl_type;
        auto stabGLSParamPressure =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterPressure() );
        if ( ( FluidMechanicsType::nOrderPressure != FluidMechanicsType::nOrderVelocity ) || !hasUpdatedTauForConvectionDiffusion )
            stabGLSParamPressure->updateTau(uconv, mu, rangeEltPressure);
        auto tau = idv(stabGLSParamPressure->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterPressure(), uconv, mu );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterPressure()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltPressure,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif
        auto stab_test = -trans(grad(p));
        if ( !fluidmec.isStationaryModel() )
        {
            auto stab_residual_bilinear_u = residualTransientJacobianExpr_u( rho,mu,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryJacobianExpr_u( rho,mu,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearForm_PatternCoupled +=
            integrate( _range=rangeEltPressure,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );

        hana::for_each( additionTermsTuple, [&rangeEltPressure,&tau,&stab_test,&fluidmec]( auto & e )
                        {
                            e.first +=
                                integrate( _range=rangeEltPressure,
                                           _expr=tau*inner(e.second,stab_test ),
                                           _geomap=fluidmec.geomap() );
                        });

    }


}


template< int StabGLSType,typename FluidMechanicsType, typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
updateResidualStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateResidual & data,
                                typename FluidMechanicsType::element_fluid_external_storage_type const& U,
                                Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                std::string const& matName, const ExprT&... exprs )
{
    fluidmec.log("FluidMechanics","updateResidualStabilizationGLS", "start" );
    static const int StabResidualType = ( FluidMechanicsType::nOrderVelocity>1 )? 1 : 0;

    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;
    auto mesh = fluidmec.mesh();
    auto Xh = fluidmec.functionSpace();
    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=fluidmec.rowStartInVector() );

    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto uconv = rho*idv(u);
    bool hasUpdatedTauForConvectionDiffusion = false;
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion( matName );
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderVelocity> stab_gls_parameter_impl_type;
        auto stabGLSParamConvectionDiffusion =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterConvectionDiffusion() );
        stabGLSParamConvectionDiffusion->updateTau(uconv, mu, rangeEltConvectionDiffusion);
        hasUpdatedTauForConvectionDiffusion = true;
        auto tau = idv(stabGLSParamConvectionDiffusion->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterConvectionDiffusion(), uconv, mu );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif
        //auto stab_test = grad(u)*uconv - mu*laplacian(u);
        auto stab_test = stabGLStestLinearExpr_u( rho,mu,idv(u),u, mpl::int_<StabGLSType>() );
        if (!fluidmec.isStationaryModel())
        {
            auto const& rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto dudt = idv(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) - idv(rhsTimeDerivative);
            auto stab_residual_u = addExpr( hana::make_tuple( residualTransientResidualExpr_u( rho,mu,idv(u),dudt,u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_u = addExpr( hana::make_tuple( residualStationaryResidualExpr_u( rho,mu,idv(u),u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? rangeEltConvectionDiffusion : intersect(markedelements(mesh,markers(d)),rangeEltConvectionDiffusion);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }
    }

    if ( fluidmec.stabilizationGLSType() == "pspg" || fluidmec.stabilizationGLSType() == "supg-pspg" || fluidmec.stabilizationGLSType() == "gls" )
    {
        auto rangeEltPressure = fluidmec.stabilizationGLSEltRangePressure( matName );
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderPressure> stab_gls_parameter_impl_type;
        auto stabGLSParamPressure =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterPressure() );
        if ( ( FluidMechanicsType::nOrderPressure != FluidMechanicsType::nOrderVelocity ) || !hasUpdatedTauForConvectionDiffusion )
            stabGLSParamPressure->updateTau(uconv, mu, rangeEltPressure);
        auto tau = idv(stabGLSParamPressure->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterPressure(), uconv, mu );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterPressure()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltPressure,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif

        auto stab_test = -trans(grad(p));
        if ( !fluidmec.isStationaryModel() )
        {
            auto const& rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto dudt = idv(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) - idv(rhsTimeDerivative);
            auto stab_residual_u = addExpr( hana::make_tuple( residualTransientResidualExpr_u( rho,mu,idv(u),dudt,u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_u = addExpr( hana::make_tuple( residualStationaryResidualExpr_u( rho,mu,idv(u),u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? rangeEltPressure : intersect(markedelements(mesh,markers(d)),rangeEltPressure);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }
    }

    fluidmec.log("FluidMechanics","updateResidualStabilizationGLS", "finish" );
}

} // namespace FluidMechanicsDetail



template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template<typename DensityExprType, typename ViscosityExprType, typename AdditionalRhsType, typename AdditionalMatType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateLinearPDEStabilisationGLS( DataUpdateLinear & data,
                                                                                                             Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                             std::string const& matName,
                                                                                                             AdditionalRhsType const& addRhsTuple, AdditionalMatType const& addMatTuple ) const
{
    if ( this->modelName() == "Navier-Stokes")
    {

        element_fluid_ptrtype fielCurrentPicardSolution;
        if ( this->solverName() == "Picard" )
        {
            auto Xh = this->functionSpace();
            const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
            fielCurrentPicardSolution = Xh->elementPtr();
            *fielCurrentPicardSolution = *Xh->elementPtr(*vecCurrentPicardSolution, this->rowStartInVector() );
        }
        auto const& BetaU = ( this->solverName() == "Oseen" )? this->timeStepBDF()->poly() : *fielCurrentPicardSolution;
        auto betaU = BetaU.template element<0>();
        auto uconv = idv( betaU );

        if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
             ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
        {
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<FluidMechanicsDetail::FModel::NavierStokes,0>( *this, data, rho, mu, uconv, matName, addRhsTuple, addMatTuple );
        }
        else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
        {
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<FluidMechanicsDetail::FModel::NavierStokes,1>( *this, data, rho, mu, uconv, matName, addRhsTuple, addMatTuple );
        }
    }
    else if ( this->modelName() == "Stokes" || this->modelName() == "StokesTransient" )
    {
        auto uconv = zero<nRealDim,1>();
        if ( this->stabilizationGLSType() == "pspg" || ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<FluidMechanicsDetail::FModel::Stokes,0>( *this, data, rho, mu, uconv, matName, addRhsTuple, addMatTuple );
        else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<FluidMechanicsDetail::FModel::Stokes,1>( *this, data, rho, mu, uconv, matName, addRhsTuple, addMatTuple );
    }
}



template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateJacobianStabilisationGLS( DataUpdateJacobian & data, element_fluid_external_storage_type const& U,
                                                                                                            Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                            std::string const& matName, const ExprT&... exprs ) const
{
    CHECK( this->modelName() == "Navier-Stokes") << "only implemented for Navier-Stokes";

    if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
         ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        FluidMechanicsDetail::updateJacobianStabilizationGLS<0>( *this, data, U, rho, mu, matName, exprs... );
    }
    else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
    {
        FluidMechanicsDetail::updateJacobianStabilizationGLS<1>( *this, data, U, rho, mu, matName, exprs... );
    }
}


template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateResidualStabilisationGLS( DataUpdateResidual & data, element_fluid_external_storage_type const& U,
                                                                                                            Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                            std::string const& matName, const ExprT&... exprs ) const
{
    CHECK( this->modelName() == "Navier-Stokes") << "only implemented for Navier-Stokes";

    if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
         ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        FluidMechanicsDetail::updateResidualStabilizationGLS<0>( *this, data, U, rho, mu, matName, exprs... );
    }
    else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
    {
        FluidMechanicsDetail::updateResidualStabilizationGLS<1>( *this, data, U, rho, mu, matName, exprs... );
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
