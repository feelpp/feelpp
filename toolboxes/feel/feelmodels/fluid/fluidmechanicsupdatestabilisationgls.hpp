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

template<typename T1>
auto
addExpr( T1 const& t1 )
{
    return t1;
}
template<typename T1,typename T2>
auto
addExpr( T1 const& t1, T2 const& t2 )
{
    return t1+t2;
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


template< int StabGLSType,typename FluidMechanicsType, typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
updateJacobianStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateJacobian & data,
                                typename FluidMechanicsType::element_fluid_external_storage_type const& U,
                                Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                typename FluidMechanicsType::range_elements_type const& range, const ExprT&... exprs )
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
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion();

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
        auto rangeEltPressure = fluidmec.stabilizationGLSEltRangePressure();
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
        auto stab_test = -rho*trans(grad(p));
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
                                typename FluidMechanicsType::range_elements_type const& range, const ExprT&... exprs )
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
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion();
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
            auto stab_residual_u = addExpr( residualTransientResidualExpr_u( rho,mu,idv(u),dudt,u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... );
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_u = addExpr( residualStationaryResidualExpr_u( rho,mu,idv(u),u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... );
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( marker(d).empty() )? rangeEltConvectionDiffusion : intersect(markedelements(mesh,marker(d)),rangeEltConvectionDiffusion);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }
    }

    if ( fluidmec.stabilizationGLSType() == "pspg" || fluidmec.stabilizationGLSType() == "supg-pspg" || fluidmec.stabilizationGLSType() == "gls" )
    {
        auto rangeEltPressure = fluidmec.stabilizationGLSEltRangePressure();

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

        auto stab_test = -rho*trans(grad(p));
        if ( !fluidmec.isStationaryModel() )
        {
            auto const& rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto dudt = idv(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) - idv(rhsTimeDerivative);
            auto stab_residual_u = addExpr( residualTransientResidualExpr_u( rho,mu,idv(u),dudt,u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... );
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_u = addExpr( residualStationaryResidualExpr_u( rho,mu,idv(u),u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... );
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( marker(d).empty() )? rangeEltPressure : intersect(markedelements(mesh,marker(d)),rangeEltPressure);
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
template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateJacobianStabilisationGLS( DataUpdateJacobian & data, element_fluid_external_storage_type const& U,
                                                                                                            Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                            range_elements_type const& range, const ExprT&... exprs ) const
{
    CHECK( this->modelName() == "Navier-Stokes") << "only implemented for Navier-Stokes";

    if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
         ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        FluidMechanicsDetail::updateJacobianStabilizationGLS<0>( *this, data, U, rho, mu, range, exprs... );
    }
    else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
    {
        FluidMechanicsDetail::updateJacobianStabilizationGLS<1>( *this, data, U, rho, mu, range, exprs... );
    }
}


template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateResidualStabilisationGLS( DataUpdateResidual & data, element_fluid_external_storage_type const& U,
                                                                                                            Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                            range_elements_type const& range, const ExprT&... exprs ) const
{
    CHECK( this->modelName() == "Navier-Stokes") << "only implemented for Navier-Stokes";

    if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
         ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        FluidMechanicsDetail::updateResidualStabilizationGLS<0>( *this, data, U, rho, mu, range, exprs... );
    }
    else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
    {
        FluidMechanicsDetail::updateResidualStabilizationGLS<1>( *this, data, U, rho, mu, range, exprs... );
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
