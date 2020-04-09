/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_UPDATESTABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_UPDATESTABILIZATIONGLS_HPP 1

#include <feel/feelvf/vf.hpp>
//#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmesh/intersect.hpp>

#include <feel/feelmodels/modelvf/exproperations.hpp>


namespace Feel
{
namespace FeelModels
{
namespace FluidMechanicsDetail
{

enum FModel { Stokes=0, NavierStokes=1 };



template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                               VelocityFieldType const& u, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + timeSteppingScaling*gradt(u)*beta_u );
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientLinearExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                               VelocityFieldType const& u, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + timeSteppingScaling*gradt(u)*beta_u ) - timeSteppingScaling*mu*laplaciant(u);
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
                                 VelocityFieldType const& u, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + timeSteppingScaling*(gradt(u)*idv(u) + gradv(u)*idt(u)) );
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientJacobianExpr_u( RhoExprType const& rho, ViscosityExprType const& mu,
                                 VelocityFieldType const& u, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + timeSteppingScaling*(gradt(u)*idv(u) + gradv(u)*idt(u)) ) - timeSteppingScaling*mu*laplaciant(u);
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
                                 VelocityFieldType const& u, PressureFieldType const& p, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return rho*( dudt + timeSteppingScaling*gradv(u)*beta_u ) + trans(gradv(p));
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename TimeDerivativeExprType,
         typename VelocityFieldType, typename PressureFieldType, typename FluidMechanicsType >
auto
residualTransientResidualExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u, TimeDerivativeExprType const& dudt,
                                 VelocityFieldType const& u, PressureFieldType const& p, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return rho*( dudt + timeSteppingScaling*gradv(u)*beta_u ) - timeSteppingScaling*mu*laplacianv(u) + trans(gradv(p));
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

template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientResidualWithoutTimeDerivativeExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                                                      VelocityFieldType const& u, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<0> )
{
    return timeSteppingScaling*rho*gradv(u)*beta_u /*+ trans(gradv(p))*/;
}
template<typename RhoExprType, typename ViscosityExprType, typename VelocityConvectionExprType, typename VelocityFieldType, typename FluidMechanicsType >
auto
residualTransientResidualWithoutTimeDerivativeExpr_u( RhoExprType const& rho, ViscosityExprType const& mu, VelocityConvectionExprType const& beta_u,
                                                      VelocityFieldType const& u, double timeSteppingScaling, FluidMechanicsType const& fluidmec, mpl::int_<1> )
{
    return timeSteppingScaling*rho*gradv(u)*beta_u - timeSteppingScaling*mu*laplacianv(u) /*+ trans(gradv(p))*/;
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

    double timeSteppingScaling = 1.;
    if ( !fluidmec.isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(fluidmec.prefix(),"time-stepping.scaling") );

    auto mesh = fluidmec.mesh();
    auto XhV = fluidmec.functionSpaceVelocity();
    auto XhP = fluidmec.functionSpacePressure();
    auto const& u = fluidmec.fieldVelocity();
    auto const& p = fluidmec.fieldPressure();
    auto const& v = u;
    auto const& q = p;

    auto bilinearFormVV_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                                _pattern=size_type(Pattern::COUPLED),
                                                _rowstart=fluidmec.rowStartInMatrix(),
                                                _colstart=fluidmec.colStartInMatrix() );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=fluidmec.rowStartInMatrix(),
                                 _colstart=fluidmec.colStartInMatrix()+1 );
    auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=fluidmec.rowStartInMatrix()+1,
                                 _colstart=fluidmec.colStartInMatrix()+0 );
    auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=fluidmec.rowStartInMatrix()+1,
                                 _colstart=fluidmec.colStartInMatrix()+1 );

    auto myLinearFormV = form1( _test=XhV, _vector=F,
                                _rowstart=fluidmec.rowStartInVector() );
    auto myLinearFormP = form1( _test=XhP, _vector=F,
                                _rowstart=fluidmec.rowStartInVector()+1 );

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
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( rho,mu,uconv,u,timeSteppingScaling, fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormVV_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeDerivative = fluidmec.timeStepBDF()->polyDeriv();
            auto stab_residual_linear = rho*idv(rhsTimeDerivative);
            myLinearFormV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( rho,mu,uconv,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormVV_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearFormVP +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );

        if ( hana::length( addRhsTuple ).value > 0 )
        {
            auto stab_residual_linear = Feel::FeelModels::vfdetail::addExpr( addRhsTuple );
            myLinearFormV +=
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
            myLinearFormV +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=timeSteppingScaling*tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }

        if ( fluidmec.timeStepping() ==  "Theta" )
        {
            auto previousSol = data.vectorInfo( prefixvm( fluidmec.prefix(),"time-stepping.previous-solution") );
            auto uOld = XhV->element( previousSol, fluidmec.rowStartInVector() );
            CHECK( fluidmec.solverName() == "Oseen" ) << "TODO";
            auto previousConv = data.vectorInfo( prefixvm( fluidmec.prefix(),"time-stepping.previous-convection-field-extrapolated") );
            auto uOldConv = XhV->element( previousConv );

            //auto stab_residual_u_old = addExpr( hana::make_tuple( residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOldConv),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            auto stab_residual_u_old = -residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(/*uOld*/uOldConv),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() );
            myLinearFormV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u_old,stab_test),
                           _geomap=fluidmec.geomap() );

            // TODO body forces
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
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( rho,mu,uconv,u,timeSteppingScaling, fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormPV +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeDerivative = fluidmec.timeStepBDF()->polyDeriv();
            auto stab_residual_linear = rho*idv(rhsTimeDerivative);
            myLinearFormP +=
                integrate( _range=rangeEltPressure,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( rho,mu,uconv,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormPV +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearFormPP +=
            integrate( _range=rangeEltPressure,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );

        if ( hana::length( addRhsTuple ).value > 0 )
        {
            auto stab_residual_linear = Feel::FeelModels::vfdetail::addExpr( addRhsTuple );
            myLinearFormP +=
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
            myLinearFormP +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=timeSteppingScaling*tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }

        if ( fluidmec.timeStepping() ==  "Theta" )
        {
            auto previousSol = data.vectorInfo( prefixvm( fluidmec.prefix(),"time-stepping.previous-solution") );
            auto uOld = XhV->element( previousSol, fluidmec.rowStartInVector() );
            CHECK( fluidmec.solverName() == "Oseen" ) << "TODO";
            auto previousConv = data.vectorInfo( prefixvm( fluidmec.prefix(),"time-stepping.previous-convection-field-extrapolated") );
            auto uOldConv = XhV->element( previousConv );

            //auto stab_residual_u_old = addExpr( hana::make_tuple( residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOldConv),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            auto stab_residual_u_old = -residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOldConv),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() );
            myLinearFormP +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u_old,stab_test),
                           _geomap=fluidmec.geomap() );

            // TODO body forces
        }

    }
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLS", "finish" );
}


template< int StabGLSType,typename FluidMechanicsType, typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
updateJacobianStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateJacobian & data,
                                typename FluidMechanicsType::element_velocity_external_storage_type const& u,
                                typename FluidMechanicsType::element_pressure_external_storage_type const& p,
                                Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                std::string const& matName, const ExprT&... exprs )
{
    static const int StabResidualType = ( FluidMechanicsType::nOrderVelocity>1 )? 1 : 0;

    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    double timeSteppingScaling = 1.;
    if ( !fluidmec.isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(fluidmec.prefix(),"time-stepping.scaling") );

    auto mesh = fluidmec.mesh();
    auto XhV = fluidmec.functionSpaceVelocity();
    auto XhP = fluidmec.functionSpacePressure();
    auto bilinearFormVV_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=J,
                                                _pattern=size_type(Pattern::COUPLED),
                                                _rowstart=fluidmec.rowStartInMatrix(),
                                                _colstart=fluidmec.colStartInMatrix() );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=fluidmec.rowStartInMatrix(),
                                 _colstart=fluidmec.colStartInMatrix()+1 );
    auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=fluidmec.rowStartInMatrix()+1,
                                 _colstart=fluidmec.colStartInMatrix()+0 );
    auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=fluidmec.rowStartInMatrix()+1,
                                 _colstart=fluidmec.colStartInMatrix()+1 );

    auto uconv = rho*idv(u);

    auto additionTermsTuple = hana::make_tuple(exprs...);

    bool hasUpdatedTauForConvectionDiffusion = false;
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion( matName );
#if 0
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
#else
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        auto tau = idv(tauFieldPtr);
#endif
        auto stab_test = stabGLStestLinearExpr_u( rho,mu,idv(u),u, mpl::int_<StabGLSType>() );
        if (!fluidmec.isStationaryModel())
        {
            auto stab_residual_bilinear_u = residualTransientJacobianExpr_u( rho,mu,u,timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormVV_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryJacobianExpr_u( rho,mu,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormVV_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearFormVP +=
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
#else
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterPressure()->fieldTauPtr();
        auto tau = idv(tauFieldPtr);
#endif
        auto stab_test = -trans(grad(p));
        if ( !fluidmec.isStationaryModel() )
        {
            auto stab_residual_bilinear_u = residualTransientJacobianExpr_u( rho,mu,u,timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormPV +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryJacobianExpr_u( rho,mu,u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearFormPV +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearFormPP +=
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
                                typename FluidMechanicsType::element_velocity_external_storage_type const& u,
                                typename FluidMechanicsType::element_pressure_external_storage_type const& p,
                                Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                std::string const& matName, const ExprT&... exprs )
{
    fluidmec.log("FluidMechanics","updateResidualStabilizationGLS", "start" );
    static const int StabResidualType = ( FluidMechanicsType::nOrderVelocity>1 )? 1 : 0;

    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !fluidmec.isStationaryModel() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(fluidmec.prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        timeSteppingScaling = data.doubleInfo( prefixvm(fluidmec.prefix(),"time-stepping.scaling") );
    }
    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
        return;

    auto mesh = fluidmec.mesh();
    auto XhV = fluidmec.functionSpaceVelocity();
    auto XhP = fluidmec.functionSpacePressure();
    auto myLinearFormV = form1( _test=XhV, _vector=R,
                               _rowstart=fluidmec.rowStartInVector() );
    auto myLinearFormP = form1( _test=XhP, _vector=R,
                               _rowstart=fluidmec.rowStartInVector()+1 );

    auto uconv = rho*idv(u);
    bool hasUpdatedTauForConvectionDiffusion = false;
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion( matName );
        //auto stab_test = grad(u)*uconv - mu*laplacian(u);
        auto stab_test = stabGLStestLinearExpr_u( rho,mu,idv(u),u, mpl::int_<StabGLSType>() );

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

        if ( !fluidmec.isStationaryModel() )
        {
            auto const& rhsTimeDerivative = fluidmec.timeStepBDF()->polyDeriv();
            auto dudt = idv(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) - idv(rhsTimeDerivative);
            auto stab_residual_u = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( residualTransientResidualExpr_u( rho,mu,idv(u),dudt,u,p,timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearFormV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_u = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( residualStationaryResidualExpr_u( rho,mu,idv(u),u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearFormV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? rangeEltConvectionDiffusion : intersect(markedelements(mesh,markers(d)),rangeEltConvectionDiffusion);
            myLinearFormV +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-timeSteppingScaling*tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }

        if ( fluidmec.timeStepping() ==  "Theta" )
        {
            auto previousSol = data.vectorInfo( prefixvm( fluidmec.prefix(),"time-stepping.previous-solution") );
            auto uOld = XhV->element( previousSol, fluidmec.rowStartInVector() );

            //auto stab_residual_u_old = addExpr( hana::make_tuple( residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOldConv),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            auto stab_residual_u_old = residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOld),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() );
            myLinearFormV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_u_old,stab_test),
                           _geomap=fluidmec.geomap() );

            // TODO body forces
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
            auto const& rhsTimeDerivative = fluidmec.timeStepBDF()->polyDeriv();
            auto dudt = idv(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) - idv(rhsTimeDerivative);
            auto stab_residual_u = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( residualTransientResidualExpr_u( rho,mu,idv(u),dudt,u,p,timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearFormP +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_u = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( residualStationaryResidualExpr_u( rho,mu,idv(u),u,p,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            myLinearFormP +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? rangeEltPressure : intersect(markedelements(mesh,markers(d)),rangeEltPressure);
            myLinearFormP +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-timeSteppingScaling*tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }

        if ( fluidmec.timeStepping() ==  "Theta" )
        {
            auto previousSol = data.vectorInfo( prefixvm( fluidmec.prefix(),"time-stepping.previous-solution") );
            auto uOld = XhV->element( previousSol, fluidmec.rowStartInVector() );

            //auto stab_residual_u_old = addExpr( hana::make_tuple( residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOldConv),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() ), exprs... ) );
            auto stab_residual_u_old = residualTransientResidualWithoutTimeDerivativeExpr_u( rho,mu,idv(uOld),uOld,1.0 - timeSteppingScaling,fluidmec, mpl::int_<StabResidualType>() );
            myLinearFormP +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_u_old,stab_test),
                           _geomap=fluidmec.geomap() );

            // TODO body forces
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
    CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(this->physicsFromCurrentType().begin()->second);

    if ( physicFluidData->equation() == "Navier-Stokes")
    {

        element_velocity_ptrtype fielCurrentPicardSolution;
        if ( this->solverName() == "Picard" )
        {
            const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
            fielCurrentPicardSolution = M_XhVelocity->elementPtr();
            *fielCurrentPicardSolution = *M_XhVelocity->elementPtr(*vecCurrentPicardSolution, this->rowStartInVector() );
        }
        auto const& betaU = ( this->solverName() == "Oseen" )? *M_fieldConvectionVelocityExtrapolated/*this->timeStepBDF()->poly()*/ : *fielCurrentPicardSolution;
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
    else if ( physicFluidData->equation() == "Stokes" || physicFluidData->equation() == "StokesTransient" )
    {
        auto uconv = vf::zero<nRealDim,1>();
        if ( this->stabilizationGLSType() == "pspg" || ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<FluidMechanicsDetail::FModel::Stokes,0>( *this, data, rho, mu, uconv, matName, addRhsTuple, addMatTuple );
        else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<FluidMechanicsDetail::FModel::Stokes,1>( *this, data, rho, mu, uconv, matName, addRhsTuple, addMatTuple );
    }
}



template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateJacobianStabilisationGLS( DataUpdateJacobian & data,
                                                                                                            element_velocity_external_storage_type const& u,
                                                                                                            element_pressure_external_storage_type const& p,
                                                                                                            Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                            std::string const& matName, const ExprT&... exprs ) const
{
    CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(this->physicsFromCurrentType().begin()->second);
    CHECK( physicFluidData->equation() == "Navier-Stokes") << "only implemented for Navier-Stokes";

    if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
         ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        FluidMechanicsDetail::updateJacobianStabilizationGLS<0>( *this, data, u, p, rho, mu, matName, exprs... );
    }
    else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
    {
        FluidMechanicsDetail::updateJacobianStabilizationGLS<1>( *this, data, u, p, rho, mu, matName, exprs... );
    }
}


template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::updateResidualStabilisationGLS( DataUpdateResidual & data,
                                                                                                            element_velocity_external_storage_type const& u,
                                                                                                            element_pressure_external_storage_type const& p,
                                                                                                            Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                                                                                            std::string const& matName, const ExprT&... exprs ) const
{
    CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(this->physicsFromCurrentType().begin()->second);
    CHECK( physicFluidData->equation() == "Navier-Stokes") << "only implemented for Navier-Stokes";

    if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
         ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        FluidMechanicsDetail::updateResidualStabilizationGLS<0>( *this, data, u, p, rho, mu, matName, exprs... );
    }
    else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
    {
        FluidMechanicsDetail::updateResidualStabilizationGLS<1>( *this, data, u, p, rho, mu, matName, exprs... );
    }
}

} // namespace FeelModels
} // namespace Feel

#endif
