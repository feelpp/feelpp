#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_ASSEMBLY_STABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_ASSEMBLY_STABILIZATIONGLS_HPP 1

#include <feel/feelvf/exproptionalconcat.hpp>
#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelvf/fluidmecdivstresstensor.hpp>

#include <feel/feelmodels/modelvf/exproperations.hpp>

namespace Feel
{
namespace FeelModels
{

namespace FluidMechanics_detail
{

template<typename FluidMechanicsType,typename ModelContextType, typename ... ExprAddedType>
auto
exprResidualImpl( FluidMechanicsType const& fm, ModelPhysicFluid<FluidMechanicsType::nDim> const& physicFluidData, MaterialProperties const& matProps,
                  ModelContextType const& mctx, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    auto const& u = mctx.field( FluidMechanicsType::FieldTag::velocity(&fm), "velocity" );
    auto const& p = mctx.field( FluidMechanicsType::FieldTag::pressure(&fm), "pressure" );
    auto const& beta_u = fm.useSemiImplicitTimeScheme() ? mctx.field( FluidMechanicsType::FieldTag::velocity_extrapolated(&fm), "velocity_extrapolated" ) : u;
    auto const& se = mctx.symbolsExpr();

    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;

    using dynamic_viscosity_law_type = DynamicViscosityLaw;
    auto dynamicViscosityLawPtr = std::static_pointer_cast<dynamic_viscosity_law_type>( matProps.law( "dynamic-viscosity" ) );
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(u), *dynamicViscosityLawPtr, matProps, se);

    using expr_convection_residual_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*gradv(u)*idv(beta_u) )>;
    using expr_viscous_stress_residual_type = std::decay_t<decltype( -timeSteppingScaling*fluidMecDivViscousStressTensor(u,physicFluidData,matProps,fm.worldComm(),fm.repository().expr(),se) )>;


    using expr_pressure_stress_redisual_type = std::decay_t<decltype( trans(gradv(p)) )>;
    using expr_time_derivative_residual_type = std::decay_t<decltype( expr_density_type{}*(fm.timeStepBDF()->polyDerivCoefficient(0)*idv(u)-idv(fm.timeStepBDF()->polyDeriv())) )>;
    using expr_gravity_force_type = std::decay_t<decltype( -timeSteppingScaling*expr_density_type{}*physicFluidData.gravityForceExpr() )>;

    auto residual_full = exprOptionalConcat<expr_convection_residual_type,
                                            expr_viscous_stress_residual_type,
                                            expr_pressure_stress_redisual_type,expr_time_derivative_residual_type,expr_gravity_force_type, ExprAddedType...>();

    if ( physicFluidData.equation() == "Navier-Stokes")
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        residual_full.expression().add( timeSteppingScaling*densityExpr*gradv(u)*idv(beta_u) );
    }
    if (!fm.isStationaryModel() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        residual_full.expression().add( densityExpr*(fm.timeStepBDF()->polyDerivCoefficient(0)*idv(u)-idv(fm.timeStepBDF()->polyDeriv())) );
    }

    // important detail, the pressure term is not include in theta scheme
    if ( !timeSteppingEvaluateResidualWithoutTimeDerivative )
        residual_full.expression().add( trans(gradv(p)) );

    auto divSigmaViscous = fluidMecDivViscousStressTensor(u,physicFluidData,matProps,fm.worldComm(),fm.repository().expr(),se,fm.stabilizationGLS_checkViscosityDependencyOnCoordinates());
    if ( !divSigmaViscous.expression().isZero() )
    {
        residual_full.expression().add( -timeSteppingScaling*divSigmaViscous );
    }

    if ( physicFluidData.gravityForceEnabled() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto const& gravityForce = physicFluidData.gravityForceExpr();
        residual_full.expression().add( -timeSteppingScaling*densityExpr*gravityForce );
    }

    return residual_full;
}

template<typename FluidMechanicsType,typename ModelContextType>
auto
exprResidual( FluidMechanicsType const& fm, ModelPhysicFluid<FluidMechanicsType::nDim> const& physicFluidData, MaterialProperties const& matProps,
              ModelContextType const& mctx, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    return exprResidualImpl( fm,physicFluidData,matProps,mctx,timeSteppingScaling,timeSteppingEvaluateResidualWithoutTimeDerivative );
}

template<typename ExprAddedType,typename FluidMechanicsType,typename ModelContextType>
auto
exprResidual( FluidMechanicsType const& fm, ModelPhysicFluid<FluidMechanicsType::nDim> const& physicFluidData, MaterialProperties const& matProps,
              ModelContextType const& mctx, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    return exprResidualImpl<FluidMechanicsType,ModelContextType,ExprAddedType>( fm,physicFluidData,matProps,mctx,timeSteppingScaling,timeSteppingEvaluateResidualWithoutTimeDerivative );
}

} // FluidMechanics_detail

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType,typename RangeType,typename ExprAddedRhsType,typename ExprAddedLhsType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateLinearPDEStabilizationGLS( DataUpdateLinear & data, ModelContextType const& mctx,
                                                                                                 ModelPhysicFluid<nDim> const& physicFluidData,
                                                                                                 MaterialProperties const& matProps, RangeType const& range,
                                                                                                 ExprAddedRhsType const& exprsAddedInResidualRhsTuple,
                                                                                                 ExprAddedLhsType const& exprsAddedInResidualLhsTuple ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );


    std::string const& matName = matProps.materialName();
    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto const& u = mctx.field( FieldTag::velocity(this), "velocity" );
    auto const& v = this->fieldVelocity();
    auto const& p = mctx.field( FieldTag::pressure(this), "pressure" );
    auto const& q = this->fieldPressure();
    auto const& beta_u = this->useSemiImplicitTimeScheme()? mctx.field( FieldTag::velocity_extrapolated(this), "velocity_extrapolated" ) : u;
    auto const& se = mctx.symbolsExpr();

    auto bilinearFormVV = form2( _test=XhV,_trial=XhV,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix() );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix()+1 );
    auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix()+1,
                                 _colstart=this->colStartInMatrix()+0 );
    auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix()+1,
                                 _colstart=this->colStartInMatrix()+1 );

    auto myLinearFormV = form1( _test=XhV, _vector=F,
                                _rowstart=this->rowStartInVector() );
    auto myLinearFormP = form1( _test=XhP, _vector=F,
                                _rowstart=this->rowStartInVector()+1 );


    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;

    using expr_density_uconv_type = std::decay_t<decltype( expr_density_type{}*idv(beta_u) )>;

    auto dynamicViscosityLawPtr = std::static_pointer_cast<dynamic_viscosity_law_type>( matProps.law( "dynamic-viscosity" ) );
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(u),*dynamicViscosityLawPtr,matProps,se);
    using expr_viscosity_type = std::decay_t<decltype( muExpr )>;
    using expr_coeff_null_type = std::decay_t<decltype( cst(0.) )>;

    // stab gls parameter
    using expr_stab_gls_parameter_type = std::decay_t<decltype( Feel::FeelModels::stabilizationGLSParameterExpr<expr_density_uconv_type,expr_viscosity_type,expr_coeff_null_type>( *this->stabilizationGLSParameterPressure() ) )>;
    std::optional<expr_stab_gls_parameter_type> tauExprSUPG, tauExprPSPG;
    if ( this->stabilizationGLSType() != "pspg" )
        tauExprSUPG.emplace( Feel::FeelModels::stabilizationGLSParameterExpr<expr_density_uconv_type,expr_viscosity_type,expr_coeff_null_type>( *this->stabilizationGLSParameterConvectionDiffusion() ) );
    if ( this->stabilizationGLSType() == "pspg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "gls" )
        tauExprPSPG.emplace( Feel::FeelModels::stabilizationGLSParameterExpr<expr_density_uconv_type,expr_viscosity_type,expr_coeff_null_type>( *this->stabilizationGLSParameterPressure() ) );

    // residual lhs
    using expr_convection_residual_lhs_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*gradt(u)*idv(beta_u) )>;
#if 0
    using expr_viscoous_stress_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*muExpr*laplaciant(u) )>;
#else
    using expr_viscoous_stress_residual_lhs_type = std::decay_t<decltype( -timeSteppingScaling*fluidMecDivViscousStressTensorLinearTrial(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se) )>;
#endif
    using expr_time_derivative_residual_lhs_type = std::decay_t<decltype( expr_density_type{}*this->timeStepBDF()->polyDerivCoefficient(0)*idt(u) )>;
    auto residual_lhs = exprOptionalConcat<expr_convection_residual_lhs_type,
                                           expr_viscoous_stress_residual_lhs_type,
                                           expr_time_derivative_residual_lhs_type>();

    // residual rhs
    using expr_time_derivative_residual_rhs_type = std::decay_t<decltype( expr_density_type{}*idv(this->timeStepBDF()->polyDeriv()))>;
    using expr_previous_residual_type = std::decay_t<decltype( -FluidMechanics_detail::exprResidual( *this, physicFluidData, matProps, mctx, 1.0-timeSteppingScaling, true ) )>;
    auto residual_rhs = exprOptionalConcat<expr_time_derivative_residual_rhs_type,expr_previous_residual_type>();


    if ( physicFluidData.equation() == "Navier-Stokes")
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        residual_lhs.expression().add( timeSteppingScaling*densityExpr*gradt(u)*idv(beta_u) );
        if ( tauExprSUPG )
            tauExprSUPG->expression().setConvection( densityExpr*idv(beta_u) );
        if ( tauExprPSPG )
            tauExprPSPG->expression().setConvection( densityExpr*idv(beta_u) );
    }
    if (!this->isStationaryModel())
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        residual_lhs.expression().add( densityExpr*this->timeStepBDF()->polyDerivCoefficient(0)*idt(u) );
        residual_rhs.expression().add( densityExpr*idv(this->timeStepBDF()->polyDeriv()) );
    }
#if 0
    if constexpr( nOrderVelocity>1 )
    {
        residual_lhs.expression().add( -timeSteppingScaling*muExpr*laplaciant(u) );
    }
#else
    auto divSigmaViscous = fluidMecDivViscousStressTensorLinearTrial(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se,this->stabilizationGLS_checkViscosityDependencyOnCoordinates());
    if ( !divSigmaViscous.expression().isZero() )
        residual_lhs.expression().add( -timeSteppingScaling*divSigmaViscous );
#endif


    if ( !this->isStationaryModel() && this->timeStepping() == "Theta" )
    {
        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
        auto previousResidualExpr = -FluidMechanics_detail::exprResidual( *this, physicFluidData, matProps, mctxPrevious, 1.0-timeSteppingScaling, true );
        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
            previousResidualExpr.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
        residual_rhs.expression().add( previousResidualExpr );
    }

    auto residual_rhs_full = Feel::FeelModels::vfdetail::addExpr( hana::append( exprsAddedInResidualRhsTuple, residual_rhs ) );

    if ( tauExprSUPG )
    {
        tauExprSUPG->expression().setDiffusion( muExpr );

        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );

#if 1
        int coeffNatureStabilization = ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" )? 0 : (this->stabilizationGLSType() == "gls")? 1 : -1;
        // test functions
        using expr_convection_test_type = std::decay_t<decltype( densityExpr*grad(u)*idv(beta_u) )>;
        auto divSigmaViscousTest = fluidMecDivViscousStressTensorLinearTest(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se,this->stabilizationGLS_checkViscosityDependencyOnCoordinates());
        using expr_viscous_stress_test_type = std::decay_t<decltype( -coeffNatureStabilization*divSigmaViscousTest )>;
        auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_viscous_stress_test_type>();
        stab_test.expression().add( densityExpr*grad(u)*idv(beta_u) );
        if ( coeffNatureStabilization != 0 && !divSigmaViscousTest.expression().isZero() )
            stab_test.expression().add( -coeffNatureStabilization*divSigmaViscousTest );
#else
        static constexpr bool velocityOrderGreaterThan1 = nOrderVelocity>1;
        auto stab_test = hana::eval_if( hana::bool_c<velocityOrderGreaterThan1>,
                                        [&u,&beta_u,&densityExpr] { return densityExpr*grad(u)*idv(beta_u); },
                                        [&u,&beta_u,&densityExpr,&muExpr] { return densityExpr*grad(u)*idv(beta_u) - muExpr*laplacian(u); } );
#endif

        auto rangeEltConvectionDiffusion = this->stabilizationGLSEltRangeConvectionDiffusion( matName );
        auto tauFieldPtr = this->stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=*tauExprSUPG);
        auto tau = idv(tauFieldPtr);

        if ( residual_lhs.expression().hasExpr() )
        {
            bilinearFormVV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner( residual_lhs,stab_test ),
                           _geomap=this->geomap() );
        }
        if ( residual_rhs.expression().hasExpr() || decltype(hana::size( exprsAddedInResidualRhsTuple ))::value > 0 )
        {
            myLinearFormV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr= tau*inner( residual_rhs_full,stab_test ),
                           _geomap=this->geomap() );
        }

        auto residual_lhs_p = trans(gradt(p));
        bilinearFormVP +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=tau*inner( residual_lhs_p,stab_test ),
                       _geomap=this->geomap() );

        hana::for_each( exprsAddedInResidualLhsTuple, [this,&XhV,&A,&rangeEltConvectionDiffusion,&tau,&stab_test]( auto & e )
                        {
                            form2( _test=XhV,_trial=std::get<0>(e),_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=std::get<1>(e) ) +=
                                integrate( _range=rangeEltConvectionDiffusion,
                                           _expr=tau*inner(std::get<2>(e),stab_test ),
                                           _geomap=this->geomap() );
                        });

        // divergence stab
        double lambdaCoeff = 1.;
        bilinearFormVV +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=lambdaCoeff*2.0*inner(idv(beta_u))*tau*divt(u)*div(v),
                       _geomap=this->geomap() );
    }

    if ( tauExprPSPG )
    {
        tauExprPSPG->expression().setDiffusion( muExpr );

        auto stab_test = -trans(grad(p));
        auto rangeEltPressure = this->stabilizationGLSEltRangePressure( matName );
        auto tauFieldPtr = this->stabilizationGLSParameterPressure()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltPressure,_expr=*tauExprPSPG);
        auto tau = idv(tauFieldPtr);

        if ( residual_lhs.expression().hasExpr() )
        {
            bilinearFormPV +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(residual_lhs,stab_test ),
                           _geomap=this->geomap() );
        }
        if ( residual_rhs.expression().hasExpr() || decltype(hana::size( exprsAddedInResidualRhsTuple ))::value > 0 )
        {
            myLinearFormP +=
                integrate( _range=rangeEltPressure,
                           _expr= tau*inner(residual_rhs_full,stab_test ),
                           _geomap=this->geomap() );
        }

        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearFormPP +=
            integrate( _range=rangeEltPressure,
                       _expr=tau*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=this->geomap() );

        hana::for_each( exprsAddedInResidualLhsTuple, [this,&XhP,&A,&rangeEltPressure,&tau,&stab_test]( auto & e )
                        {
                            form2( _test=XhP,_trial=std::get<0>(e),_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix()+1,
                                   _colstart=std::get<1>(e) ) +=
                                integrate( _range=rangeEltPressure,
                                           _expr=tau*inner(std::get<2>(e),stab_test ),
                                           _geomap=this->geomap() );
                        });
    }
}

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType,typename RangeType,typename... ExprAddedType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateJacobianStabilizationGLS( 
        DataUpdateJacobian & data, ModelContextType const& mctx,
        ModelPhysicFluid<nDim> const& physicFluidData,
        MaterialProperties const& matProps, RangeType const& range,
        const ExprAddedType&... exprsAddedInResidual ) const
{
    sparse_matrix_ptrtype& J = data.jacobian();

    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    double timeSteppingScaling = 1.;
    if ( !this->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    std::string const& matName = matProps.materialName();
    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto const& u = mctx.field( FieldTag::velocity(this), "velocity" );
    auto const& p = mctx.field( FieldTag::pressure(this), "pressure" );
    auto const& v = u;
    auto const& q = p;
    auto const& beta_u = this->useSemiImplicitTimeScheme()? mctx.field( FieldTag::velocity_extrapolated(this), "velocity_extrapolated" ) : u;
    auto const& se = mctx.symbolsExpr();

    auto bilinearFormVV = form2( _test=XhV,_trial=XhV,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix() );
    auto bilinearFormVP = form2( _test=XhV,_trial=XhP,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix()+1 );
    auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix()+1,
                                 _colstart=this->colStartInMatrix()+0 );
    auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix()+1,
                                 _colstart=this->colStartInMatrix()+1 );


    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;
    auto dynamicViscosityLawPtr = std::static_pointer_cast<dynamic_viscosity_law_type>( matProps.law( "dynamic-viscosity" ) );
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(u), *dynamicViscosityLawPtr, matProps, se);

    // jacobian of residual
    using expr_convection_residual_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*(gradt(u)*idv(u)+gradv(u)*idt(u)) )>;
    using expr_convection_extrapolated_residual_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*(gradt(u)*idv(beta_u)) )>;
    using expr_viscous_stress_residual_type = std::decay_t<decltype( -timeSteppingScaling*fluidMecDivViscousStressTensorJacobian(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se) )>;

    using expr_time_derivative_residual_type = std::decay_t<decltype( expr_density_type{}*this->timeStepBDF()->polyDerivCoefficient(0)*idt(u) )>;
    auto residual_jac = exprOptionalConcat<expr_convection_residual_type,expr_convection_extrapolated_residual_type,
                                           expr_viscous_stress_residual_type,
                                           expr_time_derivative_residual_type>();

    auto additionTermsTuple = hana::make_tuple(exprsAddedInResidual...);

    if ( physicFluidData.equation() == "Navier-Stokes")
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        if ( this->useSemiImplicitTimeScheme() )
            residual_jac.expression().add( timeSteppingScaling*densityExpr*(gradt(u)*idv(beta_u)) );
        else
            residual_jac.expression().add( timeSteppingScaling*densityExpr*(gradt(u)*idv(u)+gradv(u)*idt(u)) );
    }
    if (!this->isStationaryModel())
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        residual_jac.expression().add( densityExpr*this->timeStepBDF()->polyDerivCoefficient(0)*idt(u) );
    }

    auto divSigmaViscous = fluidMecDivViscousStressTensorJacobian(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se,this->stabilizationGLS_checkViscosityDependencyOnCoordinates());
    if ( !divSigmaViscous.expression().isZero() )
        residual_jac.expression().add( -timeSteppingScaling*divSigmaViscous );

    if ( this->stabilizationGLSType() != "pspg" )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
#if 1
        int coeffNatureStabilization = ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" )? 0 : (this->stabilizationGLSType() == "gls")? 1 : -1;
        // test functions
        using expr_convection_test_type = std::decay_t<decltype( densityExpr*grad(u)*idv(beta_u) )>;
        auto divSigmaViscousTest = fluidMecDivViscousStressTensorLinearTest(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se,this->stabilizationGLS_checkViscosityDependencyOnCoordinates());
        using expr_viscous_stress_test_type = std::decay_t<decltype( -coeffNatureStabilization*divSigmaViscousTest )>;
        auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_viscous_stress_test_type>();
        stab_test.expression().add( densityExpr*grad(u)*idv(beta_u) );
        if ( coeffNatureStabilization != 0 && !divSigmaViscousTest.expression().isZero() )
            stab_test.expression().add( -coeffNatureStabilization*divSigmaViscousTest );
#else
        static constexpr bool velocityOrderGreaterThan1 = nOrderVelocity>1;
        auto stab_test = hana::eval_if( hana::bool_c<velocityOrderGreaterThan1>,
                                        [&u,&beta_u,&densityExpr] { return densityExpr*grad(u)*idv(beta_u); },
                                        [&u,&beta_u,&densityExpr,&muExpr] { return densityExpr*grad(u)*idv(beta_u) - muExpr*laplacian(u); } );
#endif
        auto rangeEltConvectionDiffusion = this->stabilizationGLSEltRangeConvectionDiffusion( matName );
        auto tauFieldPtr = this->stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        //tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=*tauExprSUPG);
        auto tau = idv(tauFieldPtr);

        if ( residual_jac.expression().hasExpr() )
        {
            bilinearFormVV +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(residual_jac,stab_test ),
                           _geomap=this->geomap() );
        }
        auto residual_jac_p = trans(gradt(p));
        bilinearFormVP +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=tau*inner(residual_jac_p,stab_test ),
                       _geomap=this->geomap() );

        hana::for_each( additionTermsTuple, [this,&XhV,&J,&rangeEltConvectionDiffusion,&tau,&stab_test]( auto & e )
                        {
                            form2( _test=XhV,_trial=std::get<0>(e),_matrix=J,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=std::get<1>(e) ) +=
                                integrate( _range=rangeEltConvectionDiffusion,
                                           _expr=tau*inner(std::get<2>(e),stab_test ),
                                           _geomap=this->geomap() );
                        });

        // divergence stab
        double lambdaCoeff = 1.;
        bilinearFormVV +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=lambdaCoeff*2.0*inner(idv(u))*tau*divt(u)*div(v),
                       _geomap=this->geomap() );
    }

    if ( this->stabilizationGLSType() == "pspg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "gls" )
    {
        auto stab_test = -trans(grad(p));
        auto rangeEltPressure = this->stabilizationGLSEltRangePressure( matName );
        auto tauFieldPtr = this->stabilizationGLSParameterPressure()->fieldTauPtr();
        //tauFieldPtr->on(_range=rangeEltPressure,_expr=*tauExprPSPG);
        auto tau = idv(tauFieldPtr);

        if ( residual_jac.expression().hasExpr() )
        {
            bilinearFormPV +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner( residual_jac,stab_test ),
                           _geomap=this->geomap() );
        }

        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearFormPP +=
            integrate( _range=rangeEltPressure,
                       _expr=tau*inner( stab_residual_bilinear_p,stab_test ),
                       _geomap=this->geomap() );

        hana::for_each( additionTermsTuple, [this,&XhP,&J,&rangeEltPressure,&tau,&stab_test]( auto & e )
                        {
                            form2( _test=XhP,_trial=std::get<0>(e),_matrix=J,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=this->rowStartInMatrix()+1,
                                   _colstart=std::get<1>(e) ) +=
                                integrate( _range=rangeEltPressure,
                                           _expr=tau*inner( std::get<2>(e),stab_test ),
                                           _geomap=this->geomap() );
                        });
    }

}

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType,typename RangeType,typename... ExprAddedType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateResidualStabilizationGLS( 
        DataUpdateResidual & data, ModelContextType const& mctx,
        ModelPhysicFluid<nDim> const& physicFluidData,
        MaterialProperties const& matProps, RangeType const& range,
        const ExprAddedType&... exprsAddedInResidual ) const
{
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationaryModel() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );
    }
    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
        return;

    std::string const& matName = matProps.materialName();
    auto mesh = this->mesh();
    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    auto const& u = mctx.field( FieldTag::velocity(this), "velocity" );
    auto const& p = mctx.field( FieldTag::pressure(this), "pressure" );
    auto const& v = u;
    auto const& q = p;
    auto const& beta_u = this->useSemiImplicitTimeScheme()? mctx.field( FieldTag::velocity_extrapolated(this), "velocity_extrapolated" ) : u;
    auto const& se = mctx.symbolsExpr();

    auto myLinearFormV = form1( _test=XhV, _vector=R,
                                _rowstart=this->rowStartInVector() );
    auto myLinearFormP = form1( _test=XhP, _vector=R,
                                _rowstart=this->rowStartInVector()+1 );

    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;
    using expr_density_uconv_type = std::decay_t<decltype( expr_density_type{}*idv(beta_u) )>;

    auto dynamicViscosityLawPtr = std::static_pointer_cast<dynamic_viscosity_law_type>( matProps.law( "dynamic-viscosity" ) );
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(u), *dynamicViscosityLawPtr, matProps, se);
    using expr_viscosity_type = std::decay_t<decltype( muExpr )>;
    using expr_coeff_null_type = std::decay_t<decltype( cst(0.) )>;

    // stab gls parameter
    using expr_stab_gls_parameter_type = std::decay_t<decltype( Feel::FeelModels::stabilizationGLSParameterExpr<expr_density_uconv_type,expr_viscosity_type,expr_coeff_null_type>( *this->stabilizationGLSParameterPressure() ) )>;
    std::optional<expr_stab_gls_parameter_type> tauExprSUPG, tauExprPSPG;
    if ( this->stabilizationGLSType() != "pspg" )
        tauExprSUPG.emplace( Feel::FeelModels::stabilizationGLSParameterExpr<expr_density_uconv_type,expr_viscosity_type,expr_coeff_null_type>( *this->stabilizationGLSParameterConvectionDiffusion() ) );
    if ( this->stabilizationGLSType() == "pspg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "gls" )
        tauExprPSPG.emplace( Feel::FeelModels::stabilizationGLSParameterExpr<expr_density_uconv_type,expr_viscosity_type,expr_coeff_null_type>( *this->stabilizationGLSParameterPressure() ) );


    using expr_previous_residual_type = std::decay_t<decltype( FluidMechanics_detail::exprResidual( *this, physicFluidData, matProps, mctx, 1.0-timeSteppingScaling, true ) )>;
    auto residual_base = FluidMechanics_detail::exprResidual<expr_previous_residual_type>( *this, physicFluidData, matProps, mctx, timeSteppingScaling, timeSteppingEvaluateResidualWithoutTimeDerivative );
    if ( !this->isStationaryModel() && this->timeStepping() == "Theta" )
    {
        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
        auto previousResidualExpr = FluidMechanics_detail::exprResidual( *this, physicFluidData, matProps, mctxPrevious, 1.0-timeSteppingScaling, true );
        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
            previousResidualExpr.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
        residual_base.expression().add( previousResidualExpr );
    }

    auto residual_full = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( residual_base, exprsAddedInResidual... ) );

    if ( physicFluidData.equation() == "Navier-Stokes")
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        if ( tauExprSUPG )
            tauExprSUPG->expression().setConvection( densityExpr*idv(beta_u) );
        if ( tauExprPSPG )
            tauExprPSPG->expression().setConvection( densityExpr*idv(beta_u) );
    }
    if ( tauExprSUPG )
    {
        tauExprSUPG->expression().setDiffusion( muExpr );

        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
#if 1
        int coeffNatureStabilization = ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" )? 0 : (this->stabilizationGLSType() == "gls")? 1 : -1;
        // test functions
        using expr_convection_test_type = std::decay_t<decltype( densityExpr*grad(u)*idv(beta_u) )>;
        auto divSigmaViscousTest = fluidMecDivViscousStressTensorLinearTest(u,physicFluidData,matProps,this->worldComm(),this->repository().expr(),se,this->stabilizationGLS_checkViscosityDependencyOnCoordinates());
        using expr_viscous_stress_test_type = std::decay_t<decltype( -coeffNatureStabilization*divSigmaViscousTest )>;
        auto stab_test = exprOptionalConcat<expr_convection_test_type,expr_viscous_stress_test_type>();
        stab_test.expression().add( densityExpr*grad(u)*idv(beta_u) );
        if ( coeffNatureStabilization != 0 && !divSigmaViscousTest.expression().isZero() )
            stab_test.expression().add( -coeffNatureStabilization*divSigmaViscousTest );
#else
        static constexpr bool velocityOrderGreaterThan1 = nOrderVelocity>1;
        auto stab_test = hana::eval_if( hana::bool_c<velocityOrderGreaterThan1>,
                                        [&u,&beta_u,&densityExpr] { return densityExpr*grad(u)*idv(beta_u); },
                                        [&u,&beta_u,&densityExpr,&muExpr] { return densityExpr*grad(u)*idv(beta_u) - muExpr*laplacian(u); } );
#endif
        auto rangeEltConvectionDiffusion = this->stabilizationGLSEltRangeConvectionDiffusion( matName );
        auto tauFieldPtr = this->stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=*tauExprSUPG);
        auto tau = idv(tauFieldPtr);

        myLinearFormV +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=tau*inner( residual_full,stab_test ),
                       _geomap=this->geomap() );

        // divergence stab
        double lambdaCoeff = 1.;
        myLinearFormV +=
            integrate( _range=rangeEltConvectionDiffusion,
                       _expr=lambdaCoeff*2.0*inner(idv(u))*tau*divv(u)*div(v),
                       _geomap=this->geomap() );
    }

    if ( tauExprPSPG )
    {
        tauExprPSPG->expression().setDiffusion( muExpr );

        auto stab_test = -trans(grad(p));
        auto rangeEltPressure = this->stabilizationGLSEltRangePressure( matName );
        auto tauFieldPtr = this->stabilizationGLSParameterPressure()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltPressure,_expr=*tauExprPSPG);
        auto tau = idv(tauFieldPtr);

        myLinearFormP +=
            integrate( _range=rangeEltPressure,
                       _expr=tau*inner(residual_full,stab_test ),
                       _geomap=this->geomap() );
    }

}

} // namespace Feel
} // namespace FeelModels

#endif
