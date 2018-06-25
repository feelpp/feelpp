
#include <feel/feelmodels/fluid/fluidmechanics.hpp>

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

template< int StabGLSType,typename FluidMechanicsType>
void
updateLinearPDEStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateLinear & data )
{
    static const int StabResidualType = ( FluidMechanicsType::nOrderVelocity>1 )? 1 : 0;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLS", "start" );
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();

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

    //CHECK( fluidmec.solverName() == "Oseen" || fluidmec.solverName() == "Picard" ) << "invalid solver name " << fluidmec.solverName();

    typename FluidMechanicsType::element_fluid_ptrtype fielCurrentPicardSolution;
    if ( fluidmec.solverName() == "Picard" )
    {
        fielCurrentPicardSolution = Xh->elementPtr();
        *fielCurrentPicardSolution = *Xh->elementPtr(*vecCurrentPicardSolution, fluidmec.rowStartInVector() );
    }

    auto const& BetaU = ( fluidmec.solverName() == "Oseen" )? fluidmec.timeStepBDF()->poly() : *fielCurrentPicardSolution;
    auto betaU = BetaU.template element<0>();

    auto const& rho = fluidmec.materialProperties()->fieldRho();
    //auto myViscosity = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(betaU,p,*fluidmec.materialProperties());
    auto myViscosity = idv(fluidmec.materialProperties()->fieldMu());

    auto uconv = idv(rho)*idv(betaU);
#if 0
#if 1
    auto tau = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameter(), uconv, myViscosity );
#else
    static const uint16_type nStabGlsOrderPoly = (FluidMechanicsType::nOrderVelocity>1)? FluidMechanicsType::nOrderVelocity : 2;
    typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, nStabGlsOrderPoly> stab_gls_parameter_impl_type;
    auto tau =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameter() )->tau( uconv, myViscosity, mpl::int_<0/*StabParamType*/>() );
    std::cout << "decltype(tau)::imorder=" << decltype(tau)::imorder << "\n";
    std::cout << "decltype(tau)::imIsPoly=" << decltype(tau)::imIsPoly << "\n";
#endif
#endif
    bool hasUpdatedTauForConvectionDiffusion = false;
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        auto rangeEltConvectionDiffusion = fluidmec.stabilizationGLSEltRangeConvectionDiffusion();
#if 0
        typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderVelocity> stab_gls_parameter_impl_type;
        auto stabGLSParamConvectionDiffusion =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterConvectionDiffusion() );
        stabGLSParamConvectionDiffusion->updateTau(uconv, myViscosity, rangeEltConvectionDiffusion);
        hasUpdatedTauForConvectionDiffusion = true;
        auto tau = idv(stabGLSParamConvectionDiffusion->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterConvectionDiffusion(), uconv, myViscosity );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterConvectionDiffusion()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltConvectionDiffusion,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif
        auto stab_test = stabGLStestLinearExpr_u( idv(rho),myViscosity,idv(betaU),u, mpl::int_<StabGLSType>() );
        if (!fluidmec.isStationaryModel())
        {
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto stab_residual_linear = idv(rho)*idv(rhsTimeDerivative);
            myLinearForm +=
                integrate( _range=rangeEltConvectionDiffusion,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
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
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( marker(d).empty() )? rangeEltConvectionDiffusion : intersect(markedelements(mesh,marker(d)),rangeEltConvectionDiffusion);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=tau*inner(expression(d),stab_test),
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
            stabGLSParamPressure->updateTau(uconv, myViscosity, rangeEltPressure);
        auto tau = idv(stabGLSParamPressure->fieldTau());
#else
        auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterPressure(), uconv, myViscosity );
        auto tauFieldPtr = fluidmec.stabilizationGLSParameterPressure()->fieldTauPtr();
        tauFieldPtr->on(_range=rangeEltPressure,_expr=tauExpr);
        auto tau = idv(tauFieldPtr);
#endif

        auto stab_test = -idv(rho)*trans(grad(p));
        if ( !fluidmec.isStationaryModel() )
        {
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
            auto stab_residual_linear = idv(rho)*idv(rhsTimeDerivative);
            myLinearForm +=
                integrate( _range=rangeEltPressure,
                           _expr= tau*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
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
        // bodyforces
        for( auto const& d : fluidmec.bodyForces() )
        {
            auto rangeBodyForceUsed = ( marker(d).empty() )? rangeEltPressure : intersect(markedelements(mesh,marker(d)),rangeEltPressure);
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=tau*inner(expression(d),stab_test),
                           _geomap=fluidmec.geomap() );
        }
    }
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLS", "finish" );
}

template< typename FluidMechanicsType>
void
updateLinearPDEStabilizationGLSStokes( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateLinear & data )
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    if ( !buildCstPart ) // stab term are constant if viscosity/density are constant
        return;
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLSStokes", "start" );

    auto mesh = fluidmec.mesh();
    auto Xh = fluidmec.functionSpace();
    auto const& U = fluidmec.fieldVelocityPressure();
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=fluidmec.rowStartInMatrix(),
                                              _colstart=fluidmec.colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=fluidmec.rowStartInVector() );


    auto rho = idv(fluidmec.materialProperties()->fieldRho());
    //auto mu = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
    auto mu = idv(fluidmec.materialProperties()->fieldMu());
    auto uconv = zero<FluidMechanicsType::nRealDim,1>();
    //auto tau = Feel::vf::FeelModels::stabilizationGLSParameterExpr<false>( *fluidmec.stabilizationGLSParameter(), uconv, mu );
    auto rangeEltPressure = fluidmec.stabilizationGLSEltRangePressure();

#if 0
    typedef StabilizationGLSParameter<typename FluidMechanicsType::mesh_type, FluidMechanicsType::nOrderPressure> stab_gls_parameter_impl_type;
    auto stabGLSParamPressure =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( fluidmec.stabilizationGLSParameterPressure() );
    stabGLSParamPressure->template updateTau<false>(uconv, mu, rangeEltPressure);
    auto tau = idv(stabGLSParamPressure->fieldTau());
#else
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr<false>( *fluidmec.stabilizationGLSParameterPressure(), uconv, mu );
    auto tauFieldPtr = fluidmec.stabilizationGLSParameterPressure()->fieldTauPtr();
    tauFieldPtr->on(_range=rangeEltPressure,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);
#endif

    auto stab_test = -rho*trans(grad(p));
    if ( !fluidmec.isStationaryModel() )
    {
        if ( FluidMechanicsType::nOrderVelocity>1 )
        {
            auto stab_residual_bilinear_u = rho*idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) - mu*laplaciant(u);
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            auto stab_residual_bilinear_u = rho*idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0);
            bilinearForm_PatternCoupled +=
                integrate( _range=rangeEltPressure,
                           _expr=tau*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto rhsTimeDerivativeVP = fluidmec.timeStepBDF()->polyDeriv();
        auto rhsTimeDerivative = rhsTimeDerivativeVP.template element<0>();
        auto stab_residual_linear = rho*idv(rhsTimeDerivative);
        myLinearForm +=
            integrate( _range=rangeEltPressure,
                       _expr= tau*inner(stab_residual_linear,stab_test ),
                       _geomap=fluidmec.geomap() );
    }
    else if ( FluidMechanicsType::nOrderVelocity>1 )
    {
        auto stab_residual_bilinear_u =  -mu*laplaciant(u);
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
    // bodyforces
    for( auto const& d : fluidmec.bodyForces() )
    {
        auto rangeBodyForceUsed = ( marker(d).empty() )? rangeEltPressure : intersect(markedelements(mesh,marker(d)),rangeEltPressure);
        myLinearForm +=
            integrate( _range=rangeBodyForceUsed,
                       _expr=tau*inner(expression(d),stab_test),
                       _geomap=fluidmec.geomap() );
    }

    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLSStokes", "finish" );
}



} // namespace FluidMechanicsDetail

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilisationGLS( DataUpdateLinear & data ) const
{
    if ( M_stabilizationGLS )
    {
        if ( this->modelName() == "Navier-Stokes")
        {
            if ( ( this->stabilizationGLSType() == "supg" || this->stabilizationGLSType() == "supg-pspg" || this->stabilizationGLSType() == "pspg" ) ||
                 ( nOrderVelocity<=1 && ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" ) ) )
            {
                FluidMechanicsDetail::updateLinearPDEStabilizationGLS<0>( *this, data );
            }
            else if ( this->stabilizationGLSType() == "gls" || this->stabilizationGLSType() == "gls-no-pspg" )
            {
                FluidMechanicsDetail::updateLinearPDEStabilizationGLS<1>( *this, data );
            }
        }
        else if ( this->modelName() == "Stokes" || this->modelName() == "StokesTransient" )
        {
            CHECK( this->stabilizationGLSType() == "pspg" ) << "only pspg stab for Stokes";
            FluidMechanicsDetail::updateLinearPDEStabilizationGLSStokes( *this, data );
        }
    }
}




} // end namespace FeelModels
} // end namespace Feel



