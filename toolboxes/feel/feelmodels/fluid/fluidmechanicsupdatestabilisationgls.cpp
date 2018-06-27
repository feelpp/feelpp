
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
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *fluidmec.stabilizationGLSParameterPressure(), uconv, mu, false );
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





} // end namespace FeelModels
} // end namespace Feel



