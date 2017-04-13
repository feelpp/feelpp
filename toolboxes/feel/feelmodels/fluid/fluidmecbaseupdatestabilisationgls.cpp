
#include <feel/feelmodels/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{

namespace FluidMechanicsDetail
{
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



template< int StabParamType,int StabGLSType,typename FluidMechanicsType>
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

    auto const& rho = fluidmec.densityViscosityModel()->fieldRho();
    auto myViscosity = Feel::vf::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(betaU,p,*fluidmec.densityViscosityModel());
    //auto myViscosity = idv(fluidmec.densityViscosityModel()->fieldMu());

    auto uconv = idv(rho)*idv(betaU);
    auto tau = fluidmec.stabilizationGLSParameter()->tau( uconv,myViscosity, mpl::int_<StabParamType>() );
    if ( fluidmec.stabilizationGLSType() != "pspg" )
    {
        //auto stab_test = grad(u)*uconv - myViscosity*laplacian(u);
        auto stab_test = stabGLStestLinearExpr_u( idv(rho),myViscosity,idv(betaU),u, mpl::int_<StabGLSType>() );
        if (!fluidmec.isStationary())
        {
            //auto stab_residual_bilinear_u = idv(rho)*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + gradt(u)*idv(betaU) ) - myViscosity*laplaciant(u);
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr=val(tau)*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeStepVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeStep = rhsTimeStepVP.template element<0>();
            auto stab_residual_linear = idv(rho)*idv(rhsTimeStep);
            myLinearForm +=
                integrate( _range=elements(mesh),
                           _expr= val(tau)*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            //auto stab_residual_bilinear_u = idv(rho)*gradt(u)*idv(betaU) - myViscosity*laplaciant(u);
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr=val(tau)*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr=val(tau)*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );
    }

    if ( fluidmec.stabilizationGLSType() == "pspg" || fluidmec.stabilizationGLSType() == "supg-pspg" || fluidmec.stabilizationGLSType() == "gls" )
    {
        auto stab_test = -idv(rho)*trans(grad(p));
        if ( !fluidmec.isStationary() )
        {
            //auto stab_residual_bilinear_u = idv(rho)*(idt(u)*fluidmec.timeStepBDF()->polyDerivCoefficient(0) + gradt(u)*idv(betaU) ) - myViscosity*laplaciant(u);
            auto stab_residual_bilinear_u = residualTransientLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr=val(tau)*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
            auto rhsTimeStepVP = fluidmec.timeStepBDF()->polyDeriv();
            auto rhsTimeStep = rhsTimeStepVP.template element<0>();
            auto stab_residual_linear = idv(rho)*idv(rhsTimeStep);
            myLinearForm +=
                integrate( _range=elements(mesh),
                           _expr= val(tau)*inner(stab_residual_linear,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        else
        {
            //auto stab_residual_bilinear_u = idv(rho)*gradt(u)*idv(betaU) - myViscosity*laplaciant(u);
            auto stab_residual_bilinear_u = residualStationaryLinearExpr_u( idv(rho),myViscosity,idv(betaU),u,fluidmec, mpl::int_<StabResidualType>() );
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr=val(tau)*inner(stab_residual_bilinear_u,stab_test ),
                           _geomap=fluidmec.geomap() );
        }
        auto stab_residual_bilinear_p = trans(gradt(p));
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr=val(tau)*inner(stab_residual_bilinear_p,stab_test ),
                       _geomap=fluidmec.geomap() );
    }
    fluidmec.log("FluidMechanics","updateLinearPDEStabilizationGLS", "finish" );
}

template< int StabParamType,typename FluidMechanicsType>
void
updateLinearPDEStabilizationGLS( FluidMechanicsType const& fluidmec, ModelAlgebraic::DataUpdateLinear & data )
{
    if ( ( fluidmec.stabilizationGLSType() == "supg" || fluidmec.stabilizationGLSType() == "supg-pspg" || fluidmec.stabilizationGLSType() == "pspg" ) ||
         ( FluidMechanicsType::nOrderVelocity<=1 && ( fluidmec.stabilizationGLSType() == "gls" || fluidmec.stabilizationGLSType() == "gls-no-pspg" ) ) )
    {
        updateLinearPDEStabilizationGLS<StabParamType,0>( fluidmec, data );
    }
    else if ( fluidmec.stabilizationGLSType() == "gls" || fluidmec.stabilizationGLSType() == "gls-no-pspg" )
    {
        updateLinearPDEStabilizationGLS<StabParamType,1>( fluidmec, data );
    }
}

} // namespace FluidMechanicsDetail

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilisationGLS( DataUpdateLinear & data ) const
{
    // update stabilization gls
    if ( M_stabilizationGLS )
    {
        if ( M_stabilizationGLSParameter->method() == "eigenvalue" )
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<0>( *this, data );
        else if ( M_stabilizationGLSParameter->method() == "doubly-asymptotic-approximation" )
            FluidMechanicsDetail::updateLinearPDEStabilizationGLS<2>( *this, data );
    }
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateResidualStabilisationGLS( DataUpdateResidual & data, element_fluid_external_storage_type const& U ) const
{}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateJacobianStabilisationGLS( DataUpdateJacobian & data, element_fluid_external_storage_type const& U ) const
{}


} // end namespace FeelModels
} // end namespace Feel



