/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_TOOLBOXES_HEAT_UPDATESTABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_HEAT_UPDATESTABILIZATIONGLS_HPP 1

#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feeldiscr/pchv.hpp>
//#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelmesh/intersect.hpp>

namespace Feel
{
namespace FeelModels
{
template< typename ConvexType, typename BasisTemperatureType >
template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateLinearPDEStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                                                        Expr<ConductivityExprType> const& kappa,
                                                                        Expr<ConvectionExprType> const& uconv,
                                                                        RangeType const& range,
                                                                        DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateLinearPDEStabilizationGLS", "start"+sc);

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto const& v = this->fieldTemperature();

    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector() );



    auto rhocpuconv = rhocp*uconv;
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),rhocpuconv/*rhocp*uconv*/, kappa );
    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
    tauFieldPtr->on(_range=range,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);

    if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg"  )
    {
        auto stab_test = rhocp*grad(u)*uconv;
        if (!this->isStationary())
        {
            auto stab_residual_bilinear = rhocp*(idt(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(u)*uconv );
            bilinearForm_PatternCoupled +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
            auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
            auto stab_residual_linear = rhocp*idv( rhsTimeStep );
            myLinearForm +=
                integrate( _range=range,
                           _expr= val(tau*stab_residual_linear)*stab_test,
                           _geomap=this->geomap() );
        }
        else
        {
            auto stab_residual_bilinear = rhocp*gradt(u)*uconv;
            bilinearForm_PatternCoupled +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? range : intersect( markedelements(mesh,markers(d)),range );
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=tau*expression(d,this->symbolsExpr())*stab_test,
                           _geomap=this->geomap() );
        }
    }
    else if ( ( this->stabilizationGLSType() == "gls" ) || ( this->stabilizationGLSType() == "unusual-gls" ) )
    {
        int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? -1 : 1;
        auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
        if (!this->isStationary())
        {
            auto stab_residual_bilinear = rhocp*(idt(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(u)*uconv ) - kappa*laplaciant(u);
            bilinearForm_PatternCoupled +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
            auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
            auto stab_residual_linear = rhocp*idv( rhsTimeStep );
            myLinearForm +=
                integrate( _range=range,
                           _expr= val(tau*stab_residual_linear)*stab_test,
                           //_expr= tau*stab_residual_linear*stab_test,
                           _geomap=this->geomap() );
        }
        else
        {
            auto stab_residual_bilinear = rhocp*gradt(u)*uconv - kappa*laplaciant(u);
            bilinearForm_PatternCoupled +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? range : intersect( markedelements(mesh,markers(d)), range );
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=tau*expression(d,this->symbolsExpr())*stab_test,
                           _geomap=this->geomap() );
        }
    }

    this->log("Heat","updateLinearPDEStabilizationGLS", "finish"+sc);
}

template< typename ConvexType, typename BasisTemperatureType >
template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateJacobianStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                                                       Expr<ConductivityExprType> const& kappa,
                                                                       Expr<ConvectionExprType> const& uconv,
                                                                       RangeType const& range,
                                                                       DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateJacobianStabilizationGLS", "start"+sc);

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto u = Xh->element( XVec, this->rowStartInVector() );
    auto const& v = this->fieldTemperature();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
#if 0
    auto rhocpuconv = rhocp*uconv;
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),rhocpuconv/*rhocp*uconv*/, kappa );
    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
    tauFieldPtr->on(_range=range,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);
#else
    auto tau = idv( this->stabilizationGLSParameter()->fieldTauPtr() );
#endif

    if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg"  )
    {
        auto stab_test = rhocp*grad(u)*uconv;
        if (!this->isStationary())
        {
            auto stab_residual_bilinear = rhocp*(idt(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(u)*uconv );
            bilinearForm +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
        }
        else
        {
            auto stab_residual_bilinear = rhocp*gradt(u)*uconv;
            bilinearForm +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
        }
    }
    else if ( ( this->stabilizationGLSType() == "gls" ) || ( this->stabilizationGLSType() == "unusual-gls" ) )
    {
        int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? -1 : 1;
        auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
        if (!this->isStationary())
        {
            auto stab_residual_bilinear = rhocp*(idt(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(u)*uconv ) - kappa*laplaciant(u);
            bilinearForm +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
        }
        else
        {
            auto stab_residual_bilinear = rhocp*gradt(u)*uconv - kappa*laplaciant(u);
            bilinearForm +=
                integrate( _range=range,
                           _expr=tau*stab_residual_bilinear*stab_test,
                           _geomap=this->geomap() );
        }
    }

    this->log("Heat","updateJacobianStabilizationGLS", "finish"+sc);

}

template< typename ConvexType, typename BasisTemperatureType >
template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateResidualStabilizationGLS( Expr<RhoCpExprType> const& rhocp,
                                                                       Expr<ConductivityExprType> const& kappa,
                                                                       Expr<ConvectionExprType> const& uconv,
                                                                       RangeType const& range,
                                                                       DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    //bool useJacobianLinearTerms = data.useJacobianLinearTerms();
    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateResidualStabilizationGLS", "start"+sc);

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto u = Xh->element( XVec, this->rowStartInVector() );
    auto const& v = this->fieldTemperature();

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() );
#if 1
    auto rhocpuconv = rhocp*uconv;
    auto tauExpr = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),rhocpuconv/*rhocp*uconv*/, kappa );
    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
    tauFieldPtr->on(_range=range,_expr=tauExpr);
    auto tau = idv(tauFieldPtr);
#else
    auto tau = idv( this->stabilizationGLSParameter()->fieldTauPtr() );
#endif

    if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg"  )
    {
        auto stab_test = rhocp*grad(u)*uconv;
        if (!this->isStationary())
        {
            auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
            auto stab_residual = rhocp*(idv(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) - idv( rhsTimeStep ) + gradv(u)*uconv );
            myLinearForm +=
                integrate( _range=range,
                           _expr=val(tau*stab_residual)*stab_test,
                           _geomap=this->geomap() );
        }
        else
        {
            auto stab_residual = rhocp*gradv(u)*uconv;
            myLinearForm +=
                integrate( _range=range,
                           _expr=val(tau*stab_residual)*stab_test,
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? range : intersect( markedelements(mesh,markers(d)), range );
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-tau*expression(d,this->symbolsExpr())*stab_test,
                           _geomap=this->geomap() );
        }
    }
    else if ( ( this->stabilizationGLSType() == "gls" ) || ( this->stabilizationGLSType() == "unusual-gls" ) )
    {
        int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? -1 : 1;
        auto stab_test = rhocp*grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
        if (!this->isStationary())
        {
            auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
            auto stab_residual = rhocp*(idv(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) - idv( rhsTimeStep ) + gradv(u)*uconv ) - kappa*laplacianv(u);
            myLinearForm +=
                integrate( _range=range,
                           _expr=val(tau*stab_residual)*stab_test,
                           _geomap=this->geomap() );
        }
        else
        {
            auto stab_residual = rhocp*gradv(u)*uconv - kappa*laplacianv(u);
            myLinearForm +=
                integrate( _range=range,
                           _expr=val(tau*stab_residual)*stab_test,
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->bodyForces() )
        {
            auto rangeBodyForceUsed = ( markers(d).empty() )? range : intersect( markedelements(mesh,markers(d)), range );
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr=-tau*expression(d,this->symbolsExpr())*stab_test,
                           _geomap=this->geomap() );
        }
    }

    this->log("Heat","updateResidualStabilizationGLS", "finish"+sc);

}

} // namespace Feel
} // namespace FeelModels

#endif
