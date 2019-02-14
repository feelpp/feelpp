/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/heat/heat.hpp>

#include <feel/feelvf/vf.hpp>

//#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>

namespace Feel
{
namespace FeelModels
{

#if 0
HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilizationGLS( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    if ( !this->fieldVelocityConvectionIsUsedAndOperational() || buildCstPart )
        return;

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

    //if ( this->fieldVelocityConvectionIsUsedAndOperational() && !buildCstPart )
    //{
        auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
        auto thecoeff = idv(this->thermalProperties()->fieldRho())*idv(this->thermalProperties()->fieldHeatCapacity());
        auto uconv=thecoeff*idv(this->fieldVelocityConvection());
        auto rangeStabUsed = M_rangeMeshElements;
#if 0
        typedef StabilizationGLSParameter<mesh_type, nOrderTemperature> stab_gls_parameter_impl_type;
        auto stabGLSParam =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( this->stabilizationGLSParameter() );
        stabGLSParam->updateTau(uconv, kappa, rangeStabUsed);
        auto tau = idv(stabGLSParam->fieldTau());
#elif 1
        auto tau = Feel::FeelModels::stabilizationGLSParameterExpr( *this->stabilizationGLSParameter(),uconv, kappa );
#else
        static const uint16_type nStabGlsOrderPoly = (nOrderTemperature>1)? nOrderTemperature : 2;
        typedef StabilizationGLSParameter<mesh_type, nStabGlsOrderPoly> stab_gls_parameter_impl_type;
        auto tau =  std::dynamic_pointer_cast<stab_gls_parameter_impl_type>( this->stabilizationGLSParameter() )->tau( uconv, kappa, mpl::int_<0/*StabParamType*/>() );
        std::cout << "decltype(tau)::imorder=" << decltype(tau)::imorder << "\n";
        std::cout << "decltype(tau)::imIsPoly=" << decltype(tau)::imIsPoly << "\n";

#endif
        if ( nOrderTemperature <= 1 || this->stabilizationGLSType() == "supg"  )
        {
            auto stab_test = grad(u)*uconv;
            if (!this->isStationary())
            {
                auto stab_residual_bilinear = thecoeff*(idt(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(u)*idv(this->fieldVelocityConvection()) );
                bilinearForm_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr=tau*stab_residual_bilinear*stab_test,
                               _geomap=this->geomap() );
                auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                auto stab_residual_linear = thecoeff*idv( rhsTimeStep );
                myLinearForm +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= val(tau*stab_residual_linear)*stab_test,
                               _geomap=this->geomap() );
            }
            else
            {
                auto stab_residual_bilinear = thecoeff*gradt(u)*idv(this->fieldVelocityConvection());
                bilinearForm_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr=tau*stab_residual_bilinear*stab_test,
                               _geomap=this->geomap() );
            }
            for( auto const& d : this->bodyForces() )
            {
                auto theExpr = expr( expression(d),this->symbolsExpr() );
                auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(mesh,markers(d));
                myLinearForm +=
                    integrate( _range=rangeBodyForceUsed,
                               _expr=tau*theExpr*stab_test,
                               _geomap=this->geomap() );
            }
        }
        else if ( ( this->stabilizationGLSType() == "gls" ) || ( this->stabilizationGLSType() == "unusual-gls" ) )
        {
            int stabCoeffDiffusion = (this->stabilizationGLSType() == "gls")? 1 : -1;
            auto stab_test = grad(u)*uconv + stabCoeffDiffusion*kappa*laplacian(u);
            if (!this->isStationary())
            {
                auto stab_residual_bilinear = thecoeff*(idt(u)*this->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(u)*idv(this->fieldVelocityConvection()) ) - kappa*laplaciant(u);
                bilinearForm_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr=tau*stab_residual_bilinear*stab_test,
                               _geomap=this->geomap() );
                auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                auto stab_residual_linear = thecoeff*idv( rhsTimeStep );
                myLinearForm +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= val(tau*stab_residual_linear)*stab_test,
                               _geomap=this->geomap() );
            }
            else
            {
                auto stab_residual_bilinear = thecoeff*gradt(u)*idv(this->fieldVelocityConvection()) - kappa*laplaciant(u);
                bilinearForm_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr=tau*stab_residual_bilinear*stab_test,
                               _geomap=this->geomap() );
            }
            for( auto const& d : this->bodyForces() )
            {
                auto theExpr = expr( expression(d),this->symbolsExpr() );
                auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(mesh,markers(d));
                myLinearForm +=
                    integrate( _range=rangeBodyForceUsed,
                               _expr=tau*theExpr*stab_test,
                               _geomap=this->geomap() );
            }
        }
        //}

}
#endif





HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateLinearPDE", "start"+sc);
    this->timerTool("Solve").start();

    bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
    bool BuildNonCstPart_Form1TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        BuildNonCstPart_Form2TransientTerm = buildCstPart;

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

    //--------------------------------------------------------------------------------------------------//

    for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isMatrix() )
        {
            if ( buildCstPart )
            {
                auto kappa = thermalConductivity.template exprMatrix<nDim,nDim>();
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               //_expr= inner(trans(kappa*trans(gradt(u))),grad(v)),
                               _expr= grad(v)*(kappa*trans(gradt(u))),
                               _geomap=this->geomap() );
            }
        }
        else if ( thermalConductivity.isConstant() )
        {
            if ( buildCstPart )
            {
                double kappa = thermalConductivity.value();
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= kappa*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            if ( buildCstPart )
            {
                auto kappa = thermalConductivity.expr();
                //auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= kappa*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }

        auto const& rhoHeatCapacity = this->thermalProperties()->rhoHeatCapacity( matName );
        if ( rhoHeatCapacity.isConstant() )
        {
            double rhoHeatCapacityValue = rhoHeatCapacity.value();
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityValue*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
            }

            if ( !this->isStationary() )
            {
                if ( BuildNonCstPart_Form2TransientTerm )
                {
                    double thecoeff = rhoHeatCapacityValue*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= thecoeff*idt(u)*id(v),
                                   _geomap=this->geomap() );
                }
                if ( BuildNonCstPart_Form1TransientTerm )
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= rhoHeatCapacityValue*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }
        else
        {
            auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
            //auto rhoHeatCapacityExpr = idv(this->thermalProperties()->fieldRho())*idv(this->thermalProperties()->fieldHeatCapacity());
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityExpr*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
            }

            if ( !this->isStationary() )
            {
                if ( BuildNonCstPart_Form2TransientTerm )
                {
                    auto thecoeff = rhoHeatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= thecoeff*idt(u)*id(v),
                                   _geomap=this->geomap() );
                }
                if ( BuildNonCstPart_Form1TransientTerm )
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= rhoHeatCapacityExpr*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }


        // update stabilization gls
        if ( M_stabilizationGLS && buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
        {
            CHECK( !thermalConductivity.isMatrix() && thermalConductivity.isConstant() && rhoHeatCapacity.isConstant() ) << "NotImplemented";
            this->updateLinearPDEStabilizationGLS(  cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(this->fieldVelocityConvection()),range,data );
        }

    }

#if 0
    if ( this->fieldVelocityConvectionIsUsedAndOperational() && !buildCstPart )
    {
        //viscous dissipation
        if ( false/*true*/ )
        {
            double mu = 1.;
            auto defv = sym(gradv( this->fieldVelocityConvection() ) );
            auto defv2 = inner(defv,defv);

            if ( !this->fieldVelocityConvectionIsIncompressible() )
            {
#if 0
                bilinearForm_PatternCoupled +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= thecoeff*(idt(u)*divv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
#endif
                myLinearForm +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= 2*mu*defv2*id(v),
                               _geomap=this->geomap() );
            }
            else
            {
                auto incomp2 = pow( divv( this->fieldVelocityConvection() ),2 );
                myLinearForm +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= 2*mu*(defv2-(1./3.)*incomp2)*id(v),
                               _geomap=this->geomap() );
            }
        }
    }
#endif

    //--------------------------------------------------------------------------------------------------//

    // update source term
    this->updateLinearPDESourceTerm( F, buildCstPart );

    // update bc
    this->updateLinearPDEWeakBC( A,F,buildCstPart );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Heat","updateLinearPDE",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("Heat","updateNewtonInitialGuess","start" );

    auto mesh = this->mesh();
    auto u = this->spaceTemperature()->element( U, this->rowStartInVector() );

    for( auto const& d : M_bcDirichlet )
    {
        auto theExpr = expr( expression(d),this->symbolsExpr() );
        u.on(_range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
             _expr=theExpr );
    }
    // synchronize temperature dof on interprocess
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("temperature");
    if ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )
        sync( u, "=", itFindDofsWithValueImposed->second );

    this->log("Heat","updateNewtonInitialGuess","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool _BuildCstPart = data.buildCstPart();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateJacobian", "start"+sc);

    bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        BuildNonCstPart_Form2TransientTerm = buildCstPart;


    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    //auto const& u = this->fieldTemperature();
    auto const u = Xh->element(XVec, this->rowStartInVector());
    auto const& v = this->fieldTemperature();


    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isConstant() )
        {
            if ( buildCstPart )
            {
                double kappa = thermalConductivity.value();
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= kappa*inner(gradt(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto kappa = thermalConductivity.expr();
            std::string symbolStr = "heat_T";
            if ( kappa.expression().hasSymbol( symbolStr ) )
            {
                if ( buildNonCstPart )
                {
                    auto kappaEval = expr( kappa, symbolExpr(symbolStr, idv(u)) );
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= kappaEval*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                    auto kappaDiff = diff( kappa,symbolStr,1,"",this->worldComm(),this->repository().expr());
                    auto kappaDiffEval = expr( kappaDiff, symbolExpr( symbolStr, idv(u) ) );
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= kappaDiffEval*idt(u)*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                if ( buildCstPart )
                {
                    //auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= kappa*inner(gradt(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }

        auto const& rhoHeatCapacity = this->thermalProperties()->rhoHeatCapacity( matName );
        if ( rhoHeatCapacity.isConstant() )
        {
            double rhoHeatCapacityValue = rhoHeatCapacity.value();
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityValue*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
            }

            if ( !this->isStationary() )
            {
                if ( BuildNonCstPart_Form2TransientTerm )
                {
                    double thecoeff = rhoHeatCapacityValue*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= thecoeff*idt(u)*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }
        else
        {
            auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
            //auto rhoHeatCapacityExpr = idv(this->thermalProperties()->fieldRho())*idv(this->thermalProperties()->fieldHeatCapacity());

            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityExpr*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
            }

            if ( !this->isStationary() )
            {
                if ( BuildNonCstPart_Form2TransientTerm )
                {
                    auto thecoeff = rhoHeatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= thecoeff*idt(u)*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }

        // update stabilization gls
        if ( M_stabilizationGLS && buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
        {
            CHECK( !thermalConductivity.isMatrix() && thermalConductivity.isConstant() && rhoHeatCapacity.isConstant() ) << "NotImplemented";
            this->updateJacobianStabilizationGLS(  cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(this->fieldVelocityConvection()),range,data );
        }

    }

    this->updateJacobianRobinBC( J,buildCstPart );
}
HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _BuildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildNonCstPart = !_BuildCstPart;
    bool buildCstPart = _BuildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateResidual", "start"+sc);

    bool Build_TransientTerm = buildNonCstPart;
    if ( !this->isStationary() && this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        Build_TransientTerm=buildNonCstPart && !UseJacobianLinearTerms;

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& v = this->fieldTemperature();
    auto const u = Xh->element(XVec, this->rowStartInVector());

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() );

    for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( matName );
        if ( thermalConductivity.isConstant() )
        {
            if ( buildNonCstPart && !UseJacobianLinearTerms )
            {
                double kappa = thermalConductivity.value();
                myLinearForm +=
                    integrate( _range=range,
                               _expr= kappa*inner(gradv(u),grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            auto kappa = thermalConductivity.expr();
            std::string symbolStr = "heat_T";
            if ( kappa.expression().hasSymbol( symbolStr ) )
            {
                if ( buildNonCstPart )
                {
                    auto kappaEval = expr( kappa, symbolExpr( symbolStr, idv(u) ) );
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= kappaEval*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
            else
            {
                if ( buildNonCstPart && !UseJacobianLinearTerms )
                {
                    //auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= kappa*inner(gradv(u),grad(v)),
                                   _geomap=this->geomap() );
                }
            }
        }

        auto const& rhoHeatCapacity = this->thermalProperties()->rhoHeatCapacity( matName );
        if ( rhoHeatCapacity.isConstant() )
        {
            double rhoHeatCapacityValue = rhoHeatCapacity.value();
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                myLinearForm +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityValue*(gradv(u)*idv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
            }
            if ( !this->isStationary() )
            {
                if ( Build_TransientTerm )
                {
                    double thecoeff = rhoHeatCapacityValue*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= thecoeff*idv(u)*id(v),
                                   _geomap=this->geomap() );
                }
                if (buildCstPart)
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= -rhoHeatCapacityValue*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }
        else
        {
            auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
            if ( buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
            {
                myLinearForm +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityExpr*(gradv(u)*idv(this->fieldVelocityConvection()))*id(v),
                               _geomap=this->geomap() );
            }
            if ( !this->isStationary() )
            {
                if ( Build_TransientTerm )
                {
                    auto thecoeff = rhoHeatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= thecoeff*idv(u)*id(v),
                                   _geomap=this->geomap() );
                }
                if (buildCstPart)
                {
                    auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
                    myLinearForm +=
                        integrate( _range=range,
                                   _expr= -rhoHeatCapacityExpr*idv(rhsTimeStep)*id(v),
                                   _geomap=this->geomap() );
                }
            }
        }

        // update stabilization gls
        if ( M_stabilizationGLS && buildNonCstPart && this->fieldVelocityConvectionIsUsedAndOperational() )
        {
            CHECK( !thermalConductivity.isMatrix() && thermalConductivity.isConstant() && rhoHeatCapacity.isConstant() ) << "NotImplemented";
            this->updateResidualStabilizationGLS(  cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(this->fieldVelocityConvection()),range,data );
        }

    }


    this->updateResidualSourceTerm( R,buildCstPart ) ;
    this->updateResidualNeumannBC( R,buildCstPart );
    this->updateResidualRobinBC( u,R,buildCstPart );

    this->log("Heat","updateResidual", "finish");
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateResidualDofElimination","start" );

    vector_ptrtype& R = data.residual();
    auto mesh = this->mesh();
    auto u = this->spaceTemperature()->element( R,this->rowStartInVector() );
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("temperature");
    auto const& dofsWithValueImposedTemperature = ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposed->second : std::set<size_type>();
    for ( size_type thedof : dofsWithValueImposedTemperature )
        u.set( thedof,0. );
    sync( u, "=", dofsWithValueImposedTemperature );

    this->log("Heat","updateResidualDofElimination","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualNeumannBC( vector_ptrtype& R, bool buildCstPart ) const
{
    if ( this->M_bcNeumann.empty() ) return;

    if ( buildCstPart )
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = this->fieldTemperature();

        auto myLinearForm = form1( _test=Xh, _vector=R,
                                   _rowstart=this->rowStartInVector() );
        for( auto const& d : this->M_bcNeumann )
        {
            auto theExpr = expr( expression(d),this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualRobinBC( element_temperature_external_storage_type const& u, vector_ptrtype& R, bool buildCstPart ) const
{
    if ( this->M_bcRobin.empty() ) return;

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& v = this->fieldTemperature();

    auto myLinearForm = form1( _test=Xh, _vector=R,
                               _rowstart=this->rowStartInVector() );
    for( auto const& d : this->M_bcRobin )
    {
        auto theExpr1 = expr( expression1(d),this->symbolsExpr() );
        if ( !buildCstPart )
        {
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= theExpr1*idv(u)*id(v),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            auto theExpr2 = expr( expression2(d),this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= -theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateResidualSourceTerm( vector_ptrtype& R, bool buildCstPart ) const
{
    if ( this->M_volumicForcesProperties.empty() ) return;

    if ( buildCstPart )
    {
        auto myLinearForm = form1( _test=this->spaceTemperature(), _vector=R,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldTemperature();

        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto theExpr = expr( expression(d),this->symbolsExpr() );
            auto rangeBodyForceUsed = ( markers(d).empty() )? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= -theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateJacobianDofElimination","start" );

    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=RBis,_expr=cst(0.) );
    }

    this->log("Heat","updateJacobianDofElimination","finish" );

}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateJacobianRobinBC( sparse_matrix_ptrtype& J, bool buildCstPart ) const
{
    if ( this->M_bcRobin.empty() ) return;

    if ( !buildCstPart )
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = this->fieldTemperature();

        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=J,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );
        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expr( expression1(d),this->symbolsExpr() );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const 
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;
    //if ( this->M_bcDirichlet.empty() ) return;

    this->log("Heat","updateLinearPDEDofElimination","start" );

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcDirichlet )
    {
        auto theExpr = expr( expression(d),this->symbolsExpr() );
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                _element=u,_rhs=F,_expr=theExpr );
    }

    this->log("Heat","updateLinearPDEDofElimination","finish" );
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDESourceTerm( vector_ptrtype& F, bool buildCstPart ) const
{
    if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        this->M_overwritemethod_updateSourceTermLinearPDE(F,buildCstPart);
        return;
    }

    if ( this->M_volumicForcesProperties.empty() ) return;

    if ( !buildCstPart )
    {
        auto myLinearForm = form1( _test=this->spaceTemperature(), _vector=F,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = this->fieldTemperature();

        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto theExpr = expr( expression(d),this->symbolsExpr() );
            auto rangeBodyForceUsed = ( markers(d).empty() )? this->rangeMeshElements() : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= theExpr*id(v),
                           _geomap=this->geomap() );
        }
    }
}

HEAT_CLASS_TEMPLATE_DECLARATIONS
void
HEAT_CLASS_TEMPLATE_TYPE::updateLinearPDEWeakBC(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{
    if ( this->M_bcNeumann.empty() && this->M_bcRobin.empty() ) return;

    if ( !buildCstPart )
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = this->fieldTemperature();

        auto myLinearForm = form1( _test=Xh, _vector=F,
                                   _rowstart=this->rowStartInVector() );
        for( auto const& d : this->M_bcNeumann )
        {
            auto theExpr = expr( expression(d),this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= theExpr*id(v),
                           _geomap=this->geomap() );
        }

        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );
        for( auto const& d : this->M_bcRobin )
        {
            auto theExpr1 = expr( expression1(d),this->symbolsExpr() );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= theExpr1*idt(v)*id(v),
                           _geomap=this->geomap() );
            auto theExpr2 = expr( expression2(d),this->symbolsExpr() );
            myLinearForm +=
                integrate( _range=markedfaces(mesh,this->markerRobinBC( name(d) ) ),
                           _expr= theExpr1*theExpr2*id(v),
                           _geomap=this->geomap() );
        }

    }
}




} // end namespace FeelModels
} // end namespace Feel
