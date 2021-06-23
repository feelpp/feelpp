#include <feel/feelvf/vf.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

namespace Feel {
namespace FeelModels {

namespace ADRDetails {

//--------------------------------------------------------------------//
//-------------------------- Trial operators -------------------------//
//--------------------------------------------------------------------//
// Advection
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::Advection )
{
    return gradt(phi)*u;
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::Advection )
{
    return gradt(phi)*u + polyDerivCoeff*idt(phi);
}
//--------------------------------------------------------------------//
// Advection-Diffusion
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::AdvectionDiffusion )
{
    return gradt(phi)*u - D*laplaciant(phi);
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::AdvectionDiffusion )
{
    return gradt(phi)*u + polyDerivCoeff*idt(phi) - D*laplaciant(phi);
}
//--------------------------------------------------------------------//
// Advection-Reaction
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::AdvectionReaction )
{
    return gradt(phi)*u + R*idt(phi);
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::AdvectionReaction )
{
    return gradt(phi)*u + polyDerivCoeff*idt(phi) + R*idt(phi);
}
//--------------------------------------------------------------------//
// Advection-Diffusion-Reaction
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::AdvectionDiffusionReaction )
{
    return gradt(phi)*u + R*idt(phi) - D*laplaciant(phi);
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTrial( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::AdvectionDiffusionReaction )
{
    return gradt(phi)*u + polyDerivCoeff*idt(phi) + R*idt(phi) - D*laplaciant(phi);
}

//--------------------------------------------------------------------//
//--------------------------- Test operators -------------------------//
//--------------------------------------------------------------------//
// Advection
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::Advection )
{
    return grad(phi)*u;
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::Advection )
{
    return grad(phi)*u + polyDerivCoeff*id(phi);
}
//--------------------------------------------------------------------//
// Advection-Diffusion
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::AdvectionDiffusion )
{
    return grad(phi)*u - D*laplacian(phi);
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::AdvectionDiffusion )
{
    return grad(phi)*u + polyDerivCoeff*id(phi) - D*laplacian(phi);
}
//--------------------------------------------------------------------//
// Advection-Reaction
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::AdvectionReaction )
{
    return grad(phi)*u + R*id(phi);
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::AdvectionReaction )
{
    return grad(phi)*u + polyDerivCoeff*id(phi) + R*id(phi);
}
//--------------------------------------------------------------------//
// Advection-Diffusion-Reaction
template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType>
auto
opStabStationaryTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, ADRTypes::AdvectionDiffusionReaction )
{
    return grad(phi)*u + R*id(phi) - D*laplacian(phi);
}

template<typename AdvectionVelocityExprType, typename DiffusionCoeffExprType, typename ReactionCoeffExprType, typename PhiFieldType, typename TransientCoeffExprType>
auto
opStabTransientTest( AdvectionVelocityExprType const& u, DiffusionCoeffExprType const& D, ReactionCoeffExprType const& R, PhiFieldType const& phi, TransientCoeffExprType const& polyDerivCoeff, ADRTypes::AdvectionDiffusionReaction )
{
    return grad(phi)*u + polyDerivCoeff*id(phi) + R*id(phi) - D*laplacian(phi);
}

} // ADRDetails

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
template< 
    typename FunctionSpaceType,
    typename FunctionSpaceAdvectionVelocityType,
    typename BasisDiffusionCoeffType,
    typename BasisReactionCoeffType
        >
template<typename ExprT>
void
AdvDiffReac<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType>::updateLinearPDEStabilization( DataUpdateLinear & data, vf::Expr<ExprT> const& advectionVelocity )
{
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    std::string sc=(BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("AdvDiffReac","updateLinearPDEStabilization", "start"+sc );
    this->timerTool("Solve").start();

    uint16_type adrt = ADREnum::Advection;
    if( this->hasReaction() ) adrt |= ADREnum::Reaction;
    if( this->hasDiffusion() && nOrder >= 2 ) adrt |= ADREnum::Diffusion;
    ADREnum adrtype = static_cast<ADREnum>( adrt );
    
    switch ( this->stabilizationMethod() )
    {
        case AdvectionStabMethod::NONE : { break; } // remove -Wswitch warning

        case AdvectionStabMethod::GALS :
        {
            this->log("AdvDiffReac","updateLinearPDEStabilization", "assembly GLS stabilization" );
            if( adrtype == Advection ) this->updateLinearPDEStabilizationGLS<ADRTypes::Advection>( data, advectionVelocity );
            if( adrtype == AdvectionDiffusion) this->updateLinearPDEStabilizationGLS<ADRTypes::AdvectionDiffusion>( data, advectionVelocity );
            if( adrtype == AdvectionReaction) this->updateLinearPDEStabilizationGLS<ADRTypes::AdvectionReaction>( data, advectionVelocity );
            if( adrtype == AdvectionDiffusionReaction) this->updateLinearPDEStabilizationGLS<ADRTypes::AdvectionDiffusionReaction>( data, advectionVelocity );
        } //GALS
        break;

        case AdvectionStabMethod::SUPG :
        {
            this->log("AdvDiffReac","updateLinearPDEStabilization", "assembly SUPG stabilization" );
            if( adrtype == Advection ) this->updateLinearPDEStabilizationSUPG<ADRTypes::Advection>( data, advectionVelocity );
            if( adrtype == AdvectionDiffusion) this->updateLinearPDEStabilizationSUPG<ADRTypes::AdvectionDiffusion>( data, advectionVelocity );
            if( adrtype == AdvectionReaction) this->updateLinearPDEStabilizationSUPG<ADRTypes::AdvectionReaction>( data, advectionVelocity );
            if( adrtype == AdvectionDiffusionReaction) this->updateLinearPDEStabilizationSUPG<ADRTypes::AdvectionDiffusionReaction>( data, advectionVelocity );
        } //SUPG
        break;

        case AdvectionStabMethod::SGS :
        {
            this->log("AdvDiffReac","updateLinearPDEStabilization", "assembly SGS stabilization" );
            if( adrtype == Advection ) this->updateLinearPDEStabilizationSGS<ADRTypes::Advection>( data, advectionVelocity );
            if( adrtype == AdvectionDiffusion) this->updateLinearPDEStabilizationSGS<ADRTypes::AdvectionDiffusion>( data, advectionVelocity );
            if( adrtype == AdvectionReaction) this->updateLinearPDEStabilizationSGS<ADRTypes::AdvectionReaction>( data, advectionVelocity );
            if( adrtype == AdvectionDiffusionReaction) this->updateLinearPDEStabilizationSGS<ADRTypes::AdvectionDiffusionReaction>( data, advectionVelocity );
        } //SGS
        break;

        case AdvectionStabMethod::CIP :
        {
            this->log("AdvDiffReac","updateLinearPDEStabilization", "assembly CIP stabilization" );
            if( BuildNonCstPart )
            {
                auto mesh = this->mesh();
                auto space = this->functionSpace();
                auto const& phi = this->fieldSolution();
                sparse_matrix_ptrtype & A = data.matrix();

                auto beta = idv(this->fieldAdvectionVelocity());
                auto beta_norm = vf::sqrt(trans(beta)*beta);
                double stabCoeff = this->M_stabilizationCIPCoefficient;
                auto coeff = stabCoeff * M_gamma1 * hFace() * hFace() * beta_norm;

                auto bilinearForm_PatternExtended = form2( 
                        _test=space, _trial=space, _matrix=A, _pattern=size_type(Pattern::EXTENDED) 
                        );
                bilinearForm_PatternExtended += integrate(
                        _range=internalfaces(mesh),
                        _expr=coeff * inner(jumpt(gradt(phi)), jump(grad(phi)))
                        );
            }

        } //CIP
        break;
    } //switch

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("AdvDiffReac","updateLinearPDEStabilization","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}
//--------------------------------------------------------------------//
// Dispatch
template< 
    typename FunctionSpaceType,
    typename FunctionSpaceAdvectionVelocityType,
    typename BasisDiffusionCoeffType,
    typename BasisReactionCoeffType
        >
template<typename ADRT, typename ExprT>
void
AdvDiffReac<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType>::updateLinearPDEStabilizationGLS( DataUpdateLinear & data, vf::Expr<ExprT> const& advectionVelocity )
{
    using namespace ADRDetails;

    bool BuildCstPart = data.buildCstPart();
    //bool BuildNonCstPart = !BuildCstPart;
    if( BuildCstPart )
        return;

    auto mesh = this->mesh();
    auto space = this->functionSpace();

    sparse_matrix_ptrtype & A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& phi = this->fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=this->rowStartInVector() );

    ADRT adrType;

    double sigma = this->isStationary() ? 0: this->timeStepBDF()->polyDerivCoefficient(0);
    auto const u = advectionVelocity;
    auto const D = idv(this->diffusionReactionModel()->fieldDiffusionCoeff());
    auto const R = idv(this->diffusionReactionModel()->fieldReactionCoeff());
    auto const uNorm = vf::sqrt(trans(u)*u);
#if 0
    auto coeff  = val( 1/( 2*uNorm*AdvectionType::nOrder/h() + std::abs(sigma) ));
#else
    //auto coeff/*tau*/ = M_stabilizationGLSParameter->tau( uconv, kappa, mpl::int_<0/*StabParamType*/>() );
    auto coeff_expr = Feel::FeelModels::stabilizationGLSParameterExpr( *(this->stabilizationGLSParameter()), u, D, true, this->hasDiffusion() );
    auto coeffP0d = vf::project(
            _space=this->functionSpaceP0d(),
            _range=elements(mesh),
            _expr=coeff_expr
            );
    auto const coeff = idv(coeffP0d);
#endif

    if( this->isStationary() )
    {
        auto L_op = opStabStationaryTest( u, D, R, phi, adrType);
        auto L_opt = opStabStationaryTrial( u, D, R, phi, adrType);
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=this->geomap() 
                );

        for( auto const& d : this->volumicSources() )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( 
                    _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=this->geomap() 
                    );
        }
        if( this->hasSourceAdded() )
        {
            linearForm += integrate( 
                    _range=this->rangeMeshElements(),
                    _expr=coeff * inner(L_op, idv(this->sourceAdded())),
                    _geomap=this->geomap() );
        }
    }
    else
    {
        double polyDerivCoeff = this->timeStepBDF()->polyDerivCoefficient(0);
        auto polyDeriv = this->timeStepBDF()->polyDeriv();

        auto L_op = opStabTransientTest( u, D, R, phi, polyDerivCoeff, adrType);
        auto L_opt = opStabTransientTrial( u, D, R, phi, polyDerivCoeff, adrType);
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=this->geomap() 
                );


        linearForm += integrate( 
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, idv(polyDeriv)),
                _geomap=this->geomap() 
                );

        for( auto const& d : this->volumicSources() )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( 
                    _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=this->geomap() 
                    );
        }
        if( this->hasSourceAdded() )
        {
            linearForm += integrate( 
                    _range=this->rangeMeshElements(),
                    _expr=coeff * inner(L_op, idv(this->sourceAdded())),
                    _geomap=this->geomap() 
                    );
        }
    }
}

template< 
    typename FunctionSpaceType,
    typename FunctionSpaceAdvectionVelocityType,
    typename BasisDiffusionCoeffType,
    typename BasisReactionCoeffType
        >
template<typename ADRT, typename ExprT>
void
AdvDiffReac<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType>::updateLinearPDEStabilizationSUPG( DataUpdateLinear & data, vf::Expr<ExprT> const& advectionVelocity )
{
    using namespace ADRDetails;

    bool BuildCstPart = data.buildCstPart();
    //bool BuildNonCstPart = !BuildCstPart;
    if( BuildCstPart )
        return;

    auto mesh = this->mesh();
    auto space = this->functionSpace();

    sparse_matrix_ptrtype & A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& phi = this->fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=this->rowStartInVector() );

    ADRT adrType;

    auto const u = advectionVelocity;
    auto const D = idv(this->diffusionReactionModel()->fieldDiffusionCoeff());
    auto const R = idv(this->diffusionReactionModel()->fieldReactionCoeff());
    auto const uNorm = vf::sqrt(trans(u)*u);
#if 0
    auto coeff = val(vf::h() / (2 * uNorm + 0.001));
#else
    //auto coeff/*tau*/ = M_stabilizationGLSParameter->tau( uconv, kappa, mpl::int_<0/*StabParamType*/>() );
    auto coeff_expr = Feel::FeelModels::stabilizationGLSParameterExpr( *(this->stabilizationGLSParameter()), u, D, true, this->hasDiffusion() );
    auto coeffP0d = vf::project(
            _space=this->functionSpaceP0d(),
            _range=this->rangeMeshElements(),
            _expr=coeff_expr
            );
    auto const coeff = idv(coeffP0d);
#endif

    if( this->isStationary() )
    {
        auto L_op = opStabStationaryTest( u, D, R, phi, ADRTypes::Advection{} );
        auto L_opt = opStabStationaryTrial( u, D, R, phi, adrType);
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=this->geomap() 
                );
        for( auto const& d : this->volumicSources() )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( 
                    _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=this->geomap() 
                    );
        }
        if( this->hasSourceAdded() )
        {
            linearForm += integrate( 
                    _range=this->rangeMeshElements(),
                    _expr=coeff * inner(L_op, idv(this->sourceAdded())),
                    _geomap=this->geomap() 
                    );
        }
    }
    else
    {
        double polyDerivCoeff = this->timeStepBDF()->polyDerivCoefficient(0);
        auto polyDeriv = this->timeStepBDF()->polyDeriv();

        auto L_op = opStabTransientTest( u, D, R, phi, polyDerivCoeff, ADRTypes::Advection{} );
        auto L_opt = opStabTransientTrial( u, D, R, phi, polyDerivCoeff, adrType);
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=this->geomap() 
                );

        linearForm += integrate( 
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, idv(polyDeriv)),
                _geomap=this->geomap() 
                );

        for( auto const& d : this->volumicSources() )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( 
                    _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=this->geomap() 
                    );
        }
        if( this->hasSourceAdded() )
        {
            linearForm += integrate( 
                    _range=this->rangeMeshElements(),
                    _expr=coeff * inner(L_op, idv(this->sourceAdded())),
                    _geomap=this->geomap() 
                    );
        }
    }
}

template< 
    typename FunctionSpaceType,
    typename FunctionSpaceAdvectionVelocityType,
    typename BasisDiffusionCoeffType,
    typename BasisReactionCoeffType
        >
template<typename ADRT, typename ExprT>
void
AdvDiffReac<FunctionSpaceType, FunctionSpaceAdvectionVelocityType, BasisDiffusionCoeffType, BasisReactionCoeffType>::updateLinearPDEStabilizationSGS( DataUpdateLinear & data, vf::Expr<ExprT> const& advectionVelocity )
{
    using namespace ADRDetails;

    bool BuildCstPart = data.buildCstPart();
    //bool BuildNonCstPart = !BuildCstPart;
    if( BuildCstPart )
        return;

    auto mesh = this->mesh();
    auto space = this->functionSpace();

    sparse_matrix_ptrtype & A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& phi = this->fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=this->rowStartInVector() );

    ADRT adrType;

    double sigma = this->isStationary() ? 0: this->timeStepBDF()->polyDerivCoefficient(0);
    auto const u = advectionVelocity;
    auto const D = idv(this->diffusionReactionModel()->fieldDiffusionCoeff());
    auto const R = idv(this->diffusionReactionModel()->fieldReactionCoeff());
    auto const uNorm = vf::sqrt(trans(u)*u);

    auto coeff = val(1. / ( (2 * uNorm) / vf::h()  + vf::abs(sigma) ) );

    if( this->isStationary() )
    {
        auto L_op = grad(phi) * u // advection term
            - (!this->isStationary()) * this->timeStepBDF()->polyDerivCoefficient(0)*id(phi) // transient term
            - (this->hasReaction()) * R *id(phi) // reaction term
            - (this->hasDiffusion() && nOrder >= 2) * (-D)*laplacian(phi) // diffusion term
            ;
        auto L_opt = opStabStationaryTrial( u, D, R, phi, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=this->geomap() );

        for( auto const& d : this->volumicSources() )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( 
                    _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=this->geomap() 
                    );
        }
        if( this->hasSourceAdded() )
        {
            linearForm += integrate( 
                    _range=this->rangeMeshElements(),
                    _expr=coeff * inner(L_op, idv(this->sourceAdded())),
                    _geomap=this->geomap() 
                    );
        }
    }
    else
    {
        double polyDerivCoeff = this->timeStepBDF()->polyDerivCoefficient(0);
        auto polyDeriv = this->timeStepBDF()->polyDeriv();

        auto L_op = grad(phi) * u // advection term
            - (!this->isStationary()) * this->timeStepBDF()->polyDerivCoefficient(0)*id(phi) // transient term
            - (this->hasReaction()) * R *id(phi) // reaction term
            - (this->hasDiffusion() && nOrder >= 2) * (-D)*laplacian(phi) // diffusion term
            ;
        auto L_opt = opStabTransientTrial( u, D, R, phi, polyDerivCoeff, adrType);
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=this->geomap() 
                );


        linearForm += integrate( 
                _range=this->rangeMeshElements(),
                _expr=coeff * inner(L_op, idv(polyDeriv)),
                _geomap=this->geomap() 
                );

        for( auto const& d : this->volumicSources() )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( 
                    _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=this->geomap() 
                    );
        }
        if( this->hasSourceAdded() )
        {
            linearForm += integrate( 
                    _range=this->rangeMeshElements(),
                    _expr=coeff * inner(L_op, idv(this->sourceAdded())),
                    _geomap=this->geomap() 
                    );
        }
    }
}

}
}
