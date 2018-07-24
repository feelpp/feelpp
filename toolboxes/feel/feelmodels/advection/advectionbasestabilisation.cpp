#include <feel/feelvf/vf.hpp>

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

//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
//--------------------------------------------------------------------//
// Dispatch
template<typename ADRT, typename AdvectionType>
void
updateLinearPDEStabilizationGLS( AdvectionType const& adv, ModelAlgebraic::DataUpdateLinear & data )
{
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    if( BuildCstPart )
        return;

    auto mesh = adv.mesh();
    auto space = adv.functionSpace();

    sparse_matrix_ptrtype & A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& phi = adv.fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=adv.rowStartInVector() );

    ADRT adrType;

    double sigma = adv.isStationary() ? 0: adv.timeStepBDF()->polyDerivCoefficient(0);
    auto u = idv(adv.fieldAdvectionVelocity());
    auto D = idv(adv.diffusionReactionModel()->fieldDiffusionCoeff());
    auto R = idv(adv.diffusionReactionModel()->fieldReactionCoeff());
    auto uNorm = vf::sqrt(trans(u)*u);
#if 0
    auto coeff  = val( 1/( 2*uNorm*AdvectionType::nOrder/h() + std::abs(sigma) ));
#else
    //auto coeff/*tau*/ = M_stabilizationGLSParameter->tau( uconv, kappa, mpl::int_<0/*StabParamType*/>() );
    auto coeff = Feel::FeelModels::stabilizationGLSParameterExpr( *(adv.stabilizationGLSParameter()), u, D, true, adv.hasDiffusion() );
#endif

    if( adv.isStationary() )
    {
        auto L_op = opStabStationaryTest( u, D, R, phi, adrType);
        auto L_opt = opStabStationaryTrial( u, D, R, phi, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=adv.geomap() );

        for( auto const& d : adv.volumicSources() )
        {
            auto rangeUsed = ( marker(d).empty() )? elements(mesh) : markedelements(mesh,marker(d));
            linearForm += integrate( _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=adv.geomap() );
        }
        if( adv.hasSourceAdded() )
        {
            linearForm += integrate( _range=elements(mesh),
                    _expr=coeff * inner(L_op, idv(adv.sourceAdded())),
                    _geomap=adv.geomap() );
        }
    }
    else
    {
        double polyDerivCoeff = adv.timeStepBDF()->polyDerivCoefficient(0);
        auto polyDeriv = adv.timeStepBDF()->polyDeriv();

        auto L_op = opStabTransientTest( u, D, R, phi, polyDerivCoeff, adrType);
        auto L_opt = opStabTransientTrial( u, D, R, phi, polyDerivCoeff, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=adv.geomap() );


        linearForm += integrate( 
                _range=elements(mesh),
                _expr=coeff * inner(L_op, idv(polyDeriv)),
                _geomap=adv.geomap() );

        for( auto const& d : adv.volumicSources() )
        {
            auto rangeUsed = ( marker(d).empty() )? elements(mesh) : markedelements(mesh,marker(d));
            linearForm += integrate( _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=adv.geomap() );
        }
        if( adv.hasSourceAdded() )
        {
            linearForm += integrate( _range=elements(mesh),
                    _expr=coeff * inner(L_op, idv(adv.sourceAdded())),
                    _geomap=adv.geomap() );
        }
    }
}

template<typename ADRT, typename AdvectionType>
void
updateLinearPDEStabilizationSUPG( AdvectionType const& adv, ModelAlgebraic::DataUpdateLinear & data )
{
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    if( BuildCstPart )
        return;

    auto mesh = adv.mesh();
    auto space = adv.functionSpace();

    sparse_matrix_ptrtype & A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& phi = adv.fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=adv.rowStartInVector() );

    ADRT adrType;

    double sigma = adv.isStationary() ? 0: adv.timeStepBDF()->polyDerivCoefficient(0);
    auto u = idv(adv.fieldAdvectionVelocity());
    auto D = idv(adv.diffusionReactionModel()->fieldDiffusionCoeff());
    auto R = idv(adv.diffusionReactionModel()->fieldReactionCoeff());
    auto uNorm = vf::sqrt(trans(u)*u);
#if 0
    auto coeff = val(vf::h() / (2 * uNorm + 0.001));
#else
    //auto coeff/*tau*/ = M_stabilizationGLSParameter->tau( uconv, kappa, mpl::int_<0/*StabParamType*/>() );
    auto coeff = Feel::FeelModels::stabilizationGLSParameterExpr( *(adv.stabilizationGLSParameter()), u, D, true, adv.hasDiffusion() );
#endif

    if( adv.isStationary() )
    {
        auto L_op = opStabStationaryTest( u, D, R, phi, ADRTypes::Advection{} );
        auto L_opt = opStabStationaryTrial( u, D, R, phi, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=adv.geomap() );
        for( auto const& d : adv.volumicSources() )
        {
            auto rangeUsed = ( marker(d).empty() )? elements(mesh) : markedelements(mesh,marker(d));
            linearForm += integrate( _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=adv.geomap() );
        }
        if( adv.hasSourceAdded() )
        {
            linearForm += integrate( _range=elements(mesh),
                    _expr=coeff * inner(L_op, idv(adv.sourceAdded())),
                    _geomap=adv.geomap() );
        }
    }
    else
    {
        double polyDerivCoeff = adv.timeStepBDF()->polyDerivCoefficient(0);
        auto polyDeriv = adv.timeStepBDF()->polyDeriv();

        auto L_op = opStabTransientTest( u, D, R, phi, polyDerivCoeff, ADRTypes::Advection{} );
        auto L_opt = opStabTransientTrial( u, D, R, phi, polyDerivCoeff, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=adv.geomap() );


        linearForm += integrate( 
                _range=elements(mesh),
                _expr=coeff * inner(L_op, idv(polyDeriv)),
                _geomap=adv.geomap() );

        for( auto const& d : adv.volumicSources() )
        {
            auto rangeUsed = ( marker(d).empty() )? elements(mesh) : markedelements(mesh,marker(d));
            linearForm += integrate( _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=adv.geomap() );
        }
        if( adv.hasSourceAdded() )
        {
            linearForm += integrate( _range=elements(mesh),
                    _expr=coeff * inner(L_op, idv(adv.sourceAdded())),
                    _geomap=adv.geomap() );
        }
    }
}

template<typename ADRT, typename AdvectionType>
void
updateLinearPDEStabilizationSGS( AdvectionType const& adv, ModelAlgebraic::DataUpdateLinear & data )
{
    bool BuildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !BuildCstPart;
    if( BuildCstPart )
        return;

    auto mesh = adv.mesh();
    auto space = adv.functionSpace();

    sparse_matrix_ptrtype & A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto const& phi = adv.fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F, _rowstart=adv.rowStartInVector() );

    ADRT adrType;

    double sigma = adv.isStationary() ? 0: adv.timeStepBDF()->polyDerivCoefficient(0);
    auto u = idv(adv.fieldAdvectionVelocity());
    auto D = idv(adv.diffusionReactionModel()->fieldDiffusionCoeff());
    auto R = idv(adv.diffusionReactionModel()->fieldReactionCoeff());
    auto uNorm = vf::sqrt(trans(u)*u);

    auto coeff = val(1. / ( (2 * uNorm) / vf::h()  + vf::abs(sigma) ) );

    if( adv.isStationary() )
    {
        auto L_op = grad(phi) * u // advection term
            - (!adv.isStationary()) * adv.timeStepBDF()->polyDerivCoefficient(0)*id(phi) // transient term
            - (adv.hasReaction()) * R *id(phi) // reaction term
            - (adv.hasDiffusion() && AdvectionType::nOrder >= 2) * (-D)*laplacian(phi) // diffusion term
            ;
        auto L_opt = opStabStationaryTrial( u, D, R, phi, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=adv.geomap() );

        for( auto const& d : adv.volumicSources() )
        {
            auto rangeUsed = ( marker(d).empty() )? elements(mesh) : markedelements(mesh,marker(d));
            linearForm += integrate( _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=adv.geomap() );
        }
        if( adv.hasSourceAdded() )
        {
            linearForm += integrate( _range=elements(mesh),
                    _expr=coeff * inner(L_op, idv(adv.sourceAdded())),
                    _geomap=adv.geomap() );
        }
    }
    else
    {
        double polyDerivCoeff = adv.timeStepBDF()->polyDerivCoefficient(0);
        auto polyDeriv = adv.timeStepBDF()->polyDeriv();

        auto L_op = grad(phi) * u // advection term
            - (!adv.isStationary()) * adv.timeStepBDF()->polyDerivCoefficient(0)*id(phi) // transient term
            - (adv.hasReaction()) * R *id(phi) // reaction term
            - (adv.hasDiffusion() && AdvectionType::nOrder >= 2) * (-D)*laplacian(phi) // diffusion term
            ;
        auto L_opt = opStabTransientTrial( u, D, R, phi, polyDerivCoeff, adrType);
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=coeff * inner(L_op, L_opt),
                _geomap=adv.geomap() );


        linearForm += integrate( 
                _range=elements(mesh),
                _expr=coeff * inner(L_op, idv(polyDeriv)),
                _geomap=adv.geomap() );

        for( auto const& d : adv.volumicSources() )
        {
            auto rangeUsed = ( marker(d).empty() )? elements(mesh) : markedelements(mesh,marker(d));
            linearForm += integrate( _range=rangeUsed,
                    _expr=coeff * inner(L_op, expression(d)),
                    _geomap=adv.geomap() );
        }
        if( adv.hasSourceAdded() )
        {
            linearForm += integrate( _range=elements(mesh),
                    _expr=coeff * inner(L_op, idv(adv.sourceAdded())),
                    _geomap=adv.geomap() );
        }
    }
}

}
}
}
