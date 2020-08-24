/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinear_Turbulence( DataUpdateLinear & data ) const
{
    const vector_ptrtype& sol = data.currentSolution();
    auto mfields_turbulence = M_turbulenceModelType->template modelFields<FilterBasisUnknownTurbulenceModel>( sol );
    auto mfields_fluid = this->modelFields().exclude( mfields_turbulence ); ;
    auto se = Feel::vf::symbolsExpr( M_turbulenceModelType->symbolsExpr( mfields_turbulence ), this->symbolsExpr( mfields_fluid ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields_turbulence), std::move( se ) );
#if 0
    auto mfieldsPrevious_fluid = M_fm->modelFields( M_previousSolFluid );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mfieldsPrevious_cfpdes = M_cfpdes->template modelFields<cfpdess_filter_basis_type>( previousSol );
        auto sePrevious = Feel::vf::symbolsExpr( M_cfpdes->symbolsExpr( mfieldsPrevious_cfpdes ), mfieldsPrevious_fluid.symbolsExpr() );
        auto mctxPrevious = Feel::FeelModels::modelContext( std::move( mfieldsPrevious_cfpdes ), std::move( sePrevious ) );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }
#endif

    M_turbulenceModelType->template updateLinearPDE<FilterBasisUnknownTurbulenceModel>( data,mctx );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearDofElimination_Turbulence( DataUpdateLinear & data ) const
{
    auto mfields_turbulence = M_turbulenceModelType->template modelFields<FilterBasisUnknownTurbulenceModel>();
    auto mfields_fluid = this->modelFields().exclude( mfields_turbulence ); ;
    auto se = Feel::vf::symbolsExpr( M_turbulenceModelType->symbolsExpr( mfields_turbulence ), this->symbolsExpr( mfields_fluid ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields_turbulence), std::move( se ) );
    M_turbulenceModelType->template updateLinearPDEDofElimination<FilterBasisUnknownTurbulenceModel>( data, mctx );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess_Turbulence( DataNewtonInitialGuess & data ) const
{
    vector_ptrtype& U = data.initialGuess();
    auto mfields_turbulence = M_turbulenceModelType->template modelFields<FilterBasisUnknownTurbulenceModel>( U );
    auto mfields_fluid = this->modelFields().exclude( mfields_turbulence ); ;
    auto se = Feel::vf::symbolsExpr( M_turbulenceModelType->symbolsExpr( mfields_turbulence ), this->symbolsExpr( mfields_fluid ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields_turbulence), std::move( se ) );
    M_turbulenceModelType->template updateNewtonInitialGuess<FilterBasisUnknownTurbulenceModel>( data, mctx );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidual_Turbulence( DataUpdateResidual & data ) const
{
    const vector_ptrtype& sol = data.currentSolution();
    auto mfields_turbulence = M_turbulenceModelType->template modelFields<FilterBasisUnknownTurbulenceModel>( sol );
    auto mfields_fluid = this->modelFields().exclude( mfields_turbulence ); ;
    auto se = Feel::vf::symbolsExpr( M_turbulenceModelType->symbolsExpr( mfields_turbulence ), this->symbolsExpr( mfields_fluid ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields_turbulence), std::move( se ) );
    M_turbulenceModelType->template updateResidual<FilterBasisUnknownTurbulenceModel>( data, mctx );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobian_Turbulence( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& sol = data.currentSolution();
    auto mfields_turbulence = M_turbulenceModelType->template modelFields<FilterBasisUnknownTurbulenceModel>( sol );
    auto mfields_fluid = this->modelFields().exclude( mfields_turbulence ); ;
    auto se = Feel::vf::symbolsExpr( M_turbulenceModelType->symbolsExpr( mfields_turbulence ), this->symbolsExpr( mfields_fluid ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields_turbulence), std::move( se ) );
    M_turbulenceModelType->template updateJacobian<FilterBasisUnknownTurbulenceModel>( data, mctx );
}


} // end namespace FeelModels
} // end namespace Feel


