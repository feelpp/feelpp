#ifndef FEELPP_TOOLBOXES_ADR_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_ADR_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

ADVDIFFREAC_CLASS_TEMPLATE_DECLARATIONS
template <typename ModelContextType>
void
ADVDIFFREAC_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool build_AdvectiveTerm = BuildNonCstPart;
    bool build_DiffusionTerm = BuildNonCstPart;
    bool build_ReactionTerm = BuildNonCstPart;
    //bool build_SourceTerm = BuildNonCstPart;
    //bool build_BoundaryNeumannTerm = BuildNonCstPart;

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("AdvDiffReac","updateLinearPDE", "start"+sc );
    this->timerTool("Solve").start();

    auto const& se = mctx.symbolsExpr();

    auto mesh = this->mesh();
    auto space = this->functionSpace();

    auto const& advection_velocity = this->fieldAdvectionVelocity();

    auto const& phi = this->fieldSolution();
    auto const& psi = this->fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F );

    // Advection
    if(this->hasAdvection() && build_AdvectiveTerm)
    {
        this->timerTool("Solve").start();
        
        //bilinearForm += integrate(
                //_range=this->rangeMeshElements(),
                //_expr=inner((gradt(phi)*idv(advection_velocity)), id(psi)),
                //_geomap=this->geomap()
                //);
        this->M_functionAssemblyLinearAdvection( data );

        double timeElapsedAdvection = this->timerTool("Solve").stop();
        this->log("AdvDiffReac","updateLinearPDE","assembly advection terms in "+(boost::format("%1% s") %timeElapsedAdvection).str() );
    }

    // Diffusion
    if( this->hasDiffusion() && build_DiffusionTerm )
    {
        this->timerTool("Solve").start();

        auto const& D = this->diffusionReactionModel()->fieldDiffusionCoeff();

        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=inner(trans(gradt(phi)),idv(D)*trans(grad(psi))),
                _geomap=this->geomap()
                );

        double timeElapsedDiffusion = this->timerTool("Solve").stop();
        this->log("AdvDiffReac","updateLinearPDE","assembly diffusion terms in "+(boost::format("%1% s") %timeElapsedDiffusion).str() );
    }

    // Reaction
    if( this->hasReaction() && build_ReactionTerm )
    {
        this->timerTool("Solve").start();

        auto const& R = this->diffusionReactionModel()->fieldReactionCoeff();
        
        bilinearForm += integrate(
                _range=this->rangeMeshElements(),
                _expr=inner(idv(R)*idt(phi), id(psi)),
                _geomap=this->geomap()
                );

        double timeElapsedReaction = this->timerTool("Solve").stop();
        this->log("AdvDiffReac","updateLinearPDE","assembly reaction terms in "+(boost::format("%1% s") %timeElapsedReaction).str() );
    }

    // Transient terms
    if (!this->isStationary())
    {
        this->timerTool("Solve").start();

        bool build_Form2TransientTerm = BuildNonCstPart;
        bool build_Form1TransientTerm = BuildNonCstPart;
        if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        {
            build_Form2TransientTerm = BuildCstPart;
        }

        if (build_Form2TransientTerm)
        {
            bilinearForm +=
                integrate( _range=this->rangeMeshElements(),
                           _expr=M_bdf->polyDerivCoefficient(0)*inner(idt(phi),id(psi)),
                           _geomap=this->geomap() );
        }
        if (build_Form1TransientTerm)
        {
            linearForm +=
                integrate( _range=this->rangeMeshElements(),
                           _expr=inner(idv(M_bdf->polyDeriv()),id(psi)),
                           _geomap=this->geomap() );
        }

        double timeElapsedTransient = this->timerTool("Solve").stop();
        this->log("AdvDiffReac","updateLinearPDE","assembly transient terms in "+(boost::format("%1% s") %timeElapsedTransient).str() );
    }

    // Source term
    if ( BuildNonCstPart )
    {
        for( auto const& d : this->M_sources )
        {
            auto rangeUsed = ( markers(d).empty() )? elements(mesh) : markedelements(mesh,markers(d));
            linearForm += integrate( _range=rangeUsed,
                                     _expr= inner(expression(d,se),id(psi)),
                                     _geomap=this->geomap() );

        }
    }

    if( this->hasSourceAdded() && !BuildCstPart )
    {
        linearForm +=
            integrate( _range=this->rangeMeshElements(),
                       _expr= inner(idv(M_fieldSource),id(psi)),
                       _geomap=this->geomap() );
    }

    // User-defined additional terms
    this->updateLinearPDEAdditional( A, F, BuildCstPart );

    // Stabilization
    if ( this->hasAdvection() )
        this->updateLinearPDEStabilization( data );

    // Boundary conditions
    this->updateWeakBCLinearPDE(A, F, BuildCstPart);

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("AdvDiffReac","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );

}

} // namespace FeelModels
} // namespace Feel

#endif
