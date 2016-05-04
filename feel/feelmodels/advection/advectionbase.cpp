#include <feel/feelmodels/advection/advectionbase.hpp>

namespace Feel {
namespace FeelModels {

//----------------------------------------------------------------------------//
// Static member initialization
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
const std::map<std::string, AdvectionStabMethod> 
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::AdvectionStabMethodIdMap = {
    {"NONE", AdvectionStabMethod::NONE},
    {"GALS", AdvectionStabMethod::GALS},
    {"CIP", AdvectionStabMethod::CIP},
    {"SUPG", AdvectionStabMethod::SUPG},
    {"SGS", AdvectionStabMethod::SGS}
};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Construction
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::AdvectionBase( 
        space_advection_ptrtype const& space,
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
:
    super_type( prefix, worldComm, subPrefix, rootRepository ),
    M_isUpdatedForUse(false),
    M_mesh(space->mesh()),
    M_space(space),
    M_gamma1(std::pow(nOrder, -3.5))
{
    this->loadParametersFromOptionsVm();

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".AdvectionConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".AdvectionSolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
}
    
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
typename ADVECTIONBASE_CLASS_TEMPLATE_TYPE::self_ptrtype 
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::New( 
        space_advection_ptrtype const& space,
        std::string const& prefix,
        WorldComm const& _worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
{
    self_ptrtype adv ( new self_type( space, prefix, _worldComm, subPrefix, rootRepository ) );
    return adv;
}

//----------------------------------------------------------------------------//
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::build()
{
    this->log("Advection","build", "start");
    
    // Algebraic data
    this->createAlgebraicData();
    // Bdf time scheme
    this->createTimeDiscretization();

    this->log("Advection","build", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::init(bool buildModelAlgebraicFactory, model_algebraic_factory_type::appli_ptrtype const& app )
{
    if ( M_isUpdatedForUse ) return;

    this->log("Advection","init", "start" );
    this->timerTool("Constructor").start();

    this->build();

    // Vector solution
    M_blockVectorSolution.resize(1);
    M_blockVectorSolution(0) = this->fieldSolutionPtr();
    M_blockVectorSolution.buildVector( this->backend() );
    
    // Algebraic factory
    if( buildModelAlgebraicFactory )
    {
        // matrix sparsity graph
        auto graph = this->buildMatrixGraph();
        
        M_algebraicFactory.reset( new model_algebraic_factory_type( app, this->backend(), graph, graph->mapRow().indexSplit() ) );
    }

    M_isUpdatedForUse = true;

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Advection","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::loadParametersFromOptionsVm()
{
    this->log("Advection","loadParametersFromOptionsVm", "start");
    
    // Stabilization method
    const std::string stabmeth = soption( _name="advec-stab-method", _prefix=this->prefix() );
    CHECK(AdvectionStabMethodIdMap.count(stabmeth)) << stabmeth <<" is not in the list of possible stabilization methods\n";
    M_stabMethod = AdvectionStabMethodIdMap.at(stabmeth);
    
    this->log("Advection","loadParametersFromOptionsVm", "finish");
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createAlgebraicData()
{
    this->log("Advection","createAlgebraicData","start");
    this->timerTool("Constructor").start();
    
    // Solution
    M_solution.reset( new element_advection_type(M_space, "phi") );
    // Backend
    M_backend = backend_type::build( soption(_name="backend"), this->prefix(), M_space->worldComm() );

    double tElapsed = this->timerTool("Constructor").stop("create");
    this->log("Advection","createAlgebraicData", (boost::format("finish in %1% s") %tElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretization()
{
    this->log("Advection","createTimeDiscretization", "start" );
    this->timerTool("Constructor").start();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdf = bdf( _vm=Environment::vm(), _space=M_space,
                       _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"phi"+suffixName)),
                       _prefix=this->prefix(),
                       // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                       _initial_time=this->timeInitial(),
                       _final_time=this->timeFinal(),
                       _time_step=this->timeStep(),
                       _restart=this->doRestart(),
                       _restart_path=this->restartPath(),
                       _restart_at_last_save=this->restartAtLastSave(),
                       _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_bdf->setfileFormat( myFileFormat );
    M_bdf->setPathSave( (fs::path(this->rootRepository()) /
                               fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%M_bdf->bdfOrder()%this->timeStep() ).str() ) ) ).string() );

    double tElapsed = this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("Advection","createTimeDiscretization", (boost::format("finish in %1% s") %tElapsed).str() );
}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Algebraic data
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::matrixPattern() const
{
    size_type pat = size_type(Pattern::COUPLED);

    if( this->stabilizationMethod() == AdvectionStabMethod::CIP )
    {
        pat = size_type(Pattern::EXTENDED);
    }

    return pat;
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
typename ADVECTIONBASE_CLASS_TEMPLATE_TYPE::graph_ptrtype
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    return stencil( 
            _test=this->functionSpace(), _trial=this->functionSpace(),
            _pattern=this->matrixPattern(),
            _diag_is_nonzero=true,
            _close=true
            )->graph();
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = this->functionSpace()->nLocalDofWithGhost();
    return res;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Time scheme
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateTimeStepBDF() 
{
    this->log("Advection","updateTimeStepBDF", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder = this->timeStepBDF()->timeOrder();

    M_bdf->next( *M_solution );

    int currentTimeOrder = this->timeStepBDF()->timeOrder();

    this->updateTime( M_bdf->time() );

    // maybe rebuild linear
    if ( M_algebraicFactory &&
         previousTimeOrder!=currentTimeOrder &&
         this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        if (!this->rebuildCstPartInLinearSystem())
        {
            this->log("Advection","updateTimeStepBDF", "do rebuildCstLinearPDE" );
            M_algebraicFactory->rebuildCstLinearPDE(M_solution);
        }
    }


    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("Advection","updateTimeStepBDF", "finish" );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Algebraic model updates
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();
    sparse_matrix_ptrtype& A_extended = data.matrixExtended();
    bool _BuildExtendedPart = data.buildExtendedPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool build_AdvectiveTerm = BuildNonCstPart;
    bool build_Form2TransientTerm = BuildNonCstPart;
    bool build_Form1TransientTerm = BuildNonCstPart;
    //bool build_SourceTerm = BuildNonCstPart;
    //bool build_BoundaryNeumannTerm = BuildNonCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        build_Form2TransientTerm=BuildCstPart;
    }

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Advection","updateLinearPDE", "start"+sc );
    this->timerTool("Solve").start();

    auto mesh = this->mesh();
    auto space = this->functionSpace();

    auto const& advection_velocity = this->fieldAdvectionVelocity();

    auto const& phi = this->fieldSolution();
    auto const& psi = this->fieldSolution();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto linearForm = form1( _test=space, _vector=F );

    // Advection
    if(build_AdvectiveTerm)
    {
        this->timerTool("Solve").start();
        
        bilinearForm += integrate(
                _range=elements(mesh),
                _expr=(gradt(phi)*idv(advection_velocity))*id(psi),
                _geomap=this->geomap()
                );

        double timeElapsedAdvection = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly advection terms in "+(boost::format("%1% s") %timeElapsedAdvection).str() );
    }

    // Transient terms
    if (!this->isStationary())
    {
        this->timerTool("Solve").start();
        
        if (build_Form2TransientTerm)
        {
            bilinearForm += integrate( 
                    _range=elements(mesh),
                    _expr=M_bdf->polyDerivCoefficient(0)*idt(phi)*id(psi),
                    _geomap=this->geomap() 
                    );

        }
        if (build_Form1TransientTerm)
        {
            linearForm += integrate(
                    _range=elements(mesh),
                    _expr=idv(M_bdf->polyDeriv())*id(psi),
                    _geomap=this->geomap() 
                    );
        }
        
        double timeElapsedTransient = this->timerTool("Solve").stop();
        this->log("Advection","updateLinearPDE","assembly transient terms in "+(boost::format("%1% s") %timeElapsedTransient).str() );
    }

    // Stabilization
    this->updateLinearPDEStabilization( data );

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Advection","updateLinearPDE","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateLinearPDEStabilization( DataUpdateLinear & data ) const
{
    using namespace Feel::vf;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _BuildCstPart = data.buildCstPart();
    sparse_matrix_ptrtype& A_extended = data.matrixExtended();
    bool _BuildExtendedPart = data.buildExtendedPart();
    
    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;

    bool Build_StabTerm = BuildNonCstPart;

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Advection","updateLinearPDEStabilization", "start"+sc );
    this->timerTool("Solve").start();
    
    auto mesh = this->mesh();
    auto space = this->functionSpace();

    auto phi = space->elementPtr();
    auto psi = space->elementPtr();

    double stabCoeff = this->stabilizationCoefficient();

    // Forms
    auto bilinearForm = form2( _test=space, _trial=space, _matrix=A );
    auto bilinearForm_PatternExtended = form2( _test=space, _trial=space, _matrix=A, _pattern=size_type(Pattern::EXTENDED) );
    auto linearForm = form1( _test=space, _vector=F );

    if(Build_StabTerm)
    {
        auto sigma = M_bdf->polyDerivCoefficient(0);
        auto beta = idv(this->fieldAdvectionVelocity());
        auto beta_norm = vf::sqrt(trans(beta)*beta);
        auto f = idv(M_bdf->polyDeriv());
        
        switch ( this->stabilizationMethod() )
        {
            case AdvectionStabMethod::NONE : { break; } // remove -Wswitch warning

            case AdvectionStabMethod::GALS :
            {
                auto coeff  = val( stabCoeff /
                        ( 2*beta_norm*nOrder/h() + std::abs(sigma) ));

                auto L_op = grad(psi)*val(beta) + val(sigma)*id(psi);
                auto L_opt = gradt(phi)*val(beta) + val(sigma)*idt(phi);

                bilinearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * L_op * L_opt );

                linearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff*L_op*f );

                break ;
            } //GALS

            case AdvectionStabMethod::SUPG :
            {
                auto coeff = val(vf::h() / (2 * beta_norm + 0.001));

                auto L_op = grad(psi) * val(beta);
                auto L_opt = gradt(phi) * val(beta) + val(sigma) * idt(phi);

                bilinearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * L_op * L_opt );

                linearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff*L_op*f );

                break;
            } //SUPG

            case AdvectionStabMethod::SGS :
            {
                auto coeff = val(1. / ( (2 * beta_norm) / vf::h()  + vf::abs(sigma) ) );

                auto L_op = grad(psi) * val(beta) - val(sigma) * id(psi);
                auto L_opt = gradt(phi) * val(beta) + val(sigma) * idt(phi);
                
                bilinearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff * L_op * L_opt );

                linearForm += integrate(
                        _range=elements(M_mesh),
                        _expr=coeff*L_op*f );

                break;
            } //SGS

            case AdvectionStabMethod::CIP :
            {
                auto coeff = stabCoeff * M_gamma1 * hFace() * hFace() * beta_norm;

                bilinearForm_PatternExtended += integrate(
                        _range=internalfaces(M_mesh),
                        _expr=coeff * jumpt(gradt(phi)) * jump(grad(phi)) );

                break;
            } //CIP

        } //switch
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("Advection","updateLinearPDEStabilization","finish in "+(boost::format("%1% s") %timeElapsed).str() );
}

ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Advection velocity update
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
template<typename ExprT>
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::updateAdvectionVelocity(vf::Expr<ExprT> const& v_expr)
{
    M_advectionVelocity->on(_range=elements(this->mesh()), _expr=v_expr );
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Solve
ADVECTIONBASE_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTIONBASE_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Advection","solve", "start" );
    this->timerTool("Solve").start();

    std::string algebraicSolverType = "LinearSystem";
    M_algebraicFactory->solve( algebraicSolverType, M_blockVectorSolution.vector() ); 
    
    M_blockVectorSolution.localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    this->log("Advection","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

} // namespace FeelModels
} // namespace Feel

/*
//instantiation

// if one application uses an advec() method which as not been instantiated
// before, it has to include <levelsetcore/advect.hpp> which includes this file
// if compiler makes this file when comming from advect.hpp, do not re-instantiate levelset class, thus skip this part
#if !defined( FEELPP_CALLED_FROM_ADVECT_HPP )

#define FEELPP_INSTANTIATE_LEVELSET 1
#define FEELPP_INSTANTIATE_LEVELSET_TEMPLATE_FUNCTION 1
#include "levelset_instance.hpp"
#undef FEELPP_INSTANTIATE_LEVELSET_TEMPLATE_FUNCTION

#endif // FEELPP_CALLED_FROM_ADVECT_HPP

#endif // LEVELSET_ADVECTION_CPP*/
