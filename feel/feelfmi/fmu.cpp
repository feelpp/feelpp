#include "feel/feelfmi/fmu.hpp"

namespace Feel
{

FMU::FMU( std::string prefix ) :
    M_prefix( prefix ),
    M_verbose( boption(_name="fmu.verbose", _prefix=M_prefix) ),
    M_callbacks( new callbacks_type )
{
    M_callbacks->malloc = malloc;
    M_callbacks->calloc = calloc;
    M_callbacks->realloc = realloc;
    M_callbacks->free = free;
    M_callbacks->logger = importlogger;
    M_callbacks->log_level = jm_log_level_all;//FLAGS_v;
    M_callbacks->context = 0;

    M_context = fmi_import_allocate_context( M_callbacks.get() );

    // create temporary repository to extract fmu
    M_tmp_dir = Environment::appRepository() + "/fmu_tmpdir/";
    if ( M_prefix!="" )
        M_tmp_dir += M_prefix+"/";
    boost::filesystem::path dir( M_tmp_dir );
    if ( !(boost::filesystem::exists(dir)) )
        boost::filesystem::create_directory(dir);
}


FMU::~FMU()
{
    fmi_import_free_context( M_context );
    fmi_import_rmdir( M_callbacks.get(), M_tmp_dir.c_str() );
}


int FMU::load( std::string _path )
{
    std::string path = _path;
    if ( path=="" )
        path = Environment::expand( soption( _name="fmu.filename", _prefix=M_prefix ) );
    CHECK( path!="" ) << "No filename specified to load FMU. Either pass argument to the function load() or set the option fmu.filename\n";

    // Read fmu version
    auto version = fmi_import_get_fmi_version( M_context, path.c_str(), M_tmp_dir.c_str() );
    CHECK( version != fmi_version_unknown_enu ) << "The FMU version could not be determined. "
                                                << jm_get_last_error( M_callbacks.get() )
                                                << std::endl;

    if ( version == fmi_version_2_0_enu )
        M_model = fmumodel_ptrtype( new FmuModel2( M_context, M_tmp_dir, M_callbacks ) );
    else
    {
        Feel::cerr << "Error only version 2.0 FMU is supported for now\n";
        return 0;
    }

    if ( boption(_name="fmu.display-variables-info", _prefix=M_prefix) )
        M_model->printVariablesInfo();

    return 1;
}


void FMU::simulate( double t_init, double t_final, double tolerance )
{
    this->initialize( t_init, t_final, tolerance );
    M_solver->simulate();
}

void FMU::initialize( double t_init, double t_final, double tolerance )
{
    if ( !M_solver )
        initSolver();
    if ( t_init==-1 )
        t_init = M_model->defaultStartTime();
    if ( t_final==-1 )
        t_final = M_model->defaultFinalTime();
    if ( tolerance==-1 )
        tolerance = M_model->defaultTolerance();

    M_solver->initialize( t_init, t_final, tolerance );
}

void FMU::doSteps( double t_stop )
{
    M_solver->doSteps( t_stop );
}

void FMU::reset()
{
    CHECK( M_model ) <<"FMU : error in reset(). No FMU loaded\n";
    M_model->reset();
}

void FMU::initSolver()
{
    CHECK( M_model ) <<"FMU : error in initSolver. No FMU loaded\n";
    auto kind = M_model->kind();
    if ( kind=="CS" )
        M_solver = solver_ptrtype( new CSSolver( M_model ) );

    M_solver->setTimeStep( doption( _name="fmu.solver.time-step", _prefix=M_prefix ) );
    M_solver->setRelativeTol( doption( _name="fmu.solver.rtol", _prefix=M_prefix) );
}


void FMU::setSolverTimeStep( double const& step )
{
    if ( !M_solver )
        initSolver();
    M_solver->setTimeStep( step );
}


void FMU::printModelInfo()
{
    CHECK( M_model ) <<"FMU : error in printModelInfo. No FMU loaded\n";
    M_model->printInfo();
    M_model->printVariablesInfo();
}

double FMU::currentTime()
{
    if ( !M_solver )
        initSolver();
    return M_solver->currentTime();
}

} //namespace Feel
