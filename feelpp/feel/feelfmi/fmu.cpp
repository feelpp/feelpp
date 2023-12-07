#include "feel/feelfmi/fmu.hpp"

namespace Feel
{

FMU::FMU( std::string prefix ) :
    M_prefix( prefix ),
    M_verbose( boption(_name="fmu.verbose", _prefix=M_prefix) ),
    M_callbacks( new callbacks_type ),
    M_tinit( doption( _name="fmu.time-initial", _prefix=M_prefix ) ),
    M_tfinal( doption( _name="fmu.time-final", _prefix=M_prefix ) )
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
    fs::path dir( M_tmp_dir );
    if ( !(fs::exists(dir)) )
        fs::create_directory(dir);

    M_export_directory = soption( _name="fmu.export-directory", _prefix=M_prefix);
    if ( M_export_directory!="" )
        M_export_directory = Environment::expand( M_export_directory );

    if ( Environment::vm().count(prefixvm( M_prefix,"fmu.exported-variables" )) )
    {
        M_export_list = std::make_shared<var_list_type>( option( _name="fmu.exported-variables", _prefix=M_prefix).template as<std::vector<std::string> >() );
    }
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
    {
        std::string fmuname = Environment::findFileRemotely(soption( _name="fmu.filename", _prefix=M_prefix ), "fmus"  );
        path = Environment::expand( fmuname );
    }

    CHECK( path!="" ) << "No filename specified to load FMU. Either pass argument to the function load() or set the option fmu.filename\n";

    // Read fmu version
    auto version = fmi_import_get_fmi_version( M_context, path.c_str(), M_tmp_dir.c_str() );
    CHECK( version != fmi_version_unknown_enu ) << "The FMU version could not be determined. "
                                                << jm_get_last_error( M_callbacks.get() )
                                                << std::endl;

    if ( version == fmi_version_2_0_enu )
        M_model = std::make_shared<FmuModel2>( M_context, M_tmp_dir, M_callbacks );
    else
    {
        Feel::cerr << "Error only version 2.0 FMU is supported for now\n";
        return 0;
    }

    if ( boption(_name="fmu.display-variables-info", _prefix=M_prefix) )
        M_model->printVariablesInfo();
    if ( M_export_list )
    {
        M_model->setExportList( M_export_list );
        M_model->setExportDirectory( M_export_directory );
    }

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
    if ( t_init>=0 )
        M_tinit = t_init;
    else if ( M_tinit<0 )
        M_tinit = M_model->defaultStartTime();

    if ( t_final>0 )
        M_tfinal = t_final;
    else if ( M_tfinal<0 )
        M_tfinal = M_model->defaultFinalTime();

    if ( tolerance==-1 )
        tolerance = M_model->defaultTolerance();

    M_solver->initialize( M_tinit, M_tfinal, tolerance );
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

void FMU::addExportedVariables( std::vector<std::string> const& var_list )
{
    if ( !M_export_list )
        setExportedVariables( var_list );
    else
    {
        M_export_list->insert( M_export_list->end(), var_list.begin(), var_list.end() );
        M_model->setExportList( M_export_list );
        M_model->setExportDirectory( M_export_directory );
    }
}

void FMU::setExportedVariables( std::vector<std::string> const& var_list )
{
    M_export_list = std::make_shared<var_list_type>( var_list );
    M_model->setExportList( M_export_list );
    M_model->setExportDirectory( M_export_directory );
}

void FMU::setExportDirectory( std::string const& path )
{
    M_export_directory = path;
    M_model->setExportDirectory( M_export_directory );
}

double FMU::currentTime()
{
    if ( !M_solver )
        initSolver();
    return M_solver->currentTime();
}

} //namespace Feel
