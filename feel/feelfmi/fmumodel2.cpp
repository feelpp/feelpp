#include "feel/feelfmi/fmumodel2.hpp"

namespace Feel
{

FmuModel2::FmuModel2( fmi_import_context_t* context, std::string tmp_dir,
                      callbacks_ptrtype callbacks ) :
    super_type( callbacks )
{
    fmi2_type_t type;

    M_version = 2;
    M_callbackfunctions.logger = fmi2logger;
    M_callbackfunctions.allocateMemory = calloc;
    M_callbackfunctions.freeMemory = free;
    //M_callbackfunctions.stepFinished = stepFinished;
    M_callbackfunctions.componentEnvironment = 0;

    // Read informations in xml file for the model
    M_fmu = fmi2_import_parse_xml( context, tmp_dir.c_str(), 0 );
    CHECK( M_fmu) << "Error parsing XML in tmp dir " << tmp_dir <<", "
                  << jm_get_last_error( callbacks.get() ) << std::endl;
    M_allocated_xml=true;

    // Read Model name and ID
    M_name = fmi2_import_get_model_name( M_fmu );
    M_guid = fmi2_import_get_GUID( M_fmu );

    // Read model kind : CS or ME
    auto kind = fmi2_import_get_fmu_kind( M_fmu );
    if ( kind == fmi2_fmu_kind_me )
    {
        M_kind = "ME";
        M_id = fmi2_import_get_model_identifier_ME( M_fmu );
        type = fmi2_model_exchange;
    }
    else if ( kind == fmi2_fmu_kind_cs )
    {
        M_kind = "CS";
        M_id = fmi2_import_get_model_identifier_CS( M_fmu );
        type = fmi2_cosimulation;
    }
    else
        CHECK( false ) << "Unxepected FMU kind, exiting. " << jm_get_last_error( callbacks.get() ) << std::endl;

    // Load DLL
    auto status = fmi2_import_create_dllfmu( M_fmu, kind, &M_callbackfunctions );
    CHECK( status!=jm_status_error ) << "Could not create the DLL loading mechanism(C-API). " << std::endl;
    M_allocated_dll = true;

    // Instatiate the model
    std::string instance_name = "FMU model " + M_name;
    status = fmi2_import_instantiate( M_fmu, instance_name.c_str(), type, 0, 0 );
    CHECK( status!=jm_status_error ) << "FMUModel2 : fmi2_import_instantiate failed\n";


    LOG(INFO) << "FMU Model loaded : name=" << M_name <<", guid=" << M_guid
              << ", kind="<< M_kind << ", id=" << M_id <<std::endl
              << "Version returned from FMU = " << fmi2_import_get_version(M_fmu) << std::endl;
}

FmuModel2::~FmuModel2()
{
    if ( M_allocated_fmu )
    {
        fmi2_import_terminate(M_fmu);
        fmi2_import_free_instance(M_fmu);
    }
    if ( M_allocated_xml )
        fmi2_import_free( M_fmu );
    if ( M_allocated_dll )
        fmi2_import_destroy_dllfmu(M_fmu);
}


void FmuModel2::initialize()
{
    auto status = fmi2_import_enter_initialization_mode(M_fmu);
    CHECK( status==fmi2_status_ok ) << "FMUModel2 : fmi2_import_enter_initialization_mode failed\n";
    status = fmi2_import_exit_initialization_mode(M_fmu);
    CHECK( status==fmi2_status_ok ) << "FMUModel2 : fmi2_import_exit_initialization_mode failed\n";
}

void FmuModel2::setupExperiment( double const& t_init, double const& t_final, double const& tol )
{
    auto status = fmi2_import_setup_experiment( M_fmu, fmi2_true, tol, t_init, fmi2_false, t_final );
    CHECK( status==fmi2_status_ok )<< "FMUModel2 : fmi2_import_setup_experiment failed\n";
}

void FmuModel2::doStep( double t_cur, double step, bool newStep )
{
    fmi2_import_do_step( M_fmu, t_cur, step, newStep );
}


double FmuModel2::defaultStartTime()
{
    return fmi2_import_get_default_experiment_start( M_fmu );
}

double FmuModel2::defaultFinalTime()
{
    return fmi2_import_get_default_experiment_stop(M_fmu);
}

double FmuModel2::defaultTolerance()
{
    return fmi2_import_get_default_experiment_tolerance(M_fmu);
}


void FmuModel2::fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...)
{
    char buffer[1000];
    va_list args;
    va_start( args, message );
    vsnprintf( buffer, 256, message, args );
    Feel::cout << "FMUModel2 : status="<<fmi2_status_to_string(status)
               << "; " << instanceName <<"("<<category<<") : " << buffer <<std::endl;
    va_end (args);
}

void FmuModel2::stepFinished(fmi2_component_environment_t env, fmi2_status_t status)
{
    Feel::cout << "FMUModel2 : stepFinished called with fmiStatus : "<< fmi2_status_to_string(status) << std::endl;
}

} // namespace Feel
