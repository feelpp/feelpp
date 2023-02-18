#include "feel/feelfmi/fmumodel2.hpp"

namespace Feel
{

struct Fmi2Variable
{
    Fmi2Variable( fmi2_import_variable_t* v )
    {
        name = fmi2_import_get_variable_name( v );
        is_alias = fmi2_import_get_variable_alias_kind( v );
        ref = fmi2_import_get_variable_vr( v );

        auto base_type = fmi2_import_get_variable_base_type( v );
        if ( base_type==fmi2_base_type_real )
            type="double";
        else if ( base_type==fmi2_base_type_int )
            type="int";
        else if ( base_type==fmi2_base_type_bool )
            type="bool";
        else if ( base_type==fmi2_base_type_str )
            type="string";
        else if ( base_type==fmi2_base_type_enum )
            type="enum";
        auto s = fmi2_import_get_variable_description( v );
        if ( s )
            desc = s;
    }

    void print( bool with_header=false ) const
    {
        if ( with_header )
            Feel::cout << std::left<< std::setw(30) << "Name" <<std::setw(6) <<"Alias"
                       << std::setw(10) <<"Type" << "Description" << std::endl;
        Feel::cout <<std::left<< std::setw(30) <<  name <<std::setw(6) <<is_alias << std::setw(10)
                   <<type << desc << std::endl;
    }

    std::string name, value;
    bool is_alias;
    fmi2_value_reference_t ref;
    std::string type, desc;
};



FmuModel2::FmuModel2( fmi_import_context_t* context, std::string tmp_dir,
                      callbacks_ptrtype callbacks ) :
    super_type( callbacks )
{

    M_version = 2;
    M_callbackfunctions.logger = fmi2logger;
    M_callbackfunctions.allocateMemory = calloc;
    M_callbackfunctions.freeMemory = free;
    //M_callbackfunctions.stepFinished = stepFinished;
    M_callbackfunctions.componentEnvironment = 0;

    // Read information in xml file for the model
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
        M_type = fmi2_model_exchange;
    }
    else if ( kind == fmi2_fmu_kind_cs )
    {
        M_kind = "CS";
        M_id = fmi2_import_get_model_identifier_CS( M_fmu );
        M_type = fmi2_cosimulation;
    }
    else
        CHECK( false ) << "Unexpected FMU kind, exiting. " << jm_get_last_error( callbacks.get() ) << std::endl;

    // Load DLL
    auto status = fmi2_import_create_dllfmu( M_fmu, kind, &M_callbackfunctions );
    CHECK( status!=jm_status_error ) << "Could not create the DLL loading mechanism(C-API). " << std::endl;
    M_allocated_dll = true;

    // Build the map of variables
    auto v_list = fmi2_import_get_variable_list( M_fmu, 1 );
    int N = fmi2_import_get_variable_list_size( v_list );
    for ( int n=0; n<N; n++ )
    {
        auto v = fmi2_import_get_variable( v_list, n );
        auto var = std::make_shared<variable_type>(v);
        M_v_map[var->name] = var;
    }

    LOG(INFO) << "FMU Model loaded : name=" << M_name <<", guid=" << M_guid
              << ", kind="<< M_kind << ", id=" << M_id <<std::endl
              << "Version returned from FMU = " << fmi2_import_get_version(M_fmu) << std::endl;
}

FmuModel2::~FmuModel2()
{
    if ( M_setup )
        terminate();
    if ( M_allocated_fmu )
        fmi2_import_free_instance(M_fmu);
    if ( M_allocated_xml )
        fmi2_import_free( M_fmu );
    if ( M_allocated_dll )
        fmi2_import_destroy_dllfmu(M_fmu);
}

void FmuModel2::reset()
{
    terminate();
    auto status = fmi2_import_reset(M_fmu);
    CHECK( status==fmi2_status_ok ) << "FMUModel2 : fmi2_import_reset failed\n";
    M_allocated_fmu=false;
    if( M_export_list )
        M_values.clear();
    M_setup=false;
}

void FmuModel2::initialize( double t_init )
{
    // Instatiate the model
    auto status = fmi2_import_enter_initialization_mode(M_fmu);
    CHECK( status==fmi2_status_ok ) << "FMUModel2 : fmi2_import_enter_initialization_mode failed\n";
    status = fmi2_import_exit_initialization_mode(M_fmu);
    CHECK( status==fmi2_status_ok ) << "FMUModel2 : fmi2_import_exit_initialization_mode failed\n";

    if( M_export_list )
    {
        std::vector<std::string> t = { std::to_string(t_init) };
        M_values["t"] = t;
        for( auto const& var : *M_export_list )
        {
            auto it = M_v_map.find(var);
            if ( it==M_v_map.end() )
                Feel::cout << "WARNING: you asked to export the variable "<< var <<", but there is no corresponding variable.\n";

            std::vector<std::string> t = { strValue(var) };
            M_values[var] = t;
        }
    }
}

void FmuModel2::terminate()
{
    auto status = fmi2_import_terminate(M_fmu);
    CHECK( status==fmi2_status_ok ) << "FMUModel2 : fmi2_import_terminate failed\n";
}

void FmuModel2::setupExperiment( double const& t_init, double const& t_final, double const& tol )
{
    std::string instance_name = "FMU model " + M_name;
    auto status1 = fmi2_import_instantiate( M_fmu, instance_name.c_str(), M_type, 0, 0 );
    CHECK( status1!=jm_status_error ) << "FMUModel2 : fmi2_import_instantiate failed\n";
    M_allocated_fmu=true;

    auto status = fmi2_import_setup_experiment( M_fmu, fmi2_true, tol, t_init, false, t_final );
    CHECK( status==fmi2_status_ok )<< "FMUModel2 : fmi2_import_setup_experiment failed\n";
    M_setup=true;
}

void FmuModel2::doStep( double t_cur, double step, bool newStep )
{
    fmi2_import_do_step( M_fmu, t_cur, step, newStep );

    if ( M_export_list )
    {
        M_values["t"].push_back( std::to_string(t_cur+step) );
        for( auto const& var : *M_export_list )
            M_values[var].push_back( strValue(var) );
    }
}

void FmuModel2::printInfo() const
{
    Feel::cout << "FMU Model loaded : name=" << M_name <<", guid=" << M_guid
               << ", kind="<< M_kind << ", id=" << M_id <<std::endl
               << "Version returned from FMU = " << fmi2_import_get_version(M_fmu) << std::endl;
}

void FmuModel2::printVariablesInfo() const
{
    Feel::cout << "=================================================================\n"
               << "FMU Model "<< name() << " : list of variables\n"
               << "=================================================================\n";
    Feel::cout << std::left<< std::setw(30) << "Name" <<std::setw(6) <<"Alias" << std::setw(10)
               <<"Type" << "Description" << std::endl;
    for ( auto const& v : M_v_map )
        v.second->print();
    Feel::cout << "=================================================================\n";
}


void FmuModel2::exportValues()  const
{
    if ( Environment::isMasterRank() && M_export_list )
    {
        std::string path = M_export_directory=="" ? "fmu_values.csv": M_export_directory + "/fmu_values.csv";
        std::ofstream data;
        data.open( path , std::ios::trunc );
        data << "t";
        for( auto const& var : *M_export_list )
            data << "," << var;
        data << std::endl;

        auto n_step = M_values.at("t").size();
        for ( int i=0; i<n_step; i++ )
        {
            data << M_values.at("t")[i];
            for( auto const& var : *M_export_list )
            {
                if ( M_values.at(var).size()>i )
                    data<<","<<M_values.at(var)[i];
                else
                    data <<",NA";
            }
            data << std::endl;
        }
        data.close();
    }

}

double FmuModel2::defaultStartTime() const
{
    return fmi2_import_get_default_experiment_start( M_fmu );
}

double FmuModel2::defaultFinalTime() const
{
    return fmi2_import_get_default_experiment_stop(M_fmu);
}

double FmuModel2::defaultTolerance() const
{
    return fmi2_import_get_default_experiment_tolerance(M_fmu);
}

void FmuModel2::setValue( std::string name, double value )
{
    CHECK( M_v_map[name]->type=="double" ) <<"FMIModel2 : setValue() called with wrong type double when type of the variable is set as " << M_v_map[name]->type <<std::endl;
    auto ref = M_v_map[name]->ref;
    auto status = fmi2_import_set_real( M_fmu, &ref, 1, &value );
}
void FmuModel2::setValue( std::string name, int value )
{
    CHECK( M_v_map[name]->type=="int" ) <<"FMIModel2 : setValue() called with wrong type int when type of the variable is set as " << M_v_map[name]->type <<std::endl;
    auto ref = M_v_map[name]->ref;
    auto status = fmi2_import_set_integer( M_fmu, &ref, 1, &value );
}
void FmuModel2::setValue( std::string name, std::string value )
{
    CHECK( M_v_map[name]->type=="string" ) <<"FMIModel2 : setValue() called with wrong type string when type of the variable is set as " << M_v_map[name]->type <<std::endl;
    auto ref = M_v_map[name]->ref;
    fmi2_string_t str = value.c_str();
    auto status = fmi2_import_set_string( M_fmu, &ref, 1, &str );
}
void FmuModel2::setValue( std::string name, bool value )
{
    CHECK( M_v_map[name]->type=="bool" ) <<"FMIModel2 : setValue() called with wrong type bool when type of the variable is set as " << M_v_map[name]->type <<std::endl;
    auto ref = M_v_map[name]->ref;
    fmi2_boolean_t b = value;
    auto status = fmi2_import_set_boolean( M_fmu, &ref, 1, &b );
}

void FmuModel2::getValue( std::string name, double& value ) const
{
    auto ref = M_v_map.at(name)->ref;
    auto status = fmi2_import_get_real( M_fmu, &ref, 1, &value );
}

void FmuModel2::getValue( std::string name, int& value ) const
{
    auto ref = M_v_map.at(name)->ref;
    auto status = fmi2_import_get_integer( M_fmu, &ref, 1, &value );
}

void FmuModel2::getValue( std::string name, std::string& value ) const
{
    auto ref = M_v_map.at(name)->ref;
    fmi2_string_t str;
    auto status = fmi2_import_get_string( M_fmu, &ref, 1, &str );
    value = str;
}

void FmuModel2::getValue( std::string name, bool& value ) const
{
    auto ref = M_v_map.at(name)->ref;
    fmi2_boolean_t b;
    auto status = fmi2_import_get_boolean( M_fmu, &ref, 1, &b );
    value = b;
}



void FmuModel2::fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...)
{
    char buffer[1000];
    va_list args;
    va_start( args, message );
    vsnprintf( buffer, 256, message, args );
    LOG(INFO) << "FMUModel2 : status="<<fmi2_status_to_string(status)
               << "; " << instanceName <<"("<<category<<") : " << buffer <<std::endl;
    va_end (args);
}

void FmuModel2::stepFinished(fmi2_component_environment_t env, fmi2_status_t status)
{
    LOG(INFO) << "FMUModel2 : stepFinished called with fmiStatus : "<< fmi2_status_to_string(status) << std::endl;
}

std::string FmuModel2::strValue( std::string const& var )
{
    std::stringstream ss;
    auto it = M_v_map.find(var);
    if ( it!=M_v_map.end() )
    {
        auto type = M_v_map[var]->type;

        if ( type=="int" )
            ss << getValue<int>( var );
        else if ( type=="double" )
            ss << getValue<double>( var );
        else if ( type=="bool" )
            ss << getValue<bool>( var );
        else if ( type=="string" )
            ss << getValue<std::string>( var );
    }
    else
        ss << "NA";

    return ss.str();
}


} // namespace Feel
