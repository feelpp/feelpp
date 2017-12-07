#ifndef FMIWRAPPER_HPP
#define FMIWRAPPER_HPP

#include "feel/feelcore/environment.hpp"
#include <boost/filesystem.hpp>
#include "fmilib.h"
#define BUFFER 1000

namespace Feel
{


class FMU
{
public :

    FMU( std::string prefix="" ) :
        M_prefix( prefix ),
        M_verbose( true/*boption(_name="fmu.verbose", _prefix=M_prefix)*/ ),
        M_allocated_dll( false ),
        M_allocated_xml( false )
    {
        M_callbacks.malloc = malloc;
        M_callbacks.calloc = calloc;
        M_callbacks.realloc = realloc;
        M_callbacks.free = free;
        M_callbacks.logger = importlogger;
        M_callbacks.log_level = jm_log_level_all;
        M_callbacks.context = 0;

        M_context = fmi_import_allocate_context( &M_callbacks );
    }

    ~FMU()
    {
        fmi_import_free_context( M_context );
        fmi_import_rmdir( &M_callbacks, tmpDir().c_str() );
        if ( M_allocated_dll )
            fmi2_import_destroy_dllfmu(M_fmu);
        if ( M_allocated_xml )
            fmi2_import_free( M_fmu );
    }

    int load( std::string path )
    {
        // Read fmu version
        M_version = fmi_import_get_fmi_version( M_context, path.c_str(), tmpDir().c_str() );
        if ( M_version == fmi_version_unknown_enu )
        {
            Feel::cerr << "The FMU version could not be determined. "
                       << jm_get_last_error( &M_callbacks ) << std::endl;
            return 0;
        }
        else if ( M_version == fmi_version_2_0_enu )
        {
            // Read informations in xml file for the model
            M_fmu = fmi2_import_parse_xml( M_context, tmpDir().c_str(), 0 );
            if( !M_fmu )
            {
                Feel::cerr << "Error parsing XML. " << jm_get_last_error( &M_callbacks )
                           << std::endl;
                return 0;
            }
            M_allocated_xml=true;

            // Read Model name and ID
            M_model_name = fmi2_import_get_model_name( M_fmu );
            M_guid = fmi2_import_get_GUID( M_fmu );
            printf("Model name : %s with guid %s\n", M_model_name.c_str(), M_guid.c_str());

            // Read model kind : CS or ME
            M_fmukind = fmi2_import_get_fmu_kind( M_fmu );
            if ( M_fmukind == fmi2_fmu_kind_me )
            {
                M_modelID = fmi2_import_get_model_identifier_ME( M_fmu );
                printf("Model identifier for ME: %s\n", M_modelID.c_str() );
            }
            else if ( M_fmukind == fmi2_fmu_kind_cs )
            {
                M_modelID = fmi2_import_get_model_identifier_CS( M_fmu );
                printf("Model identifier for CS: %s\n", M_modelID.c_str() );
            }
            else
            {
                Feel::cerr << "Unxepected FMU kind, exiting. " << jm_get_last_error( &M_callbacks )
                           << std::endl;
                return 0;
            }

            M_callbackfunctions.logger = fmi2logger;
            M_callbackfunctions.allocateMemory = calloc;
            M_callbackfunctions.freeMemory = free;
            M_callbackfunctions.stepFinished = stepFinished;
            M_callbackfunctions.componentEnvironment = 0;

            auto status = fmi2_import_create_dllfmu( M_fmu, M_fmukind, &M_callbackfunctions );
            if ( status==jm_status_error )
            {
                Feel::cerr << "Could not create the DLL loading mechanism(C-API). " << std::endl;
                return 0;
            }
            M_allocated_dll = true;
            printf("Version returned from FMU:   %s\n", fmi2_import_get_version(M_fmu));
        }
        else
        {
            Feel::cerr << "Error only version 2.0 FMU is supported for now\n";
            return 0;
        }

        return 1;
    }

private :

    //! Create Tempory Directory to exctract the FMU in
    std::string tmpDir()
    {
        std::string dir_name = "/ssd/wahl/omcbuild/testsuite/feelom/"+M_prefix+"_fmu";
        boost::filesystem::path dir(dir_name);
        if ( !(boost::filesystem::exists(dir)) )
            boost::filesystem::create_directory(dir);
        return dir_name;
    }

    static void importlogger(jm_callbacks* c,
                             jm_string module, jm_log_level_enu_t log_level, jm_string message)
    {
        printf("FeelFMI module = %s, log level = %d: %s\n", module, log_level, message);
    }


    static void fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...)
    {
        char msg[BUFFER];
        va_list argp;
        va_start(argp, message);

        printf("fmiStatus = %s;  %s (%s): %s\n", fmi2_status_to_string(status), instanceName, category, msg);
    }

    static void stepFinished(fmi2_component_environment_t env, fmi2_status_t status)
    {
        printf("stepFinished is called wiht fmiStatus = %s\n", fmi2_status_to_string(status));

    }


private :
    std::string M_prefix, M_model_name, M_guid, M_modelID;
    bool M_verbose, M_allocated_dll, M_allocated_xml;
    jm_callbacks M_callbacks;
    fmi2_callback_functions_t M_callbackfunctions;
    fmi_import_context_t* M_context;
    fmi2_import_t* M_fmu;
    fmi2_fmu_kind_enu_t M_fmukind;
    fmi_version_enu_t M_version;



};


} //namespace Feel
#endif
