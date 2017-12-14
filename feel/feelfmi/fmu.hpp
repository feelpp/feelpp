#ifndef FMIWRAPPER_HPP
#define FMIWRAPPER_HPP


#define BUFFER 1000

#include "feel/feelcore/environment.hpp"
#include "fmilib.h"
#include "feel/feelfmi/fmumodel2.hpp"

#include <boost/filesystem.hpp>




namespace Feel
{


class FMU :
        public boost::enable_shared_from_this<FMU>
{
public :
    typedef jm_callbacks callbacks_type;
    typedef boost::shared_ptr<jm_callbacks> callbacks_ptrtype;

    typedef FmuModelBase fmumodel_type;
    typedef boost::shared_ptr<fmumodel_type> fmumodel_ptrtype;

    FMU( std::string prefix="" ) :
        M_prefix( prefix ),
        M_verbose( true/*boption(_name="fmu.verbose", _prefix=M_prefix)*/ ),
        M_allocated_dll( false ),
        M_allocated_xml( false ),
        M_callbacks( new callbacks_type )
    {
        M_callbacks->malloc = malloc;
        M_callbacks->calloc = calloc;
        M_callbacks->realloc = realloc;
        M_callbacks->free = free;
        M_callbacks->logger = importlogger;
        M_callbacks->log_level = jm_log_level_all;
        M_callbacks->context = 0;

        M_context = fmi_import_allocate_context( M_callbacks.get() );
    }

    ~FMU()
    {
        fmi_import_free_context( M_context );
        fmi_import_rmdir( M_callbacks.get(), tmpDir().c_str() );
    }

    int load( std::string path )
    {
        // Read fmu version
        auto version = fmi_import_get_fmi_version( M_context, path.c_str(), tmpDir().c_str() );
        if ( version == fmi_version_unknown_enu )
        {
            Feel::cerr << "The FMU version could not be determined. "
                       << jm_get_last_error( M_callbacks.get() ) << std::endl;
            return 0;
        }
        else if ( version == fmi_version_2_0_enu )
        {
            M_model = fmumodel_ptrtype( new FmuModel2( M_context, tmpDir(), M_callbacks ) );

        }
        else
        {
            Feel::cerr << "Error only version 2.0 FMU is supported for now\n";
            return 0;
        }

        return 1;
    }


    void simulate( double t_init=-1, double t_final=-1 )
    {
        // if ( t_init==-1 )
        //     t_init = defaultStartTime();
        // if ( t_final==-1 )
        //     t_final = defaultFinalTime();
    }

    // //! Get the default initial time of the model
    // double defaultStartTime() { return fmi2_import_get_default_experiment_start(M_fmu); }
    // //! Get the default final time of the model
    // double defaultFinalTime() { return fmi2_import_get_default_experiment_stop(M_fmu); }
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




private :
    std::string M_prefix;
    bool M_verbose, M_allocated_dll, M_allocated_xml;
    callbacks_ptrtype M_callbacks;
    fmi_import_context_t* M_context;

    fmumodel_ptrtype M_model;
};


} //namespace Feel
#endif
