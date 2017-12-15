#ifndef FMIWRAPPER_HPP
#define FMIWRAPPER_HPP


#define BUFFER 1000

#include "feel/feelcore/environment.hpp"

#include "fmilib.h"
#include <boost/filesystem.hpp>


#include "feel/feelfmi/fmumodel2.hpp"
#include "feel/feelfmi/cssolver.hpp"

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

    typedef SolverBase solver_type;
    typedef boost::shared_ptr<solver_type> solver_ptrtype;


    FMU( std::string prefix="" ) :
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
        M_tmp_dir = Environment::rootRepository() + "/fmu_tmpdir/";
        if ( M_prefix!="" )
            M_tmp_dir += M_prefix+"/";
        boost::filesystem::path dir( M_tmp_dir );
        if ( !(boost::filesystem::exists(dir)) )
            boost::filesystem::create_directory(dir);
    }

    ~FMU()
    {
        fmi_import_free_context( M_context );
        fmi_import_rmdir( M_callbacks.get(), M_tmp_dir.c_str() );
    }

    int load( std::string _path="" )
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

        return 1;
    }


    void simulate( double t_init=-1, double t_final=-1, double tolerance=-1 )
    {
        if ( !M_solver )
            initSolver();
        if ( t_init==-1 )
            t_init = M_model->defaultStartTime();
        if ( t_final==-1 )
            t_final = M_model->defaultFinalTime();
        if ( tolerance==-1 )
            tolerance = M_model->defaultTolerance();

        M_solver->simulate( t_init, t_final, tolerance );
    }

    void initSolver()
    {
        auto kind = M_model->kind();
        if ( kind=="CS" )
            M_solver = solver_ptrtype( new CSSolver( M_model ) );

        M_solver->setTimeStep( doption( _name="fmu.solver.time-step", _prefix=M_prefix ) );
        M_solver->setRelativeTol( doption( _name="fmu.solver.rtol", _prefix=M_prefix) );
    }

    void setSolverTimeStep( double const& step )
    {
        if ( !M_solver )
            initSolver();
        M_solver->setTimeStep( step );
    }

private :
    static void importlogger(jm_callbacks* c,
                             jm_string module, jm_log_level_enu_t log_level, jm_string message)
    {
        Feel::cout << "FeelFMI["<< module <<"] : " << message << std::endl;
    }



private :
    std::string M_prefix, M_tmp_dir;
    bool M_verbose;
    callbacks_ptrtype M_callbacks;
    fmi_import_context_t* M_context;

    fmumodel_ptrtype M_model;
    solver_ptrtype M_solver;
};



} //namespace Feel
#endif
