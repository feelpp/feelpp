#ifndef FMIWRAPPER_HPP
#define FMIWRAPPER_HPP

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

    FMU( std::string prefix="" );
    ~FMU();

    int load( std::string _path="" );
    void simulate( double t_init=-1, double t_final=-1, double tolerance=-1 );
    void reset();
    void initialize( double t_init=-1, double t_final=-1, double tolerance=-1 );
    void doSteps( double t_stop );

    void initSolver();
    void setSolverTimeStep( double const& step );

    void printModelInfo();

    double currentTime();


    template <typename VariableType>
    void setValue( std::string var_name, VariableType value )
    {
        CHECK( M_model ) <<"FMU trying to access variable without model\n";
        M_model->setValue( var_name, value );
    }

    template <typename VariableType>
    VariableType getValue( std::string var_name )
    {
        VariableType value;
        CHECK( M_model ) <<"FMU trying to access variable without model\n";
        M_model->getValue( var_name, value );
        return value;
    }

private :
    static void importlogger(jm_callbacks* c,
                             jm_string module, jm_log_level_enu_t log_level, jm_string message)
    {
        LOG(INFO) << "FeelFMI["<< module <<"] : " << message << std::endl;
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
