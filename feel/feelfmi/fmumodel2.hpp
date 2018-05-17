#ifndef FMUMODEL2
#define FMUMODEL2

#include "fmilib.h"

#include <boost/filesystem.hpp>

#include "feel/feelcore/environment.hpp"
#include "feel/feelfmi/fmumodelbase.hpp"

namespace Feel
{
struct Fmi2Variable;

class FmuModel2 :
        public FmuModelBase
{
public :
    typedef FmuModel2 self_type;
    typedef FmuModelBase super_type;

    typedef typename super_type::callbacks_ptrtype callbacks_ptrtype;
    typedef Fmi2Variable variable_type;
    typedef boost::shared_ptr<variable_type> variable_ptrtype;
    typedef std::map<std::string,variable_ptrtype> var_map_type;

    FmuModel2( fmi_import_context_t* context, std::string tmp_dir, callbacks_ptrtype callbacks );
    ~FmuModel2();

    void reset() override;
    void setupExperiment( double const& t_init, double const& t_final, double const& tol ) override;
    void initialize( double t_init ) override;
    void terminate() override;
    void doStep( double t_cur, double step, bool newStep ) override;

    void printVariablesInfo() override;
    void printInfo() override;
    void exportValues() override;

    double defaultStartTime() override;
    double defaultFinalTime() override;
    double defaultTolerance() override;

    void setValue( std::string name, double value ) override;
    void setValue( std::string name, int value ) override;
    void setValue( std::string name, std::string value ) override;
    void setValue( std::string name, bool value ) override;

    void getValue( std::string name, double& value ) override;
    void getValue( std::string name, int& value ) override;
    void getValue( std::string name, std::string& value ) override;
    void getValue( std::string name, bool& value ) override;

    template <typename VariableType>
    VariableType getValue( std::string name )
    {
        VariableType value;
        getValue( name, value );
        return value;
    }

private :
    static void fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...);
    static void stepFinished(fmi2_component_environment_t env, fmi2_status_t status);

    std::string strValue( std::string var );

private :
    fmi2_import_t* M_fmu;
    fmi2_callback_functions_t M_callbackfunctions;
    var_map_type M_v_map;
    std::map< std::string,std::vector<std::string> > M_values;
    fmi2_type_t M_type;

};


}



#endif
