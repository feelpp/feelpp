#ifndef FMUMODEL2
#define FMUMODEL2

#include "fmilib.h"

#include <boost/filesystem.hpp>

#include "feel/feelcore/environment.hpp"
#include "feel/feelfmi/fmumodelbase.hpp"

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

    std::string name;
    bool is_alias;
    fmi2_value_reference_t ref;
    std::string type, desc;
};


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
    void initialize() override;
    void terminate() override;
    void doStep( double t_cur, double step, bool newStep ) override;

    void printVariablesInfo() override;
    void printInfo() override;

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

private :
    static void fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...);
    static void stepFinished(fmi2_component_environment_t env, fmi2_status_t status);


private :
    fmi2_import_t* M_fmu;
    fmi2_callback_functions_t M_callbackfunctions;
    var_map_type M_v_map;
};


}



#endif
