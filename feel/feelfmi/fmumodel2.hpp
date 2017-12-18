#ifndef FMUMODEL2
#define FMUMODEL2

#include "fmilib.h"

#include <boost/filesystem.hpp>

#include "feel/feelcore/environment.hpp"
#include "feel/feelfmi/fmumodelbase.hpp"

namespace Feel
{

class FmuModel2 :
        public FmuModelBase
{
public :
    typedef FmuModel2 self_type;
    typedef FmuModelBase super_type;

    typedef typename super_type::callbacks_ptrtype callbacks_ptrtype;

    FmuModel2( fmi_import_context_t* context, std::string tmp_dir, callbacks_ptrtype callbacks );
    ~FmuModel2();

    void setupExperiment( double const& t_init, double const& t_final, double const& tol ) override;
    void initialize() override;
    void doStep( double t_cur, double step, bool newStep ) override;

    double defaultStartTime() override;
    double defaultFinalTime() override;
    double defaultTolerance() override;

private :
    static void fmi2logger(fmi2_component_environment_t env, fmi2_string_t instanceName,
                           fmi2_status_t status, fmi2_string_t category, fmi2_string_t message, ...);
    static void stepFinished(fmi2_component_environment_t env, fmi2_status_t status);


private :
    fmi2_import_t* M_fmu;
    fmi2_callback_functions_t M_callbackfunctions;
};


}



#endif
