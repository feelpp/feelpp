#ifndef FMUMODELBASE_HPP
#define FMUMODELBASE_HPP

#include "fmilib.h"

namespace Feel
{
class FmuModelBase
{
public :
    typedef jm_callbacks callbacks_type;
    typedef boost::shared_ptr<jm_callbacks> callbacks_ptrtype;

    FmuModelBase( callbacks_ptrtype callbacks ) :
        M_callbacks( callbacks ),
        M_allocated_xml( false ),
        M_allocated_dll( false ),
        M_allocated_fmu( false )
    {}

    virtual ~FmuModelBase()
    {}

    virtual void setupExperiment( double const& t_init, double const& t_final, double const& tol )=0;
    virtual void initialize()=0;
    virtual void doStep( double t_cur, double step, bool newStep )=0;

    int version() { return M_version; }
    std::string name() { return M_name; }
    std::string guid() { return M_guid; }
    std::string kind() { return M_kind; }

    virtual double defaultStartTime()=0;
    virtual double defaultFinalTime()=0;
    virtual double defaultTolerance()=0;

protected :
    callbacks_ptrtype M_callbacks;
    bool M_allocated_xml, M_allocated_dll, M_allocated_fmu;
    std::string M_name, M_guid, M_id, M_kind;
    int M_version;
};

}

#endif
