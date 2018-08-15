#ifndef FMUMODELBASE_HPP
#define FMUMODELBASE_HPP

#include "fmilib.h"

namespace Feel
{
class FmuModelBase
{
public :
    typedef jm_callbacks callbacks_type;
    typedef std::shared_ptr<jm_callbacks> callbacks_ptrtype;

    typedef std::vector<std::string> var_list_type;
    typedef std::shared_ptr<var_list_type> var_list_ptrtype;

    FmuModelBase( callbacks_ptrtype callbacks ) :
        M_callbacks( callbacks ),
        M_allocated_xml( false ),
        M_allocated_dll( false ),
        M_allocated_fmu( false ),
        M_setup( false )
    {}

    virtual ~FmuModelBase()
    {}

    void setExportList( var_list_ptrtype exp_list ) { M_export_list = exp_list; }
    void setExportDirectory( std::string const& path ) { M_export_directory=path; }

    int version() { return M_version; }
    std::string name() { return M_name; }
    std::string guid() { return M_guid; }
    std::string kind() { return M_kind; }
    bool isSetup() { return M_setup; }

    virtual void reset()=0;
    virtual void setupExperiment( double const& t_init, double const& t_final, double const& tol )=0;
    virtual void initialize( double t_init )=0;
    virtual void terminate()=0;
    virtual void doStep( double t_cur, double step, bool newStep )=0;

    virtual void printInfo()=0;
    virtual void printVariablesInfo()=0;
    virtual void exportValues()=0;

    virtual double defaultStartTime()=0;
    virtual double defaultFinalTime()=0;
    virtual double defaultTolerance()=0;

    virtual void setValue( std::string name, double value )=0;
    virtual void setValue( std::string name, int value )=0;
    virtual void setValue( std::string name, std::string value )=0;
    virtual void setValue( std::string name, bool value )=0;

    virtual void getValue( std::string name, double& value )=0;
    virtual void getValue( std::string name, int& value )=0;
    virtual void getValue( std::string name, std::string& value )=0;
    virtual void getValue( std::string name, bool& value )=0;

protected :
    callbacks_ptrtype M_callbacks;
    bool M_allocated_xml, M_allocated_dll, M_allocated_fmu, M_setup;
    std::string M_name, M_guid, M_id, M_kind, M_export_directory;
    int M_version;

    var_list_ptrtype M_export_list;
};

}

#endif
