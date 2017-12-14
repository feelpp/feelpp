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
        M_allocated_dll( false )
        {}

    int version() { return M_version; }
    std::string name() { return M_name; }
    std::string guid() { return M_guid; }

protected :
    callbacks_ptrtype M_callbacks;
    bool M_allocated_xml, M_allocated_dll;
    std::string M_name, M_guid, M_id, M_kind;
    int M_version;
};

}

#endif
