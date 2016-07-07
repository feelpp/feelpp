
#include "multifluidconfig.h"

#include <feel/feelmodels/multifluid/multifluid.cpp>
#include <feel/feelmodels/multifluid/helfrichforcemodel.hpp>

namespace Feel {
namespace FeelModels {

//template class MultiFluid< FLUIDMECHANICS_CLASS_TYPE, LEVELSET_CLASS_TYPE >;
template class MULTIFLUID_CLASS_INSTANTIATION;

// Register interface forces models
const bool helfrich_interfaceforcesmodel = 
    MULTIFLUID_CLASS_INSTANTIATION::interfaceforces_factory_type::instance().registerProduct( 
            "helfrich", 
            &detail::createInterfaceForcesModel<HelfrichForceModel, typename MULTIFLUID_CLASS_INSTANTIATION::levelset_type> );

}
}

