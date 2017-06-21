
#include <opusheat_impl.cpp>

#include <feel/feelcrb/crbplugin.hpp>

namespace Feel
{
FEELPP_CRB_PLUGIN_TEMPLATE( OpusHeatStationary, OpusHeat<true>, "opusheat-stationary" )

template class OpusHeat<true>;
}

