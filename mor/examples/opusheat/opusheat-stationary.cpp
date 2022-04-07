
#include <opusheat_impl.cpp>

#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{
FEELPP_CRB_PLUGIN_TEMPLATE( OpusHeatStationary, OpusHeat<true>, opusheat_stationary )

template class OpusHeat<true>;
}

