#include <opusheat_impl.cpp>

#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{
FEELPP_CRB_PLUGIN_TEMPLATE( OpusHeatStationary, OpusHeat<true>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,P1G1) )

template class OpusHeat<true>;
}
