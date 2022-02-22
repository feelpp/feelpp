
#include <opusheat_impl.cpp>

#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{
FEELPP_CRB_PLUGIN_TEMPLATE( OpusHeat, OpusHeat<false>, opusheat )

template class OpusHeat<false>;
}

