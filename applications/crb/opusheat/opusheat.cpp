
#include <opusheat_impl.cpp>

#include <feel/feelcrb/crbplugin.hpp>

namespace Feel
{
FEELPP_CRB_PLUGIN_TEMPLATE( OpusHeat, OpusHeat<false>, opusheat )

template class OpusHeat<false>;
}

