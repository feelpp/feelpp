#include <feel/feelmor/crbplugin.hpp>

#include <geoheatdiffusion.hpp>

namespace Feel
{

    GeoHeatDiffusion::GeoHeatDiffusion()
        :
        super_type( "GeoHeatDiffusion" )
    {
        this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME) );
        this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
    }

    GeoHeatDiffusion::GeoHeatDiffusion( po::variables_map const& vm )
        :
        super_type( "GeoHeatDiffusion" )
    {
        this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME) );
        this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
    }
    
   FEELPP_CRB_PLUGIN( GeoHeatDiffusion, FEELPP_MOR_PLUGIN_NAME )

}