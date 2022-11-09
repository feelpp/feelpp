#include "toolboxmor_heat.hpp"

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{

template<typename ToolboxType, int Options>
ToolboxMorHeat<ToolboxType,Options>::ToolboxMorHeat( std::string const& name, std::string const& prefix )
    :
    super_type( name, prefix )
{
    this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME) + fmt::format("_{}dP{}G{}",toolbox_type::nDim, toolbox_type::nOrderTemperature, toolbox_type::nOrderGeo) );
    this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
}

template<typename ToolboxType, int Options>
void
ToolboxMorHeat<ToolboxType,Options>::initModel()
{
    M_offlineToolbox = toolbox_type::New(_prefix="heat");
    M_offlineToolbox->init();
    M_offlineToolbox->printAndSaveInfo();
    if ( false )
        M_offlineToolbox->solve();

    this->setFunctionSpaces( M_offlineToolbox->spaceTemperature() );

    if ( M_offlineToolbox->hasModelProperties() )
        this->addModelData( "toolbox_json_setup", M_offlineToolbox->modelProperties().jsonData(), "toolbox_model/setup.json" );

    auto heatBoxModel = DeimMorModelToolbox<toolbox_type>::New( M_offlineToolbox );
    this->initOfflineToolbox( heatBoxModel );

    //this->initOnline();// maybe give here heatBoxModel
}

template<typename ToolboxType, int Options>
void
ToolboxMorHeat<ToolboxType,Options>::initOnlineToolbox( std::shared_ptr<DeimMorModelBase<typename super_type::mesh_type>> heatBoxModel )
{
    //auto heatBoxModel = DeimMorModelToolbox<toolbox_type>::New("heat");

    M_onlineModelProperties = std::make_shared<ModelProperties>();// this->repository().expr(),this->worldCommPtr(), this->prefix(), this->clovm() );
    if( this->hasModelData("toolbox_json_setup") )
    {
        auto & mdata = this->additionalModelData("toolbox_json_setup");
        auto const& jsonData = mdata.template fetch_data<nl::json>( this->crbModelDb().dbRepository() );
        M_onlineModelProperties->setup( jsonData );
    }

    std::dynamic_pointer_cast<DeimMorModelToolbox<toolbox_type>>( heatBoxModel )->setToolboxInitFunction(
        [this](   /*auto*/ typename toolbox_type::mesh_ptrtype  mesh ) {
            auto tbDeim = toolbox_type::New( _prefix="heat"/*M_prefix*/);
            tbDeim->setModelProperties( M_onlineModelProperties );
            tbDeim->setMesh(mesh);
            tbDeim->init();
            //tbDeim->printAndSaveInfo();
            return tbDeim;
        });

    super_type::initOnlineToolbox(heatBoxModel);
}

template<typename ToolboxType, int Options>
void
ToolboxMorHeat<ToolboxType,Options>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    super_type::setupSpecificityModel( ptree,dbDir );
    auto heatBoxModel = DeimMorModelToolbox<toolbox_type>::New("heat");
    this->initOnlineToolbox( heatBoxModel );
    //this->initOnline();
}


template <int Dim,int Order>
struct heat_type
{
    using convex_type = Simplex<Dim>;
    using base_type = Lagrange<Order, Scalar, Continuous, PointSetFekete>;
    using type = FeelModels::Heat<convex_type, base_type>;
};

using heat_2dP1G1_t = typename heat_type<2,1>::type;
using heat_3dP1G1_t = typename heat_type<3,1>::type;

template class ToolboxMorHeat<heat_2dP1G1_t>;
template class ToolboxMorHeat<heat_3dP1G1_t>;

// FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor_heat_2dP1, ToolboxMorHeat<heat_2dP1G1_t>, toolboxmor_heat_2dP1 )
// FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor_heat_3dP1, ToolboxMorHeat<heat_3dP1G1_t>, toolboxmor_heat_3dP1 )

FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor_heat_2dP1G1, ToolboxMorHeat<heat_2dP1G1_t>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,_2dP1G1) )
FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor_heat_3dP1G1, ToolboxMorHeat<heat_3dP1G1_t>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,_3dP1G1) )


}
