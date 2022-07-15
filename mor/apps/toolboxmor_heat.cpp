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
std::shared_ptr<typename ToolboxMorHeat<ToolboxType,Options>::self_type>
ToolboxMorHeat<ToolboxType,Options>::createReducedBasisModel()
{
    auto heatBox = toolbox_type::New(_prefix="heat");
    heatBox->init();
    heatBox->printAndSaveInfo();

    auto model = std::make_shared<self_type>(soption("toolboxmor.name"));
    model->initOffline( heatBox );

    return model;
}

template<typename ToolboxType, int Options>
void
ToolboxMorHeat<ToolboxType,Options>::initOffline( toolbox_ptrtype toolbox )
{
    M_offlineToolbox = toolbox;
    this->setFunctionSpaces( toolbox->spaceTemperature() );
    auto heatBoxModel = DeimMorModelToolbox<toolbox_type>::New( toolbox );

    if ( M_offlineToolbox->hasModelProperties() )
        this->addModelData( "toolbox_json_setup", M_offlineToolbox->modelProperties().jsonData(), "toolbox_model/setup.json" );

    this->initToolbox(heatBoxModel);

    this->initOnline();// maybe give here heatBoxModel
}

template<typename ToolboxType, int Options>
void
ToolboxMorHeat<ToolboxType,Options>::initOnline()
{
    auto heatBoxModel = DeimMorModelToolbox<toolbox_type>::New("heat");

    auto modelProps = std::make_shared<ModelProperties>();// this->repository().expr(),this->worldCommPtr(), this->prefix(), this->clovm() );
    if( this->hasModelData("toolbox_json_setup") )
    {
        auto & mdata = this->additionalModelData("toolbox_json_setup");
        auto const& jsonData = mdata.template fetch_data<nl::json>( this->crbModelDb().dbRepository() );
        modelProps->setup( jsonData );
    }

    heatBoxModel->setToolboxInitFunction(
        [modelProps](   /*auto*/ typename toolbox_type::mesh_ptrtype  mesh ) {
            auto tbDeim = toolbox_type::New( _prefix="heat"/*M_prefix*/);
            tbDeim->setModelProperties( modelProps );
            tbDeim->setMesh(mesh);
            tbDeim->init();
            //tbDeim->printAndSaveInfo();
            return tbDeim;
        });

    this->initOnlineToolbox(heatBoxModel);
}

template<typename ToolboxType, int Options>
void
ToolboxMorHeat<ToolboxType,Options>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    super_type::setupSpecificityModel( ptree,dbDir );
    this->initOnline();
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
