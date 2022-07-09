#include "toolboxmor_heat_plugin.hpp"

//#include <feel/feelmor/toolboxmor.hpp>
#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{

template<typename ToolboxType, int Options>
ToolboxMorPlugin<ToolboxType,Options>::ToolboxMorPlugin( std::string const& name, std::string const& prefix )
    :
    super_type( name, prefix )
{}

template<typename ToolboxType, int Options>
std::shared_ptr<typename ToolboxMorPlugin<ToolboxType,Options>::self_type>
ToolboxMorPlugin<ToolboxType,Options>::createReducedBasisModel()
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
ToolboxMorPlugin<ToolboxType,Options>::initOffline( toolbox_ptrtype toolbox )
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
ToolboxMorPlugin<ToolboxType,Options>::initOnline()
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
            auto tbDeim = std::make_shared<toolbox_type>( "heat"/*M_prefix*/);
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
ToolboxMorPlugin<ToolboxType,Options>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
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

FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor_heat_2dP1, ToolboxMorPlugin<heat_2dP1G1_t>, toolboxmor_heat_2dP1 )
FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor_heat_3dP1, ToolboxMorPlugin<heat_3dP1G1_t>, toolboxmor_heat_3dP1 )

template class ToolboxMorPlugin<heat_2dP1G1_t>;

}
