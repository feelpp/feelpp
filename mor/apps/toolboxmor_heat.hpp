/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_MOR_TOOLBOXMOR_HEAT_HPP
#define FEELPP_MOR_TOOLBOXMOR_HEAT_HPP

#include <feel/feelmor/toolboxmor.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

template<typename ToolboxType, int Options = 0>
class FEELPP_EXPORT ToolboxMorHeat : public ToolboxMor<typename ToolboxType::space_temperature_type,Options>
{
    using super_type = ToolboxMor<typename ToolboxType::space_temperature_type,Options>;
    using self_type = ToolboxMorHeat<ToolboxType,Options>;
  public:
    using toolbox_type = ToolboxType;
    using toolbox_ptrtype = std::shared_ptr<ToolboxType>;

    explicit ToolboxMorHeat( std::string const& name = "ToolboxMor"/*"ToolboxMorHeat"*/, std::string const& prefix = "" );

    toolbox_ptrtype offlineToolbox() const { return M_offlineToolbox; }

    void initModel() override;
    void initOnlineToolbox( std::shared_ptr<DeimMorModelBase<typename super_type::mesh_type>> heatBoxModel ) override;
    //void initOnline();
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

  private :
    toolbox_ptrtype M_offlineToolbox;
    std::shared_ptr<ModelProperties> M_onlineModelProperties;
};

}

#endif
