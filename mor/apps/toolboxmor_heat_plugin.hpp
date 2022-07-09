/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_MOR_TOOLBOXMOR_PLUGIN_HPP
#define FEELPP_MOR_TOOLBOXMOR_PLUGIN_HPP

#include <feel/feelmor/toolboxmor.hpp>

namespace Feel
{

template<typename ToolboxType, int Options = 0>
class FEELPP_EXPORT ToolboxMorPlugin : public ToolboxMor<typename ToolboxType::space_temperature_type,Options>
{
    using super_type = ToolboxMor<typename ToolboxType::space_temperature_type,Options>;
    using self_type = ToolboxMorPlugin<ToolboxType,Options>;
  public:
    using toolbox_type = ToolboxType;
    using toolbox_ptrtype = std::shared_ptr<ToolboxType>;

    explicit ToolboxMorPlugin( std::string const& name = "ToolboxMor"/*"ToolboxMorHeat"*/, std::string const& prefix = "" );

    void initOffline( toolbox_ptrtype toolbox );
    void initOnline();

    toolbox_ptrtype offlineToolbox() const { return M_offlineToolbox; }

    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

    static std::shared_ptr<self_type> createReducedBasisModel();
  private :
    toolbox_ptrtype M_offlineToolbox;
};

}

#endif
