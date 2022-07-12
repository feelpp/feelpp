
#include <feel/feelmor/toolboxmor.hpp>
#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{
using base_type = bases<Lagrange<1, Scalar, Continuous, PointSetFekete> >;
using space_type = FunctionSpace<Mesh<Simplex<2> >, base_type>;
// FEELPP_CRB_PLUGIN_TEMPLATE( ToolboxMor, ToolboxMor<FunctionSpace<Mesh<Simplex<2> >, base_type>>, toolboxmor )

template class ToolboxMor<space_type>;


po::options_description
makeToolboxMorOptions()
{
    po::options_description options( "ToolboxMor options" );
    options.add_options()
        ( "toolboxmor.filename", po::value<std::string>()->default_value("toolbox.json"), "json file defining the parameters" )
        ( "toolboxmor.trainset-deim-size", po::value<int>()->default_value(40), "size of the deim trainset" )
        ( "toolboxmor.trainset-mdeim-size", po::value<int>()->default_value(40), "size of the mdeim trainset" )
        ( "toolboxmor.test-deim", po::value<bool>()->default_value(false), "test deim decomposition" )
        ;
    options.add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        // .add(deimOptions("vec"))
        // .add(deimOptions("mat"))
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"))
        .add(bdf_options("OpusHeatTB"));
    return options;
}
AboutData
makeToolboxMorAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "toolboxmor",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2009-2014 Feel++ Consortium");
    return about;
}

}
