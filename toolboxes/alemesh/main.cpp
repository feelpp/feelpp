#include <feel/feelmodels/modelmesh/meshale.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmodels/modelcore/options.hpp>

template <Feel::uint16_type OrderGeo>
void
runALEMesh()
{
    using namespace Feel;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,OrderGeo>>);

    auto alemesh = FeelModels::meshale( _mesh=mesh );

    for ( std::string const& bctype : std::vector<std::string>({ "moving","fixed","free" }) )
    {
        std::string opt = (boost::format( "markers.%1%" ) %bctype).str();
        if ( Environment::vm().count( opt ) )
            alemesh->addBoundaryFlags( bctype, Environment::vm()[opt].as<std::vector<std::string> >() );
    }

    alemesh->init();
    alemesh->printAndSaveInfo();

    if ( Environment::vm().count( "mesh-adaptation-function" ) )
    {
        auto eScal = expr( soption(_name="mesh-adaptation-function") );
        //auto uScal = alemesh->functionSpace()->compSpace()->element( eScal );
        //alemesh->updateMetricMeshAdaptation( idv(uScal) );
        alemesh->updateMetricMeshAdaptation( eScal );
    }


    auto disp = alemesh->functionSpace()->element();
    if ( Environment::vm().count( "displacement-imposed" ) )
    {
        auto dispExpr = expr<FEELPP_DIM,1>( soption(_name="displacement-imposed") );
        disp.on(_range=elements(mesh),_expr=dispExpr);
    }

    alemesh->update( disp );

    alemesh->exportResults();

}


int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description myoptions( "application alemesh options" );
    myoptions.add( alemesh_options("") );
    myoptions.add_options()
        ("geo-order", Feel::po::value<int>()->default_value( 1 ), "geo-order")
        ("displacement-imposed", Feel::po::value<std::string>(), "displacement-imposed")
        ("markers.moving", po::value<std::vector<std::string> >()->multitoken(), "list of markers on moving boundary" )
        ("markers.fixed", po::value<std::vector<std::string> >()->multitoken(), "list of markers on fixed boundary" )
        ("markers.free", po::value<std::vector<std::string> >()->multitoken(), "list of markers on free boundary" )
        ("mesh-adaptation-function", Feel::po::value<std::string>(), "mesh-adaptation-function")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
                     _about=about(_name="application_alemesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runALEMesh<1>();
    return 0;
}
