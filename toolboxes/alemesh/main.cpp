#include <feel/feelmodels/modelmesh/meshale.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmodels/modelcore/options.hpp>

template <Feel::uint16_type OrderGeo>
void
runALEMesh()
{
    using namespace Feel;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);

    auto alemesh = FeelModels::meshale( _mesh=mesh );

    alemesh->addBoundaryFlags("moving","Moving");
    alemesh->addBoundaryFlags("fixed","Fixed");
    alemesh->init();
    alemesh->printAndSaveInfo();

    auto dispExpr = expr<FEELPP_DIM,1>( soption(_name="functions.f") );
    auto disp = alemesh->functionSpace()->element( dispExpr );
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
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
                     _about=about(_name="application_alemesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runALEMesh<1>();
    return 0;
}
