#include <feel/feelmodels/modelmesh/meshale.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmesh/remesh.hpp>
#include <feel/feeldiscr/quality.hpp>

template <Feel::uint16_type OrderGeo>
void
runALEMesh()
{
    using namespace Feel;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,OrderGeo>>);

    auto alemesh = FeelModels::meshale( _mesh=mesh );
    //alemesh->setComputationalDomain( "fluid", markedelements(mesh,"Omega") );
    //alemesh->setComputationalDomain( "fluid", elements(mesh) );
    alemesh->init();

    for ( std::string const& bctype : std::vector<std::string>({ "moving","fixed","free" }) )
    {
        std::string opt = (boost::format( "markers.%1%" ) %bctype).str();
        if ( Environment::vm().count( opt ) )
            alemesh->addMarkersInBoundaryCondition( bctype, Environment::vm()[opt].as<std::vector<std::string> >() );
    }

    alemesh->printAndSaveInfo();

    if ( Environment::vm().count( "mesh-adaptation-function" ) )
    {
        auto eScal = expr( soption(_name="mesh-adaptation-function") );
        //auto uScal = alemesh->functionSpace()->compSpace()->element( eScal );
        //alemesh->updateMetricMeshAdaptation( idv(uScal) );
        alemesh->updateMetricMeshAdaptation( eScal );
    }


    if ( Environment::vm().count( "displacement-imposed" ) )
    {
        auto dispExpr = expr<FEELPP_DIM,1>( soption(_name="displacement-imposed") );
        alemesh->updateDisplacementImposed( dispExpr,markedfaces(mesh,alemesh->markers("moving")) );
    }

    alemesh->updateMovingMesh();

    alemesh->exportResults();

    if ( boption( "remesh" ) )
    {
        auto const& mesh = alemesh->movingMesh();
        auto Xh = Pch<1>( mesh );
        auto met = Xh->element();
        if ( !soption( "remesh.metric" ).empty() )
            met.on( _range=elements(mesh), _expr=expr(soption("remesh.metric")) );
        else
        {
            auto [ havg, hmin, hmax ] = hMeasures( mesh );
            met.on( _range=elements(mesh), _expr=cst(hmax)  );
        }


        
        auto r =  remesher( mesh );
    
        r.setMetric( met );
        auto out = r.execute();
        out->updateForUse();
        auto ein = exporter( _mesh=mesh, _name="meshin" );
        ein->add( "etaQ", etaQ(mesh) );
        ein->add( "nsrQ", nsrQ(mesh) );
        ein->save();
        auto eout = exporter( _mesh=out, _name="remeshed" );
        eout->add( "etaQ", etaQ( out ) );
        eout->add( "nsrQ", nsrQ( out ) );
        eout->save();
    }
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
        ( "remesh", po::value<bool>()->default_value( 0 ), "remesh " )
        ( "remesh.metric", po::value<std::string>()->default_value( "" ), "remesh metric expression" )
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
                     _about=about(_name="application_alemesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runALEMesh<1>();
    return 0;
}
