#include <feel/feelmodels/modelmesh/meshale.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmesh/remesh.hpp>
#include <feel/feeldiscr/quality.hpp>
#include <feel/feells/disttoentityrange.hpp>

using namespace Feel;

/**
 * metric function equal to 
 * - 1 if the distance \p d is less than \p e
 * - exp(n*d) if the distance \p d is greater or equal than \p e 
 */
template<typename DistFieldT, typename DistToRangeExpr, typename = std::enable_if_t<is_functionspace_element_v<DistFieldT>>>
auto flatThenIncreaseAroundEntityRange( DistFieldT const& d, DistToRangeExpr const& e, int coefexp = 3 )
{
    return ( idv(d) < e ) + ( idv(d) > e-cst(1e-8) )*exp(coefexp*idv(d)/h());
}
template <Feel::uint16_type OrderGeo>
void
runALEMesh()
{

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,OrderGeo>>);

    auto ex = exporter( _mesh=mesh, _name="move", _geo="change" );
    auto genALEMesh = []( auto mesh ) {
       auto alemesh = FeelModels::meshale( _mesh=mesh );
       for ( std::string const& bctype : std::vector<std::string>({ "moving","fixed","free" }) )
       {
           std::string opt = (boost::format( "markers.%1%" ) %bctype).str();
           if ( Environment::vm().count( opt ) )
           {
               for ( auto s : Environment::vm()[opt].as<std::vector<std::string> >() )
                   std::cout << "Registering boundary flag " << s  << " to " << bctype << std::endl;
               alemesh->addMarkersInBoundaryCondition( bctype,  Environment::vm()[opt].as<std::vector<std::string> >() );
           } 
       }
       alemesh->init();
       alemesh->printAndSaveInfo();
       return alemesh;
    };
    auto alemesh = genALEMesh( mesh );
    auto disp = alemesh->functionSpace()->elementPtr();
    auto exSave = [&ex]( double t, auto alemesh, auto disp ) {

                      auto Xh = Pch<1>( alemesh->movingMesh() );
                      auto phi = distToEntityRange( Xh, markedfaces( Xh->mesh(), alemesh->markers("moving")) );

                      ex->step( t )->setMesh( alemesh->movingMesh() );
                      ex->step( t )->add( "disp", idv(*disp) );
                      ex->step( t )->add( "etaQ", etaQ( alemesh->movingMesh() ) );
                      ex->step( t )->add( "nsrQ", nsrQ( alemesh->movingMesh() ) );
                      ex->step( t )->add( "dist", phi );
                      auto nlayers = ioption(_name="remesh.metric.layers");
                      ex->step( t )->add( "met", flatThenIncreaseAroundEntityRange(phi,nlayers*h(),nlayers)*expr(soption("remesh.metric")) ); 
                      ex->save();
                  };
    exSave( 0, alemesh, disp );

    auto dispExpr = expr<FEELPP_DIM,1>( soption(_name="displacement-imposed") );
    auto updateDisp = []( auto alemesh, auto disp, auto dexpr, double t, double T0, double dt )
                          {
                              dexpr.setParameterValues( { {"t", t }, {"T0", T0 }, {"dt",dt } } );
                              alemesh->updateDisplacementImposed( dexpr, elements(alemesh->movingMesh()) );
                              alemesh->updateMovingMesh();
                          };
    double dt = doption("dt");
    double T = doption("Tfinal");
    for( double t = dt, T0 = 0; t < T; t += dt )
    {
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "Time: " << t << std::endl;
        
        if ( Environment::vm().count( "mesh-adaptation-function" ) )
        {
            auto eScal = expr( soption(_name="mesh-adaptation-function") );
            //auto uScal = alemesh->functionSpace()->compSpace()->element( eScal );
            //alemesh->updateMetricMeshAdaptation( idv(uScal) );
            alemesh->updateMetricMeshAdaptation( eScal );
        }

        updateDisp( alemesh, disp, dispExpr, t, T0, dt );

        if ( boption( "remesh" ) )
        {
            auto const& moving_mesh = alemesh->movingMesh();
            auto Xh = Pch<1>( moving_mesh );
            auto met = Xh->element();
            auto etaqmin = etaQ( moving_mesh ).min();
            
            if ( etaqmin < doption("etaqtol") )
            {
                // revert to t-dt
                updateDisp( alemesh, disp, -dispExpr, t, T0, dt );
                updateDisp( alemesh, disp, dispExpr, t-dt, T0, dt );

                // apply remesh

                // define distance function from moving boundary to adapt the metric
                auto phi = distToEntityRange( Xh, markedfaces( Xh->mesh(), alemesh->markers("moving")) );
                auto [ havg, hmin, hmax ] = hMeasures( moving_mesh, markedfaces( Xh->mesh(), alemesh->markers("moving"))  );

                if ( !soption( "remesh.metric" ).empty() )
                {
                    auto nlayers = ioption(_name="remesh.metric.layers");
                    met.on( _range=elements(moving_mesh), _expr=flatThenIncreaseAroundEntityRange(phi,nlayers*h(),nlayers)*havg);
                }
                else
                {
                    met.on( _range=elements(moving_mesh), _expr=cst(havg)  );
                }
                    auto r =  remesher( moving_mesh, std::vector<int>{}, alemesh->markers("moving") );
    
                r.setMetric( met );
                auto out = r.execute();
                out->updateForUse();
                out->setMarkerNames( mesh->markerNames() );
                mesh = out;
                alemesh = genALEMesh( mesh );
                disp = alemesh->functionSpace()->elementPtr();
                // revert to t-dt
                T0=t-dt;
                updateDisp( alemesh, disp, dispExpr, t, T0, dt );
            }
        }
        exSave( t, alemesh, disp );
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
        ( "remesh.metric.layers", po::value<int>()->default_value( 3 ), "remesh metric layers" )
        ( "remesh.metric.exp", po::value<std::string>()->default_value( "1" ), "remesh metric expression for exponent" )
        ( "remesh.metric", po::value<std::string>()->default_value( "" ), "remesh metric expression" )
        ( "dt", po::value<double>()->default_value( 0.1 ), "dt" )
        ( "Tfinal", po::value<double>()->default_value( 1 ), "T" )
        ( "etaqtol", po::value<double>()->default_value( 0.1 ), "etaqtol" )
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
                     _about=about(_name="application_alemesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runALEMesh<1>();
    return 0;
}
