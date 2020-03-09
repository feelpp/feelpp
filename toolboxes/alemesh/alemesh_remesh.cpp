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
               alemesh->addBoundaryFlags( bctype,  Environment::vm()[opt].as<std::vector<std::string> >() );
           } 
       }
       alemesh->init();
       alemesh->printAndSaveInfo();
       return alemesh;
    };
    auto alemesh = genALEMesh( mesh );
    auto disp = alemesh->functionSpace()->elementPtr();
    auto exSave = [&ex]( double t, auto alemesh, auto disp ) {
                      ex->step( t )->setMesh( alemesh->movingMesh() );
                      ex->step( t )->add( "disp", idv(*disp) );
                      ex->step( t )->add( "etaQ", etaQ( alemesh->movingMesh() ) );
                      ex->step( t )->add( "nsrQ", nsrQ( alemesh->movingMesh() ) );
                      ex->save();
                  };
    exSave( 0, alemesh, disp );
    
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


        if ( Environment::vm().count( "displacement-imposed" ) )
        {
            auto dispExpr = expr<FEELPP_DIM,1>( soption(_name="displacement-imposed") );
            dispExpr.setParameterValues( { {"t", t },{"T0", T0 } } );
            disp->on(_range=elements(alemesh->referenceMesh()),_expr=dispExpr);
        }

        alemesh->update( *disp );


        if ( boption( "remesh" ) )
        {
            auto const& moving_mesh = alemesh->movingMesh();
            auto Xh = Pch<1>( moving_mesh );
            auto met = Xh->element();
            auto etaqmin = etaQ( moving_mesh ).min();

            if ( etaqmin < doption("etaqtol") )
            {
                auto dispExpr = expr<FEELPP_DIM,1>( soption(_name="displacement-imposed") );
                dispExpr.setParameterValues( { {"t", t }, {"T0", T0 } } );
                disp->on(_range=elements(alemesh->referenceMesh()),_expr=-dispExpr);
                alemesh->update( *disp );
                dispExpr.setParameterValues( { {"t", t-dt }, {"T0", T0 } } );
                disp->on(_range=elements(alemesh->referenceMesh()),_expr=dispExpr);
                alemesh->update( *disp );
                if ( !soption( "remesh.metric" ).empty() )
                    met.on( _range=elements(moving_mesh), _expr=expr(soption("remesh.metric")) );
                else
                {
                    auto [ havg, hmin, hmax ] = hMeasures( moving_mesh );
                    met.on( _range=elements(moving_mesh), _expr=cst(hmax)  );
                }
                auto r =  remesher( moving_mesh );
    
                r.setMetric( met );
                auto out = r.execute();
                out->updateForUse();
                out->setMarkerNames( mesh->markerNames() );
                mesh = out;
                alemesh = genALEMesh( mesh );
                disp = alemesh->functionSpace()->elementPtr();
                T0 = t - dt;
                dispExpr.setParameterValues( { {"t", t }, {"T0", T0 } } );
                disp->on(_range=elements(alemesh->referenceMesh()),_expr=dispExpr);
                alemesh->update( *disp );
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
