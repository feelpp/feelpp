
#include <feel/feelmodels/levelset/levelset.hpp>

namespace Feel
{

template <uint16_type OrderLevelset, uint16_type OrderLevelsetPN = LEVELSET_PN_ORDER >
void
runLevelsetApplication()
{
    using namespace Feel;

    typedef Simplex<FEELPP_DIM,1> convex_type;
    typedef Lagrange<OrderLevelset, Scalar, Continuous, PointSetFekete> basis_type;
    //typedef FunctionSpace<Mesh<convex_type>, bases<typename FeelModels::detail::ChangeBasisPolySet<Vectorial, basis_type>::type>[>, Periodicity<NoPeriodicity><] > space_advection_velocity_type;
    //typedef FunctionSpace<Mesh<convex_type>, Feel::detail::bases<Lagrange<OrderLevelset, Vectorial, Continuous, PointSetFekete>>, double, Periodicity<NoPeriodicity>, mortars<NoMortar>> space_advection_velocity_type;
    typedef Lagrange<OrderLevelsetPN, Scalar, Continuous, PointSetFekete> basis_PN_type;

    typedef FeelModels::LevelSet<
        convex_type,
        basis_type,
        NoPeriodicity,
        basis_PN_type
        > model_type;
    
    auto LS = model_type::New("levelset");

    double ls_x0 = 0.5, ls_y0 = 0.5, ls_z0 = 0.5;
    double ls_radius = 0.25;

    //auto phi_init = sqrt( 
            //(vf::Px()-ls_x0)*(vf::Px()-ls_x0) 
            //+ (vf::Py()-ls_y0) * (vf::Py()-ls_y0) 
            ////+ (vf::Pz()-ls_z0) * (vf::Pz()-ls_z0) 
            //) - ls_radius;
    //LS->setInitialValue( phi_init );

    LS->init();
    LS->printAndSaveInfo();
    if( !LS->doRestart() )
        LS->exportResults(0.);

    bool exportDistToBoundary = boption( _name="export-dist-to-boundary" );

    std::shared_ptr<Exporter<typename model_type::mesh_type>> myExporter;

    if( exportDistToBoundary )
    {
        myExporter = exporter( 
                _mesh=LS->mesh(),
                _name="DistExport",
                _path=LS->exporterPath()
                );
    }

    if ( LS->isStationary() )
    {
        LS->solve();
        LS->exportResults();
    }
    else
    {
        int redist_every = ioption( _name="levelset.redist-every" );

        for ( int iter = 1; !LS->timeStepBase()->isFinished(); LS->updateTimeStep(), ++iter )
        {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << LS->time() << "s \n";
            Feel::cout << "============================================================\n";

            Feel::cout << "Iter since redist: " << LS->iterSinceRedistanciation() << std::endl;
            Feel::cout << "Levelset BDF order: " << LS->timeStepBDF()->timeOrder() << std::endl;

            LS->solve();
            if( redist_every > 0 && iter%redist_every == 0 )
            {
                Feel::cout << "Reinitializing... ";
                LS->redistanciate();
                Feel::cout << "done\n";
            }

            LS->exportResults();
            if( exportDistToBoundary )
            {
                auto distToBoundary = LS->distToBoundary();
                myExporter->step(iter)->add("distToBoundary", *distToBoundary );
                myExporter->save();
            }
        }
    }
}

}  // namespace Feel

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description levelsetoptions( "application levelset options" );
    levelsetoptions.add( toolboxes_options( "levelset" ) );
    levelsetoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P2, P1" )
        ("levelset.redist-every", Feel::po::value<int>()->default_value( -1 ), "redistantiate levelset every n iterations" )
        ("export-dist-to-boundary", Feel::po::value<bool>()->default_value( false ), "compute and export the distance to the boundary" )
        ;

    Environment env( 
            _argc=argc, _argv=argv,
            _desc=levelsetoptions,
            _about=about(_name="application_levelset",
                _author="Feel++ Consortium",
                _email="feelpp-devel@feelpp.org")
            );

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P2" )
        runLevelsetApplication<2>();
    else if ( feapprox == "P1" )
        runLevelsetApplication<1>();
    
    else CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
