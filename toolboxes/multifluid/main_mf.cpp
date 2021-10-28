#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {

// Warning: OrderPNLevelset must equal the FEELPP_MODELS_LEVELSET_PN_ORDER cmake variable or you must instantiate the appropriate lib
template< int Dim, uint16_type OrderVelocity, uint16_type OrderPressure, uint16_type OrderLevelset = 1, uint16_type OrderPNLevelset = LEVELSET_PN_ORDER>
void
runApplicationMultiFluid()
{
    using namespace Feel;

    typedef Simplex< Dim, FEELPP_GEO_ORDER > convex_type;
    typedef FeelModels::FluidMechanics< convex_type,
                                        Lagrange<OrderVelocity, Vectorial, Continuous, PointSetFekete>,
                                        Lagrange<OrderPressure, Scalar, Continuous, PointSetFekete>
                                            > model_fluid_type;
    typedef FeelModels::LevelSet< convex_type, 
                                  Lagrange<OrderLevelset, Scalar, Continuous, PointSetFekete>,
                                  NoPeriodicity,
                                  Lagrange<OrderPNLevelset, Scalar, Continuous, PointSetFekete>
                                  > model_levelset_type;

    typedef FeelModels::MultiFluid< model_fluid_type, model_levelset_type > model_multifluid_type;

    std::shared_ptr<model_multifluid_type> MFModel( new model_multifluid_type( "multifluid" ) );

    MFModel->init();
    MFModel->printAndSaveInfo();
    if( !MFModel->doRestart() )
        MFModel->exportResults(0);

    typedef FunctionSpace< Mesh<convex_type>,
                           bases<Lagrange<OrderVelocity-1, Scalar, Continuous, PointSetFekete>>
                               > space_velocity_surfdiv_type;
    std::shared_ptr<space_velocity_surfdiv_type> spaceSurfDivU;
    std::shared_ptr<Exporter<typename model_multifluid_type::mesh_type>> surfDivUExporter;

    for( MFModel->startTimeStep(); !MFModel->timeStepBase()->isFinished(); MFModel->updateTimeStep() )
    {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << MFModel->time() << "s \n";
            Feel::cout << "============================================================\n";

            MFModel->solve();
            MFModel->exportResults();
    }
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description multifluidoptions( "application multifluid options" );
    multifluidoptions.add( Feel::toolboxes_options("multifluid") );
    multifluidoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 2 ), "dimension")
        ("case.fm-discretization", Feel::po::value<std::string>()->default_value( "P2P1" ), "fluid discretization : P2P1")
        ("case.ls-discretization", Feel::po::value<std::string>()->default_value( "P1" ), "levelset discretization : P1")
        ("export-velocity-surface-div", Feel::po::value<bool>()->default_value( false ), "compute and export the velocity surface divergence" )
        ;

    Environment env( _argc=argc, _argv=argv,
            _desc=multifluidoptions,
            _about=about(_name="application_multifluid",
                _author="Feel++ Consortium",
                _email="feelpp-devel@feelpp.org"));

    int dim = ioption( _name="case.dimension" );
    std::string fmDiscretization = soption(_name="case.fm-discretization");
    std::string lsDiscretization = soption(_name="case.ls-discretization");

    auto dim_t = hana::make_tuple( hana::int_c<2> );

    auto fmDiscretization_t = hana::make_tuple( 
            //hana::make_tuple( "P1P1", hana::make_tuple( hana::int_c<1>, hana::int_c<1>) ),
            hana::make_tuple( "P2P1", hana::make_tuple( hana::int_c<2>, hana::int_c<1> ) ) 
            );
    auto lsDiscretization_t = hana::make_tuple( 
            hana::make_tuple( "P1", hana::make_tuple( hana::int_c<1> ) ) 
            );

    hana::for_each( 
            hana::cartesian_product( hana::make_tuple( dim_t, fmDiscretization_t, lsDiscretization_t ) ), 
            [&fmDiscretization, &lsDiscretization, &dim]( auto const& d )
                {
                    constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                    std::string const& _fmDiscretization = hana::at_c<0>( hana::at_c<1>(d) );
                    std::string const& _lsDiscretization = hana::at_c<0>( hana::at_c<2>(d) );

                    constexpr int _uorder = std::decay_t<decltype(hana::at_c<0>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                    constexpr int _porder = std::decay_t<decltype(hana::at_c<1>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                    constexpr int _lsorder = std::decay_t<decltype(hana::at_c<0>(hana::at_c<1>( hana::at_c<2>(d)) ))>::value;
                    if ( dim == _dim && fmDiscretization == _fmDiscretization && lsDiscretization == _lsDiscretization )
                        runApplicationMultiFluid<_dim, _uorder, _porder, _lsorder>();
                } 
            );

    return 0;
}
