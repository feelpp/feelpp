#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {

// Warning: OrderPNLevelset must equal the FEELPP_MODELS_LEVELSET_PN_ORDER cmake variable or you must instantiate the appropriate lib
template< uint16_type OrderVelocity, uint16_type OrderPressure, uint16_type OrderLevelset = 1, uint16_type OrderPNLevelset = LEVELSET_PN_ORDER>
void
runApplicationMultiFluid()
{
    using namespace Feel;

    typedef Simplex<FEELPP_DIM,FEELPP_GEO_ORDER> convex_type;
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
        ("fluid-fe-approximation", Feel::po::value<std::string>()->default_value( "P2P1" ), "fluid-fe-approximation : P2P1,P1P1 ")
        ("levelset-fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "levelset-fe-approximation : P1,P2 ")
        ("export-velocity-surface-div", Feel::po::value<bool>()->default_value( false ), "compute and export the velocity surface divergence" )
        ;

    Environment env( _argc=argc, _argv=argv,
            _desc=multifluidoptions,
            _about=about(_name="application_multifluid",
                _author="Feel++ Consortium",
                _email="feelpp-devel@feelpp.org"));


    std::string fluid_feapprox = soption(_name="fluid-fe-approximation");
    std::string levelset_feapprox = soption(_name="levelset-fe-approximation");
    if ( fluid_feapprox == "P2P1" )
    {
        if( levelset_feapprox == "P1" )
            runApplicationMultiFluid<2,1,1>();
        //else if( levelset_feapprox == "P2" )
            //runApplicationMultiFluid<2,1,2>();
        //else if( levelset_feapprox == "P3" )
            //runApplicationMultiFluid<2,1,3>();
        else CHECK( false ) << "invalid levelset-feapprox " << levelset_feapprox;
    }
    else CHECK( false ) << "invalid fluid-feapprox " << fluid_feapprox;

    return 0;
}
