#include <feel/feelmodels/multifluid/multifluid.hpp>

namespace Feel {

template< uint16_type OrderVelocity, uint16_type OrderPressure, uint16_type OrderLevelset = 1>
void
runApplicationMultiFluid()
{
    using namespace Feel;

    typedef FeelModels::FluidMechanics< Simplex<FEELPP_DIM,FEELPP_GEO_ORDER>,
                                        Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::LevelSet< Simplex<FEELPP_DIM, FEELPP_GEO_ORDER>,
                                  OrderLevelset > model_levelset_type;

    typedef FeelModels::MultiFluid< model_fluid_type, model_levelset_type > model_multifluid_type;

    auto MFModel = model_multifluid_type::New( "multifluid" );
    MFModel->build();

    MFModel->init();

    for( ; !MFModel->timeStepBase()->isFinished(); MFModel->updateTimeStep() )
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
    multifluidoptions.add( Feel::feelmodels_options("multifluid") );
    multifluidoptions.add_options()
        ("fluid-fe-approximation", Feel::po::value<std::string>()->default_value( "P2P1" ), "fluid-fe-approximation : P2P1,P1P1 ")
        ("levelset-fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "levelset-fe-approximation : P1,P2 ")
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
        else if( levelset_feapprox == "P2" )
            runApplicationMultiFluid<2,1,2>();
        else CHECK( false ) << "invalid levelset-feapprox " << levelset_feapprox;
    }
    /*#if FEELPP_DIM == 2
      else if ( feapprox == "P1P1" )
      runApplicationFSI<1,1>();
#endif*/
    else CHECK( false ) << "invalid fluid-feapprox " << fluid_feapprox;

    return 0;
}
