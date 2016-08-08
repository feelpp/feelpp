
#include <feel/feelmodels/levelset/levelset.hpp>

namespace Feel
{

template <uint16_type OrderLevelset>
void
runLevelsetApplication()
{
    using namespace Feel;

    typedef FeelModels::LevelSet< 
        Simplex<FEELPP_DIM,1>,
        OrderLevelset
        > model_type;
    
    auto LS = model_type::New("levelset");
    LS->build();

    double ls_x0 = 0.5, ls_y0 = 0.5, ls_z0 = 0.5;
    double ls_radius = 0.25;

    auto phi_init = sqrt( 
            (vf::Px()-ls_x0)*(vf::Px()-ls_x0) 
            + (vf::Py()-ls_y0) * (vf::Py()-ls_y0) 
            //+ (vf::Pz()-ls_z0) * (vf::Pz()-ls_z0) 
            ) - ls_radius;
    LS->setInitialValue( phi_init );

    LS->init();
    LS->printAndSaveInfo();
    //LS->exportResults();

    if ( LS->isStationary() )
    {
        LS->solve();
        LS->exportResults();
    }
    else
    {
        int reinit_every = ioption( _name="levelset.reinit-every" );

        for ( int iter = 0; !LS->timeStepBase()->isFinished(); LS->updateTimeStep(), ++iter )
        {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << LS->time() << "s \n";
            Feel::cout << "============================================================\n";

            Feel::cout << "Iter since reinit: " << LS->iterSinceReinit() << std::endl;
            Feel::cout << "Levelset BDF order: " << LS->timeStepBDF()->timeOrder() << std::endl;

            LS->solve();
            if( reinit_every > 0 && iter%reinit_every == 0 )
            {
                Feel::cout << "Reinitializing... ";
                LS->reinitialize();
                Feel::cout << "done\n";
            }

            LS->exportResults();
        }
    }
}

}  // namespace Feel

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description levelsetoptions( "application levelset options" );
    levelsetoptions.add( feelmodels_options( "levelset" ) );
    levelsetoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P2, P1" )
        ("levelset.reinit-every", Feel::po::value<int>()->default_value( -1 ), "reinitialize levelset every n iterations" )
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
