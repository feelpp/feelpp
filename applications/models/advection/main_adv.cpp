
#include <feel/feelmodels/advection/advection.hpp>

namespace Feel
{

template <uint16_type OrderAdvection>
void
runAdvectionApplication()
{
    using namespace Feel;

    typedef FeelModels::Advection< 
        Simplex<FEELPP_DIM,1>,
        Lagrange<OrderAdvection, Scalar,Continuous,PointSetFekete> > model_type;
    
    auto Adv = model_type::New("advection");

    Adv->init();
    Adv->printAndSaveInfo();

    if ( Adv->isStationary() )
    {
        Adv->solve();
        Adv->exportResults();
    }
    else
    {
        for ( ; !Adv->timeStepBase()->isFinished(); Adv->updateTimeStep() )
        {
            if (Adv->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << Adv->time() << "s \n";
                std::cout << "============================================================\n";
            }

            Adv->solve();
            Adv->exportResults();
        }
    }

}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description advectionoptions( "application advection options" );
    advectionoptions.add( feelmodels_options("advection") );
    advectionoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P2,P1 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=advectionoptions,
                   _about=about(_name="application_advection",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P2" )
        runAdvectionApplication<2>();
    else if ( feapprox == "P1" )
        runAdvectionApplication<1>();
    
    else CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
