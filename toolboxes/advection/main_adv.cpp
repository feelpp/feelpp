
#include <feel/feelmodels/advection/advection.hpp>

namespace Feel
{

template <uint16_type OrderAdvection>
void
runScalarAdvectionApplication()
{
    using namespace Feel;

    typedef FeelModels::AdvDiffReac< 
        FunctionSpace< 
            Mesh<Simplex<FEELPP_DIM,1>>, 
            bases<Lagrange<OrderAdvection, Scalar,Continuous,PointSetFekete>>
            >
        > model_type;
    
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

template <uint16_type OrderAdvection>
void
runVectorialAdvectionApplication()
{
    using namespace Feel;

    typedef FeelModels::AdvDiffReac< 
        FunctionSpace< 
            Mesh<Simplex<FEELPP_DIM,1>>, 
            bases<Lagrange<OrderAdvection, Vectorial,Continuous,PointSetFekete>>
            >
        > model_type;
    
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
    advectionoptions.add( toolboxes_options("advection") );
    advectionoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P2,P1 ")
        ("adv-type", Feel::po::value<std::string>()->default_value( "scalar" ), " advected field type : scalar, vectorial ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=advectionoptions,
                   _about=about(_name="application_advection",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    std::string advtype = soption(_name="adv-type");
    if( advtype == "scalar" )
    {
        if ( feapprox == "P2" )
            runScalarAdvectionApplication<2>();
        else if ( feapprox == "P1" )
            runScalarAdvectionApplication<1>();

        else CHECK( false ) << "invalid feapprox " << feapprox;
    }
    else if( advtype == "vectorial" )
    {
        if ( feapprox == "P2" )
            runVectorialAdvectionApplication<2>();
        else if ( feapprox == "P1" )
            runVectorialAdvectionApplication<1>();

        else CHECK( false ) << "invalid feapprox " << feapprox;
    }
    else CHECK( false ) << "invalid adv-type " << advtype;

    return 0;
}
