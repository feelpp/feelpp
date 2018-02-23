/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

namespace Feel
{

template <uint16_type OrderVelocity,uint16_type OrderPressure, uint16_type OrderGeo = 1>
void
runApplicationFluid()
{
    using namespace Feel;

    typedef FeelModels::FluidMechanics< Simplex<FEELPP_DIM,OrderGeo>,
                                        Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_type;
    auto FM = model_type::New("fluid");

    FM->init();
    FM->printAndSaveInfo();

    if ( FM->isStationary() )
    {
        FM->solve();
        FM->exportResults();
    }
    else
    {
        if ( !FM->doRestart() )
            FM->exportResults(FM->timeInitial());

        for ( ; !FM->timeStepBase()->isFinished(); FM->updateTimeStep() )
        {
            if (FM->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << FM->time() << "s \n";
                std::cout << "============================================================\n";
            }

            FM->solve();
            FM->exportResults();
        }
    }

}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description fluidmecoptions( "application fluid-mechanics options" );
    fluidmecoptions.add( toolboxes_options("fluid") );
    fluidmecoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P2P1" ), "fe-approximation : P2P1,P1P1 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=fluidmecoptions,
                   _about=about(_name="application_fluid",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P2P1" || feapprox == "P2P1G1" )
        runApplicationFluid<2,1>();
#if 0// FEELPP_DIM == 2
    else if ( feapprox == "P1P1" || feapprox == "P1P1G1" )
        runApplicationFluid<1,1>();
#endif
#if 1
    else if ( feapprox == "P2P1G2" )
        runApplicationFluid<2,1,2>();
#endif
    else CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
