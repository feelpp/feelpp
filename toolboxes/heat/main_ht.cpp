/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/heat/heat.hpp>

template <int OrderT>
void
runApplicationHeat()
{
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<FEELPP_DIM,1>,
                                      Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_type;
    boost::shared_ptr<model_type> heat( new model_type("heat") );
    heat->init();
    heat->printAndSaveInfo();

    if ( heat->isStationary() )
    {
        heat->solve();
        heat->exportResults();
    }
    else
    {
        for ( ; !heat->timeStepBase()->isFinished(); heat->updateTimeStep() )
        {
            if (heat->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << heat->time() << "s \n";
                std::cout << "============================================================\n";
            }

            heat->solve();
            heat->exportResults();
        }
    }
}

int
main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description heatoptions( "heat options" );
    heatoptions.add( toolboxes_options("heat") );
    heatoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2,P3 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatoptions,
                     _about=about(_name="toolboxes_heat",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1" )
        runApplicationHeat<1>();
#if FEELPP_INSTANTIATION_ORDER_MAX >= 2
    else if ( feapprox == "P2" )
        runApplicationHeat<2>();
#endif
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    else if ( feapprox == "P3" )
        runApplicationHeat<3>();
#endif
    else
        CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
