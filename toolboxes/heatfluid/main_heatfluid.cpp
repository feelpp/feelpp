/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

template <int OrderT,int OrderV, int OrderP>
void
runApplicationHeatFluid()
{
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<FEELPP_DIM,1>,
                                      Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::FluidMechanics< Simplex<FEELPP_DIM,1>,
                                        Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderP, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::HeatFluid< model_heat_type,model_fluid_type> model_type;
    boost::shared_ptr<model_type> heatFluid( new model_type("heat-fluid") );
    heatFluid->init();
    heatFluid->printAndSaveInfo();

    if (heatFluid->isStationary() )
    {
        heatFluid->solve();
        heatFluid->exportResults();
    }
    else
    {
        for ( ; !heatFluid->timeStepBase()->isFinished(); heatFluid->updateTimeStep() )
        {
            if (heatFluid->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << heatFluid->time() << "s \n";
                std::cout << "============================================================\n";
            }

            heatFluid->solve();
            heatFluid->exportResults();
        }
    }
}

int
main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description heatfluidoptions( "application heat-fluid options" );
    heatfluidoptions.add( toolboxes_options("heat-fluid") );
    heatfluidoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1-P2P1" ), "fe-approximation")
     ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatfluidoptions,
                     _about=about(_name="application_heatfluid",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1-P2P1" )
        runApplicationHeatFluid<1,2,1>();
    else
        CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
