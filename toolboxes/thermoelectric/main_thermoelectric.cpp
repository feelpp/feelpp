/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

template <int OrderT>
void
runApplicationThermoElectric()
{
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<FEELPP_DIM,1>,
                                      Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::Electric< Simplex<FEELPP_DIM,1>,
                                  Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_electric_type;
    typedef FeelModels::ThermoElectric< model_heat_type,model_electric_type> model_type;
    boost::shared_ptr<model_type> thermoElectric( new model_type("thermo-electric") );
    thermoElectric->init();
    thermoElectric->printAndSaveInfo();

    if (thermoElectric->isStationary() )
    {
        thermoElectric->solve();
        thermoElectric->exportResults();
    }
    else
    {
        for ( ; !thermoElectric->timeStepBase()->isFinished(); thermoElectric->updateTimeStep() )
        {
            if (thermoElectric->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << thermoElectric->time() << "s \n";
                std::cout << "============================================================\n";
            }

            thermoElectric->solve();
            thermoElectric->exportResults();
        }
    }
}

int
main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description thermoelectricoptions( "application thermo-electric options" );
    thermoelectricoptions.add( toolboxes_options("thermo-electric") );
    thermoelectricoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2,P3 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=thermoelectricoptions,
                     _about=about(_name="application_thermoelectric",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1" )
        runApplicationThermoElectric<1>();
#if FEELPP_INSTANTIATION_ORDER_MAX >= 2
    else if ( feapprox == "P2" )
        runApplicationThermoElectric<2>();
#endif
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    else if ( feapprox == "P3" )
        runApplicationThermoElectric<3>();
#endif
    else
        CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
