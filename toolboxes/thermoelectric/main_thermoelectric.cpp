/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

template <int nDim,int OrderT>
int
runApplicationThermoElectric()
{
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<nDim,1>,
                              Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::Electric< Simplex<nDim,1>,
                                  Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_electric_type;
    typedef FeelModels::ThermoElectric< model_heat_type,model_electric_type> model_type;
    std::shared_ptr<model_type> thermoElectric( new model_type("thermo-electric") );
    thermoElectric->init();
    thermoElectric->printAndSaveInfo();

    if (thermoElectric->isStationary() )
    {
        thermoElectric->solve();
        thermoElectric->exportResults();
    }
    else
    {
        for ( thermoElectric->startTimeStep() ; !thermoElectric->timeStepBase()->isFinished(); thermoElectric->updateTimeStep() )
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
    return !thermoElectric->checkResults();
}

int
main( int argc, char** argv )
{
    using namespace Feel;
    try
    {
        po::options_description thermoelectricoptions( "application thermo-electric options" );
        thermoelectricoptions.add( toolboxes_options("thermo-electric") );
        thermoelectricoptions.add_options()
            ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
            ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
            ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=thermoelectricoptions,
                        _about=about(_name="feelpp_toolbox_thermoelectric",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        Feel::FeelModels::printToolboxApplication( "thermo-electric" );

        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");

        int status = 0;
        hana::for_each( Pc_t<>, [&discretization, &dimension, &status]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype( hana::at_c<0>( d ) )>::value;
                            constexpr int _torder = std::decay_t<decltype( hana::at_c<1>( d ) )>::value;
                            std::string const& _discretization = hana::at_c<2>( d );
                            if ( dimension == _dim && discretization == _discretization )
                                status = runApplicationThermoElectric<_dim,_torder>(); } );
        return status;
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
