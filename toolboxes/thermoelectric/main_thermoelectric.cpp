/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2016-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>
#include <feel/feelmodels/modelcore/scalardiscretization.hpp>

template <int nDim,int OrderT>
int
runApplicationThermoElectric( Feel::FeelModels::ScalarDiscretization const& spec )
{
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<nDim,OrderT>,
                              Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::Electric< Simplex<nDim,OrderT>,
                                  Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_electric_type;
    typedef FeelModels::ThermoElectric< model_heat_type,model_electric_type> model_type;
    std::shared_ptr<model_type> thermoElectric( new model_type("thermo-electric", "thermoelectric", Environment::worldCommPtr(), "",
                                                            FeelModels::ModelBaseRepository(), RuntimeOrder{spec.polynomialOrder},
                                                            RuntimeOrder{spec.polynomialOrder}, RuntimeOrder{spec.geometryOrder}) );
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
            ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P<k> or P<k>G<g>")
            ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=thermoelectricoptions,
                        _about=about(_name="feelpp_toolbox_thermoelectric",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        Feel::FeelModels::printToolboxApplication( "thermo-electric" );

        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");

        auto spec = FeelModels::parseScalarDiscretization( discretization );
        if ( !spec || ( dimension != 2 && dimension != 3 ) )
        {
            if ( Environment::isMasterRank() )
                std::cerr << "Invalid thermoelectric case: expected dimension 2 or 3 and P<k>[G<g>] with positive orders\n";
            return EXIT_FAILURE;
        }
        if ( spec->polynomialOrder == 1 && spec->geometryOrder == 1 )
            return dimension == 2 ? runApplicationThermoElectric<2,1>( *spec ) : runApplicationThermoElectric<3,1>( *spec );
        return dimension == 2 ? runApplicationThermoElectric<2,Dynamic>( *spec ) : runApplicationThermoElectric<3,Dynamic>( *spec );
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
