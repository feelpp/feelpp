/*
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 * SPDX-FileContributor: GQMP <gwenael.pallares@gmail.com>
 * SPDX-FileContributor: Romain Hild <romain.hild@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2018-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmodels/modelcore/convergencemode.hpp>
#include <feel/feelmodels/modelcore/scalardiscretization.hpp>

template <typename ToolboxType>
int
runToolboxSimulation( std::shared_ptr<ToolboxType> toolbox )
{
    toolbox->init();
    toolbox->printAndSaveInfo();

    if ( toolbox->isStationary() )
    {
        toolbox->solve();
        toolbox->exportResults();
    }
    else
    {
        if ( !toolbox->doRestart() )
            toolbox->exportResults(toolbox->timeInitial());

        for ( toolbox->startTimeStep() ; !toolbox->timeStepBase()->isFinished(); toolbox->updateTimeStep() )
        {
            if (toolbox->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << toolbox->time() << "s \n";
                std::cout << "============================================================\n";
            }

            toolbox->solve();
            toolbox->exportResults();
        }
    }
    return !toolbox->checkResults();
}


template <int Dim, int Order>
int runHeatSimulation( Feel::FeelModels::ScalarDiscretization const& spec, std::string const& mode )
{
    using namespace Feel;
    using model_type = FeelModels::Heat<Simplex<Dim,Order>, Lagrange<Order,Scalar,Continuous,PointSetFekete>>;
    if ( mode == "h-convergence" )
        return runHConvergence<model_type>( "heat", runToolboxSimulation<model_type>,
            [spec]( std::string const& prefix, FeelModels::ModelBaseRepository const& repository, po::variables_map const& vm ) {
                return model_type::New( _prefix=prefix, _repository=repository, _vm=vm,
                                        _runtime_order=RuntimeOrder{spec.polynomialOrder}, _order=RuntimeOrder{spec.geometryOrder} );
            } );
    auto heat = model_type::New( _prefix="heat", _runtime_order=RuntimeOrder{spec.polynomialOrder}, _order=RuntimeOrder{spec.geometryOrder} );
    return runToolboxSimulation( heat );
}

int
main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description heatoptions( "heat options" );
    heatoptions.add( toolboxes_options("heat") );
    heatoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P<k> or P<k>G<g>")
        ("case.mode", Feel::po::value<std::string>()->default_value( "simulation" ), "mode : simulation, h-convergence")
        ("case.mode.h-convergence.hsize", po::value<std::vector<double> >()->multitoken(), "mesh hsize used in h-convergence" )
        ("case.mode.h-convergence.measures", po::value<std::vector<std::string> >()->multitoken(), "measures names used in fit checker" )
        ("case.mode.h-convergence.slopes", po::value<std::vector<double> >()->multitoken(), "reference slope for the fit checker" )
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatoptions,
                     _about=about(_name="toolboxes_heat",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    Feel::FeelModels::printToolboxApplication( "heat" );

    std::string mode = soption(_name="case.mode");
    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto spec = FeelModels::parseScalarDiscretization( discretization );
    if ( !spec || ( dimension != 2 && dimension != 3 ) ||
         ( mode != "simulation" && mode != "h-convergence" ) )
    {
        if ( Environment::isMasterRank() )
            std::cerr << "Invalid heat case: expected dimension 2 or 3, P<k>[G<g>] with positive orders, and simulation or h-convergence mode\n";
        return EXIT_FAILURE;
    }
    if ( spec->polynomialOrder == 1 && spec->geometryOrder == 1 )
        return dimension == 2 ? runHeatSimulation<2,1>( *spec, mode ) : runHeatSimulation<3,1>( *spec, mode );
    return dimension == 2 ? runHeatSimulation<2,Dynamic>( *spec, mode ) : runHeatSimulation<3,Dynamic>( *spec, mode );
}
