/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Romain Hild <romain.hild@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2019-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelmodels/electric/electric.hpp>

#include <feel/feelmodels/modelcore/scalardiscretization.hpp>

template <int nDim,int OrderPotential,int OrderGeo>
int
runApplicationElectricStatic()
{
    using namespace Feel;

    using model_type = FeelModels::Electric< Simplex<nDim,OrderGeo>,
                                             Lagrange<OrderPotential, Scalar,Continuous,PointSetFekete> >;
    std::shared_ptr<model_type> electric( new model_type("electric") );
    electric->init();
    electric->printAndSaveInfo();
    electric->solve();
    electric->exportResults();
    return !electric->checkResults();
}

template <int nDim>
int
runApplicationElectricDynamic( int orderPotential, int orderGeometry )
{
    using namespace Feel;

    using model_type = FeelModels::Electric< Simplex<nDim,Dynamic>,
                                             Lagrange<Dynamic, Scalar,Continuous,PointSetFekete> >;
    std::shared_ptr<model_type> electric( new model_type( "electric", "electric",
                                                           Environment::worldCommPtr(), "", FeelModels::ModelBaseRepository(),
                                                           RuntimeOrder{ orderPotential }, RuntimeOrder{ orderGeometry } ) );
    electric->init();
    electric->printAndSaveInfo();
    electric->solve();
    electric->exportResults();
    return !electric->checkResults();
}

int
main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description electricoptions( "electric options" );
        electricoptions.add( toolboxes_options("electric") );
        electricoptions.add_options()
            ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
            ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P<k> or P<k>G<g>")
            ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=electricoptions,
                        _about=about(_name="toolboxes_electric",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");
        auto discretizationSpec = FeelModels::parseScalarDiscretization( discretization );
        if ( !discretizationSpec )
        {
            if ( Environment::isMasterRank() )
                std::cerr << "Invalid case.discretization='" << discretization << "'. Expected P<k> or P<k>G<g> with k,g >= 1\n";
            return EXIT_FAILURE;
        }
        int status = EXIT_FAILURE;
        if ( discretizationSpec->polynomialOrder == 1 && discretizationSpec->geometryOrder == 1 )
        {
            if ( dimension == 2 )
                status = runApplicationElectricStatic<2,1,1>();
            else if ( dimension == 3 )
                status = runApplicationElectricStatic<3,1,1>();
        }
        else
        {
            if ( dimension == 2 )
                status = runApplicationElectricDynamic<2>( discretizationSpec->polynomialOrder, discretizationSpec->geometryOrder );
            else if ( dimension == 3 )
                status = runApplicationElectricDynamic<3>( discretizationSpec->polynomialOrder, discretizationSpec->geometryOrder );
        }

        if ( status == EXIT_FAILURE && Environment::isMasterRank() && ( dimension != 2 && dimension != 3 ) )
            std::cerr << "Invalid case.dimension='" << dimension << "'. Expected 2 or 3\n";
        return status;
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
