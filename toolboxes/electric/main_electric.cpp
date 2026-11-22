/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/electric/electric.hpp>

#include <optional>
#include <regex>

namespace
{
struct DiscretizationSpec
{
    int potentialOrder = 1;
    int geometryOrder = 1;
};

std::optional<DiscretizationSpec>
parseDiscretization( std::string const& discretization )
{
    static const std::regex re( R"(^P([0-9]+)(?:G([0-9]+))?$)" );
    std::smatch m;
    if ( !std::regex_match( discretization, m, re ) )
        return std::nullopt;

    DiscretizationSpec spec;
    spec.potentialOrder = std::stoi( m[1].str() );
    spec.geometryOrder = m[2].matched ? std::stoi( m[2].str() ) : 1;
    if ( spec.potentialOrder < 1 || spec.geometryOrder < 1 )
        return std::nullopt;
    return spec;
}
} // namespace

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
runApplicationElectricDynamic( int orderPotential )
{
    using namespace Feel;

    using model_type = FeelModels::Electric< Simplex<nDim,1>,
                                             Lagrange<Dynamic, Scalar,Continuous,PointSetFekete> >;
    std::shared_ptr<model_type> electric( new model_type( "electric", "electric",
                                                           Environment::worldCommPtr(), "", FeelModels::ModelBaseRepository(),
                                                           RuntimeOrder{ orderPotential } ) );
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
        auto discretizationSpec = parseDiscretization( discretization );
        if ( !discretizationSpec )
        {
            if ( Environment::isMasterRank() )
                std::cerr << "Invalid case.discretization='" << discretization << "'. Expected P<k> or P<k>G<g> with k,g >= 1\n";
            return EXIT_FAILURE;
        }
        if ( discretizationSpec->geometryOrder != 1 )
        {
            if ( Environment::isMasterRank() )
                std::cerr << "Unsupported geometry order G" << discretizationSpec->geometryOrder
                          << " for electric toolbox. Use G1.\n";
            return EXIT_FAILURE;
        }

        int status = EXIT_FAILURE;
        if ( discretizationSpec->potentialOrder == 1 && discretizationSpec->geometryOrder == 1 )
        {
            if ( dimension == 2 )
                status = runApplicationElectricStatic<2,1,1>();
            else if ( dimension == 3 )
                status = runApplicationElectricStatic<3,1,1>();
        }
        else
        {
            if ( dimension == 2 )
                status = runApplicationElectricDynamic<2>( discretizationSpec->potentialOrder );
            else if ( dimension == 3 )
                status = runApplicationElectricDynamic<3>( discretizationSpec->potentialOrder );
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
