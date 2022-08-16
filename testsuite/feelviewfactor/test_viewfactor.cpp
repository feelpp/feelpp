/**
 * @file test_viewfactor.cpp
 * @author Christophe Prud'homme (christophe.prudhomme@cemosis.fr)
 * @brief various tests for view factors computations
 * @version 0.1
 * @date 2022-07-29
 * 
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright copyright (c) 2022 Université de Strasbourg
 * 
 */

#define BOOST_TEST_MODULE test_viewfactor
#include <feel/feelcore/testsuite.hpp>

#include <fmt/ostream.h>
#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelviewfactor/unobstructedplanarviewfactor.hpp>

/** use Feel namespace */
using namespace Feel;
using Feel::project;

inline
AboutData
makeAbout()
{
    AboutData about( "test_vf" ,
                     "test_vf" ,
                     "0.2",
                     "nD(n=2,3) test vf",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@cemosis.fr", "" );
    return about;

}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "N", Feel::po::value<int>()->default_value( 2 ), "Number of problems to solve successively" )
    ;
    return opts.add( Feel::feel_options(  ) ).add( Feel::feel_options( "square" ) ).add( feel_options( "cube" ) );
}

template<typename MeshType>
void checkViewFactorEnclosure(std::string const& prefix)
{
    auto mesh = loadMesh( _mesh=new MeshType, _filename=Environment::expand(soption(fmt::format("{}.gmsh.filename",prefix) ) ) );
    auto jsons = vsoption( fmt::format("{}.json.filename", prefix ) );
    BOOST_TEST_MESSAGE( fmt::format( "reading json: {}", Environment::expand(jsons[0]) ) );
    nl::json j;
    std::ifstream f(  Environment::expand(jsons[0]) );
    f >> j;
    BOOST_TEST_MESSAGE(fmt::format("j: {}",j.dump(1)));
 
    UnobstructedPlanarViewFactor<MeshType> upvf( mesh, j );
    upvf.compute();
    BOOST_TEST_MESSAGE( fmt::format("{}", upvf.viewFactors() ) );
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( viewfactor )


BOOST_AUTO_TEST_CASE( test_square )
{
    checkViewFactorEnclosure<Mesh<Simplex<2>>>("square");
}
BOOST_AUTO_TEST_CASE( test_cube )
{
    checkViewFactorEnclosure<Mesh<Simplex<3>>>( "cube" );
}

BOOST_AUTO_TEST_SUITE_END()