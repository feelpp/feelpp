/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Thomas Lantz
       Date: 2015-04-27

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/// [all]

#define BOOST_TEST_MODULE test_integrateQuadra
#include <boost/test/data/test_case.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelpoly/multiscalequadrature.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_integrateQuadra" ,
                     "test_integrateQuadra" ,
                     "0.2",
                     "test integrate Quadra",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015-2026 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "christophe.prudhomme@feelpp.org", "", "" );
    return about;
}

namespace bdata = boost::unit_test::data;

class IntegrateQuadraFixture
{
public:
    using mesh_type = Mesh<Hypercube<2>>;
    using mesh_ptrtype = typename mesh_type::mesh_ptrtype;

    IntegrateQuadraFixture()
        : mesh( createGMSHMesh( _mesh=new mesh_type,
                                 _desc=domain( _name="polymere", _xmax=1, _ymax=1 ) ) )
    {
        // Ensure repository exists for GiNaC-generated sources
        Feel::Environment::changeRepository( _directory=boost::format( "test_integrateQuadra" ), _subdir=false );
        const auto repoDir = Feel::Environment::repository().directory();
        fs::create_directories( repoDir / "exprs" );
        fs::create_directories( fs::path( Feel::Environment::exprRepository() ) );
    }

    mesh_ptrtype mesh;
};

namespace
{
constexpr double kTolerancePct = 5.0;

const std::vector<std::string> kRunExprs = {
    "x:x:y",
    "x+y:x:y",
    "cos(x)*sin(y):x:y",
    "y*exp(x):x:y"};

const std::vector<std::string> kResolExprs = {
    "sin(x):x:y",
    "x+y:x:y",
    "x*y:x:y",
    "cos(x*y):x:y"};
} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( integrQuadra_suite )

BOOST_DATA_TEST_CASE_F( IntegrateQuadraFixture, integrals_match_default_and_multiscale,
                        bdata::make( kRunExprs ), exprString )
{
    auto g = expr( exprString );

    auto int_volume_msq = integrate( _range=elements( mesh ),
                                     _expr=g,
                                     _quad=_Q<1,MultiScaleQuadrature>() ).evaluate();
    auto int_volume_std = integrate( _range=elements( mesh ),
                                     _expr=g ).evaluate();

    auto int_boundary_msq = integrate( _range=boundaryfaces( mesh ),
                                       _expr=g,
                                       _quad=_Q<1,MultiScaleQuadrature>() ).evaluate();
    auto int_boundary_std = integrate( _range=boundaryfaces( mesh ),
                                       _expr=g ).evaluate();

    BOOST_TEST_CONTEXT( "expr=" << exprString )
    {
        BOOST_CHECK_CLOSE( int_volume_msq( 0, 0 ), int_volume_std( 0, 0 ), kTolerancePct );
        BOOST_CHECK_CLOSE( int_boundary_msq( 0, 0 ), int_boundary_std( 0, 0 ), kTolerancePct );
    }
}

BOOST_DATA_TEST_CASE_F( IntegrateQuadraFixture, projection_matches_quadrature,
                        bdata::make( kResolExprs ), exprString )
{
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto g = expr( exprString );
    auto gProj = vf::project( _space=Vh, _range=elements( mesh ), _expr=g );

    auto a = form2( _trial=Vh, _test=Vh );
    // Use a higher-order multi-scale quadrature to reduce projection error
    a = integrate( _range=elements( mesh ),
                   _expr=idt( u )*id( v ),
                   _quad=_Q<3,MultiScaleQuadrature>() );

    auto l = form1( _test=Vh );
    l = integrate( _range=elements( mesh ),
                   _expr=g*id( v ),
                   _quad=_Q<3,MultiScaleQuadrature>() );

    a.solve( _rhs=l, _solution=u );

    const auto diff = normL2( _range=elements( mesh ), _expr=idv( u )-idv( gProj ) );
    BOOST_TEST_CONTEXT( "expr=" << exprString )
    {
        BOOST_CHECK_SMALL( diff, 1e-3 );
    }
}

BOOST_AUTO_TEST_SUITE_END()
