/*
 * SPDX-FileCopyrightText: 2015-2026 Feel++ Consortium
 * SPDX-FileContributor: Thomas Lantz
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
/** @file test_integrateQuadra.cpp
 * @brief Data-driven quadrature and projection checks sharing one suite-owned mesh.
 */
/// [all]

#define BOOST_TEST_MODULE test_integrateQuadra
#include <boost/test/data/test_case.hpp>
#include <boost/filesystem.hpp>
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
#include <memory>
#include <stdexcept>

using namespace Feel;

/** @brief Describe the quadrature test and retain its original authorship. */
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

/** @brief Own the common mesh for one suite, within the Feel++ environment lifetime. */
class QuadratureMeshFixture
{
public:
    /** Reference quadrilateral mesh used by every expression dataset. */
    using MeshType = Mesh<Hypercube<2>>;
    /** Owning mesh handle; data fixtures retain it only during their test case. */
    using MeshPtr = typename MeshType::mesh_ptrtype;

    /** @brief Prepare this test's expression directory and generate the mesh once. */
    QuadratureMeshFixture()
    {
        // Preserve CTest's distinct serial/MPI repositories, including --rm cleanup.
        if ( Feel::Environment::isMasterRank() )
            fs::create_directories( fs::path( Feel::Environment::exprRepository() ) );
        Feel::Environment::worldComm().barrier();
        M_mesh = createGMSHMesh( _mesh=new MeshType,
                                _desc=domain( _name="polymere", _xmax=1, _ymax=1 ) );
        s_mesh = M_mesh;
    }

    /** @brief Release the suite owner before the global Feel++ environment tears down. */
    ~QuadratureMeshFixture() = default;

    /** @brief Acquire the active suite mesh without giving a global variable ownership. */
    static MeshPtr sharedMesh()
    {
        auto mesh = s_mesh.lock();
        if ( !mesh ) throw std::logic_error( "quadrature mesh suite fixture is not active" );
        return mesh;
    }

private:
    /** Sole mesh owner across datasets; destroyed with this suite fixture. */
    MeshPtr M_mesh;
    /** Non-owning access for data fixtures; safe after the Feel++ environment shuts down. */
    inline static std::weak_ptr<MeshType> s_mesh;
};

/** @brief Give each expression case access to the suite-owned mesh. */
class IntegrateQuadraFixture
{
public:
    /** @brief Retain the mesh for the current case, preserving independent expression state. */
    IntegrateQuadraFixture() : mesh( QuadratureMeshFixture::sharedMesh() ) {}

    /** Mesh shared across cases; each projection still owns its function-space state. */
    QuadratureMeshFixture::MeshPtr mesh;
};

namespace
{
/** Relative tolerance for the default/multiscale quadrature comparison. */
constexpr double kTolerancePct = 5.0;

/** Independent polynomial and transcendental integrands. */
const std::vector<std::string> kRunExprs = {
    "x:x:y",
    "x+y:x:y",
    "cos(x)*sin(y):x:y",
    "y*exp(x):x:y"};

/** Projection datasets covering polynomial and transcendental behavior. */
const std::vector<std::string> kResolExprs = {
    "sin(x):x:y",
    "x+y:x:y",
    "x*y:x:y",
    "cos(x*y):x:y"};
} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( integrQuadraSuite, *boost::unit_test::fixture<QuadratureMeshFixture>() )

/** @test Default and multiscale rules agree on volume and boundary integrals. */
BOOST_DATA_TEST_CASE_F( IntegrateQuadraFixture, integralsMatchDefaultAndMultiscale,
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

/** @test Projection and quadrature produce the same L2 approximation. */
BOOST_DATA_TEST_CASE_F( IntegrateQuadraFixture, projectionMatchesQuadrature,
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
