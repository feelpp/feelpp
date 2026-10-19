/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
   SPDX-FileCopyrightText: 2026 University of Strasbourg
   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   SPDX-License-Identifier: LGPL-2.1-or-later
*/
#define BOOST_TEST_MODULE entity_measure
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feelmesh/geoelement.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <array>
#include <cmath>

FEELPP_ENVIRONMENT_NO_OPTIONS

using namespace Feel;

namespace
{
template <typename Shape, int RealDim>
void checkMeasure( double expected )
{
    using entity_type = GeoND<RealDim, Shape>;
    typename entity_type::reference_convex_type reference;
    entity_type entity;
    std::array<typename entity_type::point_type, entity_type::numPoints> points;
    for ( uint16_type p = 0; p < entity_type::numPoints; ++p )
    {
        for ( int d = 0; d < RealDim; ++d )
            points[p](d) = 0.;
        for ( int d = 0; d < Shape::nDim; ++d )
            points[p](d) = reference.vertex(p)(d);
        // Tilt embedded entities so an axis-aligned special case cannot hide
        // a wrong reference measure or ambient/topological dimension mix-up.
        if constexpr ( RealDim > Shape::nDim )
            points[p](RealDim-1) = reference.vertex(p)(0);
        entity.setPoint( p, points[p] );
    }
    BOOST_CHECK_SMALL( entity.measure()-expected, 1e-12 );
    auto quadrature = entity_type::imMeasure();
    static_assert( entity_type::quad_meas_type::nDim == Shape::nDim );
    BOOST_CHECK_SMALL( quadrature.weightsSum()-reference.measure(), 1e-12 );
}

template <typename Shape, int RealDim>
void checkPhysicalMeasure( std::array<std::array<double,RealDim>,Shape::numPoints> const& coordinates,
                           double expected )
{
    using entity_type = GeoND<RealDim,Shape>;
    entity_type entity;
    std::array<typename entity_type::point_type,Shape::numPoints> points;
    for ( uint16_type p=0; p<Shape::numPoints; ++p )
    {
        for ( int d=0; d<RealDim; ++d )
            points[p](d)=coordinates[p][d];
        entity.setPoint(p,points[p]);
    }
    BOOST_CHECK_SMALL( entity.measure()-expected, 1e-12 );
}

template <typename Entity>
double triangleArea( Entity const& entity )
{
    std::array<double,3> a, b;
    for ( int d=0; d<3; ++d )
    {
        a[d]=entity.point(1)(d)-entity.point(0)(d);
        b[d]=entity.point(2)(d)-entity.point(0)(d);
    }
    double x=a[1]*b[2]-a[2]*b[1];
    double y=a[2]*b[0]-a[0]*b[2];
    double z=a[0]*b[1]-a[1]*b[0];
    return 0.5*std::sqrt(x*x+y*y+z*z);
}
}

BOOST_AUTO_TEST_CASE( simplex_measures )
{
    checkMeasure<Simplex<1,1,1>,1>( 2. );
    checkMeasure<Simplex<1,1,2>,2>( 2.*std::sqrt(2.) );
    checkMeasure<Simplex<1,1,3>,3>( 2.*std::sqrt(2.) );
    checkMeasure<Simplex<2,1,2>,2>( 2. );
    checkMeasure<Simplex<2,1,3>,3>( 2.*std::sqrt(2.) );
    checkMeasure<Simplex<3,1,3>,3>( 4./3. );
}

BOOST_AUTO_TEST_CASE( hypercube_measures )
{
    checkMeasure<Hypercube<1,1,1>,1>( 2. );
    checkMeasure<Hypercube<1,1,2>,2>( 2.*std::sqrt(2.) );
    checkMeasure<Hypercube<2,1,2>,2>( 4. );
    checkMeasure<Hypercube<2,1,3>,3>( 4.*std::sqrt(2.) );
    checkMeasure<Hypercube<3,1,3>,3>( 8. );
}

BOOST_AUTO_TEST_CASE( physical_entity_measures )
{
    // These coordinates are independent of Reference<>. Lengths come from
    // Euclidean distances, areas from determinants/cross products.
    checkPhysicalMeasure<Simplex<1,1,2>,2>( {{{10.,-7.},{13.,-3.}}}, 5. );
    checkPhysicalMeasure<Simplex<1,1,3>,3>( {{{10.,-7.,5.},{12.,-4.,11.}}}, 7. );
    checkPhysicalMeasure<Simplex<2,1,2>,2>( {{{10.,-7.},{12.,-7.},{10.,-4.}}}, 3. );
    checkPhysicalMeasure<Simplex<2,1,3>,3>( {{{10.,-7.,5.},{12.,-7.,6.},{10.,-4.,7.}}}, std::sqrt(61.)/2. );
    checkPhysicalMeasure<Hypercube<2,1,2>,2>( {{{10.,-7.},{12.,-7.},{12.,-4.},{10.,-4.}}}, 6. );
    checkPhysicalMeasure<Hypercube<2,1,3>,3>( {{{10.,-7.,5.},{12.,-7.,6.},{12.,-4.,8.},{10.,-4.,7.}}}, std::sqrt(61.) );
    checkPhysicalMeasure<Simplex<3,1,3>,3>( {{{10.,20.,30.},{12.,20.,30.},{10.,23.,30.},{10.,20.,34.}}}, 4. );
}

BOOST_AUTO_TEST_CASE( extracted_physical_surface_measures )
{
    using mesh_type=Mesh<Simplex<3>>;
    auto mesh=std::make_shared<mesh_type>(Environment::worldCommPtr());
    int rank=mesh->worldComm().globalRank();
    std::array<std::array<double,3>,4> coordinates={{{10.,20.,30.},{12.,20.,30.},{10.,23.,30.},{10.,20.,34.}}};
    for ( int p=0; p<4; ++p )
    {
        mesh_type::point_type point(4*rank+p);
        for ( int d=0; d<3; ++d )
            point(d)=coordinates[p][d]+(d==0 ? 10.*rank : 0.);
        point.setProcessId(rank);
        point.setProcessIdInPartition(rank);
        mesh->addPoint(std::move(point));
    }
    mesh_type::element_type element;
    element.setId(rank);
    element.setProcessId(rank);
    element.setProcessIdInPartition(rank);
    for ( int p=0; p<4; ++p )
        element.setPoint(p,mesh->point(4*rank+p));
    mesh->addElement(std::move(element));
    mesh->components().set(MESH_UPDATE_FACES|MESH_UPDATE_EDGES);
    mesh->updateForUse();
    // updateForUse may renumber element IDs; inspect the owned range rather
    // than assuming the input ID still equals the MPI rank.
    int cellCount=0;
    for ( auto const& wrapped : elements(mesh) )
    {
        BOOST_CHECK_SMALL( unwrap_ref(wrapped).measure()-4., 1e-12 );
        ++cellCount;
    }
    BOOST_CHECK_EQUAL(cellCount,1);
    for ( auto const& wrapped : boundaryfaces(mesh) )
    {
        auto const& face=unwrap_ref(wrapped);
        double expected=triangleArea(face);
        BOOST_CHECK_SMALL( face.measure()-expected, 1e-12 );
        BOOST_CHECK_SMALL( face.element(0).faceMeasure(face.idInElement0())-expected, 1e-12 );
    }
    auto surface=createSubmesh(_range=boundaryfaces(mesh),_view=true);
    double total=0.;
    for ( auto const& wrapped : elements(surface) )
    {
        auto const& face=unwrap_ref(wrapped);
        double expected=triangleArea(face);
        BOOST_CHECK_SMALL( face.measure()-expected, 1e-12 );
        total+=face.measure();
    }
    // Three orthogonal faces have areas 3, 4 and 6; the oblique face sqrt(61).
    BOOST_CHECK_SMALL( total-(13.+std::sqrt(61.)), 1e-11 );
}
