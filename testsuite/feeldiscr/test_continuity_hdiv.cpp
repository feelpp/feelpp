/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-20

  Copyright (C) 2014-2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_continuity.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-20
 */
#include <sstream>

#define BOOST_TEST_MODULE continuity testsuite
#include <feel/feelcore/testsuite.hpp>

#include <boost/mpl/list.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/doflayout.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelpoly/brezzidouglasmarini.hpp>

#include <algorithm>
#include <map>
#include <vector>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( continuitysuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
    typedef boost::mpl::list<boost::mpl::pair<boost::mpl::int_<2>,boost::mpl::int_<0> >,
                             // boost::mpl::pair<boost::mpl::int_<2>,boost::mpl::int_<1> >,
                             // boost::mpl::pair<boost::mpl::int_<2>,boost::mpl::int_<2> >,
                             boost::mpl::pair<boost::mpl::int_<3>,boost::mpl::int_<0> >
                             > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( HDivRT0, T, dim_types )
{
    constexpr uint16_type nDim = T::first::value;
    constexpr uint16_type nOrder = T::second::value;

    BOOST_TEST_MESSAGE( "check continuity for HDivRT in  " << nDim << "D P" << nOrder << "\n" );

    typedef Mesh<Simplex<nDim,1> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Dh<nOrder>( mesh );
    auto u = Xh->element();
    auto a1 = form1( _test=Xh );
    a1  = integrate( _range=internalfaces( mesh ), _expr=(leftface(trans(id(u))*N())+rightface(-trans(id(u))*N())) );
    std::string exprstr = (nDim==2)? "{x*y,x+y}:x:y" : "{x*y,x+z,-z*y}:x:y:z";
    u.on(  _range=elements(mesh), _expr=expr<nDim,1>(exprstr) );
    //a1.vector().printMatlab("HDivRT.m");
    BOOST_CHECK_SMALL( a1( u ), 1e-10 );
    //u.printMatlab("uRT0.m");

    BOOST_TEST_MESSAGE( "HDivRT, a1(u)=" << a1(u)  );
    BOOST_TEST_MESSAGE( "check continuity for HDivRT in  " << nDim << "D P" << nOrder << " done\n" );
}

template<typename SpacePtrType>
void
checkInternalFaceDofConsistency( SpacePtrType const& Xh,
                                 typename SpacePtrType::element_type::mesh_ptrtype const& mesh,
                                 bool requireNonIdentityPermutation )
{
    using space_type = typename SpacePtrType::element_type;
    using fe_type = typename space_type::fe_type;
    using size_type = typename space_type::size_type;
    using mesh_type = typename space_type::mesh_type;
    using face_permutation_type = typename mesh_type::element_type::face_permutation_type;

    auto const& fe = *Xh->fe();
    auto const localDofCount = fe.localDofPerComponent();
    auto dof = Xh->dof();

    int checkedFaces = 0;
    bool hasNonIdentityPermutation = false;

    for ( auto fit = mesh->beginFace(), fend = mesh->endFace(); fit != fend; ++fit )
    {
        auto const& face = fit->second;
        if ( face.isOnBoundary() )
            continue;

        const size_type e0 = face.ad_first();
        const size_type e1 = face.ad_second();
        const uint16_type lf0 = face.pos_first();
        const uint16_type lf1 = face.pos_second();
        if ( e0 == invalid_v<size_type> || e1 == invalid_v<size_type> ||
             lf0 == invalid_uint16_type_value || lf1 == invalid_uint16_type_value )
            continue;

        auto const& elt0 = mesh->element( e0 );
        auto const& elt1 = mesh->element( e1 );
        if constexpr ( fe_type::nDim == 3 )
        {
            hasNonIdentityPermutation =
                hasNonIdentityPermutation ||
                ( elt0.facePermutation( lf0 ) != face_permutation_type( face_permutation_type::IDENTITY ) ) ||
                ( elt1.facePermutation( lf1 ) != face_permutation_type( face_permutation_type::IDENTITY ) );
        }

        std::map<size_type, int> g0ToSign;
        std::map<size_type, int> g1ToSign;

        auto const& signs0 = dof->localToGlobalSigns( e0 );
        auto const& signs1 = dof->localToGlobalSigns( e1 );
        auto const& transforms0 = dof->localToGlobalTransforms( e0 );
        auto const& transforms1 = dof->localToGlobalTransforms( e1 );
        BOOST_REQUIRE_GE( transforms0.size(), static_cast<std::size_t>( localDofCount ) );
        BOOST_REQUIRE_GE( transforms1.size(), static_cast<std::size_t>( localDofCount ) );

        for ( uint16_type ldof = 0; ldof < localDofCount; ++ldof )
        {
            auto const attachment = fe.dofAttachment( ldof );
            if ( !attachment.isValid() || attachment.entityDim != 2 )
                continue;

            if ( attachment.entityId == lf0 )
            {
                auto const& transform = transforms0[ldof];
                BOOST_CHECK_EQUAL( dof->dofTransformSignProjection( transform ), signs0( ldof ) );
                if constexpr ( FiniteElementDofTransformProvider<fe_type, typename mesh_type::element_type> )
                {
                    auto const expectedTransform = fe.dofTransform( elt0, ldof );
                    BOOST_CHECK_EQUAL( static_cast<int>( transform.kind ), static_cast<int>( expectedTransform.kind ) );
                    BOOST_CHECK_EQUAL( transform.sign, expectedTransform.sign );
                }
                const size_type gdof = dof->localToGlobal( e0, ldof ).index();
                g0ToSign[gdof] = signs0( ldof );
            }
            if ( attachment.entityId == lf1 )
            {
                auto const& transform = transforms1[ldof];
                BOOST_CHECK_EQUAL( dof->dofTransformSignProjection( transform ), signs1( ldof ) );
                if constexpr ( FiniteElementDofTransformProvider<fe_type, typename mesh_type::element_type> )
                {
                    auto const expectedTransform = fe.dofTransform( elt1, ldof );
                    BOOST_CHECK_EQUAL( static_cast<int>( transform.kind ), static_cast<int>( expectedTransform.kind ) );
                    BOOST_CHECK_EQUAL( transform.sign, expectedTransform.sign );
                }
                const size_type gdof = dof->localToGlobal( e1, ldof ).index();
                g1ToSign[gdof] = signs1( ldof );
            }
        }

        BOOST_REQUIRE_EQUAL( g0ToSign.size(), g1ToSign.size() );
        BOOST_REQUIRE( !g0ToSign.empty() );
        for ( auto const& [gdof, sign0] : g0ToSign )
        {
            auto const it = g1ToSign.find( gdof );
            BOOST_REQUIRE( it != g1ToSign.end() );
            BOOST_CHECK_EQUAL( sign0, -it->second );
        }
        ++checkedFaces;
    }

    BOOST_CHECK_GT( checkedFaces, 0 );
    if ( requireNonIdentityPermutation )
        BOOST_CHECK( hasNonIdentityPermutation );
}

BOOST_AUTO_TEST_CASE( HDivRT0_FacePermutationConsistency3D )
{
    using mesh_type = Mesh<Simplex<3,1>>;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = Dh<0>( mesh );
    checkInternalFaceDofConsistency( Xh, mesh, false );
}

BOOST_AUTO_TEST_CASE( HDivBDM1_FacePermutationConsistency3D )
{
    using mesh_type = Mesh<Simplex<3,1>>;
    using basis_type = bases<BrezziDouglasMarini<0>>;
    using space_type = FunctionSpace<mesh_type, basis_type>;

    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = space_type::New( mesh );
    checkInternalFaceDofConsistency( Xh, mesh, true );
}


BOOST_AUTO_TEST_SUITE_END()
