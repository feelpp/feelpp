/*
  This file is part of the Feel library

  Copyright (C) 2007,2008 EPFL
  Copyright (C) 2010 University of Coimbra

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/**
   \file ale.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2008-04-17
*/
#ifndef __ALE_CPP
#define __ALE_CPP

#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <feel/feeldiscr/ale.hpp>

namespace Feel
{

template < class Convex >
ALE<Convex>::ALE( interval_type const& intX,
                  mesh_ptrtype& mesh,
                  po::variables_map const& vm )
    :
    b( backend_type::build( vm ) ),
    intervalX( intX ),
    reference_mesh( mesh ),
    p1_fspace( p1_functionspace_type::New( reference_mesh ) ),
    pN_fspace( pN_functionspace_type::New( reference_mesh ) ),
    new_mesh( new new_mesh_type ),
    harmonic( b->newMatrix( p1_fspace, p1_fspace ) ),
    pN_ale( pN_fspace, "pN_ale" ),
    pN_displacement( pN_fspace, "pN_displacement" ),
    pN_identity( pN_fspace, "pN_identity" ),
    ho_mesh( new ho_mesh_type( reference_mesh ) )
{
    using namespace Feel::vf;

    // Define harmonic extension operator
    p1_element_type u( p1_fspace, "u" );
    p1_element_type v( p1_fspace, "v" );

    M_timer.restart();
    size_type pattern = Pattern::COUPLED;
    form2( p1_fspace, p1_fspace, harmonic, _init=true, _pattern = pattern ) =
        integrate( elements(reference_mesh),
                   trace( trans(gradt(u))*grad(v) ) );

    form2( p1_fspace, p1_fspace, harmonic ) +=
        integrate( boundaryfaces(reference_mesh),
                   - trans((gradt(u)*N()))*id(v)
                   );

    harmonic->close();

    Log() << "[ALE] Time to generate harmonic extension operator: " << M_timer.elapsed() << "\n";

    pN_ale.container() = ublas::scalar_vector<double>( pN_ale.size(), 0.0 );

    //generate high order mesh in reference domain
    M_timer.restart();
    new_mesh = ho_mesh->getMesh();
    //Log() << "[ALE] Time to generate high order reference mesh: " << M_timer.elapsed() << "\n";

    pN_identity = vf::project( pN_fspace, elements( reference_mesh ), P() );
}




template < class Convex >
ALE<Convex>::ALE( ALE const& tc )
    :
    b( tc.b ),
    intervalX( tc.intervalX ),
    reference_mesh( tc.reference_mesh ),
    p1_fspace( tc.p1_fspace ),
    pN_fspace( tc.pN_fspace ),
    new_mesh( tc.new_mesh ),
    harmonic( tc.harmonic ),
    pN_ale( tc.pN_ale ),
    pN_displacement( tc.pN_displacement ),
    pN_identity( tc.pN_identity ),
    ho_mesh( tc.ho_mesh )
{
}


template < class Convex >
ALE<Convex>::~ALE()
{}




template < class Convex >
typename ALE<Convex>::pN_element_type&
ALE<Convex>::getMap()
{
    return pN_ale;
}

template < class Convex >
typename ALE<Convex>::pN_element_type&
ALE<Convex>::getDisplacement()
{
    return pN_displacement;
}

template < class Convex >
typename ALE<Convex>::pN_element_type&
ALE<Convex>::getIdentity()
{
    return pN_identity;
}

template < class Convex >
typename ALE<Convex>::pN_functionspace_ptrtype
ALE<Convex>::functionSpace()
{
    return pN_fspace;
}


template < class Convex >
typename ALE<Convex>::mesh_ptrtype
ALE<Convex>::getReferenceMesh()
{
    return reference_mesh;
}


template < class Convex >
void
ALE<Convex>::displacement2Map( pN_element_type const& disp, pN_element_type& u )
{
    u = getIdentity();

    u += disp;
    u.updateGlobalValues();
}


template < class Convex >
void
ALE<Convex>::map2Displacement( pN_element_type const& u,
                               pN_element_type& disp )
{
    disp = u;
    disp -= getIdentity();
    disp.updateGlobalValues();
}


template < class Convex >
void
ALE<Convex>::generateP1Map( p1_element_type& p )
{
    using namespace Feel::vf;
    ExporterQuick<mesh_type> exp( "harmonic", "ensight" );
    p1_element_type v( p1_fspace, "v");
    exp.save( 0, p );
    vector_ptrtype rhs( b->newVector( p1_fspace ) );
    rhs->zero();
    rhs->close();

    form2( p1_fspace, p1_fspace, harmonic ) += on( boundaryfaces(reference_mesh), v, rhs, idv(p),
                                                   ON_ELIMINATION );

    vector_ptrtype U( b->newVector( p1_fspace ) );
    b->solve(harmonic, harmonic, U, rhs);
    p = *U;

    exp.save( 1, p );
}

#if defined( FEEL_INSTANTIATION_MODE )
template class ALE< Simplex<2,1> >;
#if BOOST_PP_GREATER_EQUAL( FEEL_MESH_MAX_ORDER, 2 )
template class ALE< Simplex<2,2> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEEL_MESH_MAX_ORDER, 3 )
template class ALE< Simplex<2,3> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEEL_MESH_MAX_ORDER, 4 )
template class ALE< Simplex<2,4> >;
#endif // FEEL_MESH_MAX_ORDER
#if BOOST_PP_GREATER_EQUAL( FEEL_MESH_MAX_ORDER, 5 )
template class ALE< Simplex<2,5> >;
#endif // FEEL_MESH_MAX_ORDER
#endif // FEEL_INSTANTIATION_MODE
}

#endif // __ALE_HPP
