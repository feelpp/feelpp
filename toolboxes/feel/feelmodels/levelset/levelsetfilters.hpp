//! -*- mode: c++; coding: utf-9; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file levelsetfilters.hpp
//! @author Thibaut Metivet <thibaut.metivet@gmail.com>
//! @date 25 Jul 2019
//! @copyright 2019 Feel++ Consortium
//!

#ifndef _LEVELSET_FILTERS_HPP
#define _LEVELSET_FILTERS_HPP 1

namespace Feel {


template< typename ElementType >
auto
levelsetInterfaceElements( ElementType const& phi )
{
    typedef ElementType element_type;
    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename element_type::mesh_type mesh_type;

    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef Range<mesh_type,MESH_ELEMENTS> range_elements_type;

    auto const& mesh = phi.mesh();
    const rank_type pid = mesh->worldCommElements().localRank();
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    // Find interface elements
    auto it = mesh->beginOrderedElement();
    auto en = mesh->endOrderedElement();

    elements_reference_wrapper_ptrtype interfaceElts( new elements_reference_wrapper_type );

    for ( ; it!=en ; it++ )
    {
        auto const& elt = boost::unwrap_ref( *it );
        if ( elt.processId() != pid )
            continue;
        bool mark_elt = false;
        bool hasPositivePhi = false;
        bool hasNegativePhi = false;
        for (int j=0; j<ndofv; j++)
        {
            if ( phi.localToGlobal(elt.id(), j, 0) < 0. )
            {
                // phi < 0
                if( hasPositivePhi )
                {
                    mark_elt = true;
                    break; //don't need to do the others dof
                }
                hasNegativePhi = true;
            }
            if ( phi.localToGlobal(elt.id(), j, 0) > 0. )
            {
                // phi > 0
                if( hasNegativePhi )
                {
                    mark_elt = true;
                    break; //don't need to do the others dof
                }
                hasPositivePhi = true;
            }
        }
        if( mark_elt )
            interfaceElts->push_back( boost::cref(elt) );
    }

    return boost::make_tuple( 
            mpl::size_t<MESH_ELEMENTS>(),
            interfaceElts->begin(),
            interfaceElts->end(),
            interfaceElts
            );
}

template< typename ElementType >
auto
levelsetDiracElementsGeneric( ElementType const& phi, double eps )
{
    typedef ElementType element_type;
    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename element_type::mesh_type mesh_type;

    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef Range<mesh_type,MESH_ELEMENTS> range_elements_type;

    auto const& mesh = phi.mesh();
    const rank_type pid = mesh->worldCommElements().localRank();
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    // Find interface elements
    auto it = mesh->beginOrderedElement();
    auto en = mesh->endOrderedElement();

    elements_reference_wrapper_ptrtype diracElts( new elements_reference_wrapper_type );

    for ( ; it!=en ; it++ )
    {
        auto const& elt = boost::unwrap_ref( *it );
        if ( elt.processId() != pid )
            continue;
        bool mark_elt = false;
        for ( int j = 0; j < ndofv; j++ )
        {
            if ( std::abs( phi.localToGlobal(elt.id(), j, 0) ) <= eps )
            {
                mark_elt = true;
                break;
            }
        }
        if( mark_elt )
            diracElts->push_back( boost::cref(elt) );
    }

    return boost::make_tuple( 
            mpl::size_t<MESH_ELEMENTS>(),
            diracElts->begin(),
            diracElts->end(),
            diracElts
            );
}

}

#endif // _LEVELSET_FILTERS_HPP
