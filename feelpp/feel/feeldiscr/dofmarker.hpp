/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 28 Dec 2019

 Copyright (C) 2019 Feel++ Consortium

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
#ifndef FEELPP_DOFMARKER_HPP
#define FEELPP_DOFMARKER_HPP 1

#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

#include <feel/feelcore/feel.hpp>

namespace Feel {

namespace bimaps = boost::bimaps;

/**
 * 
 */
template<typename IndexT = uint32_type>
struct NoDofMarkerPolicy
{
    using index_t = IndexT;
    typedef boost::bimap<index_t,boost::bimaps::multiset_of<size_type> > dof_marker_type;
    typedef typename dof_marker_type::value_type dof2marker;

    //NoDofMarkerPolicy(  
    typename dof_marker_type::right_range_type
    markerToDof( boost::any const& marker )
        {
            typename dof_marker_type::right_range_type x;
            return x;
        }

    typename dof_marker_type::right_range_type
    markerToDofLessThan( boost::any const& marker )
        {
            typename dof_marker_type::right_range_type x;
            return x;
        }
    typename dof_marker_type::right_range_type
    markerToDofGreaterThan( boost::any const& marker )
        {
            typename dof_marker_type::right_range_type x;
            return x;
        }

    void printDofMarker(std::string const& filename )
        {
        }
    void insertDofMarker( int, int )
        {
            //M_dof_marker.insert( dof2marker( itdof->second+shift+c,  marker.value() ) );
        }
    template<typename T>
    void  renumberDofMarker( std::vector<T> const& previousGlobalIdToNewGlobalId )
        {
#if 0
            dof_marker_type newDofMarker;
            for ( auto it = M_dof_marker.left.begin(), en = M_dof_marker.left.end(); it != en; ++it )
                newDofMarker.insert( dof2marker(previousGlobalIdToNewGlobalId[it->first],it->second) );
            M_dof_marker.clear();
            M_dof_marker.swap( newDofMarker );
#endif
        }
    dof_marker_type M_dof_marker;

};

template<typename IndexT = uint32_type>
struct WithDofMarkerPolicy
{
    using index_t = IndexT;
    typedef boost::bimap<index_t,boost::bimaps::multiset_of<size_type> > dof_marker_type;
    typedef typename dof_marker_type::value_type dof2marker;

    WithDofMarkerPolicy() = default;
    WithDofMarkerPolicy( std::shared_ptr<MeshBase<index_t>> const& m  )
        : M_mesh( m  )
        {}
    typename dof_marker_type::right_range_type
    markerToDof( boost::any const& marker )
        {
            using namespace boost::bimaps;
            int id = M_mesh->markerId( marker );
            return M_dof_marker.right.range( id <= _key, _key<id+1 );
        }

    typename dof_marker_type::right_range_type
    markerToDofLessThan( boost::any const& marker )
        {
            using namespace boost::bimaps;
            int id = M_mesh->markerId( marker );
            return M_dof_marker.right.range( unbounded, _key<id );
        }
    typename dof_marker_type::right_range_type
    markerToDofGreaterThan( boost::any const& marker )
        {
            using namespace boost::bimaps;
            int id = M_mesh->markerId( marker );
            return M_dof_marker.right.range( id<_key, unbounded );
        }

    void printDofMarker(std::string const& filename )
        {
            // std::ofstream ofs( filename.c_str() );
            // BOOST_FOREACH( auto dof, _M_dof_marker )
            // {
            //     //ofs << dof.first << " " << dof.second << "\n";
            // }
            std::ofstream ofs( filename.c_str() );
            for( auto dofleft : M_dof_marker.left )
            {
                ofs << dofleft.first << " " << dofleft.second << "\n";
            }
        }
    void insertDofMarker( int key, int v )
        {
            M_dof_marker.insert( dof2marker( key,  v ) );
        }
    template<typename T>
    void  renumberDofMarker( std::vector<T> const& previousGlobalIdToNewGlobalId )
        {
            dof_marker_type newDofMarker;
            for ( auto it = M_dof_marker.left.begin(), en = M_dof_marker.left.end(); it != en; ++it )
                newDofMarker.insert( dof2marker(previousGlobalIdToNewGlobalId[it->first],it->second) );
            M_dof_marker.clear();
            M_dof_marker.swap( newDofMarker );
        }
    dof_marker_type M_dof_marker;
    std::shared_ptr<MeshBase<index_t>> M_mesh;

};



}

#endif
