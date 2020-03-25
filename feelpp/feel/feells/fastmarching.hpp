//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Author(s) : 
//!     Thibaut Metivet <thibaut.metivet@inria.fr>
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
//! @file fastmarching.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _FASTMARCHING_HPP
#define _FASTMARCHING_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include "fastmarchingdofstatus.hpp"
#include "simplexeikonalsolver.hpp"

namespace Feel {

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver = SimplexEikonalSolver >
class FastMarching: private LocalEikonalSolver< FunctionSpaceType >
{
    public:
        typedef FastMarching< FunctionSpaceType, LocalEikonalSolver > self_type;
        typedef std::shared_ptr< self_type > self_ptrtype;
        typedef LocalEikonalSolver< FunctionSpaceType > eikonal_solver_type;

        //--------------------------------------------------------------------//
        // Functionspace and mesh
        typedef FunctionSpaceType functionspace_type;
        typedef std::shared_ptr< functionspace_type > functionspace_ptrtype;
        typedef typename functionspace_type::element_type element_type;
        typedef typename functionspace_type::element_ptrtype element_ptrtype;

        static const uint16_type nDofPerElt = functionspace_type::fe_type::nDof;

        typedef typename functionspace_type::mesh_type mesh_type;
        typedef typename functionspace_type::mesh_ptrtype mesh_ptrtype;

        typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
        typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
        typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

        //--------------------------------------------------------------------//
        static constexpr uint16_type nRealDim = functionspace_type::nRealDim;
        using size_type = typename functionspace_type::size_type;
        typedef typename functionspace_type::value_type value_type;
        typedef typename node<value_type>::type node_type;
        typedef typename matrix_node<value_type>::type matrix_node_type;

        //--------------------------------------------------------------------//
        // Utility comparators
        template< typename T >
        struct greater_abs
        {
            constexpr bool operator()( T const& lhs, T const& rhs ) 
            {
                return std::abs( lhs ) > std::abs( rhs );
            }
        };
        template< typename Key, typename Val >
        struct greater_abs< std::pair< Key, Val > >
        {
            constexpr bool operator()( std::pair< Key, Val > const& lhs, std::pair< Key, Val > const& rhs ) 
            {
                return std::abs( lhs.second ) > std::abs( rhs.second );
            }
        };
        template< typename T >
        struct less_abs
        {
            constexpr bool operator()( T const& lhs, T const& rhs ) 
            {
                return std::abs( lhs ) < std::abs( rhs );
            }
        };
        template< typename Key, typename Val >
        struct less_abs< std::pair< Key, Val > >
        {
            constexpr bool operator()( std::pair< Key, Val > const& lhs, std::pair< Key, Val > const& rhs ) 
            {
                return std::abs( lhs.second ) < std::abs( rhs.second );
            }
        };
        //--------------------------------------------------------------------//
        template< typename Key, typename Val, 
                  typename Cmp = greater_abs< std::pair< Key, Val > >
                >
        class HeapMap
        {
            public:
                typedef Key key_type;
                typedef Val value_type;
                typedef std::pair< key_type, value_type > data_type;

            public:
                HeapMap( Cmp cmp = Cmp() ) : M_cmp( cmp ) {}

                auto size() const { return M_data.size(); }

                void insert( data_type const& data ) 
                {
                    M_data.push_back( data );
                    std::push_heap( M_data.begin(), M_data.end(), M_cmp );
                }
                template< typename InputIt >
                void insert( InputIt first, InputIt last )
                {
                    M_data.insert( M_data.end(), first, last );
                    std::make_heap( M_data.begin(), M_data.end(), M_cmp );
                }

                void insert_or_assign( data_type const& data )
                {
                    auto dataIt = std::find_if( std::begin( M_data ), std::end( M_data ),
                            [&data]( data_type const& d ) { return d.first == data.first; }
                            );
                    if( dataIt == std::end( M_data ) )
                    {
                        this->insert( data );
                    }
                    else
                    {
                        dataIt->second = data.second;
                        std::make_heap( M_data.begin(), M_data.end(), M_cmp );
                    }
                }

                template< typename Predicate >
                bool insert_or_assign_if( data_type const& data, Predicate p )
                {
                    auto dataIt = std::find_if( std::begin( M_data ), std::end( M_data ),
                            [&data]( data_type const& d ) { return d.first == data.first; }
                            );
                    if( dataIt == std::end( M_data ) )
                    {
                        this->insert( data );
                        return true;
                    }
                    else
                    {
                        if( p(*dataIt) )
                        {
                            dataIt->second = data.second;
                            std::make_heap( M_data.begin(), M_data.end(), M_cmp );
                            return true;
                        }
                        else
                            return false;
                    }
                }

                data_type front() const 
                {
                    return M_data.front();
                }
                data_type pop_front()
                {
                    data_type front = M_data.front();
                    std::pop_heap( M_data.begin(), M_data.end(), M_cmp );
                    M_data.pop_back();
                    return front;
                }

                void clear()
                {
                    M_data.clear();
                }

                bool empty() const { return M_data.empty(); }

                Cmp const& comp() const { return M_cmp; }

            private:
                std::vector< data_type > M_data;
                Cmp M_cmp;
        };

        typedef std::pair< size_type, value_type > pair_dof_value_type;
        typedef HeapMap< size_type, value_type > heap_type;

    public:
        //--------------------------------------------------------------------//
        // Constructor
        FastMarching( functionspace_ptrtype const& space );

        //--------------------------------------------------------------------//
        // Accessors
        functionspace_ptrtype const& functionSpace() const { return M_space; }
        mesh_ptrtype const& mesh() const { return this->functionSpace()->mesh(); }

        //--------------------------------------------------------------------//
        // Utility
        static bool greaterAbsValue( value_type const& a, value_type const& b )
        {
            return std::abs( a ) > std::abs( b );
        }
        static bool greaterAbsPairDofValue( pair_dof_value_type const& a, pair_dof_value_type const& b ) 
        {
            return self_type::greaterAbsValue( a.second, b.second );
        }

        //--------------------------------------------------------------------//
        // Result
        element_type run( element_type const& phi, range_elements_type const& rangeDone );

    private:
        using eikonal_solver_type::solveEikonal;
        
        void updateNeighborDofs( size_type dofId, element_type & sol );
        void updateCloseDofs( std::vector< size_type > const& dofCloseIds, std::vector< size_type > const& dofDoneIds, size_type eltId, element_type & sol );

        void marchNarrowBand( element_type & sol );

        void syncDofs( element_type & sol );

        element_type runImpl( element_type const& phi, range_elements_type const& rangeDone );

    private:
        functionspace_ptrtype M_space;

        std::map< size_type, std::set< rank_type > > M_dofSharedOnCluster;
        std::map< size_type, size_type > M_mapSharedDofGlobalClusterToGlobalProcess;

        std::vector< FastMarchingDofStatus > M_dofStatus;
        heap_type M_positiveCloseDofHeap;
        heap_type M_negativeCloseDofHeap;
        int M_nNewDofs;
};

} // namespace Feel


#endif // _FASTMARCHING_HPP
