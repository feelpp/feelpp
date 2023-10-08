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

        static inline const uint16_type nDofPerElt = functionspace_type::fe_type::nDof;

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
                // We also add a (boolean) tag to track deleted entries
                typedef std::pair< data_type, bool > container_data_type;
                typedef std::unique_ptr< container_data_type > container_data_ptrtype;
                typedef std::vector< container_data_ptrtype > container_type;

                template< typename C >
                class CmpAdaptor
                {
                public:
                    CmpAdaptor( C c ) : M_c( c ) {}
                    constexpr bool operator()( 
                            container_data_ptrtype const& lhs, 
                            container_data_ptrtype const& rhs ) 
                    {
                        return M_c( lhs->first, rhs->first );
                    }
                private:
                    C M_c;
                };
                typedef CmpAdaptor<Cmp> cmp_type;

            public:
                HeapMap( Cmp cmp = Cmp() ) : 
                    M_data(), M_cmp( cmp ), M_validEntriesPtr(), M_size( 0 ) 
                {}

                void insert( data_type const& data ) 
                {
                    M_data.emplace_back( new container_data_type( data, true ) );
                    M_validEntriesPtr[data.first] = M_data.back().get();
                    std::push_heap( M_data.begin(), M_data.end(), M_cmp );
                    ++M_size;
                }
                template< typename InputIt >
                void insert( InputIt first, InputIt last )
                {
                    for( auto it = first; it != last; ++it )
                    {
                        M_data.emplace_back( new container_data_type( *it, true ) );
                        M_validEntriesPtr[it->first] = M_data.back().get();
                    }
                    std::make_heap( M_data.begin(), M_data.end(), M_cmp );
                    M_size += std::distance( first, last );
                }

                void insert_or_assign( data_type const& data )
                {
                    auto dataIt = M_validEntriesPtr.find( data.first );
                    if( dataIt == M_validEntriesPtr.end() )
                    {
                        this->insert( data );
                    }
                    else
                    {
                        dataIt->second->second = false;
                        --M_size;
                        this->insert( data );
                    }
                }

                template< typename Predicate >
                bool insert_or_assign_if( data_type const& data, Predicate p )
                {
                    auto dataIt = M_validEntriesPtr.find( data.first );
                    if( dataIt == M_validEntriesPtr.end() )
                    {
                        this->insert( data );
                        return true;
                    }
                    else
                    {
                        if( p(*dataIt) )
                        {
                            dataIt->second->second = false;
                            --M_size;
                            this->insert( data );
                            return true;
                        }
                        else
                            return false;
                    }
                }

                data_type front() const 
                {
                    if( M_data.front()->second ) // Entry is valid
                    {
                        return M_data.front()->first;
                    }
                    else // Entry was deleted, remove it
                    {
                        container_type & data = const_cast<container_type &>( M_data );
                        std::pop_heap( data.begin(), data.end(), M_cmp );
                        data.pop_back();
                        return this->front();
                    }
                }
                data_type pop_front()
                {
                    container_data_type const f = *M_data.front();
                    std::pop_heap( M_data.begin(), M_data.end(), M_cmp );
                    M_data.pop_back();
                    if( f.second ) // Entry is valid
                    {
                        M_validEntriesPtr.erase( f.first.first );
                        --M_size;
                        return f.first;
                    }
                    else // Entry was deleted, return next
                        return this->pop_front();
                }

                void clear()
                {
                    M_data.clear();
                    M_validEntriesPtr.clear();
                    M_size = 0;
                }

                auto size() const { return M_size; }
                bool empty() const { return M_size == 0; }

                cmp_type const& comp() const { return M_cmp; }

            private:
                container_type M_data;
                cmp_type M_cmp;
                std::unordered_map< key_type, container_data_type * > M_validEntriesPtr;
                // Track size accounting only valid entries
                size_type M_size;
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
        // Options
        /*
         * Positive narrow band width (negative for infinite width)
         */
        value_type positiveNarrowBandWidth() const { return M_positiveNarrowBandWidth; }
        void setPositiveNarrowBandWidth( const value_type & width ) { M_positiveNarrowBandWidth = width; }
        /*
         * Negative narrow band width (negative for infinite width)
         */
        value_type negativeNarrowBandWidth() const { return M_negativeNarrowBandWidth; }
        void setNegativeNarrowBandWidth( const value_type & width ) { M_negativeNarrowBandWidth = width; }
        /*
         * Set positive and negative widths
         */
        void setNarrowBandWidth( const value_type & width ) { this->setPositiveNarrowBandWidth( width ); this->setNegativeNarrowBandWidth( width ); }
        /*
         * Local marching stride before parallel update (negative for infinite stride, ie full local marching)
         */
        value_type stride() const { return M_stride; }
        void setStride( const value_type & s ) { M_stride = s; }

        //--------------------------------------------------------------------//
        // Result
        element_type run( element_type const& phi, range_elements_type const& rangeDone );

    private:
        using eikonal_solver_type::solveEikonal;
        
        void updateNeighborDofs( size_type dofId, element_type & sol );
        void updateCloseDofs( std::vector< size_type > const& dofCloseIds, std::vector< size_type > const& dofDoneIds, size_type eltId, element_type & sol );

        void marchLocalNarrowBand( element_type & sol, const value_type & positiveBound, const value_type & negativeBound );
        template< bool HasBound >
        void marchLocalSignedNarrowBand( element_type & sol, heap_type & heap, const value_type & bound );

        void syncDofs( element_type & sol, const value_type & positiveBound = -1., const value_type & negativeBound = -1. );

        element_type runImpl( element_type const& phi, range_elements_type const& rangeDone );

    private:
        functionspace_ptrtype M_space;

        value_type M_positiveNarrowBandWidth = -1.;
        value_type M_negativeNarrowBandWidth = -1.;

        value_type M_stride = -1.;

        std::map< size_type, std::set< rank_type > > M_dofSharedOnCluster;
        std::map< size_type, size_type > M_mapSharedDofGlobalClusterToGlobalProcess;

        std::vector< FastMarchingDofStatus > M_dofStatus;
        heap_type M_positiveCloseDofHeap;
        heap_type M_negativeCloseDofHeap;
        int M_nNewDofs = 0;
};

} // namespace Feel


#endif // _FASTMARCHING_HPP
