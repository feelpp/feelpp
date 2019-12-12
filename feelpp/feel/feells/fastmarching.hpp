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
#include <feel/feells/reinit_fms_impl.hpp>

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

        std::vector< FastMarchingDofStatus > M_dofStatus;
        heap_type M_positiveCloseDofHeap;
        heap_type M_negativeCloseDofHeap;
};


template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
FastMarching< FunctionSpaceType, LocalEikonalSolver >::FastMarching(
        functionspace_ptrtype const& space ) :
    eikonal_solver_type( space ),
    M_space( space ),
    M_dofStatus( space->dof()->nLocalDofWithGhost() ),
    M_positiveCloseDofHeap(),
    M_negativeCloseDofHeap()
{
    auto const& dofTable = this->functionSpace()->dof();
    // Init shared dof map
    M_dofSharedOnCluster = dofTable->activeDofSharedOnCluster();
    for( size_type dofId = 0; dofId < dofTable->nLocalDofWithGhost(); ++dofId )
    {
        if( dofTable->dofGlobalProcessIsGhost( dofId ) )
        {
            size_type dofGCId = dofTable->mapGlobalProcessToGlobalCluster( dofId );
            rank_type procId = dofTable->procOnGlobalCluster( dofGCId );
            M_dofSharedOnCluster[dofId].insert( procId );
        }
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
typename FastMarching< FunctionSpaceType, LocalEikonalSolver >::element_type
FastMarching< FunctionSpaceType, LocalEikonalSolver >::run( element_type const& phi, range_elements_type const& rangeDone )
{
    return this->runImpl( phi, rangeDone );
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
typename FastMarching< FunctionSpaceType, LocalEikonalSolver >::element_type
FastMarching< FunctionSpaceType, LocalEikonalSolver >::runImpl( element_type const& phi, range_elements_type const& rangeDone )
{
    // Initialize return
    element_type sol = this->functionSpace()->element();
    sol.on( _range=elements( this->mesh() ), _expr=idv(phi) );

    // Initialize helper structures
    std::fill( M_dofStatus.begin(), M_dofStatus.end(), FastMarchingDofStatus::FAR );
    M_positiveCloseDofHeap.clear();
    M_negativeCloseDofHeap.clear();

    // Initialize DONE dofs
    auto itEltDone = boost::get<1>( rangeDone );
    auto enEltDone = boost::get<2>( rangeDone );
    for( ; itEltDone != enEltDone; ++itEltDone )
    {
        auto const eltDone = boost::unwrap_ref( *itEltDone );
        size_type const eltDoneId = eltDone.id();
        for( uint16_type j = 0; j < nDofPerElt; ++j )
        {
            size_type dofId = this->functionSpace()->dof()->localToGlobalId( eltDoneId, j );
            M_dofStatus[dofId] = FastMarchingDofStatus::DONE_FIX;
        }
    }
    // Initialize CLOSE (M_closeDofHeap) dofs
    itEltDone = boost::get<1>( rangeDone );
    for( ; itEltDone != enEltDone; ++itEltDone )
    {
        auto const eltDone = boost::unwrap_ref( *itEltDone );
        size_type const eltDoneId = eltDone.id();
        for( uint16_type j = 0; j < nDofPerElt; ++j )
        {
            size_type dofId = this->functionSpace()->dof()->localToGlobalId( eltDoneId, j );
            this->updateNeighborDofs( dofId, sol );
        }
    }

    // Perform parallel fast-marching
    bool closeDofHeapsAreEmptyOnAllProc = false;
    while( !closeDofHeapsAreEmptyOnAllProc )
    {
        this->marchNarrowBand( sol );
        this->syncDofs( sol );
        this->marchNarrowBand( sol );
        // Update closeDofHeapsAreEmptyOnAllProc
        bool closeDofHeapsAreEmpty = M_positiveCloseDofHeap.empty() && M_negativeCloseDofHeap.empty();
        closeDofHeapsAreEmptyOnAllProc = mpi::all_reduce(
                this->functionSpace()->worldComm(),
                closeDofHeapsAreEmpty,
                std::logical_and<bool>()
                );
    }

    return sol;
}


template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::updateNeighborDofs( size_type dofDoneId, element_type & sol )
{
    value_type const dofDoneVal = sol(dofDoneId);
    // Process elements containing dofDoneId
    auto const& [beg,end] = this->functionSpace()->dof()->globalDof( dofDoneId );
    for( auto eltIt = beg; eltIt != end; ++eltIt )
    {
        size_type const eltId = eltIt->second.elementId();
        auto const& elt = this->mesh()->element( eltId );
        rank_type const eltPid = elt.processId();

        // process element dofs
        std::vector< size_type > dofDoneIds, dofIdsToProcess;
        dofDoneIds.emplace_back( dofDoneId );
        bool hasDofToProcess = false;
        for( auto const& [lid,gid]: this->functionSpace()->dof()->localDof( eltId ) )
        {
            size_type const dofId = gid.index();
            if( dofId == dofDoneId )
                continue;
            if( M_dofStatus[dofId] != FastMarchingDofStatus::DONE_FIX 
                    && less_abs<value_type>()( dofDoneVal, sol(dofId) ) )
            {
                dofIdsToProcess.emplace_back( dofId );
                hasDofToProcess = true;
            }
            else if( M_dofStatus[dofId] & FastMarchingDofStatus::DONE )
            {
                dofDoneIds.emplace_back( dofId );
            }
        }

        if( hasDofToProcess )
        {
            this->updateCloseDofs( dofIdsToProcess, dofDoneIds, eltId, sol );
        }
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::updateCloseDofs( std::vector< size_type > const& dofCloseIds, std::vector< size_type > const& dofDoneIds, size_type eltId, element_type & sol )
{
    auto const closeDofsIdVals = this->solveEikonal( dofCloseIds, dofDoneIds, eltId, sol );
    for( auto const& [closeDofId,closeDofVal]: closeDofsIdVals )
    {
        if( less_abs<value_type>()( closeDofVal, sol(closeDofId) ) )
        {
            if( sol(dofDoneIds[0]) > 0. )
            {
                //value_type signedCloseDofVal = ( sol(dofDoneIds[0]) > 0. ) ? std::abs( closeDofVal ) : -std::abs( closeDofVal );
                sol(closeDofId) = closeDofVal;
                M_positiveCloseDofHeap.insert_or_assign( pair_dof_value_type(closeDofId, closeDofVal) );
                M_dofStatus[closeDofId] = FastMarchingDofStatus::CLOSE_NEW;

                //std::cout << "["<<this->mesh()->worldCommPtr()->localRank()<<"]" 
                //<< "updating dof(" << closeDofId << "," << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( closeDofId ) << ")"
                //<< " with value " << sol(closeDofId)
                //<< " using dofs {";
                //for( auto dofId: dofDoneIds )
                //std::cout << "(" << dofId << "," << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster(dofId) << "); ";
                //std::cout << "}" << std::endl;
            }
            else
            {
                //value_type signedCloseDofVal = ( sol(dofDoneIds[0]) > 0. ) ? std::abs( closeDofVal ) : -std::abs( closeDofVal );
                sol(closeDofId) = -closeDofVal;
                M_negativeCloseDofHeap.insert_or_assign( pair_dof_value_type(closeDofId, -closeDofVal) );
                M_dofStatus[closeDofId] = FastMarchingDofStatus::CLOSE_NEW;

                //std::cout << "["<<this->mesh()->worldCommPtr()->localRank()<<"]" 
                //<< "updating dof(" << closeDofId << "," << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( closeDofId ) << ")"
                //<< " with value " << sol(closeDofId)
                //<< " using dofs {";
                //for( auto dofId: dofDoneIds )
                //std::cout << "(" << dofId << "," << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster(dofId) << "); ";
                //std::cout << "}" << std::endl;
            }
        }
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::marchNarrowBand( element_type & sol )
{
    // Perform local fast-marching
    while( !M_positiveCloseDofHeap.empty() )
    {
        auto [nextDofDoneId, nextDofDoneValue] = M_positiveCloseDofHeap.pop_front();
        sol(nextDofDoneId) = nextDofDoneValue;
        if( M_dofStatus[nextDofDoneId] != FastMarchingDofStatus::DONE_OLD )
            M_dofStatus[nextDofDoneId] = FastMarchingDofStatus::DONE_NEW;
        this->updateNeighborDofs( nextDofDoneId, sol );
    }
    while( !M_negativeCloseDofHeap.empty() )
    {
        auto [nextDofDoneId, nextDofDoneValue] = M_negativeCloseDofHeap.pop_front();
        sol(nextDofDoneId) = nextDofDoneValue;
        if( M_dofStatus[nextDofDoneId] != FastMarchingDofStatus::DONE_OLD )
            M_dofStatus[nextDofDoneId] = FastMarchingDofStatus::DONE_NEW;
        this->updateNeighborDofs( nextDofDoneId, sol );
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::syncDofs( element_type & sol )
{
    // Dof table
    auto const& dofTable = sol.dof();
    // If sequential return
    if( dofTable->worldCommPtr()->localSize() == 1 )
        return;
    // MPI comm
    rank_type const localPid = dofTable->worldCommPtr()->localRank();
    int nRequests = 2 * dofTable->neighborSubdomains().size();
    mpi::request * mpiRequests = new mpi::request[nRequests];
    int cntRequests = 0;

    // Send ghost dofs to the unique active dof
    std::map< rank_type, std::vector< std::pair<size_type, value_type> > > dataToSend, dataToRecv;
    for( auto const& dofShared: M_dofSharedOnCluster )
    {
        size_type const dofId = dofShared.first;
        if( M_dofStatus[dofId] & FastMarchingDofStatus::NEW )
        {
            // Set status to OLD
            M_dofStatus[dofId] &= ~FastMarchingDofStatus::NEW;
            M_dofStatus[dofId] |= FastMarchingDofStatus::OLD;

            size_type dofGCId = dofTable->mapGlobalProcessToGlobalCluster( dofId );
            value_type dofVal = sol(dofId);
            for( rank_type p: dofShared.second )
                dataToSend[p].emplace_back( dofGCId, dofVal );
        }
    }
    // Perform communications
    cntRequests = 0;
    for( rank_type p: dofTable->neighborSubdomains() )
    {
        //std::cout << "["<<localPid<<"]" << "sending to " << p << ": ";
        //for( auto const& dataS: dataToSend[p] )
        //{
            //std::cout << " (" << dataS.first << "," << dataS.second << "), ";
        //}
        //std::cout << std::endl;
        mpiRequests[cntRequests++] = dofTable->worldCommPtr()->localComm().isend( p, 0, dataToSend[p] );
        mpiRequests[cntRequests++] = dofTable->worldCommPtr()->localComm().irecv( p, 0, dataToRecv[p] );
    }
    mpi::wait_all( mpiRequests, mpiRequests + nRequests );
    // Update the value of active dofs if needed
    for( auto const& dataR: dataToRecv )
    {
        for( auto const& dofGCIdVal: dataR.second )
        {
            // Find received dof global process id
            size_type const dofGCId = dofGCIdVal.first;
            auto resSearchDof = dofTable->searchGlobalProcessDof( dofGCId );
            DCHECK( boost::get<0>( resSearchDof ) ) << "[" << localPid << "]" << " dof " << dofGCId << " not found\n";
            size_type const dofId = boost::get<1>( resSearchDof );
            //size_type const dofId = dofGCId - dofTable->firstDofGlobalCluster();
            // Update current value with received ghost values if needed
            value_type valCurrent = sol(dofId);
            value_type valNew = dofGCIdVal.second;
            if( less_abs<value_type>()( valNew, valCurrent ) )
            {
                heap_type & closeDofHeap = ( valNew > 0. ) ? M_positiveCloseDofHeap : M_negativeCloseDofHeap;
                sol.set( dofId, valNew );
                closeDofHeap.insert_or_assign( pair_dof_value_type( dofId, valNew ) );
                M_dofStatus[dofId] = FastMarchingDofStatus::DONE_OLD;
                //std::cout << "["<<localPid<<"]" << "updating active dof(" << dofId << "," << dofGCId << ")" 
                    //<< std::endl;
            }
        }
    }

    // Clean
    delete [] mpiRequests;
}

} // namespace Feel


#endif // _FASTMARCHING_HPP
