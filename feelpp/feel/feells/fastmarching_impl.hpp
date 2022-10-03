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


#include <feel/feeldiscr/syncdofs.hpp>

#include "fastmarching.hpp"

//#define DEBUG_FM_COUT

namespace Feel {


template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
FastMarching< FunctionSpaceType, LocalEikonalSolver >::FastMarching(
        functionspace_ptrtype const& space ) :
    eikonal_solver_type( space ),
    M_space( space ),
    M_dofStatus( space->dof()->nLocalDofWithGhost() ),
    M_positiveCloseDofHeap(),
    M_negativeCloseDofHeap(),
    M_nNewDofs( 0 )
{
    auto const& dofTable = this->functionSpace()->dof();
    // Init shared dof map
    // First find active dofs with sharing procs
    M_dofSharedOnCluster = dofTable->activeDofSharedOnCluster();
    // Then send this information with sharing procs so that
    // they also know which are the other ghost-owning procs
    rank_type const localPid = dofTable->worldCommPtr()->localRank();
    int nRequests = 2 * dofTable->neighborSubdomains().size();
    mpi::request * mpiRequests = new mpi::request[nRequests];
    std::map< rank_type, std::vector< std::pair< size_type, std::set<rank_type> > > > dataToSend, dataToRecv;
    for( auto const& activeDofShared: M_dofSharedOnCluster )
    {
        size_type const dofId = activeDofShared.first;
        size_type dofGCId = dofTable->mapGlobalProcessToGlobalCluster( dofId );
        for( rank_type const p : activeDofShared.second )
            dataToSend[p].emplace_back( dofGCId, activeDofShared.second );
        // Update active dofs globalClusterToGlobalProcess map
        M_mapSharedDofGlobalClusterToGlobalProcess[dofGCId] = dofId;
    }
    // Update ghost dofs globalClusterToGlobalProcess map
    size_type const ghostDofIdStart = dofTable->nLocalDofWithoutGhost();
    size_type const ghostDofIdEnd = dofTable->nLocalDofWithGhost();
    for( size_type k = ghostDofIdStart ; k < ghostDofIdEnd ; ++k )
    {
        M_mapSharedDofGlobalClusterToGlobalProcess[dofTable->mapGlobalProcessToGlobalCluster(k)] = k;
    }

    // Send shared dofs
    int cntRequests = 0;
    for( rank_type const p: dofTable->neighborSubdomains() )
    {
        mpiRequests[cntRequests++] = dofTable->worldCommPtr()->localComm().isend( p, 0, dataToSend[p] );
        mpiRequests[cntRequests++] = dofTable->worldCommPtr()->localComm().irecv( p, 0, dataToRecv[p] );
    }
    mpi::wait_all( mpiRequests, mpiRequests + nRequests );
    for( auto const& dataR: dataToRecv )
    {
        for( auto const& dofGCIdPids: dataR.second )
        {
            size_type const dofGCId = dofGCIdPids.first;
            size_type const dofId = M_mapSharedDofGlobalClusterToGlobalProcess[dofGCId];
#if !defined(NDEBUG)
            auto resSearchDof = dofTable->searchGlobalProcessDof( dofGCId );
            DCHECK( boost::get<0>( resSearchDof ) ) << "[" << localPid << "]" << " dof " << dofGCId << " not found\n";
            size_type const dofId2 = boost::get<1>( resSearchDof );
            CHECK( dofId == dofId2 ) << "[" << localPid << "]" << "dof id must be the same : " << dofId << " vs " << dofId2 << "(" << dofGCId << "," << dofTable->firstDofGlobalCluster() << ")" << std::endl;
#endif
            M_dofSharedOnCluster[dofId].insert( dataR.first );
            for( rank_type const p : dofGCIdPids.second )
                if( p != localPid )
                    M_dofSharedOnCluster[dofId].insert( p );
        }
    }
#ifdef DEBUG_FM_COUT
    std::cout << "["<<dofTable->worldCommPtr()->localRank()<<"]" << "dofSharedOnCluster = ";
    for( auto const& dofShared: M_dofSharedOnCluster )
    {
        std::cout << " (" << dofShared.first << ", {" ;
        for( rank_type p: dofShared.second )
            std::cout << p << ",";
        std::cout << "} )";
    }
    std::cout << std::endl;
#endif
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
    sol = phi;

    // Initialize helper structures
    std::fill( M_dofStatus.begin(), M_dofStatus.end(), FastMarchingDofStatus::FAR );
    M_positiveCloseDofHeap.clear();
    M_negativeCloseDofHeap.clear();

    // Initialize DONE dofs
    auto itEltDone = rangeDone.template get<1>();
    auto enEltDone = rangeDone.template get<2>();
    for( ; itEltDone != enEltDone; ++itEltDone )
    {
        auto const eltDone = boost::unwrap_ref( *itEltDone );
        size_type const eltDoneId = eltDone.id();
#ifdef DEBUG_FM_COUT
        std::cout << "["<<this->mesh()->worldCommPtr()->localRank()<<"]" 
            << "fixing elt " << eltDoneId << " with dofs { ";
#endif
        for( auto const& lDof: this->functionSpace()->dof()->localDof( eltDoneId ) )
        {
            size_type dofId = lDof.second.index();
            M_dofStatus[dofId] = FastMarchingDofStatus::DONE_FIX;
#ifdef DEBUG_FM_COUT
            std::cout << "(" 
                << lDof.first.localDof() << ","
                << dofId << ","
                << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( dofId ) << "; "
                << sol( dofId ) << "); ";
#endif
        }
#ifdef DEBUG_FM_COUT
        std::cout << " }" << std::endl;
#endif
    }
    Feel::syncDofs( M_dofStatus, *sol.dof(), rangeDone,
            []( FastMarchingDofStatus curStatus, std::set<FastMarchingDofStatus> ghostVals ) 
            -> FastMarchingDofStatus
            {
                if( curStatus == FastMarchingDofStatus::DONE_FIX )
                    return FastMarchingDofStatus::DONE_FIX;
                for( auto const& ghostStatus: ghostVals )
                {
                    if( ghostStatus == FastMarchingDofStatus::DONE_FIX )
                    {
                        return FastMarchingDofStatus::DONE_FIX;
                    }
                }
                return curStatus;
            }
            );
    // Initialize CLOSE (M_closeDofHeap) dofs
    for( size_type dofId = 0; dofId < this->functionSpace()->nLocalDofWithGhost(); ++dofId )
    {
        if( M_dofStatus[dofId] & FastMarchingDofStatus::DONE )
        {
            this->updateNeighborDofs( dofId, sol );
        }
    }

    bool hasPositiveBound = M_positiveNarrowBandWidth > 0.;
    value_type positiveBound = -1.;
    bool hasNegativeBound = M_negativeNarrowBandWidth > 0.;
    value_type negativeBound = -1.;
    bool hasStride = M_stride > 0.;

    // Perform parallel fast-marching
    value_type minPositiveAbsGlobalValue = -1.;
    value_type minNegativeAbsGlobalValue = -1.;
    while( true )
    {
        if( hasPositiveBound )
        {
            // Compute global absolute positive min values
            value_type minPositiveAbsLocalValue = ( !M_positiveCloseDofHeap.empty() ) ? std::abs( M_positiveCloseDofHeap.front().second ) : M_positiveNarrowBandWidth;
            minPositiveAbsGlobalValue = mpi::all_reduce(
                this->functionSpace()->worldComm(),
                minPositiveAbsLocalValue,
                mpi::minimum<value_type>()
                );
            positiveBound = hasStride ? std::min( minPositiveAbsGlobalValue + M_stride, M_positiveNarrowBandWidth ) : M_positiveNarrowBandWidth;
        }
        if( hasNegativeBound )
        {
            // Compute global absolute negative min values
            value_type minNegativeAbsLocalValue = ( !M_negativeCloseDofHeap.empty() ) ? std::abs( M_negativeCloseDofHeap.front().second ) : M_negativeNarrowBandWidth;
            minNegativeAbsGlobalValue = mpi::all_reduce(
                this->functionSpace()->worldComm(),
                minNegativeAbsLocalValue,
                mpi::minimum<value_type>()
                );
            negativeBound = hasStride ? std::min( minNegativeAbsGlobalValue + M_stride, M_negativeNarrowBandWidth ) : M_negativeNarrowBandWidth;
        }

        // Update maxNumberOfNewDofsOnAllProc
        int maxNumberOfNewDofsOnAllProc = mpi::all_reduce(
                this->functionSpace()->worldComm(),
                M_nNewDofs,
                mpi::maximum<int>()
                );
        if( !hasPositiveBound )
        {
            if( !hasNegativeBound )
            {
                // Check if both positive and negative heaps are empty on all proc
                bool closeDofHeapsAreEmpty = M_positiveCloseDofHeap.empty() && M_negativeCloseDofHeap.empty();
                bool closeDofHeapsAreEmptyOnAllProc = mpi::all_reduce(
                        this->functionSpace()->worldComm(),
                        closeDofHeapsAreEmpty,
                        std::logical_and<bool>()
                        );
                if( closeDofHeapsAreEmptyOnAllProc && maxNumberOfNewDofsOnAllProc == 0 )
                    break;
            }
            else
            {
                // Check positive heap is empty on all proc
                bool closeDofPositiveHeapIsEmpty = M_positiveCloseDofHeap.empty();
                bool closeDofPositiveHeapIsEmptyOnAllProc = mpi::all_reduce(
                        this->functionSpace()->worldComm(),
                        closeDofPositiveHeapIsEmpty,
                        std::logical_and<bool>()
                        );
                if( minNegativeAbsGlobalValue >= M_negativeNarrowBandWidth
                        && closeDofPositiveHeapIsEmptyOnAllProc
                        && maxNumberOfNewDofsOnAllProc == 0 )
                    break;
            }
        }
        else
        {
            if( !hasNegativeBound )
            {
                // Check negative heap is empty on all proc
                bool closeDofNegativeHeapIsEmpty = M_negativeCloseDofHeap.empty();
                bool closeDofNegativeHeapIsEmptyOnAllProc = mpi::all_reduce(
                        this->functionSpace()->worldComm(),
                        closeDofNegativeHeapIsEmpty,
                        std::logical_and<bool>()
                        );
                if( minPositiveAbsGlobalValue >= M_positiveNarrowBandWidth
                        && closeDofNegativeHeapIsEmptyOnAllProc
                        && maxNumberOfNewDofsOnAllProc == 0 )
                    break;
            }
            else
            {
                if( ( minPositiveAbsGlobalValue >= M_positiveNarrowBandWidth ) 
                        && ( minNegativeAbsGlobalValue >= M_negativeNarrowBandWidth ) 
                        && maxNumberOfNewDofsOnAllProc == 0 )
                    break;
            }
        }

        this->marchLocalNarrowBand( sol, positiveBound, negativeBound );
        this->syncDofs( sol, positiveBound, negativeBound );
        // The second marchLocalNarrowBand is only needed when the stride is finite
        if( hasStride )
            this->marchLocalNarrowBand( sol, positiveBound, negativeBound );
    }

    return sol;
}


template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::updateNeighborDofs( size_type dofDoneId, element_type & sol )
{
    value_type const dofDoneVal = sol(dofDoneId);
    int dofDoneSgn = (0. < dofDoneVal) - (dofDoneVal < 0.);
    // Process elements containing dofDoneId
#ifdef DEBUG_FM_COUT
    std::cout << "["<<this->mesh()->worldCommPtr()->localRank()<<"]" 
        << "processing neighbors of dof(" << dofDoneId << "," 
        << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( dofDoneId ) << "; "
        << dofDoneVal
        << ") : ";
#endif
    auto const& [beg,end] = this->functionSpace()->dof()->globalDof( dofDoneId );
    for( auto eltIt = beg; eltIt != end; ++eltIt )
    {
        size_type const eltId = eltIt->second.elementId();
        auto const& elt = this->mesh()->element( eltId );
        rank_type const eltPid = elt.processId();

#ifdef DEBUG_FM_COUT
        std::cout << "elt " << eltId << " with dofs { ";
#endif
        // process element dofs
        std::vector< size_type > dofDoneIds, dofCloseIdsToProcess, dofDoneIdsToProcess;
        dofDoneIds.emplace_back( dofDoneId );
        for( auto const& [lid,gid]: this->functionSpace()->dof()->localDof( eltId ) )
        {
            size_type const dofId = gid.index();
            if( dofId == dofDoneId )
                continue;
            // CLOSE or DONE neighbors are only considered
            // if they have the same sign than the dofDone
            int dofSgn = (0. < sol(dofId)) - (sol(dofId) < 0.);
            if( dofSgn != dofDoneSgn && dofDoneSgn != 0.
                    && !(M_dofStatus[dofId] & FastMarchingDofStatus::FAR) )
                continue;
#ifdef DEBUG_FM_COUT
            std::cout << "(" << dofId << "," 
                << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster(dofId) << "; "
                << int(M_dofStatus[dofId]) << ", "
                << sol(dofId) << ", "
                << int(less_abs<value_type>()( dofDoneVal, sol(dofId) ))
                << "); ";
#endif
            if( M_dofStatus[dofId] != FastMarchingDofStatus::DONE_FIX 
                    && less_abs<value_type>()( dofDoneVal, sol(dofId) ) )
            {
                if( M_dofStatus[dofId] & FastMarchingDofStatus::DONE )
                {
                    dofDoneIdsToProcess.emplace_back( dofId );
                }
                else
                {
                    dofCloseIdsToProcess.emplace_back( dofId );
                }
            }
            if( M_dofStatus[dofId] & FastMarchingDofStatus::DONE )
            {
                dofDoneIds.emplace_back( dofId );
            }
        }

#ifdef DEBUG_FM_COUT
        std::cout << " };" << std::endl;
#endif

        if( dofCloseIdsToProcess.size() != 0 )
            this->updateCloseDofs( dofCloseIdsToProcess, dofDoneIds, eltId, sol );
        for( size_type const dofDoneIdToProcess : dofDoneIdsToProcess )
        {
            std::vector< size_type > dofIdsToProcess = { dofDoneIdToProcess };
            std::vector< size_type > dofDoneIdsWithoutIdToProcess( dofDoneIds.size() - 1 );
            std::remove_copy( dofDoneIds.begin(), dofDoneIds.end(), dofDoneIdsWithoutIdToProcess.begin(), dofDoneIdToProcess );
            this->updateCloseDofs( dofIdsToProcess, dofDoneIdsWithoutIdToProcess, eltId, sol );
        }
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::updateCloseDofs( std::vector< size_type > const& dofIdsToProcess, std::vector< size_type > const& dofDoneIds, size_type eltId, element_type & sol )
{
    auto const closeDofsIdVals = this->solveEikonal( dofIdsToProcess, dofDoneIds, eltId, sol );
    for( auto const& [closeDofId,closeDofVal]: closeDofsIdVals )
    {
#ifdef DEBUG_FM_COUT
        int dofDoneSgn = (0. < sol(dofDoneIds[0])) - (sol(dofDoneIds[0]) < 0.);
        std::cout << "["<<this->mesh()->worldCommPtr()->localRank()<<"]" 
            << "updating dof(" << closeDofId << "," 
            << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( closeDofId ) << ")"
            << " with value " << sol(closeDofId)
            //<< " to value " << closeDofVal << "(" << ( (sol(dofDoneIds[0]) > 0.) ? "+":"-" ) << ")"
            << " to value " << closeDofVal << "(" << dofDoneSgn << ")"
            << " using dofs {";
        for( auto dofId: dofDoneIds )
            std::cout << "(" << dofId << "," 
                << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster(dofId) << "; "
                << sol(dofId) << "); ";
        std::cout << "} ";
        if( less_abs<value_type>()( closeDofVal, sol(closeDofId) ) )
            std::cout << "accept";
        std::cout << std::endl;
#endif
        if( less_abs<value_type>()( closeDofVal, sol(closeDofId) ) )
        {
            // we look for the first non-zero sign of done dofs
            int dofSgn = 0;
            for( size_type dofDoneId: dofDoneIds )
            {
                if( sol(dofDoneId) > 0. )
                {
                    dofSgn = 1;
                    break;
                }
                else if( sol(dofDoneId) < 0. )
                {
                    dofSgn = -1;
                    break;
                }
            }
            if( dofSgn >= 0 ) // positive or no signed done dof (default to positive)
            {
                sol(closeDofId) = closeDofVal;
                M_positiveCloseDofHeap.insert_or_assign( pair_dof_value_type(closeDofId, closeDofVal) );
            }
            else if( dofSgn < 0 )
            {
                sol(closeDofId) = -closeDofVal;
                M_negativeCloseDofHeap.insert_or_assign( pair_dof_value_type(closeDofId, -closeDofVal) );
            }

            M_dofStatus[closeDofId] = FastMarchingDofStatus::CLOSE_NEW;
        }
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::marchLocalNarrowBand( element_type & sol, const value_type & positiveBound, const value_type & negativeBound )
{
    // Perform local positive and negative fast-marching
    if( positiveBound > 0. )
        this->marchLocalSignedNarrowBand<true>( sol, M_positiveCloseDofHeap, positiveBound );
    else
        this->marchLocalSignedNarrowBand<false>( sol, M_positiveCloseDofHeap, positiveBound );

    if( negativeBound > 0. )
        this->marchLocalSignedNarrowBand<true>( sol, M_negativeCloseDofHeap, negativeBound );
    else
        this->marchLocalSignedNarrowBand<false>( sol, M_negativeCloseDofHeap, negativeBound );
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
template< bool HasBound >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::marchLocalSignedNarrowBand( element_type & sol, heap_type & heap, const value_type & bound )
{
    if constexpr( HasBound )
        DCHECK( bound > 0. ) << "bound should be positive";

    // Perform local fast-marching
    while( !heap.empty() )
    {
        auto [nextDofDoneId, nextDofDoneValue] = heap.pop_front();
        sol(nextDofDoneId) = nextDofDoneValue;
        if constexpr( HasBound )
        {
            if( std::abs( nextDofDoneValue ) > bound )
                break;
        }
        if( M_dofStatus[nextDofDoneId] != FastMarchingDofStatus::DONE_OLD )
            M_dofStatus[nextDofDoneId] = FastMarchingDofStatus::DONE_NEW;
#ifdef DEBUG_FM_COUT
        std::cout << "["<<this->mesh()->worldCommPtr()->localRank()<<"]" 
            << "updating dof(" << nextDofDoneId << "," 
            << this->functionSpace()->dof()->mapGlobalProcessToGlobalCluster( nextDofDoneId ) << ")"
            << " with value " << sol(nextDofDoneId)
            << " and status " << int(M_dofStatus[nextDofDoneId])
            << std::endl;
#endif
        this->updateNeighborDofs( nextDofDoneId, sol );
    }
}

template< typename FunctionSpaceType, template < typename > class LocalEikonalSolver >
void
FastMarching< FunctionSpaceType, LocalEikonalSolver >::syncDofs( element_type & sol, const value_type & positiveBound, const value_type & negativeBound )
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

    // Reset local number of NEW dofs
    M_nNewDofs = 0;

    // Send ghost dofs to the unique active dof
    std::map< rank_type, std::vector< std::pair<size_type, value_type> > > dataToSend, dataToRecv;
    for( auto const& dofShared: M_dofSharedOnCluster )
    {
        size_type const dofId = dofShared.first;
        if( M_dofStatus[dofId] & FastMarchingDofStatus::NEW )
        {
            // Update local number of NEW dofs
            M_nNewDofs += 1;
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
#ifdef DEBUG_FM_COUT 
        std::cout << "["<<localPid<<"]" << "sending to " << p << ": ";
        for( auto const& dataS: dataToSend[p] )
        {
            std::cout << " (" << dataS.first << "," << dataS.second << "), ";
        }
        std::cout << std::endl;
#endif
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
            size_type const dofId = M_mapSharedDofGlobalClusterToGlobalProcess[dofGCId];
#if !defined(NDEBUG)
            auto resSearchDof = dofTable->searchGlobalProcessDof( dofGCId );
            DCHECK( boost::get<0>( resSearchDof ) ) << "[" << localPid << "]" << " dof " << dofGCId << " not found\n";
            size_type const dofId2 = boost::get<1>( resSearchDof );
            CHECK( dofId == dofId2 ) << "[" << localPid << "]" << "dof id must be the same : " << dofId << " vs " << dofId2 << "(" << dofGCId << "," << dofTable->firstDofGlobalCluster() << ")" << std::endl;
#endif
            // Update current value with received ghost values if needed
            value_type valCurrent = sol(dofId);
            value_type valNew = dofGCIdVal.second;
            if( less_abs<value_type>()( valNew, valCurrent ) )
            {
                heap_type & closeDofHeap = ( valNew > 0. ) ? M_positiveCloseDofHeap : M_negativeCloseDofHeap;
                sol.set( dofId, valNew );
                closeDofHeap.insert_or_assign( pair_dof_value_type( dofId, valNew ) );
                if( ( positiveBound > 0. && valNew > positiveBound )
                 || ( negativeBound > 0. && valNew < 0. && std::abs(valNew) > negativeBound ) )
                    M_dofStatus[dofId] = FastMarchingDofStatus::CLOSE_OLD;
                else
                    M_dofStatus[dofId] = FastMarchingDofStatus::DONE_OLD;
#ifdef DEBUG_FM_COUT
                std::cout << "["<<localPid<<"]" << "updating active dof(" << dofId << "," << dofGCId << ")" 
                    << std::endl;
#endif
            }
        }
    }

    // Clean
    delete [] mpiRequests;
}

} // namespace Feel

