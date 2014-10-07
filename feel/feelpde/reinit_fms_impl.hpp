/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Winkelmann
              Vincent Doyeux
   Date     : Mon Feb 24 15:18:07 2014

   Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_REINIT_FMS_IMPL_HPP
#define FEELPP_REINIT_FMS_IMPL_HPP 1

#define FEELPP_INSTANTIATE_FMS 1

#include <feel/feelpde/reinit_fms.hpp>

#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operators.hpp>

//#define FM_EXPORT 1

#if defined(FM_EXPORT)
#include <feel/feelfilters/exporter.hpp>
#endif

namespace Feel
{

template<typename FunctionSpaceType, typename periodicity_type>
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                  periodicity_type __periodicity)
    :
    M_functionspace( __functionspace ),
    checkStatus(backend()->newVector(M_functionspace)),
    valueAtClose(backend()->newVector(M_functionspace)),
    M_periodicity(__periodicity),
    M_neighbors(),
    M_coords( __functionspace->dof()->nDof() ),
    M_translation( __periodicity.translation() ),
    firstDof( M_functionspace->dof()->firstDofGlobalCluster() ),
    M_nbDofTag1(0)
{

    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    auto it = M_functionspace->mesh()->beginElementWithProcessId();
    auto en = M_functionspace->mesh()->endElementWithProcessId();

    // create the first neighbours data structure (M_neighbors)
    for ( ; it!=en ; ++it )
        {
            std::vector<size_type> indices( ndofv );
            for ( uint16_type j = 0; j < ndofv;++j )
                {
                  size_type index = M_functionspace->dof()->localToGlobal(*it, j, 0).index();

                  if (M_functionspace->dof()->dofGlobalProcessIsGhost( index ))
                    M_ghostClusterToProc[ processorToCluster( index ) ] = index;

                    indices[j] = index;
                    for ( uint16_type c = 0; c < Dim; ++c )
                      M_coords[index][c] = M_functionspace->dof()->dofPoint( index ).template get<0>()[c];
                }
            for ( uint16_type j = 0; j < ndofv;++j )
                for ( uint16_type k = j+1; k < ndofv;++k )
                    {
                      M_neighbors[indices[j]].insert( indices[k] );
                      M_neighbors[indices[k]].insert( indices[j] );
                    }
        }

    // periodic tables if necessary
    if ( M_periodicity.isPeriodic() )
        createPeriodicCorrespondanceTable();

}


template<typename FunctionSpaceType, typename periodicity_type>
void
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
createPeriodicCorrespondanceTable()
{
    /* create a data structure periodic,
       a bi-directionnal map containing :
       idGlobal of dofs having tag1, corresponding by trans idGlobal of dofs having tag2
     */
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    // store idGlobTag1 and coords
    typedef std::vector< std::pair<size_type, node_type> > idCoord_type;
    idCoord_type idsTag1;
    std::set< size_type > dofTag1Done;
    auto rg1 = markedfaces( M_functionspace->mesh(), M_periodicity.tag1() );
    for ( auto it = rg1.template get<1>(), en = rg1.template get<2>(); it!=en; ++it)
        for (size_type k=0; k<M_functionspace->dof()->nLocalDofOnFace() ; ++k)
            {
                const size_type index = boost::get<0>( M_functionspace->dof()->localToGlobal( *it, k, 0 ) );
                const size_type idOnCluster = processorToCluster(index);
                const node_type coordPointPlusTrans = get<0>( M_functionspace->dof()->dofPoint( index ) ) + M_translation;

                if ( ! dofTag1Done.count(idOnCluster) && (! M_functionspace->dof()->dofGlobalProcessIsGhost(index)) )
                    {
                        idsTag1.push_back( std::make_pair(idOnCluster, coordPointPlusTrans) );
                        dofTag1Done.insert(idOnCluster);
                    }
            }

    dofTag1Done.clear(); // no need anymore

    M_nbDofTag1 = mpi::all_reduce(Environment::worldComm().globalComm(),
                                          idsTag1.size(),
                                          std::plus<size_type>() );

    std::vector< idCoord_type > allIdsTag1;
    mpi::all_gather( Environment::worldComm().globalComm(),
                     idsTag1,
                     allIdsTag1 );

    typedef std::vector< std::pair< size_type, size_type > > idstag1tag2_type;
    idstag1tag2_type idsTag1_idsTag2;

    auto distNodes = []( const node_type& a, const node_type& b )->double
        {
            const node_type c = a-b;
            double dist=0;
            for (int d=0; d<Dim; ++d)
                dist += std::abs( c[d] );
            return dist;
        };

    // for each point on the periodic face
    // search if there is a correspondant point on
    // the other face
    std::set< size_type > dofTag2Done;

    auto rg2 = markedfaces( M_functionspace->mesh(), M_periodicity.tag2() );
    for (auto it = rg2.template get<1>(), en=rg2.template get<2>(); it!=en; ++it)
        for (size_type k=0; k<M_functionspace->dof()->nLocalDofOnFace() ; ++k)
            {
                const size_type indexTag2 = boost::get<0>( M_functionspace->dof()->localToGlobal( *it, k, 0 ) );
                if ( ! dofTag2Done.count( indexTag2 ) )
                    {
                        const size_type globalIndexTag2 = processorToCluster( indexTag2 );
                        const node_type coordPointTag2 = get<0>( M_functionspace->dof()->dofPoint( indexTag2 ) );
                        bool dofFound = false;

                        for (auto const& allIdsVec : allIdsTag1)
                            {
                                for (auto const& idCoordPointTag1 : allIdsVec )
                                    {
                                        const double distPoints = distNodes( coordPointTag2, idCoordPointTag1.second );

                                        if (distPoints < 1e-10)
                                            {
                                                idsTag1_idsTag2.push_back( std::make_pair( idCoordPointTag1.first,
                                                                                           globalIndexTag2 ) );
                                                dofFound = true;
                                                dofTag2Done.insert(indexTag2);
                                                break;
                                            }
                                    }

                                if (dofFound)
                                    break;
                            }
                    }
            }

    dofTag2Done.clear();

    // only one proc gather all the informations in a bi-map
    std::vector< idstag1tag2_type > all_idsTag1_idsTag2;
    mpi::gather( Environment::worldComm().globalComm(),
                 idsTag1_idsTag2,
                 all_idsTag1_idsTag2,
                 0);

    if (Environment::worldComm().rank() == 0)
        for (auto const& idTag1Tag2 : all_idsTag1_idsTag2)
            for (auto const& id1_id2 : idTag1Tag2 )
                M_idTag1_idTag2.left.insert( id1_id2 );

    mpi::broadcast( Environment::worldComm(),
                    M_idTag1_idTag2,
                    0 );

    CHECK( M_idTag1_idTag2.left.size() == M_nbDofTag1 ) <<"problem in periodicity table of fast marching\n"
                                                              <<"all nodes tagged 1 should have a conter part in tag2\n"
                                                              <<"M_idTag1_idTag2.left.size() = "<< M_idTag1_idTag2.left.size() <<"\n"
                                                              <<"nbDofTag1 = "<< M_nbDofTag1 <<"\n";

}// createPeriodicCorrespondanceTable


template<typename FunctionSpaceType, typename periodicity_type>
void
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
reduceDonePoints(element_type const& __v, element_type& status, std::set<size_type>& done )
{
  /*  Communicate the DONE points across all the processes  */

  if (Environment::worldComm().size() == 1)
    {
      nbTotalDone = done.size();
      return ;
    }

  /* make the sum of the DONE(=2) + FAR(=0)
     the sum is then communicated through all the processors ( close() ) */
  checkStatus->zero();

  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
      checkStatus->add(k, status(k) == DONE ? 1 : 0);

  // communications are done here. The sum of the real values and ghost values is done
  // thus if a dof has been set to DONE, even on a ghost value, checkStatus is 1
  checkStatus->close();

  int nbDoneWithoutGhost=0;
  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
    if ( (*checkStatus)(k) )
      {
        done.insert( k );
        status.set(k, DONE);
        if (! M_functionspace->dof()->dofGlobalProcessIsGhost( k ) )
          ++nbDoneWithoutGhost;
      }

  nbTotalDone = mpi::all_reduce(Environment::worldComm().globalComm(),
                                nbDoneWithoutGhost,
                                std::plus<size_type>() );
}


template<typename FunctionSpaceType, typename periodicity_type>
void
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
reduceClosePoints(heap_type& theHeap, element_type& status )
{
  if (Environment::worldComm().size() == 1)
    return ;

  //need to communicate the list of all the CLOSE elements and their values to all processors
  checkStatus->zero();
  valueAtClose->zero();

  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
    if ( status(k) == CLOSE )
      {
        checkStatus->add(k, 1);
        valueAtClose->add(k, theHeap.valueAtIndex( k ) );
      }

  checkStatus->close();
  valueAtClose->close();

  // store global id of every CLOSE value and its phi associated value
  std::map< size_type, value_type > phiValueAtGlobalId;

  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
    {

      // if a dof has been marked as CLOSE on an other processor, mark it as CLOSE and take its value
      if ( ((*checkStatus)(k) == 1) && (status(k) != CLOSE) )
        {
          status(k) = CLOSE;
          heap_entry_type newEntry = { (*valueAtClose)(k), k };
          theHeap.push( newEntry );
        }

      // if a point is marked as CLOSE from at least two different processors, store its value of phi in order to be compared with the others
      else if ( ( (*checkStatus)(k) > 1 ) && (status(k) == CLOSE) )
        phiValueAtGlobalId[ processorToCluster( k ) ] = theHeap.valueAtIndex( k );
    }


  std::vector< std::map< size_type, value_type > > all_phiValueAtGlobalId;

  mpi::all_gather( Environment::worldComm().globalComm(),
                   phiValueAtGlobalId,
                   all_phiValueAtGlobalId );

  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
      if ( (*checkStatus)(k) > 1 )
        {
          const size_type idOnCluster = processorToCluster( k );

          for ( auto const & phiAtGlobId : all_phiValueAtGlobalId )
              if ( phiAtGlobId.count( idOnCluster ) )
                {
                  heap_entry_type newEntry = { phiAtGlobId.at(idOnCluster) , k};
                  if ( theHeap.checkExistingEntry( k ) )
                    theHeap.change( newEntry ); //change only if phiAtGlobId[idOnCluster] < phi already stored
                  else
                    theHeap.push( newEntry );
                }

          status(k) = CLOSE;

        }
}



template<typename FunctionSpaceType, typename periodicity_type>
typename ReinitializerFMS<FunctionSpaceType, periodicity_type>::element_type
ReinitializerFMS<FunctionSpaceType, periodicity_type>::operator()
  ( element_type const& phi, bool useMarker2AsDoneMarker ) /*const*/
{

#if defined( FM_EXPORT )
    auto ex = exporter(_mesh=M_functionspace->mesh(), _name="fastmarchin");
#endif

    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    auto __v = vf::project(M_functionspace, elements(M_functionspace->mesh()), idv(phi) );

    // acquire interface (=done) cells
    // the chosen approach assumes monotonicity of the FE-function in the
    // element, thus valid only for P1 elements

    std::set<size_type> done;
    auto status = vf::project( M_functionspace, elements(M_functionspace->mesh()), vf::cst(FAR) );

    /* VD, sometime, I need to give myself the elements which are already done, thus use marker2*/
    if (useMarker2AsDoneMarker)
      {
        auto it_marked = M_functionspace->mesh()->elementsWithMarker2(1, M_functionspace->mesh()->worldComm().localRank()).first;
        auto en_marked = M_functionspace->mesh()->elementsWithMarker2(1, M_functionspace->mesh()->worldComm().localRank()).second;

        for ( ; it_marked!=en_marked ; ++it_marked )
          for ( uint16_type j = 0; j < ndofv; ++j )
            {
              size_type index = M_functionspace->dof()->localToGlobal(*it_marked, j, 0).index();
              done.insert( index );
              status[index] = DONE;
            }
      }

    else
        {
          auto it = M_functionspace->mesh()->beginElementWithProcessId();
          auto en = M_functionspace->mesh()->endElementWithProcessId();

            for ( ; it!=en ; ++it )
                {
                    uint16_type nPlus = 0;
                    uint16_type nMinus = 0;
                    std::vector<size_type> indices( ndofv );
                    for ( uint16_type j = 0; j < ndofv; ++j )
                        {
                            size_type index = M_functionspace->dof()->localToGlobal(*it, j, 0).index();
                            indices[j] = index;
                            if ( phi[index] < 0.0 )
                                ++nMinus;
                            else if ( phi[index] > 0.0 )
                                ++nPlus;
                        }
                    //if the element has nodes with positive and negative values
                    //mark as done
                    if ( nPlus != ndofv && nMinus != ndofv )
                            for ( uint16_type j = 0; j < ndofv; ++j )
                                {
                                    done.insert( indices[j] );
                                    status[indices[j]] = DONE;
                                }
                }
        }

    // //communicate the DONE list between all the proc
    reduceDonePoints(__v, status, done );

    heap_type theHeap;
    // initialize close distances in heap and mark close points in status array
    for ( auto dit = done.begin(); dit != done.end(); ++dit )
      fmsHeapUpdate( *dit, __v, status, theHeap );

    done.clear(); // not needed any more, save memory...

    // if there is no interface, do nothing
    if (nbTotalDone == 0)
        return __v;

    // for debug purpuses only
    auto checkHeap = [&] (std::string msg)
        {
#if 0
        for (auto entry : theHeap )
          CHECK( status( entry.second ) == CLOSE )
            << std::endl << msg << std::endl
            << "on proc " << Environment::worldComm().rank()
            << " the entry at ldof "<< entry.second
            << " and entry at cdof "<< processorToCluster( entry.second )
            << " is marked as " << status( entry.first ) <<std::endl ;
#endif
        };


    checkHeap( "before reduce close points" );
    reduceClosePoints( theHeap, status );
    checkHeap( "after reduce close points" );

    auto minNewEntry = [](std::pair<heap_entry_type, size_type> a,
                          std::pair<heap_entry_type, size_type> b)
      { // return the entry having the minimum abs(phi) value
        return std::abs(a.first.first) < std::abs(b.first.first) ? a : b; };


    const int nbTotalIterFM = M_functionspace->dof()->nDof() - nbTotalDone - M_nbDofTag1;

    unsigned int count_iteration=0;

    for (int i=0; i<nbTotalIterFM; ++i)
      {
        /*VD : the heap is sorted with the smallest phi value
          at the top.
          thus, the new accepted element is always the top of the heap */

        // contains entry type and global index on Cluster
        std::pair<heap_entry_type, size_type> newAccepted = { theHeap.front(), processorToCluster(theHeap.front().second) };

        // the real new accepted value is the min of all the phi computed in the heaps
        newAccepted = mpi::all_reduce(Environment::worldComm().globalComm(),
                                      newAccepted,
                                      minNewEntry);

        size_type newIdOnCluster = newAccepted.second;

        bool dofIsPresentOnProcess = false;

        //        size_type newIdOnProc;
        std::vector< size_type > newIdsOnProc;
        if ( M_functionspace->dof()->dofGlobalClusterIsOnProc( newIdOnCluster ) )
              newIdsOnProc.push_back( clusterToProcessor(newIdOnCluster) );

        else if ( M_ghostClusterToProc.count( newIdOnCluster ) )
            newIdsOnProc.push_back( M_ghostClusterToProc[newIdOnCluster] );

        if (M_periodicity.isPeriodic())
            {
                size_type idOnClusterPeriodicPoint = invalid_size_type_value;
                // look from tag1 to tag 2
                if (M_idTag1_idTag2.left.count( newIdOnCluster ) )
                    idOnClusterPeriodicPoint = M_idTag1_idTag2.left.at( newIdOnCluster );
                else if (M_idTag1_idTag2.right.count( newIdOnCluster ) )
                    idOnClusterPeriodicPoint = M_idTag1_idTag2.right.at( newIdOnCluster );

                if ( idOnClusterPeriodicPoint != invalid_size_type_value )
                    {
                        if ( M_functionspace->dof()->dofGlobalClusterIsOnProc( idOnClusterPeriodicPoint ) )
                            newIdsOnProc.push_back( clusterToProcessor(idOnClusterPeriodicPoint) );

                        else if ( M_ghostClusterToProc.count( idOnClusterPeriodicPoint ) )
                            newIdsOnProc.push_back( M_ghostClusterToProc[idOnClusterPeriodicPoint] );
                    }
            }

        for( size_type newIdOnProc : newIdsOnProc )
          {
            theHeap.removeFromHeap( newIdOnProc );

            status[newIdOnProc] = DONE;

            __v[newIdOnProc] = newAccepted.first.first;

            fmsHeapUpdate( newIdOnProc, __v, status, theHeap );
          }

        checkHeap("in loop");

#if defined( FM_EXPORT )
        count_iteration++;
        ex->step(count_iteration)->add("v", __v);
        ex->step(count_iteration)->add("status", status);
        ex->save();
#endif

      } // marching loop

    return __v;

} // operator()


template<typename FunctionSpaceType, typename periodicity_type>
void ReinitializerFMS<FunctionSpaceType, periodicity_type>::
fmsHeapUpdate( size_type idDone,
               element_type const& __v,
               element_type& status,
               heap_type& theHeap ) const
{
    // nbr : neighbours of the node just done
    std::set<size_type> const & nbrs = M_neighbors.find(idDone)->second;

    std::vector<size_type> ids( 1, idDone );
    for ( auto n0it = nbrs.begin(); n0it != nbrs.end(); ++n0it )
        {
            if ( status[*n0it] == FAR )
                status[*n0it] = CLOSE;
            if ( status[*n0it] == CLOSE )
                {

                  /* to give a reference, compute phi with only one DONE neighbours
                    it is sure that *n0it is a DONE neighbours */
                    // one neighbor
                    ids.push_back(*n0it);
                    value_type phiNew = fmsDistN( ids, __v );
                    ids.pop_back();

                    /*compute all the phi possible with all the neighbors around and returns the smallest one*/
                    phiNew = fmsDistRec( ids, *n0it, __v, phiNew, status );

                    theHeap.change( std::make_pair( phiNew, *n0it ) );
                } // if CLOSE
        } // loop over neighbor 0
} // fmsHeapUpdate


template<typename FunctionSpaceType, typename periodicity_type>
typename ReinitializerFMS<FunctionSpaceType, periodicity_type>::value_type
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
fmsDistRec( std::vector<size_type> & ids,
            size_type idClose,
            element_type const& __v,
            value_type phiOld,
            element_type const& status ) const
{
    /*search for all DONE neighbors
      recalculate phi with all the neighbors
      returns the smallest phi  */

    /* the method fmsDistN computes the distance function for a given CLOSE
       using the DONE neighbors in the same element
       thus, first it has to search for the ids of the neighbors DONE
       being in the same element */

    // only allows for getting at maximum 2 nodes in 2d and 3 nodes in 3d
    if ( ids.size() >= Dim )
        return phiOld;

    value_type phiNew(phiOld);

    std::set<size_type> const & nbrs = M_neighbors.find(idClose)->second;
    for (auto nit = nbrs.begin(); nit != nbrs.end(); ++nit )
        {

          // if the next neighbors is not a DONE neighbors
            if ( status[*nit] != DONE )
                continue;

            // if this node has already been accepted, stop
            bool unique = true;
            for ( auto idsit=ids.begin(); idsit != ids.end(); ++idsit )
                unique &= ( *idsit != *nit );
            if ( !unique ) // points must be unique
                continue;

            if ( M_neighbors.find(ids[0])->second.find(*nit) ==
                 M_neighbors.find(ids[0])->second.end() )
                continue;

            // one neighbor more
            ids.push_back(*nit);

            ids.push_back(idClose);
            value_type phiCand = fmsDistN( ids, __v );
            ids.pop_back();
            if ( phiCand != 0.0 )
                phiNew = closerOne( phiCand, phiNew );

            phiNew = fmsDistRec( ids, idClose, __v, phiNew, status );

            ids.pop_back();
        }
    return phiNew;
} // fmsDistRec


template<typename FunctionSpaceType, typename periodicity_type>
typename ReinitializerFMS<FunctionSpaceType, periodicity_type>::value_type
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
fmsDistN( std::vector<size_type> const & ids,
          element_type const & __v ) const
{
    /*
      VD :
      return value of phi calculated from the KNOWN points having id contained in ids
      first entry if ids is the index of the point to calculate
      values of phi at ids are stored in __v (calculated previously)
    */

    uint32_type nPts = ids.size()-1; // number of KNOWN points

    std::vector<point_type> basis( nPts );
    std::vector<value_type> wNorm( nPts );
    point_type grad;
    value_type n_rest = 1.0;
    std::vector<std::vector<value_type> >
        q( nPts, std::vector<value_type>( nPts, 0.0 ) );
    std::vector<value_type> n( nPts, 0.0 );
    value_type eps = type_traits<value_type>::epsilon();
    for ( uint32_type i=0; i<nPts; ++i )
        {
            basis[i] = M_coords[ids[i+1]] - M_coords[ids[0]];
            for ( uint32_type j=0; j<i; j++ )
                {
                    q[i][j] = dot( basis[i], basis[j] );
                    basis[i] -= q[i][j]*basis[j];
                }
            wNorm[i] = 1.0/norm(basis[i]);
            basis[i] *= wNorm[i];
            if ( i<nPts-1 )
                {
                    n[i] = __v[ids[i+1]]-__v[ids[0]];
                    for ( uint32_type k=0; k<i; ++k )
                        n[i] -= ( __v[ids[k+1]] - __v[ids[0]] ) * q[i][k] * wNorm[k];
                    n[i] *= wNorm[i];
                    n_rest -= n[i]*n[i];
                }
            else
                {
                    if ( n_rest < -10.0*eps )
                        return 0.0;
                    else if ( n_rest < 0.0 )
                        n_rest = 0.0;
                    n[i] = std::sqrt( n_rest );
                    if ( __v[ids[0]] < 0.0 )
                        n[i] *= -1.0;
                }
            grad += n[i]*basis[i];
        }

    // verify gradient to new value comes through convex hull of points used
    // to construct it
    point_type dx( M_coords[ids[nPts]] - M_coords[ids[0]] );
    std::vector<value_type> lambda( nPts, 0.0 );
    lambda[nPts-1] = dot( dx , basis[nPts-1] ) / dot( grad, basis[nPts-1] );
    value_type lambdaTot = 0.0;
    bool inside = true;
    uint32_type j=nPts-1;
    while ( j > 0 )
        {
            --j;
            lambda[j] = dot( dx, basis[j] ) - lambda[nPts-1]*n[j];
            for ( uint32_type i=j+1; i<nPts-1; ++i )
                lambda[j] -= lambda[i] * q[i][j];
            lambda[j] *= wNorm[j];
            lambdaTot += lambda[j];
            inside &= lambda[j] >= 0.0;
            // inside &= lambda[j] <= 1.0 // for line products
        }
    inside &= lambdaTot <= 1.0; // simplex assumed!

    return inside ? __v[ids[0]] + dot( grad, dx ) : 0.0;
} // fmsDistN



} // Feel
#endif
