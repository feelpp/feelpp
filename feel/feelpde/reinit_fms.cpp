#ifndef REINIT_FMS_CPP
#define REINIT_FMS_CPP 1

#define FEELPP_INSTANTIATE_FMS 1

#include "reinit_fms.hpp"

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
    ndofOnCluster( __functionspace->dof()->nDof() )
{
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    auto it = M_functionspace->mesh()->beginElementWithProcessId();
    auto en = M_functionspace->mesh()->endElementWithProcessId();

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

}

template<typename FunctionSpaceType, typename periodicity_type>
void
ReinitializerFMS<FunctionSpaceType, periodicity_type>::
reduceDonePoints(element_type const& __v, heap_type& theHeap, element_type& status, std::set<size_type>& done )
{
  /*  Communicate the DONE points across all the processes  */

  if (Environment::worldComm().size() == 1)
    return ;

  /* make the sum of the DONE(=2) + FAR(=0)
     the sum is then communicated through all the processors ( close() ) */
  checkStatus->zero();

  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
    checkStatus->add(k, status(k) == DONE ? 1 : 0);

  // communications are done here. The sum of the real values and ghost values is done
  // thus if a dof has been set to DONE, even on a ghost value, checkStatus is 1
  checkStatus->close();

  for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
    if ( (*checkStatus)(k) )
      {
        done.insert( k );
        status.set(k, DONE);
      }

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

    heap_type theHeap;

    // //communicate the DONE list between all the proc
    if ( ! useMarker2AsDoneMarker )
      reduceDonePoints(__v, theHeap, status, done );

    // initialize close distances in heap and mark close points in status array
    for ( auto dit = done.begin(); dit != done.end(); ++dit )
      fmsHeapUpdate( *dit, __v, status, theHeap );

    done.clear(); // not needed any more, save memory...

    reduceClosePoints( theHeap, status );

#if 0
    // for debug purpuses only
    auto checkHeap = [&] ()
        {
        for (auto entry : theHeap )
          CHECK( status( entry.second ) == CLOSE )
            << "on proc " << Environment::worldComm().rank()
            << " the entry at ldof "<< entry.second
            << " and entry at cdof "<< processorToCluster( entry.second )
            << " is marked as " << status( entry.first ) <<std::endl ;
      };
#endif

    size_type sumAllSizes= mpi::all_reduce(Environment::worldComm().globalComm(),
                                           theHeap.size(),
                                           std::plus<size_type>() );

    // marching loop
    // continue since all the Heap are not at 0
    int count_iteration = 0;

    while ( sumAllSizes )
      {

        /*VD : the heap is sorted with the smallest phi value
          at the top.
          thus, the new accepted element is always the top of the heap */

        auto newAccepted = theHeap.front();

        // the real new accepted value is the min of all the phi computed in the heaps
        newAccepted = mpi::all_reduce(Environment::worldComm().globalComm(),
                                      newAccepted,
                                      theHeap.min);

        size_type newIdOnCluster = 0;
        // extract the global on cluster index of the point accepted
        if ( theHeap.front() == newAccepted )
          newIdOnCluster = processorToCluster(newAccepted.second);


        // propagate the global id on cluster of the new accepted to all the proc
        newIdOnCluster = mpi::all_reduce(Environment::worldComm().globalComm(),
                                         newIdOnCluster,
                                         mpi::maximum<size_type>() );

        bool dofIsPresentOnProcess = false;

        size_type newIdOnProc;
        if ( M_functionspace->dof()->dofGlobalClusterIsOnProc( newIdOnCluster ) )
          {
            newIdOnProc = clusterToProcessor(newIdOnCluster);
            dofIsPresentOnProcess=true;
          }
        else if ( M_ghostClusterToProc.count( newIdOnCluster ) )
          {
            newIdOnProc = M_ghostClusterToProc[newIdOnCluster];
            dofIsPresentOnProcess=true;
          }

        if (dofIsPresentOnProcess)
          {

            theHeap.removeFromHeap( newIdOnProc );

            status[newIdOnProc] = DONE;

            __v[newIdOnProc] = newAccepted.first;

            fmsHeapUpdate( newIdOnProc, __v, status, theHeap );

            if (M_periodicity.isPeriodic())
              updatePeriodicPoint(newAccepted, __v, status, theHeap);

          }

        sumAllSizes = mpi::all_reduce(Environment::worldComm().globalComm(),
                                      theHeap.size(),
                                      std::plus<size_type>() );

        // avoid infinite loop if a heap is not cleaned correctly for any reason
        // CHECK(++count_iteration < ndofOnCluster)<<"something is wrong in fastmarching loop. The march should converge in less than ndofOnCluster = "
        //                                         <<ndofOnCluster<<" iterations and still has not converged in "<<count_iteration<<" iterations"<<std::endl;

#if defined( FM_EXPORT )
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
void ReinitializerFMS<FunctionSpaceType, periodicity_type>::
updatePeriodicPoint(heap_entry_type const& newAccepted, element_type& __v, element_type& status, heap_type& theHeap) const
{

  /* ++++++++++++++++  This function has not been ported to parallel yet +++++++++++++ */

    // search if the point at newIdx has a corresponding point by periodic translation
    // if yes, mark the point as done and update the heap in consequence

    vf::node_type pt_tofindPlus(Dim);
    vf::node_type pt_tofindMinus(Dim);
    vf::node_type pt_found(Dim);

    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    uint16_type newIdx = newAccepted.second;

    for (int c = 0 ; c<Dim; c++)
        {
            pt_tofindPlus [c] = M_coords[newIdx][c] + M_translation[c];
            pt_tofindMinus[c] = M_coords[newIdx][c] - M_translation[c];
        }

    int foundPlusMinus = 0;
    size_type id_found;

    if (M_functionspace->findPoint(pt_tofindPlus, id_found, pt_found))
        foundPlusMinus = 1;
    else if (M_functionspace->findPoint(pt_tofindMinus, id_found, pt_found))
        foundPlusMinus = -1;

    if (foundPlusMinus)
      for (int k = 0; k < ndofv; ++k)
        {
          size_type index_k = M_functionspace->dof()->localToGlobal(id_found, k, 0).index();
          value_type diff_position=0;
          for (int c = 0; c < Dim; c++)
            diff_position += std::abs(M_coords[index_k][c] - (M_coords[newIdx][c] + foundPlusMinus * M_translation[c] ) );
          if ( diff_position < 1e-9 )
            {
              //set the value of phi at the point
              __v[index_k] = newAccepted.first;

              //remove it from the heap
              theHeap.removeFromHeap(index_k);

              // mark as DONE
              status[index_k]=DONE;

              //add its neigbours in CLOSE
              fmsHeapUpdate(index_k, __v, status, theHeap);
              break ; // don't need to search for other corresp
            }
        }

} //updatePeriodicDof


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


template class ReinitializerFMS< spaceP1LS_type, Periodic<> > ;
template class ReinitializerFMS< spaceP1LS_type, NoPeriodicity  > ;
template class ReinitializerFMS< spaceP1LS_3d_type, NoPeriodicity > ;
template class ReinitializerFMS< spaceP1LS_3d_type, Periodic<> > ;

} // Feel
#endif
