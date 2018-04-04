#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operators.hpp>

#define FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS \
    template< typename FunctionSpaceType, typename PeriodicityType, typename HeapDataType > \
    /**/ 
#define FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE \
    FastMarchingBase<FunctionSpaceType, PeriodicityType, HeapDataType> \
    /**/ 

namespace Feel {

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::FastMarchingBase( 
        functionspace_ptrtype const& space,
        periodicity_type periodicity )
: 
    M_functionspace( space ),
    M_status( M_functionspace->elementPtr() ),
    M_checkStatus(backend()->newVector(M_functionspace)),
    M_valueAtClose(backend()->newVector(M_functionspace)),
    M_neighbors(),
    M_distance( M_functionspace->elementPtr() ),
    M_periodicity(periodicity),
    M_coords( M_functionspace->dof()->nDof() ),
    M_translation( periodicity.translation() ),
    M_firstDof( M_functionspace->dof()->firstDofGlobalCluster() ),
    M_nbDofTag1(0)
{
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    // create the first neighbours data structure (M_neighbors)
    for ( auto const& eltWrap : elements( M_functionspace->mesh() ) )
    {
        auto const& elt = boost::unwrap_ref( eltWrap );
        std::vector<size_type> indices( ndofv );
        for ( uint16_type j = 0; j < ndofv;++j )
        {
            size_type index = M_functionspace->dof()->localToGlobal(elt, j, 0).index();

            if (M_functionspace->dof()->dofGlobalProcessIsGhost( index ))
            {
                M_ghostClusterToProc[ processorToCluster( index ) ] = index;
            }

            indices[j] = index;
            for ( uint16_type c = 0; c < Dim; ++c )
            {
                M_coords[index][c] = M_functionspace->dof()->dofPoint( index ).template get<0>()[c];
            }
        }
        for ( uint16_type j = 0; j < ndofv;++j )
        {
            for ( uint16_type k = j+1; k < ndofv;++k )
            {
                M_neighbors[indices[j]].insert( indices[k] );
                M_neighbors[indices[k]].insert( indices[j] );
            }
        }
    }

    // periodic tables if necessary
    if ( M_periodicity.isPeriodic() )
    {
        createPeriodicCorrespondanceTable();
    }
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
void
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::run( 
        element_type const& phi, 
        bool useMarker2AsMarkerDone )
{
    FEELPP_ASSERT( FunctionSpaceType::fe_type::nOrder == 1 )("FunctionSpaceType needs to be a finite element space of order 1");

    this->initMarch( phi, useMarker2AsMarkerDone );

    // if there is no interface, do nothing
    if (M_nbTotalDone == 0)
        return;

    //// for debug purposes only
    //auto checkHeap = [&] (std::string msg)
        //{
//#if 0
        //for (auto entry : M_heap )
          //CHECK( (*M_status)( entry.second ) == CLOSE )
            //<< std::endl << msg << std::endl
            //<< "on proc " << Environment::worldComm().rank()
            //<< " the entry at ldof "<< entry.second
            //<< " and entry at cdof "<< processorToCluster( entry.second )
            //<< " is marked as " << (*M_status)( entry.first ) <<std::endl ;
//#endif
        //};

    //checkHeap( "before reduce close points" );
    reduceClosePoints();
    //checkHeap( "after reduce close points" );

    auto minNewEntry = [](
            std::pair<heap_entry_data_type, size_type> a,
            std::pair<heap_entry_data_type, size_type> b )
    { // return the entry having the minimum abs(phi) value
        return std::abs(a.first.value()) < std::abs(b.first.value()) ? a : b; 
    };

    const int nbTotalIterFM = M_functionspace->dof()->nDof() - M_nbTotalDone - M_nbDofTag1;
    //for (int i=0; i<nbTotalIterFM; ++i)

    bool allHeapAreEmpty = false;
    //while( M_heap.size() != 0 )
    while( !allHeapAreEmpty )
    {
        /* The heap is sorted with the smallest phi value at the top.
         * Thus, the new accepted element is always the top of the heap 
         */

        // contains entry type and global index on Cluster
        auto heap_front = M_heap.front();
        std::pair<heap_entry_data_type, size_type> newAccepted = 
        { heap_front, processorToCluster(heap_front.index()) };

        // the real new accepted value is the min of all the phi computed in the heaps
        newAccepted = mpi::all_reduce(Environment::worldComm().globalComm(),
                                      newAccepted,
                                      minNewEntry);

        size_type newIdOnCluster = newAccepted.second;

        bool dofIsPresentOnProcess = false;

        std::vector< size_type > newIdsOnProc;
        if ( M_functionspace->dof()->dofGlobalClusterIsOnProc( newIdOnCluster ) )
              newIdsOnProc.push_back( clusterToProcessor(newIdOnCluster) );

        else if ( M_ghostClusterToProc.count( newIdOnCluster ) )
            newIdsOnProc.push_back( M_ghostClusterToProc[newIdOnCluster] );

        for( size_type newIdOnProc : newIdsOnProc )
        {
            M_heap.removeFromHeap( newIdOnProc );

            this->processDof( newIdOnProc, newAccepted.first.value(), newAccepted.first.data() );

            this->updateHeap( newIdOnProc );
        }

        const bool heapIsEmpty = M_heap.size() == 0;
        allHeapAreEmpty = mpi::all_reduce(M_functionspace->mesh()->worldComm(),
                                          heapIsEmpty,
                                          std::logical_and<bool>() );
    }

    M_heap.clear();
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::clusterToProcessor( size_type dof ) const
{
    return dof - M_firstDof; 
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::processorToCluster( size_type dof ) const
{
    return M_functionspace->dof()->mapGlobalProcessToGlobalCluster( dof ); 
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
void 
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::createPeriodicCorrespondanceTable()
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
            const size_type index = boost::get<0>( M_functionspace->dof()->localToGlobal( boost::unwrap_ref( *it ), k, 0 ) );
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
            const size_type indexTag2 = boost::get<0>( M_functionspace->dof()->localToGlobal( boost::unwrap_ref( *it ), k, 0 ) );
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

    CHECK( M_idTag1_idTag2.left.size() == M_nbDofTag1 ) 
        <<"problem in periodicity table of fast marching\n"
        <<"all nodes tagged 1 should have a conter part in tag2\n"
        <<"M_idTag1_idTag2.left.size() = "<< M_idTag1_idTag2.left.size() <<"\n"
        <<"nbDofTag1 = "<< M_nbDofTag1 <<"\n";

}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
void 
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::reduceDonePoints(
        std::set<size_type>& doneIds )
{
  /*  Broadcast the DONE points across all the processes  */
    if (Environment::worldComm().size() == 1)
    {
        M_nbTotalDone = doneIds.size();
        return ;
    }

    /* make the sum of the DONE(=2) + FAR(=0)
       the sum is then communicated through all the processors ( close() ) */
    M_checkStatus->zero();

    for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
        M_checkStatus->add(k, (*M_status)(k) == DONE ? 1 : 0);

    // communications are done here. The sum of the real values and ghost values is done
    // thus if a dof has been set to DONE, even on a ghost value, M_checkStatus is 1
    M_checkStatus->close();

    int nbDoneWithoutGhost=0;
    for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
        if ( (*M_checkStatus)(k) )
        {
            doneIds.insert( k );
            (*M_status).set(k, DONE);
            if (! M_functionspace->dof()->dofGlobalProcessIsGhost( k ) )
                ++nbDoneWithoutGhost;
        }

    M_nbTotalDone = mpi::all_reduce(Environment::worldComm().globalComm(),
            nbDoneWithoutGhost,
            std::plus<size_type>() );
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
void 
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::reduceClosePoints()
{
    if (Environment::worldComm().size() == 1)
        return ;

    //need to communicate the list of all the CLOSE elements and their values to all processors
    M_checkStatus->zero();
    //M_valueAtClose->zero();

    for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
        if ( (*M_status)(k) == CLOSE )
        {
            M_checkStatus->add(k, 1);
            //M_valueAtClose->add(k, M_heap.valueAtIndex( k ) );
        }

    M_checkStatus->close();
    //M_valueAtClose->close();

    // store global id of every CLOSE value and its phi associated value
    std::map< size_type, value_data_type > phiValueAtGlobalId;

    for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
    {
        // if a dof has been marked as CLOSE on an other processor, mark it as CLOSE and take its value
        //if ( ((*M_checkStatus)(k) == 1) && ((*M_status)(k) != CLOSE) )
        //{
            //(*M_status)(k) = CLOSE;
            //std::vector<value_type> datak = (this->hasOptionalData())? (*M_dataAtClose)(k): std::vector<value_type>();
            //heap_entry_data_type newEntry = { (*M_valueAtClose)(k), k, datak };
            //M_heap.push( newEntry );
        //}

        // if a point is marked as CLOSE from at least two different processors, store its value of phi in order to be compared with the others
        //else if ( ( (*M_checkStatus)(k) > 1 ) && ((*M_status)(k) == CLOSE) )
        if ( ( (*M_checkStatus)(k) > 0 ) && ((*M_status)(k) == CLOSE) )
        {
            phiValueAtGlobalId[ processorToCluster( k ) ] = value_data_type(
                    M_heap.valueAtIndex(k), M_heap.dataAtIndex(k) );
        }
    }


    std::vector< std::map< size_type, value_data_type > > all_phiValueAtGlobalId;

    mpi::all_gather( Environment::worldComm().globalComm(),
            phiValueAtGlobalId,
            all_phiValueAtGlobalId );

    for (size_type k = 0 ; k < M_functionspace->nLocalDof() ; ++k)
        //if ( (*M_checkStatus)(k) > 1 )
        if ( (*M_checkStatus)(k) > 0 )
        {
            const size_type idOnCluster = processorToCluster( k );

            for ( auto const & phiAtGlobId : all_phiValueAtGlobalId )
                if ( phiAtGlobId.count( idOnCluster ) )
                {
                    heap_entry_type newEntry = { phiAtGlobId.at(idOnCluster).value() , k };
                    M_heap.change( newEntry, phiAtGlobId.at(idOnCluster).data()  ); //change only if phiAtGlobId[idOnCluster] < phi already stored
                }

            (*M_status)(k) = CLOSE;
        }
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
void
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::initMarch( 
        element_type const& phi, 
        bool useMarker2AsMarkerDone )
{
    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    *M_distance = vf::project(M_functionspace, elements(M_functionspace->mesh()), idv(phi) );

    // acquire interface (=doneIds) cells
    // the chosen approach assumes monotonicity of the FE-function in the
    // element, thus valid only for P1 elements

    std::set<size_type> doneIds;
    (*M_status) = vf::project( M_functionspace, elements(M_functionspace->mesh()), vf::cst(FAR) );

    /* VD, sometime, I need to give myself the elements which are already done, thus use marker2*/
    if (useMarker2AsMarkerDone)
    {
        auto rangeMarked2Elements = M_functionspace->mesh()->elementsWithMarker2(1, M_functionspace->mesh()->worldComm().localRank());
        auto it_marked = std::get<0>( rangeMarked2Elements );
        auto en_marked = std::get<1>( rangeMarked2Elements );

        for ( ; it_marked!=en_marked ; ++it_marked )
        {
            auto const& elt = boost::unwrap_ref( *it_marked );
            for ( uint16_type j = 0; j < ndofv; ++j )
            {
                size_type index = M_functionspace->dof()->localToGlobal( elt, j, 0).index();
                doneIds.insert( index );
                (*M_status)[index] = DONE;
            }
        }
    }

    else
    {
        auto rangeElements = M_functionspace->mesh()->elementsWithProcessId();
        auto it = std::get<0>( rangeElements );
        auto en = std::get<1>( rangeElements );

        for ( ; it!=en ; ++it )
        {
            auto const& elt = boost::unwrap_ref( *it );
            uint16_type nPlus = 0;
            uint16_type nMinus = 0;
            std::vector<size_type> indices( ndofv );
            for ( uint16_type j = 0; j < ndofv; ++j )
            {
                size_type index = M_functionspace->dof()->localToGlobal(elt, j, 0).index();
                indices[j] = index;
                if ( phi[index] < 0.0 )
                    ++nMinus;
                else if ( phi[index] > 0.0 )
                    ++nPlus;
            }
            //if the element has nodes with positive and negative values
            //mark as done
            if ( nPlus != ndofv && nMinus != ndofv )
            {
                for ( uint16_type j = 0; j < ndofv; ++j )
                {
                    doneIds.insert( indices[j] );
                    (*M_status)[indices[j]] = DONE;
                }
            }
        }
    }

    // communicate the DONE list between all the proc
    reduceDonePoints( doneIds );

    // initialize close distances in heap and mark close points in (*M_status) array
    for ( auto dit = doneIds.begin(); dit != doneIds.end(); ++dit )
        this->updateHeap( *dit );

    //doneIds.clear(); // not needed any more, save memory...
    //return doneIds;
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
void 
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::updateHeap(
        size_type idDone )
{
    // nbr : neighbours of the node just done
    std::set<size_type> const & nbrs = M_neighbors.find(idDone)->second;

    std::vector<size_type> ids( 1, idDone );
    for ( auto n0it = nbrs.begin(); n0it != nbrs.end(); ++n0it )
    {
        if ((*M_status)[*n0it] == FAR )
           (*M_status)[*n0it] = CLOSE;
        if ((*M_status)[*n0it] == CLOSE )
        {
            /* to give a reference, compute phi with only one DONE neighbours
               it is sure that *n0it is a DONE neighbours */
            // one neighbor
            ids.push_back(*n0it);
            value_type phiNew = fmsDistN( ids, *M_distance );
            ids.pop_back();

            /*compute all the phi possible with all the neighbors around and returns the smallest one*/
            phiNew = fmsDistRec( ids, *n0it, phiNew );

            M_heap.change( std::make_pair( phiNew, *n0it ) );
        } // if CLOSE
    } // loop over neighbor 0
}

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
typename FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::value_type
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::fmsDistN( 
        std::vector<size_type> const & ids,
        element_type const & __v
        ) const
{
    /*VD :
     * return value of phi calculated from the KNOWN points having id contained in ids
     * first entry if ids is the index of the point to calculate
     * values of phi at ids are stored in __v (calculated previously)
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

FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
typename FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::value_type
FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE::fmsDistRec( 
        std::vector<size_type> & ids,
        size_type idClose,
        value_type phiOld ) const
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
        if ((*M_status)[*nit] != DONE )
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
        value_type phiCand = fmsDistN( ids, *M_distance );
        ids.pop_back();
        if ( phiCand != 0.0 )
            phiNew = closerOne( phiCand, phiNew );

        phiNew = fmsDistRec( ids, idClose, phiNew );

        ids.pop_back();
    }
    return phiNew;
} // fmsDistRec

} // namespace Feel

#undef FASTMARCHINGBASE_CLASS_TEMPLATE_DECLARATIONS
#undef FASTMARCHINGBASE_CLASS_TEMPLATE_TYPE
