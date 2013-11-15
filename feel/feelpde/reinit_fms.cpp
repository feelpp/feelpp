#ifndef REINIT_FMS_CPP
#define REINIT_FMS_CPP 1

#define FEELPP_INSTANTIATE_FMS 1

#include <feel/feelpde/reinit_fms.hpp>

namespace Feel
{
template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type>
ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::
ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                  IteratorRange const& r ,
                  periodicity_type __periodicity)
    :
    M_functionspace( __functionspace ),
    M_range( r ),
    M_periodicity(__periodicity),
    M_neighbors(),
    M_coords( __functionspace->dof()->nDof() )
{
    using namespace Feel;

    LOG(INFO) << "ReinitializerFMS constructor from space and iterator range";

    fe_type* __fe = __functionspace->fe().get();

    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    typename gm_type::precompute_ptrtype
        __geopc( new typename gm_type::precompute_type(__functionspace->gm(),
                                                       __fe->points() ) );

    const uint16_type ndofv = functionspace_type::fe_type::nDof;
    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = M_range;

    gm_context_ptrtype __c( new gm_context_type( __functionspace->gm(),
                                                 *it,
                                                 __geopc ) );
    // acquire neighborship and node coordinates
    for ( ; it!=en ; ++it )
        {
            __c->update( *it, __geopc );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            t_expr_type tensor_expr( vf::P(), mapgmc );
            std::vector<size_type> indices( ndofv );
            for ( uint16_type __j = 0; __j < ndofv;++__j )
                {
                    size_type index = boost::get<0>(M_functionspace->dof()->localToGlobal( it->id(), __j, 0 ));
                    indices[__j] = index;
                    for ( uint16_type c = 0; c < Dim; ++c )
                        M_coords[index][c] = tensor_expr.evalq( c, 0, __j );
                }
            for ( uint16_type __j = 0; __j < ndofv;++__j )
                for ( uint16_type __k = __j+1; __k < ndofv;++__k )
                    {
                        M_neighbors[indices[__j]].insert( indices[__k] );
                        M_neighbors[indices[__k]].insert( indices[__j] );
                    }
        }

    // get the translation
    M_translation = __periodicity.translation();

}


template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type>
typename ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::element_type
ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::operator()
    ( element_type const& phi ) const
{
    using namespace Feel;

    //     Debug() << "[ReinitFMS] operator()\n";

    element_type __v( M_functionspace );

    __v.clear();
    fe_type* __fe = __v.functionSpace()->fe().get();

    // we should manipulate the same type of functions on the left and
    // on the right
    //BOOST_STATIC_ASSERT(( boost::is_same<return_value_type, typename functionspace_type::return_value_type>::value ));

    const uint16_type ndofv = functionspace_type::fe_type::nDof;
    CHECK( __v.size() == M_functionspace->dof()->nDof() )
        << "Invalid size : " <<  __v.size() << "!=" <<  M_functionspace->dof()->nDof();
    // assert functionspace_type::nComponents == 1
    __v.resize( M_functionspace->dof()->nDof() );
    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = M_range;

    // acquire interface (=done) cells
    // Note: Might be based on indicatorGamma from ReinitializerILP
    //       but this might not be more efficient, at the cost of an unneeded
    //       dependency. So (re-)do it here for the moment.
    // doing also initialization of done distances
    // the chosen approach assumes monotonicity of the FE-function in the
    // element, thus valid only for P1 elements
    std::set<size_type> done;
    std::vector<status_type> status( __v.size(), FAR );

    for ( ; it!=en ; ++it )
        {
            uint16_type nPlus = 0;
            uint16_type nMinus = 0;
            std::vector<size_type> indices( ndofv );
            for ( uint16_type __j = 0; __j < ndofv;++__j )
                {
                    size_type index = phi.start() + boost::get<0>(M_functionspace->dof()->localToGlobal( it->id(), __j, 0 ));
                    indices[__j] = index;
                    // std::cout<<"indices["<<__j<<"]="<<indices[__j]<<std::endl;
                    // std::cout<<"phi["<<index<<"]="<<phi[index]<<std::endl;
                    if ( phi[index] < 0.0 )
                        ++nMinus;
                    else if ( phi[index] > 0.0 )
                        ++nPlus;
                }
            //if the element has nodes with positive and negative values
            //mark as done
            if ( nPlus != ndofv && nMinus != ndofv )
                {
                    for ( uint16_type __j = 0; __j < ndofv; ++__j )
                        {
                            done.insert( indices[__j] );
                            status[indices[__j]] = DONE;
                            __v[indices[__j]] = phi[indices[__j]];
                        }
                }
        }

    Feel::details::FmsHeap<value_type> theHeap;

    // initialize close distances in heap and mark close points in status array
    for ( std::set<size_type>::iterator dit = done.begin();
          dit != done.end(); ++dit )
        {
            fmsHeapUpdate( *dit, __v, status, theHeap );
        } // loop over done points
    done.clear(); // not needed any more, save memory...

    // marching loop
    while ( theHeap.size() > 0 )
        {

            /*VD : the heap is sorted with the smallest phi value
             at the top.
             thus, the new accepted element is always the top of the heap */

            // mark closest done and save distance
            typename Feel::details::FmsHeap<value_type>::heap_entry_type
                newAccepted = theHeap.pop();
            uint16_type newIdx = newAccepted.second;
            status[newIdx] = DONE;
            __v[newIdx] = newAccepted.first;
            // update heap
            fmsHeapUpdate( newIdx, __v, status, theHeap );

            if (M_periodicity.isPeriodic())
                updatePeriodicPoint(newAccepted, __v, status, theHeap);


        } // marching loop

    return __v;

} // operator()


template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type>
void ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::
updatePeriodicPoint(typename Feel::details::FmsHeap<value_type>::heap_entry_type const& newAccepted, element_type& __v, std::vector<status_type>& status, Feel::details::FmsHeap<value_type>& theHeap) const
{
    // search if the point at newIdx has a corresponding point by periodic translation
    // if yes, mark the point as done and update the heap in consequence

    vf::node_type pt_tofindPlus(Dim);
    vf::node_type pt_tofindMinus(Dim);
    vf::node_type pt_found(Dim);

    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    uint16_type newIdx = newAccepted.second;

    for (int __c = 0 ; __c<Dim; __c++)
        {
            pt_tofindPlus [__c] = M_coords[newIdx][__c] + M_translation[__c];
            pt_tofindMinus[__c] = M_coords[newIdx][__c] - M_translation[__c];
        }

    int foundPlusMinus = 0;
    size_type id_found;

    if (M_functionspace->findPoint(pt_tofindPlus, id_found, pt_found))
        foundPlusMinus = 1;
    else if (M_functionspace->findPoint(pt_tofindMinus, id_found, pt_found))
        foundPlusMinus = -1;

    if (foundPlusMinus)
        for (int __k = 0; __k < ndofv; ++__k)
            {
                value_type index_k = boost::get<0>(M_functionspace->dof()->localToGlobal( id_found, __k, 0 ));
                value_type diff_position=0;
                for (int __c = 0; __c < Dim; __c++)
                    diff_position += std::abs(M_coords[index_k][__c] - (M_coords[newIdx][__c] + foundPlusMinus * M_translation[__c] ) );
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



template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type>
void ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::
fmsHeapUpdate( size_type idDone,
               element_type const& __v,
               std::vector<status_type>& status,
               Feel::details::FmsHeap<value_type>& theHeap ) const
{
    using namespace Feel;

    // nbr : neighbours of the node just done
    std::set<size_type> const & nbrs = M_neighbors.find(idDone)->second;
    std::set<size_type>::const_iterator n0it;
    std::vector<size_type> ids( 1, idDone );
    for ( n0it = nbrs.begin(); n0it != nbrs.end(); ++n0it )
        {
            if ( status[*n0it] == FAR )
                status[*n0it] = CLOSE;
            if ( status[*n0it] == CLOSE )
                {
                    // one neighbor
                    ids.push_back(*n0it);
                    value_type phiNew = fmsDistN( ids, __v );
                    ids.pop_back();

                    phiNew = fmsDistRec( ids, *n0it, __v, phiNew, status );

                    theHeap.change( std::make_pair( phiNew, *n0it ) );
                } // if CLOSE
        } // loop over neighbor 0
} // fmsHeapUpdate


template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type>
typename ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::value_type
ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::
fmsDistRec( std::vector<size_type> & ids,
            size_type idClose,
            element_type const& __v,
            value_type phiOld,
            std::vector<status_type> const& status ) const
{
    using namespace Feel;

    /*
      VD :
      search for all DONE neighbours
      recalculate phi with all the neighbours
      returns the smallest phi
    */

    if ( ids.size() >= Dim )
        return phiOld;

    value_type phiNew(phiOld);

    std::set<size_type> const & nbrs = M_neighbors.find(idClose)->second;
    std::set<size_type>::const_iterator nit;
    for ( nit = nbrs.begin(); nit != nbrs.end(); ++nit )
        {
            if ( status[*nit] != DONE )
                continue;
            bool unique = true;
            for ( std::vector<size_type>::const_iterator idsit=ids.begin();
                  idsit != ids.end(); ++idsit )
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


template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type>
typename ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::value_type
ReinitializerFMS<FunctionSpaceType, IteratorRange, periodicity_type>::
fmsDistN( std::vector<size_type> const & ids,
          element_type const & __v ) const
{
    using namespace Feel;

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


// 2d
template class ReinitializerFMS< spaceP1LS_type, itRangeLS, Periodic<> > ;
template class ReinitializerFMS< spaceP1LS_type, itRangeLS,  NoPeriodicity  > ;
// 3d
template class ReinitializerFMS< spaceP1LS_3d_type, itRange_3d_LS, NoPeriodicity > ;
template class ReinitializerFMS< spaceP1LS_3d_type, itRange_3d_LS, Periodic<> > ;

} // Feel
#endif
