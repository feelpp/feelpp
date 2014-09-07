/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-18

  Copyright (C) 2005,2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file reinit_fms.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-18
 */
#ifndef __Reinit_Fms_H
#define __Reinit_Fms_H 1

#include <feel/feelvf/vf.hpp>

#include "fmsheap.hpp"
#include "fmspoint.hpp"

namespace Feel
{

template<typename FunctionSpaceType, typename IteratorRange>
class ReinitializerFMS
{
public:


    /** @name Typedefs
     */
    //@{

    enum status_type { DONE=0, CLOSE=1, FAR=2 };

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    static const uint16_type Dim = geoelement_type::nDim;

    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef IteratorRange range_iterator;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                      IteratorRange const& r )
        :
        M_functionspace( __functionspace ),
        M_range( r ),
        M_neighbors(),
        M_coords( __functionspace->dof()->nDof() )
    {
        DVLOG(2) << "ReinitializerFMS constructor from space and iterator range\n";

        fe_type* __fe = __functionspace->fe().get();

        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        typename gm_type::precompute_ptrtype
        __geopc( new typename gm_type::precompute_type( __functionspace->gm(),
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

            for ( uint16_type __j = 0; __j < ndofv; ++__j )
            {
                size_type index = boost::get<0>( M_functionspace->dof()->localToGlobal( it->id(), __j, 0 ) );
                indices[__j] = index;

                for ( uint16_type c = 0; c < Dim; ++c )
                    M_coords[index][c] = tensor_expr.evalq( c, 0, __j );
            }

            for ( uint16_type __j = 0; __j < ndofv; ++__j )
                for ( uint16_type __k = __j+1; __k < ndofv; ++__k )
                {
                    M_neighbors[indices[__j]].insert( indices[__k] );
                    M_neighbors[indices[__k]].insert( indices[__j] );
                }
        }
    }

    ReinitializerFMS( ReinitializerFMS const& __vfi )
        :
        M_functionspace( __vfi.M_functionspace ),
        M_range( __vfi.M_range ),
        M_neighbors( __vfi.M_neighbors ),
        M_coords( __vfi.M_coords )
    {
        DVLOG(2) << "ReinitializerFMS copy constructor\n";
    }

    virtual ~ReinitializerFMS() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    element_type operator() ( element_type const& phi ) const;

    //@}

private:

    typedef typename functionspace_type::fe_type fe_type;

    typedef typename functionspace_type::gm_type gm_type;
    typedef typename gm_type::template Context<vm::POINT, geoelement_type>
    gm_context_type;
    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, geoelement_type> fecontext_type;

    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef __typeof__( vf::P() ) expression_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;

    typedef details::FmsPoint<value_type, Dim> point_type;
    typedef std::map<size_type, std::set<size_type> > neighbors_type;


    void fmsHeapUpdate( size_type idDone,
                        element_type const& __v,
                        std::vector<status_type>& status,
                        details::FmsHeap<value_type>& theHeap ) const;

    value_type fmsDistN( std::vector<size_type> const& ids,
                         element_type const& __v ) const;

    value_type fmsDistRec( std::vector<size_type> & ids,
                           size_type idClose,
                           element_type const& __v,
                           value_type phiOld,
                           std::vector<status_type> const& status ) const;

    value_type closerOne( value_type a, value_type b ) const
    {
        return a*a < b*b ? a : b;
    }

    functionspace_ptrtype const& M_functionspace;
    range_iterator M_range;
    neighbors_type M_neighbors;
    std::vector<point_type> M_coords;

};

template<typename FunctionSpaceType, typename Iterator>
typename ReinitializerFMS<FunctionSpaceType, Iterator>::element_type
ReinitializerFMS<FunctionSpaceType, Iterator>::operator()
( element_type const& phi ) const
{

    //     VLOG(1) << "[ReinitFMS] operator()\n";

    element_type __v( M_functionspace );
    __v.clear();

    fe_type* __fe = __v.functionSpace()->fe().get();

    // we should manipulate the same type of functions on the left and
    // on the right
    //BOOST_STATIC_ASSERT(( boost::is_same<return_value_type, typename functionspace_type::return_value_type>::value ));

    const uint16_type ndofv = functionspace_type::fe_type::nDof;
    FEELPP_ASSERT( __v.size() == M_functionspace->dof()->nDof() )( __v.size() )( M_functionspace->dof()->nDof() ).warn( "invalid size" );
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

        for ( uint16_type __j = 0; __j < ndofv; ++__j )
        {
            size_type index = phi.start() + boost::get<0>( M_functionspace->dof()->localToGlobal( it->id(), __j, 0 ) );
            indices[__j] = index;

            if ( phi[index] < 0.0 )
                ++nMinus;

            else if ( phi[index] > 0.0 )
                ++nPlus;
        }

        if ( nPlus != ndofv && nMinus != ndofv )
            for ( uint16_type __j = 0; __j < ndofv; ++__j )
            {
                done.insert( indices[__j] );
                status[indices[__j]] = DONE;
                __v[indices[__j]] = phi[indices[__j]];
            }
    }

    details::FmsHeap<value_type> theHeap;

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
        // mark closest done and save distance
        typename details::FmsHeap<value_type>::heap_entry_type
        newAccepted = theHeap.pop();
        uint16_type newIdx = newAccepted.second;
        status[newIdx] = DONE;
        __v[newIdx] = newAccepted.first;

        // update heap
        fmsHeapUpdate( newIdx, __v, status, theHeap );
    } // marching loop

    //std::cout << "v = " << __v << "\n";
    return __v;

} // operator()


template<typename FunctionSpaceType, typename Iterator>
void ReinitializerFMS<FunctionSpaceType, Iterator>::
fmsHeapUpdate( size_type idDone,
               element_type const& __v,
               std::vector<status_type>& status,
               details::FmsHeap<value_type>& theHeap ) const
{
    std::set<size_type> const & nbrs = M_neighbors.find( idDone )->second;
    std::set<size_type>::const_iterator n0it;
    std::vector<size_type> ids( 1, idDone );

    for ( n0it = nbrs.begin(); n0it != nbrs.end(); ++n0it )
    {
        if ( status[*n0it] == FAR )
            status[*n0it] = CLOSE;

        if ( status[*n0it] == CLOSE )
        {
            // one neighbor
            ids.push_back( *n0it );
            value_type phiNew = fmsDistN( ids, __v );
            ids.pop_back();

            phiNew = fmsDistRec( ids, *n0it, __v, phiNew, status );

            theHeap.change( std::make_pair( phiNew, *n0it ) );
        } // if CLOSE
    } // loop over neighbor 0
} // fmsHeapUpdate


template<typename FunctionSpaceType, typename Iterator>
typename ReinitializerFMS<FunctionSpaceType, Iterator>::value_type
ReinitializerFMS<FunctionSpaceType, Iterator>::
fmsDistRec( std::vector<size_type> & ids,
            size_type idClose,
            element_type const& __v,
            value_type phiOld,
            std::vector<status_type> const& status ) const
{
    if ( ids.size() >= Dim )
        return phiOld;

    value_type phiNew( phiOld );

    std::set<size_type> const & nbrs = M_neighbors.find( idClose )->second;
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

        if ( M_neighbors.find( ids[0] )->second.find( *nit ) ==
                M_neighbors.find( ids[0] )->second.end() )
            continue;

        // one neighbor more
        ids.push_back( *nit );

        ids.push_back( idClose );
        value_type phiCand = fmsDistN( ids, __v );
        ids.pop_back();

        if ( phiCand != 0.0 )
            phiNew = closerOne( phiCand, phiNew );

        phiNew = fmsDistRec( ids, idClose, __v, phiNew, status );

        ids.pop_back();
    }

    return phiNew;
} // fmsDistRec


template<typename FunctionSpaceType, typename Iterator>
typename ReinitializerFMS<FunctionSpaceType, Iterator>::value_type
ReinitializerFMS<FunctionSpaceType, Iterator>::
fmsDistN( std::vector<size_type> const & ids,
          element_type const & __v ) const
{
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

        wNorm[i] = 1.0/norm( basis[i] );
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

            n[i] = sqrt( n_rest );

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

} // namespace Feel

#endif /* __Reinit_Fms_H */
