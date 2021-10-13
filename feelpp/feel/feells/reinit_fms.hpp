/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2012-2016 Feel++ Consortium

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
/**
   big changes : handle parallel and a lot of code simplifications
   Vincent Doyeux
   2014-01-10
 */

#ifndef __Reinit_Fms_H
#define __Reinit_Fms_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelmesh/filters.hpp>

#include <feel/feells/fmsheap.hpp>
#include <feel/feells/fmspoint.hpp>
#include <feel/feells/lstypes.hpp>
#include <boost/serialization/set.hpp>

namespace Feel
{
/**
 * @brief Fast Marching algorithm
 * @details implements the fast marching algorithm
 *
 * @tparam FunctionSpaceType function space type
 * @tparam periodicity_type = NoPeriodicity periodicity of the distance function
 */
template<typename FunctionSpaceType, typename periodicity_type = NoPeriodicity>
class ReinitializerFMS
{
public:

    static_assert( FunctionSpaceType::fe_type::nOrder == 1, "FunctionSpaceType needs to be a finite element space of order 1");
    static_assert( ! FunctionSpaceType::is_periodic , "Space for fast marching must be non periodic, but periodicity can be given as second template argument");

    /** @name Typedefs
     */
    //@{

    enum status_type {FAR=0, CLOSE=1, DONE=2};

    typedef ReinitializerFMS<FunctionSpaceType, periodicity_type> self_type;
    typedef std::shared_ptr< self_type > self_ptrtype;

    typedef Backend<double> backend_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    using size_type = typename functionspace_type::size_type;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    static const uint16_type Dim = geoelement_type::nDim;

    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    typedef typename periodicity_type::node_type node_type;

    ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                      periodicity_type __periodicity=NoPeriodicity());


    static self_ptrtype New(  functionspace_ptrtype const& __functionspace,
                              periodicity_type __periodicity=NoPeriodicity() )
    {
        self_ptrtype fm( new self_type(__functionspace,  __periodicity) );
        return fm;
    }

    virtual ~ReinitializerFMS() {}


    //element_type operator() ( element_type const& phi, bool useMarker2AsDoneMarker=false );
    element_type operator() ( element_type const& phi, range_elements_type const& rangeInitialElts );
    element_type operator() ( element_type const& phi );

    // same as operator()
    //element_type march ( element_type const& phi, bool useMarker2AsDoneMarker=false )
    //{ return this->operator()( phi, useMarker2AsDoneMarker) ; }
    element_type march ( element_type const& phi, range_elements_type const& rangeInitialElts )
    { return this->operator()( phi, rangeInitialElts) ; }
    element_type march ( element_ptrtype phi, range_elements_type const& rangeInitialElts )
    { return this->operator()( *phi, rangeInitialElts) ; }
    element_type march ( element_type const& phi )
    { return this->operator()( phi ) ; }
    element_type march ( element_ptrtype phi )
    { return this->operator()( *phi ) ; }


private:

    typedef Feel::details::FmsHeap<value_type> heap_type;
    typedef typename heap_type::heap_entry_type heap_entry_type;
    typedef Feel::details::FmsPoint<value_type, Dim> point_type;

    /* The map indexes are the global indexes on the PROC and the set contains the global indexes on the CLUSTER of its neigbors */
    typedef std::map<size_type, std::set<size_type> > neighbors_type;

    inline size_type clusterToProcessor( size_type dof )
    { return dof - firstDof; }

    inline size_type processorToCluster( size_type dof )
    { return M_functionspace->dof()->mapGlobalProcessToGlobalCluster( dof ); }

    void createPeriodicCorrespondanceTable();

    void reduceDonePoints(element_type const& __v, element_type& status, std::set<size_type>& done );

    void reduceClosePoints(heap_type& theHeap, element_type& status );

    void fmsHeapUpdate( size_type idDone,
                        element_type const& __v,
                        element_type& status,
                        heap_type& theHeap ) const;

    value_type fmsDistN( std::vector<size_type> const& ids,
                         element_type const& __v ) const;

    value_type fmsDistRec( std::vector<size_type> & ids,
                           size_type idClose,
                           element_type const& __v,
                           value_type phiOld,
                           element_type const& status ) const;

    value_type closerOne( value_type a, value_type b ) const
    {
        return a*a < b*b ? a : b;
    }

    functionspace_ptrtype const& M_functionspace;
    vector_ptrtype checkStatus;
    vector_ptrtype valueAtClose;
    periodicity_type M_periodicity;
    neighbors_type M_neighbors;
    std::map< size_type, size_type> M_ghostClusterToProc;
    boost::bimap< size_type, size_type > M_idTag1_idTag2;
    std::vector<point_type> M_coords;
    node_type M_translation;
    const size_type firstDof;
    int M_nbDofTag1;
    int nbTotalDone;
};



#if !defined( FEELPP_INSTANTIATE_FMS )
extern template class Feel::ReinitializerFMS< ls_space_type<2>, Feel::Periodic<> > ;
extern template class Feel::ReinitializerFMS< ls_space_type<2>, Feel::NoPeriodicity  > ;
extern template class Feel::ReinitializerFMS< ls_space_type<2,2>, Feel::NoPeriodicity  > ;
extern template class Feel::ReinitializerFMS< ls_space_type<3>, Feel::Periodic<> > ;
extern template class Feel::ReinitializerFMS< ls_space_type<3>, Feel::NoPeriodicity > ;
extern template class Feel::ReinitializerFMS< ls_space_type<3,2>, Feel::NoPeriodicity > ;
// Hypercube
extern template class Feel::ReinitializerFMS< lsh_space_type<2>, Feel::NoPeriodicity  > ;
extern template class Feel::ReinitializerFMS< lsh_space_type<3>, Feel::NoPeriodicity > ;
#endif


template<typename FunctionSpaceType, typename periodicity_type = NoPeriodicity>
std::shared_ptr< ReinitializerFMS< FunctionSpaceType, periodicity_type > >
fms ( std::shared_ptr<FunctionSpaceType> const& Xh, periodicity_type p = NoPeriodicity() )
{
    auto fm = ReinitializerFMS< FunctionSpaceType, periodicity_type >::New( Xh, p );
    return fm;
}

} // namespace Feel

#endif /* __Reinit_Fms_H */
