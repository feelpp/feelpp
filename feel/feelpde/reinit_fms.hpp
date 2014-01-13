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
/**
   big changes : handle parallel and a lot of code simplifications
   Vincent Doyeux
   2014-01-10
 */

#ifndef __Reinit_Fms_H
#define __Reinit_Fms_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelvf/vf.hpp>

#include "fmsheap.hpp"
#include "fmspoint.hpp"
#include <boost/serialization/set.hpp>

namespace Feel
{

template<typename FunctionSpaceType, typename periodicity_type = NoPeriodicity>
class ReinitializerFMS
{
public:

    /** @name Typedefs
     */
    //@{

    enum status_type { FAR=0, CLOSE=1, DONE=2 };

    typedef ReinitializerFMS<FunctionSpaceType, periodicity_type> self_type;
    typedef boost::shared_ptr< self_type > self_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    static const uint16_type Dim = geoelement_type::nDim;


    ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                      periodicity_type __periodicity=NoPeriodicity());


    static self_ptrtype New(  functionspace_ptrtype const& __functionspace,
                              periodicity_type __periodicity=NoPeriodicity() )
    {
        self_ptrtype fm( new self_type(__functionspace,  __periodicity) );
        return fm;
    }

    virtual ~ReinitializerFMS() {}


    element_type operator() ( element_type const& phi, bool useMarker2AsDoneMarker=false );

    // same as operator()
    element_type march ( element_type const& phi, bool useMarker2AsDoneMarker=false )
    { return this->operator()( phi, useMarker2AsDoneMarker) ; }


private:

    typedef Feel::details::FmsPoint<value_type, Dim> point_type;

    /* The map indexes are the global indexes on the PROC and the set contains the global indexes on the CLUSTER of its neigbors */
    typedef std::map<size_type, std::set<size_type> > neighbors_type;

    inline size_type clusterToProcessor( size_type dof )
    { return M_functionspace->dof()->mapGlobalClusterToGlobalProcess( dof - firstDof ); }

    inline size_type processorToCluster( size_type dof )
    { return M_functionspace->dof()->mapGlobalProcessToGlobalCluster( dof ); }

    void reduceDonePoints(element_type const& __v, Feel::details::FmsHeap<value_type>& theHeap, element_type& status, std::set<size_type>& done );

    void updatePeriodicPoint(typename Feel::details::FmsHeap<value_type>::heap_entry_type const& newAccepted, element_type& __v, element_type& status, Feel::details::FmsHeap<value_type>& theHeap) const;

    void fmsHeapUpdate( size_type idDone,
                        element_type const& __v,
                        element_type& status,
                        Feel::details::FmsHeap<value_type>& theHeap ) const;

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
    periodicity_type M_periodicity;
    neighbors_type M_neighbors;
    std::map< size_type, size_type> M_ghostClusterToProc;
    std::vector<point_type> M_coords;
    vf::node_type M_translation;
    const size_type firstDof;

};

// Instantiate 2d and 3d only when compiling reinit_fms.cpp (where FEELPP_INSTANTIATE_FMS is defined), else "extern" avoid useless instantiation
typedef Feel::bases<Feel::Lagrange<1, Feel::Scalar> > basisP1LS_type;

//2d
typedef Feel::Mesh< Feel::Simplex<2> > mesh_typeLS;
typedef Feel::FunctionSpace<mesh_typeLS, basisP1LS_type, Feel::Periodicity <Feel::NoPeriodicity> > spaceP1LS_type;

// 3d
typedef Feel::Mesh< Feel::Simplex<3> > mesh_3d_typeLS;
typedef Feel::FunctionSpace<mesh_3d_typeLS, basisP1LS_type, Feel::Periodicity <Feel::NoPeriodicity> > spaceP1LS_3d_type;

#if !defined( FEELPP_INSTANTIATE_FMS )
extern template class Feel::ReinitializerFMS< spaceP1LS_type, Feel::Periodic<> > ;
extern template class Feel::ReinitializerFMS< spaceP1LS_type, Feel::NoPeriodicity  > ;
extern template class Feel::ReinitializerFMS< spaceP1LS_3d_type, Feel::NoPeriodicity > ;
extern template class Feel::ReinitializerFMS< spaceP1LS_3d_type, Feel::Periodic<> > ;
#endif


template<typename FunctionSpaceType, typename periodicity_type = NoPeriodicity>
boost::shared_ptr< ReinitializerFMS< FunctionSpaceType, periodicity_type > >
fms ( boost::shared_ptr<FunctionSpaceType> const& Xh, periodicity_type p = NoPeriodicity() )
{
    auto fm = ReinitializerFMS< FunctionSpaceType, periodicity_type >::New( Xh, p );
    return fm;
}

} // namespace Feel

#endif /* __Reinit_Fms_H */
