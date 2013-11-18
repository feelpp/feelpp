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

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelpde/fmsheap.hpp>
#include <feel/feelpde/fmspoint.hpp>


namespace Feel
{

template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type = NoPeriodicity>
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
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    static const uint16_type Dim = geoelement_type::nDim;

    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef IteratorRange range_iterator;

    ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                      IteratorRange const& r ,
                      periodicity_type __periodicity=NoPeriodicity());


    // ReinitializerFMS( ReinitializerFMS const& __vfi )
    //     :
    //     M_functionspace( __vfi.M_functionspace ),
    //     M_range( __vfi.M_range ),
    //     M_neighbors( __vfi.M_neighbors ),
    //     M_coords( __vfi.M_coords )
    //     {
    //         Debug( 5065 ) << "ReinitializerFMS copy constructor\n";
    //     }

    virtual ~ReinitializerFMS() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    element_type operator() ( element_type const& phi) const;
    //@}

private:

    typedef typename functionspace_type::fe_type fe_type;

    typedef typename functionspace_type::gm_type gm_type;
    typedef typename gm_type::template Context<vm::POINT, geoelement_type>
    gm_context_type;
    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, geoelement_type> fecontext_type;

    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef decltype( vf::P() ) expression_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;

    typedef Feel::details::FmsPoint<value_type, Dim> point_type;
    typedef std::map<size_type, std::set<size_type> > neighbors_type;

    void updatePeriodicPoint(typename Feel::details::FmsHeap<value_type>::heap_entry_type const& newAccepted,
                             element_type& __v,
                             std::vector<status_type>& status,
                             Feel::details::FmsHeap<value_type>& theHeap) const;

    void fmsHeapUpdate( size_type idDone,
                        element_type const& __v,
                        std::vector<status_type>& status,
                        Feel::details::FmsHeap<value_type>& theHeap ) const;

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
    periodicity_type M_periodicity;
    neighbors_type M_neighbors;
    std::vector<point_type> M_coords;
    vf::node_type M_translation;
};

// needed for instantiation (copied from levelset.hpp)
typedef Mesh< Simplex<2> > mesh_typeLS;
typedef boost::shared_ptr<mesh_typeLS> mesh_ptrtypeLS;
typedef bases<Lagrange<1, Scalar> > basisP1LS_type;
typedef FunctionSpace<mesh_typeLS, basisP1LS_type, Periodicity <NoPeriodicity> > spaceP1LS_type;
mesh_ptrtypeLS ___mesh_LS;
typedef decltype( elements( ___mesh_LS ) ) itRangeLS;

typedef ReinitializerFMS< spaceP1LS_type, itRangeLS, Periodic<> > fms_2d_type;
typedef ReinitializerFMS< spaceP1LS_type, itRangeLS, NoPeriodicity > fms_2d_periodic_type;

#if !defined( FEELPP_INSTANTIATE_FMS )
extern template class ReinitializerFMS< spaceP1LS_type, itRangeLS, Periodic<> > ;
extern template class ReinitializerFMS< spaceP1LS_type, itRangeLS,  NoPeriodicity  > ;
#endif

// 3d
typedef Mesh< Simplex<3> > mesh_3d_typeLS;
typedef boost::shared_ptr<mesh_3d_typeLS> mesh_3d_ptrtypeLS;
typedef FunctionSpace<mesh_3d_typeLS, basisP1LS_type, Periodicity <NoPeriodicity> > spaceP1LS_3d_type;
mesh_3d_ptrtypeLS ___mesh_3d_LS;
typedef decltype( elements( ___mesh_3d_LS ) ) itRange_3d_LS;

typedef ReinitializerFMS< spaceP1LS_3d_type, itRange_3d_LS, NoPeriodicity > fms_3d_type;
typedef  ReinitializerFMS< spaceP1LS_3d_type, itRange_3d_LS, Periodic<> > fms_3d_periodic_type;

#if !defined( FEELPP_INSTANTIATE_FMS )
extern template class ReinitializerFMS< spaceP1LS_3d_type, itRange_3d_LS, NoPeriodicity > ;
extern template class ReinitializerFMS< spaceP1LS_3d_type, itRange_3d_LS, Periodic<> > ;
#endif

template<typename FunctionSpaceType, typename IteratorRange, typename periodicity_type = NoPeriodicity>
boost::shared_ptr<ReinitializerFMS< FunctionSpaceType, IteratorRange, periodicity_type>>
fms( boost::shared_ptr<FunctionSpaceType> const& Xh,
     IteratorRange r,
     periodicity_type p = NoPeriodicity() )
{
    boost::shared_ptr<ReinitializerFMS< FunctionSpaceType, IteratorRange, periodicity_type>> ret( new ReinitializerFMS< FunctionSpaceType, IteratorRange, periodicity_type>( Xh, r ) );
    return ret;
}

} // namespace Feel

#endif /* __Reinit_Fms_H */
