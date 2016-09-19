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

#include <feel/feells/fastmarchingbase.hpp>

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
class ReinitializerFMS : public FastMarchingBase<FunctionSpaceType, periodicity_type>
{
public:
    typedef FastMarchingBase<FunctionSpaceType, periodicity_type> super_type;
    typedef ReinitializerFMS<FunctionSpaceType, periodicity_type> self_type;
    typedef boost::shared_ptr< self_type > self_ptrtype;

    typedef Backend<double> backend_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;
    typedef typename functionspace_type::value_type value_type;

    typedef typename super_type::heap_data_type heap_data_type;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type geoelement_type;
    static const uint16_type Dim = geoelement_type::nDim;

    typedef typename periodicity_type::node_type node_type;

    ReinitializerFMS( functionspace_ptrtype const& __functionspace,
                      periodicity_type __periodicity=NoPeriodicity());


    static self_ptrtype New(  functionspace_ptrtype const& __functionspace,
                              periodicity_type __periodicity=NoPeriodicity() )
    {
        self_ptrtype fm( new self_type(__functionspace,  __periodicity) );
        return fm;
    }

    element_type operator() ( element_type const& phi, bool useMarker2AsDoneMarker=false )
    {
        this->run( phi, useMarker2AsDoneMarker );
        return *(this->getDistance());
    }

    // same as operator()
    element_type march ( element_type const& phi, bool useMarker2AsDoneMarker=false )
    { return this->operator()( phi, useMarker2AsDoneMarker) ; }

    element_type march ( element_ptrtype phi, bool useMarker2AsDoneMarker=false )
    { return this->operator()( *phi, useMarker2AsDoneMarker) ; }


private:
    void processDof( size_type idOnProc, value_type val, heap_data_type const& opt_data );
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
boost::shared_ptr< ReinitializerFMS< FunctionSpaceType, periodicity_type > >
fms ( boost::shared_ptr<FunctionSpaceType> const& Xh, periodicity_type p = NoPeriodicity() )
{
    auto fm = ReinitializerFMS< FunctionSpaceType, periodicity_type >::New( Xh, p );
    return fm;
}

} // namespace Feel

#endif /* __Reinit_Fms_H */
