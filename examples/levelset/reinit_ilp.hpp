/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-23

  Copyright (C) 2006 EPFL

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
/**
   \file reinit_ilp.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-23
 */

#ifndef _REINIT_ILP_HPP_
#define _REINIT_ILP_HPP_

#include <life/lifediscr/functionspace.hpp>
#include <life/lifepoly/im.hpp>
#include <life/lifevf/vf.hpp>

#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifediscr/fsfunctionallinear.hpp>

namespace Life
{

// An interface local projection reinitialization solver

template<class Space,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class ReinitializerILP
{
public:

    // -- TYPEDEFS --
    typedef Space space_type;
    typedef typename space_type::value_type value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    static const uint16_type order = space_type::basis_type::nOrder;
    // assert nOrder == 1
    static const uint16_type imOrder = 2;

    static const uint16_type Dim = space_type::nDim;


    /* matrix */
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    /* mesh */
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    /* space */
    typedef fusion::vector<fem::Lagrange<Dim, 0,
                                         Scalar, Discontinuous, value_type> >
    basis_i_type;
    typedef FunctionSpace<mesh_type, basis_i_type, value_type> space_i_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<space_i_type> space_i_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_i_type::element_type element_i_type;

    /* quadrature */
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    ReinitializerILP( const space_ptrtype& space,
                      const backend_ptrtype& backend )
        :
        M_mesh( space->mesh() ),
        M_spaceLS( space ),
        M_massLSinterf( space, space, backend ),
        M_im()
    {}

    // reinitialize
    element_type operator()( element_type const& phi,
                             element_i_type const& indicatorGamma );

private:

    mesh_ptrtype M_mesh;

    space_ptrtype M_spaceLS;

    OperatorLinear<space_type, space_type> M_massLSinterf;

    im_type M_im;

}; // class ReinitializerILP



template<class Space,
         template<uint16_type,uint16_type,uint16_type> class Entity>
typename ReinitializerILP<Space, Entity>::element_type
ReinitializerILP<Space, Entity>::operator()
    ( element_type const& phi,
      element_i_type const& indicatorGamma )
{
    element_type phiNew;
    phiNew = phi;

    using namespace Life::vf;

    M_mesh->updateMarker2( indicatorGamma );

    // matrix
//     value_type eps = type_traits<value_type>::epsilon();
    M_massLSinterf =
        integrate( marked2elements(M_mesh, 1), M_im,
                   idt(phiNew) * id(phiNew)
                   );
//         integrate( marked2elements(M_mesh, 0), M_im,
//                    eps * idt(phiNew) * id(phiNew)
//                    );
    M_massLSinterf.mat().close();

    // rhs
    FsFunctionalLinear<space_type> rhs( M_spaceLS );
    rhs =
        integrate( marked2elements(M_mesh, 1), M_im,
                   idv(phi) / sqrt( gradv(phi)*trans(gradv(phi)) )
                   * id(phiNew)
                   );
//         integrate( marked2elements(M_mesh, 0), M_im,
//                    eps * idv(phi) * id(phiNew)
//                    );

    // projection
    M_massLSinterf.applyInverse( phiNew, rhs );

    return phiNew;

} // operator()

} // Life

#endif /* _REINIT_ILP_HPP_ */
