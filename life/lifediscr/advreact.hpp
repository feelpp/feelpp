/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-04

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
   \file advreact.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-04
 */

#ifndef _ADVREACT_HPP_
#define _ADVREACT_HPP_

#include <life/lifediscr/functionspace.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifevf/vf.hpp>

#include <life/lifediscr/operatorlinear.hpp>

namespace Life
{

/**
 * \class AdvReact
 * \brief Advection-Reaction solver
 *
 * An interior penalty stabilized advection-reaction solver.
 * It solves the following equation:
 * sigma * phi + beta^T grad phi = f on the domain
 *                           phi = g on the inflow boundary
 * the inflow boundary is detected through beta
 *
 */
template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class AdvReact
{
public:

    // -- TYPEDEFS --
    typedef Space space_type;
    typedef typename space_type::value_type value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    static const uint16_type Dim = space_type::nDim;


    static const value_type polyOrder = space_type::basis_type::nOrder;

    /* mesh */
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;

    /* space */
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /* quadrature */
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    AdvReact( const space_ptrtype& space,
              const backend_ptrtype& backend )
        :
        M_mesh( space->mesh() ),
        M_space( space ),
        M_phi( space, "phi" ),
        M_operator( space, space, backend ),
        M_operatorStab( space, space, backend ),
        M_rhs( space ),
        M_updated( false ),
        M_stabcoeff( 0.1 * std::pow( polyOrder, -3.5 ) ),
        M_im()
    {
    }

    // setting of options
    void set_stabcoeff( double stabcoeff )
    {
        M_stabcoeff = stabcoeff * std::pow( polyOrder, -3.5 );
    }

    // update operator and rhs with given expressions
    template<typename Esigma, typename Ebeta,
             typename Ef, typename Eg>
    void update( const Esigma& sigma,
                 const Ebeta& beta,
                 const Ef& f,
                 const Eg& g,
                 bool updateStabilization = true );

    // solve system, call not needed, but possible
    void solve();

    // result
    const element_type phi()
    {
        solve();
        return M_phi;
    }

private:

    mesh_ptrtype M_mesh;

    space_ptrtype M_space;

    element_type M_phi;

    OperatorLinear<space_type, space_type> M_operator;
    OperatorLinear<space_type, space_type> M_operatorStab;
    FsFunctionalLinear<space_type> M_rhs;

    bool M_updated;

    double M_stabcoeff;

    im_type M_im;

}; // class AdvReact



template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename Esigma, typename Ebeta,
         typename Ef, typename Eg>
void AdvReact<Space, imOrder, Entity>::update(const Esigma& sigma,
                                              const Ebeta& beta,
                                              const Ef& f,
                                              const Eg& g,
                                              bool updateStabilization)
{
    M_updated = true;

    using namespace Life::vf;

    if ( M_stabcoeff != 0.0 )
    {
        if ( updateStabilization )
        {
            M_operatorStab =
                integrate( internalfaces(M_mesh), M_im,
                           val( M_stabcoeff*vf::pow(hFace(),2.0) *
                                abs(trans(beta)*N()) ) *
                           (jumpt(gradt(M_phi)) * jump(grad(M_phi)))
                           );
        }
        M_operator = M_operatorStab;
    }

    M_operator +=
        integrate( elements(M_mesh), M_im,
                   ( val(sigma)*idt(M_phi) + gradt(M_phi)*val(beta) ) * id(M_phi)
                   ) +
        integrate( boundaryfaces(M_mesh), M_im,
                   - chi( trans(beta)*N() < 0 ) * /* inflow */
                   val(trans(beta)*N()) * idt(M_phi) * id(M_phi)
                   );
    M_operator.mat().close();
//     M_operator.mat().printMatlab("M_advReact.m");

    M_rhs =
        integrate( elements(M_mesh), M_im,
                   val(f) * id(M_phi)
                   ) +
        integrate( boundaryfaces(M_mesh), M_im,
                   - chi( trans(beta)*N() < 0 ) * /* inflow */
                   val( (trans(beta)*N())*(g) ) * id(M_phi)
                   );

//     M_rhs.container().printMatlab("F_advReact.m");

} // update

template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void AdvReact<Space, imOrder, Entity>::solve()
{
    // -- make sure solve is needed
    if ( !M_updated )
    {
        return;
    }
    M_updated = false;

    // -- solve
    M_operator.applyInverse( M_phi, M_rhs );

} // AdvReact::solve

} // Life

#endif /* _ADVREACT_HPP_ */
