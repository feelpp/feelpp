/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-04

  Copyright (C) 2006 EPFL

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
   \file advreact.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-04
 */

#ifndef _ADVREACT_HPP_
#define _ADVREACT_HPP_

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/operatorlinear.hpp>

namespace Feel
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
        M_rhsStab( space ),
        M_updated( false ),
        M_stabcoeff( 0.1 * std::pow( polyOrder, -3.5 ) ),
        M_im()
    {
        M_StabMethod="SGS";
        DoStabilize=true;
    }

    // setting of options
    void set_stabcoeff( double stabcoeff )
    {
        M_stabcoeff = stabcoeff * std::pow( polyOrder, -3.5 );
    }

    void set_StabMethod( std::string Method )
    {
        //possible : CIP, GLS, SGS, SUPG
        M_StabMethod=Method;
    }

    // update operator and rhs with given expressions
    template<typename Esigma, typename Ebeta,
             typename Ef, typename Eg>
    void update( const Esigma& sigma,
                 const Ebeta& beta,
                 const Ef& f,
                 const Eg& g,
                 bool updateStabilization = true
                 );

    // solve system, call not needed, but possible
    void solve();

    // result
    const element_type phi()
    {
        solve();
        return M_phi;
    }

    bool DoStabilize;

private:

    mesh_ptrtype M_mesh;

    space_ptrtype M_space;

    element_type M_phi;

    OperatorLinear<space_type, space_type> M_operator;
    OperatorLinear<space_type, space_type> M_operatorStab;
    FsFunctionalLinear<space_type> M_rhs;
    FsFunctionalLinear<space_type> M_rhsStab;

    bool M_updated;

    double M_stabcoeff;

    im_type M_im;

    std::string M_StabMethod;

}; // class AdvReact


template<class Space, uint16_type imOrder,
         template<uint16_type,uint16_type,uint16_type> class Entity>
template<typename Esigma, typename Ebeta,
         typename Ef, typename Eg>
void AdvReact<Space, imOrder, Entity>::update(const Esigma& sigma,
                                              const Ebeta& beta,
                                              const Ef& f,
                                              const Eg& g,
                                              bool updateStabilization
                                              )
{
    M_updated = true;

    using namespace Feel::vf;

    M_operator =
        integrate( elements(M_mesh), M_im,
                   ( val(sigma)*idt(M_phi) + gradt(M_phi)*val(beta) ) * id(M_phi)
                   ) +
        integrate( boundaryfaces(M_mesh), M_im,
                   - chi( trans(beta)*N() < 0 ) * /* inflow */
                   val(trans(beta)*N()) * idt(M_phi) * id(M_phi)
                   );

    if (DoStabilize)
        {
            if ( updateStabilization )
                {
                    //good review of stabilization methods in [Chaple 2006]
                    if (M_StabMethod=="CIP" && (M_stabcoeff != 0.0) )
                                {
                                    /* don't work properly in 3D (because of internalfaces) */
                                    M_operatorStab=
                                        integrate( internalfaces(M_mesh), M_im,
                                                   val( M_stabcoeff*vf::pow(hFace(),2.0) *
                                                        abs(trans(beta)*N()) ) *
                                                   (jumpt(gradt(M_phi)) * jump(grad(M_phi)))
                                                   );
                                }//Continuous Interior Penalty

                    else if (M_StabMethod== "SUPG")
                        {
                            AUTO(coeff, vf::h()/(2*vf::sqrt(val(trans(beta))*val(beta))));
                            AUTO(L_op, ( grad(M_phi)*val(beta) ) );
                            AUTO(L_opt, ( gradt(M_phi)*val(beta) + val(sigma)*idt(M_phi) ) );

                            M_operatorStab=
                                integrate(elements(M_mesh), M_im,
                                          coeff
                                          * L_op
                                          * L_opt );

                            M_rhsStab=
                                integrate(elements(M_mesh), M_im,
                                          coeff*L_op*val(f));
                        }//Streamline Upwind Petrov Galerkin

                    else if (M_StabMethod== "GLS")
                        {
                            AUTO(coeff, 1.0 / (2*vf::sqrt(val(trans(beta))*val(beta))/vf::h()+vf::abs(val(sigma))));
                            AUTO(L_op, ( grad(M_phi)*val(beta) + val(sigma)*id(M_phi) ) );
                            AUTO(L_opt, ( gradt(M_phi)*val(beta) + val(sigma)*idt(M_phi) ) );

                            M_operatorStab=
                                integrate( elements(M_mesh), M_im,
                                           coeff
                                           * L_op
                                           * L_opt );

                            M_rhsStab=
                                integrate(elements(M_mesh), M_im,
                                          coeff*L_op*val(f));

                        }//Galerkin Least Square

                    else if (M_StabMethod== "SGS")
                        {
                            AUTO(coeff, 1.0 / (2*vf::sqrt(val(trans(beta))*val(beta))/vf::h()+vf::abs(val(sigma))));
                            AUTO(L_op, (grad(M_phi)* val(beta) - val(sigma) * id(M_phi) ) );
                            AUTO(L_opt, ( gradt(M_phi)*val(beta) + val(sigma)*idt(M_phi) ) );

                            M_operatorStab=
                                integrate(elements(M_mesh), M_im,
                                          coeff
                                          * L_op
                                          * L_opt );
                            M_rhsStab=
                                integrate(elements(M_mesh), M_im,
                                          coeff*L_op*val(f));
                        }//Subgrid Scale method
                }//update stabilization
            M_operator.add(1.0, M_operatorStab);
        }//DoStabilize
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
    if (DoStabilize)
        M_rhs.add(M_rhsStab);

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

} // Feel

#endif /* _ADVREACT_HPP_ */
