/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-01-09

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
   \file indicator.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-01-09
 */

#ifndef _INDICATOR_HPP_
#define _INDICATOR_HPP_

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feeldiscr/fsfunctionallinear.hpp>

namespace Feel
{

// Class collecting

template<class Space,
         template<uint16_type,uint16_type,uint16_type> class Entity>
class Indicator
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

    Indicator( const space_ptrtype& space,
               const space_i_ptrtype& indicatorSpace,
               const backend_ptrtype& backend )
        :
        M_mesh( space->mesh() ),
        M_spaceLS( space ),
        M_spaceIndic( indicatorSpace ),
        M_kappa( indicatorSpace, "kappa" ),
        M_indicatorGamma( indicatorSpace, "indicatorGamma" ),
        M_indicatorChange( indicatorSpace, "indicatorChange" ),
        M_massIndic( indicatorSpace, indicatorSpace, backend ),
        M_im()
    {
        using namespace Feel::vf;

        // build P0 mass matrix once only, at construction
        M_massIndic =
            integrate( elements( M_mesh ), M_im,
                       idt( M_kappa )*id( M_kappa )
                     );
        M_massIndic.mat().close();
    }

    // update indicator
    void update( element_type const& phiNew,
                 element_type const& phiOld );

    // update indicator
    void update( element_type const& phi );

    // return P0 indicator function of interface region
    // (=elements cut by interface)
    element_i_type indicatorGamma() const
    {
        return M_indicatorGamma;
    }

    // return P0 indicator function kappa with following properties:
    // * kappa = -1 on "-" elements
    // * -1 < kappa < +1 on elements cut by interface
    // * kappa = +1 on "+" elements
    element_i_type kappa() const
    {
        return M_kappa;
    }

    // return P0 indicator function of regions where sign of phi has changed
    element_i_type indicatorChange() const
    {
        return M_indicatorChange;
    }

private:

    void updateKappa( element_type const& phi );

    mesh_ptrtype M_mesh;

    space_ptrtype M_spaceLS;
    space_i_ptrtype M_spaceIndic;

    element_i_type M_kappa;
    element_i_type M_indicatorGamma;
    element_i_type M_indicatorChange;

    OperatorLinear<space_i_type, space_i_type> M_massIndic;

    im_type M_im;

}; // class Indicator



template<class Space,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void
Indicator<Space, Entity>::updateKappa( element_type const& phi )
{
    using namespace Feel::vf;

    // --- sigma: node-wise sign of phi as P1 function
    element_type sigma(  M_spaceLS, "sigma" );
    sigma = vf::project( M_spaceLS,
                         elements( M_mesh ),
                         sign( idv( phi ) )
                       );

    // --- kappa: element-wise mean of sigma as P0-indicator function
    FsFunctionalLinear<space_i_type> rhsKappa( M_spaceIndic );
    rhsKappa =
        integrate( elements( M_mesh ), M_im,
                   idv( sigma )*id( M_kappa )
                 );
    M_massIndic.applyInverse( M_kappa, rhsKappa );
}


template<class Space,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void
Indicator<Space, Entity>::update( element_type const& phi )
{
    updateKappa( phi );

    using namespace Feel::vf;

    // --- indicatorGamma: 1 in elements crossed by interface, 0 elsewhere
    M_indicatorGamma = vf::project( M_spaceIndic,
                                    elements( M_mesh ),
                                    chi( ( idv( M_kappa ) <  1.0-0.5/( Dim+1 ) ) &&
                                         ( idv( M_kappa ) > -1.0+0.5/( Dim+1 ) ) ) );
} // Indicator::update


template<class Space,
         template<uint16_type,uint16_type,uint16_type> class Entity>
void
Indicator<Space, Entity>::update( element_type const& phiNew,
                                  element_type const& phiOld )
{
    updateKappa( phiNew );

    using namespace Feel::vf;

    // --- sigmaDelta: node-wise sign of phiOld*phiNew as P1 function
    element_type sigmaDelta(  M_spaceLS, "sigmaDelta" );
    sigmaDelta = vf::project( M_spaceLS,
                              elements( M_mesh ),
                              sign( idv( phiOld )*idv( phiNew ) )
                            );

    // --- kappaDelta: element-wise mean of sigmaDelta as P0-indicator function
    element_i_type kappaDelta( M_spaceIndic, "kappaDelta" );
    FsFunctionalLinear<space_i_type> rhsKappa( M_spaceIndic );
    rhsKappa =
        integrate( elements( M_mesh ), M_im,
                   idv( sigmaDelta )*id( kappaDelta )
                 );
    M_massIndic.applyInverse( kappaDelta, rhsKappa );

    // --- indicatorChange: 1 in elements where phi changed sign, 0 elsewhere
    M_indicatorChange =
        vf::project( M_spaceIndic,
                     elements( M_mesh ),
                     chi( ( ( idv( M_kappa ) <  1.0-0.5/( Dim+1 ) ) &&
                            ( idv( M_kappa ) > -1.0+0.5/( Dim+1 ) ) ) ||
                          ( idv( kappaDelta ) < 1.0-0.5/( Dim+1 ) )   )
                   );

} // Indicator::update(...)

} // Feel

#endif /* _INDICATOR_HPP_ */
