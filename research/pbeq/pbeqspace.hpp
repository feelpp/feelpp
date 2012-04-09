/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-07-11

  Copyright (C) 2007 Unil

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
   \file pbeqspace.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-07-11
 */

#ifndef _PBEQSPACE_HPP
#define _PBEQSPACE_HPP

#include <boost/numeric/ublas/vector.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/vector.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>

#include "pbeqapplication.hpp"
#include "molecule.hpp"


namespace Feel
{

/**
 * Function spaces for the Posisson-Bolzmann equations
 */
class PbeqSpace
{
public:

    /** @name Typedefs
     */
    //@{

    static const uint16_type Dim = 3;
    static const uint16_type Order   = 2;
    static const uint16_type OrderHeavyside = 1;
    static const uint16_type imOrder = 2*Order + OrderHeavyside;
    static const uint16_type imOrderHeavyside = OrderHeavyside;


    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<Dim, 1,Dim> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /*basis*/
    typedef fusion::vector<fem::Lagrange<Dim, OrderHeavyside, Scalar, Discontinuous, value_type, Simplex> >  heavyside_basis_type;
    typedef fusion::vector<fem::Lagrange<Dim, Order, Scalar, Continuous, value_type, Simplex> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, heavyside_basis_type, value_type> heavyside_space_type;
    typedef boost::shared_ptr<heavyside_space_type>                    heavyside_space_ptrtype;
    // element_type is the type of an element of $V_h$
    typedef heavyside_space_type::element_type                heavyside_element_type;

    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type>                    space_ptrtype;
    // element_type is the type of an element of $V_h$
    typedef space_type::element_type                element_type;

    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Simplex> im_type;
    typedef im_type::nodes_type nodes_type;

    typedef IM<Dim, imOrderHeavyside, value_type, Simplex> heavyside_im_type;

    //typedef ublas::matrix<value_type, ublas::column_major> nodes_type;

    /*function type, see also file operations.hpp::403 */
    /*
      ublas::vector<value_type> f( nodes_type const& ) );
    */

    /* derived from molecule */
    typedef Molecule molecule_type;
    typedef molecule_type::node_type node_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{

    PbeqSpace( );

    PbeqSpace( value_type const meshSize,
               value_type const farFactor,
               value_type const farBnd,
               node_type& stretch,
               node_type& translation );

    PbeqSpace( value_type const meshSize,
               value_type const farFactor,
               value_type const farBnd );

    PbeqSpace( PbeqSpace const& tc );


    // ~PbeqSpace() {}

    //@}

    /** @name Operator overloads
     */
    //@{
    //@}

    void setUpInvStretch();

    heavyside_element_type heavyside( molecule_type const& molecule ) const;

    heavyside_element_type fastHeavyside( molecule_type const& molecule ) const;

    void intvrho( molecule_type const& molecule, vector_ptrtype rhs ) const;

    heavyside_element_type chargeDensity( molecule_type const& molecule ) const;

    mesh_ptr_type mesh()
    {
        return M_mesh;
    }

    space_ptrtype Xh()
    {
        return M_Space;   // Space for the FE function
    }

    heavyside_space_ptrtype HSpace()
    {
        return M_HeavysideSpace;   // Space for the Heavyside function
    }


    /**
     * load the mesh  as well as create the finite element spaces
     */
    bool loadMesh( std::string const& meshname = entity_type::name().append( ".msh" ) );

    /**
     * create the mesh using mesh size \c meshSize as well as the finite element spaces
     */
    void createMesh( bool const geoOnly=false );

    /** @name setters and getters
     */
    //@{
    void setMeshSize( value_type const meshSize )
    {
        M_meshSize = meshSize;
    }
    void setFarFactor( value_type const farFactor )
    {
        M_farFactor = farFactor;
    }
    void setFarBnd( value_type const farBnd )
    {
        M_farBnd = farBnd;
    }
    void setStretch    ( node_type const& stretch )
    {
        M_stretch = stretch;
        setUpInvStretch();
    }
    void setTranslation( node_type const& translation )
    {
        M_translation = translation;
    }
    void setSmoothWindow( value_type const& sW )
    {
        M_sW = sW;
    }


    value_type meshSize() const
    {
        return M_meshSize;
    }
    value_type farFactor() const
    {
        return M_farFactor;
    }
    value_type farBnd() const
    {
        return M_farBnd;
    }

    node_type const& stretch()     const
    {
        return M_stretch;
    }
    node_type const& translation() const
    {
        return M_translation;
    }

    //@}

private:

    /**
     * create the finite Element spaces
     */
    void createSpaces( );

    /**
     * returns a finite element representation of the heavyside function
     * centered in center with radius radius
     */
    heavyside_element_type heavyside( value_type const& radius,
                                      node_type  const& center ) const;

    heavyside_element_type chargeDensity( value_type const& radius,
                                          value_type const& charge,
                                          node_type  const& center ) const;

private:

    value_type M_meshSize;
    value_type M_farFactor;
    value_type M_farBnd;
    node_type  M_stretch, M_invStretch, M_translation;

    value_type M_sW;

    bool M_meshSetted;
    bool M_spacesSetted;

    mesh_ptr_type M_mesh;

    heavyside_space_ptrtype M_HeavysideSpace; // Space for the Heavyside function
    space_ptrtype           M_Space;          // Space for the FE function

    std::map<std::string,std::pair<boost::timer,value_type> > timers;
    std::map<std::string,value_type> stats;
}; // PbeqSpace





} // end namespace Feel

#endif
