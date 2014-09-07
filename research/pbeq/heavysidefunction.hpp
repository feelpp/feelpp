/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-08-24

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
   \file heavysidefunction.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-08-24
 */

#ifndef HEAVYSIDEFUNCTION_H
#define HEAVYSIDEFUNCTION_H

#include "pbeqspace.hpp"

namespace Feel
{

class heavysideFunction
{
public:

    typedef PbeqSpace                            pbeqspace_type;
    typedef pbeqspace_type::nodes_type           nodes_type;
    typedef pbeqspace_type::molecule_type        molecule_type;
    typedef pbeqspace_type::heavyside_space_type space_type;
    typedef pbeqspace_type::mesh_type            mesh_type;
    typedef pbeqspace_type::node_type            node_type;

    static const uint16_type Dim;

    typedef double value_type;

    // element_type is the type of an element of $V_h$
    typedef mesh_type::element_type          element_type;

    typedef element_type::gm_type             gm_type;
    typedef element_type::gm_ptrtype          gm_ptrtype;

    heavysideFunction() : M_elt( 0 ), M_molecule( 0 ), M_stretch( 0 ), M_translation( 0 ) {}
    ~heavysideFunction()
    {
        M_elt=0;
        M_molecule=0;
        M_stretch=0;
        M_translation=0;
    }

    ublas::vector<double> operator()( nodes_type const& pointsOnRef ) const;

    value_type operator()( node_type const& pointOnRef ) const;

    void setSmoothWindow( value_type const& sW );

    void setMolecule( molecule_type const* molecule )
    {
        M_molecule = molecule;
    }
    void unSetMolecule( )
    {
        M_molecule = 0;
    }

    void setElement( element_type const* elt )
    {
        M_elt = elt;
    }
    void unSetElement( )
    {
        M_elt = 0;
    }

    void setStretch( node_type const*  stretch )
    {
        M_stretch = stretch;
    }
    void unSetStretch( )
    {
        M_stretch = 0;
    }

    void setTranslation( node_type const*  translation )
    {
        M_translation = translation;
    }
    void unSetTranslation( )
    {
        M_translation = 0;
    }

private:

    nodes_type transformToReal( nodes_type const& Gt ) const;


private:
    element_type const* M_elt;
    molecule_type const* M_molecule;
    node_type const*     M_stretch;
    node_type const*     M_translation;
    value_type           M_sW;
    value_type           M_sW2;
    value_type           M_2sW;
    value_type           M_4sW3;
    value_type           M_34sW2;

};

} // namespace Feel

#endif /* __HEAVYSIDEFUNCTION_H */
