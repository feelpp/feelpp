/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
      Date: 2005-08-10

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier Grenoble 1

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
  \file refentity.hpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2005-08-10
*/
#ifndef __RefEntity_H
#define __RefEntity_H 1

#include <stdexcept>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>

namespace Feel
{
template<size_type ShapeE, typename T = double>
class Entity
{
};

/**
 * \class Reference
 * \brief Reference convex
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename Geo, uint16_type Dim = 1, uint16_type Order = 1, uint16_type RDim = Dim, typename T = double>
class Reference
{};

template<typename Geo, uint16_type Dim, uint16_type Order, uint16_type RDim, typename T>
std::ostream&
operator<<( std::ostream& os,
            Reference<Geo, Dim, Order, RDim, T> const& ref )
{
    typedef Reference<Geo, Dim, Order, RDim, T> ref_type;
    os << "     Dimension: " << ref_type::nDim << "\n"
       << "         Order: " << ref_type::nOrder << "\n"
       << "Real dimension: " << ref_type::nRealDim << "\n";
    os << " Vertices: " << ref.vertices() << "\n";
    os << "  Normals: " << ref.normals() << "\n";
    return os;
}
template<typename RefEntity>
void toPython( RefEntity const& e, std::string str = "simplex" )
{
    typedef typename RefEntity::value_type value_type;
    typedef typename RefEntity::node_type node_type;
    std::ostringstream ostr;
    ostr << str
         << "_" << RefEntity::nDim
         << "_" << RefEntity::nOrder
         << "_" << RefEntity::nRealDim
         << ".py";
    std::ofstream ofs( ostr.str().c_str() );

    ofs << "from pyx import *\n";
    ofs << "p=path.path(";

    for ( int i = 0; i < RefEntity::numEdges; ++i )
    {
        for ( int j = 0; j < 2; ++j )
        {
            node_type x( 2 );

            if ( RefEntity::nRealDim == 1 )
            {
                x( 0 ) = e.edgeVertex( i,j )( 0 );
                x( 1 ) = value_type( 0 );
            }

            if ( RefEntity::nRealDim == 2 )
            {
                x = e.edgeVertex( i, j );
            }

            if ( RefEntity::nRealDim == 3 )
            {
                x( 0 ) = e.edgeVertex( i, j )( 0 )+e.edgeVertex( i, j )( 1 )*std::cos( M_PI/4 );
                x( 1 ) = e.edgeVertex( i, j )( 2 )+e.edgeVertex( i, j )( 1 )*std::sin( M_PI/4 );
            }

            if ( j == 0 )
                ofs << "path.moveto(" << double( x( 0 ) )<< "," << double( x( 1 ) ) << "),\n";

            else if ( j == 1 )
                ofs << "path.lineto(" << double( x( 0 ) )<< "," << double( x( 1 ) ) << "),\n";
        }
    }

    ofs << "path.closepath() )\n";
    ofs << "c = canvas.canvas()\n"
        << "c.stroke(p, [style.linewidth.Thin])\n";

    for ( int i = 0; i < RefEntity::numPoints; ++i )
    {
        node_type x( 2 );

        if ( RefEntity::nRealDim == 1 )
        {
            x( 0 ) = e.point( i )( 0 );
            x( 1 ) = value_type( 0 );
        }

        if ( RefEntity::nRealDim == 2 )
        {
            x = e.point( i );
        }

        if ( RefEntity::nRealDim == 3 )
        {
            x( 0 ) = e.point( i )( 0 )+e.point( i )( 1 )*std::cos( M_PI/4 );
            x( 1 ) = e.point( i )( 2 )+e.point( i )( 1 )*std::sin( M_PI/4 );
        }

        ofs << "c.fill ( path.circle(" << double( x( 0 ) ) << "," << double( x( 1 ) )<< ", 0.05 ),[deco.filled([color.grey.black])])\n";
        ofs << "c.text(" << double( x( 0 ) ) << "," << double( x( 1 ) ) << ", \"" << i << "\")\n";
    }

    ofs << "c.writePDFfile(\"" << str << "_" << RefEntity::nDim
        << "_" << RefEntity::nOrder
        << "_" << RefEntity::nRealDim
        << "\", paperformat=\"a4\")\n";
}

} // Feel

#include <feel/feelmesh/refsimplex.hpp>
#include <feel/feelmesh/refhypercube.hpp>
#endif /* __RefEntity_H */
