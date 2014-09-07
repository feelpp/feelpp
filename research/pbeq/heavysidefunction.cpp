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
   \file heavysidefunction.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-08-24
 */


#include "heavysidefunction.hpp"

namespace Feel
{

void
heavysideFunction::setSmoothWindow( value_type const& sW )
{
    M_sW = sW;
    M_sW2 = sW*sW;
    M_2sW = 2*sW;
    M_4sW3 = 1./( 4*sW*sW*sW );
    M_34sW2 = 3./( 4*sW*sW );
}

ublas::vector<value_type>
heavysideFunction::operator()( nodes_type const& pointsOnRef ) const
{
    ublas::vector<value_type> result( pointsOnRef.size2() );
    std::fill( result.begin(), result.end(), 1.0 );

    nodes_type pointsHat( transformToReal( pointsOnRef ) );

    molecule_type::atoms_const_iterator_type atom( M_molecule->begin() );
    int i;

    if ( atom == M_molecule->end() ) return result;

    node_type Ellipse( pbeqspace_type::Dim );
    node_type point( pbeqspace_type::Dim );

    value_type r, dr, dr2;


    for ( i = 0; i < pointsOnRef.size2(); ++i )
    {
        point = element_prod( *M_stretch,column( pointsHat,i ) ) + ( *M_translation );

        node_type mine( column( pointsHat,i ) );
        //std::cout << "mine = " << mine(0)  << " "<< mine(1) << " " << mine(2) << std::endl;
        //std::cout << "point = " << point(0)  << " "<< point(1) << " " << point(2) << std::endl;

        for ( atom = M_molecule->begin(); atom != M_molecule->end(); ++atom )
        {
            Ellipse =  point  - atom->center();

            r = norm_2( Ellipse );

            if  ( r <= atom->radius() - M_sW )
            {
                result( i ) = 0;
                break;
            }

            if  ( M_sW == 0 || r >= atom->radius() + M_sW )
                continue;

            dr  = r  - atom->radius() + M_sW;
            dr2 = dr*dr;

            result( i ) *= - M_4sW3 * ( dr*dr2 ) + M_34sW2 * dr2;

        }

        //std::cout << "result(" << i << ") = " << result(i) << std::endl;
    }


    return result;
}


heavysideFunction::nodes_type
heavysideFunction::transformToReal( nodes_type const& Gt ) const
{
    gm_ptrtype _gm_ptr( new gm_type );

    gm_type::precompute_ptrtype __geopc( new gm_type::precompute_type( _gm_ptr, Gt ) );

    gm_type::Context<vm::POINT, element_type> gmc( _gm_ptr, *M_elt, __geopc );

    return gmc.xReal();
}


heavysideFunction::value_type
heavysideFunction::operator()( node_type const& pointHat ) const
{
    node_type Ellipse( pbeqspace_type::Dim );
    node_type point( pbeqspace_type::Dim );

    point = element_prod( *M_stretch,pointHat ) + ( *M_translation );

    molecule_type::atoms_const_iterator_type atom( M_molecule->begin() );

    for ( ; atom != M_molecule->end(); ++atom )
    {
        Ellipse =  point  - atom->center();

        //if ( norm_inf(Ellipse) >  atom->radius() ) continue;

        if  ( norm_2( Ellipse ) < atom->radius2() )
        {
            return 0;
        }
    }

    return 1;
}

} // namespace Feel
