/*
This file is part of the Feel library
Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#include <feel/feelmesh/bareitems.hpp>

namespace Feel
{
/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from three Point size_type's
  \return  bool is false if orientation  has been changed.
  \param i is a Point size_type
  \param j is a Point size_type
  \param k is a Point size_type

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

*/
std::pair<BareFace, bool>
makeBareFace( size_type const i, size_type const j, size_type const k )
{
    if ( i < j && i < k )
    {
        if ( j < k )
        {
            return std::make_pair( BareFace( i, j, k ), true );
        }

        else
        {
            return std::make_pair( BareFace( i, k, j ), false );
        }
    }

    else if ( j < k && j < i )
    {
        if ( k < i )
        {
            return std::make_pair( BareFace( j, k, i ), true );
        }

        else
        {
            return std::make_pair( BareFace( j, i, k ), false );
        }
    }

    else
    {
        if ( i < j )
        {
            return std::make_pair( BareFace( k, i, j ), true );
        }

        else
        {
            return std::make_pair( BareFace( k, j, i ), false );
        }
    }
}

/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from four Point size_type's
  \return is false if orientation  has been changed.
  \param i is a Point size_type
  \param j is a Point size_type
  \param k is a Point size_type
  \param l is a Point size_type
  \pre i, j, k and l >0. i!=j!=k!=l

  To be used with Quad faces.

\remarks For quad faces the construction process is more complex. We start from
  the smallest vertex and we take the first three vertices in the
   sequence. We then procede as for the triangles.
*/

std::pair<BareFace, bool>
makeBareFace( size_type const i, size_type const j, size_type const k, size_type const l )
{
    std::vector<size_type> helper( 4 );
    helper[ 0 ] = i;
    helper[ 1 ] = j;
    helper[ 2 ] = k;
    helper[ 3 ] = l;
    std::vector<size_type>::iterator vi = std::max_element( helper.begin(), helper.end() );
    std::rotate( helper.begin(), vi, helper.end() );
    return makeBareFace( helper[ 1 ], helper[ 2 ], helper[ 3 ] );
}
}
