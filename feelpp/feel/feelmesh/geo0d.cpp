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
#include <feel/feelmesh/geo0d.hpp>

namespace Feel
{
/*--------------------------------------------------------------
                                 Geo0D
  ---------------------------------------------------------------*/
template<uint16_type Dim>
Geo0D<Dim>::Geo0D()
    :
    MeshEntityWithBoundary( 0 ),
    _coor( Dim )
{
    _coor = ZeroVector( _coor.size() );
}

template<uint16_type Dim>
Geo0D<Dim>::Geo0D( uint16_type id, bool boundary )
    :
    MeshEntityWithBoundary( id, boundary ),
    _coor( Dim )
{
    _coor = ZeroVector( _coor.size() );
}

template<uint16_type Dim>
Geo0D<Dim>::Geo0D( uint16_type id, Real x, Real y, Real z, bool boundary )
    :
    MeshEntityWithBoundary( id, boundary ),
    _coor()
{
    if ( Dim < 2 )
        _coor[ 0 ] = x;

    if (  Dim < 3 )
        _coor[ 1 ] = y;

    if ( Dim == 3 )
        _coor[ 2 ] = z;
}

template<uint16_type Dim>
Geo0D<Dim>::Geo0D( Geo0D const & G )
    :
    MeshEntityWithBoundary( G._id, G._boundary ),
    _coor( G._coor )
{
}

template<uint16_type Dim>
Geo0D<Dim> &
Geo0D<Dim>::operator=( Geo0D<Dim> const & G )
{
    if (  this == &G )
        return *this;

    _id = G._id;
    _boundary = G._boundary;
    _coor = G._coor;
    return *this;
}

template<uint16_type Dim>
std::ostream &
Geo0D<Dim>::showMe( bool verbose, std::ostream & out ) const
{
    out.setf( std::ios::scientific, std::ios::floatfield );
    out << " Geo0D object " << std::endl;

    if ( verbose )
    {
        unsigned i;
        out << " node:" << std::endl;
        out << node() << "\n";
    }

    out << "id = " << id() << "  ";
    out << "----- END OF Geo0D data ---" << std::endl << std::endl;
    return out;
}

}
