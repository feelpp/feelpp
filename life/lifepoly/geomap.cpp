/*
 This file is part of the Life library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#include <cmath>

#include <boost/lambda/lambda.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>

#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/blas/blas.hpp>

#include <life/lifealg/glas.hpp>
#include <life/lifemesh/geoentity.hpp>

#include <life/lifepoly/geomap.hpp>
#include <life/lifemesh/geond.hpp>
#include <life/lifemesh/geo0d.hpp>



namespace Life
{

template<typename Elem, typename T>
RealToReference<Elem,T>::RealToReference( Elem const& elem )
    :
    _M_gm( new gm_type ),
    _M_igm( _M_gm, elem )
{}


template<typename Elem, typename T>
typename RealToReference<Elem,T>::points_type
RealToReference<Elem,T>::operator()( points_type const& pts ) const
{
    return _M_igm( pts );
}

template class RealToReference< GeoND<2,GeoEntity<Simplex<2,1> >, double>, double >;
template class RealToReference< GeoND<3,GeoEntity<Simplex<3,1> >, double>, double >;

#if defined( LIFE_HAVE_QD_REAL)
template class RealToReference< GeoND<2,GeoEntity<Simplex<2,1> >, qd_real>, qd_real >;
template class RealToReference< GeoND<2,GeoEntity<Simplex<2,1> >, dd_real>, dd_real >;
#endif // LIFE_HAVE_QD_REAL

} // Life
