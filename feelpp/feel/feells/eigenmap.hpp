//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Author(s) : 
//!     Thibaut Metivet <thibaut.metivet@inria.fr>
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file eigenmap.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _EIGEN_MAP_HPP
#define _EIGEN_MAP_HPP 1

#include <Eigen/Core>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

namespace Feel {

template< int Rows, typename T >
auto
eigenMap( ublas::vector<T> const& v, int nrows = Rows)
{
    typedef Eigen::Matrix<T, Rows, 1, Eigen::ColMajor> EigenVectorType;
    return Eigen::Map< EigenVectorType >(
            const_cast<T*>( v.data().begin() ), nrows, 1
            );
}

template< int Rows, int Cols = Eigen::Dynamic, typename T >
auto
eigenMap( ublas::matrix<T, ublas::column_major> const& m, int nRows = Rows, int cols = Cols )
{
    typedef Eigen::Matrix<T, Rows, Cols, Eigen::ColMajor> EigenMatrixType;
    int nCols = cols;
    if( nCols == Eigen::Dynamic )
        nCols = m.size2();
    return Eigen::Map< EigenMatrixType >( 
            const_cast<T*>( m.data().begin() ), nRows, nCols
            );
}

} // namespace Feel

#endif // _EIGEN_MAP_HPP
