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
//! @file distancepointtoface.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _DISTANCE_POINT_TO_FACE_HPP
#define _DISTANCE_POINT_TO_FACE_HPP 1

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "eigenmap.hpp"

namespace Feel {

namespace detail {
namespace geometry {

template< typename Derived1, typename Derived2, typename Derived3 >
auto 
distancePointToSegment( 
        // Point vertex
        Eigen::MatrixBase<Derived1> const& pt, 
        // Segment vertices
        Eigen::MatrixBase<Derived2> const& p1,
        Eigen::MatrixBase<Derived3> const& p2
        );

template< typename Derived1, typename Derived2 >
auto 
distancePointToSegment( 
        Eigen::MatrixBase<Derived1> const& pt, 
        Eigen::MatrixBase<Derived2> const& seg 
        );

template< typename Derived1, typename Derived2, typename Derived3, typename Derived4 >
auto 
distancePointToTriangle( 
        // Point vertex
        Eigen::MatrixBase<Derived1> const& pt, 
        // Triangles vertices
        Eigen::MatrixBase<Derived2> const& p1,
        Eigen::MatrixBase<Derived3> const& p2,
        Eigen::MatrixBase<Derived4> const& p3
        );

template< typename Derived1, typename Derived2 >
auto 
distancePointToTriangle( 
        Eigen::MatrixBase<Derived1> const& pt, 
        Eigen::MatrixBase<Derived2> const& tri 
        );

template< uint16_t RealDim, typename PointType, typename MatrixType >
auto
distancePointToFace( PointType const& pt, MatrixType const& face )
{
    static_assert( RealDim == 2 || RealDim == 3, "Dimension must be 2 or 3" );

    auto P = eigenMap<RealDim>( pt );
    auto F = eigenMap<RealDim>( face );

    if constexpr ( RealDim == 2 )
        return distancePointToSegment( P, F );
    else if constexpr ( RealDim == 3 )
        return distancePointToTriangle( P, F );

}

//--------------------------------------------------------------------//

template< typename Derived1, typename Derived2, typename Derived3 >
auto 
distancePointToSegment( 
        // Point vertex
        Eigen::MatrixBase<Derived1> const& P, 
        // Segment vertices
        Eigen::MatrixBase<Derived2> const& P1,
        Eigen::MatrixBase<Derived3> const& P2
        )
{
    typedef typename Derived1::Scalar value_type;

    auto const P1P2 = P2 - P1;
    auto const P1P = P - P1;

    value_type r1 = P1P.dot(P1P2);
    if( r1 <= 0. )
    {
        return P1P.norm();
    }
    value_type r2 = P1P2.squaredNorm();
    if( r2 <= r1 )
    {
        return (P-P2).norm();
    }

    auto Proj = P1 + (r1/r2)*P1P2;
    return (P-Proj).norm();
} 

template< typename Derived1, typename Derived2 >
auto 
distancePointToSegment( 
        Eigen::MatrixBase<Derived1> const& pt, 
        Eigen::MatrixBase<Derived2> const& seg 
        )
{
    return distancePointToSegment( pt, seg.col(0), seg.col(1) );
}


template< typename Derived1, typename Derived2, typename Derived3 >
auto
triangleFrameTransformationMatrix(
        Eigen::MatrixBase<Derived1> const& P1,
        Eigen::MatrixBase<Derived2> const& P2,
        Eigen::MatrixBase<Derived3> const& P3
        )
{
    auto const P1P2 = P2 - P1;
    auto const P1P3 = P3 - P1;

    typedef typename Derived1::Scalar value_type;

    // Change frame to (P1P2, P1P3, unit(P1P2xP1P3))
    Eigen::Matrix<value_type, 3, 3> R;
    R.col(0) = P1P2;
    R.col(1) = P1P3;
    R.col(2) = (R.col(0).cross( R.col(1) )).normalized();
    return R.inverse().eval();
}

template< typename Derived0, typename Derived1, typename Derived2, typename Derived3, typename Derived4 >
auto 
distancePointToTriangleWithTransformationMatrix( 
        // Point vertex
        Eigen::MatrixBase<Derived0> const& P, 
        // Triangles vertices
        Eigen::MatrixBase<Derived1> const& P1,
        Eigen::MatrixBase<Derived2> const& P2,
        Eigen::MatrixBase<Derived3> const& P3,
        // Triangle frame transformation matrix
        Eigen::MatrixBase<Derived4> const& Rinv
        )
{
    Eigen::Matrix<typename Derived0::Scalar, 3, 1> Ptild = Rinv * (P - P1);

    if( Ptild(0) < 0. )
    {
        // Xtild < 0: dist = dist(P,[P1P3])
        return distancePointToSegment( P, P1, P3 );
    }
    if( Ptild(1) < 0. )
    {
        // Ytild < 0: dist = dist(P,[P1P2])
        return distancePointToSegment( P, P1, P2 );
    }
    if( 1. < Ptild(0)+Ptild(1) )
    {
        // 1-Xtild < Ytild: dist = dist(P,[P2P3])
        return distancePointToSegment( P, P2, P3 );
    }
    // Otherwise Ptild projects in Tri: dist = |Ptild_z|
    return std::abs( Ptild(2) );
}

template< typename Derived0, typename Derived1, typename Derived2, typename Derived3 >
auto 
distancePointToTriangle( 
        // Point vertex
        Eigen::MatrixBase<Derived0> const& P, 
        // Triangles vertices
        Eigen::MatrixBase<Derived1> const& P1,
        Eigen::MatrixBase<Derived2> const& P2,
        Eigen::MatrixBase<Derived3> const& P3
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived0, 3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1, 3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2, 3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived3, 3)

    //auto const P1P2 = P2 - P1;
    //auto const P1P3 = P3 - P1;

    //typedef typename Derived1::Scalar value_type;

    //// Change frame to (P1P2, P1P3, unit(P1P2xP1P3))
    //Eigen::Matrix<value_type, 3, 3> R;
    //R.col(0) = P1P2;
    //R.col(1) = P1P3;
    //R.col(2) = (R.col(0).cross( R.col(1) )).normalized();
    //auto Rinv = R.inverse();
    auto Rinv = triangleFrameTransformationMatrix( P1, P2, P3 );
    return distancePointToTriangleWithTransformationMatrix( P, P1, P2, P3, Rinv );
} 

template< typename Derived1, typename Derived2 >
auto 
distancePointToTriangle( 
        Eigen::MatrixBase<Derived1> const& pt, 
        Eigen::MatrixBase<Derived2> const& tri 
        )
{
    return distancePointToTriangle( pt, tri.col(0), tri.col(1), tri.col(2) );
}

} // namespace geometry
} // namespace detail
} // namespace Feel

#endif // _DISTANCE_POINT_TO_FACE_HPP
