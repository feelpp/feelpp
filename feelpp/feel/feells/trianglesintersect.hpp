/*
 * 2D and 3D triangle-triangle intersection test.
 *
 * This file contains a C++ "Eigen-based" implementation 
 * of the triangle-triangle overlap test described in
 *
 * "Fast and Robust Triangle-Triangle Overlap Test
 *  Using Orientation Predicates"  P. Guigue - O. Devillers
 *
 *  Journal of Graphics Tools, 8(1), 2003
 *
 * Author: Thibaut Metivet <thibaut.metivet@inria.fr>
 * Copyright (C) 2019 Thibaut Metivet
 *
 */

#ifndef _TRIANGLES_INTERSECT_HPP
#define _TRIANGLES_INTERSECT_HPP 1

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Feel {
namespace detail {
namespace geometry {

template< typename Derived1, typename Derived2 >
bool triangles2DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        );

template< typename Derived1, typename Derived2 >
bool triangles3DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        );

//--------------------------------------------------------------------//
// 2D intersection test

template< typename Derived >
auto counterClockwiseTest( 
        Eigen::MatrixBase<Derived> const& a,
        Eigen::MatrixBase<Derived> const& b,
        Eigen::MatrixBase<Derived> const& c
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived,2)

    return (a(0)-c(0))*(b(1)-c(1))-(a(1)-c(1))*(b(0)-c(0));
}

template< typename Derived1, typename Derived2 >
bool vertex2DIntersectionTest( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,2)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,2)

    if (counterClockwiseTest(r2,p2,q1) >= 0.)
    {
        if (counterClockwiseTest(r2,q2,q1) <= 0.)
        {
            if (counterClockwiseTest(p1,p2,q1) > 0.) 
            {
                if (counterClockwiseTest(p1,q2,q1) <= 0.) 
                    return true; 
                else 
                    return false;
            } 
            else 
            {
                if (counterClockwiseTest(p1,p2,r1) >= 0.)
                    if (counterClockwiseTest(q1,r1,p2) >= 0.) 
                        return true; 
                    else 
                        return false;
                else 
                    return false;
            }
        }
        else 
        {
            if (counterClockwiseTest(p1,q2,q1) <= 0.)
                if (counterClockwiseTest(r2,q2,r1) <= 0.)
                    if (counterClockwiseTest(q1,r1,q2) >= 0.) 
                        return true; 
                    else 
                        return false;
                else 
                    return false;
            else 
                return false;
        }
    }
    else
    {
        if (counterClockwiseTest(r2,p2,r1) >= 0.) 
            if (counterClockwiseTest(q1,r1,r2) >= 0.)
                if (counterClockwiseTest(p1,p2,r1) >= 0.) 
                    return true;
                else return false;
            else 
                if (counterClockwiseTest(q1,r1,q2) >= 0.) 
                {
                    if (counterClockwiseTest(r2,r1,q2) >= 0.) 
                        return true; 
                    else 
                        return false; 
                }
                else 
                    return false;
        else 
            return false;
    }
}

template< typename Derived1, typename Derived2 >
bool edge2DIntersectionTest( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,2)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,2)

    if (counterClockwiseTest(r2,p2,q1) >= 0.) 
    {
        if (counterClockwiseTest(p1,p2,q1) >= 0.) 
        { 
            if (counterClockwiseTest(p1,q1,r2) >= 0.) 
                return true; 
            else 
                return false;
        } 
        else 
        { 
            if (counterClockwiseTest(q1,r1,p2) >= 0.)
            { 
                if (counterClockwiseTest(r1,p1,p2) >= 0.) 
                    return true; 
                else 
                    return false;
            } 
            else 
                return false; 
        } 
    } 
    else 
    {
        if (counterClockwiseTest(r2,p2,r1) >= 0.) 
        {
            if (counterClockwiseTest(p1,p2,r1) >= 0.) 
            {
                if (counterClockwiseTest(p1,r1,r2) >= 0.) 
                    return true;  
                else 
                {
                    if (counterClockwiseTest(q1,r1,r2) >= 0.) 
                        return true; 
                    else 
                        return false;
                }
            }
            else  
                return false; 
        }
        else 
            return false; 
    }
}

template< typename Derived1, typename Derived2 >
bool ccwTriangles2DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,2)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,2)

    if ( counterClockwiseTest(p2,q2,p1) >= 0. ) 
    {
        if ( counterClockwiseTest(q2,r2,p1) >= 0. ) 
        {
            if ( counterClockwiseTest(r2,p2,p1) >= 0. ) 
                return true;
            else 
                return edge2DIntersectionTest(p1,q1,r1,p2,q2,r2);
        } 
        else 
        {  
            if ( counterClockwiseTest(r2,p2,p1) >= 0. ) 
                return edge2DIntersectionTest(p1,q1,r1,r2,p2,q2);
            else 
                return vertex2DIntersectionTest(p1,q1,r1,p2,q2,r2);
        }
    }
    else 
    {
        if ( counterClockwiseTest(q2,r2,p1) >= 0. ) 
        {
            if ( counterClockwiseTest(r2,p2,p1) >= 0. ) 
                return edge2DIntersectionTest(p1,q1,r1,q2,r2,p2);
            else  
                return vertex2DIntersectionTest(p1,q1,r1,q2,r2,p2);
        }
        else 
            return vertex2DIntersectionTest(p1,q1,r1,r2,p2,q2);
    }
}

template< typename Derived1, typename Derived2 >
bool triangles2DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,2)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,2)

    if ( counterClockwiseTest(p1,q1,r1) < 0. )
        if ( counterClockwiseTest(p2,q2,r2) < 0. )
            return ccwTriangles2DIntersect(p1,r1,q1,p2,r2,q2);
        else
            return ccwTriangles2DIntersect(p1,r1,q1,p2,q2,r2);
    else
        if ( counterClockwiseTest(p2,q2,r2) < 0. )
            return ccwTriangles2DIntersect(p1,q1,r1,p2,r2,q2);
        else
            return ccwTriangles2DIntersect(p1,q1,r1,p2,q2,r2);
}

//--------------------------------------------------------------------//
// 3D intersection test

template< typename Derived1, typename Derived2 >
bool checkMinMaxInterval3DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,3)

    auto const p2mq1 = p2 - q1;
    auto const p1mq1 = p1 - q1;
    auto const N1 = p2mq1.cross(p1mq1);
    auto const q2mq1 = q2 - q1;
    if( q2mq1.dot(N1) > 0. )
        return false;

    auto const p2mp1 = p2 - p1;
    auto const r1mp1 = r1 - p1;
    auto const N2 = p2mp1.cross(r1mp1);
    auto const r2mp1 = r2 - p1;
    if( r2mp1.dot(N2) > 0. )
        return false;
    else
        return true;
}

template< typename Derived1, typename Derived2, typename Derived3 >
bool coplanarTriangles3DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2,
        // Normal
        Eigen::MatrixBase<Derived3> const& N
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived3,3)

    typedef Eigen::Matrix<typename Derived1::Scalar, 2, 1> vertex1_type;
    typedef Eigen::Matrix<typename Derived2::Scalar, 2, 1> vertex2_type;

    vertex1_type P1, Q1, R1;
    vertex2_type P2, Q2, R2;
    Eigen::PermutationMatrix<3> permutationMatrix;

    auto const Nabs = N.cwiseAbs().eval();

    // Project the 3D triangles onto 2D such that the area of the projection is maximized.
    if (( Nabs(0) > Nabs(2) ) && ( Nabs(0) >= Nabs(1) )) 
    {
        // Project onto plane YZ
        permutationMatrix.indices() << 2, 1, 0;

        P1 = ( permutationMatrix * q1 ).template head<2>();
        Q1 = ( permutationMatrix * p1 ).template head<2>();
        R1 = ( permutationMatrix * r1 ).template head<2>();

        P2 = ( permutationMatrix * q2 ).template head<2>();
        Q2 = ( permutationMatrix * p2 ).template head<2>();
        R2 = ( permutationMatrix * r2 ).template head<2>();
    } 
    else if (( Nabs(1) > Nabs(2) ) && ( Nabs(1) >= Nabs(0) )) 
    {
        // Project onto plane XZ
        permutationMatrix.indices() << 0, 2, 1;

        P1 = ( permutationMatrix * q1 ).template head<2>();
        Q1 = ( permutationMatrix * p1 ).template head<2>();
        R1 = ( permutationMatrix * r1 ).template head<2>();

        P2 = ( permutationMatrix * q2 ).template head<2>();
        Q2 = ( permutationMatrix * p2 ).template head<2>();
        R2 = ( permutationMatrix * r2 ).template head<2>();
    } 
    else 
    {
        // Project onto plane XY
        P1 = p1.template head<2>();
        Q1 = q1.template head<2>();
        R1 = r1.template head<2>();

        P2 = p2.template head<2>();
        Q2 = q2.template head<2>();
        R2 = r2.template head<2>();
    }

    return triangles2DIntersect(P1,Q1,R1,P2,Q2,R2);
}

template< typename Derived1, typename Derived2, typename Derived3 >
bool triangles3DIntersect_testTriangle2VerticesPermutations( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2,
        // Triangle 1 normal
        Eigen::MatrixBase<Derived3> const& N1,
        // Triangle 2 vertices distance sign to triangle 1 plane
        typename Derived2::Scalar dp2,
        typename Derived2::Scalar dq2,
        typename Derived2::Scalar dr2
        )
{
    if (dp2 > 0.) 
    { 
        if (dq2 > 0.) 
            return checkMinMaxInterval3DIntersect(p1,r1,q1,r2,p2,q2);
        else if (dr2 > 0.) 
            return checkMinMaxInterval3DIntersect(p1,r1,q1,q2,r2,p2);
        else 
            return checkMinMaxInterval3DIntersect(p1,q1,r1,p2,q2,r2);
    }
    else if (dp2 < 0.) 
    { 
        if (dq2 < 0.) 
            return checkMinMaxInterval3DIntersect(p1,q1,r1,r2,p2,q2);
        else if (dr2 < 0.) 
            return checkMinMaxInterval3DIntersect(p1,q1,r1,q2,r2,p2);
        else 
            return checkMinMaxInterval3DIntersect(p1,r1,q1,p2,q2,r2);
    } 
    else 
    { 
        if (dq2 < 0.) 
        { 
            if (dr2 >= 0.)  
                return checkMinMaxInterval3DIntersect(p1,r1,q1,q2,r2,p2);
            else 
                return checkMinMaxInterval3DIntersect(p1,q1,r1,p2,q2,r2);
        } 
        else if (dq2 > 0.) 
        { 
            if (dr2 > 0.) 
                return checkMinMaxInterval3DIntersect(p1,r1,q1,p2,q2,r2);
            else  
                return checkMinMaxInterval3DIntersect(p1,q1,r1,q2,r2,p2);
        } 
        else  
        { 
            if (dr2 > 0.) 
                return checkMinMaxInterval3DIntersect(p1,q1,r1,r2,p2,q2);
            else if (dr2 < 0.) 
                return checkMinMaxInterval3DIntersect(p1,r1,q1,r2,p2,q2);
            else 
                return coplanarTriangles3DIntersect(p1,q1,r1,p2,q2,r2,N1);
        }
    }
}

template< typename Derived1, typename Derived2, typename Derived3 >
bool triangles3DIntersect_testTrianglesVerticesPermutations( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2,
        // Triangle 1 normal
        Eigen::MatrixBase<Derived3> const& N1,
        // Triangle 2 normal
        Eigen::MatrixBase<Derived3> const& N2,
        // Triangle 1 vertices distance sign to triangle 2 plane
        Eigen::MatrixBase<Derived3> const& d1,
        // Triangle 2 vertices distance sign to triangle 1 plane
        Eigen::MatrixBase<Derived3> const& d2
        )
{
    typename Derived3::Scalar dp1 = d1(0);
    typename Derived3::Scalar dq1 = d1(1);
    typename Derived3::Scalar dr1 = d1(2);

    typename Derived3::Scalar dp2 = d2(0);
    typename Derived3::Scalar dq2 = d2(1);
    typename Derived3::Scalar dr2 = d2(2);

    // Test triangle 1 vertices permutations
    if (dp1 > 0.) 
    {
        if (dq1 > 0.) 
            return triangles3DIntersect_testTriangle2VerticesPermutations(r1,p1,q1,p2,r2,q2,N1,dp2,dr2,dq2);
        else if (dr1 > 0.) 
            return triangles3DIntersect_testTriangle2VerticesPermutations(q1,r1,p1,p2,r2,q2,N1,dp2,dr2,dq2);
        else 
            return triangles3DIntersect_testTriangle2VerticesPermutations(p1,q1,r1,p2,q2,r2,N1,dp2,dq2,dr2);
    } 
    else if (dp1 < 0.) 
    {
        if (dq1 < 0.) 
            return triangles3DIntersect_testTriangle2VerticesPermutations(r1,p1,q1,p2,q2,r2,N1,dp2,dq2,dr2);
        else if (dr1 < 0.) 
            return triangles3DIntersect_testTriangle2VerticesPermutations(q1,r1,p1,p2,q2,r2,N1,dp2,dq2,dr2);
        else 
            return triangles3DIntersect_testTriangle2VerticesPermutations(p1,q1,r1,p2,r2,q2,N1,dp2,dr2,dq2);
    } 
    else 
    {
        if (dq1 < 0.) 
        {
            if (dr1 >= 0.) 
                return triangles3DIntersect_testTriangle2VerticesPermutations(q1,r1,p1,p2,r2,q2,N1,dp2,dr2,dq2);
            else 
                return triangles3DIntersect_testTriangle2VerticesPermutations(p1,q1,r1,p2,q2,r2,N1,dp2,dq2,dr2);
        }
        else if (dq1 > 0.) 
        {
            if (dr1 > 0.) 
                return triangles3DIntersect_testTriangle2VerticesPermutations(p1,q1,r1,p2,r2,q2,N1,dp2,dr2,dq2);
            else 
                return triangles3DIntersect_testTriangle2VerticesPermutations(q1,r1,p1,p2,q2,r2,N1,dp2,dq2,dr2);
        }
        else  
        {
            if (dr1 > 0.) 
                return triangles3DIntersect_testTriangle2VerticesPermutations(r1,p1,q1,p2,q2,r2,N1,dp2,dq2,dr2);
            else if (dr1 < 0.) 
                return triangles3DIntersect_testTriangle2VerticesPermutations(r1,p1,q1,p2,r2,q2,N1,dp2,dr2,dq2);
            else 
                return coplanarTriangles3DIntersect(p1,q1,r1,p2,q2,r2,N1);
        }
    }
}

template< typename Derived1, typename Derived2 >
bool triangles3DIntersect( 
        // Triangle 1 vertices
        Eigen::MatrixBase<Derived1> const& p1,
        Eigen::MatrixBase<Derived1> const& q1,
        Eigen::MatrixBase<Derived1> const& r1,
        // Triangle 2 vertices
        Eigen::MatrixBase<Derived2> const& p2,
        Eigen::MatrixBase<Derived2> const& q2,
        Eigen::MatrixBase<Derived2> const& r2 
        )
{
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived1,3)
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived2,3)

    // Define sgn lambda to compute sign of value (-1, 0, +1)
    auto const sgn = [](auto val)->int
    {
        auto zero = decltype(val)(0);
        return (zero < val) - (val < zero);
    };

    // Compute distance signs of p1, q1 and r1 to the plane of triangle(p2,q2,r2)
    auto const p2mr2 = p2 - r2;
    auto const q2mr2 = q2 - r2;
    auto const N2 = p2mr2.cross(q2mr2).eval();

    Eigen::Matrix<typename Derived1::Scalar, 3, 1> d1;
    auto const p1mr2 = p1 - r2;
    d1(0) = p1mr2.dot(N2);
    auto const q1mr2 = q1 - r2;
    d1(1) = q1mr2.dot(N2);
    auto const r1mr2 = r1 - r2;
    d1(2) = r1mr2.dot(N2);

    // Test planes intersection
    int dp1sgn = sgn(d1(0));
    int dq1sgn = sgn(d1(1));
    int dr1sgn = sgn(d1(2));
    if( (dp1sgn == dq1sgn) && (dp1sgn == dr1sgn) && dp1sgn != 0 )
        return false;

    // Compute distance signs of p2, q2 and r2 to the plane of triangle(p1,q1,r1)
    auto const q1mp1 = q1 - p1;
    auto const r1mp1 = r1 - p1;
    auto const N1 = q1mp1.cross(r1mp1).eval();

    Eigen::Matrix<typename Derived2::Scalar, 3, 1> d2;
    auto const p2mr1 = p2 - r1;
    d2(0) = p2mr1.dot(N1);
    auto const q2mr1 = q2 - r1;
    d2(1) = q2mr1.dot(N1);
    auto const r2mr1 = r2 - r1;
    d2(2) = r2mr1.dot(N1);

    // Test planes intersection
    int dp2sgn = sgn(d2(0));
    int dq2sgn = sgn(d2(1));
    int dr2sgn = sgn(d2(2));
    if( (dp2sgn == dq2sgn) && (dp2sgn == dr2sgn) && dp2sgn != 0 )
        return false;

    // Test the different triangles vertices permutations
    return triangles3DIntersect_testTrianglesVerticesPermutations(p1,q1,r1,p2,q2,r2,N1,N2,d1,d2);
}


} // namespace geometry
} // namespace detail
} // namespace Feel

#endif // _TRIANGLES_INTERSECT_HPP
