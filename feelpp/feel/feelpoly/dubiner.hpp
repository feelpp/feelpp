/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-06

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

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
   \file dubiner.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-06
 */
#ifndef __Dubiner_H
#define __Dubiner_H 1

#include <vector>

#include <boost/lambda/if.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelpoly/expansions.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelmesh/pointset.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/warpblend.hpp>
#include <feel/feelpoly/expansiontypes.hpp>

namespace Feel
{

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
class Dubiner;


template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct DubinerTraits
{
    static inline const uint16_type nDim = Dim;
    static inline const uint16_type nRealDim = RealDim;
    static inline const uint16_type nOrder = Degree;
    static inline const uint16_type nConvexOrderDiff = nDim+nOrder+1;
    static inline const bool is_normalized = NormalizationPolicy::is_normalized;

    /** @name Typedefs
     */
    //@{

    /*
     * numerical type
     */
    typedef T value_type;

    template<uint16_type order, typename V = value_type>
    struct Convex
    {
        typedef Simplex<nDim, order, nDim/*nRealDim*/> type;
        typedef Reference<Simplex<nDim, order, nDim/*nRealDim*/>, nDim, order, nDim/*nRealDim*/, V>  reference_type;
    };

    template<typename NewT>
    struct ChangeValueType
    {
        typedef Dubiner<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy> type;
        typedef DubinerTraits<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy> traits_type;
    };

    template<uint16_type NewOrder>
    struct ChangeOrder
    {
        typedef Dubiner<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy> type;
        typedef DubinerTraits<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy> traits_type;
    };

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef typename Convex<nOrder>::type convex_type;
    typedef typename Convex<nOrder>::reference_type reference_convex_type;

    typedef typename Convex<nConvexOrderDiff>::type diff_convex_type;
    typedef typename Convex<nConvexOrderDiff>::reference_type diff_reference_convex_type;

    using diff_pointset_type = if_t<nDim == 2,
                                    PointSetWarpBlend<diff_convex_type, nConvexOrderDiff, value_type>,
                                    PointSetEquiSpaced<diff_convex_type, nConvexOrderDiff, value_type>>;

    /*
     * storage policy
     */
    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::points_type points_type;
    typedef typename storage_policy::node_type node_type;
}; // class DubinerTraits

template<int D, int O>
struct DubinerTag
{
    static const int Dim = D;
    static const int Order = O;
};
/**
 * \class Dubiner
 * \brief Dubiner polynomial orthonormal basis
 *
 * This class represents the Dubiner polynomials up to degree \c
 * Degree on a simplex in dimension \c Dim.
 *
 *
 * The dubiner polynomials in 1D, the segment \f$[-1;1]\f$ are defined
 * using Jacobi polynomials as follows:
 * \f$ \phi_i(x) = P_i^{0,0}(x) \f$
 *where \f$P_i^{0,0}(x)\f$ is the i-th Jacobi polynomial evaluated at
 * \f$x \in [-1;1]\f$ with weights \f$(0,0)\f$.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 *
 * @see Robert C. Kirby Algorithm 839: FIAT - A New Paradigm for
 * Computing Finite Element Basis Functions, ACM Trans. Math. Software
 * Vol. 30 No. 4 pp 502-516
 *
 * @see M. Dubiner. Spectral methods on triangles and other
 * domains. J. Sci. Comput., 6:345–390, 1991.
 *
 * @see G.E. Karniadakis and S.J. Sherwin, ''Spectral/hp Element
 * Methods for CFD,'' Oxford University Press, March 1999.
 *
 */
template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Dubiner

{

public:
    typedef DubinerTraits<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> traits_type;

    static inline const uint16_type nDim = traits_type::nDim;
    static inline const uint16_type nRealDim = traits_type::nRealDim;
    static inline const uint16_type nOrder = traits_type::nOrder;
    static inline const uint16_type nConvexOrderDiff = traits_type::nConvexOrderDiff;
    static inline const bool is_normalized = traits_type::is_normalized;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = false;
    static inline const bool is_product = true;
    typedef Discontinuous continuity_type;

    /** @name Typedefs
     */
    //@{

    typedef Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> self_type;

    /*
     * for now can be used only as a basis but we might be interested
     * to have then expressed in other basis like in the moment basis
     */
    typedef self_type basis_type;

    typedef typename traits_type::value_type value_type;

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;

    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    /*
     * storage policy
     */
    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Dubiner();
    Dubiner( Dubiner const & d )
        :
        M_refconvex(),
        M_pts( d.M_pts ),
        M_D( d.M_D )
    {

    }
    Dubiner( Dubiner && d ) = default;
    
    ~Dubiner() = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& d )
    {
        if ( this != &d )
        {
            M_pts = d.M_pts;
            M_D = d.M_D;
        }

        return *this;
    }
    self_type& operator=( self_type && d ) = default;

    //ublas::matrix_column<matrix_type const> operator()( node_type const& pt ) const
    matrix_type operator()( node_type const& pt ) const
    {
        points_type pts( pt.size(), 1 );
        ublas::column( pts, 0 ) = pt;
        return evaluate( pts );
    }

    matrix_type operator()( points_type const& pts ) const
    {
        return evaluate( pts );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * Number of polynomials in set
     */
    size_type size() const
    {
        return convex_type::polyDims( nOrder );
    }

    /**
     * \return the maximum degree of the Dubiner polynomial to be
     * constructed
     */
    uint16_type degree() const
    {
        return nOrder;
    }

    /**
     * \return self as a basis
     */
    self_type const& basis() const
    {
        return *this;
    }

    /**
     * \return true if the Dubiner polynomials are normalized, false
     * otherwise
     */
    bool isNormalized() const
    {
        return is_normalized;
    }

    /**
     * \return the \c familyName()
     */
    std::string familyName() const
    {
        return "dubiner";
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{



    /**
     * Dubiner polynomials is an orthonormal basis, the coefficients
     * of the polynomials of the basis are the canonical vectors and
     * represented by the identity matrix (lines are polynomials and
     * columns are the polynomial basis )
     *
     * This function is correct only if we use the Dubiner polynomials
     * as a basis
     */
    matrix_type coeff() const
    {
#if 0
        std::cout << "[Dubiner::coeff] coeff = "
                  << ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() )
                  << "\n";
#endif
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() );
    }


    /**
     * evaluate the Dubiner polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    static matrix_type evaluate( points_type const& __pts )
    {
        return evaluate( __pts, int_c<nDim>{} );
    }

    /**
     * @brief Evaluate Dubiner polynomials up to given order
     *
     * For Dynamic order support: when order > nOrder (compile-time),
     * falls back to runtime evaluation using dyna::JacobiBatchEvaluation.
     *
     * @param __pts Points to evaluate at
     * @param order Polynomial order (can be > nOrder for dynamic case)
     * @return Evaluation matrix
     */
    static matrix_type evaluate( points_type const& __pts, uint16_type order )
    {
        // For Dynamic order support: use runtime evaluation when order > nOrder
        if ( order > nOrder )
        {
            return evaluateRuntime( __pts, order, int_c<nDim>{} );
        }

        // Static path: evaluate at compile-time order and truncate
        auto full = evaluate( __pts );
        const size_type nrows = convex_type::polyDims( order );
        if ( nrows == full.size1() )
            return full;
        matrix_type out( nrows, full.size2() );
        ublas::project( out, ublas::range( 0, nrows ), ublas::range( 0, full.size2() ) ) =
            ublas::project( full, ublas::range( 0, nrows ), ublas::range( 0, full.size2() ) );
        return out;
    }

    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts )
    {
        return derivate( __pts, int_c<nDim>{} );
    }

    /**
     * @brief Derivate Dubiner polynomials up to given order
     *
     * For Dynamic order support: when order > nOrder (compile-time),
     * falls back to runtime derivation using dyna::JacobiBatchDerivation.
     *
     * @param __pts Points to evaluate at
     * @param order Polynomial order (can be > nOrder for dynamic case)
     * @return Vector of derivation matrices (one per spatial dimension)
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, uint16_type order )
    {
        // For Dynamic order support: use runtime derivation when order > nOrder
        if ( order > nOrder )
        {
            return derivateRuntime( __pts, order, int_c<nDim>{} );
        }

        // Static path: evaluate at compile-time order and truncate
        auto full = derivate( __pts );
        const size_type nrows = convex_type::polyDims( order );
        if ( full.size() == 0 || nrows == full[0].size1() )
            return full;
        vector_matrix_type out( full.size() );
        for ( size_type i = 0; i < full.size(); ++i )
        {
            out[i].resize( nrows, full[i].size2() );
            ublas::project( out[i], ublas::range( 0, nrows ), ublas::range( 0, full[i].size2() ) ) =
                ublas::project( full[i], ublas::range( 0, nrows ), ublas::range( 0, full[i].size2() ) );
        }
        return out;
    }

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i ) const
    {
        return M_D[i];
    }

    /**
     * @brief Compute derivation matrix at runtime for given order
     *
     * Used for dynamic polynomial order support. Computes the derivation
     * matrix D_i such that d/dx_i f = D_i * f for polynomial coefficients f.
     *
     * @param i Derivative direction (0: x, 1: y, 2: z)
     * @param order Polynomial order
     * @return Derivation matrix of size (nBasis x nBasis) for the given order
     */
    matrix_type d( uint16_type i, uint16_type order ) const
    {
        // For compile-time order, just return precomputed matrix
        if ( order == nOrder )
            return M_D[i];

        // Compute derivation matrix at runtime for given order
        // Following the same algorithm as the constructor
        // The pointset must have exactly nBasis points for the matrix to be square
        // For Simplex in nDim dimensions: nBasis = (order+1)*(order+2)/2 for 2D, etc.
        using diff_convex_type = Simplex<nDim, 1, nDim>;

        // Use dynamic order pointset to get exactly nBasis points for the given order
        // Use brace initialization to avoid most vexing parse
        PointSetEquiSpaced<diff_convex_type, Dynamic, value_type> diff_pts_gen{ RuntimeOrder{ order } };
        points_type diff_pts = diff_pts_gen.points();

        // Evaluate basis at differentiation points - A should be square (nBasis x nBasis)
        matrix_type A = evaluate( diff_pts, order );

        // Invert A to get interpolation matrix
        matrix_type D_mat = ublas::identity_matrix<value_type>( A.size1(), A.size2() );
        LU<matrix_type> lu( A );
        matrix_type C = lu.solve( D_mat );

        // Compute derivatives at the points
        vector_matrix_type d_vec = derivate( diff_pts, order );

        // Derivation matrix: D_i = d_vec[i] * C
        matrix_type result = ublas::prod( d_vec[i], C );
        glas::clean( result );
        return result;
    }

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& derivate( uint16_type i ) const
    {
        return M_D[i];
    }

    //@}

private:
private:


    static matrix_type
    evaluate( points_type const& __pts, int_c<0> )
        {
            matrix_type m(1,1);
            m(0,0)=1;
            return m;
        }
    /**
     * Evaluation at a set of points of the expansion basis in 1D on
     * the line
     */
    static matrix_type
    evaluate( points_type const& __pts, int_c<1> )
    {
        // Delegate to unified implementation with compile-time order
        return evaluateRuntime( __pts, nOrder, int_c<1>{} );
    }

    template<typename AE>
    static vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts, int_c<0> )
        {
            vector_matrix_type m(1);
            m[0].resize(1,1);
            m[0](0,0)=0;
            return m;
        }
    /**
     * derivation at a set of points of the expansion basis in 1D on
     * the line
     */
    template<typename AE>
    static vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts, int_c<1> )
    {
        // Delegate to unified implementation with compile-time order
        return derivateRuntime( __pts, nOrder, int_c<1>{} );
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type evaluate( points_type const& __pts, int_c<2> );

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, int_c<2> );

    /**
     * Evaluation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    static matrix_type evaluate( points_type const& __pts, int_c<3> );

    /**
     * derivation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, int_c<3> );

    //
    // Runtime evaluation methods for Dynamic order support
    //

    /**
     * @brief Runtime evaluation in 0D (point)
     */
    static matrix_type evaluateRuntime( points_type const& __pts, uint16_type order, int_c<0> )
    {
        (void)order;
        matrix_type m( 1, 1 );
        m( 0, 0 ) = 1;
        return m;
    }

    /**
     * @brief Runtime evaluation in 1D (line)
     */
    static matrix_type evaluateRuntime( points_type const& __pts, uint16_type order, int_c<1> )
    {
        // Convert matrix_row to vector for JacobiBatchEvaluation
        ublas::vector<value_type> pts( __pts.size2() );
        for ( size_t i = 0; i < __pts.size2(); ++i )
            pts( i ) = __pts( 0, i );

        matrix_type m( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), pts ) );

        if ( is_normalized )
        {
            for ( uint16_type i = 0; i < m.size1(); ++i )
                ublas::row( m, i ) *= math::sqrt( value_type( i ) + value_type( 0.5 ) );
        }

        return m;
    }

    /**
     * @brief Runtime evaluation in 2D (triangle)
     */
    static matrix_type evaluateRuntime( points_type const& __pts, uint16_type order, int_c<2> );

    /**
     * @brief Runtime evaluation in 3D (tetrahedron)
     */
    static matrix_type evaluateRuntime( points_type const& __pts, uint16_type order, int_c<3> );

    //
    // Runtime derivation methods for Dynamic order support
    //

    /**
     * @brief Runtime derivation in 0D (point)
     */
    template<typename AE>
    static vector_matrix_type derivateRuntime( ublas::matrix_expression<AE> const& __pts,
                                                uint16_type order, int_c<0> )
    {
        (void)order;
        vector_matrix_type m( 1 );
        m[0].resize( 1, 1 );
        m[0]( 0, 0 ) = 0;
        return m;
    }

    /**
     * @brief Runtime derivation in 1D (line)
     */
    template<typename AE>
    static vector_matrix_type derivateRuntime( ublas::matrix_expression<AE> const& __pts,
                                                uint16_type order, int_c<1> )
    {
        vector_matrix_type D( 1 );
        D[0].resize( order + 1, __pts().size2() );

        // Copy matrix row to vector for JacobiBatchDerivation
        ublas::vector<value_type> pts_vec( __pts().size2() );
        for ( size_type k = 0; k < __pts().size2(); ++k )
            pts_vec( k ) = __pts()( 0, k );

        D[0] = dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), pts_vec );

        if ( is_normalized )
            for ( uint16_type i = 0; i <= order; ++i )
                ublas::row( D[0], i ) *= math::sqrt( value_type( i ) + value_type( 0.5 ) );

        return D;
    }

    /**
     * @brief Runtime derivation in 2D (triangle)
     */
    template<typename AE>
    static vector_matrix_type derivateRuntime( ublas::matrix_expression<AE> const& __pts,
                                                uint16_type order, int_c<2> );

    /**
     * @brief Runtime derivation in 3D (tetrahedron)
     */
    template<typename AE>
    static vector_matrix_type derivateRuntime( ublas::matrix_expression<AE> const& __pts,
                                                uint16_type order, int_c<3> );

private:
    reference_convex_type M_refconvex;
    points_type M_pts;

    /**
     * Derivation matrix
     * \note construct it only once per dubiner polynomials
     */
    std::vector<matrix_type> M_D;

}; // class Dubiner

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::Dubiner()
    :
    M_refconvex(),
    M_pts( M_refconvex.makePoints( Dim, 0 ) ),
    M_D( Dim )
{
    reference_convex_type refconvex;
    // constructor pointset for differentiation only in
    // the interior(1)
    diff_pointset_type diff_pts( 1 );
    matrix_type A( evaluate( diff_pts.points() ) );
    
    matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
    LU<matrix_type> lu( A );
    matrix_type C = lu.solve( D );
    
    vector_matrix_type d ( derivate( diff_pts.points() ) );
    for ( size_type i = 0; i < d.size(); ++i )
    {
        M_D[i] = ublas::prod( d[i], C );
        glas::clean( M_D[i] );
    }
    
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, int_c<2> )
{
    // Delegate to unified implementation with compile-time order
    return evaluateRuntime( __pts, nOrder, int_c<2>{} );
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Dubiner<Dim, RealDim,  Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, int_c<2> )
{
    // Delegate to unified implementation with compile-time order
    return derivateRuntime( __pts, nOrder, int_c<2>{} );
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, int_c<3> )
{
    // Delegate to unified implementation with compile-time order
    return evaluateRuntime( __pts, nOrder, int_c<3>{} );
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, int_c<3> )
{
    // Delegate to unified implementation with compile-time order
    return derivateRuntime( __pts, nOrder, int_c<3>{} );
}

//
// Runtime evaluation implementations for Dynamic order support
//

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluateRuntime(
    points_type const& __pts, uint16_type order, int_c<2> )
{
    // Number of DOFs for simplex of order 'order' in 2D: (order+1)(order+2)/2
    const size_type ndof = ( order + 1 ) * ( order + 2 ) / 2;
    matrix_type res( ndof, __pts.size2() );

    details::etas<TRIANGLE, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );

    // Use runtime Jacobi evaluation
    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    std::vector<matrix_type> bs( order + 1 );

    for ( uint16_type i = 0; i <= order; ++i )
    {
        bs[i].resize( order - i + 1, eta2s.size() );
        bs[i] = dyna::JacobiBatchEvaluation( order - i, value_type( 2 * i + 1 ), value_type( 0 ), eta2s );
    }

    // Use runtime scalings
    matrix_type scalings = details::scalingsRuntime( order, eta2s );

    for ( uint16_type cur = 0, k = 0; k <= order; ++k )
    {
        for ( uint16_type i = 0; i <= k; ++i, ++cur )
        {
            uint16_type ii = k - i;
            uint16_type jj = i;

            if ( is_normalized )
            {
                value_type normalization = math::sqrt( ( value_type( ii ) + value_type( 0.5 ) ) *
                                                        ( value_type( ii + jj ) + value_type( 1 ) ) );

                for ( uint16_type l = 0; l < as.size2(); ++l )
                    res( cur, l ) = normalization * as( ii, l ) * scalings( ii, l ) * bs[ii]( jj, l );
            }
            else
            {
                for ( uint16_type l = 0; l < as.size2(); ++l )
                    res( cur, l ) = as( ii, l ) * scalings( ii, l ) * bs[ii]( jj, l );
            }
        }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluateRuntime(
    points_type const& __pts, uint16_type order, int_c<3> )
{
    // Number of DOFs for simplex of order 'order' in 3D: (order+1)(order+2)(order+3)/6
    const size_type ndof = ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
    matrix_type res( ndof, __pts.size2() );

    DCHECK( __pts.size1() == 3 ) << "invalid space dimension for point set, should be 3, it is currently "
                                  << __pts.size1();

    details::etas<TETRAHEDRON, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );
    ublas::vector<value_type> eta3s = ublas::row( etas(), 2 );

    // Use runtime Jacobi evaluation
    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    std::vector<matrix_type> bs( order + 1 );
    ublas::matrix<matrix_type> cs( order + 1, order + 1 );

    for ( uint16_type i = 0; i <= order; ++i )
    {
        bs[i].resize( order - i + 1, eta2s.size() );
        bs[i] = dyna::JacobiBatchEvaluation( order - i, value_type( 2 * i + 1 ), value_type( 0 ), eta2s );

        for ( uint16_type j = 0; j <= order - i; ++j )
        {
            cs( i, j ).resize( order - i - j + 1, eta3s.size() );
            cs( i, j ) = dyna::JacobiBatchEvaluation( order - i - j,
                                                       value_type( 2 * ( i + j + 1 ) ), value_type( 0 ), eta3s );
        }
    }

    // Use runtime scalings
    matrix_type scalings2 = details::scalingsRuntime( order, eta2s );
    matrix_type scalings3 = details::scalingsRuntime( order, eta3s );

    for ( uint16_type cur = 0, k = 0; k <= order; ++k )
    {
        for ( uint16_type i = 0; i <= k; ++i )
        {
            for ( uint16_type j = 0; j <= k - i; ++j, ++cur )
            {
                uint16_type ii = k - i - j;
                uint16_type jj = j;
                uint16_type kk = i;

                if ( is_normalized )
                {
                    value_type normalization = math::sqrt( ( value_type( ii ) + value_type( 0.5 ) ) *
                                                            ( value_type( ii + jj ) + value_type( 1 ) ) *
                                                            ( value_type( ii + jj + kk ) + value_type( 1.5 ) ) );

                    for ( uint16_type l = 0; l < as.size2(); ++l )
                        res( cur, l ) = normalization * ( as( ii, l ) *
                                                           scalings2( ii, l ) * bs[ii]( jj, l ) *
                                                           scalings3( ii + jj, l ) * cs( ii, jj )( kk, l ) );
                }
                else
                {
                    for ( uint16_type l = 0; l < as.size2(); ++l )
                        res( cur, l ) = as( ii, l ) *
                                        scalings2( ii, l ) * bs[ii]( jj, l ) *
                                        scalings3( ii + jj, l ) * cs( ii, jj )( kk, l );
                }
            }
        }
    }

    return res;
}

//
// Runtime derivation implementations for Dynamic order support
//

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivateRuntime(
    ublas::matrix_expression<AE> const& __pts, uint16_type order, int_c<2> )
{
    const size_type ndof = ( order + 1 ) * ( order + 2 ) / 2;
    vector_matrix_type res( 2 );
    res[0].resize( ndof, __pts().size2() );
    res[1].resize( ndof, __pts().size2() );

    details::etas<TRIANGLE, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );

    // Runtime Jacobi polynomials
    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    matrix_type das( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    std::vector<matrix_type> bs( order + 1 );
    std::vector<matrix_type> dbs( order + 1 );

    for ( uint16_type i = 0; i <= order; ++i )
    {
        bs[i].resize( order - i + 1, eta2s.size() );
        dbs[i].resize( order - i + 1, eta2s.size() );
        bs[i] = dyna::JacobiBatchEvaluation( order - i, value_type( 2 * i + 1 ), value_type( 0 ), eta2s );
        dbs[i] = dyna::JacobiBatchDerivation( order - i, value_type( 2 * i + 1 ), value_type( 0 ), eta2s );
    }

    matrix_type scalings = details::scalingsRuntime( order, eta2s );
    ublas::vector<value_type> one( ublas::scalar_vector<value_type>( eta1s.size(), value_type( 1 ) ) );
    ublas::vector<value_type> tmp( ublas::scalar_vector<value_type>( eta1s.size(), value_type( 1 ) ) );

    for ( uint16_type k = 0, cur = 0; k <= order; ++k )
    {
        for ( uint16_type i = 0; i <= k; ++i, ++cur )
        {
            uint16_type ii = k - i;
            uint16_type jj = i;

            // x derivation
            ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                                              ublas::row( bs[ii], jj ) );

            if ( ii > 0 )
                ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( res[0], cur ),
                                                                  ublas::row( scalings, ii - 1 ) );

            // y derivation
            ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                                              ublas::row( bs[ii], jj ) );
            ublas::row( res[1], cur ) = value_type( 0.5 ) * ublas::element_prod( ublas::row( res[1], cur ),
                                                                                   ( one + eta1s ) );

            if ( ii > 0 )
                ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( res[1], cur ),
                                                                  ublas::row( scalings, ii - 1 ) );

            // derivate (1-x)^ii
            tmp = ublas::element_prod( ublas::row( scalings, ii ),
                                        ublas::row( dbs[ii], jj ) );

            if ( ii > 0 )
                tmp -= value_type( 0.5 ) * value_type( ii ) *
                       ublas::element_prod( ublas::row( scalings, ii - 1 ),
                                             ublas::row( bs[ii], jj ) );

            // add contrib to y derivation
            ublas::row( res[1], cur ) += ublas::element_prod( ublas::row( as, ii ), tmp );

            // orthonormalize if required
            if ( is_normalized )
            {
                value_type normalization = math::sqrt( ( value_type( ii ) + value_type( 0.5 ) ) *
                                                        ( value_type( ii + jj ) + value_type( 1 ) ) );
                ublas::row( res[0], cur ) *= normalization;
                ublas::row( res[1], cur ) *= normalization;
            }
        }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivateRuntime(
    ublas::matrix_expression<AE> const& __pts, uint16_type order, int_c<3> )
{
    const size_type ndof = ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
    vector_matrix_type res( 3 );
    res[0].resize( ndof, __pts().size2() );
    res[1].resize( ndof, __pts().size2() );
    res[2].resize( ndof, __pts().size2() );

    FEELPP_ASSERT( __pts().size1() == 3 )( __pts().size1() ).error( "invalid space dimension" );

    details::etas<TETRAHEDRON, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );
    ublas::vector<value_type> eta3s = ublas::row( etas(), 2 );

    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    matrix_type das( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    std::vector<matrix_type> bs( order + 1 );
    std::vector<matrix_type> dbs( order + 1 );
    ublas::matrix<matrix_type> cs( order + 1, order + 1 );
    ublas::matrix<matrix_type> dcs( order + 1, order + 1 );

    for ( uint16_type i = 0; i <= order; ++i )
    {
        bs[i].resize( order - i + 1, eta2s.size() );
        dbs[i].resize( order - i + 1, eta2s.size() );
        bs[i] = dyna::JacobiBatchEvaluation( order - i, value_type( 2 * i + 1 ), value_type( 0 ), eta2s );
        dbs[i] = dyna::JacobiBatchDerivation( order - i, value_type( 2 * i + 1 ), value_type( 0 ), eta2s );

        for ( uint16_type j = 0; j <= order - i; ++j )
        {
            cs( i, j ).resize( order - i - j + 1, eta3s.size() );
            dcs( i, j ).resize( order - i - j + 1, eta3s.size() );
            cs( i, j ) = dyna::JacobiBatchEvaluation( order - i - j,
                                                       value_type( 2 * ( i + j + 1 ) ), value_type( 0 ), eta3s );
            dcs( i, j ) = dyna::JacobiBatchDerivation( order - i - j,
                                                        value_type( 2 * ( i + j + 1 ) ), value_type( 0 ), eta3s );
        }
    }

    matrix_type scalings2 = details::scalingsRuntime( order, eta2s );
    matrix_type scalings3 = details::scalingsRuntime( order, eta3s );

    ublas::vector<value_type> one( ublas::scalar_vector<value_type>( eta1s.size(), value_type( 1 ) ) );
    ublas::vector<value_type> tmp( ublas::scalar_vector<value_type>( eta1s.size(), value_type( 1 ) ) );

    for ( uint16_type cur = 0, k = 0; k <= order; ++k )
    {
        for ( uint16_type i = 0; i <= k; ++i )
        {
            for ( uint16_type j = 0; j <= k - i; ++j, ++cur )
            {
                uint16_type ii = k - i - j;
                uint16_type jj = j;
                uint16_type kk = i;

                // x derivation
                ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                                                  ublas::row( bs[ii], jj ) );
                ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( res[0], cur ),
                                                                  ublas::row( cs( ii, jj ), kk ) );

                if ( ii > 0 )
                    ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( res[0], cur ),
                                                                      ublas::row( scalings2, ii - 1 ) );

                if ( ii + jj > 0 )
                    ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( res[0], cur ),
                                                                      ublas::row( scalings3, ii + jj - 1 ) );

                // y derivation
                ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                                                  ublas::row( bs[ii], jj ) );
                ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( res[1], cur ),
                                                                  ublas::row( cs( ii, jj ), kk ) );
                ublas::row( res[1], cur ) = value_type( 0.5 ) * ublas::element_prod( ublas::row( res[1], cur ),
                                                                                       ( one + eta1s ) );

                if ( ii > 0 )
                    ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( res[1], cur ),
                                                                      ublas::row( scalings2, ii - 1 ) );

                if ( ii + jj > 0 )
                    ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( res[1], cur ),
                                                                      ublas::row( scalings3, ii + jj - 1 ) );

                // derivate (1-x)^ii for y
                tmp = ublas::element_prod( ublas::row( scalings2, ii ),
                                            ublas::row( dbs[ii], jj ) );

                if ( ii > 0 )
                    tmp -= value_type( 0.5 ) * value_type( ii ) *
                           ublas::element_prod( ublas::row( scalings2, ii - 1 ),
                                                 ublas::row( bs[ii], jj ) );

                tmp = ublas::element_prod( tmp, ublas::row( as, ii ) );
                tmp = ublas::element_prod( tmp, ublas::row( cs( ii, jj ), kk ) );

                if ( ii + jj > 0 )
                    tmp = ublas::element_prod( tmp, ublas::row( scalings3, ii + jj - 1 ) );

                ublas::row( res[1], cur ) += tmp;

                // z derivation
                ublas::row( res[2], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                                                  ublas::row( bs[ii], jj ) );
                ublas::row( res[2], cur ) = ublas::element_prod( ublas::row( res[2], cur ),
                                                                  ublas::row( cs( ii, jj ), kk ) );
                ublas::row( res[2], cur ) = value_type( 0.5 ) * ublas::element_prod( ublas::row( res[2], cur ),
                                                                                       ( one + eta1s ) );

                if ( ii > 0 )
                    ublas::row( res[2], cur ) = ublas::element_prod( ublas::row( res[2], cur ),
                                                                      ublas::row( scalings2, ii - 1 ) );

                if ( ii + jj > 0 )
                    ublas::row( res[2], cur ) = ublas::element_prod( ublas::row( res[2], cur ),
                                                                      ublas::row( scalings3, ii + jj - 1 ) );

                // derivate (1-x)^ii for z (same as y)
                tmp = ublas::element_prod( ublas::row( scalings2, ii ),
                                            ublas::row( dbs[ii], jj ) );

                if ( ii > 0 )
                    tmp -= value_type( 0.5 ) * value_type( ii ) *
                           ublas::element_prod( ublas::row( scalings2, ii - 1 ),
                                                 ublas::row( bs[ii], jj ) );

                tmp = ublas::element_prod( tmp, ublas::row( as, ii ) );
                tmp = ublas::element_prod( tmp, ublas::row( cs( ii, jj ), kk ) );
                tmp = value_type( 0.5 ) * ublas::element_prod( tmp, ( one + eta2s ) );

                if ( ii + jj > 0 )
                    tmp = ublas::element_prod( tmp, ublas::row( scalings3, ii + jj - 1 ) );

                ublas::row( res[2], cur ) += tmp;

                // derivate scaling for z
                tmp = ublas::element_prod( ublas::row( scalings3, ii + jj ),
                                            ublas::row( dcs( ii, jj ), kk ) );

                if ( ii + jj > 0 )
                    tmp -= value_type( 0.5 ) * value_type( ii + jj ) *
                           ublas::element_prod( ublas::row( cs( ii, jj ), kk ),
                                                 ublas::row( scalings3, ii + jj - 1 ) );

                tmp = ublas::element_prod( tmp, ublas::row( as, ii ) );
                tmp = ublas::element_prod( tmp, ublas::row( bs[ii], jj ) );
                tmp = ublas::element_prod( tmp, ublas::row( scalings2, ii ) );

                ublas::row( res[2], cur ) += tmp;

                if ( is_normalized )
                {
                    value_type normalization = math::sqrt( ( value_type( ii ) + value_type( 0.5 ) ) *
                                                            ( value_type( ii + jj ) + value_type( 1 ) ) *
                                                            ( value_type( ii + jj + kk ) + value_type( 1.5 ) ) );

                    ublas::row( res[0], cur ) *= normalization;
                    ublas::row( res[1], cur ) *= normalization;
                    ublas::row( res[2], cur ) *= normalization;
                }
            }
        }
    }

    return res;
}

}
#endif /* __Dubiner_H */
