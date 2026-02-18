/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-02-20

  Copyright (C) 2006 EPFL
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
   \file legendre.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-02-20
 */
#ifndef __Legendre_H
#define __Legendre_H 1

#include <vector>

#include <boost/lambda/if.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelpoly/expansions.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/gausslobatto.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/expansiontypes.hpp>

namespace Feel
{
template< class Convex,
          int Order,
          typename T >
class PointSetGaussLobatto;

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
class Legendre;


template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct LegendreTraits
{
    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = RealDim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrderDiff = nOrder+2;
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
        typedef Hypercube<nDim, order, /*nRealDim*/nDim> type;
        typedef Reference<Hypercube<nDim, order, /*nRealDim*/nDim>, nDim, order, nDim/*nRealDim*/, V>  reference_type;
    };

    template<typename NewT>
    struct ChangeValueType
    {
        typedef Legendre<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy> type;
        typedef LegendreTraits<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy> traits_type;
    };

    template<uint16_type NewOrder>
    struct ChangeOrder
    {
        typedef Legendre<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy> type;
        typedef LegendreTraits<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy> traits_type;
    };

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef typename Convex<nOrder>::type convex_type;
    typedef typename Convex<nOrder>::reference_type reference_convex_type;

    typedef typename Convex<nConvexOrderDiff>::type diff_convex_type;
    typedef typename Convex<nConvexOrderDiff>::reference_type diff_reference_convex_type;

    typedef PointSetGaussLobatto<diff_convex_type, nConvexOrderDiff, value_type> diff_pointset_type;

    /*
     * storage policy
     */
    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::points_type points_type;
    typedef typename storage_policy::node_type node_type;
};

template<int D, int O>
struct LegendreTag
{
    static const int Dim = D;
    static const int Order = O;
};
/**
 * \class Legendre
 * \brief Legendre polynomial orthonormal basis
 *
 * This class represents the Legendre polynomials up to degree \c
 * Degree on a simplex in dimension \c Dim.
 *
 *
 * The legendre polynomials in 1D, the segment \f$[-1;1]\f$ are defined
 * using Jacobi polynomials as follows:
 * \f$ \phi_i(x) = P_i^{0,0}(x) \f$
 *where \f$P_i^{0,0}(x)\f$ is the i-th Jacobi polynomial evaluated at
 * \f$x \in [-1;1]\f$ with weights \f$(0,0)\f$.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
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
class Legendre
{
public:

    typedef LegendreTraits<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> traits_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = RealDim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrder = nOrder+2;
    static inline const bool is_normalized = NormalizationPolicy::is_normalized;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = false;
    static inline const bool is_product = true;
    typedef Discontinuous continuity_type;

    /** @name Typedefs
     */
    //@{

    typedef Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> self_type;

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

    Legendre()
        :
        M_refconvex(),
        M_pts( M_refconvex.makePoints( nDim, 0 ) )
    {
        this->initDerivation();
    }
    Legendre( Legendre const & d )
        :
        M_refconvex(),
        M_pts( d.M_pts )
    {

    }

    ~Legendre()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& d )
    {
        if ( this != &d )
        {
            M_pts = d.M_pts;
            _S_D = d._S_D;
        }

        return *this;
    }

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
     * \return the maximum degree of the Legendre polynomial to be
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
     * \return true if the Legendre polynomials are normalized, false
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
        return "legendre";
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
     * Legendre polynomials is an orthonormal basis, the coefficients
     * of the polynomials of the basis are the canonical vectors and
     * represented by the identity matrix (lines are polynomials and
     * columns are the polynomial basis )
     *
     * This function is correct only if we use the Legendre polynomials
     * as a basis
     */
    matrix_type coeff() const
    {
#if 0
        std::cout << "[Legendre::coeff] coeff = "
                  << ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() )
                  << "\n";
#endif
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), M_pts.size2() );
    }


    /**
     * evaluate the Legendre polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    static matrix_type evaluate( points_type const& __pts )
    {
        return evaluate( __pts, int_c<nDim>{} );
    }

    /**
     * @brief Evaluate Legendre polynomials up to given order
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
     * @brief Derivate Legendre polynomials up to given order
     *
     * For Dynamic order support: when order > nOrder (compile-time),
     * falls back to runtime derivation using dyna::JacobiBatchDerivation.
     *
     * @param __pts Points to evaluate at
     * @param order Polynomial order (can be > nOrder for dynamic case)
     * @return Vector of derivation matrices (one per spatial dimension)
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts, uint16_type order )
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
     * \brief derivatives of Legendre polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    static matrix_type const& d( uint16_type i )
    {
        return _S_D[i];
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
    static matrix_type d( uint16_type i, uint16_type order )
    {
        // For compile-time order, just return precomputed matrix
        if ( order == nOrder )
            return _S_D[i];

        // Compute derivation matrix at runtime for given order
        // The pointset must have exactly nBasis points for the matrix to be square
        // For Hypercube: nBasis = (order+1)^nDim
        using diff_convex_type = Hypercube<nDim, 1, nDim>;

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
     * \brief derivatives of Legendre polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    static matrix_type const& derivate( uint16_type i )
    {
        return _S_D[i];
    }


    //@}

private:
private:

    static value_type normalization( int i )
    {
        return ( is_normalized?math::sqrt( value_type( i ) + 0.5 ) : value_type( 1 ) );
    }
    static value_type normalization( int i, int j )
    {
        return ( is_normalized?math::sqrt( ( value_type( i ) + 0.5 ) *
                                           ( value_type( j ) + 0.5 ) ) : value_type( 1 ) );
    }
    static value_type normalization( int i, int j, int k )
    {
        return ( is_normalized?math::sqrt( ( value_type( i ) + 0.5 ) *
                                           ( value_type( j ) + 0.5 ) *
                                           ( value_type( k ) + 0.5 ) ) : value_type( 1 ) );
    }
    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type
    evaluate( points_type const& __pts, int_c<1> )
    {
        // Delegate to unified implementation with compile-time order
        return evaluateRuntime( __pts, nOrder, int_c<1>{} );
    }

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
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
                ublas::row( m, i ) *= normalization( i );
        }

        return m;
    }

    /**
     * @brief Runtime evaluation in 2D (quad)
     */
    static matrix_type evaluateRuntime( points_type const& __pts, uint16_type order, int_c<2> );

    /**
     * @brief Runtime evaluation in 3D (hexahedron)
     */
    static matrix_type evaluateRuntime( points_type const& __pts, uint16_type order, int_c<3> );

    //
    // Runtime derivation methods for Dynamic order support
    //

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
                ublas::row( D[0], i ) *= normalization( i );

        return D;
    }

    /**
     * @brief Runtime derivation in 2D (quad)
     */
    template<typename AE>
    static vector_matrix_type derivateRuntime( ublas::matrix_expression<AE> const& __pts,
                                                uint16_type order, int_c<2> );

    /**
     * @brief Runtime derivation in 3D (hexahedron)
     */
    template<typename AE>
    static vector_matrix_type derivateRuntime( ublas::matrix_expression<AE> const& __pts,
                                                uint16_type order, int_c<3> );

    static void initDerivation();
private:
    reference_convex_type M_refconvex;
    points_type M_pts;
    /**
     * \c true if differentation matrix initialized, \c false
     * otherwise
     */
    static bool _S_has_derivation;

    /**
     * Derivation matrix
     * \note construct it only once per dubiner polynomials
     */
    static std::vector<matrix_type> _S_D;

}; // class Dubiner

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
bool Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::_S_has_derivation = false;

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
std::vector<typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type>
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::_S_D;

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
void
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::initDerivation()
{
#if 0
    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;

    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;
#endif // 0

    if ( _S_has_derivation == false )
    {
        _S_has_derivation = true;

        reference_convex_type refconvex;
        // constructor pointset for differentiation only in
        // the interior(1)
        diff_pointset_type diff_pts( 1 );
        matrix_type A( evaluate( diff_pts.points() ) );

#if 1
        matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
        LU<matrix_type> lu( A );
        matrix_type C = lu.solve( D );

        vector_matrix_type d ( derivate( diff_pts.points() ) );
        _S_D.resize( d.size() );

        for ( size_type i = 0; i < d.size(); ++i )
        {
            _S_D[i] = ublas::prod( d[i], C );
            glas::clean( _S_D[i] );
        }

#endif
    }
}

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, int_c<2> )
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, int_c<2> )
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, int_c<3> )
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, int_c<3> )
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluateRuntime(
    points_type const& __pts, uint16_type order, int_c<2> )
{
    // Number of DOFs for hypercube of order 'order' in 2D: (order+1)^2
    const size_type ndof = ( order + 1 ) * ( order + 1 );
    matrix_type res( ndof, __pts.size2() );

    // Copy matrix rows to vectors for JacobiBatchEvaluation
    ublas::vector<value_type> eta1s( __pts.size2() );
    ublas::vector<value_type> eta2s( __pts.size2() );
    for ( size_type k = 0; k < __pts.size2(); ++k )
    {
        eta1s( k ) = __pts( 0, k );
        eta2s( k ) = __pts( 1, k );
    }

    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    matrix_type bs( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta2s ) );

    for ( uint16_type cur = 0, i = 0; i <= order; ++i )
    {
        for ( uint16_type j = 0; j <= order; ++j, ++cur )
        {
            ublas::row( res, cur ) = normalization( i, j ) *
                                      ublas::element_prod( ublas::row( as, i ), ublas::row( bs, j ) );
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluateRuntime(
    points_type const& __pts, uint16_type order, int_c<3> )
{
    // Number of DOFs for hypercube of order 'order' in 3D: (order+1)^3
    const size_type ndof = ( order + 1 ) * ( order + 1 ) * ( order + 1 );
    matrix_type res( ndof, __pts.size2() );

    // Copy matrix rows to vectors for JacobiBatchEvaluation
    ublas::vector<value_type> eta1s( __pts.size2() );
    ublas::vector<value_type> eta2s( __pts.size2() );
    ublas::vector<value_type> eta3s( __pts.size2() );
    for ( size_type l = 0; l < __pts.size2(); ++l )
    {
        eta1s( l ) = __pts( 0, l );
        eta2s( l ) = __pts( 1, l );
        eta3s( l ) = __pts( 2, l );
    }

    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta1s ) );
    matrix_type bs( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta2s ) );
    matrix_type cs( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), eta3s ) );

    for ( uint16_type cur = 0, i = 0; i <= order; ++i )
    {
        for ( uint16_type j = 0; j <= order; ++j )
        {
            for ( uint16_type k = 0; k <= order; ++k, ++cur )
            {
                ublas::row( res, cur ) = normalization( i, j, k ) *
                                          ublas::element_prod( ublas::element_prod( ublas::row( as, i ),
                                                                                     ublas::row( bs, j ) ),
                                                               ublas::row( cs, k ) );
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivateRuntime(
    ublas::matrix_expression<AE> const& __pts, uint16_type order, int_c<2> )
{
    const size_type ndof = ( order + 1 ) * ( order + 1 );
    vector_matrix_type res( 2 );
    res[0].resize( ndof, __pts().size2() );
    res[1].resize( ndof, __pts().size2() );

    // Copy matrix rows to vectors for JacobiBatchEvaluation/Derivation
    ublas::vector<value_type> pts_x( __pts().size2() );
    ublas::vector<value_type> pts_y( __pts().size2() );
    for ( size_type k = 0; k < __pts().size2(); ++k )
    {
        pts_x( k ) = __pts()( 0, k );
        pts_y( k ) = __pts()( 1, k );
    }

    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), pts_x ) );
    matrix_type das( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), pts_x ) );
    matrix_type bs( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), pts_y ) );
    matrix_type dbs( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), pts_y ) );

    for ( uint16_type cur = 0, i = 0; i <= order; ++i )
    {
        for ( uint16_type j = 0; j <= order; ++j, ++cur )
        {
            ublas::row( res[0], cur ) = normalization( i, j ) *
                                         ublas::element_prod( ublas::row( das, i ), ublas::row( bs, j ) );
            ublas::row( res[1], cur ) = normalization( i, j ) *
                                         ublas::element_prod( ublas::row( as, i ), ublas::row( dbs, j ) );
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
typename Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivateRuntime(
    ublas::matrix_expression<AE> const& __pts, uint16_type order, int_c<3> )
{
    const size_type ndof = ( order + 1 ) * ( order + 1 ) * ( order + 1 );
    vector_matrix_type res( 3 );
    res[0].resize( ndof, __pts().size2() );
    res[1].resize( ndof, __pts().size2() );
    res[2].resize( ndof, __pts().size2() );

    // Copy matrix rows to vectors for JacobiBatchEvaluation/Derivation
    ublas::vector<value_type> pts_x( __pts().size2() );
    ublas::vector<value_type> pts_y( __pts().size2() );
    ublas::vector<value_type> pts_z( __pts().size2() );
    for ( size_type l = 0; l < __pts().size2(); ++l )
    {
        pts_x( l ) = __pts()( 0, l );
        pts_y( l ) = __pts()( 1, l );
        pts_z( l ) = __pts()( 2, l );
    }

    matrix_type as( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), pts_x ) );
    matrix_type das( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), pts_x ) );
    matrix_type bs( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), pts_y ) );
    matrix_type dbs( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), pts_y ) );
    matrix_type cs( dyna::JacobiBatchEvaluation( order, value_type( 0 ), value_type( 0 ), pts_z ) );
    matrix_type dcs( dyna::JacobiBatchDerivation( order, value_type( 0 ), value_type( 0 ), pts_z ) );

    for ( uint16_type cur = 0, i = 0; i <= order; ++i )
    {
        for ( uint16_type j = 0; j <= order; ++j )
        {
            for ( uint16_type k = 0; k <= order; ++k, ++cur )
            {
                ublas::row( res[0], cur ) = ( normalization( i, j, k ) *
                                               ublas::element_prod( ublas::element_prod( ublas::row( das, i ),
                                                                                          ublas::row( bs, j ) ),
                                                                    ublas::row( cs, k ) ) );

                ublas::row( res[1], cur ) = ( normalization( i, j, k ) *
                                               ublas::element_prod( ublas::element_prod( ublas::row( as, i ),
                                                                                          ublas::row( dbs, j ) ),
                                                                    ublas::row( cs, k ) ) );

                ublas::row( res[2], cur ) = ( normalization( i, j, k ) *
                                               ublas::element_prod( ublas::element_prod( ublas::row( as, i ),
                                                                                          ublas::row( bs, j ) ),
                                                                    ublas::row( dcs, k ) ) );
            }
        }
    }

    return res;
}

}
#endif /* __Legendre_H */
