/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-02-20
 */
#ifndef __Legendre_H
#define __Legendre_H 1

#include <vector>

#include <boost/lambda/if.hpp>

#include <life/lifemesh/refentity.hpp>
#include <life/lifealg/glas.hpp>
#include <life/lifealg/lu.hpp>
#include <life/lifepoly/expansions.hpp>
#include <life/lifepoly/policy.hpp>
#include <life/lifepoly/gausslobatto.hpp>
#include <life/lifepoly/equispaced.hpp>
#include <life/lifepoly/basis.hpp>
#include <life/lifepoly/expansiontypes.hpp>

namespace Life
{
template< class Convex,
          uint16_type Order,
          typename T >
class PointSetGaussLobatto;

template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
struct Legendre;


template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct LegendreTraits
{
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrderDiff = nOrder+2;
    static const bool is_normalized = NormalizationPolicy::is_normalized;

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
        typedef SimplexProduct<nDim, order, nDim> type;
        typedef Reference<SimplexProduct<nDim, order, nDim>, nDim, order, nDim, V>  reference_type;
    };

    template<typename NewT>
    struct ChangeValueType
    {
        typedef Legendre<Dim, Degree, NormalizationPolicy, NewT, StoragePolicy> type;
        typedef LegendreTraits<Dim, Degree, NormalizationPolicy, NewT, StoragePolicy> traits_type;
    };

    template<uint16_type NewOrder>
    struct ChangeOrder
    {
        typedef Legendre<Dim, NewOrder, NormalizationPolicy, T, StoragePolicy> type;
        typedef LegendreTraits<Dim, NewOrder, NormalizationPolicy, T, StoragePolicy> traits_type;
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
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Legendre : public Basis<LegendreTag<Dim, Degree>, T>
{
    typedef Basis<LegendreTag<Dim, Degree>, T> super;
public:

    typedef LegendreTraits<Dim, Degree, NormalizationPolicy, T, StoragePolicy> traits_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrder = nOrder+2;
    static const bool is_normalized = NormalizationPolicy::is_normalized;
    static const bool isTransformationEquivalent = true;
    static const bool isContinuous = false;
    typedef Discontinuous continuity_type;

    /** @name Typedefs
     */
    //@{

    typedef Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy> self_type;

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
        super( *this ),
        _M_refconvex(),
        _M_pts( _M_refconvex.makePoints( nDim, 0 ) )
    {
    }
    Legendre( Legendre const & d )
        :
        super( *this ),
        _M_refconvex(),
        _M_pts( d._M_pts )
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
                _M_pts = d._M_pts;
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
    size_type size() const { return convex_type::polyDims( nOrder ); }

    /**
     * \return the maximum degree of the Legendre polynomial to be
     * constructed
     */
    uint16_type degree() const { return nOrder; }

    /**
     * \return self as a basis
     */
    self_type const& basis() const { return *this; }

    /**
     * \return true if the Legendre polynomials are normalized, false
     * otherwise
     */
    bool isNormalized() const { return is_normalized; }

    /**
     * \return the \c familyName()
     */
    std::string familyName() const { return "legendre"; }

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
                  << ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), _M_pts.size2() )
                  << "\n";
#endif
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( nOrder ), _M_pts.size2() );
    }


    /**
     * evaluate the Legendre polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    static matrix_type evaluate( points_type const& __pts )
    {
        return evaluate( __pts, mpl::int_<nDim>() );
    }

    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts )
    {
        return derivate( __pts, mpl::int_<nDim>() );
    }


    //@}

private:
private:

    static value_type normalization( int i ) { return (is_normalized?math::sqrt( value_type( i ) + 0.5 ) : value_type(1)); }
    static value_type normalization( int i, int j ) { return (is_normalized?math::sqrt( ( value_type( i ) + 0.5 ) *
                                                                                        ( value_type( j ) + 0.5 ) ) : value_type(1)); }
    static value_type normalization( int i, int j, int k ) { return (is_normalized?math::sqrt( ( value_type( i ) + 0.5 ) *
                                                                                               ( value_type( j ) + 0.5 ) *
                                                                                               ( value_type( k ) + 0.5 ) ) : value_type(1)); }
    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type
    evaluate( points_type const& __pts, mpl::int_<1> )
    {
        matrix_type m ( JacobiBatchEvaluation<nOrder,value_type>( 0.0, 0.0, ublas::row(__pts, 0) ) );
        if ( is_normalized )
            {
                for ( uint16_type i = 0;i < m.size1(); ++i )
                    ublas::row( m, i ) *= normalization( i );
            }
        return m;
    }

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    static vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<1> )
    {
        LIFE_ASSERT( __pts().size1() == 1 )( __pts().size1() )( __pts().size2() ).error("invalid points");
        // Debug() << "Expansion::derivate<1>] number of points " << __pts().size2() << "\n";

        vector_matrix_type D( 1 );
        D[0].resize( nOrder+1, __pts().size2() );
        D[0] = JacobiBatchDerivation<nOrder,value_type>( 0.0, 0.0, ublas::row(__pts(),0) );
        if ( is_normalized )
            for ( uint16_type i = 0; i < nOrder+1; ++i )
                ublas::row( D[0], i ) *= normalization( i );
        return D;
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type evaluate( points_type const& __pts, mpl::int_<2> );

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> );

    /**
     * Evaluation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    static matrix_type evaluate( points_type const& __pts, mpl::int_<3> );

    /**
     * derivation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> );

private:
    reference_convex_type _M_refconvex;
    points_type _M_pts;

};
template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<2> )
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    ublas::vector<value_type> eta1s = ublas::row( __pts, 0 );
    ublas::vector<value_type> eta2s = ublas::row( __pts, 1 );

    matrix_type as( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, eta1s ) );
    matrix_type bs( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, eta2s ) );

    for ( uint16_type cur = 0, i = 0; i < nOrder+1;++i )
        {
            for ( uint16_type j = 0; j < nOrder+1; ++j,++cur )
                {
                    ublas::row( res, cur ) = normalization( i, j ) * ublas::element_prod( ublas::row( as, i),
                                                                                          ublas::row( bs, j) );
                }
        }
    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<2> )
{
    vector_matrix_type res( 2 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    // evaluate Legendre polynomials components
    matrix_type as( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 0 ) ) );
    matrix_type das( JacobiBatchDerivation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 0 ) ) );
    matrix_type bs( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 1 ) ) );
    matrix_type dbs( JacobiBatchDerivation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 1 ) ) );

    for ( uint16_type cur = 0, i = 0; i < nOrder+1;++i )
        {
            for ( uint16_type j = 0; j < nOrder+1; ++j,++cur )
                {
                    ublas::row( res[0], cur ) = normalization( i, j ) * ublas::element_prod( ublas::row( das, i),
                                                                                             ublas::row( bs, j) );
                    ublas::row( res[1], cur ) = normalization( i, j ) * ublas::element_prod( ublas::row( as, i),
                                                                                             ublas::row( dbs, j) );
                }
        }

    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, mpl::int_<3> )
{
    matrix_type res( convex_type::polyDims( nOrder ), __pts.size2() );

    ublas::vector<value_type> eta1s = ublas::row( __pts, 0 );
    ublas::vector<value_type> eta2s = ublas::row( __pts, 1 );
    ublas::vector<value_type> eta3s = ublas::row( __pts, 2 );

    matrix_type as( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, eta1s ) );
    matrix_type bs( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, eta2s ) );
    matrix_type cs( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, eta3s ) );

    for ( uint16_type cur = 0, i = 0; i < nOrder+1;++i )
        {
            for ( uint16_type j = 0; j < nOrder+1; ++j )
                {
                    for ( uint16_type k = 0; k < nOrder+1; ++k,++cur )
                        {
                            ublas::row( res, cur ) =
                                normalization( i, j, k )*
                                ublas::element_prod(ublas::element_prod( ublas::row( as, i),
                                                                         ublas::row( bs, j) ),
                                                    ublas::row( cs, k ) );
                        }
                }
        }
    return res;
}

template<uint16_type Dim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Legendre<Dim, Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, mpl::int_<3> )
{
    vector_matrix_type res( 3 );
    res[0].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[1].resize( convex_type::polyDims( nOrder ), __pts().size2() );
    res[2].resize( convex_type::polyDims( nOrder ), __pts().size2() );

    // evaluate Legendre polynomials components
    matrix_type as( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 0 ) ) );
    matrix_type das( JacobiBatchDerivation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 0 ) ) );
    matrix_type bs( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 1 ) ) );
    matrix_type dbs( JacobiBatchDerivation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 1 ) ) );
    matrix_type cs( JacobiBatchEvaluation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 2 ) ) );
    matrix_type dcs( JacobiBatchDerivation<nOrder, value_type>( 0.0, 0.0, ublas::row( __pts(), 2 ) ) );

    for ( uint16_type cur = 0, i = 0; i < nOrder+1;++i )
        {
            for ( uint16_type j = 0; j < nOrder+1; ++j )
                {
                    for ( uint16_type k = 0; k < nOrder+1; ++k,++cur )
                        {
                            ublas::row( res[0], cur ) = ( normalization( i, j, k ) *
                                                          ublas::element_prod( ublas::element_prod( ublas::row( das, i),
                                                                                                    ublas::row( bs, j) ),
                                                                               ublas::row( cs, k) ) );

                            ublas::row( res[1], cur ) = ( normalization( i, j, k ) *
                                                          ublas::element_prod( ublas::element_prod( ublas::row( as, i),
                                                                                                    ublas::row( dbs, j) ),
                                                                               ublas::row( cs, k) ) );

                            ublas::row( res[2], cur ) = ( normalization( i, j, k ) *
                                                          ublas::element_prod( ublas::element_prod( ublas::row( as, i),
                                                                                                    ublas::row( bs, j) ),
                                                                               ublas::row( dcs, k) ) );


                        }
                }
        }

    return res;
}


}
#endif /* __Legendre_H */
