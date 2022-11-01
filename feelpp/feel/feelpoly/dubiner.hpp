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
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelpoly/expansions.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelmesh/pointset.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/warpblend.hpp>
#include <feel/feelpoly/expansiontypes.hpp>

namespace Feel
{

template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
class Dubiner;


template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct DubinerTraits
{
    inline static const uint16_type nDim = Dim;
    inline static const uint16_type nRealDim = RealDim;
    inline static const int nOrder = Degree;
    inline static const bool is_order_dynamic = ( Degree == Dynamic );
    inline static const int nConvexOrderDiff = (is_order_dynamic)?Dynamic:nDim+nOrder+1;
    inline static const bool is_normalized = NormalizationPolicy::is_normalized;

    /** @name Typedefs
     */
    //@{

    /*
     * numerical type
     */
    typedef T value_type;

    template<int order, typename V = value_type>
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

    template<int NewOrder>
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

    using diff_pointset_type = boost::mp11::mp_if_c<nDim==2,
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
         int Degree = Dynamic,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Dubiner

{
    
public:
    typedef DubinerTraits<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> traits_type;

    inline static const int nDim = traits_type::nDim;
    inline static const int nRealDim = traits_type::nRealDim;
    inline static const int nOrder = traits_type::nOrder;
    inline static const int nConvexOrderDiff = traits_type::nConvexOrderDiff;
    inline static const bool is_order_dynamic = ( Degree == Dynamic );
    inline static const bool is_normalized = traits_type::is_normalized;
    inline static const bool isTransformationEquivalent = true;
    inline static const bool isContinuous = false;
    inline static const bool is_product = true;
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

    /**
     * @brief Construct a new Dubiner polynomial of dynamic degree \p Degree
     * 
     * @param degree polynomial degree
     */
    explicit Dubiner( int degree );
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
    size_type size() const noexcept 
    {
        if constexpr ( is_order_dynamic )
            return convex_type::polyDims( M_order );
        else
            return convex_type::polyDims( nOrder );

    }

    /**
     * \return the maximum degree of the Dubiner polynomial to be
     * constructed
     */
    uint16_type degree() const noexcept
    {
        return M_order;
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
                  << ublas::identity_matrix<value_type>( reference_convex_type::polyDims( M_order ), M_pts.size2() )
                  << "\n";
#endif
        return ublas::identity_matrix<value_type>( reference_convex_type::polyDims( M_order ), M_pts.size2() );
    }


    /**
     * @brief evaluate the Dubiner polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    matrix_type evaluate( points_type const& __pts ) const
    {
        return evaluate( __pts, is_order_dynamic?M_order:nOrder ); 
    }


    /**
     * evaluate the Dubiner polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    matrix_type evaluate( points_type const& __pts, int order  ) const
    {
         if constexpr ( is_order_dynamic )
            return evaluate( __pts, M_order, mpl::int_<nDim>() );
        else
            return evaluate( __pts, nOrder, mpl::int_<nDim>() );
    }

    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts )
    {
        //static_assert(!is_order_dynamic, "Polynomial order must be known at compile time, it cannot be Dynamic");
        return derivate( __pts, nOrder, mpl::int_<nDim>() );
    }

    template<typename AE>
    vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts, int order = Dynamic )
    {
        if ( is_order_dynamic )
            return derivate( __pts, M_order, mpl::int_<nDim>() );
        else
            return derivate( __pts, nOrder, mpl::int_<nDim>() );
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
    evaluate( points_type const& __pts, int /*order*/, mpl::int_<0> )
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
    evaluate( points_type const& __pts, int order, mpl::int_<1> )
    {
        matrix_type m ( JacobiBatchEvaluation<value_type>( order, 0.0, 0.0, ublas::row( __pts, 0 ) ) );

        if constexpr ( is_normalized )
        {
            for ( uint16_type i = 0; i < m.size1(); ++i )
                ublas::row( m, i ) *= math::sqrt( value_type( i )+0.5 );
        }

        return m;
    }

    template<typename AE>
    static vector_matrix_type
    derivate( ublas::matrix_expression<AE> const& __pts, int /*order*/, mpl::int_<0> )
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
    derivate( ublas::matrix_expression<AE> const& __pts, int order, mpl::int_<1> )
    {
        vector_matrix_type D( 1 );
        D[0].resize( order+1, __pts().size2() );
        D[0] = JacobiBatchDerivation<value_type>( order, 0.0, 0.0, ublas::row( __pts(),0 ) );

        if constexpr ( is_normalized )
            for ( uint16_type i = 0; i < order+1; ++i )
                ublas::row( D[0], i ) *= math::sqrt( value_type( i )+0.5 );

        return D;
    }

    /**
     * Evaluation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    static matrix_type evaluate( points_type const& __pts, int order, mpl::int_<2> );

    /**
     * derivation at a set of points of the expansion basis in 2D on
     * the triangle
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, int order, mpl::int_<2> );

    /**
     * Evaluation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    static matrix_type evaluate( points_type const& __pts, int order, mpl::int_<3> );

    /**
     * derivation at a set of points of the expansion basis in 3D on
     * the tetrahedron
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE> const& __pts, int order, mpl::int_<3> );

private:

    int M_order;
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
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::Dubiner()
    :
    Dubiner( Degree )
{}

template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::Dubiner(int order)
    :
    M_order( order ),
    M_refconvex(order),
    M_pts( M_refconvex.makePoints( Dim, 0 ) ),
    M_D( Dim )
{
    reference_convex_type refconvex( order );
    // constructor pointset for differentiation only in
    // the interior(1)
    DVLOG( 4 ) << fmt::format( "[Dubiner({})] pts:{}\n", order, em( M_pts ) );
    diff_pointset_type  diff_pts(nDim+M_order+1, 1 );
    em_cmatrix_col_type<T> Ae = em(diff_pts.points());
    DVLOG( 4 ) << fmt::format( "diff_pts: {}\n, diff_order: {}\n ", Ae, nDim + M_order + 1 );
    matrix_type A( this->evaluate( diff_pts.points(), order ) );
    DVLOG( 4 ) << fmt::format( "A: {}\n ", em_matrix_row_type<T>(A.data().begin(),A.size1(),A.size2()) );
    matrix_type D = ublas::identity_matrix<value_type>( A.size1(), A.size2()  );
    LU<matrix_type> lu( A );
    matrix_type C = lu.solve( D );
    
    vector_matrix_type d ( this->derivate( diff_pts.points(), order ) );
    for ( size_type i = 0; i < d.size(); ++i )
    {
        M_D[i] = ublas::prod( d[i], C );
        glas::clean( M_D[i] );
    }
    
}

template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, int order, mpl::int_<2> )
{
    DVLOG( 4 ) << fmt::format( "[Dubiner::evaluate] order: {} compile time order: {}, polynomialset size:{} \n", order, Degree, convex_type::polyDims( order ) );
    matrix_type res( convex_type::polyDims( order ), __pts.size2() );

    details::etas<TRIANGLE, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );

    //std::cout << "xis = " << __pts << "\n";
    //std::cout << "etas = " << etas() << "\n";

    matrix_type as( JacobiBatchEvaluation<value_type>( order, 0.0, 0.0, eta1s ) );
    std::vector<matrix_type> bs( order+1 );

    for ( int i = 0; i < order+1; ++i )
    {
        bs[ i ].resize( order-i, eta2s.size() );
        bs[ i ] = dyna::JacobiBatchEvaluation( order-i, value_type( 2*i+1 ), value_type( 0.0 ), eta2s );
    }


    details::scalings<nOrder, T> scalings( order, eta2s );


    for ( uint16_type cur = 0, k = 0; k < order+1; ++k )
    {
        for ( uint16_type i = 0; i < k+1; ++i,++cur )
        {
            uint16_type ii = k-i;
            uint16_type jj = i;

            if ( is_normalized )
            {
                value_type normalization = math::sqrt( ( value_type( ii )+0.5 )*( value_type( ii+jj )+1.0 ) );

                for ( uint16_type l = 0; l < as.size2(); ++l )
                    res( cur, l ) = normalization*as( ii,l )*scalings()( ii,l )*bs[ii]( jj,l );
            }

            else
            {
                for ( uint16_type l = 0; l < as.size2(); ++l )
                    res( cur, l ) = as( ii,l )*scalings()( ii,l )*bs[ii]( jj,l );
            }
        }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Dubiner<Dim, RealDim,  Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, int order, mpl::int_<2> )
{
    vector_matrix_type res( 2 );
    res[0].resize( convex_type::polyDims( order ), __pts().size2() );
    res[1].resize( convex_type::polyDims( order ), __pts().size2() );

    // transform wrapped coordinates in cartesian coordinates
    details::etas<TRIANGLE, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );

    //std::cout << "xis = " << __pts << "\n";
    //std::cout << "etas = " << etas() << "\n";

    // evaluate Dubiner polynomials components
    matrix_type as(  JacobiBatchEvaluation<value_type>( order, 0.0, 0.0, eta1s ) );
    matrix_type das( JacobiBatchDerivation<value_type>( order, 0.0, 0.0, eta1s ) );
    //std::cout << "das= " << das <<  "\n";
    std::vector<matrix_type> bs( order+1 );
    std::vector<matrix_type> dbs( order+1 );

    for ( uint16_type i = 0; i < order+1; ++i )
    {
        bs[ i ].resize( order-i, eta2s.size() );
        dbs[ i ].resize( order-i, eta2s.size() );
        bs[ i ] = dyna::JacobiBatchEvaluation( order-i, value_type( 2*i+1 ), value_type( 0.0 ), eta2s );
        dbs[ i ] = dyna::JacobiBatchDerivation( order-i, value_type( 2*i+1 ), value_type( 0.0 ), eta2s );

        //std::cout << "dbs["<< i << "]= " << dbs[i] <<  "\n";
    }

    details::scalings<nOrder, T> scalings( order, eta2s );
    //std::cout << "scalings = " << scalings() << "\n";
    ublas::vector<value_type> one( ublas::scalar_vector<value_type>( eta1s.size(), 1.0 ) );
    ublas::vector<value_type> tmp( ublas::scalar_vector<value_type>( eta1s.size(), 1.0 ) );

    // assemble Dubiner polynomials components
    for ( uint16_type k = 0, cur = 0; k < order+1; ++k )
    {
        for ( uint16_type i = 0; i < k+1; ++i, ++cur )
        {
            uint16_type ii = k-i;
            uint16_type jj = i;


            // x derivation
            ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                        ublas::row( bs[ii], jj ) );

            if ( ii > 0 )
                ublas::row( res[0], cur ) = element_prod( ublas::row( res[0], cur ),
                                            ublas::row( scalings(), ii-1 ) );

            // y derivation
            ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                        ublas::row( bs[ii], jj ) );
            ublas::row( res[1], cur ) = 0.5 * element_prod( ublas::row( res[1], cur ), ( one+eta1s ) );

            if ( ii > 0 )
                ublas::row( res[1], cur ) = element_prod( ublas::row( res[1], cur ),
                                            ublas::row( scalings(), ii-1 ) );

            // derivate (1-x)^ii
            tmp = ublas::element_prod( ublas::row( scalings(), ii ),
                                       ublas::row( dbs[ii], jj ) );

            if ( ii > 0 )
                tmp -= 0.5 * ii * ublas::element_prod( ublas::row( scalings(), ii-1 ),
                                                       ublas::row( bs[ii], jj ) );

            // add contrib to y derivation
            ublas::row( res[1], cur ) += ublas::element_prod( ublas::row( as, ii ), tmp );

            // orthonormalize if required
            if ( is_normalized )
            {
                value_type normalization = math::sqrt( ( value_type( ii )+0.5 )*( value_type( ii+jj )+1.0 ) );
                ublas::row( res[0], cur ) *= normalization;
                ublas::row( res[1], cur ) *= normalization;
            }
        }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::evaluate( points_type const& __pts, int order, mpl::int_<3> )
{
    matrix_type res( convex_type::polyDims( order ), __pts.size2() );

    DCHECK( __pts.size1() == 3 ) << "invalid space dimension for point set, should be 3, it is currently " << __pts.size1();

    details::etas<TETRAHEDRON, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );
    ublas::vector<value_type> eta3s = ublas::row( etas(), 2 );

    //std::cout << "xis = " << __pts << "\n";
    //std::cout << "etas = " << etas() << "\n";

    matrix_type as( JacobiBatchEvaluation<value_type>( order, 0.0, 0.0, eta1s ) );
    std::vector<matrix_type> bs( order+1 );
    ublas::matrix<matrix_type> cs( order+1, order+1 );

    for ( int i = 0; i < order+1; ++i )
    {
        bs[ i ].resize( order-i, eta2s.size() );
        bs[ i ] = dyna::JacobiBatchEvaluation( order-i, value_type( 2*i+1 ), value_type( 0.0 ), eta2s );

        for ( int j = 0; j < order+1-i; ++j )
        {
            cs( i, j ).resize( order-i-j, eta3s.size() );
            cs( i, j ) = dyna::JacobiBatchEvaluation( order-i-j,
                         value_type( 2*( i+j+1 ) ), value_type( 0.0 ), eta3s );
        }
    }


    details::scalings<nOrder, T> scalings2( order, eta2s );
    details::scalings<nOrder, T> scalings3( order, eta3s );

    for ( uint16_type cur = 0, k = 0; k < order+1; ++k )
    {
        for ( uint16_type i = 0; i < k+1; ++i )
        {
            for ( uint16_type j = 0; j < k+1-i; ++j,++cur )
            {
                uint16_type ii = k-i-j;
                uint16_type jj = j;
                uint16_type kk = i;


                if ( is_normalized )
                {
                    value_type normalization = math::sqrt( ( value_type( ii )+0.5 )*
                                                           ( value_type( ii+jj )+1.0 )*
                                                           ( value_type( ii+jj+kk )+1.5 ) );

                    for ( uint16_type l = 0; l < as.size2(); ++l )
                        res( cur, l ) = normalization*( as( ii,l )*
                                                        scalings2()( ii,l )*bs[ii]( jj,l )*
                                                        scalings3()( ii+jj,l )*cs( ii, jj )( kk,l ) );
                }

                else
                {
                    for ( uint16_type l = 0; l < as.size2(); ++l )
                        res( cur, l ) = as( ii,l )*
                                        scalings2()( ii,l )*bs[ii]( jj,l )*
                                        scalings3()( ii+jj,l )*cs( ii, jj )( kk,l );
                }
            }
        }
    }

    return res;
}

template<uint16_type Dim,
         uint16_type RealDim,
         int Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
template<typename AE>
typename Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::vector_matrix_type
Dubiner<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>::derivate( ublas::matrix_expression<AE> const& __pts, int order, mpl::int_<3> )
{
    vector_matrix_type res( 3 );
    res[0].resize( convex_type::polyDims( order ), __pts().size2() );
    res[1].resize( convex_type::polyDims( order ), __pts().size2() );
    res[2].resize( convex_type::polyDims( order ), __pts().size2() );

    FEELPP_ASSERT( __pts().size1() == 3 )( __pts().size1() ).error( "invalid space dimension" );

    details::etas<TETRAHEDRON, value_type> etas( __pts );
    ublas::vector<value_type> eta1s = ublas::row( etas(), 0 );
    ublas::vector<value_type> eta2s = ublas::row( etas(), 1 );
    ublas::vector<value_type> eta3s = ublas::row( etas(), 2 );

    //std::cout << "xis = " << __pts << "\n";
    //std::cout << "etas = " << etas() << "\n";

    matrix_type as(  JacobiBatchEvaluation<value_type>( order, 0.0, 0.0, eta1s ) );
    matrix_type das( JacobiBatchDerivation<value_type>( order, 0.0, 0.0, eta1s ) );
    std::vector<matrix_type> bs( order+1 );
    std::vector<matrix_type> dbs( order+1 );
    ublas::matrix<matrix_type> cs( order+1, order+1 );
    ublas::matrix<matrix_type> dcs( order+1, order+1 );

    for ( int i = 0; i < order+1; ++i )
    {
        bs[ i ].resize( order-i, eta2s.size() );
        dbs[ i ].resize( order-i, eta2s.size() );
        bs[ i ] = dyna::JacobiBatchEvaluation( order-i, value_type( 2*i+1 ), value_type( 0.0 ), eta2s );
        dbs[ i ] = dyna::JacobiBatchDerivation( order-i, value_type( 2*i+1 ), value_type( 0.0 ), eta2s );

        for ( int j = 0; j < order+1-i; ++j )
        {
            cs( i, j ).resize( order-i-j, eta3s.size() );
            dcs( i, j ).resize( order-i-j, eta3s.size() );
            cs( i, j ) = dyna::JacobiBatchEvaluation( order-i-j,
                         value_type( 2*( i+j+1 ) ), value_type( 0.0 ), eta3s );
            dcs( i, j ) = dyna::JacobiBatchDerivation( order-i-j,
                          value_type( 2*( i+j+1 ) ), value_type( 0.0 ), eta3s );
        }
    }


    details::scalings<nOrder, T> scalings2( order, eta2s );
    details::scalings<nOrder, T> scalings3( order, eta3s );

    //std::cout << "scalings = " << scalings() << "\n";
    ublas::vector<value_type> one( ublas::scalar_vector<value_type>( eta1s.size(), 1.0 ) );
    ublas::vector<value_type> tmp( ublas::scalar_vector<value_type>( eta1s.size(), 1.0 ) );


    for ( uint16_type cur = 0, k = 0; k < order+1; ++k )
    {
        for ( uint16_type i = 0; i < k+1; ++i )
        {
            for ( uint16_type j = 0; j < k+1-i; ++j,++cur )
            {
                uint16_type ii = k-i-j;
                uint16_type jj = j;
                uint16_type kk = i;


                // x derivation
                ublas::row( res[0], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                            ublas::row( bs[ii], jj ) );
                ublas::row( res[0], cur ) = element_prod( ublas::row( res[0], cur ),
                                            ublas::row( cs( ii, jj ), kk ) );

                if ( ii > 0 )
                    ublas::row( res[0], cur ) = element_prod( ublas::row( res[0], cur ),
                                                ublas::row( scalings2(), ii-1 ) );

                if ( ii+jj > 0 )
                    ublas::row( res[0], cur ) = element_prod( ublas::row( res[0], cur ),
                                                ublas::row( scalings3(), ii+jj-1 ) );

                // y derivation
                ublas::row( res[1], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                            ublas::row( bs[ii], jj ) );
                ublas::row( res[1], cur ) = element_prod( ublas::row( res[1], cur ),
                                            ublas::row( cs( ii, jj ), kk ) );
                ublas::row( res[1], cur ) = 0.5 * element_prod( ublas::row( res[1], cur ),
                                            ( one+eta1s ) );

                if ( ii > 0 )
                    ublas::row( res[1], cur ) = element_prod( ublas::row( res[1], cur ),
                                                ublas::row( scalings2(), ii-1 ) );

                if ( ii+jj > 0 )
                    ublas::row( res[1], cur ) = element_prod( ublas::row( res[1], cur ),
                                                ublas::row( scalings3(), ii+jj-1 ) );

                // derivate (1-x)^ii
                tmp = ublas::element_prod( ublas::row( scalings2(), ii ),
                                           ublas::row( dbs[ii], jj ) );

                if ( ii > 0 )
                    tmp -= 0.5 * ii * ublas::element_prod( ublas::row( scalings2(), ii-1 ),
                                                           ublas::row( bs[ii], jj ) );

                tmp = ublas::element_prod( tmp,
                                           ublas::row( as, ii ) );
                tmp = ublas::element_prod( tmp,
                                           ublas::row( cs( ii, jj ), kk ) );

                if ( ii+jj > 0 )
                    tmp = ublas::element_prod( tmp,
                                               ublas::row( scalings3(), ii+jj-1 ) );

                // add contrib to y derivation
                ublas::row( res[1], cur ) += tmp;

                // z derivation
                ublas::row( res[2], cur ) = ublas::element_prod( ublas::row( das, ii ),
                                            ublas::row( bs[ii], jj ) );
                ublas::row( res[2], cur ) = element_prod( ublas::row( res[2], cur ),
                                            ublas::row( cs( ii, jj ), kk ) );
                ublas::row( res[2], cur ) = 0.5 * element_prod( ublas::row( res[2], cur ),
                                            ( one+eta1s ) );

                if ( ii > 0 )
                    ublas::row( res[2], cur ) = element_prod( ublas::row( res[2], cur ),
                                                ublas::row( scalings2(), ii-1 ) );

                if ( ii+jj > 0 )
                    ublas::row( res[2], cur ) = element_prod( ublas::row( res[2], cur ),
                                                ublas::row( scalings3(), ii+jj-1 ) );

                // derivate (1-x)^ii
                tmp = ublas::element_prod( ublas::row( scalings2(), ii ),
                                           ublas::row( dbs[ii], jj ) );

                if ( ii > 0 )
                    tmp -= 0.5 * ii * ublas::element_prod( ublas::row( scalings2(), ii-1 ),
                                                           ublas::row( bs[ii], jj ) );

                tmp = ublas::element_prod( tmp,
                                           ublas::row( as, ii ) );
                tmp = ublas::element_prod( tmp,
                                           ublas::row( cs( ii, jj ), kk ) );
                tmp = 0.5 * element_prod( tmp, ( one+eta2s ) );

                if ( ii+jj > 0 )
                    tmp = ublas::element_prod( tmp,
                                               ublas::row( scalings3(), ii+jj-1 ) );

                // add contrib to z derivation
                ublas::row( res[2], cur ) += tmp;

                // derivate (1-x)^ii
                tmp = ublas::element_prod( ublas::row( scalings3(), ii+jj ),
                                           ublas::row( dcs( ii, jj ), kk ) );

                if ( ii+jj > 0 )
                    tmp -= 0.5*( ii+jj )*ublas::element_prod( ublas::row( cs( ii, jj ), kk ),
                            ublas::row( scalings3(), ii+jj-1 ) );

                tmp = ublas::element_prod( tmp,
                                           ublas::row( as, ii ) );
                tmp = ublas::element_prod( tmp,
                                           ublas::row( bs[ ii ], jj ) );
                tmp = ublas::element_prod( tmp,
                                           ublas::row( scalings2(), ii ) );

                // add contrib to z derivation
                ublas::row( res[2], cur ) += tmp;

                if ( is_normalized )
                {
                    value_type normalization = math::sqrt( ( value_type( ii )+0.5 )*
                                                           ( value_type( ii+jj )+1.0 )*
                                                           ( value_type( ii+jj+kk )+1.5 ) );

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
