/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-30

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
/**
   \file orthonormalpolynomialset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-30
 */
#ifndef __OrthonormalPolynomialSet_H
#define __OrthonormalPolynomialSet_H 1
namespace Life
{
/** \cond DETAIL */
namespace detail
{
/**
 * \internal
 * \class OrthonormalPolynomialSet
 * \brief a set of orthonormal polynomials over a convex
 *
 * On the simplicies we use the Dubiner basis
 *
 */
template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType = Scalar,
         typename T = double,
         template<uint16_type,uint16_type,uint16_type> class Convex = Simplex>
class OrthonormalPolynomialSet
{};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class OrthonormalPolynomialSet<Dim, Order, PolySetType, T, Simplex>
    :
    public PolynomialSet<Dubiner<Dim, Order, Normalized<true>, T, StorageUBlas>, PolySetType >
{
    typedef PolynomialSet<Dubiner<Dim, Order, Normalized<true>, T, StorageUBlas>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    typedef OrthonormalPolynomialSet<Dim, Order, PolySetType, T, Simplex> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef Dubiner<Dim, Order, Normalized<true>, T, StorageUBlas> basis_type;
    typedef Simplex<Dim, Order, Dim> convex_type;
    template<int O>
    struct convex
    {
        typedef Simplex<Dim, O, Dim> type;
    };
    typedef Reference<convex_type, nDim, nOrder, nDim, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;
    typedef typename matrix_node<value_type>::type points_type;

    OrthonormalPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
        if ( is_tensor2 )
            std::cout << "[orthonormalpolynomialset] m = " << m << "\n";
        if ( !(ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                      m ) < 1e-10 ) )
            std::cout << "m1=" << m << "\n"
                      << "m2=" << polyset_type::toMatrix( polyset_type::toType( m ) ) << "\n"
                      << ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) - m ) << "\n";
        LIFE_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                             m ) < 1e-10 )( m ).warn ( "invalid transformation" );
        this->setCoefficient( polyset_type::toType( m ), true );
    }

    OrthonormalPolynomialSet<Dim, Order, Scalar,T, Simplex > toScalar() const
    {
        return OrthonormalPolynomialSet<Dim, Order, Scalar,T, Simplex >();
    }

    /**
     * \return the family name of the polynomial set
     */
    std::string familyName() const { return "dubiner"; }


    points_type points() const { return points_type(); }
    points_type points( int f ) const { return points_type(); }
};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
const uint16_type OrthonormalPolynomialSet<Dim, Order,PolySetType,T, Simplex>::nLocalDof;


template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class OrthonormalPolynomialSet<Dim, Order, PolySetType, T, SimplexProduct>
    :
    public PolynomialSet<Legendre<Dim, Order, Normalized<true>, T>, PolySetType >
{
    typedef PolynomialSet<Legendre<Dim, Order, Normalized<true>, T>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    typedef OrthonormalPolynomialSet<Dim, Order, PolySetType, T, SimplexProduct> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    typedef typename super::component_type component_type;
    typedef T value_type;
    typedef Legendre<Dim, Order, Normalized<true>, T> basis_type;
    typedef SimplexProduct<Dim, Order, Dim> convex_type;
    typedef typename matrix_node<value_type>::type points_type;

    template<int O>
    struct convex
    {
        typedef SimplexProduct<Dim, O, Dim> type;
    };
    typedef Reference<convex_type, nDim, nOrder, nDim, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

    OrthonormalPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
        if ( is_tensor2 )
            std::cout << "[orthonormalpolynomialset] m = " << m << "\n";
        LIFE_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                             m ) < 1e-10 )( m ).warn ( "invalid transformation" );
        this->setCoefficient( polyset_type::toType( m ), true );
    }

    OrthonormalPolynomialSet<Dim, Order, Scalar,T, SimplexProduct > toScalar() const
    {
        return OrthonormalPolynomialSet<Dim, Order, Scalar,T, SimplexProduct >();
    }
    std::string familyName() const { return "legendre"; }
    points_type points() const { return points_type(); }
    points_type points( int f ) const { return points_type(); }
};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
const uint16_type OrthonormalPolynomialSet<Dim, Order,PolySetType,T, SimplexProduct>::nLocalDof;
} // detail
/** \encond  */

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar>
class OrthonormalPolynomialSet
{
public:
    template<uint16_type N,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                                  mpl::identity<detail::OrthonormalPolynomialSet<N,Order,PolySetType,T,Simplex> >,
                                  mpl::identity<detail::OrthonormalPolynomialSet<N,Order,PolySetType,T,SimplexProduct> > >::type::type result_type;
        typedef result_type type;
    };

    typedef OrthonormalPolynomialSet<Order,Scalar> component_basis_type;
};

} // Life
#endif /* __OrthonormalPolynomialSet_H */
