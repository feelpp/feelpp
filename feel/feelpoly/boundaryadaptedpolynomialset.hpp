/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file boundaryadaptedpolynomialset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-30
 */
#ifndef __BoundaryAdaptedPolynomialSet_H
#define __BoundaryAdaptedPolynomialSet_H 1

namespace Feel
{
/** \cond DETAIL */
namespace detail
{
/**
 * \internal
 * \class BoundaryAdaptedPolynomialSet
 * \brief a set of boundary adapted polynomials over a convex
 *
 */

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType = Scalar,
         typename T = double,
         template<uint16_type,uint16_type,uint16_type> class Convex = Simplex>
class BoundaryAdaptedPolynomialSet
{};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class BoundaryAdaptedPolynomialSet<Dim, Order, PolySetType, T, Simplex>
    :
public PolynomialSet<BoundaryAdapted<Dim, Order, T>, PolySetType >
{
    typedef PolynomialSet<BoundaryAdapted<Dim, Order, T>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    typedef BoundaryAdaptedPolynomialSet<Dim, Order, PolySetType, T, Simplex> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;

    static const bool is_continuous = true;

    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef BoundaryAdapted<Dim, Order, T> basis_type;
    typedef typename basis_type::points_type points_type;

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

    BoundaryAdaptedPolynomialSet()
        :
        super( basis_type() )
    {

        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );

        if ( is_tensor2 )
            std::cout << "[boundaryadaptedpolynomialset] m = " << m << "\n";

        this->setCoefficient( polyset_type::toType( m ), true );
    }

    BoundaryAdaptedPolynomialSet<Dim, Order, Scalar,T, Simplex > toScalar() const
    {
        return BoundaryAdaptedPolynomialSet<Dim, Order, Scalar,T, Simplex >();
    }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const
    {
        return "dubinerba";
    }

    points_type const& points() const
    {
        return this->basis().points();
    }
    points_type const& points( int f ) const
    {
        return this->basis().points( f );
    }
};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
const uint16_type BoundaryAdaptedPolynomialSet<Dim, Order, PolySetType,T, Simplex >::nDofPerEdge;

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
const uint16_type BoundaryAdaptedPolynomialSet<Dim, Order, PolySetType,T, Simplex >::nLocalDof;


template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class BoundaryAdaptedPolynomialSet<Dim, Order, PolySetType, T, Hypercube>
    :
public PolynomialSet<TensorisedBoundaryAdapted<Dim, Order, T>, PolySetType >
{
    typedef PolynomialSet<TensorisedBoundaryAdapted<Dim, Order, T>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    typedef BoundaryAdaptedPolynomialSet<Dim, Order, PolySetType, T, Hypercube> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_continuous = true;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef TensorisedBoundaryAdapted<Dim, Order, T> basis_type;
    typedef Hypercube<Dim, Order, Dim> convex_type;
    template<int O>
    struct convex
    {
        typedef Hypercube<Dim, O, Dim> type;
    };
    typedef Reference<convex_type, nDim, nOrder, nDim, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    typedef typename basis_type::points_type points_type;


    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;


    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    /** mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
        mpl::identity<typename mpl::max<mpl::int_<0>,mpl::int_<(nOrder-1)*(nOrder-1)*(nOrder-1)> >::type >,
        mpl::identity<mpl::int_<0> >::type::type::value;
    **/

    //!< Total number of degrees of freedom
    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

    BoundaryAdaptedPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );

        if ( is_tensor2 )
            std::cout << "[boundaryadaptedpolynomialset] m = " << m << "\n";

        this->setCoefficient( polyset_type::toType( m ), true );
    }

    BoundaryAdaptedPolynomialSet<Dim, Order, Scalar,T, Hypercube > toScalar() const
    {
        return BoundaryAdaptedPolynomialSet<Dim, Order, Scalar,T, Hypercube >();
    }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const
    {
        return "tensorizedba";
    }

    points_type const& points() const
    {
        return this->basis().points();
    }
    points_type const& points( int f ) const
    {
        return this->basis().points( f );
    }
};
} // detail

/// \encond 

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar>
class BoundaryAdaptedPolynomialSet
{
public:
    template<uint16_type N,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                mpl::identity<Feel::detail::BoundaryAdaptedPolynomialSet<N,Order,PolySetType,T,Simplex> >,
                mpl::identity<Feel::detail::
BoundaryAdaptedPolynomialSet<N,Order,PolySetType,T,Hypercube> > >::type::type result_type;
        typedef result_type type;
    };

    typedef BoundaryAdaptedPolynomialSet<Order,Scalar> component_basis_type;
};

}
#endif /* __BoundaryAdaptedPolynomialSet_H */
