/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file orthogonalpolynomialset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-30
 */
#ifndef __OrthogonalPolynomialSet_H
#define __OrthogonalPolynomialSet_H 1

namespace Feel
{
/**
 * \class OrthogonalPolynomialSet
 * \brief a set of orthogonal polynomials over a convex
 *
 * On the simplicies we use the Dubiner basis
 *
 */
template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType = Scalar,
         typename T = double,
         template<uint16_type,uint16_type,uint16_type> class Convex = Simplex>
class OrthogonalPolynomialSet
{};


template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class OrthogonalPolynomialSet<Dim, Order, PolySetType, T, Simplex>
    :
public PolynomialSet<Dubiner<Dim, Order, Normalized<false>, T, StorageUBlas>, PolySetType >
{
    typedef PolynomialSet<Dubiner<Dim, Order, Normalized<false>, T, StorageUBlas>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const bool isTransformationEquivalent = true;

    typedef OrthogonalPolynomialSet<Dim, Order, PolySetType, T, Simplex> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef Dubiner<Dim, Order, Normalized<false>, T, StorageUBlas> basis_type;
    typedef Simplex<Dim, Order, Dim> convex_type;
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

    OrthogonalPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
        this->setCoefficient( polyset_type::toType( m ), true );
    }

    OrthogonalPolynomialSet<Dim, Order, Scalar,T, Simplex > toScalar() const
    {
        return OrthogonalPolynomialSet<Dim, Order, Scalar,T, Simplex >();
    }

    std::string familyName() const
    {
        return "dubiner";
    }
};
template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
const uint16_type OrthogonalPolynomialSet<Dim, Order,PolySetType,T, Simplex>::nLocalDof;

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
class OrthogonalPolynomialSet<Dim, Order, PolySetType, T, Hypercube>
    :
public PolynomialSet<Legendre<Dim, Order, Normalized<false>, T>, PolySetType >
{
    typedef PolynomialSet<Legendre<Dim, Order, Normalized<false>, T>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const bool isTransformationEquivalent = true;

    typedef OrthogonalPolynomialSet<Dim, Order, PolySetType, T, Hypercube> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef Legendre<Dim, Order, Normalized<false>, T> basis_type;
    typedef Hypercube<Dim, Order, Dim> convex_type;
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

    OrthogonalPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
        this->setCoefficient( polyset_type::toType( m ), true );
    }

    OrthogonalPolynomialSet<Dim, Order, Scalar,T, Hypercube > toScalar() const
    {
        return OrthogonalPolynomialSet<Dim, Order, Scalar,T, Hypercube >();
    }

    std::string familyName() const
    {
        return "legendre";
    }
};

template<uint16_type Dim,
         uint16_type Order,
         template<uint16_type> class PolySetType,
         typename T>
const uint16_type OrthogonalPolynomialSet<Dim, Order,PolySetType,T, Hypercube>::nLocalDof;
}
#endif /* __OrthogonalPolynomialSet_H */
