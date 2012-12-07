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
   \file constrainedpolynomialset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-30
 */
#ifndef __ConstrainedPolynomialSet_H
#define __ConstrainedPolynomialSet_H 1


namespace Feel
{

template<typename P> class Functional;
template<typename P> class FunctionalSet;

/**
 * \class ConstrainedPolynomialSet
 */
template<typename Poly>
class ConstrainedPolynomialSet
    :
public mpl::if_<mpl::bool_<Poly::is_scalar>,
    mpl::identity<PolynomialSet<Poly> >,
    mpl::identity<PolynomialSet<Poly, Vectorial> > >::type::type
{
    typedef typename mpl::if_<mpl::bool_<Poly::is_scalar>,
            mpl::identity<PolynomialSet<Poly> >,
            mpl::identity<PolynomialSet<Poly, Vectorial> > >::type::type super;
public:
    /** @name Constants
     */
    //@{

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const bool is_vectorial = super::is_vectorial;
    static const bool is_scalar = super::is_scalar;
    static const uint16_type nComponents = super::nComponents;
    //@}

    /** @name Typedefs
     */
    //@{
    typedef ConstrainedPolynomialSet<Poly> self_type;
    typedef Poly space_type;
    typedef typename super::value_type value_type;
    typedef typename super::basis_type basis_type;

    typedef typename super::polyset_type polyset_type;
    typedef Functional<space_type> constraint_type;
    typedef FunctionalSet<space_type> constraintset_type;


    typedef PolynomialSet<space_type, Scalar> component_type;
    typedef typename mpl::if_<mpl::bool_<is_scalar>,
            mpl::identity<Polynomial<space_type> >,
            mpl::identity<Polynomial<space_type, Vectorial> > >::type::type polynomial_type;

    typedef typename super::convex_type convex_type;
    typedef typename basis_type::matrix_type matrix_type;
    typedef typename basis_type::points_type points_type;


    BOOST_STATIC_ASSERT( ( boost::is_same<typename matrix_type::value_type, value_type>::value ) );
    BOOST_STATIC_ASSERT( ( boost::is_same<typename matrix_type::value_type, typename points_type::value_type>::value ) );

    //@}

    /** @name Constructors, destructor
     */
    //@{

    ConstrainedPolynomialSet()
        :
        super()
    {}

    void setConstraints( constraintset_type const& fset )
    {
        // form the matrix associated with the functional fset applied
        // to the functionSpace
        matrix_type m( fset( fset.functionSpace() ) );

        //std::cout << "[ConstrainedPolynomialSet] m = " << m << "\n";

        // apply svd to determine the intersection of the kernels of
        // the linear functionals
        if ( is_vectorial )
        {
            //m = vectorialToMatrix( m, nComponents );
            std::cout << "[ConstrainedPolynomialSet] v2m(m) = " << m << "\n";
        }

        SVD<matrix_type> svd( m );

        //extract the coefficients of V associated with the null
        //singular values
        matrix_type mv ( ublas::subrange( svd.V(), svd.S().size(), svd.V().size1(), 0, m.size2() ) );
        //std::cout << "[ConstrainedPolynomialSet] mv = " << mv << "\n";
        this->setCoefficient( polyset_type::toType( mv ), true  );

    }

    //@}
};
} // Feel
#endif /* __ConstrainedPolynomialSet_H */
