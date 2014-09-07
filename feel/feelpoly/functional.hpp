/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-06

  Copyright (C) 2005-2006 EPFL

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
   \file functional.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-06
 */
#ifndef __Functional_H
#define __Functional_H 1

#include <boost/operators.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <feel/feelcore/feel.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

/**
 * \class Functional
 * \brief represents a linear functional
 *
 * A functional is defined by a polynomial set and a set of
 * coefficients
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
*/
template<typename Space>
class Functional
    :
public boost::addable<Functional<Space> >
{
    typedef boost::addable<Functional<Space> > super;
public:


    /** @name Typedefs
     */
    //@{

    typedef typename Space::value_type value_type;
    typedef Functional<Space> self_type;
    typedef Space space_type;
    typedef Space polynomialset_type;
    typedef typename space_type::polynomial_type polynomial_type;
    typedef typename space_type::basis_type basis_type;
    typedef typename space_type::matrix_type matrix_type;

    static const uint16_type nComponents = space_type::nComponents;

    // representation type for the functionals
    typedef ublas::matrix<value_type> rep_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Functional()
        :
        super(),
        M_p(),
        M_coeff()
    {}

    Functional( space_type const& P )
        :
        super(),
        M_p( P ),
        M_coeff( M_p.coeff() )

    {

    }

    //template<class AE>
    Functional( space_type const& P,
                matrix_type const& coeff )
    //ublas::matrix_expression<AE> const& coeff )
        :
        super(),
        M_p( P ),
        M_coeff( coeff )
    {
        //FEELPP_ASSERT( M_coeff.size1() == nComponents && M_coeff.size2() == M_p.polynomialDimensionPerComponent() )
        //( M_coeff.size1() )( M_coeff.size2() )( M_p.polynomialDimension() ).error( "invalid coefficient size" );
        //         FEELPP_ASSERT( M_coeff.size2() == 1 )( M_coeff.size2() ).error( "there should be only one column" );
    }

    Functional( Functional const & __f )
        :
        M_p( __f.M_p ),
        M_coeff( __f.M_coeff )
    {
        //FEELPP_ASSERT( M_coeff.size1() == nComponents && M_coeff.size2() == M_p.polynomialDimensionPerComponent() )
        //( M_coeff.size1() )( M_coeff.size2() )( M_p.polynomialDimensionPerComponent() ).error( "invalid coefficient size" );
    }

    virtual ~Functional()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& __f )
    {
        if ( this != &__f )
        {
            M_p = __f.M_p;
            M_coeff = __f.M_coeff;
        }

        return *this;
    }

    /**
     * add to another functional
     * it generates automatically operator+ thanks to addable
     */
    self_type& operator+=( const self_type& __f )
    {
        M_coeff += __f.M_coeff;
        return *this;
    }

    /**
     * apply the functional to a polynomial
     *
     *
     * \param p polynomial
     * \return matrix resulting from the application of the functional to the polynomial
     */
    virtual matrix_type operator()( polynomial_type const& p ) const
    {
        FEELPP_ASSERT( p.coeff().size2()  == M_coeff.size2() )
        ( p.coeff() )( M_coeff ).error( "invalid polynomial" );

        return ublas::prod( p.coeff(), ublas::trans( M_coeff ) );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the dimension of the polynomial space
     */
    uint16_type size() const
    {
        return M_coeff.size2();
    }

    /**
     * \return the coefficient of the functional in the basis
     * associated with the polynomial space
     */
    rep_type const& coeff() const
    {
        return M_coeff;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //template<class AE>
    //void setCoefficient( ublas::matrix_expression<AE> const& __coeff )
    void setCoefficient( matrix_type const& __coeff )
    {
        M_coeff = __coeff;
    }


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    space_type M_p;
    rep_type M_coeff;
};

}
#endif /* __Functional_H */
