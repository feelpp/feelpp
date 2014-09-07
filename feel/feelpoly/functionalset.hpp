/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-11

  Copyright (C) 2005,2006 EPFL

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
   \file functionalset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-11
 */
#ifndef __FunctionalSet_H
#define __FunctionalSet_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <feel/feelpoly/functional.hpp>

namespace Feel
{
/**
 * \class FunctionalSet
 * \brief Set of functionals
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<typename Space>
class FunctionalSet
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Space space_type;
    typedef typename space_type::value_type value_type;


    typedef FunctionalSet<Space> functionalset_type;
    typedef functionalset_type self_type;
    typedef Functional<Space> functional_type;


    typedef typename space_type::matrix_type matrix_type;

    typedef std::vector<functional_type> fset_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    FunctionalSet()
        :
        M_space(),
        M_fset(),
        M_mat()
    {}

    FunctionalSet( space_type const& s )
        :
        M_space( s ),
        M_fset(),
        M_mat()
    {
    }
    FunctionalSet( space_type const& s, std::vector<functional_type> const& fset )
        :
        M_space( s ),
        M_fset( fset ),
        M_mat( space_type::nComponents*fset.size(), fset[0].coeff().size2() )
    {
        //std::cout << "FunctionalSet: " << fset[0].coeff() <<  "\n";
        this->setFunctionalSet( fset );
    }
    FunctionalSet( FunctionalSet const & fset )
        :
        M_space( fset.M_space ),
        M_fset( fset.M_fset ),
        M_mat( fset.M_mat )
    {}

    ~FunctionalSet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& fset )
    {
        if ( this != fset )
        {
            M_space = fset.M_space;
            M_fset = fset.M_fset;
            M_mat = fset.M_mat;
        }

        return *this;
    }

    /**
     * \return the i-th functional
     */
    functional_type const& operator()( uint16_type i ) const
    {
        return M_fset[i];
    }

    /**
     * \return the value of the functional set applied to a polynomial
     */
    matrix_type operator()( space_type const& p ) const
    {
        //FEELPP_ASSERT( M_mat.size2() == ublas::trans(p.coeff()).size1() )( M_mat.size1() )( p.coeff().size1() ).error( "incompatible dimension between functional and polynomial.\n Is the space correctly defined?" );

        return ublas::prod( space_type::polyset_type::toMatrix( M_mat ),
                            ublas::trans( space_type::polyset_type::toMatrix( p.coeff() ) ) );
    }
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the function space from which the functionals take
     * their values
     */
    space_type const& functionSpace() const
    {
        return M_space;
    }

    /**
     * This works only if the basis is orthonormal
     *
     * \return the representation of the functional set using basis
     * of the function space.
     */
    matrix_type const& rep() const
    {
        return M_mat;
    }


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the function space
     */
    void setFunctionSpace( space_type const& __space )
    {
        M_space = __space;
    }

    /**
     * set the Functional set
     */
    void setFunctionalSet( std::vector<functional_type> const& fset )
    {
        M_fset = fset;


        if ( space_type::is_scalar )
        {
            // update matrix associated with functionals applied to the
            // basis of the function space
            M_mat = ublas::zero_matrix<value_type>( fset.size(), fset[0].coeff().size2() );

            //std::cout << "mat size" << M_mat << "\n";
            for ( uint16_type i = 0; i < fset.size(); ++i )
            {
                //std::cout << "Functional " << i << "=" << fset[i].coeff() << "\n";
                ublas::row( M_mat, i ) = ublas::row( fset[i].coeff(), 0 );
            }

            //std::cout << "mat size" << M_mat << "\n";

        }

        else
        {
            // update matrix associated with functionals applied to the
            // basis of the function space
            M_mat = ublas::zero_matrix<value_type>( space_type::nComponents*fset.size(), fset[0].coeff().size2() );

            for ( uint16_type i = 0; i < fset.size(); ++i )
            {
                ublas::project( M_mat,
                                ublas::range( i*space_type::nComponents, ( i+1 )*space_type::nComponents ),
                                ublas::range( 0, M_mat.size2() ) ) = ublas::scalar_matrix<value_type>( space_type::nComponents, M_mat.size2(), -1 );
                ublas::project( M_mat,
                                ublas::range( i*space_type::nComponents, ( i+1 )*space_type::nComponents ),
                                ublas::range( 0, M_mat.size2() ) ) = fset[i].coeff();
            }
        }
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:
    space_type M_space;
    fset_type M_fset;
    matrix_type M_mat;
};
} // Feel
#endif /* __FunctionalSet_H */
