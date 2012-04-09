/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
        _M_space(),
        _M_fset(),
        _M_mat()
    {}

    FunctionalSet( space_type const& s )
        :
        _M_space( s ),
        _M_fset(),
        _M_mat()
    {
    }
    FunctionalSet( space_type const& s, std::vector<functional_type> const& fset )
        :
        _M_space( s ),
        _M_fset( fset ),
        _M_mat( space_type::nComponents*fset.size(), fset[0].coeff().size2() )
    {
        //std::cout << "FunctionalSet: " << fset[0].coeff() <<  "\n";
        this->setFunctionalSet( fset );
    }
    FunctionalSet( FunctionalSet const & fset )
        :
        _M_space( fset._M_space ),
        _M_fset( fset._M_fset ),
        _M_mat( fset._M_mat )
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
            _M_space = fset._M_space;
            _M_fset = fset._M_fset;
            _M_mat = fset._M_mat;
        }

        return *this;
    }

    /**
     * \return the i-th functional
     */
    functional_type const& operator()( uint16_type i ) const
    {
        return _M_fset[i];
    }

    /**
     * \return the value of the functional set applied to a polynomial
     */
    matrix_type operator()( space_type const& p ) const
    {
        //FEELPP_ASSERT( _M_mat.size2() == ublas::trans(p.coeff()).size1() )( _M_mat.size1() )( p.coeff().size1() ).error( "incompatible dimension between functional and polynomial.\n Is the space correctly defined?" );

        return ublas::prod( space_type::polyset_type::toMatrix( _M_mat ),
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
        return _M_space;
    }

    /**
     * This works only if the basis is orthonormal
     *
     * \return the representation of the functional set using basis
     * of the function space.
     */
    matrix_type const& rep() const
    {
        return _M_mat;
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
        _M_space = __space;
    }

    /**
     * set the Functional set
     */
    void setFunctionalSet( std::vector<functional_type> const& fset )
    {
        _M_fset = fset;


        if ( space_type::is_scalar )
        {
            // update matrix associated with functionals applied to the
            // basis of the function space
            _M_mat = ublas::zero_matrix<value_type>( fset.size(), fset[0].coeff().size2() );

            //std::cout << "mat size" << _M_mat << "\n";
            for ( uint16_type i = 0; i < fset.size(); ++i )
            {
                //std::cout << "Functional " << i << "=" << fset[i].coeff() << "\n";
                ublas::row( _M_mat, i ) = ublas::row( fset[i].coeff(), 0 );
            }

            //std::cout << "mat size" << _M_mat << "\n";

        }

        else
        {
            // update matrix associated with functionals applied to the
            // basis of the function space
            _M_mat = ublas::zero_matrix<value_type>( space_type::nComponents*fset.size(), fset[0].coeff().size2() );

            for ( uint16_type i = 0; i < fset.size(); ++i )
            {
                ublas::project( _M_mat,
                                ublas::range( i*space_type::nComponents, ( i+1 )*space_type::nComponents ),
                                ublas::range( 0, _M_mat.size2() ) ) = ublas::scalar_matrix<value_type>( space_type::nComponents, _M_mat.size2(), -1 );
                ublas::project( _M_mat,
                                ublas::range( i*space_type::nComponents, ( i+1 )*space_type::nComponents ),
                                ublas::range( 0, _M_mat.size2() ) ) = fset[i].coeff();
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
    space_type _M_space;
    fset_type _M_fset;
    matrix_type _M_mat;
};
} // Feel
#endif /* __FunctionalSet_H */
