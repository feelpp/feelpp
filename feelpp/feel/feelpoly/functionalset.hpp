/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Core>
#include <feel/feelcore/warnon.hpp>
// clang-format on


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
    using eigen_matrix_type = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using eigen_matrix_map_type = Eigen::Map<eigen_matrix_type>;
    using eigen_const_matrix_map_type = Eigen::Map<eigen_matrix_type const>;

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
        M_mat()
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
        if ( this != &fset )
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

    /** \return basis-dependent matrix representation of the ordered dual set. */
    matrix_type const& dualMatrix() const noexcept
    {
        return M_mat;
    }

    /** \return zero-copy Eigen view of dualMatrix(). */
    [[nodiscard]] eigen_const_matrix_map_type eigenDualMatrix() const noexcept
    {
        return eigen_const_matrix_map_type( M_mat.data().begin(),
                                            M_mat.size1(), M_mat.size2() );
    }

    /** \return number of mathematical functionals in the ordered set. */
    [[nodiscard]] std::size_t size() const noexcept
    {
        return M_fset.size();
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
        if ( fset.empty() )
        {
            M_mat.resize( 0, 0, false );
            return;
        }

        const std::size_t rowsPerFunctional = space_type::is_scalar ? 1 : space_type::nComponents;
        const std::size_t columns = fset.front().coeff().size2();
        M_mat.resize( rowsPerFunctional*fset.size(), columns, false );
        eigen_matrix_map_type dualMap( M_mat.data().begin(), M_mat.size1(), M_mat.size2() );
        dualMap.setZero();

        for ( std::size_t i = 0; i < fset.size(); ++i )
        {
            auto const functionalMap = fset[i].rieszRepresentation();
            CHECK_EQ( static_cast<std::size_t>( functionalMap.rows() ), rowsPerFunctional )
                << "invalid functional value shape in ordered dual set";
            CHECK_EQ( static_cast<std::size_t>( functionalMap.cols() ), columns )
                << "inconsistent functional representation width in ordered dual set";
            dualMap.middleRows( static_cast<Eigen::Index>( i*rowsPerFunctional ),
                                static_cast<Eigen::Index>( rowsPerFunctional ) ) = functionalMap;
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
