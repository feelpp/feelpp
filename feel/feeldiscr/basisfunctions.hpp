/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014 Feel++ Consortium

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
   \file basisfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_BASISFUNCTIONS_HPP
#define FEELPP_BASISFUNCTIONS_HPP 1

namespace Feel
{
/**
 * \class BasisFunctions
 * \brief base class for basis functions
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename SpaceType, bool IsTest = true>
class BasisFunctions
{
public:


    /** @name Constants
     */
    //@{

    static const bool is_test = IsTest;
    static const bool is_trial = !IsTest;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef SpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    using space_type = functionspace_type;
    using space_ptrtype = functionspace_ptrtype;
    typedef typename SpaceType::basis_type basis_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    BasisFunctions(boost::shared_ptr<space_type> const& X )
        :
        M_space( X )
        {}
    //! copy constructor
    BasisFunctions( BasisFunctions const & bf )
        :
        M_space( bf.M_space )
        {}
    //! destructor
    virtual ~BasisFunctions() {}

    //@}


    /** @name Accessors
     */
    //@{
    space_ptrtype space() const { return M_space; }
    basis_type const& basis() const { return M_space->fe(); }
    //@}

protected:

    space_ptrtype M_space;
};
template<typename SpaceType> using Test = BasisFunctions<SpaceType,0>;

template<typename SpaceType>
Test<SpaceType> test( boost::shared_ptr<SpaceType> const& X )
{
    return Test<SpaceType>( X );
}

template<typename SpaceType> using Trial = BasisFunctions<SpaceType,1>;
template<typename SpaceType>
Trial<SpaceType> trial( boost::shared_ptr<SpaceType> const& X )
{
    return Trial<SpaceType>( X );
}

}
#Endif /* FEELPP_BASISFUNCTIONS_HPP */
