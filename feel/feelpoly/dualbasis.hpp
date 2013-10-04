/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-10

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
   \file dualbasis.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-10
 */
#ifndef __DualBasis_H
#define __DualBasis_H 1

namespace Feel
{
/**
 * \class DualBasis
 * \brief basis of a space P' dual of some space P
 *
 * The element of this basis are functionals evaluated at a node set
 * It contains also the information about the dof table
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 */
template<typename Primal>
class DualBasis
{
public:


    /** @name Typedefs
     */
    //@{
    typedef DualBasis<Primal> self_type;

    typedef Primal primal_space_type;
    typedef typename primal_space_type::basis_type basis_type;

    //@}
    /** @name Constants
    */
    //@{
    static const uint16_type nDim = primal_space_type::nDim;
    static const uint16_type nOrder = primal_space_type::nOrder;
    //@}
    /** @name Constructors, destructor
     */
    //@{


    DualBasis( primal_space_type const& primal )
        :
        M_primal( primal )
    {}
    DualBasis( DualBasis const & b )
        :
        M_primal( b.M_primal )
    {}
    ~DualBasis()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& dual )
    {
        if ( this != &dual )
        {
            M_primal = dual.M_primal;
        }

        return *this;
    }

    basis_type const& operator()() const
    {
        return M_primal.basis();
    }

    //@}

    /** @name Accessors
     */
    //@{

    primal_space_type const& primalSpace() const
    {
        return M_primal;
    }

    basis_type const& basis() const
    {
        return M_primal.basis();
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{



    //@}



protected:

private:

    primal_space_type M_primal;
};
} // Feel

#endif /* __DualBasis_H */
