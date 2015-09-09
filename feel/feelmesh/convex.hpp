/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2015 Feel++ Consortium

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
   \file convex.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-30
 */
#ifndef FEELPP_MESH_CONVEX_HPP
#define FEELPP_MESH_CONVEX_HPP 1

namespace Feel
{
/**
 * \class Convex
 * \brief Convex base class
 *
 * @author Christophe Prud'homme
 * @see
 */
class Convex
{
public:


    /** @name Constants
     */
    //@{

    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    constexpr Convex( uint16_type Dim, uint16_type Order, uint16_type RDim )
        :
        M_dim( Dim ),
        M_order( Order ),
        M_rdim( RDim )
        {}

    constexpr Convex( uint16_type Dim, uint16_type RDim )
        :
        M_dim( Dim ),
        M_order( 1 ),
        M_rdim( RDim )
        {}

    constexpr Convex( Convex const & ) = default;
    Convex( Convex && ) = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Convex& operator=( Convex const & o ) = default;
    Convex& operator=( Convex && o ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    constexpr uint16_type dimension() const { return M_rdim; }
    constexpr uint16_type topologicalDimension() const { return M_dim; }
    constexpr uint16_type geometricalDimension() const { return M_rdim; }
    constexpr uint16_type realDimension() const { return M_rdim; }
    constexpr uint16_type order() const { return M_order; }


    //@}

    /** @name  Mutators
     */
    //@{

    //! set the order to \p o, allow for constexpr
    constexpr void setOrder( uint16_type o ) { M_order = o; }

    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:
    uint16_type M_dim;
    uint16_type M_order;
    uint16_type M_rdim;

private:

};

// Alias to Convex for compatibility
using ConvexBase = Convex;

} // Feel
#endif /* FEELPP_MESH_CONVEX_HPP */
