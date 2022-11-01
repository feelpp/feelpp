/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2022 Feel++ Consortium
  Copyright (C) 2011-2022 Université de Strasbourg

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
#ifndef FEELPP_CONVEX_HPP
#define FEELPP_CONVEX_HPP 1

namespace Feel
{
class ConvexBase {};
/**
 * @brief Convex base class
 */
template<uint16_type Dim, int Order, uint16_type RDim = Dim>
class Convex : public ConvexBase
{
public:


    /** @name Constants
     */
    //@{
    inline static const uint16_type nDim = Dim;
    inline static const int nOrder = Order;
    inline static const uint16_type nRealDim = RDim;


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    Convex() {}

    //! copy constructor
    Convex( Convex const & ) {}

    //! destructor
    virtual ~Convex() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Convex& operator=( Convex const & o ) = default;

    //@}

    /** @name Accessors
     */
    //@{


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

};
} // Feel
#endif /* __Convex_H */
