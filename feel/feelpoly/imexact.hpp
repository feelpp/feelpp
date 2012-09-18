/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-07-02

  Copyright (C) 2006 EPFL

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
   \file imexact.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-07-02
 */
#ifndef __IMExact_H
#define __IMExact_H 1

namespace Feel
{
/**
 * \class IMExact
 * \brief exact integration method
 *
 *  \ingroup Polynomial
 *  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<typename T = double>
class IMExact
{
public:


    /** @name Typedefs
     */
    //@{

    static const bool is_exact = true;
    static const bool is_face_im = false;

    typedef T value_type;
    typedef ublas::matrix<value_type,ublas::column_major> points_type;

    typedef IMExact<value_type> face_quadrature_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    IMExact() {}
    IMExact( IMExact const & ) {}
    ~IMExact() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * no points with exact integration
     */
    uint16_type nPoints() const
    {
        return 0;
    }

    /**
     * dummy points
     */
    points_type points() const
    {
        return points_type();
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * the integration has already been calculated elsewhere
     * @see feel/feelf/operators2.hpp
     */
    template<typename Expression>
    value_type integrate( Expression  const&  f ) const
    {
        uint32_type k = 0;
        return value_type( f( k ) );
    }

    //@}



protected:

private:

};
} // Feel
#endif /* __IMExact_H */
