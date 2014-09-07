/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-11-16

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
   \file fsfunctional.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-11-16
 */

#ifndef _FSFUNCTIONAL_HPP_
#define _FSFUNCTIONAL_HPP_

namespace Feel
{

// Functional on function space
template<class Space>
class FsFunctional
{
public:

    // -- TYPEDEFS --
    typedef Space space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef typename space_type::value_type value_type;

    FsFunctional( space_ptrtype space ) :
        M_space( space )
    {
    }

    virtual ~FsFunctional() {}

    // apply the functional
    virtual value_type
    operator()( const element_type& x ) const = 0;

    space_ptrtype space() const
    {
        return M_space;
    }

private:

    space_ptrtype M_space;

}; // class FsFunctional

} // Feel

#endif /* _FSFUNCTIONAL_HPP_ */
