/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file unaryfunctor.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-18
 */

#include <feel/feelvf/unaryfunctor.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
UnaryFunctor::UnaryFunctor( std::string const& name,
                            boost::shared_ptr<FunctorDomain> const& domain )
    : M_name( name ),
      M_domain( domain )
{
}

/// \endcond
}
}
