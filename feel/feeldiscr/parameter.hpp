/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-20

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file parameter.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-20
 */
#ifndef FEELPP_feeldiscr_parameter_H
#define FEELPP_feeldiscr_parameter_H 1

#include <boost/parameter.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/aux_/maybe.hpp>

#include <feel/feelpoly/parameter.hpp>

namespace Feel
{

namespace parameter = boost::parameter;

BOOST_PARAMETER_TEMPLATE_KEYWORD( model_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( mesh_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( expansion_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( bases_list )
BOOST_PARAMETER_TEMPLATE_KEYWORD( basis_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( domain_value_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( periodicity_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( continuity_type )
BOOST_PARAMETER_TEMPLATE_KEYWORD( mortar_type )

}
#endif /* __parameter_H */
