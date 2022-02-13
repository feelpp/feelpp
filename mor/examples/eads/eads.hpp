/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-03

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file eads.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-03
 */
#ifndef __Eads_H
#define __Eads_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>


namespace Feel
{
/**
 * Defines the \c AboutData for a Eads testcase application
 * \param appname string containing the application for the Eads testcase
 */
AboutData makeEadsAbout( std::string const& appname );

/**
 * defines the options for all Eads applications
 */
po::options_description makeEadsOptions(  );
}
#endif /* __Eads_H */
