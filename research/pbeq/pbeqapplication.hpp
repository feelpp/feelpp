/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-08-31

  Copyright (C) 2007 Unil

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
   \file pbeqapplication.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-08-31
 */
#ifndef __PBEQAPPLICATION_H
#define __PBEQAPPLICATION_H 1

#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>

namespace Feel
{

typedef double value_type;
typedef Application  application_type;
typedef Backend<value_type> backend_type;
typedef boost::shared_ptr<backend_type> backend_ptrtype;

} // end namespace Feel


#endif // __PBEQAPPLICATION_H
