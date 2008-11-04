/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-08-31

  Copyright (C) 2007 Unil

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
   \file pbeqapplication.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-08-31
 */
#ifndef __PBEQAPPLICATION_H
#define __PBEQAPPLICATION_H 1

#include <life/lifecore/application.hpp>
#include <life/lifealg/backend.hpp>

namespace Life
{

typedef double value_type;
typedef Application  application_type;
typedef Backend<value_type> backend_type;
typedef boost::shared_ptr<backend_type> backend_ptrtype;

} // end namespace Life


#endif // __PBEQAPPLICATION_H
