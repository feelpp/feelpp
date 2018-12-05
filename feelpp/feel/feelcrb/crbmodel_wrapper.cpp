/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-25

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file crbmodel_wrapper.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-25
 */
#include <boost/python.hpp>
using namespace boost::python;

#include <crbmodel.hpp>
#include <heat1d.hpp>

BOOST_PYTHON_MODULE( heat1d )
{
    using namespace Feel;
    class_<CRBModel<Heat1D> >( "Heat1D" )
    .def( init<>() )
    .def( "thetaAq", &CRBModel<Heat1D>::thetaAq )
    ;


}

