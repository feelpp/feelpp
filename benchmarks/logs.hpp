/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the FeelV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-06-22

  Copyright (C) 2005,2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file logs.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-06-22
 */
#ifndef CYLINDER_LOGS_HPP
#define CYLINDER_LOGS_HPP 1

#include <boost/log/log.hpp>

BOOST_DECLARE_LOG( app )
BOOST_DECLARE_LOG( dbg )
BOOST_DECLARE_LOG( err )
BOOST_DECLARE_LOG( warn )
BOOST_DECLARE_LOG( info )

#endif /* CYLINDER_LOGS_HPP */
