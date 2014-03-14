/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-30

  Copyright (C) 2014 Feel++ Consortium

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
   \file worldscomm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
#ifndef FEELPP_WORLDSCOMM_HPP
#define FEELPP_WORLDSCOMM_HPP 1

#include <vector>

#include <feel/feelcore/worldcomm.hpp>

namespace Feel {

/**
 * a set of worlds communicator
 *
 * @author Christophe Prud'homme
 * @see
 */
class WorldsComm : public std::vector<WorldComm>
{
public:
    typedef std::vector<WorldComm> super;

    WorldsComm( std::vector<WorldComm> const& s ) : super(s) {}
    WorldsComm( std::vector<WorldComm>&& s ) : super(s) {}
};

inline WorldsComm
worldsComm( WorldComm const& wc )
{
    return WorldsComm( std::vector<WorldComm>( 1, wc ) );
}


}
#endif /* FEELPP_WORLDSCOMM_HPP */
