/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-05-03

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
   \file debug.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-05-03
 */
#include <Eigen/Core>
#include <feel/feelcore/debug.hpp>

namespace Feel
{
template<typename T, int N1=Eigen::Dynamic, int N2=Eigen::Dynamic>
DebugStream&
operator<<( DebugStream& __os, Eigen::Matrix<T,N1,N2> const& __n )
{

    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}

template<typename T, int N1=Eigen::Dynamic, int N2=Eigen::Dynamic>
NdebugStream&
operator<<( NdebugStream& __os, Eigen::Matrix<T,N1,N2> const&  )
{
    return __os;
}


} //  Feel
