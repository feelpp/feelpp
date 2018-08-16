/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-12-30

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file feelmodelslog.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-12-30
 */


#ifndef FEELMODELS_LOG_HPP
#define FEELMODELS_LOG_HPP 1

#include <iostream>
#include <string>

#include <feel/feelcore/worldcomm.hpp>

namespace Feel {
namespace FeelModels {

    void Log(std::string const& _className,std::string const& _functionName, std::string const& _msg);
    void Log(std::string const& _className,std::string const& _functionName, std::string const& _msg, WorldComm const& worldComm, bool allproc=false);
    void Log(std::string const& _msgbefore,std::string const& _className, std::string const& _functionName, std::string const& _msg, WorldComm const& worldComm, bool allproc=false);


} // namespace FeelModels
} // namespace Feel


#endif // FEELMODELS_LOG_HPP
