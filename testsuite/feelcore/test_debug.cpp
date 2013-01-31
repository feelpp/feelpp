/*
  This file is part of the Feel library.

  Author: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

  Copyright (C) 2004 EPFL
  Copyright (C) 2013 Feel++ Consortium

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
#include <iostream>
#include <feel/feelcore/feel.hpp>

int
main()
{
    VLOG(1) << "Glog verbosity level 1 is now enabled\n";

    VLOG(2) << "Glog verbosity level 2 is now enabled\n";
}

