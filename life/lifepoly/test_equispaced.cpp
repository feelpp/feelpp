/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-09-13

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file test_equispaced.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-09-13
 */
#include <life/lifecore/life.hpp>
#include <life/lifepoly/equispaced.hpp>

int main()
{
    using namespace Life;

    PointSetEquiSpaced<Simplex<2, 3, 2> , 3, double> pset3;
    pset3.toPython();

    PointSetEquiSpaced<Simplex<2, 4, 2> , 4, double> pset4;
    pset4.toPython();

    PointSetEquiSpaced<Simplex<2, 5, 2> , 5, double> pset5;
    pset5.toPython();
}
