/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
      Date: 05/10/2015

 Copyright (C) 2015 Feel++ Consortium

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
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelts/timeset.hpp>


namespace Feel
{
void
TimeSet::print()
{
    if ( Environment::isMasterRank() )
    {
        std::cout <<    "------------------------------------------------------------\n";
        std::cout << "Time " << t() << "s, step=" << this->index() << ", k_{n+1}=" << k() << " k_n=" << kprev(1) << "\n";
    }
}
void
TimeSet::save( std::string const& fname )
{
    std::ofstream ofs(fname);
    for(auto i : range(this->size()))
    {
        ofs << std::setw( 5 ) << std::right << i << " "
            << std::setw( 11 ) << std::scientific << std::setprecision( 2 ) << std::right << this->at(i).timeStep() << " "
            << std::setw( 11 ) << std::scientific << std::setprecision( 2 ) << std::right << this->at(i).time()
            << std::setw( 11 ) << std::scientific << std::setprecision( 2 ) << std::right << this->at(i).error();
    }
}
} /* Feel */
