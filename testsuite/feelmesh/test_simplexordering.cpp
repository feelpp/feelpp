/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 12 août 2015
 
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
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/refentity.hpp>

int main()
{
    using namespace Feel;
    Simplex<2> T2(3);
    std::cout << "T2: " << T2 << std::endl;
    auto Tref21 = reference( T2 );
    std::cout << "Tref21: " << Tref21 << std::endl;
    T21.setOrder( 2 );

    Simplex<3> T3(3);
    auto Tref31 = reference( T3 );
    std::cout << "T3: " << T3 << std::endl;
    std::cout << "Tref31: " << Tref31 << std::endl;
    
    
    
}

