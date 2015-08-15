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
#include <feel/feelpoly/equispaced.hpp>

template<typename C>
void
test_convex( C && c )
{
    using namespace Feel;
    std::cout << c.name() << ":" << c << std::endl;
    auto Tref21 = reference( c, 1 );
    std::cout << "reference(C,1): " << Tref21 << std::endl;
    Tref21.setOrder( 2 );
    std::cout << "reference(C,2): " << Tref21 << std::endl;
}
template<typename C>
void
test_equispaced( C && c, int o )
{
    using namespace Feel;
    PointSetEquiSpaced<C,double> pset( o );
    std::cout << pset << std::endl;
}
int main()
{
    using namespace Feel; 
    test_convex( Simplex<1>() );
    test_convex( Simplex<1,3>() );
    test_convex( Simplex<2>() );
    test_convex( Simplex<2,3>() );
    test_convex( Simplex<3>() );

    test_convex( Hypercube<1>() );
    test_convex( Hypercube<2>() );
    test_convex( Hypercube<3>() );

    test_equispaced( Simplex<2>(), 1 );
    test_equispaced( Simplex<3>(), 1 );
    
}

