/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 26 Mar 2015
 
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
#include <vector>
#include <iostream>


#include <boost/mpi.hpp>

int main( int argc, char**argv )
{
    namespace mpi=boost::mpi;
    mpi::environment env( argc, argv );
    mpi::communicator world;
    std::vector<int> x ( 2 );
    x[0] = world.rank()%2;
    x[1] = world.rank()%2+1;
    std::vector<int> v( 2*world.size() );
    mpi::all_gather( world, x.data(), 2, v );
    std::cout << "size(): " << v.size() << std::endl;
    std::for_each( v.begin(), v.end(), []( int i ) {  std::cout << i << ","; } );
    std::vector<int> p1,p2;
    for( int i = 0;i < world.size(); ++i )
    {
        if ( v[2*i] ) p1.push_back( i );
        if ( v[2*i+1] ) p2.push_back( i );
    }
    auto g1 = world.group().include( p1.begin(), p1.end() );
    mpi::communicator g1comm ( world, g1 );
    p1.push_back( p2.front() );
    auto g2 = world.group().include( p1.begin(), p1.end() );    
    mpi::communicator g2comm ( world, g2 );
    
    std::string value;
    if ( g2comm )
    {
    
        if (g2comm.rank() == g2comm.size()-1) {
            value = "Hello, World!";
        }
        
        mpi::broadcast(g2comm, value, g2comm.size()-1);
        
        std::cout << "Process #" << g2comm.rank() << " global #" << world.rank() << " says " << value
                  << std::endl;
    }
    
}
    
