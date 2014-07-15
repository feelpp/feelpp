/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-07-15

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
   \file dofneighbors.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-07-15
 */
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description dofneighborsoptions( "dofneighbors options" );
	dofneighborsoptions.add_options()
    ( "dof", po::value<int>()->default_value( 0 ), "global dof id" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc= dofneighborsoptions,
                     _about=about(_name="dofneighbors",
                                  _desc="print neighbor dof information of a dof",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<2>( mesh );

    // loop over the multiset of global dof the dof table provides a relation
    // between local dof (element id + local dof id) and global dof. While the
    // global dof is unique, it is possibly associated to multiple local dof.
    // we can work here with the global view (global dof) of the relation or the
    // local view (local dof) of the relation.
    for( auto const& dof : Vh->dof()->globalDof( ioption( "dof" ) ) )
    {
        // first print the local dof associated to the global dof 'dof'
        std::cout << "dof element id : " << dof.second.elementId() << " local dof : " << dof.second.localDof() << "\n";
        // now print all neighbor dof including itself
        for( auto const& neighbordof : Vh->dof()->localDof( dof.second.elementId() ) )
        {
            std::cout << "  |- local dof neighbor id  : " << neighbordof.first.localDof()
                      << " -> global dof neighbor id " << neighbordof.second.index() << "\n";
        }
    }
    return 0;

}
