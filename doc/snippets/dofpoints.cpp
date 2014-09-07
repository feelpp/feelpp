/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-03-15

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
   \file dofpoints.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-03-15
 */
//# marker1 #
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<2> > mesh_type;
    typedef FunctionSpace<mesh_type,bases<Lagrange<2,Vectorial> > > fs_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="hypercube-2",
                                        _usenames=true,
                                        _shape="hypercube",
                                        _dim=2,
                                        _h=2 ) );
    auto Xh = fs_type::New( mesh );
    auto B = Xh->element();

    std::cout << "number of degees of freedom: " << Xh->nDof() << "\n";
    std::cout << "B.size: : " << B.size() << "\n";
    auto dofpt_it = Xh->dof()->dofPointBegin();
    auto dofpt_en = Xh->dof()->dofPointEnd();

    for ( ; dofpt_it != dofpt_en; ++dofpt_it )
    {
        // the data structure associated with each dot point is a tuple containing:
        //   - the dof point coordinate
        //   - the dof index
        //   - the dof component associated with the dofpoint
        auto dofpt_coord = dofpt_it->get<0>();
        auto dofpt_id = dofpt_it->get<1>();
        auto dofpt_comp = dofpt_it->get<2>();

        // do something with the coordinate and store it in the proper vector entry in B
        auto r = vec( ( Px()-dofpt_coord[0] )*( Px()-dofpt_coord[0] ),
                      ( Py()-dofpt_coord[1] )*( Py()-dofpt_coord[1] ) );
        auto I = integrate( elements( mesh ), r ).evaluate();
        std::cout << "I = " << I << "\n";
        B[dofpt_id] = I( dofpt_comp, 0 );
        std::cout << "Dof coordinate[" << dofpt_id << "]=" << dofpt_coord << ", component=" << dofpt_comp << "\n";
        std::cout << "B[" << dofpt_id << "]=" << B[dofpt_id] << "\n";

    }
}
//# endmarker1 #
