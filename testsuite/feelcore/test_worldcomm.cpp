/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Guillaume Dollé <gdolle@unistra.fr>
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

//#include <feel/feelcore/environment.hpp>
//#include <feel/feelcore/worldcomm.hpp>
//#include <feel/feelfilters/loadmesh.hpp>
//#include <feel/feelvf/form.hpp>
//#include <feel/feeldiscr/pch.hpp>
//#include <feel/feelvf/expr.hpp>
#include <feel/feel.hpp>

Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_worldcomm" ,
                           "test_worldcomm" ,
                           "0.1",
                           "MPI communicators test",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    return about;
}

// Create Npcomm communicators and solve in parallel two laplacians with
// two different boundary strong conditions.
int main( int argc, char**argv )
{
    using namespace Feel;
    Environment env(_argc=argc, _argv=argv,
                    _about=makeAbout() );

    auto world = Environment::worldComm();
    // Number of communicators to create (Np=k*Ncomm).
    const int Npcomm = 2;
    const int Np = world.size();


    if( (Np>=Npcomm) && Npcomm > 0 )
    {
        auto color = world.rank() % Npcomm;

        WorldComm w( color );

        //w.showMe();

        auto mesh = loadMesh( _mesh = new Mesh< Simplex<2> >,
                              _filename = "test_twodomains.geo",
                              _worldcomm = w.localComm(),
                              _partitions = w.localSize()
                              );

        // Create a feel++ worldcomm using local communicators.
        auto bend = backend( _worldcomm=w.localComm() );

        auto Xh = Pch<1>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();

        // Bilinear form2 create the matrix using the backend worldcomm.
        auto a = form2( _test=Xh, _trial=Xh, _backend=bend );

        a+= integrate( _range=markedelements(mesh,"omega"),
                       _expr=gradt(u)*trans(grad(v) ) );

        // Linear form1 create the vector using the backend worldcomm.
        auto l = form1( _test=Xh, _backend=bend );
        l+= integrate( _range=markedelements(mesh,"omega"),
                       _expr=id(v) );

        //
        if( color == 0 )
        {
            a += on( _range=markedfaces( mesh, "wall"), _rhs=l, _element=u, _expr=cst(0.) );
            //a += on( _range=markedfaces( mesh, "wall_left"), _rhs=l, _element=u, _expr=cst(1.) );
            //a += on( _range=markedpoints( mesh, "left"), _rhs=l, _element=u, _expr=cst(1.) );
        }
        else
        {
            a += on( _range=markedfaces( mesh, "wall"), _rhs=l, _element=u, _expr=cst(1.) );
            //a += on( _range=markedfaces( mesh, "wall_right"), _rhs=l, _element=u, _expr=cst(1.) );
            //a += on( _range=markedpoints( mesh, "right"), _rhs=l, _element=u, _expr=cst(1.) );
        }

        // solveb use the backend worldcomm
        a.solveb( _rhs=l, _solution=u, _backend=bend );


        std::cout << "global rank: " << w.localRank()
                  << " local rank:" << w.globalRank()
                  << " u size: " << u.size()
                  << " u[0]: " << u[0] << std::endl;

        // Exporter use the mesh worldcomm
        auto e = exporter( _mesh=mesh, _name=( boost::format("problem%1%") % color ).str() );
        e->add("u",u);
        e->save();
    }
}

