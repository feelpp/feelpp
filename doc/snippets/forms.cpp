/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 07 Dec 2014
 
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
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

struct OpT
{
    typedef Mesh<Simplex<2,1>> mesh_type;
    typedef typename Mesh<Simplex<2,1>>::ptrtype mesh_ptrtype;
    typedef typename Feel::meta::Pch<mesh_type,1>::type space_type;
    typedef typename space_type::element_type element_type;
    typedef typename Feel::meta::Pch<mesh_type,1>::ptrtype space_ptrtype;
    typedef typename Feel::meta::BilinearForm<space_type,space_type>::type form2_type;
    typedef typename Feel::meta::LinearForm<space_type>::type form1_type;
    typedef typename Feel::meta::Exporter<mesh_type>::ptrtype exporter_ptrtype;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    element_type u;
    form2_type a;
    form1_type l;
    exporter_ptrtype e;
    
    OpT() 
        {
            mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
            Xh = Pch<1>( mesh );
            u = Xh->element("u");
            auto v = Xh->element("g");
        
            l = form1( _test=Xh );
            l = integrate(_range=elements(mesh),
                          _expr=id(v));
        
            a = form2( _trial=Xh, _test=Xh);
            a = integrate(_range=elements(mesh),
                          _expr=gradt(u)*trans(grad(v)) );
            a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.) );
            e = exporter( _mesh=mesh );
        }
    double solve(double t = 1.)
        {
            a.solve(_rhs=l,_solution=u);
            e->step(t)->add( "u", u );
            e->save();
            // returns the mean of u : ()\int_\Omega u)/(\int_\Omega 1)
            return mean( _expr=idv(u), _range=elements(mesh))(0,0);
        }
};


int main(int argc, char** argv)
{

	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="forms",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    OpT op;
    
    double m = op.solve( 1.0 );
    if ( Environment::isMasterRank() )
        std::cout << "mean value = " << m << "\n";

    m = op.solve( 2.0 );
    if ( Environment::isMasterRank() )
        std::cout << "mean value = " << m << "\n";

    
    
    
    return 0;
}
