/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>

  Date 2013-02-19

  Copyright (C) 2013 Université de Strasbourg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <feel/feel.hpp>
using namespace Feel;

int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="stokesCerebroVeinous",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>>,
                          _filename="brainVeines20.msh",
                          _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                          _physical_are_elementary_regions=true,
                          _partitions=Environment::worldComm().localSize(),
                          _worldcomm=Environment::worldComm()
                         );



    // function space
    auto Vh = THch<2>( mesh );

    // element U=(u,p) in Vh
    auto U = Vh->element();
    auto u = U.element<0>();
    auto p = U.element<1>();

    auto hsize = mean(elements(mesh),h());
    std::cout << "   hsize = " << hsize << "\n";
    //+++++++++++++++++++++++++ Velocity at inlet +++++++++++++++++++++++++++

    auto sumAreaIN = integrate(markedfaces(mesh,52),cst(1.)).evaluate()(0,0); //52 is a physical group assembling all the entries
    std::cout << "The total sum of the inlet area = "<<sumAreaIN << "\n";
    auto u_inlet = -6000./sumAreaIN;
    std::cout << "u_inlet = "<< u_inlet << "\n";
    auto flowIn = integrate(markedfaces(mesh,52),inner(u_inlet*N(),N() ) ).evaluate()(0,0);
    std::cout<<"flow_inlet = "<<flowIn<<"\n";

    //+++++++++++++++++++++++++ Velocity at Outlet +++++++++++++++++++++++++++
    auto sumAreaOUT = integrate(markedfaces(mesh,51),cst(1.)).evaluate()(0,0);
    std::cout << "The total sum of the outlet area = "<<sumAreaOUT << "\n";
    auto u_outlet = 6000./sumAreaOUT;
    std::cout << "u_outlet = "<< u_outlet << "\n";
    auto flowOut = integrate(markedfaces(mesh,51),inner(u_outlet*N(),N() ) ).evaluate()(0,0) ;
    std::cout<<"flow_outlet = "<<flowOut<<"\n";
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // left hand side
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=trace(gradt(u)*trans(grad(u))) );

    a+= integrate(_range=elements(mesh),
                  _expr=-div(u)*idt(p)-divt(u)*id(p));

    // right hand side
    auto f = vec(cst(0.),cst(0.),cst(0.));
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh), _expr=inner(f,id(u)));
    std::cout << "Start imposing boundary conditions \n";


    //++++++++ Boundary Conditions ++++++++++++++++
    //++++ At the inlets ++++
    a+=on(_range=markedfaces(mesh,52), _rhs=l, _element=u, _expr=vec(u_inlet*Nx(),u_inlet*Ny(),u_inlet*Nz()));
    //++++ On the wall ++++
    a+=on(_range=markedfaces(mesh,53),  _rhs=l, _element=u, _expr=vec(cst(0.), cst(0.), cst(0.)) );
    std::cout << "Start solving the system \n";

    /*//++++ At the outlets ++++
    a+=on(_range=markedfaces(mesh,7),  _rhs=l, _element=u, _expr=vec(u_outlet*Nx(),u_outlet*Ny(),u_outlet*Nz()) );
    a+=on(_range=markedfaces(mesh,24), _rhs=l, _element=u, _expr=vec(u_outlet*Nx(),u_outlet*Ny(),u_outlet*Nz()) );*/

    //+++++++++ Solve +++++++++
    a.solve(_rhs=l,_solution=U);

    //Flow at inlet
    auto flowInCal = integrate(markedfaces(mesh,52),inner(idv(u),N() ) ).evaluate()(0,0);
    //Flow at the outlet
    auto flowOutCal = integrate(markedfaces(mesh,51),inner(idv(u),N() ) ).evaluate()(0,0) ;
    std::cout.precision(20);
    std::cout<<"flow_inlet calcule = "<<flowInCal<<"\n";
    std::cout<<"flow_outlet calcule = "<<flowOutCal<<"\n";
    std::cout<<"Diff flow In/Out = "<<flowInCal+flowOutCal<<"\n";
    std::cout<<"Diff flow In/In = "<<flowInCal-flowIn<<"\n";
    std::cout<<"Diff flow Out/Out = "<<flowOutCal-flowOut<<"\n";

    //auto hsize = mean(elements(mesh),h());

    LOG(INFO) << "   hsize = " << hsize << "\n";
    LOG(INFO) << "[dof]             number of dof: " << Vh->nDof() << "\n";
    LOG(INFO) << "[dof]        number of dof/proc: " << Vh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]          number of dof(U): " << Vh->functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof]     number of dof/proc(U): " << Vh->functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]          number of dof(P): " << Vh->functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof]     number of dof/proc(P): " << Vh->functionSpace<1>()->nLocalDof()  << "\n";

    std::cout << "Mesh size = "<< hsize << "\n";
    std::cout << "[dof]             number of dof: " << Vh->nDof() << "\n";
    std::cout << "[dof]        number of dof/proc: " << Vh->nLocalDof() << "\n";
    std::cout << "[dof]          number of dof(U): " << Vh->functionSpace<0>()->nDof()  << "\n";
    std::cout << "[dof]     number of dof/proc(U): " << Vh->functionSpace<0>()->nLocalDof()  << "\n";
    std::cout << "[dof]          number of dof(P): " << Vh->functionSpace<1>()->nDof()  << "\n";
    std::cout << "[dof]     number of dof/proc(P): " << Vh->functionSpace<1>()->nLocalDof()  << "\n";


    // save results
    auto e = exporter( _mesh=mesh );
   e->add( "u", u );
   e->add( "p", p );
   e->save();
}
/// [marker_main]
