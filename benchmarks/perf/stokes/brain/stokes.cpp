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

/**
 \page StokesTutorial Stokes Tutorial
\author Feel++ Consortium
\date 2013-02-19

\tableofcontents
<br>
<hr>
<br>

\section Stokes_Theory Theory
Let solve the stokes equation considering a Poiseuille flow profile. <br>
We have the following system of equations,
<br><center>\f$
\left\{
\begin{aligned}
-\mu\Delta \bf u + \nabla p & =  \bf f & \text{on}\; \Omega \;, \\
                \nabla\cdot\bf u & =  0 & \text{on}\; \Omega \;, \\
                \bf u  & =  g & \text{on}\; \Gamma \;, \\
\end{aligned}
\right
\f$</center><br>
where \f$u\in [H_g^1(\Omega)]^d\f$ denotes the flow speed, \f$p\in [L_0^2(\Omega)]\f$ the fluid pressure, \f$\mu\f$ the fluid viscosity.<br>
The last boundary condition expresses a null pressure fixed on the outlet.<br>
The Poiseuille profile on the boundary is,
<br><center>\f$
g(x,y)=\left(
\begin{aligned}
 y(1-y) \\
 0      \\
\end{aligned}
\right)
\f$</center><br>
The method used to obtain the strong formulation is closed to the one used for the laplacian (see section \ref Laplacian ).<br>
We multiply the first equation by a test function \f$v\in H^1(\Omega)\f$ and we integrate on the domain \f$\Omega\f$,
<br><center>\f$
\begin{aligned}
 \left(
\int_\Omega \mu \nabla \mathbf u : \nabla \mathbf v
-\int_{\partial\Omega} \frac{\partial \mathbf u}{\partial \mathbf n} \cdot \mathbf v
\right)
+\int_\Omega ( \nabla\cdot(p \mathbf v) - \mathbf p \nabla\cdot v )
=\int_\Omega \mathbf f \cdot \mathbf v \;.
\end{aligned}
\right)
\f$</center><br>
where \f$n\f$ denotes a normal vector on the boundary.<br>
The divergence theorem (or Gauss's theorem) gives,
<br><center>\f$
\begin{aligned}
\int_\Omega \nabla\cdot(p \mathbf v) = \int_{\partial\Omega} p \mathbf v\cdot \mathbf n \;.
\end{aligned}
\right)
\f$</center><br>
We have to add a consistency terms to the equation to guaranty the symmetry of the bilinear form.<br>
This term is provided by the second equation. We multiply this equation by a test function \f$q\in L_2(\Omega)\f$ and we integrate on the domain \f$\Omega\f$,
<br><center>\f$
\begin{aligned}
\int_{\Omega} \nabla\cdot\mathbf u q = 0 \;,
\end{aligned}
\right)
\f$</center><br>
Finally, we deduce from the equations and after rearranging the integrals the variationnal formulation,
<br><center>\f$
\begin{aligned}
\int_\Omega \mu \nabla \mathbf u :\nabla \mathbf v
+\int_\Omega \left( \nabla\cdot\mathbf u q - p \nabla\cdot\mathbf v \right)
+
    \int_{\partial\Omega} \left( p \mathbf n -
   \frac{\partial \mathbf u}{\partial \mathbf n} \right)
     \cdot \mathbf v
=\int_\Omega \mathbf f \cdot \mathbf v
\end{aligned}
\right)
\f$</center><br>
Let us assume now that \f$(\mathbf v,q) \in [H_0^1(\Omega)]^d \times L_0^2(\Omega)\f$, the variationnal formulation leads to:
Find \f$(\mathbf u,p)\in [H_g^1(\Omega)]^d\times L_0^2(\Omega) \f$ such that for all \f$(\mathbf v,q) \in [H_0^1(\Omega)]^d \times L_0^2(\Omega)\f$
<br><center>\f$
\begin{aligned}
\int_\Omega \mu \nabla \mathbf u :\nabla \mathbf v
+\int_\Omega \left( \nabla\cdot\mathbf u q - p \nabla\cdot\mathbf v \right)
=\int_\Omega \mathbf f \cdot \mathbf v
\end{aligned}
\right)
\f$</center><br>
Or equivalently:
<br><center>\f$
\begin{aligned}
  a((\mathbf u,p),(\mathbf v,q)) = l((\mathbf v,q))
\end{aligned}
\right)
\f$</center><br>
where \f$a\f$ is a bilinear form, continuous, coercive and where \f$l\f$ is a linear form.

\section Stokes_Implementation Implementation
Let's see the \feel code corresponding to this mathematical statement (source \c "doc/manual/tutorial/mystokes.cpp").<br>
We suppose for this example the viscosity \f$\mu=1\f$ and \f$\mathbf f = 0\f$.

The procedure to create the mesh is very simple.<br>
You have to provide to the command line (or via the cfg file) the gmsh.filename option.<br>
You can provide a \c .geo or a \c .msh file (created via gmsh).

As for the laplacian problem, the code is very closed to the mathematical formulation.<br>
We define the product of function spaces for the flow speed and the flow pressur using \c THch<order>()(TH stands for Taylor-Hoods) function which is \c Pch<N+1> \f$\times\f$ \c Pch<N> for respectively flow speed and pressure spaces.<br>
We take an element
\f$U=\left(
    \begin{array}{c}
        u \\
        p \\
    \end{array}
\right)
\f$
in this space. Then we define the integrals of the variationnal formulation for the left and the right hand side. Finally, we apply the Poiseuille profile on the boundary.<br>
We call the solver to resolve the problem (\ref Solver).
\snippet mystokes.cpp marker_main

*/
/// [marker_main]
int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="stokesCerebroVeinous",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    // create the mesh
    //auto mesh = loadMesh(_mesh=new Mesh<Simplex< 2 > > );


    //++++++++++++
    auto mesh = loadGMSHMesh( _mesh=new Mesh<Simplex<3>>,
                              _filename="brainVeines.msh",
                              //_h = meshSizeInit(),
                              _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                              _physical_are_elementary_regions=true,
                              _partitions=Environment::worldComm().localSize(),
                              _worldcomm=Environment::worldComm() );
    //double meshSize = option(_name="gmsh.hsize").template as<double>();
    //std::cout<<"Mesh size = "<<mesh->h()<< "\n";
    /*mesh->addMarkerName( "inlet", 52, 2 );
    mesh->addMarkerName( "outlet", 53, 2 );
    mesh->addMarkerName( "wall", 54, 2 );*/
    //++++++++++++


    // function space
    auto Vh = THch<2>( mesh );

    // element U=(u,p) in Vh
    auto U = Vh->element();
    auto u = U.element<0>();
    auto p = U.element<1>();

    int i;
    //+++++++++++++++++++++++++ Velocity at inlet +++++++++++++++++++++++++++
    auto sumAreaIN = integrate(markedfaces(mesh,1),cst(1.)).evaluate()(0,0);
    std::cout << "Area of inlet 1 = "<< sumAreaIN << "\n";
    double Area;
    std::cout.precision(20);
    for (i=2;i<=31;i++)
        {
            if( (i!=7) && (i!=24))
               {
                   Area = integrate(markedfaces(mesh,i),cst(1.)).evaluate()(0,0);
                   sumAreaIN += Area;
                   std::cout << "Area of inlet "<<i<< " = "<< Area << "\n";
               }
        }
    std::cout << "The total sum of the inlet area = "<<sumAreaIN << "\n";
    auto u_inlet = -6000./sumAreaIN;
    std::cout << "u_inlet = "<< u_inlet << "\n";

    auto flowIn = integrate(markedfaces(mesh,1),inner(u_inlet*N(),N() ) ).evaluate()(0,0);
    for (i=2;i<=31;i++)
        {
            if( (i!=7) && (i!=24))
                flowIn += integrate(markedfaces(mesh,i),inner(u_inlet*N(),N() ) ).evaluate()(0,0);
        }
    std::cout<<"flow_inlet = "<<flowIn<<"\n";

    //+++++++++++++++++++++++++ Velocity at Outlet +++++++++++++++++++++++++++
    auto sumAreaOUT = integrate(markedfaces(mesh,7),cst(1.)).evaluate()(0,0) + integrate(markedfaces(mesh,24),cst(1.)).evaluate()(0,0);
    std::cout << "The total sum of the outlet area = "<<sumAreaOUT << "\n";
    auto u_outlet = 6000./sumAreaOUT;
    std::cout << "u_outlet = "<< u_outlet << "\n";
    auto flowOut = integrate(markedfaces(mesh,7),inner(u_outlet*N(),N() ) ).evaluate()(0,0) + integrate(markedfaces(mesh,24),inner(u_outlet*N(),N() ) ).evaluate()(0,0);
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
    //++++ On the wall ++++
    a+=on(_range=markedfaces(mesh,50),  _rhs=l, _element=u, _expr=vec(cst(0.), cst(0.), cst(0.)) );
    std::cout << "Start solving the system \n";

    //++++ At the inlets ++++
    for (i=1;i<=31;i++)
        {
            if( (i!=7) && (i!=24))
                a+=on(_range=markedfaces(mesh,i), _rhs=l, _element=u, _expr=vec(u_inlet*Nx(),u_inlet*Ny(),u_inlet*Nz()));
        }
    //++++ At the outlets ++++
    a+=on(_range=markedfaces(mesh,7),  _rhs=l, _element=u, _expr=vec(u_outlet*Nx(),u_outlet*Ny(),u_outlet*Nz()) );
    a+=on(_range=markedfaces(mesh,24), _rhs=l, _element=u, _expr=vec(u_outlet*Nx(),u_outlet*Ny(),u_outlet*Nz()) );
    

    //+++++++++ Solve +++++++++
    a.solve(_rhs=l,_solution=U);

    //Flow at inlet
    auto flowInCal = integrate(markedfaces(mesh,1),inner(idv(u),N() ) ).evaluate()(0,0);
    for (i=2;i<=31;i++)
        {
            if( (i!=7) && (i!=24))
                flowInCal += integrate(markedfaces(mesh,i),inner(idv(u),N() ) ).evaluate()(0,0);
        }
 
    //Flow at the outlet
    auto flowOutCal = integrate(markedfaces(mesh,7),inner(idv(u),N() ) ).evaluate()(0,0) + integrate(markedfaces(mesh,24),inner(idv(u),N() ) ).evaluate()(0,0);
    std::cout.precision(20);
    std::cout<<"flow_inlet calcule = "<<flowInCal<<"\n";
    std::cout<<"flow_outlet calcule = "<<flowOutCal<<"\n";
    std::cout<<"Diff flow In/Out = "<<flowInCal+flowOutCal<<"\n";
    std::cout<<"Diff flow In/In = "<<flowInCal-flowIn<<"\n";
    std::cout<<"Diff flow Out/Out = "<<flowOutCal-flowOut<<"\n";

    // LOG(INFO) << "   hsize = " << meshSizeInit() << "\n";
    LOG(INFO) << "[dof]             number of dof: " << Vh->nDof() << "\n";
    LOG(INFO) << "[dof]        number of dof/proc: " << Vh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]          number of dof(U): " << Vh->template functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof]     number of dof/proc(U): " << Vh->template functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]          number of dof(P): " << Vh->template functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof]     number of dof/proc(P): " << Vh->template functionSpace<1>()->nLocalDof()  << "\n";

    //std::cout << "   hsize = " << meshSizeInit() << "\n";
    std::cout << "[dof]             number of dof: " << Vh->nDof() << "\n";
    std::cout << "[dof]        number of dof/proc: " << Vh->nLocalDof() << "\n";
    std::cout << "[dof]          number of dof(U): " << Vh->template functionSpace<0>()->nDof()  << "\n";
    std::cout << "[dof]     number of dof/proc(U): " << Vh->template functionSpace<0>()->nLocalDof()  << "\n";
    std::cout << "[dof]          number of dof(P): " << Vh->template functionSpace<1>()->nDof()  << "\n";
    std::cout << "[dof]     number of dof/proc(P): " << Vh->template functionSpace<1>()->nLocalDof()  << "\n";


    // save results
    auto e = exporter( _mesh=mesh );
   e->add( "u", u );
   e->add( "p", p );
   e->save();
}
/// [marker_main]
