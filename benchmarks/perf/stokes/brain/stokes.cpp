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
bool display_flowrates( std::ostream& os, std::map<std::string, double> const& flowrates );

int main(int argc, char**argv )
{
    using namespace Feel;
    po::options_description cerebroveinusoptions( "Cerebroveinus options" );
    cerebroveinusoptions.add_options()
    ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
    ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
    ( "Q", po::value<double>()->default_value( 0.0000052166503023447670669 ), "Flow rate" )
    ;
    
    Environment env( _argc=argc, _argv=argv,
                     _desc=cerebroveinusoptions,
                     _about=about(_name="stokesCerebroVeinous",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    /*auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>>,
                          _filename="brainVeines20.msh",
                          _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                          _physical_are_elementary_regions=true,
                          _partitions=Environment::worldComm().localSize(),
                          _worldcomm=Environment::worldComm()
                         );*/
    constexpr int dim = FEELPP_DIM;
    constexpr int p_order = FEELPP_ORDER_P;
    auto mesh = loadMesh( new Mesh<Simplex<3>> );
    size_type nbdyfaces = nelements(boundaryfaces(mesh));


    // function space
    auto Vh = THch<2>( mesh );

    // element U=(u,p) in Vh
    auto U = Vh->element();
    auto u = U.element<0>();
    auto p = U.element<1>();
    double mu = doption(_name="mu");
    double Q = doption(_name="Q");
    
    if ( Environment::isMasterRank() )
    {
        std::cout << "      number of elements : " << mesh->numGlobalElements() << std::endl;
        std::cout << "         number of faces : " << mesh->numGlobalFaces() << std::endl;
        std::cout << "number of boundary faces : " << nbdyfaces << std::endl;
        std::cout << "      number of points : " << mesh->numGlobalPoints() << std::endl;
        std::cout << "    number of vertices : " << mesh->numGlobalVertices() << std::endl;
        std::cout << " - mesh characteristic" << std::endl;
        std::cout << "                h max : " << mesh->hMax() << std::endl;
        std::cout << "                h min : " << mesh->hMin() << std::endl;
        std::cout << "                h avg : " << mesh->hAverage() << std::endl;
        std::cout << "              measure : " << mesh->measure() << std::endl;

        std::cout << "[dof]             number of dof: " << Vh->nDof() << "\n";
        std::cout << "[dof]        number of dof/proc: " << Vh->nLocalDof() << "\n";
        std::cout << "[dof]          number of dof(U): " << Vh->functionSpace<0>()->nDof()  << "\n";
        std::cout << "[dof]     number of dof/proc(U): " << Vh->functionSpace<0>()->nLocalDof()  << "\n";
        std::cout << "[dof]          number of dof(P): " << Vh->functionSpace<1>()->nDof()  << "\n";
        std::cout << "[dof]     number of dof/proc(P): " << Vh->functionSpace<1>()->nLocalDof()  << "\n";
    }

    std::vector<std::string> flowrates_str={"INLET","OUTLET","inlet1","inlet2","inlet3","inlet4","inlet5","inlet6","inlet7","inlet8","inlet9","inlet10","inlet11","inlet12","inlet13","inlet14","inlet15","inlet16","inlet17","inlet18","inlet19","inlet20","inlet21","inlet22","inlet23","inlet24","inlet25","inlet26","inlet27","inlet28","inlet29","outlet1","outlet2"};
    auto f_flowIn = form1(_test=Vh->functionSpace<0>() );
    std::map<std::string,form1_type<functionspace_type<decltype(Vh->functionSpace<0>())>>> flowrates_f;
    for( auto const& f : flowrates_str )
    {
        flowrates_f[f] = form1(_test=Vh->functionSpace<0>() );
        flowrates_f[f] = integrate(_range=markedfaces(mesh,f), _expr=inner(id(u),N()));
    }
    std::ofstream ofs( "flowrates.md" );
    std::map<std::string,double> flowrates;
    for( auto & f : flowrates_f )
    {
        flowrates[f.first]= 0;;
    }
    //display_flowrates( ofs, flowrates );

    auto hsize = mean(elements(mesh),h());
    //+++++++++++++++++++++++++ Velocity at inlet +++++++++++++++++++++++++++
    auto sumAreaIN = integrate(markedfaces(mesh,"INLET"),cst(1.)).evaluate()(0,0); //52 is a physical group assembling all the entries
    auto u_inlet = -Q/sumAreaIN;
    auto flowIn = integrate(markedfaces(mesh,"INLET"),inner(u_inlet*N(),N() ) ).evaluate()(0,0);

    //+++++++++++++++++++++++++ Velocity at Outlet +++++++++++++++++++++++++++
    auto sumAreaOUT = integrate(markedfaces(mesh,"OUTLET"),cst(1.)).evaluate()(0,0);
    auto u_outlet = Q/sumAreaOUT;
    auto flowOut = integrate(markedfaces(mesh,"OUTLET"),inner(u_outlet*N(),N() ) ).evaluate()(0,0) ;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( Environment::isMasterRank() )
    {
        std::cout << "   mean hsize = " << hsize << "\n";
        std::cout << "The total sum of the inlet area = "<<sumAreaIN << "\n";
        std::cout << "u_inlet = "<< u_inlet << "\n";
        std::cout<<"flow_inlet = "<<flowIn<<"\n";
        std::cout << "The total sum of the outlet area = "<<sumAreaOUT << "\n";
        std::cout << "u_outlet = "<< u_outlet << "\n";
        std::cout<<"flow_outlet = "<<flowOut<<"\n";
    }
    // left hand side
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=mu*trace(gradt(u)*trans(grad(u))) );

    a+= integrate(_range=elements(mesh),
                  _expr=-div(u)*idt(p)-divt(u)*id(p));

    // right hand side
    auto f = vec(cst(0.),cst(0.),cst(0.));
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh), _expr=inner(f,id(u)));
    std::cout << "Start imposing boundary conditions \n";


    //++++++++ Boundary Conditions ++++++++++++++++
    //++++ At the inlets ++++
    a+=on(_range=markedfaces(mesh,"INLET"), _rhs=l, _element=u, _expr=vec(u_inlet*Nx(),u_inlet*Ny(),u_inlet*Nz()));
    //++++ On the wall ++++
    a+=on(_range=markedfaces(mesh,"wall"),  _rhs=l, _element=u, _expr=vec(cst(0.), cst(0.), cst(0.)) );
    std::cout << "Start solving the system \n";

    /*//++++ At the outlets ++++
    a+=on(_range=markedfaces(mesh,7),  _rhs=l, _element=u, _expr=vec(u_outlet*Nx(),u_outlet*Ny(),u_outlet*Nz()) );
    a+=on(_range=markedfaces(mesh,24), _rhs=l, _element=u, _expr=vec(u_outlet*Nx(),u_outlet*Ny(),u_outlet*Nz()) );*/

    //+++++++++ Solve +++++++++
    a.solve(_rhs=l,_solution=U);
#if 0
    //Flow at inlet
    auto flowInCal = integrate(markedfaces(mesh,"INLET"),inner(idv(u),N() ) ).evaluate()(0,0);
    //Flow at the outlet
    auto flowOutCal = integrate(markedfaces(mesh,"OUTLET"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    //Flow at each of the inlet sections
    auto flowIn1 = integrate(markedfaces(mesh,"inlet1"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn2 = integrate(markedfaces(mesh,"inlet2"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn3 = integrate(markedfaces(mesh,"inlet3"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn4 = integrate(markedfaces(mesh,"inlet4"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn5 = integrate(markedfaces(mesh,"inlet5"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn6 = integrate(markedfaces(mesh,"inlet6"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn7 = integrate(markedfaces(mesh,"inlet7"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn8 = integrate(markedfaces(mesh,"inlet8"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn9 = integrate(markedfaces(mesh,"inlet9"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn10 = integrate(markedfaces(mesh,"inlet10"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn11 = integrate(markedfaces(mesh,"inlet11"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn12 = integrate(markedfaces(mesh,"inlet12"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn13 = integrate(markedfaces(mesh,"inlet13"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn14 = integrate(markedfaces(mesh,"inlet14"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn15 = integrate(markedfaces(mesh,"inlet15"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn16 = integrate(markedfaces(mesh,"inlet16"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn17 = integrate(markedfaces(mesh,"inlet17"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn18 = integrate(markedfaces(mesh,"inlet18"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn19 = integrate(markedfaces(mesh,"inlet19"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn20 = integrate(markedfaces(mesh,"inlet20"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn21 = integrate(markedfaces(mesh,"inlet21"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn22 = integrate(markedfaces(mesh,"inlet22"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn23 = integrate(markedfaces(mesh,"inlet23"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn24 = integrate(markedfaces(mesh,"inlet24"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn25 = integrate(markedfaces(mesh,"inlet25"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn26 = integrate(markedfaces(mesh,"inlet26"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn27 = integrate(markedfaces(mesh,"inlet27"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn28 = integrate(markedfaces(mesh,"inlet28"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowIn29 = integrate(markedfaces(mesh,"inlet29"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    //Flow at each of the outlet sections
    auto flowOut1 = integrate(markedfaces(mesh,"outlet1"),inner(idv(u),N() ) ).evaluate()(0,0) ;
    auto flowOut2 = integrate(markedfaces(mesh,"outlet2"),inner(idv(u),N() ) ).evaluate()(0,0) ;

    if ( Environment::isMasterRank() )
    {
        std::cout.precision(20);
        std::cout<<"flow_inlet calcule = "<<flowInCal<<"\n";
        std::cout<<"flow_outlet calcule = "<<flowOutCal<<"\n";
        std::cout<<"Diff flow In/Out = "<<flowInCal+flowOutCal<<"\n";
        std::cout<<"Diff flow In/In = "<<flowInCal-flowIn<<"\n";
        std::cout<<"Diff flow Out/Out = "<<flowOutCal-flowOut<<"\n";

        std::cout<<"flow inlet1: = "<<flowIn1<<"\n";
        std::cout<<"flow inlet2: = "<<flowIn2<<"\n";
        std::cout<<"flow inlet3: = "<<flowIn3<<"\n";
        std::cout<<"flow inlet4: = "<<flowIn4<<"\n";
        std::cout<<"flow inlet5: = "<<flowIn5<<"\n";
        std::cout<<"flow inlet6: = "<<flowIn6<<"\n";
        std::cout<<"flow inlet7: = "<<flowIn7<<"\n";
        std::cout<<"flow inlet8: = "<<flowIn8<<"\n";
        std::cout<<"flow inlet9: = "<<flowIn9<<"\n";
        std::cout<<"flow inlet10: = "<<flowIn10<<"\n";
        std::cout<<"flow inlet11: = "<<flowIn11<<"\n";
        std::cout<<"flow inlet12: = "<<flowIn12<<"\n";
        std::cout<<"flow inlet13: = "<<flowIn13<<"\n";
        std::cout<<"flow inlet14: = "<<flowIn14<<"\n";
        std::cout<<"flow inlet15: = "<<flowIn15<<"\n";
        std::cout<<"flow inlet16: = "<<flowIn16<<"\n";
        std::cout<<"flow inlet17: = "<<flowIn17<<"\n";
        std::cout<<"flow inlet18: = "<<flowIn18<<"\n";
        std::cout<<"flow inlet19: = "<<flowIn19<<"\n";
        std::cout<<"flow inlet20: = "<<flowIn20<<"\n";
        std::cout<<"flow inlet21: = "<<flowIn21<<"\n";
        std::cout<<"flow inlet22: = "<<flowIn22<<"\n";
        std::cout<<"flow inlet23: = "<<flowIn23<<"\n";
        std::cout<<"flow inlet24: = "<<flowIn24<<"\n";
        std::cout<<"flow inlet25: = "<<flowIn25<<"\n";
        std::cout<<"flow inlet26: = "<<flowIn26<<"\n";
        std::cout<<"flow inlet27: = "<<flowIn27<<"\n";
        std::cout<<"flow inlet28: = "<<flowIn28<<"\n";
        std::cout<<"flow inlet29: = "<<flowIn29<<"\n";

        std::cout<<"flow outlet1: = "<<flowOut1<<"\n";
        std::cout<<"flow outlet2: = "<<flowOut2<<"\n";
    }
#endif


    for( auto & f : flowrates_f )
        flowrates[f.first]= f.second(u);

    if (Environment::isMasterRank())
    {
        std::cout<< "\n \n \n ============= FLOW OVER RADIAL SECTIONS =================\n";
        bool ok = display_flowrates( std::cout, flowrates);
        if ( !ok )
            std::cout << "WARNING ! flow rate issue at time \n";
        display_flowrates( ofs,flowrates );
        std::cout<< "\n \n \n ==========================================================\n";
    }
    // save results
    auto e = exporter( _mesh=mesh );
   e->add( "u", u );
   e->add( "p", p );
   e->save();
}
/// [marker_main]
