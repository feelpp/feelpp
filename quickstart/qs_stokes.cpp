/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 12:13:15 2014

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
#include <feel/feel.hpp>
#include <feel/feelpde/preconditionerbtpcd.hpp>

int main(int argc, char**argv )
{
    //! [marker1]
    using namespace Feel;
	po::options_description stokesoptions( "Stokes options" );
	stokesoptions.add_options()
		( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
		( "fixpoint.tol", po::value<double>()->default_value( 1e-8 ), "tolerance" )
		( "fixpoint.maxit", po::value<double>()->default_value( 10 ), "max iteration" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc=stokesoptions,
                     _about=about(_name="qs_stokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

#if 0
    double meshSize = option(_name="gmsh.hsize").as<double>();
    GeoTool::Rectangle R( meshSize,"myRectangle",GeoTool::Node(0,0),GeoTool::Node(5,1));
    R.setMarker(_type="line",_name="inlet",_marker4=true);
    R.setMarker(_type="line",_name="outlet",_marker2=true);
    R.setMarker(_type="line",_name="wall",_marker1=true,_marker3=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);
    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2>>,_name="qs_stokes");
#else
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
#endif



    auto g = expr<2,1>( soption(_name="functions.g") );
    auto wall = expr<2,1>( soption(_name="functions.h") );
    auto Vh = THch<2>( mesh );
    auto U = Vh->element();
    auto Un = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto un = Un.element<0>();
    auto v = V.element<0>(g,"poiseuille");
    auto p = U.element<1>();
    auto pn = Un.element<1>();
    auto q = V.element<1>();

    if ( Environment::isMasterRank() )
    {
        std::cout << "FunctionSpace\tVelocity\tPressure\n";
        std::cout.width(16);
        std::cout << std::left << Vh->nDof();
        std::cout.width(16);
        std::cout << std::left << Vh->functionSpace<0>()->nDof();
        std::cout.width(16);
        std::cout << std::left << Vh->functionSpace<1>()->nDof() << "\n";
    }
    auto deft = gradt( u );
    auto def = grad( v );
    double mu = doption(_name="mu");
    double fixPtTol = doption(_name="fixpoint.tol");
    int fixPtMaxIt = doption(_name="fixpoint.maxit");
    auto l = form1( _test=Vh );
    auto r = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh);
    auto at = form2( _trial=Vh, _test=Vh);
    a += integrate( _range=elements( mesh ), _expr=mu*inner( deft,def ) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) + divt( u )*id( q ) );

    auto e = exporter( _mesh=mesh );

    std::map<std::string,std::set<flag_type>> bcs;
    bcs["Dirichlet"].insert(mesh->markerName("inlet"));
    bcs["Dirichlet"].insert(mesh->markerName("wall"));
    bcs["Neumann"].insert(mesh->markerName("outlet"));
    
    auto incru = normL2( _range=elements(mesh), _expr=idv(u)-idv(un));
    auto incrp = normL2( _range=elements(mesh), _expr=idv(p)-idv(pn));
    int fixedpt_iter = 0;
    at=a;
    at+=on(_range=markedfaces(mesh,"wall"), _rhs=l, _element=u,
           _expr=zero<2,1>() ) ;
    at+=on(_range=markedfaces(mesh,"inlet"), _rhs=l, _element=u,
           _expr=g );

    if ( Environment::isMasterRank() )
        std::cout << " - Setting up Precondition BtPCD...\n";
    auto a_btpcd = btpcd( _space=Vh, _bc=bcs, _matrix= at.matrixPtr() );
    a_btpcd->setMatrix( at.matrixPtr() );
    map_vector_field<2,1,2> m_dirichlet;
    m_dirichlet["inlet"]=g;
    m_dirichlet["wall"]=wall;

    if ( Environment::isMasterRank() )
        std::cout << " - Solving Stokes...\n";
    if ( boption("btpcd") )
    {
        a_btpcd->update( at.matrixPtr(), zero<2,1>(), m_dirichlet );
        at.solveb(_rhs=l,_solution=U,_backend=backend(),_prec=a_btpcd);
    }
    else
        at.solve(_rhs=l,_solution=U);

    e->step(0)->add( "u", u );
    e->step(0)->add( "p", p );
    e->save();
#if 1
    auto deltaU = Vh->element();
    //m_dirichlet["inlet"]=wall;
    do
    {
        if ( Environment::isMasterRank() )
            std::cout << " - Assemble nonlinear terms  ...\n";
        at = a;
        at += integrate( _range=elements(mesh),_expr=trans(id(v))*(gradt(u)*idv(u)) );
        at += integrate( _range=elements(mesh), _expr=trans(id(v))*gradv(u)*idt(u) );

        if ( Environment::isMasterRank() )
            std::cout << " - Assemble residual  ...\n";
        r = integrate( _range=elements( mesh ), _expr=mu*inner( gradv(u),def ) );
        r +=integrate( _range=elements( mesh ), _expr=-div( v )*idv( p ) + divv( u )*id( q ) );
        r += integrate( _range=elements(mesh),_expr=trans(id(v))*(gradv(u)*idv(u)) );
        if ( Environment::isMasterRank() )
            std::cout << " - Assemble BC   ...\n";
        at+=on(_range=markedfaces(mesh,"wall"), _rhs=r, _element=u,
               _expr=zero<2,1>() ) ;
        at+=on(_range=markedfaces(mesh,"inlet"), _rhs=r, _element=u,
               _expr=zero<2,1>() );
        r.vectorPtr()->scale(-1);
        if ( Environment::isMasterRank() )
        {
            std::cout << "non linear iteration " << fixedpt_iter << " \n";
        }
        if ( Environment::isMasterRank() )
            std::cout << " - Solve   ...\n";
        if ( boption("btpcd") )
        {
            a_btpcd->update( at.matrixPtr(), idv(u), m_dirichlet );
            at.solveb(_rhs=r,_solution=deltaU/*U*/,_backend=backend(),_prec=a_btpcd );
        }
        else
            at.solveb(_rhs=r,_solution=deltaU/*U*/,_backend=backend(_rebuild=true) );
        U.add(1.,deltaU);
        incru = normL2( _range=elements(mesh), _expr=idv(u)-idv(un));
        incrp = normL2( _range=elements(mesh), _expr=idv(p)-idv(pn));
        fixedpt_iter++;
        if ( Environment::isMasterRank() )
        {
            std::cout << "Iteration "  << fixedpt_iter << "\n";
            std::cout << " . ||u-un|| = " << incru << std::endl;
            std::cout << " . ||p-pn|| = " << incrp << std::endl;
        }
        if  ( fixedpt_iter < 10 )
        {
            Un = U;
            Un.close();
        }
        e->step(fixedpt_iter)->add( "u", u );
        e->step(fixedpt_iter)->add( "p", p );
        e->save();

    }
    while ( ( incru > fixPtTol && incrp > fixPtTol ) && ( fixedpt_iter < fixPtMaxIt ) );
#endif 
    
    return 0;
}
