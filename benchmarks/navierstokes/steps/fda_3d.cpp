/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date     : Tue Feb 25 12:13:15 2014

 Copyright (C) 2014-2015 Feel++ Consortium

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
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feelpde/preconditionerblockns.hpp>

int main(int argc, char**argv )
{
    constexpr int dim = FEELPP_DIM;
    constexpr int order_p= FEELPP_ORDER_P;
    
    using namespace Feel;
	po::options_description steadyfdaoptions( "Steady FDA options" );
	steadyfdaoptions.add_options()
		( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "rho", po::value<double>()->default_value( 1.0 ), "coeff" )
        ( "penaldir", po::value<double>()->default_value( 100 ), "coeff" )
        ( "sym", po::value<bool>()->default_value( 0 ), "use symmetric deformation tensor" )
        ( "stokes.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Stokes preconditioner: petsc, PM, Blockns" )
        ( "picard", po::value<bool>()->default_value( 1 ), "picard" )
        ( "picard.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Stokes preconditioner: petsc, PM, Blockns" )
		( "picard.tol", po::value<double>()->default_value( 1e-8 ), "tolerance" )
		( "picard.maxit", po::value<double>()->default_value( 10 ), "max iteration" )
        ( "newton", po::value<bool>()->default_value( 1 ), "newton" )
        ( "newton.preconditioner", po::value<std::string>()->default_value( "petsc" ), "Stokes preconditioner: petsc, PM, Blockns" )
        ( "newton.tol", po::value<double>()->default_value( 1e-8 ), "tolerance" )
		( "newton.maxit", po::value<double>()->default_value( 10 ), "max iteration" )
		;
    steadyfdaoptions.add( backend_options( "stokes" ) )
        .add( backend_options( "newton" ) )
         .add( backend_options( "picard" ) );
	Environment env( _argc=argc, _argv=argv,
                     _desc=steadyfdaoptions,
                     _about=about(_name=(boost::format("steady_fda_%1%d")%dim).str(),
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<dim>>);
    auto Vh = THch<order_p>( mesh );
    auto U = Vh->element();
    auto Un = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto un = Un.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto pn = Un.element<1>();
    auto q = V.element<1>();
    double mu = doption(_name="mu");
    double rho = doption(_name="rho");

    if ( Environment::isMasterRank() )
    {
        std::cout << "Re\tFunctionSpace\tVelocity\tPressure\n";
        std::cout.width(16);
        std::cout << std::left << 0.0461*0.012*rho/mu;
        std::cout.width(16);
        std::cout << std::left << Vh->nDof();
        std::cout.width(16);
        std::cout << std::left << Vh->functionSpace<0>()->nDof();
        std::cout.width(16);
        std::cout << std::left << Vh->functionSpace<1>()->nDof() << "\n";

        std::cout << "[blockns]\n";
        std::cout << " - cd: " << boption( "blockns.cd" ) << "\n";
        std::cout << " - pcd: " << boption( "blockns.pcd" ) << "\n";
        std::cout << " - pcd.inflow: " << soption( "blockns.pcd.inflow" ) << "\n";
        std::cout << " - pcd.outflow: " << soption( "blockns.pcd.outflow" ) << "\n";
        std::cout << " - pcd.order: " << ioption( "blockns.pcd.order" ) << "\n";

    }
    auto deft = gradt( u );
    auto def = grad( v );
    
    auto mybdf = bdf( _space=Vh, _name="mybdf" );
    
    double fixPtTol = doption(_name="picard.tol");
    int fixPtMaxIt = doption(_name="picard.maxit");
    double newtonTol = doption(_name="newton.tol");
    int newtonMaxIt = doption(_name="newton.maxit");

    
    BoundaryConditions bcs;
    map_vector_field<dim,1,2> m_dirichlet { bcs.getVectorFields<dim> ( "velocity", "Dirichlet" ) };

    auto l = form1( _test=Vh );
    auto r = form1( _test=Vh );
    auto residual = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh);
    auto at = form2( _trial=Vh, _test=Vh);
    a += integrate( _range=elements( mesh ), _expr=mu*inner( deft,def ) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );
    if ( boption("blockns.weakdir") )
    {
        for( auto const& d : m_dirichlet )
        {
            LOG(INFO) << "Setting Dirichlet condition on " << d.first << " with " << d.second;
            a += integrate( _range=markedfaces( mesh, d.first ),
                            _expr=(trans(-mu*gradt(u)*N()+idt(p)*N()))*id(v)+(trans(-mu*grad(u)*N()+id(p)*N()))*idt(v)+doption("penaldir")*trans(idt(u))*id(v)/hFace() );
        }
    }
    if ( Environment::numberOfProcessors() == 1 )
    {
        a.matrix().printMatlab( "A.m" );
    }
    auto e = exporter( _mesh=mesh );

    
    auto incru = normL2( _range=elements(mesh), _expr=idv(u)-idv(un));
    auto incrp = normL2( _range=elements(mesh), _expr=idv(p)-idv(pn));
    at+=a;
    if ( boption("blockns.weakdir") )
    {
        for( auto const& d : m_dirichlet )
        {
            LOG(INFO) << "Setting Dirichlet condition on " << d.first << " with " << d.second;
            l += integrate( _range=markedfaces( mesh, d.first ),
                            _expr=trans(d.second)*(-mu*grad(u)*N()+id(p)*N()+doption("penaldir")*id(v)/hFace() ) );
        }
    }
    else
    {
        for( auto const& c : m_dirichlet )
        {
            at+=on(_range=markedfaces(mesh,c.first), _rhs=l, _element=u,
                   _expr=c.second ) ;
        }
    }

    tic();
    auto a_blockns = blockns( _space=Vh,
                             _type=soption("stokes.preconditioner"),
                             _bc=bcs,
                             _matrix= at.matrixPtr(),
                             _alpha=rho*mybdf->polyDerivCoefficient(0),
                             _mu=mu,
                             _rho=rho,
                             _prefix="velocity" );

    
    toc(" - Setting up Precondition Blockns...");
    
    a_blockns->setMatrix( at.matrixPtr() );

    tic();
    
    
    if ( soption("stokes.preconditioner") != "petsc" )
    {
        a_blockns->update( at.matrixPtr(), zero<dim,1>(), m_dirichlet );
        at.solveb(_rhs=l,_solution=U,_backend=backend(_name="stokes"),_prec=a_blockns);
    }
    else
        at.solve(_rhs=l,_solution=U);
    toc(" - Solving Stokes...");

    tic();
    e->step(0)->add( "u", u );
    e->step(0)->add( "p", p );
    e->save();
    toc(" - Exporting Stokes results...");
    
    Un = U;
    for ( mybdf->start();  mybdf->isFinished() == false; mybdf->next(U) )
    {
        int fixedpt_iter = 0;

        double res = 0;
        
        auto deltaU = Vh->element();
        auto bdf_poly = mybdf->polyDeriv();
        auto rhsu =  bdf_poly.element<0>();
        auto extrap = mybdf->poly();
        auto extrapu = extrap.element<0>();

        tic();
        r.zero();
        at.zero();
        // add BDF term to the right hand side from previous time steps
        r = integrate( _range=elements(mesh), _expr=rho*(trans(idv(rhsu))*id(u) ) );
        
        at += a;
        at += integrate( _range=elements(mesh),_expr=rho*trans(id(v))*(gradt(u)*idv(extrapu)) );
        toc( " - Picard:: Assemble nonlinear terms  ..." );
        tic();
        
        for( auto const& c : m_dirichlet )
        {
            LOG(INFO) << "Applying Dirichlet condition " << c.second << " on " << c.first;
            if ( boption( "blockns.weakdir" ) )
            {
                at += integrate( _range=markedfaces( mesh, c.first ),
                                    _expr=(trans(-mu*gradt(u)*N()+idt(p)*N()))*id(v)+(trans(-mu*grad(u)*N()+id(p)*N()))*idt(v+doption("penaldir")*trans(idt(u))*id(v)/hFace() );
                r += integrate( _range=markedfaces( mesh, c.first ),
                                    _expr=trans(-mu*grad(u)*N()+id(p)*N()+doption("penaldir")*id(v)/hFace())*c.second );
            }
            else
            {
                at+=on(_range=markedfaces(mesh,c.first), _rhs=r, _element=u,
                           _expr=c.second ) ;
            }
        }
        toc(" - Picard:: Assemble BC   ...");

            
            
        /*u=vf::project(Vh->template functionSpace<0>(), elements(mesh), zero<3,1>());
        p=vf::project(Vh->template functionSpace<1>(), elements(mesh), constant(0.0));*/
            
            
        if ( soption("picard.preconditioner") != "petsc" )
        {
            a_blockns->update( at.matrixPtr(), idv(extrapu), m_dirichlet );
            at.solveb(_rhs=r,_solution=U,_backend=backend(_name="picard",_rebuild=true),_prec=a_blockns );
        }
        else
            at.solveb(_rhs=r,_solution=U,_backend=backend(_rebuild=true) );
        
        toc(" - Picard:: Solve   ...");
        incru = normL2( _range=elements(mesh), _expr=idv(u)-idv(un));
        incrp = normL2( _range=elements(mesh), _expr=idv(p)-idv(pn));
        fixedpt_iter++;
            
        if ( Environment::isMasterRank() )
        {
            std::cout << "Iteration "  << fixedpt_iter << "\n";
            std::cout << " . ||u-un|| = " << incru << std::endl;
            std::cout << " . ||p-pn|| = " << incrp << std::endl;
        }
        Un = U;
        Un.close();

        e->step(mybdf->time())->add( "u", u );
        e->step(mybdf->time())->add( "p", p );
        auto mean_p = mean( _range=elements(mesh), _expr=idv(p) )( 0, 0 );
        e->step(mybdf->time())->addScalar( "mean_p", mean_p );
        e->save();

    }
return 0;
}
