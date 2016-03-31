/*This file is part of the Feel library

  Author(s): Guillaume Dolle <gdolle@unistra.fr>
       Date: 2016-03-01

  Copyright (C) 2014-2016 Feel++ Consortium

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
#include <functional>

#include <feel/feel.hpp>
#include <feel/feelopt/nlopt.hpp>


constexpr static int ORDER = 1;
// Fixed: do not change.
constexpr static int DIM = 2;
constexpr static int N_UNKNOWNS = 3;


int main(int argc, char**argv )
{
    using namespace Feel;
    po::options_description nloptoptions( "NLOpt options" );
    nloptoptions.add_options()
        ( "k", po::value<double>()->default_value( 1.0 ), "Diffusion coefficient background" )
        ( "k0", po::value<double>()->default_value( 2 ), "Diffusion coefficient inclusion 1 in [0,1]" )
        ( "k1", po::value<double>()->default_value( 3 ), "Diffusion coefficient inclusion 2 in [0,1]" )
        ( "k2", po::value<double>()->default_value( 4 ), "Diffusion coefficient inclusion 3 in [0,1]" )
        ( "k0init", po::value<double>()->default_value( 1.0 ), "Diffusion coefficient inclusion 1 in [0,1]" )
        ( "k1init", po::value<double>()->default_value( 1.0 ), "Diffusion coefficient inclusion 2 in [0,1]" )
        ( "k2init", po::value<double>()->default_value( 1.0 ), "Diffusion coefficient inclusion 3 in [0,1]" )
        ( "lb", po::value<double>()->default_value( 1e-4 ), "Lowerbound value" )
        ( "ub", po::value<double>()->default_value( 1e1 ), "Upperbound value" )
        ;

    Environment env( _argc=argc, _argv=argv,
                     _desc=nloptoptions,
                     _about=about(_name="heat_nlopt",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    using vec = std::vector<double>;

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<DIM>>,
                          _filename="circle_incs.geo" );

    auto g = expr( soption(_name="functions.g") );

    auto Vh = Pch<ORDER>( mesh ) ;
    auto a = form2( _trial=Vh, _test=Vh);
    auto l = form1( _test=Vh);

    auto u = Vh->element("u");
    auto p = Vh->element("p");
    auto v = Vh->element( "v" );

    auto u_obs = Vh->element("u_obs");
    auto k_real = Vh->element("k_real");

    auto k = Vh->element();

    auto e = exporter( mesh );
    int iter=0;

    std::list<std::string> omega = {"homogeneous","inclusion0","inclusion1","inclusion2"};

    auto forwardProblem = [&]( const double* x )
    {
        k.on( _range=markedelements( mesh, omega ), _expr=cst(doption("k")) );
        k.on( _range=markedelements( mesh, "inclusion0" ), _expr=cst(x[0]) );
        k.on( _range=markedelements( mesh, "inclusion1" ), _expr=cst(x[1]) );
        k.on( _range=markedelements( mesh, "inclusion2" ), _expr=cst(x[2]) );

        LOG(INFO) << "Compute forward problem";
        a.zero();
        l.zero();
        l = integrate(_range=markedelements(mesh,omega),
                      _expr=cst(0.));
        a = integrate(_range=markedelements(mesh,omega),
                      _expr=idv(k)*gradt(u)*trans(grad(v)) );
        // Dirichlet boundary condition.
        a+=on( _range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
        a.solve( _rhs=l, _solution=u );

        iter++;
    };

    auto adjointProblem = [&]( const double* x )
    {
        LOG(INFO) << "Compute adjoint problem";
        a.zero();
        l.zero();
        l = integrate( _range=markedelements(mesh,omega),
                       _expr=(idv(u)-idv(u_obs))*id(v) );
        a = integrate( _range=markedelements(mesh,omega),
                       _expr=idv(k)*gradt(p)*trans(grad(v)) );
        // Dirichlet boundary condition.
        a+=on( _range=boundaryfaces(mesh), _rhs=l, _element=p, _expr=cst(0) );
        a.solve( _rhs=l, _solution=p );
    };

    // Objective function.
    auto J = [&]()->double
    {
        // Compute objective function 1/2 * int_U (u-u_obs)^2
        double j = 1./2 * integrate( _range=markedelements(mesh,omega),
                                     _expr=(idv(u) - idv(u_obs))*(idv(u) - idv(u_obs)) ).evaluate()(0,0);
        return j;
    };

    // Gradient
    auto gradJ = [&]( const double* x, double* grad )->double*
    {
        grad[0] = integrate( _range=markedelements(mesh,"inclusion0"),
                             _expr=-gradv(p)*trans(gradv(u)) ).evaluate()(0,0);
        grad[1] = integrate( _range=markedelements(mesh,"inclusion1"),
                             _expr=-gradv(p)*trans(gradv(u)) ).evaluate()(0,0);
        grad[2] = integrate( _range=markedelements(mesh,"inclusion2"),
                             _expr=-gradv(p)*trans(gradv(u)) ).evaluate()(0,0);
        return grad;
    };

    // We could use the vector version here (see ::NLOPT::vfunc).
    auto myfunc = [&]( unsigned n, const double* x, double* grad, void *my_func_data )->double
    {
        forwardProblem(x);
        double obj=J();

        adjointProblem(x);
        gradJ(x, grad);

        Feel::cout << "objective function = " <<  obj
                   << " for diffusion coefficients k0=" << x[0] << ", k1=" << x[1] << ", k2=" << x[2]
                   << " gradient DJk0=" << grad[0] << ", DJk1=" << grad[1] << ", DJk2=" << grad[2]
                   << "\n";
        return obj;
    };

    // We compute a numerical observation which is the solution of the forward problem
    // for real coefficients k0, k1, k2.

    // Diffusion coefficient initialization.
    k.on( _range=markedelements( mesh, omega ), _expr=cst(doption("k")) );
    vec x = { doption("k0"), doption("k1"), doption("k2") };
    forwardProblem( &x[0] );
    // Keep a copy here
    k_real=k;
    u_obs=u;

    // Set default inclusion diffusion coefficients.
    x = { doption("k0init"), doption("k1init"), doption("k2init") };

    ::nlopt::opt opt( ::nlopt::LD_LBFGS, N_UNKNOWNS );

    // We keep diffusion coefficients in [0,1].
    vec lb = { doption("lb"), doption("lb"), doption("lb") }; // lowerbound
    vec ub = { doption("ub"), doption("ub"), doption("ub") }; // upperbound

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    ::nlopt::func f = std::ref( myfunc );

    opt.set_min_objective( f, nullptr );
    opt.set_maxeval( ioption("nlopt.maxeval") );
    opt.set_xtol_rel( 1e-4 );

    double minf;

    //myfunc(0,&x[0],NULL,NULL);

    ::nlopt::result result = opt.optimize(x, minf);
    if (result < 0)
    {
        Feel::cout << "nlopt failed!\n";
    }
    else
    {
        Feel::cout << "Diffusion coefficient found!\n"
                   << "k0=" << x[0] << ", k1=" << x[1] << ", k2=" << x[2]
                   << "\n";
    }

    if(e->doExport())
    {
        e->add( "k_real", k_real );
        e->add( "u_obs", u_obs );
        e->add( "u", u );
        e->add( "k", k );
        e->save();
    }

    return 0;
}

/***********************************************************************[ MODE ]
 -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
 -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
\******************************************************************************/
