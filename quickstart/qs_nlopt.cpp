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
#include <functional>

#include <feel/feel.hpp>
#include <feel/feelopt/nlopt.hpp>

namespace Feel
{
typedef double (*fobj)(unsigned n, const double *x, double *fgrad, void *my_func_data);

struct myf
{
    decltype( loadMesh(_mesh=new Mesh<Simplex<2>>) ) mesh;
    meta::Pch<Mesh<Simplex<2>>,2>::ptrtype Vh;
    double mu;
    decltype(form2( _trial=Vh, _test=Vh)) a;
    decltype(form1( _test=Vh)) l;

    myf()
        :
        mesh( loadMesh(_mesh=new Mesh<Simplex<2>>) ),
        Vh( Pch<2>( mesh ) ),
        mu(doption(_name="mu")),
        a( form2( _trial=Vh, _test=Vh) ),
        l( form1( _test=Vh) )
        {
        }

    double operator()(unsigned n, const double *x, double *fgrad, void *my_func_data)
        {
            tic();
            auto u = Vh->element("u");
            auto g = expr( soption(_name="functions.g") );
            auto v = Vh->element( g, "g" );
            l.zero();
            l = integrate(_range=elements(mesh),
                          _expr=(Px()/(1+x[1]*x[0]))*id(v));
            a.zero();
            a = integrate(_range=elements(mesh),
                          _expr=x[0]*gradt(u)*trans(grad(v)) +idt(u)*id(u) );
            a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
            a.solve(_rhs=l,_solution=u);
            auto m = mean( _range=elements(mesh), _expr=idv(u) )(0,0);
            toc( "[call myf::operator()]" );
            if ( Environment::isMasterRank() )
                std::cout << "objective function = " << m
                          << " at (" << x[0] << "," << x[1]
                          << ") = " << "\n";
            return m;
        }

};
}

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	po::options_description nloptoptions( "NLOpt options" );
	nloptoptions.add_options()
    ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=nloptoptions,
                   _about=about(_name="qs_nlopt",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));


    ::nlopt::opt opt(::nlopt::LN_NEWUOA_BOUND, 2);

    std::vector<double> lb(2),ub(2);
    lb[0] = 1e-3;lb[1] = 1e-3;
    ub[0] = 10;ub[1] = 10;
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);


    myf F;
    // we don't want to copy myf around so use std::ref()
    ::nlopt::func f = std::ref(F);
    opt.set_min_objective( f, NULL);
    opt.set_maxeval( ioption("nlopt.maxeval") );

    opt.set_xtol_rel(1e-4);

    std::vector<double> x(2);
    x[0] = doption(_name="mu");
    x[1] = doption(_name="mu");
    double minf;

    ::nlopt::result result = opt.optimize(x, minf);
    if (result < 0) {
        std::cout << "nlopt failed!\n";
    }
    else {
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = " << minf << "\n";
    }
#if 0
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", u );
    e->add( "g", v );
    e->save();
#endif
    return 0;
    //# endmarker4 #
}
