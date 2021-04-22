/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 22 Nov 2019

 Copyright (C) 2019 Feel++ Consortium

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
#include <map>
#include <vector>
#include <feel/feelcore/environment.hpp>
#include <feel/feelvf/ginac.hpp>
#include <matplot/matplot.h>



int main( int argc, char** argv )
{
    using namespace Feel;
    using Feel::cout;
    po::options_description meshpartoptions( "Plot options" );
	meshpartoptions.add_options()
        ( "tmin", po::value<double>()->default_value( 0 ), "min value" )
        ( "tmax", po::value<double>()->default_value( 1 ), "max value" )
        ( "tstep", po::value<double>()->default_value( 100 ), "number of value between tmin and tmax" )
        ( "tsamp", po::value<std::string>()->default_value( "linear" ), "type of t sampling: linear or log" )
        ( "type", po::value<std::string>()->default_value( "plot" ), "plot, semilogx, semilogy, loglog" )
        ( "scalar_expr", po::value<std::vector<std::string>>()->default_value( {"t:t"} ), "list of scalar expressions with name and representations" )
        ( "vectorial_expr", po::value<std::vector<std::string>>()->default_value( {"gv|{sin(2*pi*x),sin(2*pi*x),sin(2*pi*x)}:x|nodal|element"} ), "list of vectorial  expressions with name and representations" )

        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=meshpartoptions,
                     _about=about( _name="plot" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" )
                     );

    using namespace matplot;
    auto tsamp = [&]( double tmin, double tmax, double tstep ) {
        if  ( soption("tsamp") == "log" )
            return logspace( tmin, tmax, tstep );
        else
            return linspace( tmin, tmax, tstep );
    };
    std::vector<double> x = tsamp(doption( "tmin" ), doption( "tmax" ), doption( "tstep" ) );
    int n = vsoption( "scalar_expr" ).size();
    int i = 0;
    for( std::string e_str : vsoption( "scalar_expr" ))
    {
        subplot(n,1,i);
        std::vector<double> y = transform(x, [&e_str](auto x) { 
            auto e = expr(e_str);
            e.setParameterValues( { { "t", x } } );
            return e.evaluate()( 0, 0 ); });
        if ( soption("type") == "plot" )
            plot(x, y, "-x");
        if ( soption("type") == "loglog" )
            loglog(x, y, "-x");
        if ( soption("type") == "semilogx" )
            semilogx(x, y, "-x");
        if ( soption( "type" ) == "semilogy" )
            semilogy( x, y, "-x" );
        title(e_str);
        //hold(matplot::on);
        ++i;
    }
    wait();

    return 0;

}

