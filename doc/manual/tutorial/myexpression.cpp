/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Vincent HUBER <vincent.huber@cemosis.fr>

   Date 2013-02-18

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

/// [marker1]
namespace Feel
{
  struct myFunctor
  {
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {
      double a = x[0];
      double b = x[1];
      return _esp*std::min(a,b); 
    }
    void setEps(double _a=0.){ _esp = _a;}
    double _esp;
  };
}
/// [marker1]

int main(int argc, char**argv )
{
/// [marker2]
  po::options_description opts ( "myOptions ");
  opts.add_options()
    ( "_ex_", po::value<std::string>()->default_value( "cos(x)*sin(y)" ), "my Expression to project" );

  using namespace Feel;

  Environment env( _argc=argc, _argv=argv,
      _desc=opts.add( feel_options() ),
      _about=about(_name="myexpression",
        _author="Feel++ Consortium",
        _email="feelpp-devel@feelpp.org"));

  // create mesh
  const int Dim = 2;
  auto mesh = unitSquare();
  //Exporter
  auto e = exporter( _mesh=mesh );

  // function space
  auto Vh = Pch<1>( mesh );
  auto u = Vh->element();
  auto v = Vh->element();

  //Define symbols, ie "x", "y"...
  std::vector<GiNaC::symbol> vars = symbols<Dim>();
  std::string _u = option(_name="_ex_").as<std::string>();

  //Create the expression from the string
  GiNaC::ex expr_u = parse(_u,vars);
  u =  vf::project( _space=Vh, _range=elements(mesh), _expr=expr(expr_u,vars) );
/// [marker2]

  // export results
  e->add( "u", u );

/// [marker3]
  Feel::myFunctor _f;
 _f.setEps(0.1); 
  v = vf::project( _space=Vh, _range=elements(mesh), _expr=idf(_f));
/// [marker3]

  // export results
  e->add( "v", v );
  e->save();
}
