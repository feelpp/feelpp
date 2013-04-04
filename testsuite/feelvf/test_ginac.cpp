/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2008-02-07

   Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
/**
   \file test_ginac.cpp
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \date 2013-04-04
*/
#include <iostream>
#include <string>

#include <feel/feel.hpp>
#include <feel/feelvf/ginac.hpp>
#include <boost/range/algorithm/for_each.hpp>
using namespace Feel;

int main( int argc, char* argv[] )
{
    using GiNaC::table;
    using GiNaC::symbol;

    ex exact_parsed;
    std::vector<symbol> vars;
    std::string exact;
    int dim;

    std::cout << "dim=" << std::flush; std::cin >> dim; std::cout << std::flush;
    std::cout << "exact=" << std::flush; std::cin >> exact; std::cout << std::flush;
    switch (dim) {
    case(1) : {
        vars = symbols<1>();
        exact_parsed = parse(exact, vars);
        auto f = expr(exact,vars);
        break;
    }
    case(2) : {
        vars = symbols<2>();
        exact_parsed = parse(exact, vars); 
        boost::for_each( ex_to<table>(exaxt_parsed), [](std::pair<std::string, ex> const& s ) {LOG(INFO) << "Symbol " << s.first << " added\n";} );
        auto f = expr(exact,vars);
        break;
    }
    case(3) : {
        vars = symbols<3>();
        exact_parsed = parse(exact, vars);
        auto f = expr(exact,vars);
        break;
    }
    default: {
        std::cerr << "wrong dimension - should be lesser or egal to 3\n";
        return 1;
    }
    }

    std::cout << exact_parsed << std::endl;
    return 0;
}
