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
#include <list>

#include <feel/feel.hpp>
#include <feel/feelvf/ginac.hpp>
#include <boost/range/algorithm/for_each.hpp>
using namespace Feel;

class my_visitor :
    public GiNaC::visitor,
    public GiNaC::symbol::visitor
{
public:
    const std::list<std::string> & get_symbols()
    {
        symbols.sort();
        symbols.unique();
        return symbols;
    }

private:
    std::list<std::string> symbols;

    void visit(const symbol & s)
    { 
        symbols.push_back(s.get_name()); 
    }

};
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description ginacoptions("Ginac options");
    ginacoptions.add_options()
        ("exact", Feel::po::value<std::string>()->default_value( "" ), "name of the input")
        ("dim", Feel::po::value<int>()->default_value( 0 ), "geometric dimension")
        ;
    return ginacoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "ginac" ,
                           "ginac" ,
                           "0.1",
                           "test ginac integration with Feelpp",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2012-2013 Laboratoire national des Champs magnetiques Intenses");

    about.addAuthor("Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "");
    return about;

}

int main( int argc, char* argv[] )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    using GiNaC::symtab;
    using GiNaC::symbol;

    ex exact_parsed;
    std::vector<symbol> vars;
    std::string exact;
    int dim;

    std::cout << "strict_ginac_parser=" << strict_ginac_parser << std::endl << std::flush;
    std::cout << "dim=" << std::flush; std::cin >> dim; std::cout << std::flush;
    std::cout << "exact=" << std::flush; std::cin >> exact; std::cout << std::flush;
    switch (dim) {
    case(1) : {
        vars = symbols<1>();
        exact_parsed = parse(exact, vars);
        break;
    }
    case(2) : {
        vars = symbols<2>();
        exact_parsed = parse(exact, vars); 
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

    // retrieve symbols - not working because nops is not "recursive"
    std::cout << "Loading symbols from : " << exact_parsed << std::endl << std::flush;
    std::cout << "Contains " << exact_parsed.GiNaC::ex::nops() << " Ginac:ex objects\n" << std::flush;

    for (GiNaC::const_iterator i=exact_parsed.begin(); i!=exact_parsed.end(); ++i)
        {
            if ( GiNaC::is_a<symbol>(*i) )
                std::cout << "Found Symbol : " << GiNaC::ex_to<symbol>(*i).get_name() << std::endl;
        }
    std::cout << std::flush;
    //not working because ex is not mutable
    //boost::for_each(exact_parsed, [](GiNaC::ex const& e) {if (GiNaC::is_a<symbol>(e)) std::cout << "Found Symbol : " <<  GiNaC::ex_to<symbol>(e).get_name() << "\n";});

    // Retrieve each symbols using viitor
    std::cout << "Loading symbols from : " << exact_parsed << " (visitor)" << std::endl << std::flush;
    my_visitor v;
    exact_parsed.traverse(v);
    std::list<std::string> symbols = v.get_symbols();
    boost::for_each(symbols, [](std::string const& s) {std::cout << "Found Symbol :" << s << std::endl;});
    return 0;
}
