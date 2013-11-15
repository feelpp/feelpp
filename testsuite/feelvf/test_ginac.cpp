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

#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>

using namespace Feel;

class MyVisitor :
    public GiNaC::visitor,
    public GiNaC::symbol::visitor
{
public:
    const std::list<symbol> & symbols()
    {
        return syms;
    }

    const std::list<std::string> & symbol_names()
    {
        sym_names.sort();
        sym_names.unique();
        return sym_names;
    }

    bool hassymbol(std::string const& s, symbol &found_symbol)
    {
        bool result = false;
        boost::for_each(syms, [&s, &result, &found_symbol](symbol const& sym)
                        {
                            if (sym.get_name() == s)
                                {
                                    result = true;
                                    found_symbol = sym;
                                }
                        });
        return result;
    }

private:
    std::list<symbol> syms;
    std::list<std::string> sym_names;

    void visit(const symbol & s)
    {
        syms.push_back(s);
        sym_names.push_back(s.get_name());
    }

};
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description ginacoptions("Ginac options");
    ginacoptions.add_options()
        ("dim", Feel::po::value<int>()->default_value( 1 ), "geometric dimension")
        ("params", Feel::po::value<std::string>()->default_value( "a;b" ), "name of parameters")
        ("exact", Feel::po::value<std::string>()->default_value( "a*x+b" ), "name of the input")
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
    std::vector<symbol> vars, parameters;
    std::vector<std::string> lst_params;
    std::string exact;
    std::string params;
    int dim = 0;

    if ( option(_name="ginac.strict-parser").as<bool>() )
        std::cout << "Strict Ginac Parser enabled\n";
    if ( !env.vm()["dim"].as<int>() )
        {
            std::cout << "dim=" << std::flush;
            std::cin >> dim;
            std::cout << std::flush;
        }
        else
        dim = env.vm()["dim"].as<int>(); 

    // Load param list
    if ( !env.vm()["params"].as<std::string>().empty() )
        {
            boost::split(lst_params, params, boost::is_any_of(";"));
            parameters = symbols(lst_params);
        }

    // load exact
    if ( env.vm()["exact"].as<std::string>().empty() )
        {
            // Load param list
            std::cout << "params (eg params=\"a;b\")=" << std::flush; std::cin >> params;  std::cout << std::flush;
            if ( !params.empty() )
                {
                    boost::split(lst_params, params, boost::is_any_of(";"));
                    parameters = symbols(lst_params);
                }
            std::cout << parameters.size() << " params read\n";

            std::cout << "exact=" << std::flush;
            std::cin >> exact;
            std::cout << std::flush;
        }
    else
      exact = env.vm()["exact"].as<std::string>();

    switch (dim) {
    case(1) : {
        vars = symbols<1>();
        break;
    }
    case(2) : {
        vars = symbols<2>();
        break;
    }
    case(3) : {
        vars = symbols<3>();
        break;
    }
    default: {
        std::cerr << "wrong dimension - should be lesser or egal to 3 - is equatl to " << dim << "\n";
        return 1;
    }
    }

#if 0
    if ( !parameters.size() )
        exact_parsed = parse(exact, vars);
    else
#endif
        exact_parsed = parse(exact, vars, parameters);

#if 0
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
    //boost::for_each(exact_parsed, [](GiNaC::ex const& e) 
    //{
    //if (GiNaC::is_a<symbol>(e)) std::cout << "Found Symbol : " <<  GiNaC::ex_to<symbol>(e).get_name() << "\n";
    //});
#endif

    // Retrieve each symbols using visitor
    std::cout << "Loading symbols from : " << exact_parsed << " (visitor)" << std::endl << std::flush;
    MyVisitor v;
    exact_parsed.traverse(v);
    std::list<std::string> symbols = v.symbol_names();
    boost::for_each(symbols, [](std::string const& s) {std::cout << "Found Symbol :" << s << std::endl;});

    // Substitute expressions
    std::cout << "Replace vars\n";
    boost::for_each(vars, [&exact_parsed](symbol const& sym)
                    {
                        ex f = exact_parsed;
                        std::cout << "Replace " << sym.get_name() << " in " << f << " :";
                        std::cout << substitute(f, sym, 1.0) << std::endl;
                    });

    std::cout << "Replace params\n";
    boost::for_each(parameters, [&exact_parsed](symbol const& sym)
                    {
                        ex f = exact_parsed;
                        if (f.has(sym))
                            {
                                std::cout << "Replace " << sym.get_name() << " in " << f << " :";
                                std::cout << substitute(f, sym, 1.0) << std::endl;
                            }
                    });

    std::cout << "Replace params using strings\n";
    boost::for_each(lst_params, [&exact_parsed](std::string const& s)
                    {
                        ex f = exact_parsed;
                        symbol sym;
                        //use MyVisitor class to find if symbol is present in f
                        MyVisitor f_v;
                        f.traverse(f_v);
                        if ( f_v.hassymbol(s, sym) )
                            {
                                std::cout << "Replace " << s << " in " << f << " :";
                                std::cout << substitute(f, sym, 1.0) << std::endl;
                            }
                    });

    std::cout << "Replace params using strings by function g()\n";
    boost::for_each(lst_params, [&exact_parsed, &exact, &vars, &parameters](std::string const& s)
                    {
                        ex f = exact_parsed;
                        ex f1;
                        symbol sym;
                        //use MyVisitor class to find if symbol is present in f
                        MyVisitor f_v;
                        f.traverse(f_v);
                        if ( f_v.hassymbol(s, sym) )
                            {
                                std::cout << "Replace " << s << " in " << f << " by :";
                                std::cin >> exact;
                                ex f1 = parse(exact, vars, parameters);
                                std::cout << std::flush;
                                std::cout << substitute(f, sym, f1) << std::endl;
                            }
                    });
		    
    exact="sin(x)";
    auto g = expr(exact, vars);
    		    
    return 0;
}
