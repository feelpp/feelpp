/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   Date: 2011-16-12

   Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
   Copyright (C) CNRS

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
   \file parsedfunction.h
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2012-06-15
*/

#ifndef __PARSEDFUNCTION_GINAC_HPP
#define __PARSEDFUNCTION_GINAC_HPP 1

#include <ginac/ginac.h>

// #include <feel/options.hpp>
// #include <feel/feelcore/environment.hpp>
// #include <feel/feelcore/feel.hpp>
// #include <feel/feelcore/application.hpp>
// #include <feel/feelalg/backend.hpp>
// #include <feel/feeldiscr/functionspace.hpp>
// #include <feel/feeldiscr/mesh.hpp>
// #include <feel/feelmesh/filters.hpp>
// #include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
// #include <feel/feelfilters/geotool.hpp>

#include <boost/assert.hpp>
#include <boost/regex.hpp>

using namespace Feel;
using  GiNaC::symbol;
using  GiNaC::lst;
using  GiNaC::ex;
using  GiNaC::ex_to;

#include <string>
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
        for(auto& sym: syms)
            {
                if (sym.get_name() == s)
                    {
                        result = true;
                        found_symbol = sym;
                    }
            };
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

template<int Dim, typename T = double >
class parsedfunc_ginac
{
public :
    std::string function_str;
    std::string params_str;
    ex function_ex;
    std::vector<symbol> vars;
    std::string filename; // ginac file associated with the parsedfunc

    parsedfunc_ginac()
    {
        DVLOG(2) << "create empty function\n";
    };

    parsedfunc_ginac(const std::string & expression, std::string name="", std::vector<std::string> const& p = std::vector<std::string>())
    :
        function_str( expression ),
        vars( symbols<Dim>() ),
        filename( name )
    {
        if(!p.empty())
            {
                std::vector<symbol> params = Symbols(p);
                // Add parameters has symbols in vars
                vars.insert( vars.end(), params.begin(), params.end() );
            }

        LOG(INFO) << "Loading function : " << function_str << std::endl;
        for(symbol const& s: vars) {LOG(INFO) << "Symbol " << s.get_name() << " found\n";};
        function_ex = parse(function_str, vars);
        LOG(INFO) << "function parsed is : " << function_ex << "\n";
    };

    std::string func_str(){return function_str;}
    GiNaC::ex func_ex(){return function_ex;}
    std::vector<symbol> func_symbols(){return vars;}

    auto getExpr()->decltype( expr( function_str, vars, filename ) )
    {
        return expr( function_str, vars, filename );
    }

    auto getExpr(std::map<std::string, T> values)->decltype( expr( function_str, vars, filename ) )
    {
        //auto myexpr = expr( function_ex, vars, filename );
        auto myexpr = expr( function_str, vars, filename );
        if(!values.empty())
            myexpr.expression().setParameterValues( values );
        return myexpr;
    }

    auto getExpr(std::map<std::string, Feel::vf::Expr<T> > values)->decltype( expr( function_str, vars, filename ) )
    {
        auto myexpr = expr( function_str, vars, filename );
        if(!values.empty())
            {
                std::map<std::string, double> Tmap;
                for( typename std::map<std::string, Feel::vf::Expr<T> >::iterator it = values.begin(); it!=values.end(); it++)
                    Tmap.insert( std::pair<std::string, double>(it->first, (it->second).evaluate()) );

                myexpr.expression().setParameterValues( Tmap );
            }
        return myexpr;
    }

    bool isConstant()
    {
        MyVisitor v;
        function_ex.traverse(v);
        std::list<std::string> symbols = v.symbol_names();

#if 0
        std::cout << "Ginac expr: " << function_str << " : symbols [" << std::flush;
        for(std::string const& s: symbols) {std::cout << "Found Symbol :" << s << std::endl;};
        std::cout << "]\n" << std::flush;

        // ex_to<numeric> is an unsafe cast, so check the type first
        if (symbols.empty())
            {
                ex f = evalf(function_ex);
                if ( GiNaC::is_a<GiNaC::numeric>(f) )
                    {
                        double d = ex_to<GiNaC::numeric>(f).to_double();
                        cout << d << endl;
                    }
            }
#endif
        return symbols.empty();
    }

    double getConstant()
    {
        double val = 0; // should be init to NAN

        //assert(this->isConstant());
        if (this->isConstant())
            {
                ex f = evalf(function_ex); // ex_to<numeric> is an unsafe cast, so check the type first
                if ( GiNaC::is_a<GiNaC::numeric>(f) )
                    {
                        val = ex_to<GiNaC::numeric>(f).to_double();
                        //std::cout << val << std::endl;
                    }
            }
        else
            {
                std::cout << "Ginac expr: " << function_str << " is not Constant\n" << std::flush;
                LOG(INFO) << "function parsed is not Constant : " << function_ex << " \n";
            }

        return val;
    }
};

inline void remove_so_files()
{
    fs::directory_iterator end;
    for( fs::directory_iterator iter( fs::current_path() ) ; iter != end ; ++iter ) {
        if ( !is_directory( *iter ) )
            {
                std::string sofile_regx = ".+(.so)+";
                boost::regex e(sofile_regx);
                boost::smatch what;

                if(boost::regex_match(iter->path().string(), what, e, boost::match_extra))
                    std::remove( iter->path().string().c_str() );
            }
    }
}

#endif /* __PARSEDFUNCTION_GINAC_HPP */
