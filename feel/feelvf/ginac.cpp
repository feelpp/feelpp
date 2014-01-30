/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
   Date: 2012-10-24

   Copyright (C) 2012 Université de Strasbourg

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
/**
   \file ginac.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-10-24
*/
#include <feel/feelcore/environment.hpp>

#include <ginac/ginac.h>
#include <boost/foreach.hpp>
#include <boost/range/algorithm/for_each.hpp>

namespace GiNaC
{
ex parse( std::string const& str, std::vector<symbol> const& syms, std::vector<symbol> const& params = std::vector<symbol>() )
{
    using namespace Feel;
    LOG(INFO) << "Parsing " << str << " using GiNaC";

    LOG(INFO) << "Number of symbols " << syms.size() << "\n";

    for(int i =0; i < syms.size();++i)
        LOG(INFO) <<" - symbol : "  << syms[i].get_name();

    LOG(INFO) << "Number of params " << params.size() << "\n";

    for(int i =0; i < params.size();++i)
        LOG(INFO) <<" - param : "  << params[i].get_name();

    using GiNaC::symbol;
    using GiNaC::symtab;
    using GiNaC::parser;
    using GiNaC::parse_error;
    symtab table;
    LOG(INFO) <<"Inserting symbols in symbol table";

    table["x"]=syms[0];
    if ( syms.size() == 2 )
    {
        table["y"]=syms[1];
    }
    if ( syms.size() == 3 )
    {
        table["y"]=syms[1];
        table["z"]=syms[2];
    }
    std::vector<symbol> total_syms;
    boost::for_each( syms, [&table, &total_syms]( symbol const& param )
                     {
                         total_syms.push_back(symbol(param));
                         LOG(INFO) << "adding param: " << param << std::endl;
                         table[param.get_name()] = param;
                     } );

    LOG(INFO) <<"Inserting params and in symbol table";

    boost::for_each( params, [&table, &total_syms]( symbol const& param )
                     {
                         total_syms.push_back(symbol(param));
                         LOG(INFO) << "adding param: " << param << std::endl;
                         table[param.get_name()] = param;
                     } );

    LOG(INFO) <<"Defining parser";
    parser reader(table ,option(_name="ginac.strict-parser").as<bool>()); // true to ensure that no more symbols are added

    LOG(INFO) <<"parse expression\n";
    ex e; // = reader(str);
    try
    {
        e = reader(str);
#if 0
        if (!reader.strict)
        {
            symtab table_symbols = reader.get_syms();
            boost::for_each( table_symbols, [](std::pair<std::string, ex> const& s )
                             {
                                 LOG(INFO) << "Symbol " << s.first << " added\n";
                             }
                );
        }
#endif
    }
    catch (std::invalid_argument& err)
    {
        reader.strict = false;
        e =reader(str);

        std::cerr << "GiNaC error parsing " << e << " : " << err.what() << std::endl;
        exit(1);
    }
    catch ( ... )
    {
        std::cerr << "Exception of unknown type!\n";
    }

    LOG(INFO) << "parsed expression :" << e << "\n";
    return e;
}

matrix
grad( ex const& f, std::vector<symbol> const& l )
{
    //std::cout << "Dim=" << Dim << "\n";
    lst g;
    std::for_each( l.begin(), l.end(), [&] ( symbol const& x ) { g.append( f.diff( x ) ); } );
    return matrix( 1, l.size(), g );
}
matrix
grad( std::string const& s, std::vector<symbol> const& l )
{
    return grad( parse( s, l ),  l );
}

matrix
grad( matrix const& f, std::vector<symbol> const& l )
{
    matrix ff = f;
    if ( f.rows() == 1 )
        ff = f.transpose();

    matrix g( ff.rows(), l.size() );
    for( int i = 0; i < ff.rows(); ++i )
    {
        for( int j = 0; j < l.size(); ++j )
        {
            g.set( i, j, ff(i,0).diff( l[j] ) );
        }
    }
    return g;
}

matrix
div( matrix const& f, std::vector<symbol> const& l )
{
    matrix ff = f;
    if ( f.cols() == 1 )
        ff = f.transpose();

    matrix g(ff.rows(), 1 );
    for( int i = 0; i < ff.rows(); i++ )
    {
        ex e;
        for( int j = 0; j < ff.cols(); j++ )
            e += ff(i,j).diff( l[j] );
        g(i,0) = e;
    }
    return g;
}

ex
laplacian( ex const& f, std::vector<symbol> const& l )
{
    ex g;
    std::for_each( l.begin(), l.end(),
                   [&] ( symbol const& x )
                   {
                       g += f.diff( x,2 );
                   } );
    return g;
}
ex
laplacian( std::string const& s, std::vector<symbol> const& l )
{
    return laplacian( parse( s, l ), l );
}
matrix
laplacian( matrix const& f, std::vector<symbol> const& l )
{
    matrix g(f.rows(),1);
    for(int i = 0; i < f.rows(); ++i )
    {
        ex lap;
        std::for_each( l.begin(), l.end(),
                       [&] ( symbol const& x )
                       {
                           lap += f(i,0).diff( x,2 );
                       } );
        g(i,0) = lap;
    }
    return g;
}

ex diff(ex const& f, symbol const& l, const int n)
{
    return f.diff( l,n );
}


matrix diff(matrix const& f, symbol const& l, const int n)
{
    matrix g(f.rows(),1);
    for(int i = 0; i < f.rows(); ++i )
    {
        g(i,0) = f(i,0).diff( symbol(l),n );
    }
    return g;
}

ex substitute(ex const &f, symbol const& s, const double val )
{
    return f.subs(GiNaC::lst(s), GiNaC::lst(GiNaC::numeric(val)));
}

ex substitute(ex const &f, symbol const& s, ex const& g )
{
    return f.subs(GiNaC::lst(s), GiNaC::lst(GiNaC::ex(g)));
}

matrix substitute(matrix const &f, symbol const& s, const double val )
{
    matrix ff(f.rows(), f.cols());

    for(int i=0; i<f.rows(); ++i)
    {
        for(int j=0; j<f.cols(); ++j)
            ff.set(i, j, f(i,j).subs(GiNaC::lst(s), GiNaC::lst(GiNaC::numeric(val))) );
    }
    return ff;
}

matrix substitute(matrix const &f, symbol const& s, ex const& g )
{
    matrix ff(f.rows(), f.cols());

    for(int i=0; i<f.rows(); ++i)
    {
        for(int j=0; j<f.cols(); ++j)
            ff.set(i, j, f(i,j).subs(GiNaC::lst(s), GiNaC::lst(GiNaC::ex(g))) );
    }
    return ff;
}

}

