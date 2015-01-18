
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

std::pair< ex, std::vector<symbol> >
parse( std::string const& str, std::string const& seps, std::vector<symbol> const& params )
{
    using namespace Feel;
    using GiNaC::symbol;
    using GiNaC::lst;
    using GiNaC::ex;
    using GiNaC::parser;
    using GiNaC::symtab;
    using GiNaC::parse_error;

    LOG(INFO) << "Parsing " << str << " using GiNaC";
    std::vector<std::string> fields;
    boost::split( fields, str, boost::is_any_of(seps), boost::token_compress_on );
    int fsize = fields.size();
    CHECK( fsize  > 0 ) << "bad expression format";
    std::string strexpr( fields[0] );
    std::vector<std::string> strsyms;
    if(fsize==1)
        strsyms.push_back("0"); // no symbols means constant expression
    else
        for( auto it=fields.begin()+1; it!=fields.end(); ++it )
            strsyms.push_back( *it );
    std::vector<symbol> syms;
    std::for_each( strsyms.begin(), strsyms.end(),
                   [&syms] ( std::string const& sym ) { syms.push_back( symbol(sym) ); } );


    LOG(INFO) << "Number of symbols " << syms.size() << "\n";
    for(int i =0; i < syms.size();++i)
        LOG(INFO) <<" - symbol : "  << syms[i].get_name();

    LOG(INFO) << "Number of params " << params.size() << "\n";
    for(int i =0; i < params.size();++i)
        LOG(INFO) <<" - param : "  << params[i].get_name();

    symtab table;
    LOG(INFO) <<"Inserting symbols in symbol table";
#if 0
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
#endif
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

    for ( auto it=table.begin(),en=table.end() ; it!=en ; ++it )
        LOG(INFO) <<" - table : "  << it->first << "\t" << it->second;


    LOG(INFO) <<"Defining parser";
    parser reader(table ,option(_name="ginac.strict-parser").as<bool>()); // true to ensure that no more symbols are added

    LOG(INFO) <<"parse expression: " << strexpr;
    ex e; // = reader(str);
    try
    {
        e = reader(strexpr);
    }
    catch (std::invalid_argument& err)
    {
        reader.strict = false;
        e =reader(strexpr);

        std::cerr << "GiNaC error parsing " << e << " : " << err.what() << std::endl;
        exit(1);
    }
    catch ( ... )
    {
        std::cerr << "Exception of unknown type!\n";
    }

    LOG(INFO) << "e=" << e << "\n";
    return std::make_pair(e,syms);
}

matrix
grad( ex const& f, std::vector<symbol> const& l )
{
    lst g;
    std::for_each( l.begin(), l.end(),
                   [&] ( symbol const& x ) { g.append( f.diff( x ) ); } );
    if ( f.info(info_flags::list) ) // is_a<lst>( f )
    {
        std::vector<ex> v;
        for( int e = 0; e < g.nops(); ++e )
        {
            if ( is_a<lst>( g.op(e) ) )
            {
                // be careful, we need to transpose because g is actually the
                // transpose of the gradient
                for( int i = 0; i < g.op(e).nops(); ++i )
                    v.push_back( g.op(i).op(e) );
            }
            else
                v.push_back( g.op(e) );
        }

        matrix h( g.nops(), g.op(0).nops(), v );
        return h;
    }
    else
    {
        std::vector<ex> v( g.nops() );
        matrix h( 1, g.nops(), g );
        return h;

    }
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
div( ex const& f, std::vector<symbol> const& l )
{
    // matricial fields are not yet supported
    CHECK( is_a<lst>( f ) ) << "Invalid expression " << f << " : cannot compute its divergence\n";
    lst g;
    std::for_each( l.begin(), l.end(),
                   [&] ( symbol const& x ) { g.append( f.diff( x ) ); } );
    std::vector<ex> v(1);
    for( int e = 0; e < g.nops(); ++e )
    {
        v[0] += g.op(e).op(e);
    }
    matrix h( 1, 1, v );
    return h;
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

matrix
curl( ex const& f, std::vector<symbol> const& l )
{
	CHECK(f.nops() > 1) << "Invalid expression for curl operator: the expression has to be a vector";
	CHECK(f.nops() <=3) << "Invalid expression for curl operator: the expression has to be a vector with max 3 components";
  if   ( is_a<lst>( f ) ){
    if   ( f.nops() == 2 ){
        std::vector<ex> v(1);
        CHECK( l[0].get_name() == "x") << "Symbol x not present in list of symbols, cannot compute curl(" << f << ")\n";
        CHECK( l[1].get_name() == "y") << "Symbol y not present in list of symbols, cannot compute curl(" << f << ")\n";
        v[0] = f.op(1).diff(l[0])-f.op(0).diff(l[1]);
        matrix h( 1, 1, v );
        return h;
    }
    if   ( f.nops() == 3){
        CHECK( l[0].get_name() == "x") << "Symbol x not present in list of symbols, cannot compute curl(" << f << ")\n";
        CHECK( l[1].get_name() == "y") << "Symbol y not present in list of symbols, cannot compute curl(" << f << ")\n";
        CHECK( l[2].get_name() == "z") << "Symbol z not present in list of symbols, cannot compute curl(" << f << ")\n";
        std::vector<ex> v(3);
        v[0]=f.op(2).diff(l[1])-f.op(1).diff(l[2]);
        v[1]=f.op(0).diff(l[2])-f.op(2).diff(l[0]);
        v[2]=f.op(1).diff(l[0])-f.op(0).diff(l[1]);
        matrix h( 3, 1, v );
        return h;
    }
	}else{ //   ( is_a<lst>( f ) )
		if   ( f.nops() == 2 ){
        CHECK( l[0].get_name() == "x") << "Symbol x not present in list of symbols, cannot compute curl(" << f << ")\n";
        CHECK( l[1].get_name() == "y") << "Symbol y not present in list of symbols, cannot compute curl(" << f << ")\n";
        std::vector<ex> v(1);
        v[0] = f[1].diff(l[0])-f[0].diff(l[1]);
        matrix h( 1, 1, v );
        return h;
		}
    if   ( f.nops() == 3 ){
			CHECK( l[0].get_name() == "x") << "Symbol x not present in list of symbols, cannot compute curl(" << f << ")\n";
      CHECK( l[1].get_name() == "y") << "Symbol y not present in list of symbols, cannot compute curl(" << f << ")\n";
      CHECK( l[2].get_name() == "z") << "Symbol z not present in list of symbols, cannot compute curl(" << f << ")\n";
      std::vector<ex> v(3);
      v[0]=f[2].diff(l[1])-f[1].diff(l[2]);
      v[1]=f[0].diff(l[2])-f[2].diff(l[0]);
      v[2]=f[1].diff(l[0])-f[0].diff(l[1]);
      matrix h( 3, 1, v );
      return h;
		}
	}
	CHECK(0) << "Invalid expression " << f << " cannot compute its curl\n";
}

matrix
curl( matrix const& f, std::vector<symbol> const& l )
{
	LOG(INFO) << "matrix version\n";
    CHECK(0) << "not implemented yet\n";
}

matrix
laplacian( ex const& f, std::vector<symbol> const& l )
{
    ex e = f.evalm();
    if ( e.is_a_matrix() ) //is_a<matrix>(e) )
    {
        matrix m( ex_to<matrix>(e) );
        matrix g( m.rows(),1 );
        for( int i = 0; i < m.rows(); ++i )
        {
            std::for_each( l.begin(), l.end(),
                           [&] ( symbol const& x )
                           {
                               g(i,0) += m(i,0).diff( x,2 );
                           } );
        }
        return g;
    }
    else if ( e.info(info_flags::list) ) //is_a<lst>(e) )
    {
        std::vector<ex> g(e.nops());
        for(int n = 0; n < e.nops(); ++n )
        {
            std::for_each( l.begin(), l.end(),
                           [&] ( symbol const& x )
                           {
                               CHECK ( !e.op(0).info(info_flags::list) ) //  !is_a<lst>(e.op(0)) )
                                   << "the matricial case (expression: " << e << ") is not implemented, please contact feelpp-devel@feelpp.org";
                               g[n] += e.op(n).diff( x,2 );
                           } );
        }
        matrix h(e.nops(),1,g);
        return h;
    }
    else
    {
        ex g;
        std::for_each( l.begin(), l.end(),
                       [&] ( symbol const& x )
                       {
                           g += e.diff( x,2 );
                       } );
        matrix h(1,1,std::vector<ex>(1,g));
        return h;

    }
}
matrix
laplacian( std::string const& s, std::vector<symbol> const& l )
{
    LOG(INFO) << "compute laplacian of " << s << std::endl;
    return laplacian( parse( s, l ), l );
}
matrix
laplacian( matrix const& f, std::vector<symbol> const& l )
{
    std::cout <<  "use matrix lap\n";
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
ex
laplacian( std::string const& s, std::vector<symbol> const& l, std::vector<symbol> const& p )
{
    return laplacian( parse( s, l, p ), l )(0,0);
}
matrix diff(ex const& f, symbol const& l, const int n)
{
    matrix ret(1,1);
    ret(0,0)=f.diff( l,n );
    return ret;
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
