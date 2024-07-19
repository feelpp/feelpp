/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-10-24
*/
#include <boost/algorithm/string/predicate.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelvf/ginac.hpp>
#include <fmt/core.h>
#include <fmt/ostream.h>

namespace GiNaC
{
std::string str( ex && f )
{
    std::ostringstream ostr;
    ostr << f;
    return ostr.str();
}
std::string str( ex const& f )
{
    std::ostringstream ostr;
    ostr << f;
    return ostr.str();
}
std::string strsymbol( std::vector<symbol> const& f )
{
    std::ostringstream ostr;
    ostr << "(";
    for(int i =0; i < f.size();++i)
    {
        ostr << f[i].get_name();
        if ( i < f.size()-1 )
            ostr << ",";
    }
    ostr << ")";
    return ostr.str();
}

ex parse( std::string const& str, std::vector<symbol> const& syms, std::vector<symbol> const& params )
{
    using namespace Feel;
    VLOG(1) << "Parsing " << str << " using GiNaC with " << syms.size() << " symbols";
    VLOG(1) <<" . symbols : "  << strsymbol(syms);
    VLOG(1) <<" . parameters : "  << strsymbol(params);

    using GiNaC::symbol;
    using GiNaC::symtab;
    using GiNaC::parser;
    using GiNaC::parse_error;
    symtab table;
    VLOG(1) <<"Inserting symbols in symbol table";

    std::vector<symbol> total_syms;
    boost::for_each( syms, [&table, &total_syms]( symbol const& param )
                     {
                         total_syms.push_back(symbol(param));
                         table[param.get_name()] = param;
                     } );

    VLOG(1) <<"Inserting params and in symbol table: " << strsymbol(total_syms);

    boost::for_each( params, [&table, &total_syms]( symbol const& param )
                     {
                         total_syms.push_back(symbol(param));
                         table[param.get_name()] = param;
                     } );
    VLOG(1) << " . table : " << table;
    VLOG(1) <<"Defining parser";
    parser reader(table ,option(_name="ginac.strict-parser").as<bool>()); // true to ensure that no more symbols are added

    VLOG(1) <<"parse expression\n";
    ex e; // = reader(str);
    try
    {
        e = reader(str);
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
    VLOG(1) << "e=" << e << "\n";
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

    VLOG(1) << "Parsing " << str << " using GiNaC";
    std::vector<std::string> fields;
    boost::split( fields, str, boost::is_any_of(seps), boost::token_compress_on );
    int fsize = fields.size();
    CHECK( fsize  > 0 ) << "bad expression format";
    std::string strexpr( fields[0] );
    std::vector<std::string> strsyms;
#if 0
    if(fsize==1)
        strsyms.push_back("0"); // no symbols means constant expression
    else
        for( auto it=fields.begin()+1; it!=fields.end(); ++it )
            strsyms.push_back( *it );
#else
    for( auto it=fields.begin()+1; it!=fields.end(); ++it )
        strsyms.push_back( *it );
#endif
    std::vector<symbol> syms;
    std::for_each( strsyms.begin(), strsyms.end(),
                   [&syms] ( std::string const& sym ) { syms.push_back( symbol(sym) ); } );


    VLOG(1) << " . Number of symbols " << syms.size() << "\n";
    VLOG(1) << " . symbols : "  << strsymbol(syms);
    VLOG(1) << " . Number of params " << params.size() << "\n";
    VLOG(1) << " . symbols : "  << strsymbol(params);

    symtab table;
    VLOG(1) <<"Inserting symbols in symbol table";
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
                         table.insert({param.get_name(), param});
                     } );

    VLOG(1) <<"Inserting params and in symbol table : " << strsymbol(total_syms);

    boost::for_each( params, [&table, &total_syms]( symbol const& param )
                     {
                         total_syms.push_back(symbol(param));
                         table.insert({param.get_name(), param});
                     } );
#if 0
    for ( auto it=table.begin(),en=table.end() ; it!=en ; ++it )
        VLOG(1) <<" - table : "  << it->first << "\t" << it->second;
#else
    VLOG(1) << " . table : " << table;
#endif

    VLOG(1) <<"Defining parser";
    parser reader(table ,option(_name="ginac.strict-parser").as<bool>()); // true to ensure that no more symbols are added

    VLOG(1) <<"parse expression: " << strexpr;
    if ( boost::algorithm::contains( strexpr, "// Not supported in C" ) )
    {
        VLOG(1) <<"invalid code: " << table;
        throw std::invalid_argument( fmt::format( "invalid code: ", table ) );
    }
    ex e; // = reader(str);
    try
    {
        VLOG(1) << "parse expression 2: " << strexpr;
        e = reader(strexpr);
        VLOG(1) << "parse expression 3: " << strexpr;
    }
    catch (std::invalid_argument& err)
    {
        throw std::invalid_argument( fmt::format( "GiNaC error parsing {}: {}", e, err.what() ) );
        reader.strict = false;
        e =reader(strexpr);
        throw std::invalid_argument( fmt::format( "GiNaC error parsing {}: {}", e, err.what() ) );
    }
    catch ( ... )
    {
        throw;
    }

    VLOG(1) << "e=" << e << "\n";
    return {e,syms};
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
    return matrix{};
}

matrix
curl( matrix const& f, std::vector<symbol> const& l )
{
	VLOG(1) << "matrix version\n";
    CHECK(0) << "not implemented yet\n";
    return matrix{};
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
    VLOG(1) << "compute laplacian of " << s << std::endl;
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
    if ( is_a<lst>(f) )
    {
        std::vector<ex> g(f.nops());
        for( int i = 0; i < f.nops(); ++i )
            g[i] = f.op(i).diff( l, n );
        matrix ret(f.nops(),1,g);
        return ret;
    }
    else
    {
        matrix ret(1,1);
        ret(0,0)=f.diff( l,n );
        return ret;
    }
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
    return f.subs({s}, {GiNaC::numeric(val)});
}

ex substitute(ex const &f, symbol const& s, ex const& g )
{
    return f.subs({s}, {GiNaC::ex(g)});
}

matrix substitute(matrix const &f, symbol const& s, const double val )
{
    matrix ff(f.rows(), f.cols());

    for(int i=0; i<f.rows(); ++i)
    {
        for(int j=0; j<f.cols(); ++j)
            ff.set(i, j, f(i,j).subs({s}, {GiNaC::numeric(val)}) );
    }
    return ff;
}

matrix substitute(matrix const &f, symbol const& s, ex const& g )
{
    matrix ff(f.rows(), f.cols());

    for(int i=0; i<f.rows(); ++i)
    {
        for(int j=0; j<f.cols(); ++j)
            ff.set(i, j, f(i,j).subs({s}, {GiNaC::ex(g)}) );
    }
    return ff;
}

#if 0
int totalDegree( GiNaC::ex const& expr, std::vector<std::pair<GiNaC::symbol,Feel::uint16_type>> const& symbolsDegree )
{
    int res = 0;
    for ( int sd=0;sd<symbolsDegree.size();++sd )
    {
        GiNaC::symbol symbDegree = symbolsDegree[sd].first;
        int unitDegree = symbolsDegree[sd].second;
        int mydegree = expr.degree( symbDegree )*unitDegree;
        std::vector<std::pair<GiNaC::symbol,Feel::uint16_type>> newSymbolsDegree;
        for (int k=0;k<symbolsDegree.size();++k)
        {
            if ( sd == k ) continue;
            newSymbolsDegree.push_back( symbolsDegree[k] );
        }
        int coeffDegree = totalDegree( expr.lcoeff( symbDegree ), newSymbolsDegree );
        res = std::max( res, mydegree + coeffDegree );
    }
    return res;
}

#else

// update total degree from a symbol (the first) and all subsymbols dependencies
// implementation based on a heap .
void totalDegreeImpl( std::vector<std::tuple<std::reference_wrapper<const GiNaC::symbol>,std::optional<GiNaC::ex>,std::vector<int>,Feel::uint16_type> > & symbolsDegree )
{
    using namespace Feel;
    GiNaC::symbol symbDegreeComputation;
    std::vector<int> heapSymbolIds = {0};
    int currentSymbolId = heapSymbolIds.back();
    while ( !heapSymbolIds.empty() )
    {
        currentSymbolId = heapSymbolIds.back();
        auto & thesymbolData = symbolsDegree.at( currentSymbolId );
        if ( std::get<3>(thesymbolData) != invalid_v<uint16_type> )
        {
            heapSymbolIds.pop_back();
            continue;
        }

        auto const& symbolDeps = std::get<2>( thesymbolData );
        std::vector<int> needSymbolIds;
        bool canBeComputed = true;
        for ( int k : symbolDeps )
        {
            auto & thesubsymbolData = symbolsDegree.at( k );
            if ( std::get<3>( thesubsymbolData ) == invalid_v<uint16_type> )
            {
                needSymbolIds.push_back( k );
                canBeComputed = false;
            }
        }

        if ( canBeComputed )
        {
            auto symbolEx = std::get<1>( thesymbolData );
            for ( int k : symbolDeps )
            {
                auto const& thesubsymbolData = symbolsDegree.at( k );
                auto const& subsymbol = std::get<0>( thesubsymbolData ).get();
                uint16_type symbolDegree = std::get<3>( thesubsymbolData );
                // susbtitute subsymbol wtih symbDegreeComputation
                *symbolEx = symbolEx->subs( subsymbol == pow(symbDegreeComputation,symbolDegree) );
            }
            // compute degree with symbDegreeComputation
            int mydegree = symbolEx->degree( symbDegreeComputation );
            std::get<3>( thesymbolData ) = mydegree;
            heapSymbolIds.pop_back();
        }
        else
        {
            for ( int k : needSymbolIds )
                heapSymbolIds.push_back( k );
        }
    }

}


int totalDegree( GiNaC::ex const& expr, std::vector<std::pair<GiNaC::symbol,Feel::uint16_type>> const& symbolsDegreeIn )
{
    using namespace Feel;

    GiNaC::symbol symbCurentExpr;
    std::vector<std::tuple<std::reference_wrapper<const GiNaC::symbol>,std::optional<GiNaC::ex>,std::vector<int>,uint16_type> > symbolsDegree;
    std::vector<int> subsymbolIds;
    for (int k=0;k<symbolsDegreeIn.size();++k)
        subsymbolIds.push_back( k+1 );

    symbolsDegree.push_back( std::make_tuple( std::cref(symbCurentExpr),
                                              std::make_optional<GiNaC::ex>( expr ),
                                              subsymbolIds, invalid_v<uint16_type> ) );
    for ( auto const& [symb,deg] : symbolsDegreeIn )
        symbolsDegree.push_back( std::make_tuple( std::cref(symb),
                                                  std::optional<GiNaC::ex>{},
                                                  std::vector<int>{}, deg ) );

    totalDegreeImpl(symbolsDegree);

    return std::get<3>( symbolsDegree.front() );
}

#endif

std::vector<std::pair< bool,double> >
toNumericValues( ex const& expr )
{
    std::vector<std::pair< bool,double> > res;
    if ( is_a<matrix>(expr) )
    {
        matrix m( ex_to<matrix>(expr) );
        int ni = m.rows();
        int nj = m.cols();
        res.resize( ni*nj, std::make_pair( false, double(0.) ) );
        for( int i = 0; i < ni; ++i )
        {
            for( int j = 0; j < nj; ++j )
            {
                GiNaC::ex funEvalf = m(i,j).evalf();
                if ( GiNaC::is_a<GiNaC::numeric>(funEvalf) )
                {
                    GiNaC::numeric funNumeric = GiNaC::ex_to<GiNaC::numeric>(funEvalf);
                    if ( funNumeric.is_real() )
                        res[i*nj+j] = std::make_pair( true, funNumeric.to_double() );
                }
            }
        }
    }
    else if ( is_a<lst>(expr) )
    {
        int no = expr.nops();
        res.resize( no, std::make_pair( false, double(0.) ) );
        for( int i = 0; i < no; ++i )
        {
            GiNaC::ex funEvalf = expr.op(i).evalf();
            if ( GiNaC::is_a<GiNaC::numeric>(funEvalf) )
            {
                GiNaC::numeric funNumeric = GiNaC::ex_to<GiNaC::numeric>(funEvalf);
                if ( funNumeric.is_real() )
                    res[i] = std::make_pair( true, funNumeric.to_double() );
            }
        }
    }
    else
    {
        res.resize( 1, std::make_pair( false, double(0.) ) );
        GiNaC::ex funEvalf = expr.evalf();
        if ( GiNaC::is_a<GiNaC::numeric>(funEvalf) )
        {
            GiNaC::numeric funNumeric = GiNaC::ex_to<GiNaC::numeric>(funEvalf);
            if ( funNumeric.is_real() )
                res[0] = std::make_pair( true, funNumeric.to_double() );
        }
    }
    return res;
}


}
