/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2012-10-15

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
   \file ginac.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-10-15
*/
#ifndef FEELPP_GINAC_HPP
#define FEELPP_GINAC_HPP 1

#include <feel/feelcore/feelmacros.hpp>

// #if defined(__clang__)
// #if !defined(__APPLE__) && FEELPP_CLANG_AT_LEAST(3,9)
// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wundefined-var-template"
// #endif
// #endif

#include <ginac/ginac.h>
extern template GiNaC::registered_class_info GiNaC::container<std::list>::reg_info;
extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;

// #if defined(__clang__)
// #if !defined(__APPLE__) && FEELPP_CLANG_AT_LEAST(3,9)
// #pragma clang diagnostic pop
// #endif
// #endif

#include <boost/fusion/container/vector.hpp>

#include <boost/parameter/preprocessor.hpp>

#include <boost/foreach.hpp>
#include <boost/range/algorithm/for_each.hpp>

#include <feel/feelvf/expr.hpp>

namespace GiNaC
{
std::string str( ex && f );
std::string str( ex const& f );
matrix grad( ex const& f, std::vector<symbol> const& l );
matrix laplacian( ex const& f, std::vector<symbol> const& l );
matrix grad( std::string const& s, std::vector<symbol> const& l );
matrix laplacian( std::string const& s, std::vector<symbol> const& l );

matrix grad( matrix const& f, std::vector<symbol> const& l );
matrix div( ex const& f, std::vector<symbol> const& l );
matrix div( matrix const& f, std::vector<symbol> const& l );
matrix curl( ex const& f, std::vector<symbol> const& l );
matrix curl( matrix const& f, std::vector<symbol> const& l );
matrix laplacian( matrix const& f, std::vector<symbol> const& l );
ex laplacian( std::string const& s, std::vector<symbol> const& l, std::vector<symbol> const& p );

matrix diff(ex const& f, symbol const& l, const int n);
matrix diff(matrix const& f, symbol const& l, const int n);

ex substitute(ex const& f, symbol const& l, const double val );
ex substitute(ex const& f, symbol const& l, ex const & g );

matrix substitute(matrix const& f, symbol const& l, const double val );
matrix substitute(matrix const& f, symbol const& l, ex const & g );

//ex parse( std::string const& str, std::vector<symbol> const& syms );
ex parse( std::string const& str, std::vector<symbol> const& syms, std::vector<symbol> const& params = std::vector<symbol>());

} // GiNaC

namespace Feel
{
using GiNaC::matrix;
using GiNaC::symbol;
using GiNaC::lst;
using GiNaC::ex;
using GiNaC::parser;
using GiNaC::diff;
using GiNaC::laplacian;
using GiNaC::grad;
using GiNaC::div;
using GiNaC::parse;

template<int Dim> inline std::vector<symbol> symbols() { return {symbol("x")}; }
template<> inline std::vector<symbol> symbols<1>() { return {symbol("x")}; }
template<> inline std::vector<symbol> symbols<2>() { return {symbol("x"),symbol("y") };}
template<> inline std::vector<symbol> symbols<3>() { return {symbol("x"),symbol("y"),symbol("z") };}

inline
std::vector<symbol>
symbols( std::initializer_list<std::string> l )
{
    std::vector<symbol> s;
    std::for_each( l.begin(), l.end(), [&s] ( std::string const& sym ) { s.push_back( symbol(sym) ); } );
    return s;
}

inline
std::vector<symbol>
symbols( std::vector<std::string> l )
{
    std::vector<symbol> s;
    std::for_each( l.begin(), l.end(), [&s] ( std::string const& sym ) { s.push_back( symbol(sym) ); } );
    return s;
}

class Symbols : public std::vector<symbol>
{
public:
    Symbols():std::vector<symbol>(symbols({"x","y","z", "t"})) {}
    Symbols(std::initializer_list<std::string> s ):std::vector<symbol>(symbols(s)) {}
    Symbols(std::vector<std::string> const& s ):std::vector<symbol>(symbols(s)) {}
};

template<typename... Args>
class Fields
    :
    public boost::fusion::vector<Args...>
{
public:
    typedef boost::fusion::vector<Args...> super;
    typedef Fields<Args...> this_type;
	static const int s = sizeof...(Args);
    Fields( super const& m) : super( m ) {}

};

} // Feel namespace



// Feel::vf
namespace Feel
{
using GiNaC::matrix;
using GiNaC::symbol;
using GiNaC::lst;
using GiNaC::ex;
using GiNaC::parser;
}

namespace GiNaC
{
/**
 * \brief Parse a string expression
 *
 * \param str the string to parse
 * \param seps symbols separator
 * \param params parameters
 *
 * ### Format
 * The string format is: "GiNaC::ex,GiNaC::symbol,GiNaC::symbol,..."
 *
 * ### example :
 * ```auto a = parse("sqrt(x*y):x:y")```
 *
 * \return a pair containing the GiNaC expression, and a vector of GiNaC symbols.
 */
std::pair< ex, std::vector<symbol> >
parse( std::string const& str, std::string const& seps=":", std::vector<symbol> const& params = std::vector<symbol>());

} // GiNaC namespace

#include <feel/feelvf/ginacbase.hpp>
#include <feel/feelvf/detail/ginacbuildlibrary.hpp>
#include <feel/feelvf/detail/ginacex.hpp>
#include <feel/feelvf/detail/ginacexvf.hpp>
#include <feel/feelvf/detail/ginacmatrix.hpp>

namespace Feel
{
namespace vf
{

// return the number of components
inline
int nbComp( std::string expression )
{
    auto parseExpr = GiNaC::parse(expression);

    auto ginacEvalm = parseExpr.first.evalm();
    bool isLst = GiNaC::is_a<GiNaC::lst>(  ginacEvalm );
    int nComp = 1;
    if ( isLst )
        nComp = ginacEvalm.nops();
    return nComp;
}

inline
Expr< GinacEx<2> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    VLOG(2)<< "expr(GiNaC::ex)\n";
    return Expr< GinacEx<2> >(  GinacEx<2>( f, lsym, std::string(""), filename, world, dirLibExpr ) );
}

inline
Expr< GinacEx<2> >
expr( std::string const& s, std::vector<GiNaC::symbol> const& lsym, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return Expr< GinacEx<2> >(  GinacEx<2>( parse(s,lsym), lsym, s, filename, world, dirLibExpr ) );
}


/**
 * \brief functor enabling ginac
 *
 */
template<int Order>
inline
Expr< GinacEx<Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return Expr< GinacEx<Order> >(  GinacEx<Order>( f, lsym, std::string(""), filename, world, dirLibExpr ));
}

template<int Order>
inline
Expr< GinacEx<Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string const& exprDesc="", std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return Expr< GinacEx<Order> >(  GinacEx<Order>( f, lsym, exprDesc, filename, world, dirLibExpr ));
}

template<int Order>
inline
Expr< GinacEx<Order> >
expr( std::string const& s, std::vector<GiNaC::symbol> const& lsym, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return Expr< GinacEx<Order> >(  GinacEx<Order>( parse(s,lsym), lsym, s, filename, world, dirLibExpr ) );
}

/**
 * @brief Create an Feel++ expression from a GiNaC expression as a string
 *
 * @tparam Order     Expression order
 * @param s          String containing the ginac expression and symbols
 * @param filename   Shared file
 *
 * @return Feel++ Expression
 */
template<int Order>
inline
Expr< GinacEx<Order> >
expr( std::string const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacEx<Order> >(  GinacEx<Order>( g.first, g.second, s, filename, world, dirLibExpr ) );
}

/**
* @brief Create an Feel++ expression from a GiNaC expression as a string
*
* @param s          String containing the ginac expression and symbols
* @param filename   Shared file
*
* @return Feel++ Expression
*/
inline
Expr< GinacEx<2> > expr( std::string const& s, std::string const& filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return expr<2>( s, filename, world, dirLibExpr );
}

/**
* @brief Create an Feel++ expression from a GiNaC expression as a string
*
* @param s          String containing the ginac expression and symbols
* @param filename   Shared file
*
* @return Feel++ Expression
*/
inline
Expr< GinacEx<2> > expr( const char* s_, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string s = s_;
    return expr( s, filename, world, dirLibExpr );
}


/**
* @brief Create an Feel++ expression from a GiNaC expression as a string
*
* @param s          String containing the ginac expression and symbols
* @param mp         Map containing cst symbol and their values to apply
* @param filename   Shared file
*
* @return Feel++ Expression
*/
inline
Expr< GinacEx<2> > expr( std::string const& s, std::map<std::string,double> const& mp, std::string filename="",
                         WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    auto ginacEx = expr(s,filename, world, dirLibExpr);
    ginacEx.setParameterValues( mp );
    return ginacEx;
}
#if 0
inline
Expr< GinacEx<2> > expr( std::string const& s, std::pair<std::string,double> const& mp, std::string filename="" )
{
    return expr( s, { { mp.first, mp.second } }, filename, world );
}
#endif

/**
 * @brief Create an Feel++ expression from a GiNaC expression as a string
 *
 * @tparam Order     Expression order
 * @param s          String containing the ginac expression and symbols
 * @param mp         Map containing cst symbol and their values to apply
 * @param filename   Shared file
 *
 * @return Feel++ Expression
 */
template<int Order>
inline
Expr< GinacEx<Order> >
expr( std::string const& s, std::map<std::string,double> const& mp, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    auto ginacEx = expr<Order>(s,filename, world, dirLibExpr);
    ginacEx.setParameterValues( mp );
    return ginacEx;
}
template<int Order>
inline
Expr< GinacEx<Order> >
expr( std::string const& s, std::pair<std::string,double> const& mp, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return expr<Order>( s, { { mp.first, mp.second } }, filename, world, dirLibExpr );
}

// ------------------------------------------------------------
// Ginac expression  with feel++ expression
// ------------------------------------------------------------

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( ex const& myexpr, std::vector<GiNaC::symbol> const & syms , GiNaC::symbol const& s, ExprT const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back( std::make_pair(s,e) );
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( myexpr, syms, std::string(""), VFmap, filename, world ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::vector<GiNaC::symbol> const& lsym, std::pair<GiNaC::symbol,ExprT> const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back( e );
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( parse(s,lsym), lsym, s, VFmap, filename, world ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( ex const& myexpr, std::vector<GiNaC::symbol> const & syms , std::initializer_list<GiNaC::symbol> const& s, std::initializer_list<ExprT> const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( s.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::initializer_list<GiNaC::symbol>::iterator it1 = s.begin();
    typename std::initializer_list<ExprT>::iterator it2 = e.begin();
    for(; it1!=s.end(); it1++)
        {
            VFmap.push_back( std::make_pair(*it1,*it2) );
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( myexpr, syms, std::string(""), VFmap, filename, world ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( ex const& myexpr, std::vector<GiNaC::symbol> const & syms , std::vector<GiNaC::symbol> const& s, std::vector<ExprT> const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( s.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::vector<GiNaC::symbol>::const_iterator it1 = s.begin();
    typename std::vector<ExprT>::const_iterator it2 = e.begin();
    for( ; it1!=s.end(); it1++)
        {
            VFmap.push_back( std::make_pair(*it1, *it2) );
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( myexpr, syms, std::string(""), VFmap, filename, world ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& myexpr, std::vector<GiNaC::symbol> const & syms , std::vector<GiNaC::symbol> const& s, std::vector<ExprT> const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( s.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::vector<GiNaC::symbol>::const_iterator it1 = s.begin();
    typename std::vector<ExprT>::const_iterator it2 = e.begin();
    for( ; it1!=s.end(); it1++)
        {
            VFmap.push_back( std::make_pair(*it1, *it2) );
            it2++;
        }
    //std::cout << "ginac ExVF : filename = " << filename << std::endl;
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( parse(myexpr,syms), syms, myexpr, VFmap, filename, world ) );
}

/**
 * @brief Create an Feel++ expression from a GiNaC expression as a string
 *
 * @tparam Order     Expression order
 * @param s          String containing the ginac expression and symbols
 * @param se         String containing the ginac symbol associated to a Feel++ expression (e.g. a finite element field)
 * @param filename   Shared file
 *
 * @return Feel++ Expression
 */

template<typename ExprT,int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::string const& se, ExprT const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    auto it = std::find_if( g.second.begin(), g.second.end(),
                            [&se]( GiNaC::symbol const& s ) { return s.get_name() == se; } );
    LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << se << " in expression " << s;
    LOG(INFO) << "g.second.size() " << g.second.size() << "\n";
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back(std::make_pair(*it, e));
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( g.first, g.second, s, VFmap, filename, world ) );
}

template<typename ExprT,int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::initializer_list<std::string> const& se, std::initializer_list<ExprT> const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( se.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::initializer_list<std::string>::iterator it1 = se.begin();
    typename std::initializer_list<ExprT>::iterator it2 = e.begin();
    for( ; it1!=se.end(); it1++)
        {
            auto it = std::find_if( g.second.begin(), g.second.end(),
                                    [&it1]( GiNaC::symbol const& s ) { return s.get_name() == *it1; } );
            LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << *it1 << " in expression " << s;
            VFmap.push_back(std::make_pair(*it, *it2));
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( g.first, g.second, s, VFmap, filename, world ) );
}

template<typename ExprT,int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::vector<std::string> const& se, std::vector<ExprT> const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( se.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::vector<std::string>::const_iterator it1 = se.begin();
    typename std::vector<ExprT>::const_iterator it2 = e.begin();
    for( ; it1!=se.end(); it1++)
        {
            auto it = std::find_if( g.second.begin(), g.second.end(),
                                    [&it1]( GiNaC::symbol const& s ) { return s.get_name() == *it1; } );
            LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << *it1 << " in expression " << s;
            VFmap.push_back(std::make_pair(*it, *it2));
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( g.first, g.second, s, VFmap, filename, world ) );
}

template<typename ExprT,int Order>
inline
Expr< GinacExVF<ExprT,Order> >
expr( Expr<GinacEx<Order>> const& s, std::string const& se, ExprT const& e )
{
    auto const& symbexpr = s.expression();
    auto it = std::find_if( symbexpr.symbols().begin(), symbexpr.symbols().end(),
                            [&se]( GiNaC::symbol const& s ) { return s.get_name() == se; } );
    CHECK( it != symbexpr.symbols().end() ) << "invalid symbol " << se << " in expression ";// << s;
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back(std::make_pair(*it, e));
    auto res = Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( symbexpr.expression(), symbexpr.symbols(), symbexpr.fun(), VFmap, symbexpr.exprDesc() ) );
    res.expression().setParameterValues( symbexpr.parameterValue() );
    return res;
}

/**
 * @brief Create an Feel++ expression from a GiNaC expression as a string
 *
 * @tparam Order     Expression order
 * @param s          String containing the ginac expression and symbols
 * @param ds         String containing the ginac symbol with respect to which \p s is differentiated
 * @param se         String containing the ginac symbol associated to a Feel++ expression (e.g. a finite element field)
 * @param filename   Shared file
 *
 * @return Feel++ Expression
 */
template<typename ExprT,int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
diff( std::string const& s, std::string const& ds, std::string const& se, ExprT const& e, std::string filename="", WorldComm const& world=Environment::worldComm() )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    auto it = std::find_if( g.second.begin(), g.second.end(),
                            [&se]( GiNaC::symbol const& s ) { return s.get_name() == se; } );
    auto diff_it = std::find_if(g.second.begin(), g.second.end(),
                                [&ds]( GiNaC::symbol const& s ) { return s.get_name() == ds; } );
    LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << se << " in expression " << s;
    LOG_IF( WARNING, (diff_it == g.second.end() ) ) << "invalid symbol " << ds << " in expression " << s << " for differentiation";
    auto diffe = diff(g.first,*diff_it);
    LOG(INFO) << "diff(" << s << "," << ds << ")=" << diffe;

    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back(std::make_pair(*it, e));
    std::string exprDesc = (boost::format("diff(%1%,%2%)")%s %ds ).str();
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( diffe, g.second, exprDesc, VFmap, filename, world ) );
}

// ------------------------------------------------------------
// Matrix expression
// ------------------------------------------------------------
/**
 * @brief Create an Feel++ expression from a GiNaC expression as a string
 *
 * @tparam Order     Expression order
 * @param s          String containing the ginac expression and symbols
 * @param filename   Shared file
 *
 * @return Feel++ Expression
 */
template<int M, int N, int Order=2>
inline
Expr< GinacMatrix<M,N,Order> >
expr( std::string const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( g.first, g.second, s, filename, world, dirLibExpr ) );
}
template<int M, int N, int Order=2>
inline
Expr< GinacMatrix<M,N,Order> >
expr( std::string const& s, std::map<std::string,double> const& mp, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    auto ginacMat = expr<M,N,Order>(s,filename, world, dirLibExpr);
    ginacMat.setParameterValues( mp );
    return ginacMat;
}
template<int M, int N, int Order=2>
inline
Expr< GinacMatrix<M,N,Order> >
expr( std::string const& s, std::pair<std::string,double> const& mp, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return expr<M,N,Order>( s, { { mp.first, mp.second } }, filename, world, dirLibExpr );
}


inline
Expr< GinacMatrix<1,1,2> >
expr( GiNaC::matrix const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    VLOG(2) << "expr(Ginac::matrix(1,1)\n";
    return Expr< GinacMatrix<1,1,2> >(  GinacMatrix<1,1,2>( f, lsym, std::string(""), filename, world, dirLibExpr ) );
}

/**
 * \brief functor enabling ginac
 *
 */
template<int M, int N, int Order>
inline
Expr< GinacMatrix<M,N,Order> >
expr( GiNaC::matrix const& f, std::vector<GiNaC::symbol> const& lsym, std::string const& exprDesc="", std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( f, lsym, exprDesc, filename, world, dirLibExpr ) );
}

template<int M, int N, int Order>
inline
Expr< GinacMatrix<M,N,Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string const& exprDesc="", std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( f, lsym, exprDesc, filename, world, dirLibExpr ) );
}

template<int M,int Order=2>
inline
Expr<GinacMatrix<1,M,Order> >
grad( Expr<GinacEx<Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("grad(%1%)")% s.expression().exprDesc() ).str();
    return expr<1,M,Order>( GiNaC::grad(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}

template<int M,int Order=2>
inline
Expr<GinacMatrix<M,M,Order> >
grad( Expr<GinacMatrix<M,1,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("grad(%1%)")% s.expression().exprDesc() ).str();
    return expr<M,M,Order>( GiNaC::grad(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
// Divergence
template<int M=1,int Order=2>
inline
Expr<GinacMatrix<M,1,Order> >
div( Expr<GinacEx<Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("div(%1%)")% s.expression().exprDesc() ).str();
    return expr<M,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
template<int M,int Order=2>
inline
Expr<GinacMatrix<1,1,Order> >
div( Expr<GinacMatrix<M,1,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("div(%1%)")% s.expression().exprDesc() ).str();
    return expr<1,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
template<int M,int Order=2>
inline
Expr<GinacMatrix<1,1,Order> >
div( Expr<GinacMatrix<1,M,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("div(%1%)")% s.expression().exprDesc() ).str();
    return expr<1,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
template<int M,int N,int Order=2>
inline
Expr<GinacMatrix<M,1,Order> >
div( Expr<GinacMatrix<M,M,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("div(%1%)")% s.expression().exprDesc() ).str();
    return expr<M,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
// Curl
template<int M=1,int Order=2>
inline
Expr<GinacMatrix<M,1,Order> >
curl( Expr<GinacEx<Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("curl(%1%)")% s.expression().exprDesc() ).str();
    return expr<M,1,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
template<int M,int Order=2>
inline
Expr<GinacMatrix<((M==2)?1:3),1,Order> >
curl( Expr<GinacMatrix<M,1,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("curl(%1%)")% s.expression().exprDesc() ).str();
    return expr<((M==2)?1:3),1,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}
template<int Order=2>
inline
Expr<GinacMatrix<2,1,Order> >
curl( Expr<GinacMatrix<1,1,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("curl(%1%)")% s.expression().exprDesc() ).str();
    return expr<2,1,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}

template<int M,int Order=2>
inline
Expr<GinacMatrix<1,M,Order> >
curl( Expr<GinacMatrix<1,M,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("curl(%1%)")% s.expression().exprDesc() ).str();
    return expr<1,M,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}

// Laplacian
template<int Order=2>
inline
Expr<GinacMatrix<1,1,Order> >
laplacian( Expr<GinacEx<Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("laplacian(%1%)")% s.expression().exprDesc() ).str();
    return expr<1,1,Order>( GiNaC::laplacian(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}

template<int M,int N=1, int Order=2>
inline
Expr<GinacMatrix<M,N,Order> >
laplacian( Expr<GinacMatrix<M,N,Order>> const& s, std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("laplacian(%1%)")% s.expression().exprDesc() ).str();
    return expr<M,N,Order>( GiNaC::laplacian(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}

template<int Order=2>
inline
Expr<GinacEx<Order> >
diff( Expr<GinacEx<Order>> const& s, std::string const& diffVariable, int diffOrder = 1, std::string filename="",
      WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    std::string exprDesc = (boost::format("diff(%1%)_%2%_o%3%")% s.expression().exprDesc() %diffVariable %diffOrder ).str();
    auto const& exprSymb = s.expression().symbols();
    auto it = std::find_if( exprSymb.begin(), exprSymb.end(),
                            [&diffVariable]( GiNaC::symbol const& s ) { return s.get_name() == diffVariable; } );
    CHECK( it != exprSymb.end() ) << "symbol " << diffVariable << "not found in expression";
    return expr<Order>( GiNaC::diff( s.expression().expression(),{ *it },diffOrder)(0,0),
                        s.expression().symbols(), exprDesc, filename, world, dirLibExpr );
}


template<int Order,int Order1,int Order2>
inline
Expr<GinacEx<Order> >
expr_mult( Expr<GinacEx<Order1>> const& e1, Expr<GinacEx<Order2>> const& e2, std::string filename="",
           WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    GiNaC::ex newex = e1.expression().expression()*e2.expression().expression();
    std::string exprDesc = str( newex );
    std::set<std::string> exprSymb;
    for ( GiNaC::symbol const& s1 : e1.expression().symbols() )
        exprSymb.insert( s1.get_name() );
    for ( GiNaC::symbol const& s2 : e2.expression().symbols() )
        exprSymb.insert( s2.get_name() );
    for ( std::string const& s : exprSymb )
        exprDesc += (":"+s);
    return expr<Order>( exprDesc, filename, world, dirLibExpr );
}

template<int Order>
inline
Expr<GinacEx<Order> >
expr_mult( Expr<GinacEx<Order>> const& e1, double e2, std::string filename="",
           WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr="" )
{
    GiNaC::ex newex = e1.expression().expression()*e2;
    std::string exprDesc = str( newex );
    std::set<std::string> exprSymb;
    for ( GiNaC::symbol const& s1 : e1.expression().symbols() )
        exprSymb.insert( s1.get_name() );
    for ( std::string const& s : exprSymb )
        exprDesc += (":"+s);
    return expr<Order>( exprDesc, filename, world, dirLibExpr );
}




template<int Order=2>
using scalar_field_expression=Expr<GinacEx<Order>>;

/**
 * defines a dictionary of scalar fields
 * 
 * this data structure creates a dictionary of scalar fields, it associates a
 * string to a Ginac Expr of rank 0.
 * 
 * \code
 * auto e = expr("x+y:x:y");
 * auto z = expr("0:x:y");
 * map_scalar_field m { { "inlet", e }, { "wall", z } };
 * \endcode
 */
template<int Order=2>
struct map_scalar_field: public std::map<std::string,scalar_field_expression<Order>>
{
    typedef std::map<std::string,scalar_field_expression<Order>> super;
    typedef super type;
    using value_type = typename super::value_type;
    map_scalar_field() = default;
    map_scalar_field(std::initializer_list<value_type> __l ) : super( __l ) {}
    map_scalar_field(map_scalar_field&& f ) = default;
    map_scalar_field(map_scalar_field const& f ) = default;
    map_scalar_field& operator=(map_scalar_field && f ) = default;
    map_scalar_field& operator=(map_scalar_field const& f ) = default;
    void setParameterValues( std::map<std::string,double> const& pv )
    {
        for( auto & f : *this )
            f.second.setParameterValues( pv );
    }
    
};

typedef std::map<std::string,Expr<GinacEx<2>>> map_scalar_field_type;

template<int Order=2>
struct map_scalar_fields: public std::map<std::string,std::vector<scalar_field_expression<Order>>>
{
    typedef std::map<std::string,std::vector<scalar_field_expression<Order>>> super;
    typedef super type;
    using value_type = typename super::value_type;
    map_scalar_fields() = default;
    map_scalar_fields(std::initializer_list<value_type> __l ) : super( __l ) {}
    map_scalar_fields(map_scalar_fields&& f ) = default;
    map_scalar_fields(map_scalar_fields const& f ) = default;
    map_scalar_fields& operator=(map_scalar_fields && f ) = default;
    map_scalar_fields& operator=(map_scalar_fields const& f ) = default;
    void setParameterValues( std::map<std::string,double> const& pv )
    {
        for( auto & f : *this )
            for ( auto & g : f.second )
                g.setParameterValues( pv );
    }
};

template<int Order>
std::string const&
marker( std::pair<const std::string, scalar_field_expression<Order>> const& p  )
{
    return p.first;
}
template<int Order>
scalar_field_expression<Order> const&
expression( std::pair<const std::string, scalar_field_expression<Order>> const& p  ) 
{
    return p.second;
}
template<int Order>
scalar_field_expression<Order>&
expression( std::pair<const std::string, scalar_field_expression<Order>> & p  ) 
{
    return p.second;
}

template<int Order>
std::string const&
marker( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> const& p  )
{
    return p.first;
}
template<int Order>
scalar_field_expression<Order> const&
expression1( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> const& p  ) 
{
    return p.second[0];
}
template<int Order>
scalar_field_expression<Order>&
expression1( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> & p  ) 
{
    return p.second[0];
}

template<int Order>
scalar_field_expression<Order> const&
expression2( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> const& p  ) 
{
    return p.second[1];
}
template<int Order>
scalar_field_expression<Order>&
expression2( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> & p  ) 
{
    return p.second[1];
}

template<int Order>
std::vector<scalar_field_expression<Order>> const&
expression( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> const& p  ) 
{
    return p.second;
}

template<int Order>
std::vector<scalar_field_expression<Order>>&
expression( std::pair<const std::string, std::vector<scalar_field_expression<Order>>> & p  ) 
{
    return p.second;
}


template<int M, int N=1, int Order=2>
using vector_field_expression=Expr<GinacMatrix<M,N,Order>>;

template<int M, int N, int Order=2>
using matrix_field_expression=vector_field_expression<M,N,Order>;

/**
 * defines a dictionary of vector fields
 * 
 * this data structure creates a dictionary of fields, it associates a string to
 * a Ginac Expr of rank 1. In the case of vector fields the size of the vector
 * field must be given
 * 
 * \code
 * auto e = expr("{x,y}:x:y");
 * auto z = expr("{0,0}:x:y");
 * map_vector_field<2> m { { "inlet", e }, { "wall", z } };
 * \endcode
 */
template<int M, int N=1, int Order=2>
struct map_vector_field: public std::map<std::string,Expr<GinacMatrix<M,N,Order>>>
{
    typedef std::map<std::string,Expr<GinacMatrix<M,N,Order>>> super;
    typedef super type;
    using value_type = typename super::value_type;
    map_vector_field() = default;
    map_vector_field(std::initializer_list<value_type> __l ) : super( __l ) {}
    map_vector_field(map_vector_field&& f ) = default;
    map_vector_field(map_vector_field const& f ) = default;
    map_vector_field& operator=(map_vector_field && f ) = default;
    map_vector_field& operator=(map_vector_field const& f ) = default;
    void setParameterValues( std::map<std::string,double> const& pv )
    {
        for( auto & f : *this )
            f.second.setParameterValues( pv );
    }    
};

template<int M, int N=1, int Order=2>
struct map_vector_fields: public std::map<std::string,std::vector<Expr<GinacMatrix<M,N,Order>>>>
{
    typedef std::map<std::string,Expr<std::vector<GinacMatrix<M,N,Order>>>> super;
    typedef super type;
    using value_type = typename super::value_type;
    map_vector_fields() = default;
    map_vector_fields(std::initializer_list<value_type> __l ) : super( __l ) {}
    map_vector_fields(map_vector_fields&& f ) = default;
    map_vector_fields(map_vector_fields const& f ) = default;
    map_vector_fields& operator=(map_vector_fields && f ) = default;
    map_vector_fields& operator=(map_vector_fields const& f ) = default;
    void setParameterValues( std::map<std::string,double> const& pv )
    {
        for( auto & f : *this )
            for ( auto & g : f.second )
                g.setParameterValues( pv );
    }
};

/**
 * define a matrix field map. providing M and N is required
 */
template<int M, int N, int Order=2>
using map_matrix_field = map_vector_field<M,N,Order>;

template<int M, int N, int Order>
std::string const&
marker( std::pair<const std::string, Expr<GinacMatrix<M,N,Order>>> const& p  )
{
    return p.first;
}
template<int M, int N, int Order>
Expr<GinacMatrix<M,N,Order>> const&
expression( std::pair<const std::string, Expr<GinacMatrix<M,N,Order>>> const& p  ) 
{
    return p.second;
}
template<int M, int N, int Order>
Expr<GinacMatrix<M,N,Order>> &
expression( std::pair<const std::string, Expr<GinacMatrix<M,N,Order>>> & p  ) 
{
    return p.second;
}

template<int M, int N, int Order>
std::string const&
marker( std::pair<const std::string, std::vector<Expr<GinacMatrix<M,N,Order>>>> const& p  )
{
    return p.first;
}
template<int M, int N, int Order>
Expr<GinacMatrix<M,N,Order>> const&
expression1( std::pair<const std::string, std::vector<Expr<GinacMatrix<M,N,Order>>>>  const& p  )
{
    return p.second[0];
}
template<int M, int N, int Order>
Expr<GinacMatrix<M,N,Order>> &
expression1( std::pair<const std::string, std::vector<Expr<GinacMatrix<M,N,Order>>>> & p  )
{
    return p.second[0];
}
template<int M, int N, int Order>
Expr<GinacMatrix<M,N,Order>> const&
expression2( std::pair<const std::string, std::vector<Expr<GinacMatrix<M,N,Order>>>> const& p  )
{
    return p.second[1];
}
template<int M, int N, int Order>
Expr<GinacMatrix<M,N,Order>> &
expression2( std::pair<const std::string, std::vector<Expr<GinacMatrix<M,N,Order>>>> & p  )
{
    return p.second[1];
}

} // vf
} // Feel

#endif /* FEELPP_GINAC_HPP */
