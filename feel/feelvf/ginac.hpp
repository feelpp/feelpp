/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
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
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-10-15
*/
#ifndef FEELPP_GINAC_HPP
#define FEELPP_GINAC_HPP 1

#include <ginac/ginac.h>
#include <boost/parameter/preprocessor.hpp>

#include <boost/foreach.hpp>
#include <boost/range/algorithm/for_each.hpp>

namespace GiNaC
{
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

ex diff(ex const& f, symbol const& l, const int n);
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

class GinacExprManagerImpl :
        public std::map<std::string, boost::shared_ptr<GiNaC::FUNCP_CUBA> > ,
        public boost::noncopyable
{
public:
    typedef boost::shared_ptr<GiNaC::FUNCP_CUBA> value_type;
    typedef std::string key_type;
    typedef std::map<key_type,value_type> ginac_expr_manager_type;
};

typedef Feel::Singleton<GinacExprManagerImpl> GinacExprManager;

struct GinacExprManagerDeleterImpl
{
    void operator()() const
        {
            VLOG(2) << "[GinacManagerDeleter] clear GinacExprManager Singleton: " << GinacExprManager::instance().size() << "\n";
            GinacExprManager::instance().clear();
            VLOG(2) << "[GinacManagerDeleter] clear GinacExprManager done\n";
        }
};
typedef Feel::Singleton<GinacExprManagerDeleterImpl> GinacExprManagerDeleter;
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
#include <feel/feelvf/detail/ginacex.hpp>
#include <feel/feelvf/detail/ginacexvf.hpp>
#include <feel/feelvf/detail/ginacmatrix.hpp>

namespace Feel
{
namespace vf
{
inline
Expr< GinacEx<2> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="" )
{
    return Expr< GinacEx<2> >(  GinacEx<2>( f, lsym, filename ) );
}

inline
Expr< GinacEx<2> >
expr( std::string const& s, std::vector<GiNaC::symbol> const& lsym, std::string filename="" )
{
    return Expr< GinacEx<2> >(  GinacEx<2>( parse(s,lsym), lsym, filename) );
}


/**
 * \brief functor enabling ginac
 *
 */
template<int Order>
inline
Expr< GinacEx<Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="" )
{
    return Expr< GinacEx<Order> >(  GinacEx<Order>( f, lsym, filename ));
}

template<int Order>
inline
Expr< GinacEx<Order> >
expr( std::string const& s, std::vector<GiNaC::symbol> const& lsym, std::string filename="" )
{
    return Expr< GinacEx<Order> >(  GinacEx<Order>( parse(s,lsym), lsym, filename) );
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
Expr< GinacEx<2> > expr( std::string const& s, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacEx<2> >(  GinacEx<2>( g.first, g.second, filename) );
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
Expr< GinacEx<2> > expr(  const char* s, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacEx<2> >(  GinacEx<2>( g.first, g.second, filename) );
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
expr( std::string const& s, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacEx<Order> >(  GinacEx<Order>( g.first, g.second, filename) );
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
Expr< GinacEx<2> > expr( std::string const& s, std::map<std::string,double> const& mp, std::string filename="" )
{
    auto ginacEx = expr(s,filename);
    ginacEx.setParameterValues( mp );
    return ginacEx;
}
#if 0
inline
Expr< GinacEx<2> > expr( std::string const& s, std::pair<std::string,double> const& mp, std::string filename="" )
{
    return expr( s, { { mp.first, mp.second } }, filename );
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
expr( std::string const& s, std::map<std::string,double> const& mp, std::string filename="" )
{
    auto ginacEx = expr<Order>(s,filename);
    ginacEx.setParameterValues( mp );
    return ginacEx;
}
template<int Order>
inline
Expr< GinacEx<Order> >
expr( std::string const& s, std::pair<std::string,double> const& mp, std::string filename="" )
{
    return expr<Order>( s, { { mp.first, mp.second } }, filename );
}

// ------------------------------------------------------------
// Ginac expression  with feel++ expression
// ------------------------------------------------------------

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( ex const& myexpr, std::vector<GiNaC::symbol> const & syms , GiNaC::symbol const& s, ExprT const& e, std::string filename="" )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back( std::make_pair(s,e) );
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( myexpr, syms, VFmap, filename ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::vector<GiNaC::symbol> const& lsym, std::pair<GiNaC::symbol,ExprT> const& e, std::string filename="" )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back( e );
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( parse(s,lsym), lsym, VFmap, filename ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( ex const& myexpr, std::vector<GiNaC::symbol> const & syms , std::initializer_list<GiNaC::symbol> const& s, std::initializer_list<ExprT> const& e, std::string filename="" )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( s.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::initializer_list<GiNaC::symbol>::iterator it1 = s.begin();
    typename std::initializer_list<ExprT>::iterator it2 = e.begin();
    for(it1; it1!=s.end(); it1++)
        {
            VFmap.push_back( std::make_pair(*it1,*it2) );
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( myexpr, syms, VFmap, filename ) );
}

template<typename ExprT, int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( ex const& myexpr, std::vector<GiNaC::symbol> const & syms , std::vector<GiNaC::symbol> const& s, std::vector<ExprT> const& e, std::string filename="" )
{
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( s.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::vector<GiNaC::symbol>::const_iterator it1 = s.begin();
    typename std::vector<ExprT>::const_iterator it2 = e.begin();
    for(it1; it1!=s.end(); it1++)
        {
            VFmap.push_back( std::make_pair(*it1, *it2) );
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( myexpr, syms, VFmap, filename ) );
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
expr( std::string const& s, std::string const& se, ExprT const& e, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    auto it = std::find_if( g.second.begin(), g.second.end(),
                            [&se]( GiNaC::symbol const& s ) { return s.get_name() == se; } );
    LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << se << " in expression " << s;

    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    VFmap.push_back(std::make_pair(*it, e));
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( g.first, g.second, VFmap, filename) );
}

template<typename ExprT,int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::initializer_list<std::string> const& se, std::initializer_list<ExprT> const& e, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( se.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::initializer_list<std::string>::iterator it1 = se.begin();
    typename std::initializer_list<ExprT>::iterator it2 = e.begin();
    for(it1; it1!=se.end(); it1++)
        {
            auto it = std::find_if( g.second.begin(), g.second.end(),
                                    [&it1]( GiNaC::symbol const& s ) { return s.get_name() == *it1; } );
            LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << *it1 << " in expression " << s;
            VFmap.push_back(std::make_pair(*it, *it2));
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( g.first, g.second, VFmap, filename) );
}

template<typename ExprT,int Order=2>
inline
Expr< GinacExVF<ExprT,Order> >
expr( std::string const& s, std::vector<std::string> const& se, std::vector<ExprT> const& e, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    std::vector< std::pair<GiNaC::symbol,ExprT> > VFmap;
    CHECK( se.size() == e.size() ) << "List of expressions and associated symbols have not the same size \n";

    typename std::vector<std::string>::const_iterator it1 = se.begin();
    typename std::vector<ExprT>::const_iterator it2 = e.begin();
    for(it1; it1!=se.end(); it1++)
        {
            auto it = std::find_if( g.second.begin(), g.second.end(),
                                    [&it1]( GiNaC::symbol const& s ) { return s.get_name() == *it1; } );
            LOG_IF( WARNING, (it == g.second.end() ) ) << "invalid symbol " << *it1 << " in expression " << s;
            VFmap.push_back(std::make_pair(*it, *it2));
            it2++;
        }
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( g.first, g.second, VFmap, filename) );
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
diff( std::string const& s, std::string const& ds, std::string const& se, ExprT const& e, std::string filename="" )
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
    return Expr< GinacExVF<ExprT,Order> >(  GinacExVF<ExprT,Order>( diffe, g.second, VFmap, filename) );
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
expr( std::string const& s, std::string filename="" )
{
    std::pair< ex, std::vector<GiNaC::symbol> > g = GiNaC::parse(s);
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( g.first, g.second, filename) );
}


inline
Expr< GinacMatrix<1,1,2> >
expr( GiNaC::matrix const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="")
{
    return Expr< GinacMatrix<1,1,2> >(  GinacMatrix<1,1,2>( f, lsym, filename ) );
}

/**
 * \brief functor enabling ginac
 *
 */
template<int M, int N, int Order>
inline
Expr< GinacMatrix<M,N,Order> >
expr( GiNaC::matrix const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="" )
{
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( f, lsym, filename) );
}

template<int M, int N, int Order>
inline
Expr< GinacMatrix<M,N,Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym, std::string filename="" )
{
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( f, lsym, filename ) );
}

template<int M,int Order=2>
inline
Expr<GinacMatrix<1,M,Order> >
grad( Expr<GinacEx<Order>> const& s, std::string filename="" )
{
    return expr<1,M,Order>( GiNaC::grad(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}

template<int M,int Order=2>
inline
Expr<GinacMatrix<M,M,Order> >
grad( Expr<GinacMatrix<M,1,Order>> const& s, std::string filename="" )
{
    return expr<M,M,Order>( GiNaC::grad(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
// Divergence
template<int M=1,int Order=2>
inline
Expr<GinacMatrix<M,1,Order> >
div( Expr<GinacEx<Order>> const& s, std::string filename="" )
{
    return expr<M,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
template<int M,int Order=2>
inline
Expr<GinacMatrix<1,1,Order> >
div( Expr<GinacMatrix<M,1,Order>> const& s, std::string filename="" )
{
    return expr<1,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
template<int M,int Order=2>
inline
Expr<GinacMatrix<1,1,Order> >
div( Expr<GinacMatrix<1,M,Order>> const& s, std::string filename="" )
{
    return expr<1,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
template<int M,int N,int Order=2>
inline
Expr<GinacMatrix<M,1,Order> >
div( Expr<GinacMatrix<M,M,Order>> const& s, std::string filename="" )
{
    return expr<M,1,Order>( GiNaC::div(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
// Curl
template<int M=1,int Order=2>
inline
Expr<GinacMatrix<M,1,Order> >
curl( Expr<GinacEx<Order>> const& s, std::string filename="" )
{
    return expr<M,1,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
template<int M,int Order=2>
inline
Expr<GinacMatrix<((M==2)?1:3),1,Order> >
curl( Expr<GinacMatrix<M,1,Order>> const& s, std::string filename="" )
{
    return expr<((M==2)?1:3),1,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}
template<int Order=2>
inline
Expr<GinacMatrix<2,1,Order> >
curl( Expr<GinacMatrix<1,1,Order>> const& s, std::string filename="" )
{
    return expr<2,1,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}

template<int M,int Order=2>
inline
Expr<GinacMatrix<1,M,Order> >
curl( Expr<GinacMatrix<1,M,Order>> const& s, std::string filename="" )
{
    return expr<1,M,Order>( GiNaC::curl(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}

// Laplacian
template<int Order=2>
inline
Expr<GinacMatrix<1,1,Order> >
laplacian( Expr<GinacEx<Order>> const& s, std::string filename="" )
{
    return expr<1,1,Order>( GiNaC::laplacian(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}

template<int M,int N=1, int Order=2>
inline
Expr<GinacMatrix<M,N,Order> >
laplacian( Expr<GinacMatrix<M,N,Order>> const& s, std::string filename="" )
{
    return expr<M,N,Order>( GiNaC::laplacian(s.expression().expression(),s.expression().symbols()), s.expression().symbols(), filename );
}

} // vf
} // Feel

#endif /* FEELPP_GINAC_HPP */
