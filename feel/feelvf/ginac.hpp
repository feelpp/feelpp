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
#ifndef __Ginac_H
#define __Ginac_H 1

#include <ginac/ginac.h>
#include <boost/parameter/preprocessor.hpp>

namespace GiNaC
{

    matrix grad( ex const& f, std::vector<symbol> const& l );
    ex laplacian( ex const& f, std::vector<symbol> const& l );
    matrix grad( std::string const& s, std::vector<symbol> const& l );
    ex laplacian( std::string const& s, std::vector<symbol> const& l );

    matrix grad( matrix const& f, std::vector<symbol> const& l );
    matrix div( matrix const& f, std::vector<symbol> const& l );
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
using  GiNaC::matrix;
using  GiNaC::symbol;
using  GiNaC::lst;
using  GiNaC::ex;
using  GiNaC::parser;

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
#if 0
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
#endif
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


namespace vf
{

/// \cond detail
/**
 * \class Ginac
 * \brief allow runtime ginac in expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order = 2>
class GinacEx
{
public:


    /** @name Typedefs
     */
    //@{


    static const size_type context = vm::POINT;
    static const bool is_terminal = false;
    static const uint16_type imorder = Order;
    static const bool imIsPoly = false;

    template<typename Funct>
    struct HasTestFunction
    {
        static const bool result = false;
    };
    template<typename Funct>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    typedef GiNaC::ex expression_type;
    typedef GinacEx<Order> this_type;
    typedef double value_type;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;

    template<typename TheExpr>
    struct Lambda
    {
        typedef this_type type;
    };
    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  ) { return *this; }

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit GinacEx( expression_type const & fun, std::vector<GiNaC::symbol> const& syms, std::string filename="",
                      WorldComm const& world=Environment::worldComm() )
        :
        M_fun( fun ),
        M_syms( syms),
        M_params( vec_type::Zero( M_syms.size() ) ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(filename.empty()?filename:(fs::current_path()/filename).string())
        {
            DVLOG(2) << "Ginac constructor with expression_type \n";
            GiNaC::lst exprs(fun);
            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            std::string filenameWithSuffix = M_filename + ".so";

            // register the GinacExprManager into Feel::Environment so that it gets the
            // GinacExprManager is cleared up when the Environment is deleted
            static bool observed=false;
            if ( !observed )
            {
                Environment::addDeleteObserver( GinacExprManagerDeleter::instance() );
                observed = true;
            }

            bool hasLinked = ( GinacExprManager::instance().find( filename ) != GinacExprManager::instance().end() )? true : false;
            if ( hasLinked )
            {
                M_cfun = GinacExprManager::instance().find( filename )->second;
            }
            else
            {
                // master rank check if the lib exist and compile this one if not done
                if ( ( world.isMasterRank() && !fs::exists( filenameWithSuffix ) ) || M_filename.empty() )
                {
                    DVLOG(2) << "GiNaC::compile_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
                    GiNaC::compile_ex(exprs, syml, *M_cfun, M_filename);
                    hasLinked=true;
                    if ( !M_filename.empty() )
                        GinacExprManager::instance().operator[]( filename ) = M_cfun;
                }
                // wait the lib compilation
                if ( !M_filename.empty() ) world.barrier();
                // link with other process
                if ( !hasLinked )
                {
                    DVLOG(2) << "GiNaC::link_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
                    GiNaC::link_ex(filenameWithSuffix, *M_cfun);
                    if ( !M_filename.empty() )
                        GinacExprManager::instance().operator[]( filename ) = M_cfun;
                }
            }
            this->setParameterFromOption();
        }

    GinacEx( GinacEx const & fun )
    :
    M_fun( fun.M_fun ),
    M_syms( fun.M_syms),
    M_params( fun.M_params ),
    M_cfun( fun.M_cfun ),
    M_filename( fun.M_filename )
        {
            if( !(M_fun==fun.M_fun && M_syms==fun.M_syms && M_filename==fun.M_filename) || M_filename.empty() )
            {
                DVLOG(2) << "Ginac copy constructor : compile object file \n";
                GiNaC::lst exprs(M_fun);
                GiNaC::lst syml;
                std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );
                GiNaC::compile_ex(exprs, syml, *M_cfun, M_filename);

            }
            else
            {
#if 0
                DVLOG(2) << "Ginac copy constructor : link with existing object file \n";
                boost::mpi::communicator world;
                // std::string pid = boost::lexical_cast<std::string>(world.rank());
                // std::string filenameWithSuffix = M_filename + pid + ".so";
                std::string filenameWithSuffix = M_filename + ".so";
                GiNaC::link_ex(filenameWithSuffix, *M_cfun);
#endif
            }
        }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    vec_type const& parameterValue() const { return M_params; }
    value_type parameterValue( int p ) const { return M_params[p]; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setParameterFromOption()
        {
            std::map<std::string,value_type> m;
            for( auto const& s : M_syms )
            {
                if ( Environment::vm().count( s.get_name() ) )
                {
                    // use try/catch in order to catch casting exception for
                    // option that do not return double. Indeed we are only
                    // collecting symbols in option database which can be cast
                    // to numerical types
                    try
                    {
                        value_type v = option( _name=s.get_name() ).template as<double>();
                        m.insert( std::make_pair( s.get_name(), v ) );
                        LOG(INFO) << "symbol " << s.get_name() << " found in option with value " << v;
                    }
                    catch(...)
                    {}
                }
            }
            this->setParameterValues( m );
        }

    void setParameterValues( vec_type const& p )
        {
            CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
            M_params = p;
        }
    void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
            for( auto p : mp )
            {
                auto it = std::find_if( M_syms.begin(), M_syms.end(),
                                        [&p]( GiNaC::symbol const& s ) { return s.get_name() == p.first; } );
                if ( it != M_syms.end() )
                {
                    M_params[it-M_syms.begin()] = p.second;
                    LOG(INFO) << "setting parameter : " << p.first << " with value: " << p.second;
                    LOG(INFO) << "parameter: " << M_params;
                }
                else
                {
                    LOG(INFO) << "Invalid parameters : " << p.first << " with value: " << p.second;
                }
            }
        }

    //@}

    /** @name  Methods
     */
    //@{

    const GiNaC::FUNCP_CUBA& fun() const
        {
            return *M_cfun;
        }

    std::vector<GiNaC::symbol> const& syms() const { return M_syms; }

    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        //typedef typename expression_type::value_type value_type;
        typedef double value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
    // change 0 into rank
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<0>,mpl::int_<0> >,
                              mpl::identity<Shape<gmc_type::nDim, Scalar, false, false> >,
                              typename mpl::if_<mpl::equal_to<mpl::int_<0>,mpl::int_<1> >,
                                                mpl::identity<Shape<gmc_type::nDim, Vectorial, false, false> >,
                                                mpl::identity<Shape<gmc_type::nDim, Tensor2, false, false> > >::type >::type::type shape;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;

    struct is_zero
    {
        static const bool value = false;
    };

    tensor( this_type const& expr,
            Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        :
        M_fun( expr.fun() ),
        M_gmc( fusion::at_key<key_type>( geom ).get() ),
        M_nsyms( expr.syms().size() ),
        M_y( vec_type::Zero(M_gmc->nPoints()) ),
        M_x( expr.parameterValue() )
        {
        }

    tensor( this_type const& expr,
            Geo_t const& geom, Basis_i_t const& /*fev*/ )
        :
        M_fun( expr.fun() ),
        M_gmc( fusion::at_key<key_type>( geom ).get() ),
        M_nsyms( expr.syms().size() ),
        M_y( vec_type::Zero(M_gmc->nPoints()) ),
        M_x(  expr.parameterValue() )
        {
        }

    tensor( this_type const& expr, Geo_t const& geom )
        :
        M_fun( expr.fun() ),
        M_gmc( fusion::at_key<key_type>( geom ).get() ),
        M_nsyms( expr.syms().size() ),
        M_y( vec_type::Zero(M_gmc->nPoints()) ),
        M_x( expr.parameterValue() )

        {
        }
    template<typename IM>
    void init( IM const& im )
        {

        }
    void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
    void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
    void update( Geo_t const& geom )
        {
            M_gmc =  fusion::at_key<key_type>( geom ).get();

            int no = 1;
            int ni = M_nsyms;///gmc_type::nDim;

            for(int q = 0; q < M_gmc->nPoints();++q )
            {
                for(int k = 0;k < gmc_type::nDim;++k )
                    M_x[k]=M_gmc->xReal( q )[k];
                M_fun(&ni,M_x.data(),&no,&M_y[q]);
            }

        }

    void update( Geo_t const& geom, uint16_type /*face*/ )
        {
            M_gmc =  fusion::at_key<key_type>( geom ).get();

            int no = 1;
            int ni = M_nsyms;//gmc_type::nDim;
            for(int q = 0; q < M_gmc->nPoints();++q )
            {
                for(int k = 0;k < gmc_type::nDim;++k )
                    M_x[k]=M_gmc->xReal( q )[k];
                M_fun(&ni,M_x.data(),&no,&M_y[q]);
            }
        }

    template<typename CTX>
    void updateContext( CTX const& ctx )
        {
            update( ctx->gmContext() );
        }

    value_type
    evalij( uint16_type i, uint16_type j ) const
        {
            return 0;
        }


    value_type
    evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q];
        }



    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q];
        }

    value_type
    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q];
        }

    GiNaC::FUNCP_CUBA M_fun;
    gmc_ptrtype M_gmc;

    int M_nsyms;
    vec_type M_y;
    vec_type M_x;

};

private:
mutable expression_type  M_fun;
std::vector<GiNaC::symbol> M_syms;
vec_type M_params;
boost::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
std::string M_filename;
};

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

template<int M=1, int N=1, int Order = 2>
class GinacMatrix
{
public:


    /** @name Typedefs
     */
    //@{


    static const size_type context = vm::POINT;
    static const bool is_terminal = false;
    static const uint16_type imorder = Order;
    static const bool imIsPoly = false;

    template<typename Funct>
    struct HasTestFunction
    {
        static const bool result = false;
    };
    template<typename Funct>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    typedef GiNaC::ex expression_type;
    typedef GinacMatrix<M,N,Order> this_type;
    typedef double value_type;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;

    template<typename TheExpr>
    struct Lambda
    {
        typedef this_type type;
    };
    template<typename TheExpr>
    typename Lambda<TheExpr>::type
    operator()( TheExpr const& e  ) { return *this; }

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit GinacMatrix( GiNaC::matrix const & fun, std::vector<GiNaC::symbol> const& syms, std::string filename="",
                          WorldComm const& world=Environment::worldComm() )
        :
        M_fun( fun.evalm() ),
        M_syms( syms),
        M_params( vec_type::Zero( M_syms.size() ) ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(filename.empty()?filename:(fs::current_path()/filename).string())
        {
            DVLOG(2) << "Ginac matrix constructor with expression_type \n";
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            std::string filenameWithSuffix = M_filename + ".so";

            // register the GinacExprManager into Feel::Environment so that it gets the
            // GinacExprManager is cleared up when the Environment is deleted
            static bool observed=false;
            if ( !observed )
            {
                Environment::addDeleteObserver( GinacExprManagerDeleter::instance() );
                observed = true;
            }

            bool hasLinked = ( GinacExprManager::instance().find( filename ) != GinacExprManager::instance().end() )? true : false;
            if ( hasLinked )
            {
                M_cfun = GinacExprManager::instance().find( filename )->second;
            }
            else
            {
                // master rank check if the lib exist and compile this one if not done
                if ( ( world.isMasterRank() && !fs::exists( filenameWithSuffix ) ) || M_filename.empty() )
                {
                    DVLOG(2) << "GiNaC::compile_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
                    GiNaC::compile_ex(exprs, syml, *M_cfun, M_filename);
                    hasLinked=true;
                    if ( !M_filename.empty() )
                        GinacExprManager::instance().operator[]( filename ) = M_cfun;
                }
                // wait the lib compilation
                if ( !M_filename.empty() ) world.barrier();
                // link with other process
                if ( !hasLinked )
                {
                    DVLOG(2) << "GiNaC::link_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
                    GiNaC::link_ex(filenameWithSuffix, *M_cfun);
                    if ( !M_filename.empty() )
                        GinacExprManager::instance().operator[]( filename ) = M_cfun;
                }
            }
        }
    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms, std::string filename="",
                          WorldComm const& world=Environment::worldComm() )
        :
        M_fun(fun.evalm()),
        M_syms( syms),
        M_params( vec_type::Zero( M_syms.size() ) ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(filename.empty()?filename:(fs::current_path()/filename).string())
        {
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            std::string filenameWithSuffix = M_filename + ".so";

            // register the GinacExprManager into Feel::Environment so that it gets the
            // GinacExprManager is cleared up when the Environment is deleted
            static bool observed=false;
            if ( !observed )
            {
                Environment::addDeleteObserver( GinacExprManagerDeleter::instance() );
                observed = true;
            }

            bool hasLinked = ( GinacExprManager::instance().find( filename ) != GinacExprManager::instance().end() )? true : false;
            if ( hasLinked )
            {
                M_cfun = GinacExprManager::instance().find( filename )->second;
            }
            else
            {
                // master rank check if the lib exist and compile this one if not done
                if ( ( world.isMasterRank() && !fs::exists( filenameWithSuffix ) ) || M_filename.empty() )
                {
                    DVLOG(2) << "GiNaC::compile_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
                    GiNaC::compile_ex(exprs, syml, *M_cfun, M_filename);
                    hasLinked=true;
                    if ( !M_filename.empty() )
                        GinacExprManager::instance().operator[]( filename ) = M_cfun;
                }
                // wait the lib compilation
                if ( !M_filename.empty() ) world.barrier();
                // link with other process
                if ( !hasLinked )
                {
                    DVLOG(2) << "GiNaC::link_ex with filenameWithSuffix " << filenameWithSuffix << "\n";
                    GiNaC::link_ex(filenameWithSuffix, *M_cfun);
                    if ( !M_filename.empty() )
                        GinacExprManager::instance().operator[]( filename ) = M_cfun;
                }
            }
        }

    GinacMatrix( GinacMatrix const & fun )
    :
    M_fun( fun.M_fun ),
    M_syms( fun.M_syms),
    M_params( fun.M_params ),
    M_cfun( fun.M_cfun ),
    M_filename( fun.M_filename )
        {
            if( !(M_fun==fun.M_fun && M_syms==fun.M_syms && M_filename==fun.M_filename) || M_filename.empty() )
            {
                DVLOG(2) << "Ginac copy constructor : compile object file \n";
                GiNaC::lst exprs;
                for( int i = 0; i < fun.M_fun.nops(); ++i ) exprs.append( fun.M_fun.op(i) );

                GiNaC::lst syml;
                std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );
                GiNaC::compile_ex(exprs, syml, *M_cfun, M_filename);
            }
            else
            {
#if 0
                DVLOG(2) << "Ginac copy constructor : link with existing object file \n";
                boost::mpi::communicator world;
                //std::string pid = boost::lexical_cast<std::string>(world.rank());
                std::string filenameWithSuffix = M_filename + ".so";
                GiNaC::link_ex(filenameWithSuffix, *M_cfun);
#endif
            }
        }


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    vec_type const& parameterValue() const { return M_params; }
    value_type parameterValue( int p ) const { return M_params[p]; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setParameterValues( vec_type const& p )
        {
            CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
            M_params = p;
        }
    void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            CHECK( M_params.size() == M_syms.size() ) << "Invalid number of parameters " << M_params.size() << " >= symbol size : " << M_syms.size();
            for( auto p : mp )
            {
                auto it = std::find_if( M_syms.begin(), M_syms.end(),
                                        [&p]( GiNaC::symbol const& s ) { return s.get_name() == p.first; } );
                if ( it != M_syms.end() )
                {
                    M_params[it-M_syms.begin()] = p.second;
                    LOG(INFO) << "setting parameter : " << p.first << " with value: " << p.second;
                    LOG(INFO) << "parameter: " << M_params;
                }
                else
                {
                    LOG(INFO) << "Invalid parameters : " << p.first << " with value: " << p.second;
                }
            }
        }

    //@}

    /** @name  Methods
     */
    //@{


    //@}

    const GiNaC::FUNCP_CUBA& fun() const
        {
            return *M_cfun;
        }

    std::vector<GiNaC::symbol> const& syms() const { return M_syms; }
    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        //typedef typename expression_type::value_type value_type;
        typedef double value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<0> >,
                                  mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        typedef typename mn_to_shape<gmc_type::nDim,M,N>::type shape;
        // be careful that the matrix passed to ginac must be Row Major,
        // however if the number of columns is 1 then eigen3 fails with
        // an assertion, so we have a special when N=1 and have the
        // matrix column major which is ok in this case
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<shape::N>, mpl::int_<1>>,
                                  mpl::identity<Eigen::Matrix<value_type,shape::M,1>>,
                                  mpl::identity<Eigen::Matrix<value_type,shape::M,shape::N,Eigen::RowMajor>>>::type::type mat_type;
        typedef std::vector<mat_type> loc_type;
        typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() )
            {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() )

            {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() )
            {
            }

        template<typename IM>
        void init( IM const& im )
            {

            }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            {
                update( geom );
            }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
            {
                update( geom );
            }
        void update( Geo_t const& geom )
            {
                M_gmc =  fusion::at_key<key_type>( geom ).get();

                int no = M*N;
                int ni = M_nsyms;//gmc_type::nDim;
                for(int q = 0; q < M_gmc->nPoints();++q )
                {
                    for(int k = 0;k < gmc_type::nDim;++k )
                        M_x[k]=M_gmc->xReal( q )[k];
                    for( int k = gmc_type::nDim; k < M_x.size(); ++k )
                        M_x[k] = 0;
                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
                    ;
                }

            }

        void update( Geo_t const& geom, uint16_type /*face*/ )
            {
                M_gmc =  fusion::at_key<key_type>( geom ).get();

                int no = M*N;
                int ni = M_nsyms;//gmc_type::nDim;
                for(int q = 0; q < M_gmc->nPoints();++q )
                {
                    for(int k = 0;k < gmc_type::nDim;++k )
                        M_x[k]=M_gmc->xReal( q )[k];
                    for( int k = gmc_type::nDim; k < M_x.size(); ++k )
                        M_x[k] = 0;
                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
                }
            }


        value_type
        evalij( uint16_type i, uint16_type j ) const
            {
                return 0;
            }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return M_y[q](c1,c2);
            }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return M_y[q](c1,c2);
            }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return M_y[q](c1,c2);
            }

        GiNaC::FUNCP_CUBA M_fun;
        gmc_ptrtype M_gmc;
        int M_nsyms;
        loc_type M_y;
        vec_type M_x;
    };

private:
    mutable expression_type  M_fun;
    std::vector<GiNaC::symbol> M_syms;
    vec_type M_params;
    boost::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
    std::string M_filename;
}; // GinacMatrix
/// \endcond

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


} // vf
    } // Feel
#endif /* __Ginac_H */
