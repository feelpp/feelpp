/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-02-13

  Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_DETAIL_GINACEX_HPP
#define FEELPP_DETAIL_GINACEX_HPP 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * @brief allow runtime ginac in expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order = 2>
class GinacEx : public Feel::vf::GiNaCBase
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
            using namespace GiNaC;
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

//                    try
//                    {
//                        expression_type e( option( _name=s.get_name() ).template as<std::string>(), 0 );
//                        if( is_a<numeric>(e) )
//                        {
//                            LOG(INFO) << "symbol " << s.get_name() << " found in option with value " << v;
//                        }
//                        else
//                        {
//                            ;
//                        }
//                    }
//                    catch(...)
//                    {}
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

    value_type
    evaluate( std::map<std::string,value_type> const& mp  )
    {
        this->setParameterValues( mp );
        int no = 1;
        int ni = M_syms.size();//gmc_type::nDim;
        value_type res;
        (*M_cfun)(&ni,M_params.data(),&no,&res);
        return res;
    }

private:
    mutable expression_type  M_fun;
    std::vector<GiNaC::symbol> M_syms;
    vec_type M_params;
    boost::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
    std::string M_filename;
};

} // vf
} // feel

#endif /* FEELPP_DETAIL_GINACEX_HPP */
