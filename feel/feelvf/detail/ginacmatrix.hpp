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
#ifndef FEELPP_DETAIL_GINACMATRIX_HPP
#define FEELPP_DETAIL_GINACMATRIX_HPP 1

namespace Feel
{
namespace vf {

/**
 * Handle Ginac matrix expression
 */
template<int M=1, int N=1, int Order = 2>
class GinacMatrix : public Feel::vf::GiNaCBase
{
public:


    /** @name Typedefs
     */
    //@{
    typedef Feel::vf::GiNaCBase super;

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

    typedef Eigen::MatrixXd evaluate_type;


    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vec_type;

    template<typename... TheExpr>
    struct Lambda
    {
        typedef this_type type;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) { return *this; }

    //@}

    /** @name Constructors, destructor
     */
    //@{

    GinacMatrix() : super() {}
    explicit GinacMatrix( GiNaC::matrix const & fun, std::vector<GiNaC::symbol> const& syms, std::string const& exprDesc,
                          std::string filename="", WorldComm const& world=Environment::worldComm() )
        :
        super( syms ),
        M_fun( fun.evalm() ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename( Environment::expand( (filename.empty() || fs::path(filename).is_absolute())? filename : (fs::current_path()/filename).string() ) ),
        M_exprDesc( exprDesc )
        {
            DVLOG(2) << "Ginac matrix matrix constructor with expression_type \n";
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            // get filename if not given
            if ( M_filename.empty() && !M_exprDesc.empty() )
            {
                M_filename = Feel::vf::detail::ginacGetDefaultFileName( M_exprDesc );
            }

            // build ginac lib and link if necessary
            Feel::vf::detail::ginacBuildLibrary( exprs, syml, M_exprDesc, M_filename, world, M_cfun );
        }
    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms, std::string const& exprDesc,
                          std::string filename="", WorldComm const& world=Environment::worldComm() )
        :
        super(syms),
        M_fun(fun.evalm()),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename( Environment::expand( (filename.empty() || fs::path(filename).is_absolute())? filename : (fs::current_path()/filename).string() ) ),
        M_exprDesc( exprDesc )
        {
            DVLOG(2) << "Ginac matrix ex constructor with expression_type \n";
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            // get filename if not given
            if ( M_filename.empty() && !M_exprDesc.empty() )
            {
                M_filename = Feel::vf::detail::ginacGetDefaultFileName( M_exprDesc );
            }

            // build ginac lib and link if necessary
            Feel::vf::detail::ginacBuildLibrary( exprs, syml, M_exprDesc, M_filename, world, M_cfun );
        }

    GinacMatrix( GinacMatrix const & fun )
    :
        super(fun),
        M_fun( fun.M_fun ),
        M_cfun( fun.M_cfun ),
        M_filename( fun.M_filename ),
        M_exprDesc( fun.M_exprDesc )
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

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}
    const expression_type& expression() const
        {
            return M_fun;
        }

    const GiNaC::FUNCP_CUBA& fun() const
        {
            return *M_cfun;
        }

    std::vector<GiNaC::symbol> const& syms() const { return M_syms; }

    std::string const& exprDesc() const { return M_exprDesc; }


    Eigen::MatrixXd
    evaluate( std::map<std::string,value_type> const& mp  )
    {
        this->setParameterValues( mp );
        int no = 1;
        int ni = M_syms.size();//gmc_type::nDim;
        Eigen::MatrixXd res(M,N);
        (*M_cfun)(&ni,M_params.data(),&no,res.data());
        return res;
    }

    Eigen::MatrixXd
    evaluate( bool parallel = true, WorldComm const& worldcomm = Environment::worldComm() ) const
    {
        int no = M*N;
        int ni = M_syms.size();
        Eigen::MatrixXd res(M,N);
        (*M_cfun)(&ni,M_params.data(),&no,res.data());
        return res;
    }

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
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() )
            {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() )

            {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
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
                    for ( auto const& comp : M_expr.indexSymbolXYZ() )
                        M_x[comp.second] = M_gmc->xReal( q )[comp.first];
                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
                }

            }

        void update( Geo_t const& geom, uint16_type /*face*/ )
            {
                M_gmc =  fusion::at_key<key_type>( geom ).get();

                int no = M*N;
                int ni = M_nsyms;//gmc_type::nDim;
                for(int q = 0; q < M_gmc->nPoints();++q )
                {
                    for ( auto const& comp : M_expr.indexSymbolXYZ() )
                        M_x[comp.second] = M_gmc->xReal( q )[comp.first];
                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
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

        this_type const& M_expr;
        GiNaC::FUNCP_CUBA M_fun;
        gmc_ptrtype M_gmc;
        int M_nsyms;
        loc_type M_y;
        vec_type M_x;
    };

private:
    mutable expression_type  M_fun;
    boost::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
    std::string M_filename;
    std::string M_exprDesc;
}; // GinacMatrix

template<int M,int N,int Order>
std::ostream&
operator<<( std::ostream& os, GinacMatrix<M,N,Order> const& e )
{
    os << e.expression();
    return os;
}

/// \endcond
/// \endcond
}} // Feel

#endif /* FEELPP_DETAIL_GINACMATRIX_HPP */
