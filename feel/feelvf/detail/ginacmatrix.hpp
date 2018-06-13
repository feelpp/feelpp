/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-02-13

  Copyright (C) 2014-2016 Feel++ Consortium

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
template<int M=1, int N=1, int Order = 2, typename SymbolsExprType = SymbolsExpr<> >
class FEELPP_EXPORT GinacMatrix : public Feel::vf::GiNaCBase
{
public:


    /** @name Typedefs
     */
    //@{
    typedef Feel::vf::GiNaCBase super;

    typedef SymbolsExprType symbols_expression_type;

    struct FunctorsVariadicExpr
    {
        struct Context
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<size_type, T1::value | T2::second_type::context >{};
                }
        };
        template<typename Funct>
        struct HasTestFunction
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || T2::second_type::template HasTestFunction<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTrialFunction
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || T2::second_type::template HasTrialFunction<Funct>::result >{};
                }
        };
    };

    static const size_type context =  std::decay_t<decltype( hana::fold( symbols_expression_type{}, hana::integral_constant<size_type,vm::POINT|vm::JACOBIAN|vm::KB|vm::NORMAL>{}, typename FunctorsVariadicExpr::Context{} ) )>::value;
    static const bool is_terminal = false;
    static const uint16_type imorder = Order;
    static const bool imIsPoly = false;

    template<typename Funct>
    struct HasTestFunction
    {
        static const bool result =  std::decay_t<decltype( hana::fold( symbols_expression_type{},
                                                                       hana::integral_constant<bool,false>{},
                                                                       typename FunctorsVariadicExpr::template HasTestFunction<Funct>{} ) )>::value;
    };
    template<typename Funct>
    struct HasTrialFunction
    {
        static const bool result =  std::decay_t<decltype( hana::fold( symbols_expression_type{},
                                                                       hana::integral_constant<bool,false>{},
                                                                       typename FunctorsVariadicExpr::template HasTrialFunction<Funct>{} ) )>::value;
    };

    template<typename Func>
    static const bool has_test_basis = false;
    template<typename Func>
    static const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    typedef GiNaC::ex expression_type;
    typedef GinacMatrix<M,N,Order,SymbolsExprType> this_type;
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
                          std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr=Environment::exprRepository() )
        :
        super( syms ),
        M_fun( fun.evalm() ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(),
        M_exprDesc( exprDesc )
        {
            std::string filenameExpanded = Environment::expand( filename );
            M_filename = (filenameExpanded.empty() || fs::path(filenameExpanded).is_absolute())? filenameExpanded : (fs::path(Environment::exprRepository())/filenameExpanded).string();

            DVLOG(2) << "Ginac matrix matrix constructor with expression_type \n";
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            // get filename if not given
            if ( M_filename.empty() && !M_exprDesc.empty() )
            {
                M_filename = Feel::vf::detail::ginacGetDefaultFileName( M_exprDesc, dirLibExpr );
            }

            // build ginac lib and link if necessary
            Feel::vf::detail::ginacBuildLibrary( exprs, syml, M_exprDesc, M_filename, world, M_cfun );
        }
    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms, std::string const& exprDesc,
                          std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr=Environment::exprRepository() )
        :
        super(syms),
        M_fun(fun.evalm()),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(),
        M_exprDesc( exprDesc )
        {
            std::string filenameExpanded = Environment::expand( filename );
            M_filename = (filenameExpanded.empty() || fs::path(filenameExpanded).is_absolute())? filenameExpanded : (fs::path(Environment::exprRepository())/filenameExpanded).string();

            DVLOG(2) << "Ginac matrix ex constructor with expression_type \n";
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );

            // get filename if not given
            if ( M_filename.empty() && !M_exprDesc.empty() )
            {
                M_filename = Feel::vf::detail::ginacGetDefaultFileName( M_exprDesc, dirLibExpr );
            }

            // build ginac lib and link if necessary
            Feel::vf::detail::ginacBuildLibrary( exprs, syml, M_exprDesc, M_filename, world, M_cfun );
        }

    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms,
                          GiNaC::FUNCP_CUBA const& cfun,  std::string const& exprDesc, symbols_expression_type const& expr )
        :
        super(syms),
        M_fun(fun.evalm()),
        M_cfun( new GiNaC::FUNCP_CUBA( cfun ) ),
        M_exprDesc( exprDesc ),
        M_expr( expr )
        {}

    GinacMatrix( GinacMatrix && fun ) = default;
    GinacMatrix( GinacMatrix const & fun ) = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    this_type& operator=( this_type const& ) = default;
    this_type& operator=( this_type && ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    symbols_expression_type const& symbolsExpression() const { return M_expr; }

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

    uint16_type index( std::string const& sname ) const
    {
        auto it = std::find_if( M_syms.begin(), M_syms.end(),
                                [=]( GiNaC::symbol const& s ) { return s.get_name() == sname; } );
        if ( it != M_syms.end() )
        {
            return it-M_syms.begin();
        }
        return invalid_uint16_type_value;
    }

    const std::vector<uint16_type> indices() const
    {
        std::vector<uint16_type> indices_vec;
        hana::for_each( M_expr, [&]( auto const& e )
                        {
                            indices_vec.push_back( this->index( e.first ) );
                        });
        return indices_vec;
    }

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
        struct TransformExprToTensor
        {
            template <typename T>
            struct apply {
                using type = typename T::second_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
            };

            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    return typename TransformExprToTensor::apply<T>::type( t.second,Geo_t{} );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) const
                {
                    return typename TransformExprToTensor::apply<T>::type( t.second,geom,fev,feu );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom, Basis_i_t const& fev ) const
                {
                    return typename TransformExprToTensor::apply<T>::type( t.second,geom,fev );
                }
            template <typename T>
            constexpr auto operator()(T const& t, Geo_t const& geom ) const
                {
                    return typename TransformExprToTensor::apply<T>::type( t.second,geom );
                }
        };

        using tuple_tensor_expr_type = std::decay_t<decltype( hana::transform( symbols_expression_type{}, TransformExprToTensor{} ) ) >;

        //typedef typename expression_type::value_type value_type;
        typedef double value_type;

        using key_type = key_t<Geo_t>;
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
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() ),
            M_t_expr( hana::transform( expr.symbolsExpression(), [&geom,&fev,&feu](auto const& t) { return TransformExprToTensor{}(t,geom,fev,feu); } ) ),
            M_t_expr_index( expr.indices() )
            {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() ),
            M_t_expr( hana::transform( expr.symbolsExpression(), [&geom,&fev](auto const& t) { return TransformExprToTensor{}(t,geom,fev); } ) ),
            M_t_expr_index( expr.indices() )
            {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), mat_type::Zero() ),
            M_x( expr.parameterValue() ),
            M_t_expr( hana::transform( expr.symbolsExpression(), [&geom](auto const& t) { return TransformExprToTensor{}(t,geom); } ) ),
            M_t_expr_index( expr.indices() )
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
                uint16_type k=0;
                hana::for_each( M_t_expr, [&k,&geom,this]( auto & e )
                                {
                                    if ( M_t_expr_index[k] != invalid_uint16_type_value )
                                        e.update( geom );
                                    ++k;
                                });

                M_gmc =  fusion::at_key<key_type>( geom ).get();

                int no = M*N;
                int ni = M_nsyms;//gmc_type::nDim;
                for(int q = 0; q < M_gmc->nPoints();++q )
                {
                    for ( auto const& comp : M_expr.indexSymbolXYZ() )
                        M_x[comp.second] = M_gmc->xReal( q )[comp.first];
                    // is it called for updates on faces? need to check that...
                    for ( auto const& comp : M_expr.indexSymbolN() )
                        M_x[comp.second] = M_gmc->unitNormal( q )[comp.first-3];
                    uint16_type k=0;
                    hana::for_each( M_t_expr, [&k,&q,this]( auto const& e )
                                    {
                                        uint16_type idx = M_t_expr_index[k];
                                        if ( idx != invalid_uint16_type_value )
                                            M_x[idx] = e.evalq( 0, 0, q );
                                        ++k;
                                    });
                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
                }

            }

        void update( Geo_t const& geom, uint16_type face )
            {
                uint16_type k=0;
                hana::for_each( M_t_expr, [&k,&geom,&face,this]( auto & e )
                                {
                                    if ( M_t_expr_index[k] != invalid_uint16_type_value )
                                        e.update( geom, face );
                                    ++k;
                                });

                M_gmc =  fusion::at_key<key_type>( geom ).get();

                int no = M*N;
                int ni = M_nsyms;//gmc_type::nDim;
                for(int q = 0; q < M_gmc->nPoints();++q )
                {
                    for ( auto const& comp : M_expr.indexSymbolXYZ() )
                        M_x[comp.second] = M_gmc->xReal( q )[comp.first];
                    for ( auto const& comp : M_expr.indexSymbolN() )
                        M_x[comp.second] = M_gmc->unitNormal( q )[comp.first-3];
                    uint16_type k=0;
                    hana::for_each( M_t_expr, [&k,&q,this]( auto const& e )
                                    {
                                        uint16_type idx = M_t_expr_index[k];
                                        if ( idx != invalid_uint16_type_value )
                                            M_x[idx] = e.evalq( 0, 0, q );
                                        ++k;
                                    });
                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
                }
            }

        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
            {
                boost::fusion::vector<CTX...> ctxvec( ctx... );
                update( boost::fusion::at_c<0>( ctxvec )->gmContext() );
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
        tuple_tensor_expr_type M_t_expr;
        const std::vector<uint16_type> M_t_expr_index;
    };

private:
    mutable expression_type  M_fun;
    boost::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
    std::string M_filename;
    std::string M_exprDesc;
    symbols_expression_type M_expr;
}; // GinacMatrix

template<int M,int N,int Order>
FEELPP_EXPORT std::ostream&
operator<<( std::ostream& os, GinacMatrix<M,N,Order> const& e )
{
    os << e.expression();
    return os;
}

template<int M, int N, int Order>
FEELPP_EXPORT std::string
str( GinacMatrix<M,N,Order> && e )
{
    return str( std::forward<GinacMatrix<M,N,Order>>(e).expression() );
}

template<int M, int N, int Order>
FEELPP_EXPORT std::string
str( GinacMatrix<M,N,Order> const& e )
{
    return str( e.expression() );
}

/// \endcond
/// \endcond
}} // Feel

#endif /* FEELPP_DETAIL_GINACMATRIX_HPP */
