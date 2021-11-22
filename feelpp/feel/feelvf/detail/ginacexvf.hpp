/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-02-14

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
#if !defined( FEELPP_VF_DETAIL_GINACEXVF_HPP )
#define FEELPP_VF_DETAIL_GINACEXVF_HPP 1

#include <feel/feelvf/detail/ginacmatrix.hpp>

#if 0
#include <any>

namespace Feel{
namespace vf{

/**
 * \class Ginac
 * \brief allow runtime ginac in expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order=2, typename SymbolsExprType = symbols_expression_empty_t >
class FEELPP_EXPORT GinacExVF : public Feel::vf::GiNaCBase
{
public:


    /** @name Typedefs
     */
    //@{
    typedef Feel::vf::GiNaCBase super;

    typedef SymbolsExprType symbols_expression_type;
    typedef typename symbols_expression_type::tuple_type symbols_expression_tuple_type;
    static const int nSymbolsExpr = std::decay_t<decltype(hana::size(symbols_expression_tuple_type{}))>::value;

    struct FunctorsVariadicExpr
    {
        template<typename Funct>
        struct HasTestFunction
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || std::tuple_element<1,typename T2::value_type>::type::template HasTestFunction<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTrialFunction
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || std::tuple_element<1,typename T2::value_type>::type::template HasTrialFunction<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTestBasis
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || std::tuple_element<1,typename T2::value_type>::type::template has_test_basis<Funct>::result >{};
                }
        };
        template<typename Funct>
        struct HasTrialBasis
        {
            template <typename T1,typename T2>
            constexpr auto operator()( T1 const& res,T2 const& e ) const
                {
                    return hana::integral_constant<bool, T1::value || std::tuple_element<1,typename T2::value_type>::type::template has_trial_basis<Funct>::result >{};
                }
        };

    };


    struct TransformSymbolsExprTupleToAny
    {
        template <typename T>
        constexpr auto operator()(T const& t) const
            {
                return std::vector<std::any>{};
            }
    };
    using tuple_expand_symbols_expr_type = std::decay_t<decltype( hana::transform( symbols_expression_tuple_type{}, TransformSymbolsExprTupleToAny{} ) ) >;

    static const size_type context = vm::DYNAMIC;

    static const bool is_terminal = false;

    template<typename Funct>
    struct HasTestFunction
    {
        static const bool result =  std::decay_t<decltype( hana::fold( symbols_expression_tuple_type{},
                                                                       hana::integral_constant<bool,false>{},
                                                                       typename FunctorsVariadicExpr::template HasTestFunction<Funct>{} ) )>::value;
    };
    template<typename Funct>
    struct HasTrialFunction
    {
        static const bool result =  std::decay_t<decltype( hana::fold( symbols_expression_tuple_type{},
                                                                       hana::integral_constant<bool,false>{},
                                                                       typename FunctorsVariadicExpr::template HasTrialFunction<Funct>{} ) )>::value;
    };

    template<typename Funct>
    static const bool has_test_basis = std::decay_t<decltype( hana::fold( symbols_expression_tuple_type{},
                                                                          hana::integral_constant<bool,false>{},
                                                                          typename FunctorsVariadicExpr::template HasTestBasis<Funct>{} ) )>::value;
    template<typename Funct>
    static const bool has_trial_basis = std::decay_t<decltype( hana::fold( symbols_expression_tuple_type{},
                                                                           hana::integral_constant<bool,false>{},
                                                                           typename FunctorsVariadicExpr::template HasTrialBasis<Funct>{} ) )>::value;
    using test_basis = std::nullptr_t;//TODO//typename expression_type::test_basis;
    using trial_basis = std::nullptr_t;//TODO//typename expression_type::trial_basis;

    typedef GiNaC::ex ginac_expression_type;
    typedef GinacExVF<Order,SymbolsExprType> this_type;
    typedef double value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;

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
    GinacExVF() : super(){}

    explicit GinacExVF( value_type value )
        :
        super(),
        M_fun( value ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_exprDesc( std::to_string( value ) ),
        M_isPolynomial( true ),
        M_polynomialOrder( 0 ),
        M_numericValue( value )
        {
            M_isNumericExpression = true ;
        }

    explicit GinacExVF( ginac_expression_type const & fun,
                        std::vector<GiNaC::symbol> const& syms,
                        std::string const& exprDesc,
                        std::string filename="",
                        WorldComm const& world=Environment::worldComm(),
                        std::string const& dirLibExpr=Environment::exprRepository(),
                        symbols_expression_type const& expr = symbols_expression_type() )
        :
        super( syms ),
        M_fun( fun ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(),
        M_exprDesc( exprDesc ),
        M_expr( expr ),
        M_isPolynomial( false ),
        M_polynomialOrder( Order ),
        M_numericValue( 0 )
        {
            this->updateNumericExpression();

            if ( !M_isNumericExpression )
            {
                std::string filenameExpanded = Environment::expand( filename );
                M_filename = (filenameExpanded.empty() || fs::path(filenameExpanded).is_absolute())? filenameExpanded : (fs::path(Environment::exprRepository())/filenameExpanded).string();

                DVLOG(2) << "Ginac constructor with expression_type \n";
                GiNaC::lst exprs({fun});
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

            this->updateForUse();
        }

    explicit GinacExVF( ginac_expression_type const & fun,
                        std::vector<GiNaC::symbol> const& syms,
                        GiNaC::FUNCP_CUBA const& cfun,
                        std::string const& exprDesc,
                        symbols_expression_type const& expr )
        :
        super( syms ),
        M_fun( fun ),
        M_cfun( new GiNaC::FUNCP_CUBA( cfun ) ),
        M_exprDesc( exprDesc ),
        M_expr( expr ),
        M_isPolynomial( false ),
        M_polynomialOrder( Order ),
        M_numericValue( 0 )
        {
            this->updateNumericExpression();
            this->updateForUse();
        }

    GinacExVF( GinacExVF const & fun ) = default;
    GinacExVF( GinacExVF && fun ) = default;

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

    //! polynomial order
    uint16_type polynomialOrder() const { return M_polynomialOrder; }

    //! expression is polynomial?
    bool isPolynomial() const { return M_isPolynomial; }

    ginac_expression_type const& expression() const
    {
        return M_fun;
    }

    const GiNaC::FUNCP_CUBA& fun() const
    {
        return *M_cfun;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! return true if the expression can be evaluated (TODO : iterate over symbols expression)
    bool isEvaluable() const
    {
        return M_isNumericExpression || ( M_indexSymbolXYZ.empty() && M_indexSymbolN.empty() && (M_syms.size() == M_symbolNameToValue.size()) );
    }
    bool isConstant() const { return this->isEvaluable(); }

    void setNumericValue( double val ) { M_isNumericExpression = true; M_numericValue = val; }


    const std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type> > > indices() const
    {
        std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type> > > indices_vec;

        hana::for_each( M_expr.tuple(), [this,&indices_vec]( auto const& evec )
                        {
                            for ( auto const& e : evec )
                            {
                                std::vector<std::tuple<uint16_type,uint16_type,uint16_type> > tmp;
                                if ( std::get<2>( e ).empty() ) // no suffix comp
                                {
                                    uint16_type idx = this->index( std::get<0>( e ) );
                                    if ( idx != invalid_v<uint16_type> )
                                        tmp.push_back( std::make_tuple( idx, 0, 0 ) );
                                }
                                else
                                {
                                    for ( auto const& [_suffix,compArray] : std::get<2>( e ) )
                                    {
                                        uint16_type idx = this->index( std::get<0>( e ) + _suffix );
                                        if ( idx != invalid_v<uint16_type> )
                                        {
                                            uint16_type c1 = compArray[0];
                                            uint16_type c2 = compArray[1];
                                            tmp.push_back( std::make_tuple( idx, c1, c2 ) );
                                        }
                                    }
                                }
                                indices_vec.push_back( tmp );
                            }
                        });
        return indices_vec;
    }
    bool isZero() const { return M_isNumericExpression? (M_numericValue == 0) : M_fun.is_zero(); }
    std::vector<GiNaC::symbol> const& syms() const { return M_syms; }

    std::string const& exprDesc() const { return M_exprDesc; }

    symbols_expression_type const& symbolsExpression() const { return M_expr; }
    symbols_expression_type & symbolsExpression() { return M_expr; }

    tuple_expand_symbols_expr_type const& expandSymbolsExpression() const { return M_expandSymbolsExpr; }
    tuple_expand_symbols_expr_type & expandSymbolsExpression() { return M_expandSymbolsExpr; }

    void setParameterValues( std::map<std::string, value_type> const& mp ) override
    {
        if ( M_isNumericExpression )
            return;

        super::setParameterValues( mp );

        auto exprIndex = this->indices();
        uint16_type k=0;
        hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&mp,&k,&exprIndex]( auto seId )
                        {
                            auto & evec = hana::at( this->symbolsExpression().tuple(), hana::int_c<seId> );
                            auto & evecExpand = hana::at( this->expandSymbolsExpression(), hana::int_c<seId> );
                            int nSubExpr = evec.size();
                            DCHECK( evecExpand.size() == nSubExpr ) << "something wrong";
                            for (int l=0;l<nSubExpr;++l,++k)
                            {
                                if ( exprIndex[k].empty() )
                                    continue;
                                auto & e = evec[l];
                                auto & theexprBase = std::get<1>( e );
                                auto & theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( this->symbolsExpression() ))>&>(evecExpand[l]);
                                theexprBase.setParameterValues( mp );
                                theexpr.setParameterValues( mp );
                            }
                        });
    }

    void updateParameterValues( std::map<std::string,double> & pv ) const
    {
        for ( auto const& [param,val] : this->symbolNameToValue() )
            pv[param] = val;

        auto exprIndex = this->indices();
        uint16_type k=0;
        hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&pv,&k,&exprIndex]( auto seId )
                        {
                            auto const& evec = hana::at( this->symbolsExpression().tuple(), hana::int_c<seId> );
                            auto const& evecExpand = hana::at( this->expandSymbolsExpression(), hana::int_c<seId> );
                            int nSubExpr = evec.size();
                            DCHECK( evecExpand.size() == nSubExpr ) << "something wrong";
                            for (int l=0;l<nSubExpr;++l,++k)
                            {
                                if ( exprIndex[k].empty() )
                                    continue;
                                auto const& e = evec[l];
                                auto const& theexprBase = std::get<1>( e );
                                auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( this->symbolsExpression() ))>const&>(evecExpand[l]);
                                theexpr.updateParameterValues( pv );
                            }
                        });
    }

    template <typename ExpandSymbolsExprType>
    auto applySymbolsExpr( ExpandSymbolsExprType const& se ) const
    {
        auto newSymbolExpr = symbolsExpr( this->symbolsExpression(), se );
        using new_expr_type = GinacExVF<Order,std::decay_t<decltype(newSymbolExpr)>>;
        new_expr_type res( this->expression(), this->symbols(), this->fun(), this->exprDesc(), newSymbolExpr );
        res.setParameterValues( this->symbolNameToValue() );
        return res;
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        if constexpr ( nSymbolsExpr > 0 )
            return this->hasSymbolDependencyImpl( symb, Feel::vf::symbolsExpr( this->symbolsExpression(), se ) );
        else
            return this->hasSymbolDependencyImpl( symb, se );

    }

    template <typename TheSymbolExprType>
        void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
    {
        if constexpr ( nSymbolsExpr > 0 )
            this->dependentSymbolsImpl( symb, res, Feel::vf::symbolsExpr( this->symbolsExpression(), se ) );
        else
            this->dependentSymbolsImpl( symb, res, se );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto newse = Feel::vf::symbolsExpr( this->symbolsExpression(), se );

        std::map<std::string,std::set<std::string>> symbolToApplyDiff;
        this->dependentSymbolsImpl( diffVariable, symbolToApplyDiff, newse );

        auto tupleDiffSymbolsExpr = hana::transform( this->symbolsExpression().tuple(), [this,&diffVariable,&symbolToApplyDiff,&newse,&world,&dirLibExpr](auto const& evec)
                                                     {
                                                         int nSubExpr = evec.size();
                                                         using _expr_type = std::decay_t<decltype( std::get<1>( evec.front() ) )>;
                                                         using _expr_diff_type =  std::decay_t<decltype(  std::get<1>( evec.front() ).template diff<diffOrder>( diffVariable,world,dirLibExpr,newse ) )>;
                                                         symbol_expression_t<_expr_diff_type> seDiffExpr;
                                                         for (int l=0;l<nSubExpr;++l)
                                                         {
                                                             auto const& e = evec[l];
                                                             std::string const& currentSymbName = std::get<0>( e );

                                                             auto itFindSymb = symbolToApplyDiff.find( currentSymbName );
                                                             if ( itFindSymb == symbolToApplyDiff.end() ) // not depend on diffVariable
                                                                 continue;
                                                             auto const& theexpr = std::get<1>( e );

                                                             auto currentDiffExpr = theexpr.template diff<diffOrder>( diffVariable,world,dirLibExpr,newse );
                                                             std::string currentDiffSymbName = (boost::format( "diff_%1%_%2%_%3%" )%currentSymbName %diffVariable %diffOrder ).str();
                                                             seDiffExpr.add( currentDiffSymbName, currentDiffExpr, std::get<2>( e ) );
                                                         }
                                                         return seDiffExpr;
                                                     });


        std::vector<GiNaC::symbol> resSymbol = this->symbols();
        ginac_expression_type res;
        if ( super::hasSymbol( diffVariable ) )
        {
            GiNaC::symbol diffSymb = Feel::symbols( std::vector<std::string>( { diffVariable } ), this->symbols() )[0];
            res += this->expression().diff( diffSymb, diffOrder );
        }

        for ( auto const& [currentSymbName,currentSymbSuffixes] : symbolToApplyDiff )
        {
            std::string currentDiffSymbName = (boost::format( "diff_%1%_%2%_%3%" )%currentSymbName %diffVariable %diffOrder ).str();
            for ( std::string const& currentSymbSuffix : currentSymbSuffixes )
            {
                std::string currentSymbNameWithSuffix = currentSymbName + currentSymbSuffix;
                auto itSym = std::find_if( this->symbols().begin(), this->symbols().end(),
                                           [&currentSymbNameWithSuffix]( GiNaC::symbol const& s ) { return s.get_name() == currentSymbNameWithSuffix; } );
                if ( itSym == this->symbols().end() )
                    continue;
                GiNaC::symbol const& currentSymb = *itSym;
                resSymbol.push_back(  GiNaC::symbol(currentDiffSymbName+currentSymbSuffix) );
                res += resSymbol.back()*this->expression().diff( currentSymb, diffOrder );
            }
        }

        auto seWithDiff = Feel::vf::symbolsExpr( this->symbolsExpression(), SymbolsExpr( std::move( tupleDiffSymbolsExpr ) ) );
        using symbols_expression_with_diff_type = std::decay_t<decltype( seWithDiff )>;
        using _expr_type = GinacExVF<Order,symbols_expression_with_diff_type>;
#if 0
        std::string exprDesc = str( res );
#else
        std::string exprDesc = (boost::format("diff(%1%)_%2%_o%3%")% this->exprDesc() %diffVariable %diffOrder ).str();
#endif
        _expr_type resExpr( res, resSymbol, exprDesc, ""/*filename*/, world, dirLibExpr, seWithDiff );

        std::map<std::string,double> pv;
        this->updateParameterValues( pv );
        resExpr.setParameterValues( pv );

        return resExpr;
    }

    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        template <typename T>
        struct TransformExprToTensor
        {
            using type = typename T::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        };

        using tuple_tensor_expr_type = std::decay_t<decltype( hana::transform( symbols_expression_tuple_type{}, TransformSymbolsExprTupleToAny{} ) ) >;

        typedef double value_type;

        using key_type = key_t<Geo_t>;
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
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_is_zero( expr.isZero() ),
            M_is_constant( expr.isConstant() ),
            M_t_expr_index( expr.indices() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( vec_type::Zero(M_gmc->nPoints()) ),
            M_x( expr.parameterValue() ),
            M_yConstant( (M_is_constant)? expr.evaluate()(0,0) : value_type(0) )
            {
                this->initSubTensor( geom, fev, feu );
            }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_is_zero( expr.isZero() ),
            M_is_constant( expr.isConstant() ),
            M_t_expr_index( expr.indices() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( vec_type::Zero(M_gmc->nPoints()) ),
            M_x( expr.parameterValue() ),
            M_yConstant( (M_is_constant)? expr.evaluate()(0,0) : value_type(0) )
            {
                this->initSubTensor( geom, fev );
            }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_is_zero( expr.isZero() ),
            M_is_constant( expr.isConstant() ),
            M_t_expr_index( expr.indices() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( vec_type::Zero(M_gmc->nPoints()) ),
            M_x( expr.parameterValue() ),
            M_yConstant( (M_is_constant)? expr.evaluate()(0,0) : value_type(0) )
            {
                this->initSubTensor( geom );
            }

        template<typename IM>
        void init( IM const& im ) {}

        FEELPP_DONT_INLINE void updateFun(Geo_t const& geom )
            {
                M_gmc =  fusion::at_key<key_type>( geom ).get();

                int no = 1;
                int ni = M_nsyms;///gmc_type::nDim;

                for(int q = 0; q < M_gmc->nPoints();++q )
                {
                    for ( auto const& comp : M_expr.indexSymbolXYZ() )
                        M_x[comp.second] = M_gmc->xReal( q )[comp.first];
                    // is it called for updates on faces? need to check that...
                    for ( auto const& comp : M_expr.indexSymbolN() )
                        M_x[comp.second] = M_gmc->unitNormal( q )[comp.first-3];

                    uint16_type k=0;
                    hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&geom,&k,&q]( auto seId ) {
                            auto const& evec = hana::at( M_expr.symbolsExpression().tupleExpr, hana::int_c<seId> );
                            int nSubExpr = evec.size();
                            auto & evecTensorExpr = hana::at(M_t_expr, hana::int_c<seId> );
                            for (int l=0;l<nSubExpr;++l,++k)
                            {
                                if ( M_t_expr_index[k].empty() )
                                    continue;
                                auto const& e = evec[l];
                                auto const& theexprBase = std::get<1>( e );
                                using subtensor_type = typename TransformExprToTensor<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr.symbolsExpression() ))>>::type;
                                auto & subTensor = std::any_cast<subtensor_type&>(evecTensorExpr[l]);
                                for ( auto const& [idx,c1,c2] : M_t_expr_index[k] )
                                    M_x[idx] = subTensor.evalq( c1, c2, q );
                            }
                        });

                    M_fun(&ni,M_x.data(),&no,&M_y[q]);
                }
            }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            {
                this->updateImpl( geom, fev, feu );
            }
        void update( Geo_t const& geom, Basis_i_t const& fev )
            {
                this->updateImpl( geom, fev );
            }
        void update( Geo_t const& geom )
            {
                this->updateImpl( geom );
            }
        void update( Geo_t const& geom, uint16_type face )
            {
                this->updateImpl( geom, face );
            }

        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
            {
                if ( M_is_zero ) return;
                boost::fusion::vector<CTX...> ctxvec( ctx... );
                update( boost::fusion::at_c<0>( ctxvec )->gmContext() );
            }

#if 0
        value_type
        evalij( uint16_type i, uint16_type j ) const
            {
                return 0;
            }
#endif
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return evalq( c1,c2,q );
            }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return evalq( c1,c2,q );
            }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return ( M_is_constant )? M_yConstant : M_y[q];
            }

    private :

        template<typename... TheArgsType>
        void initSubTensor( const TheArgsType&... theInitArgs )
            {
                if ( M_is_constant ) return;

                uint16_type k=0;
                hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&theInitArgs...,&k]( auto seId ) {
                        auto const& evec = hana::at( M_expr.symbolsExpression().tupleExpr, hana::int_c<seId> );
                        auto const & evecExpand = hana::at( M_expr.expandSymbolsExpression(), hana::int_c<seId> );
                        int nSubExpr = evec.size();
                        DCHECK( evecExpand.size() == nSubExpr ) << "something wrong";
                        auto & evecTensorExpr = hana::at(M_t_expr, hana::int_c<seId> );
                        evecTensorExpr.resize( nSubExpr );
                        for (int l=0;l<nSubExpr;++l,++k)
                        {
                            if ( M_t_expr_index[k].empty() )
                                continue;

                            auto const& e = evec[l];
                            auto const& theexprBase = std::get<1>( e );
                            auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr.symbolsExpression() ))> const&>(evecExpand[l]);
                            evecTensorExpr[l] = typename TransformExprToTensor<std::decay_t<decltype(theexpr)>>::type( theexpr,theInitArgs... );
                        }
                    });
            }


        template<typename... TheArgsType>
        void updateImpl( Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                if ( M_is_zero ) return;
                if ( M_is_constant ) return;

                uint16_type k=0;
                hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&geom,&theUpdateArgs...,&k]( auto seId ) {
                        auto const& evec = hana::at( M_expr.symbolsExpression().tupleExpr, hana::int_c<seId> );
                        int nSubExpr = evec.size();
                        auto & evecTensorExpr = hana::at(M_t_expr, hana::int_c<seId> );
                        for (int l=0;l<nSubExpr;++l,++k)
                        {
                            if ( M_t_expr_index[k].empty() )
                                continue;
                            auto const& e = evec[l];
                            auto const& theexprBase = std::get<1>( e );
                            using subtensor_type = typename TransformExprToTensor<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr.symbolsExpression() ))>>::type;
                            auto & subTensor = std::any_cast<subtensor_type&>(evecTensorExpr[l]);
                            subTensor.update( geom,theUpdateArgs... );
                        }
                    });

                updateFun( geom );
            }

    private :
        this_type const& M_expr;
        GiNaC::FUNCP_CUBA M_fun;
        const bool M_is_zero;
        const bool M_is_constant;
        tuple_tensor_expr_type M_t_expr;
        const std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> M_t_expr_index; // (id,c1,c2)
        gmc_ptrtype M_gmc;

        int M_nsyms;
        vec_type M_y;
        vec_type M_x;
        const value_type M_yConstant;
    };

    evaluate_type
    evaluate( std::map<std::string,value_type> const& mp  )
    {
        this->setParameterValues( mp );
        return this->evaluateImpl( true, Environment::worldCommPtr() );
    }

    evaluate_type
    evaluate( bool parallel = true, worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr() ) const
    {
        return this->evaluateImpl( parallel, worldcomm );
    }
private :

    void updateNumericExpression()
    {
        auto resToNum = toNumericValues( M_fun );
        CHECK( resToNum.size() == 1 )  << "invalid size " << resToNum.size() << " : must be 1";
        M_isNumericExpression = resToNum[0].first;
        if ( M_isNumericExpression )
        {
            M_numericValue = resToNum[0].second;
            M_isPolynomial = 1;
            M_polynomialOrder = 0;
        }
    }
    void updateForUse()
    {
        if ( M_isNumericExpression )
            return;

        std::vector<std::pair<GiNaC::symbol,int>> symbTotalDegree;
        for ( auto const& thesymbxyz : this->indexSymbolXYZ() )
            symbTotalDegree.push_back( std::make_pair( M_syms[thesymbxyz.second], 1 ) );
        bool symbExprArePolynomials = true;
        hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&symbTotalDegree,&symbExprArePolynomials]( auto seId )
                        {
                            auto const& evec = hana::at( M_expr.tupleExpr, hana::int_c<seId> );
                            int nSubExpr = evec.size();
                            auto & evecExpand = hana::at(M_expandSymbolsExpr, hana::int_c<seId> );
                            evecExpand.resize( nSubExpr );
                            for (int k=0;k<nSubExpr;++k)
                            {
                                auto const& e = evec[k];

                                if ( std::get<2>( e ).empty() )
                                {
                                    uint16_type idx = this->index( std::get<0>( e ) );
                                    if ( idx == invalid_uint16_type_value )
                                        continue;

                                    // call an optional update function
                                    if ( std::get<3>( e ) )
                                        std::get<3>( e )();

                                    auto const& theexprBase = std::get<1>( e );
                                    auto theexpr = theexprBase.applySymbolsExpr( M_expr );
                                    evecExpand[k] = theexpr;

                                    M_context = M_context | Feel::vf::dynamicContext( theexpr );

                                    if ( theexpr.isPolynomial() )
                                        symbTotalDegree.push_back( std::make_pair( M_syms[idx], theexpr.polynomialOrder() ) );
                                    else
                                        symbExprArePolynomials = false;
                                }
                                else
                                {
                                    if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                         continue;

                                    // call an optional update function
                                    if ( std::get<3>( e ) )
                                        std::get<3>( e )();

                                    auto const& theexprBase = std::get<1>( e );
                                    auto theexpr = theexprBase.applySymbolsExpr( M_expr );
                                    evecExpand[k] = theexpr;

                                    M_context = M_context | Feel::vf::dynamicContext( theexpr );

                                    for ( auto const& [_suffix,compArray] : std::get<2>( e ) )
                                    {
                                        uint16_type idx = this->index( std::get<0>( e ) + _suffix );
                                        if ( idx != invalid_v<uint16_type> )
                                        {
                                            if ( theexpr.isPolynomial() )
                                                symbTotalDegree.push_back( std::make_pair( M_syms[idx], theexpr.polynomialOrder() ) );
                                            else
                                                symbExprArePolynomials = false;
                                        }
                                    }
                                }
                            }
                        });

        if ( symbExprArePolynomials )
        {
            GiNaC::lst symlIsPoly;
            for ( auto const& thesymb : symbTotalDegree )
                symlIsPoly.append( thesymb.first );
            M_isPolynomial = M_fun.is_polynomial( symlIsPoly );
        }
        else
            M_isPolynomial = false;

        if ( M_isPolynomial )
            M_polynomialOrder = totalDegree( M_fun, symbTotalDegree );
        else
            M_polynomialOrder = Order;
    }

    evaluate_type
    evaluateImpl( bool parallel, worldcomm_ptr_t const& worldcomm ) const
    {
        if ( M_isNumericExpression )
            return evaluate_type::Constant( M_numericValue );
        int no = 1;
        int ni = M_syms.size();

        vec_type x( ni );
        for ( uint16_type k=0;k<ni;++k )
            x[k] = M_params[k];

        hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&x,&parallel,&worldcomm]( auto seId )
                        {
                            auto const& evec = hana::at( M_expr.tupleExpr, hana::int_c<seId> );
                            auto const& evecExpand = hana::at(M_expandSymbolsExpr, hana::int_c<seId> );
                            int nSubExpr = evec.size();
                            for (int k=0;k<nSubExpr;++k)
                            {
                                auto const& e = evec[k];

                                if ( std::get<2>( e ).empty() )
                                {
                                    uint16_type idx = this->index( std::get<0>( e ) );
                                    if ( idx == invalid_v<uint16_type> )
                                        continue;

                                    auto const& theexprBase = std::get<1>( e );
                                    auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr ))> const&>(evecExpand[k]);

                                    x[idx] = theexpr.evaluate( parallel, worldcomm )(0,0);
                                }
                                else
                                {
                                    if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                        continue;

                                    auto const& theexprBase = std::get<1>( e );
                                    auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr ))> const&>(evecExpand[k]);

                                    //auto const& theexpr = std::get<1>( e );
                                    for ( auto const& [_suffix,compArray] : std::get<2>( e ) )
                                    {
                                        uint16_type idx = this->index( std::get<0>( e ) + _suffix );
                                        if ( idx != invalid_v<uint16_type> )
                                        {
                                            uint16_type c1 = compArray[0];
                                            uint16_type c2 = compArray[1];
                                            x[idx] = theexpr.evaluate( parallel, worldcomm )(c1,c2);
                                        }
                                    }
                                }
                            }
                        });

        value_type res;
        (*M_cfun)(&ni,x.data(),&no,&res );

        return evaluate_type::Constant( res );
    }

    template <typename AnElementOfSymbolExprType>
    bool hasAtLeastOneSymbolDependency( AnElementOfSymbolExprType const& e ) const
    {
        bool res = false;
        if ( std::get<2>( e ).empty() )
        {
            res = this->index( std::get<0>( e ) ) != invalid_v<uint16_type>;
        }
        else
        {
            for ( auto const& [_suffix,compArray] : std::get<2>( e ) )
                if ( this->index( std::get<0>( e ) + _suffix ) != invalid_v<uint16_type> )
                {
                    res = true;
                    break;
                }
        }
        return res;
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependencyImpl( std::string const& symb, TheSymbolExprType const& se ) const
    {
        if ( super::hasSymbol( symb ) )
            return true;
        bool res = false;

        hana::for_each( se.tuple(), [this,&symb,&res,&se]( auto const& evec )
                        {
                            if ( res )
                                return;

                            int nSubExpr = evec.size();
                            for (int l=0;l<nSubExpr;++l)
                            {
                                auto const& e = evec[l];
                                std::string const& currentSymbName = std::get<0>( e );
                                // if ( !super::hasSymbol( currentSymbName ) )
                                //     continue;
                                if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                    continue;

                                auto const& currentExpr = std::get<1>( e );
                                res = currentExpr.hasSymbolDependency( symb, se );
                                if ( res )
                                    break;
                            }
                        });
        return res;
    }

    template <typename TheSymbolExprType>
        void dependentSymbolsImpl( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
    {
        hana::for_each( se.tuple(), [this,&symb,&res,&se]( auto const& evec )
                        {
                            int nSubExpr = evec.size();
                            for (int l=0;l<nSubExpr;++l)
                            {
                                auto const& e = evec[l];
                                std::string const& currentSymbName = std::get<0>( e );
                                if ( res.find( currentSymbName ) != res.end() )
                                    continue;

                                // if ( !super::hasSymbol( currentSymbName ) )
                                //     continue;
                                if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                    continue;

                                auto const& currentExpr = std::get<1>( e );
                                if ( !currentExpr.hasSymbolDependency( symb, se ) )
                                    continue;

                                if ( std::get<2>( e ).empty() )
                                    res[currentSymbName].insert( "" );
                                else
                                {
                                    for ( auto const& [_suffix,compArray] : std::get<2>( e ) )
                                        res[currentSymbName].insert( _suffix );
                                }

                                //if ( res.insert( std::get<0>( e ) ).second )
                                currentExpr.dependentSymbols( symb, res, se );
                            }
                        });
    }

private:
    mutable ginac_expression_type  M_fun;
    std::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
    std::string M_filename;
    std::string M_exprDesc;
    symbols_expression_type M_expr;
    tuple_expand_symbols_expr_type M_expandSymbolsExpr;
    bool M_isPolynomial;
    uint16_type M_polynomialOrder;
    bool M_isNumericExpression;
    value_type M_numericValue;
};

template<int Order,typename SymbolsExprType>
FEELPP_EXPORT std::ostream&
operator<<( std::ostream& os, GinacExVF<Order,SymbolsExprType> const& e )
{
    os << e.expression();
    return os;
}

template<int Order,typename SymbolsExprType>
FEELPP_EXPORT std::string
str( GinacExVF<Order,SymbolsExprType> && e )
{
    return str( std::forward<GinacExVF<Order,SymbolsExprType>>(e).expression() );
}
template<int Order,typename SymbolsExprType>
FEELPP_EXPORT std::string
str( GinacExVF<Order,SymbolsExprType> const& e )
{
    return str( e.expression() );
}

template<int Order = 2>
using GinacEx = GinacExVF<Order>;


}} // feel::vf
#endif

#endif
