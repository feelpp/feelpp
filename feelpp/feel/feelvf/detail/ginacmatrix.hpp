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

#include <any>

namespace Feel
{
namespace vf {

/**
 * Handle Ginac matrix expression
 */
template<int M=1, int N=1, int Order = 2, typename SymbolsExprType = symbols_expression_empty_t >
class FEELPP_EXPORT GinacMatrix : public Feel::vf::GiNaCBase
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
    using test_basis = std::nullptr_t; // TODO
    using trial_basis = std::nullptr_t; // TODO

    typedef GiNaC::ex expression_type;
    typedef GinacMatrix<M,N,Order,SymbolsExprType> this_type;
    typedef double value_type;

    //typedef Eigen::MatrixXd evaluate_type;
    // be careful that the matrix passed to ginac must be Row Major,
    // however if the number of columns is 1 then eigen3 fails with
    // an assertion, so we have a special when N=1 and have the
    // matrix column major which is ok in this case
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<N>, mpl::int_<1>>,
                              mpl::identity<Eigen::Matrix<value_type,M,1>>,
                              mpl::identity<Eigen::Matrix<value_type,M,N,Eigen::RowMajor>>>::type::type evaluate_type;



    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vec_type;

    template<typename... TheExpr>
    struct Lambda
    {
        using new_se_type = typename SymbolsExprType::template Lambda<TheExpr...>::type;
        using type = GinacMatrix<M,N,Order,new_se_type>;
    };
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  )
    {
        typename Lambda<TheExpr...>::type res( this->expression(), this->symbols(), this->fun(), this->exprDesc(), M_expr.applyLambda( e... ) );
        res.setParameterValues( this->symbolNameToValue() );
        return res;
    }

    //@}

    /** @name Constructors, destructor
     */
    //@{

    GinacMatrix() : super() {}

    explicit GinacMatrix( value_type value )
        :
        super(),
        M_fun( value ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_exprDesc( std::to_string( value ) ),
        M_isPolynomial( true ),
        M_polynomialOrder( 0 ),
        M_numericValue( evaluate_type::Constant( value ) )
        {
            M_isNumericExpression = true ;
        }

    explicit GinacMatrix( GiNaC::matrix const & fun, std::vector<GiNaC::symbol> const& syms, std::string const& exprDesc,
                          std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr=Environment::exprRepository(),
                          symbols_expression_type const& expr = symbols_expression_type() )
        :
        super( syms ),
        //M_fun( fun.evalm() ),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(),
        M_exprDesc( exprDesc ),
        M_expr( expr ),
        M_isPolynomial( false ),
        M_polynomialOrder( Order ),
        M_numericValue( evaluate_type::Zero() )
        {
            GiNaC::ex funWithEvalm = fun.evalm();
            if constexpr ( M*N == 1 )
            {
                CHECK( funWithEvalm.nops() == 1 ) << "expression " << funWithEvalm << "should have only one element, but " << funWithEvalm.nops();
                M_fun = funWithEvalm.op( 0 );
                if ( GiNaC::is_a<GiNaC::lst>( M_fun ) )
                    CHECK( false ) << "scalar expression should not be a lst";
            }
            else
                M_fun = funWithEvalm;

            this->updateNumericExpression();

            if ( !M_isNumericExpression )
            {
                std::string filenameExpanded = Environment::expand( filename );
                M_filename = (filenameExpanded.empty() || fs::path(filenameExpanded).is_absolute())? filenameExpanded : (fs::path(Environment::exprRepository())/filenameExpanded).string();

                DVLOG(2) << "Ginac matrix matrix constructor with expression_type \n";
                GiNaC::lst exprs;
                if constexpr ( M*N == 1 )
                {
                    exprs.append( M_fun );
                }
                else
                {
                    for( int i = 0; i < M_fun.nops(); ++i )
                        exprs.append( M_fun.op(i) );
                }

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
    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms, std::string const& exprDesc,
                          std::string filename="", WorldComm const& world=Environment::worldComm(), std::string const& dirLibExpr=Environment::exprRepository(),
                          symbols_expression_type const& expr = symbols_expression_type() )
        :
        super(syms),
        //M_fun(fun.evalm()),
        M_cfun( new GiNaC::FUNCP_CUBA() ),
        M_filename(),
        M_exprDesc( exprDesc ),
        M_expr( expr ),
        M_isPolynomial( false ),
        M_polynomialOrder( Order ),
        M_numericValue( evaluate_type::Zero() )
        {
            GiNaC::ex funWithEvalm = fun.evalm();
            if constexpr ( M*N == 1 )
            {
                if ( GiNaC::is_a<GiNaC::lst>( funWithEvalm ) )
                {
                    CHECK( funWithEvalm.nops() == 1 ) << "expression " << funWithEvalm << "should have only one element, but " << funWithEvalm.nops();
                    M_fun = funWithEvalm.op( 0 );
                    if ( GiNaC::is_a<GiNaC::lst>( M_fun ) )
                        CHECK( false ) << "scalar expression should not be a lst";
                }
                else
                    M_fun = funWithEvalm;
            }
            else
                M_fun = funWithEvalm;

            this->updateNumericExpression();

            if ( !M_isNumericExpression )
            {
                std::string filenameExpanded = Environment::expand( filename );
                M_filename = (filenameExpanded.empty() || fs::path(filenameExpanded).is_absolute())? filenameExpanded : (fs::path(Environment::exprRepository())/filenameExpanded).string();

                DVLOG(2) << "Ginac matrix ex constructor with expression_type \n";
                GiNaC::lst exprs;
                if constexpr ( M*N == 1 )
                {
                    exprs.append( M_fun );
                }
                else
                {
                    for( int i = 0; i < M_fun.nops(); ++i )
                        exprs.append( M_fun.op(i) );
                }

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

    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms,
                          GiNaC::FUNCP_CUBA const& cfun,  std::string const& exprDesc, symbols_expression_type const& expr )
        :
        super(syms),
        M_fun(fun.evalm()),
        M_cfun( new GiNaC::FUNCP_CUBA( cfun ) ),
        M_exprDesc( exprDesc ),
        M_expr( expr ),
        M_isPolynomial( false ),
        M_polynomialOrder( Order ),
        M_numericValue( evaluate_type::Zero() )
        {
            this->updateNumericExpression();
            this->updateForUse();
        }

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

    //! polynomial order
    uint16_type polynomialOrder() const { return M_polynomialOrder; }

    //! expression is polynomial?
    bool isPolynomial() const { return M_isPolynomial; }

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
                                auto & theexprBase = e.expr();
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
                                auto const& theexprBase = e.expr();
                                auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( this->symbolsExpression() ))>const&>(evecExpand[l]);
                                theexpr.updateParameterValues( pv );
                            }
                        });
    }

    void renameSymbols( std::map<std::string,std::string> const& old2new )
    {
        for ( auto const& [oldSymbName,newSymbName] : old2new )
        {
            uint16_type indexOldSymb = this->index( oldSymbName );
            if ( indexOldSymb == invalid_uint16_type_value )
                continue;

            GiNaC::symbol oldSymb = M_syms[indexOldSymb];
            M_syms[indexOldSymb] = GiNaC::symbol( newSymbName );
            M_fun = M_fun.subs( oldSymb == M_syms[indexOldSymb] );

            auto itFindSymbToValue = M_symbolNameToValue.find( oldSymbName );
            if ( itFindSymbToValue != M_symbolNameToValue.end() )
                M_symbolNameToValue[newSymbName] = itFindSymbToValue->second;
            else
                M_symbolNameToValue[newSymbName] = 0;
            M_symbolNameToValue.erase( oldSymbName );
        }
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

    template <typename ExpandSymbolsExprType>
    auto applySymbolsExpr( ExpandSymbolsExprType const& se ) const
    {
        auto newSymbolExpr = symbolsExpr( this->symbolsExpression(), se );
        using new_expr_type = GinacMatrix<M,N,Order,std::decay_t<decltype(newSymbolExpr)>>;
        new_expr_type res( this->expression(), this->symbols(), this->fun(), this->exprDesc(), newSymbolExpr );
        res.setParameterValues( this->symbolNameToValue() );
        return res;
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto newse = Feel::vf::symbolsExpr( this->symbolsExpression(), se );

        std::map<std::string, std::map<std::string,std::set<std::string>>> mapSymbolToApplyDiff; // diffVariable -> ( symbol to apply diff, (set of suffix) )
        this->dependentSymbolsImpl( diffVariable, mapSymbolToApplyDiff[diffVariable], newse );
        auto diff_se = this->diffSymbolsExprImpl<diffOrder>( mapSymbolToApplyDiff, world, dirLibExpr, newse );


        std::vector<GiNaC::symbol> resSymbol = this->symbols();
        std::vector<GiNaC::ex> res(M*N);
        //CHECK( GiNaC::is_a<GiNaC::lst>(this->expression()) ) << "something wrong";
        if ( super::hasSymbol( diffVariable ) )
        {
            GiNaC::symbol diffSymb = Feel::symbols( std::vector<std::string>( { diffVariable } ), this->symbols() )[0];
            if constexpr ( M*N==1 )
                res[0] += this->expression().diff( diffSymb, diffOrder );
            else
            {
                for( int i = 0; i < M*N; ++i )
                    res[i] += this->expression().op(i).diff( diffSymb, diffOrder );
            }
        }

        for ( auto const& [currentSymbName,currentSymbSuffixes] : mapSymbolToApplyDiff[diffVariable] )
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
                resSymbol.push_back( GiNaC::symbol(currentDiffSymbName+currentSymbSuffix) );
                if constexpr ( M*N==1 )
                    res[0] += resSymbol.back()*this->expression().diff( currentSymb, diffOrder );
                else
                {
                    for( int i = 0; i < M*N; ++i )
                        res[i] += resSymbol.back()*this->expression().op(i).diff( currentSymb, diffOrder );
                }
            }
        }

        auto seWithDiff = Feel::vf::symbolsExpr( this->symbolsExpression(), diff_se );
        using symbols_expression_with_diff_type = std::decay_t<decltype( seWithDiff )>;
        using _expr_type = GinacMatrix<M,N,Order,symbols_expression_with_diff_type>;
        GiNaC::matrix resmat(M,N,res);
#if 0
        std::string exprDesc = str( resmat );
#else
        std::string exprDesc = (boost::format("diff(%1%)_%2%_o%3%")% this->exprDesc() %diffVariable %diffOrder ).str();
#endif
        _expr_type resExpr( resmat, resSymbol, exprDesc, ""/*filename*/, world, dirLibExpr, seWithDiff );

        std::map<std::string,double> pv;
        this->updateParameterValues( pv );
        resExpr.setParameterValues( pv );

        return resExpr;
    }

    template <int Dim>
    auto grad( std::vector<std::string> const& diffVariables, WorldComm const& world, std::string const& dirLibExpr ) const
    {
        CHECK( N == 1 ) << "not implmented";
        CHECK( Dim == diffVariables.size() ) << "incompatible Dim with diffVariables : " << Dim << " vs " << diffVariables.size();

        auto const& newse = this->symbolsExpression();
        static const int diffOrder = 1;

        std::map<std::string, std::map<std::string,std::set<std::string>>> mapSymbolToApplyDiff; // diffVariable -> ( symbol to apply diff, (set of suffix) )
        for ( std::string const& diffVariable : diffVariables )
            this->dependentSymbolsImpl( diffVariable, mapSymbolToApplyDiff[diffVariable], newse );
        auto diff_se = this->diffSymbolsExprImpl<diffOrder>( mapSymbolToApplyDiff, world, dirLibExpr, newse );


        std::vector<GiNaC::symbol> resSymbol = this->symbols();
        std::vector<GiNaC::ex> res(M*N*Dim);
        for ( int d=0;d<Dim;++d )
        {
            std::string const& diffVariable = diffVariables[d];
            if ( super::hasSymbol( diffVariable ) )
            {
                GiNaC::symbol diffSymb = Feel::symbols( std::vector<std::string>( { diffVariable } ), this->symbols() )[0];
                if constexpr ( M*N==1 )
                    res[d] += this->expression().diff( diffSymb, diffOrder );
                else
                {
                    for( int i = 0; i < M*N; ++i )
                        res[i*Dim+d] += this->expression().op(i).diff( diffSymb, diffOrder );
                }
            }

            for ( auto const& [currentSymbName,currentSymbSuffixes] : mapSymbolToApplyDiff[diffVariable] )
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
                    resSymbol.push_back( GiNaC::symbol(currentDiffSymbName+currentSymbSuffix) );
                    if constexpr ( M*N==1 )
                                     res[d] += resSymbol.back()*this->expression().diff( currentSymb, diffOrder );
                    else
                    {
                        for( int i = 0; i < M*N; ++i )
                            res[i*Dim+d] += resSymbol.back()*this->expression().op(i).diff( currentSymb, diffOrder );
                    }
                }
            }
        }

        auto seWithDiff = Feel::vf::symbolsExpr( this->symbolsExpression(), diff_se );
        using symbols_expression_with_diff_type = std::decay_t<decltype( seWithDiff )>;
        using _expr_type = GinacMatrix<M,Dim,Order,symbols_expression_with_diff_type>;
        GiNaC::matrix resmat(M,Dim,res);
        std::string exprDesc = (boost::format("grad(%1%)")% this->exprDesc() ).str();
        for ( std::string const& diffVariable : diffVariables )
            exprDesc += "_" + diffVariable;
        _expr_type resExpr( resmat, resSymbol, exprDesc, ""/*filename*/, world, dirLibExpr, seWithDiff );

        std::map<std::string,double> pv;
        this->updateParameterValues( pv );
        resExpr.setParameterValues( pv );

        return resExpr;
    }

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

    //! return true if the expression can be evaluated (TODO : iterate over symbols expression)
    bool isEvaluable() const
    {
        return M_isNumericExpression || ( M_indexSymbolXYZ.empty() && M_indexSymbolN.empty() && (M_syms.size() == M_symbolNameToValue.size()) );
    }
    bool isConstant() const { return this->isEvaluable(); }

    void setNumericValue( double val ) { M_isNumericExpression = true; M_numericValue = evaluate_type::Constant( val ); }


    const std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type> > >  indices() const
    {
        std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type> > > indices_vec;
        hana::for_each( M_expr.tupleExpr, [this,&indices_vec]( auto const& evec )
                        {
                            for ( auto const& e : evec )
                            {
                                std::vector<std::tuple<uint16_type,uint16_type,uint16_type> > tmp;
                                if ( e.componentSuffix().empty() ) // no suffix comp
                                {
                                    uint16_type idx = this->index( e.symbol() );
                                    if ( idx != invalid_v<uint16_type> )
                                        tmp.push_back( std::make_tuple( idx, 0, 0 ) );
                                }
                                else
                                {
                                    for ( auto const& [_suffix,compArray] : e.componentSuffix() )
                                    {
                                        uint16_type idx = this->index( e.symbol() + _suffix );
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

    Eigen::MatrixXd
    evaluate( std::map<std::string,value_type> const& mp  )
    {
        this->setParameterValues( mp );
        return this->evaluateImpl( true, Environment::worldCommPtr() );
    }

    Eigen::MatrixXd
    evaluate( bool parallel = true, worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr() ) const
    {
        return this->evaluateImpl( parallel, worldcomm );
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

        //typedef typename expression_type::value_type value_type;
        typedef double value_type;

        using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        using shape = ShapeGeneric<gmc_t<Geo_t>::nDim,M,N>;
        typedef std::vector<evaluate_type> loc_type;
        //typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_is_constant( expr.isConstant() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), evaluate_type::Zero() ),
            M_x( expr.parameterValue() ),
            M_yConstant( (M_is_constant)? expr.evaluate() : evaluate_type::Zero() ),
            M_t_expr_index( expr.indices() )
            {
                this->initSubTensor( geom, fev, feu );
            }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_is_constant( expr.isConstant() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), evaluate_type::Zero() ),
            M_x( expr.parameterValue() ),
            M_yConstant( (M_is_constant)? expr.evaluate() : evaluate_type::Zero() ),
            M_t_expr_index( expr.indices() )
            {
                this->initSubTensor( geom, fev );
            }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_fun( expr.fun() ),
            M_is_constant( expr.isConstant() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_nsyms( expr.syms().size() ),
            M_y( M_gmc->nPoints(), evaluate_type::Zero() ),
            M_x( expr.parameterValue() ),
            M_yConstant( (M_is_constant)? expr.evaluate() : evaluate_type::Zero() ),
            M_t_expr_index( expr.indices() )
            {
                this->initSubTensor( geom );
            }

        template<typename IM>
        void init( IM const& im ) {}

        FEELPP_DONT_INLINE void updateFun( Geo_t const& geom )
            {
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
                    hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSymbolsExpr> ), [this,&geom,&k,&q]( auto seId ) {
                            auto const& evec = hana::at( M_expr.symbolsExpression().tupleExpr, hana::int_c<seId> );
                            int nSubExpr = evec.size();
                            auto & evecTensorExpr = hana::at(M_t_expr, hana::int_c<seId> );
                            for (int l=0;l<nSubExpr;++l,++k)
                            {
                                if ( M_t_expr_index[k].empty() )
                                    continue;
                                auto const& e = evec[l];
                                auto const& theexprBase = e.expr();
                                using subtensor_type = typename TransformExprToTensor<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr.symbolsExpression() ))>>::type;
                                auto & subTensor = std::any_cast<subtensor_type&>(evecTensorExpr[l]);
                                for ( auto const& [idx,c1,c2] : M_t_expr_index[k] )
                                    M_x[idx] = subTensor.evalq( c1, c2, q );
                            }
                        });

                    M_fun(&ni,M_x.data(),&no,M_y[q].data());
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
#if !defined( NDEBUG )
                CHECK( c1 < M ) << "invalid c1 " << c1 << " sould be less than " << M
                                << "shape : " << shape::M << " " << shape::N << " is_scalar=" << shape::is_scalar << " is_vectorial=" << shape::is_vectorial << " is_tensor2="<<shape::is_tensor2;
                CHECK( c2 < N ) << "invalid c2 " << c2 << " sould be less than " << N;
#endif
                return ( M_is_constant )? M_yConstant(c1,c2) : M_y[q](c1,c2);
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
                            auto const& theexprBase = e.expr();
                            auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr.symbolsExpression() ))> const&>(evecExpand[l]);
                            evecTensorExpr[l] = typename TransformExprToTensor<std::decay_t<decltype(theexpr)>>::type( theexpr,theInitArgs... );
                        }
                    });
            }
        template<typename... TheArgsType>
        void updateImpl( Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                //if ( M_is_zero ) return;
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
                            auto const& theexprBase = e.expr();
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
        const bool M_is_constant;
        gmc_ptrtype M_gmc;
        int M_nsyms;
        loc_type M_y;
        vec_type M_x;
        const evaluate_type M_yConstant;
        tuple_tensor_expr_type M_t_expr;
        const std::vector<std::vector<std::tuple<uint16_type,uint16_type,uint16_type>>> M_t_expr_index; // (id,c1,c2)
    };

private :

    void updateNumericExpression()
    {
        auto resToNum = toNumericValues( M_fun );
        CHECK( resToNum.size() == M*N )  << "invalid size " << resToNum.size() << " vs " << M*N;
        M_isNumericExpression = true;
        for ( auto const& bd : resToNum )
        {
            if ( !bd.first )
            {
                M_isNumericExpression = false;
                break;
            }
        }
        if ( M_isNumericExpression )
        {
            for( int i = 0; i < M; ++i )
                for( int j = 0; j < N; ++j )
                    M_numericValue(i,j) = resToNum[i*N+j].second;
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
                                if ( e.componentSuffix().empty() )
                                {
                                    uint16_type idx = this->index( e.symbol() );
                                    if ( idx == invalid_uint16_type_value )
                                        continue;

                                    // call an optional update function
                                    if ( e.updateFunction() )
                                        e.updateFunction()();

                                    auto const& theexprBase = e.expr();
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
                                    if ( e.updateFunction() )
                                        e.updateFunction()();

                                    auto const& theexprBase = e.expr();
                                    auto theexpr = theexprBase.applySymbolsExpr( M_expr );
                                    evecExpand[k] = theexpr;

                                    M_context = M_context | Feel::vf::dynamicContext( theexpr );

                                    for ( auto const& [_suffix,compArray] : e.componentSuffix() )
                                    {
                                        uint16_type idx = this->index( e.symbol() + _suffix );
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

        // int no = M_fun.nops();
        // CHECK( no == M*N )  << "invalid size " << no << " vs " << M*N << "\n" ;
        static constexpr int no = M*N;
        if ( symbExprArePolynomials )
        {
            GiNaC::lst symlIsPoly;
            for ( auto const& thesymb : symbTotalDegree )
                symlIsPoly.append( thesymb.first );

            std::vector<bool> compExprIsPoly( no );
            if constexpr ( no == 1 )
                 compExprIsPoly[0] = M_fun.is_polynomial( symlIsPoly );
            else
            {
                for( int i = 0; i < no; ++i )
                    compExprIsPoly[i] = M_fun.op(i).is_polynomial( symlIsPoly );
            }
            M_isPolynomial = std::accumulate( compExprIsPoly.begin(), compExprIsPoly.end(), true, std::logical_and<bool>() );
        }
        else
            M_isPolynomial = false;

        if ( M_isPolynomial )
        {
            std::vector<uint16_type>  compExprPolyOrder( no );
            if constexpr ( no == 1 )
                compExprPolyOrder[0] = totalDegree( M_fun, symbTotalDegree );
            else
            {
                for( int i = 0; i < no; ++i )
                    compExprPolyOrder[i] = totalDegree( M_fun.op(i), symbTotalDegree );
            }
            M_polynomialOrder = *std::max_element( compExprPolyOrder.begin(), compExprPolyOrder.end() );
        }
        else
            M_polynomialOrder = Order;
    }

    evaluate_type
    evaluateImpl( bool parallel, worldcomm_ptr_t const& worldcomm ) const
    {
        if ( M_isNumericExpression )
            return M_numericValue;
        int no = M*N;
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

                                if ( e.componentSuffix().empty() )
                                {
                                    uint16_type idx = this->index( e.symbol() );
                                    if ( idx == invalid_v<uint16_type> )
                                        continue;

                                    auto const& theexprBase = e.expr();
                                    auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr ))> const&>(evecExpand[k]);

                                    x[idx] = theexpr.evaluate( parallel, worldcomm )(0,0);
                                }
                                else
                                {
                                    if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                        continue;

                                    auto const& theexprBase = e.expr();
                                    auto const& theexpr = std::any_cast<std::decay_t<decltype(theexprBase.applySymbolsExpr( M_expr ))> const&>(evecExpand[k]);

                                    for ( auto const& [_suffix,compArray] : e.componentSuffix() )
                                    {
                                        uint16_type idx = this->index( e.symbol() + _suffix );
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

        evaluate_type res = evaluate_type::Zero();
        (*M_cfun)(&ni,x.data(),&no,res.data());
        return res;
    }

    template <typename AnElementOfSymbolExprType>
    bool hasAtLeastOneSymbolDependency( AnElementOfSymbolExprType const& e ) const
    {
        bool res = false;
        if ( e.componentSuffix().empty() )
        {
            res = this->index( e.symbol() ) != invalid_v<uint16_type>;
        }
        else
        {
            for ( auto const& [_suffix,compArray] : e.componentSuffix() )
                if ( this->index( e.symbol() + _suffix ) != invalid_v<uint16_type> )
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
                                std::string const& currentSymbName = e.symbol();
                                // if ( !super::hasSymbol( currentSymbName ) )
                                //     continue;
                                if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                    continue;

                                auto const& currentExpr = e.expr();
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
                                std::string const& currentSymbName = e.symbol();
                                if ( res.find( currentSymbName ) != res.end() )
                                    continue;

                                // if ( !super::hasSymbol( currentSymbName ) )
                                //     continue;
                                if ( !this->hasAtLeastOneSymbolDependency( e ) )
                                    continue;

                                auto const& currentExpr = e.expr();
                                if ( !currentExpr.hasSymbolDependency( symb, se ) )
                                    continue;

                                if ( e.componentSuffix().empty() )
                                    res[currentSymbName].insert( "" );
                                else
                                {
                                    for ( auto const& [_suffix,compArray] : e.componentSuffix() )
                                        res[currentSymbName].insert( _suffix );
                                }

                                currentExpr.dependentSymbols( symb, res, se );
                            }
                        });
    }


    template <int diffOrder,typename TheSymbolExprType>
        auto diffSymbolsExprImpl( std::map< std::string,std::map<std::string, std::set<std::string>>> const& mapSymbolToApplyDiff,
                                  WorldComm const& world, std::string const& dirLibExpr,
                                  TheSymbolExprType const& newse ) const
    {
        auto tupleDiffSymbolsExpr = hana::transform( this->symbolsExpression().tuple(), [this,&mapSymbolToApplyDiff,&newse,&world,&dirLibExpr](auto const& evec)
                                                     {
                                                         int nSubExpr = evec.size();
                                                         using _expr_type = std::decay_t<decltype( evec.front().expr() )>;
                                                         using _expr_diff_type =  std::decay_t<decltype( evec.front().expr().template diff<diffOrder>( "",world,dirLibExpr,newse ) )>;
                                                         symbol_expression_t<_expr_diff_type> seDiffExpr;
                                                         for (int l=0;l<nSubExpr;++l)
                                                         {
                                                             auto const& e = evec[l];
                                                             std::string const& currentSymbName = e.symbol();
                                                             auto const& theexpr = e.expr();

                                                             for ( auto const& [diffVariable,symbolToApplyDiff] : mapSymbolToApplyDiff )
                                                             {
                                                                 auto itFindSymb = symbolToApplyDiff.find( currentSymbName );
                                                                 if ( itFindSymb == symbolToApplyDiff.end() ) // not depend on diffVariable
                                                                     continue;

                                                                 auto currentDiffExpr = theexpr.template diff<diffOrder>( diffVariable,world,dirLibExpr,newse );
                                                                 std::string currentDiffSymbName = (boost::format( "diff_%1%_%2%_%3%" )%currentSymbName %diffVariable %diffOrder ).str();
                                                                 seDiffExpr.add( currentDiffSymbName, currentDiffExpr, e.componentSuffix() );
                                                             }
                                                         }
                                                         return seDiffExpr;
                                                     });
        return SymbolsExpr( std::move( tupleDiffSymbolsExpr ) );
    }

private:
    mutable expression_type  M_fun;
    std::shared_ptr<GiNaC::FUNCP_CUBA> M_cfun;
    std::string M_filename;
    std::string M_exprDesc;
    symbols_expression_type M_expr;
    tuple_expand_symbols_expr_type M_expandSymbolsExpr;
    bool M_isPolynomial;
    uint16_type M_polynomialOrder;
    evaluate_type M_numericValue;

}; // GinacMatrix

template<int M,int N,int Order,typename SymbolsExprType>
FEELPP_EXPORT std::ostream&
operator<<( std::ostream& os, GinacMatrix<M,N,Order,SymbolsExprType> const& e )
{
    os << e.expression();
    return os;
}

template<int M, int N, int Order,typename SymbolsExprType>
FEELPP_EXPORT std::string
str( GinacMatrix<M,N,Order,SymbolsExprType> && e )
{
    return str( std::forward<GinacMatrix<M,N,Order,SymbolsExprType>>(e).expression() );
}

template<int M, int N, int Order,typename SymbolsExprType>
FEELPP_EXPORT std::string
str( GinacMatrix<M,N,Order,SymbolsExprType> const& e )
{
    return str( e.expression() );
}

template<int Order=2, typename SymbolsExprType = symbols_expression_empty_t >
using GinacExVF =GinacMatrix<1,1,Order,SymbolsExprType>;

template<int Order = 2>
using GinacEx = GinacExVF<Order>;

/// \endcond
/// \endcond
}} // Feel



#endif /* FEELPP_DETAIL_GINACMATRIX_HPP */
