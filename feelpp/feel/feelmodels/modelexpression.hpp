/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 22 Feb 2018

 Copyright (C) 2018 Feel++ Consortium

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
#ifndef FEELPP_MODELEXPRESSION_HPP
#define FEELPP_MODELEXPRESSION_HPP 1

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelmodels/modelindexes.hpp>

namespace Feel {

class FEELPP_EXPORT ModelExpression
{
public :

    static const uint16_type expr_order = 2;
    typedef scalar_field_expression<expr_order> expr_scalar_type;
    typedef vector_field_expression<2,1,expr_order> expr_vectorial2_type;
    typedef vector_field_expression<3,1,expr_order> expr_vectorial3_type;
    typedef matrix_field_expression<2,2,expr_order> expr_matrix22_type;
    typedef matrix_field_expression<3,3,expr_order> expr_matrix33_type;

    static constexpr auto expr_shapes = hana::make_tuple( hana::make_tuple(hana::int_c<1>,hana::int_c<1>),
                                                          hana::make_tuple(hana::int_c<2>,hana::int_c<1>),
                                                          hana::make_tuple(hana::int_c<3>,hana::int_c<1>),
                                                          hana::make_tuple(hana::int_c<2>,hana::int_c<2>),
                                                          hana::make_tuple(hana::int_c<3>,hana::int_c<3>) );

    using value_type = typename expr_scalar_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic >;

    ModelExpression() = default;
    explicit ModelExpression( value_type val ) { this->setExprScalar( Feel::vf::expr( val ) ); }
    ModelExpression( ModelExpression const& ) = default;
    ModelExpression( ModelExpression && ) = default;
    ModelExpression& operator=( ModelExpression const& ) = default;
    ModelExpression& operator=( ModelExpression && ) = default;

    bool hasExprScalar() const { return (M_exprScalar)? true : false; }
    bool hasExprVectorial2() const { return (M_exprVectorial2)? true : false; }
    bool hasExprVectorial3() const { return (M_exprVectorial3)? true : false; }
    bool hasExprMatrix22() const { return (M_exprMatrix22)? true : false; }
    bool hasExprMatrix33() const { return (M_exprMatrix33)? true : false; }
    template <int M,int N>
    bool hasExprMatrix() const
    {
        if ( M==2 && N==2 ) return this->hasExprMatrix22();
        else if ( M==3 && N==3 ) return this->hasExprMatrix33();
        else return false;
    }
    template <int M,int N>
    bool hasExpr() const
    {
        if ( M==1 && N == 1 ) return this->hasExprScalar();
        else if ( M == 2 && N == 1 ) return this->hasExprVectorial2();
        else if ( M == 3 && N == 1 ) return this->hasExprVectorial3();
        else if ( M==2 && N==2 ) return this->hasExprMatrix22();
        else if ( M==3 && N==3 ) return this->hasExprMatrix33();
        else return false;
    }

    bool hasExpr(int M,int N) const
    {
        if ( M==1 && N == 1 ) return this->hasExprScalar();
        else if ( M == 2 && N == 1 ) return this->hasExprVectorial2();
        else if ( M == 3 && N == 1 ) return this->hasExprVectorial3();
        else if ( M==2 && N==2 ) return this->hasExprMatrix22();
        else if ( M==3 && N==3 ) return this->hasExprMatrix33();
        else return false;
    }

    bool hasAtLeastOneExpr() const { return  this->hasExprScalar() || this->hasExprVectorial2() || this->hasExprVectorial3() || this->hasExprMatrix22() || this->hasExprMatrix33(); }

    bool isScalar() const { return hasExprScalar(); }
    bool isVector() const { return hasExprVectorial2() || hasExprVectorial3(); }
    bool isMatrix() const { return hasExprMatrix22() || hasExprMatrix33(); }

    bool isEvaluable() const
    {
        if ( !this->hasAtLeastOneExpr() )
            return false;
        bool res = true;
        hana::for_each( expr_shapes, [this,&res]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                                res = res && this->expr<ni,nj>().expression().isEvaluable();
                        });
        return res;
    }
    bool isConstant() const { return this->isEvaluable(); }

    evaluate_type evaluate() const
    {
        evaluate_type res;
        hana::for_each( expr_shapes, [this,&res]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                                res = this->expr<ni,nj>().evaluate();
                        });
        return res;
    }


    value_type value() const { CHECK( this->hasExprScalar() && this->exprScalar().expression().isEvaluable() ) << "expression is not scalar or evaluable";return this->exprScalar().evaluate()(0,0); }
    void setValue( value_type val ) { CHECK( this->hasExprScalar() ) << "expr should be scalar"; this->exprScalar().expression().setNumericValue( val ); }

    expr_scalar_type const& exprScalar() const { CHECK( this->hasExprScalar() ) << "no Scalar expression"; return *M_exprScalar; }
    expr_vectorial2_type const& exprVectorial2() const { CHECK( this->hasExprVectorial2() ) << "no Vectorial2 expression"; return *M_exprVectorial2; }
    expr_vectorial3_type const& exprVectorial3() const { CHECK( this->hasExprVectorial3() ) << "no Vectorial3 expression"; return *M_exprVectorial3; }
    expr_matrix22_type const& exprMatrix22() const { CHECK( this->hasExprMatrix22() ) << "no Matrix22 expression"; return *M_exprMatrix22; }
    expr_matrix33_type const& exprMatrix33() const { CHECK( this->hasExprMatrix33() ) << "no Matrix33 expression"; return *M_exprMatrix33; }

    expr_scalar_type & exprScalar() { CHECK( this->hasExprScalar() ) << "no Scalar expression"; return *M_exprScalar; }
    expr_vectorial2_type & exprVectorial2() { CHECK( this->hasExprVectorial2() ) << "no Vectorial2 expression"; return *M_exprVectorial2; }
    expr_vectorial3_type & exprVectorial3() { CHECK( this->hasExprVectorial3() ) << "no Vectorial3 expression"; return *M_exprVectorial3; }
    expr_matrix22_type & exprMatrix22() { CHECK( this->hasExprMatrix22() ) << "no Matrix22 expression"; return *M_exprMatrix22; }
    expr_matrix33_type & exprMatrix33() { CHECK( this->hasExprMatrix33() ) << "no Matrix33 expression"; return *M_exprMatrix33; }

    template <int M=1,int N=1>
        expr_scalar_type const& expr( typename std::enable_if< M==1 && N == 1>::type* = nullptr ) const { return this->exprScalar(); }
    template <int M=1,int N=1>
        expr_vectorial2_type const& expr( typename std::enable_if< M==2 && N == 1>::type* = nullptr ) const { return this->exprVectorial2(); }
    template <int M=1,int N=1>
        expr_vectorial3_type const& expr( typename std::enable_if< M==3 && N == 1>::type* = nullptr ) const { return this->exprVectorial3(); }
    template <int M=1,int N=1>
        expr_matrix22_type const& expr( typename std::enable_if< M==2 && N == 2>::type* = nullptr ) const { return this->exprMatrix22(); }
    template <int M=1,int N=1>
        expr_matrix33_type const& expr( typename std::enable_if< M==3 && N == 3>::type* = nullptr ) const { return this->exprMatrix33(); }

    template <int M=1,int N=1>
        expr_scalar_type & expr( typename std::enable_if< M==1 && N == 1>::type* = nullptr ) { return this->exprScalar(); }
    template <int M=1,int N=1>
        expr_vectorial2_type & expr( typename std::enable_if< M==2 && N == 1>::type* = nullptr ) { return this->exprVectorial2(); }
    template <int M=1,int N=1>
        expr_vectorial3_type & expr( typename std::enable_if< M==3 && N == 1>::type* = nullptr ) { return this->exprVectorial3(); }
    template <int M=1,int N=1>
        expr_matrix22_type & expr( typename std::enable_if< M==2 && N == 2>::type* = nullptr ) { return this->exprMatrix22(); }
    template <int M=1,int N=1>
        expr_matrix33_type & expr( typename std::enable_if< M==3 && N == 3>::type* = nullptr ) { return this->exprMatrix33(); }

    template <int M,int N>
        expr_matrix22_type const& exprMatrix( typename std::enable_if< M==2 && N == 2>::type* = nullptr ) const { return this->exprMatrix22(); }
    template <int M,int N>
        expr_matrix33_type const& exprMatrix( typename std::enable_if< M==3 && N == 3>::type* = nullptr ) const { return this->exprMatrix33(); }

    void setExprScalar( expr_scalar_type const& expr ) { M_exprScalar = boost::optional<expr_scalar_type>( expr ); }
    void setExprVectorial2( expr_vectorial2_type const& expr ) { M_exprVectorial2 = boost::optional<expr_vectorial2_type>( expr ); }
    void setExprVectorial3( expr_vectorial3_type const& expr ) { M_exprVectorial3 = boost::optional<expr_vectorial3_type>( expr ); }
    void setExprMatrix22( expr_matrix22_type const& expr ) { M_exprMatrix22 = boost::optional<expr_matrix22_type>( expr ); }
    void setExprMatrix33( expr_matrix33_type const& expr ) { M_exprMatrix33 = boost::optional<expr_matrix33_type>( expr ); }

    //! set an expression from a key of a ptree p
    void setExpr( std::string const& key, pt::ptree const& p, WorldComm const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes = ModelIndexes() );
    void setExpr( std::string const& expr, WorldComm const& worldComm = Environment::worldComm(), std::string const& directoryLibExpr = "" );

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        hana::for_each( expr_shapes, [this,&mp]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                                this->expr<ni,nj>().setParameterValues( mp );
                        });
    }

    void eraseParameterValues( std::set<std::string> const& symbNames )
    {
        hana::for_each( expr_shapes, [this,&symbNames]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                                this->expr<ni,nj>().expression().eraseParameterValues( symbNames );
                        });
    }

    void updateParameterValues( std::set<std::string> const& symbNames, std::map<std::string,double> & pv ) const
    {
        hana::for_each( expr_shapes, [this,&symbNames,&pv]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                            {
                                auto const& theexpr = this->expr<ni,nj>();
                                if ( theexpr.expression().isEvaluable() )
                                {
                                    auto theEvalExpr = theexpr.evaluate();
                                    if constexpr ( ni == 1 && nj == 1 )
                                                 {
                                                     for ( std::string const& s : symbNames )
                                                         pv[s] = theEvalExpr(0,0);
                                                 }
                                    else
                                    {
                                        for ( auto const& [_suffix,compArray] : SymbolExprComponentSuffix(ni,nj ) )
                                        {
                                            uint16_type c1 = compArray[0];
                                            uint16_type c2 = compArray[1];
                                            for ( std::string const& s : symbNames )
                                                pv[ s + _suffix ] = theEvalExpr(c1,c2);
                                        }
                                    }
                                }
                            }
                        });
    }

    void updateParameterValues( std::string const& symbName, std::map<std::string,double> & pv ) const { this->updateParameterValues( std::set<std::string>({ symbName }), pv );}

    void updateSymbolNames( std::set<std::string> const& symbNames, std::set<std::string> & outputSymbNames ) const
    {
        hana::for_each( expr_shapes, [this,&symbNames,&outputSymbNames]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                            {
                                if constexpr ( ni == 1 && nj == 1 )
                                {
                                    outputSymbNames.insert( symbNames.begin(), symbNames.end() );
                                }
                                else
                                {
                                    for ( auto const& [_suffix,compArray] : SymbolExprComponentSuffix(ni,nj ) )
                                    {
                                        for ( std::string const& s : symbNames )
                                            outputSymbNames.insert( s + _suffix );
                                    }
                                }
                            }
                        });
    }
    void updateSymbolNames( std::string const& symbName, std::set<std::string> & outputSymbNames ) const { this->updateSymbolNames( std::set<std::string>({ symbName }), outputSymbNames ); }

    //! rename symbols in symbolics expr with mapping \old2new
    void renameSymbols( std::map<std::string,std::string> const& old2new )
    {
        hana::for_each( expr_shapes, [this,&old2new]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                                this->expr<ni,nj>().expression().renameSymbols( old2new );
                        });
    }


    template<typename ExprT>
    constexpr auto
    exprScalar( std::string const& symb, ExprT const& e ) const
    {
        return Feel::vf::expr( this->exprScalar(), symbolExpr(symb, e) );
    }

    std::string exprToString() const
    {
        std::string res;
        hana::for_each( expr_shapes, [this,&res]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                                res = str( this->expr<ni,nj>().expression() ); // we guess that we have only one expression
                        });
        return res;
    }

    std::tuple<std::string,SymbolExprComponentSuffix> exprInformations() const
    {
        std::string exprStr;
        SymbolExprComponentSuffix comps;
        hana::for_each( expr_shapes, [this,&exprStr,&comps]( auto const& e_ij )
                        {
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )  // we guess that we have only one expression
                            {
                                exprStr = str( this->expr<ni,nj>().expression() );
                                comps = SymbolExprComponentSuffix( ni,nj );
                            }
                        });
        return std::make_tuple( std::move(exprStr), std::move(comps) );
    }

    template <typename TheSymbolExprType = symbols_expression_empty_t>
    bool hasSymbolDependency( std::set<std::string> const& symbolsStr, TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
    {
        if ( this->isConstant() )
            return false;

        if ( symbolsStr.empty() )
            return false;

        bool res = false;
        hana::for_each( expr_shapes, [this,&symbolsStr,&se,&res]( auto const& e_ij )
                        {
                            if ( res )
                                return;
                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
                            if ( this->hasExpr<ni,nj>() )
                            {
                                for ( std::string const& symbolStr : symbolsStr )
                                {
                                    if ( this->expr<ni,nj>().hasSymbolDependency( symbolStr, se ) )
                                    {
                                        res = true;
                                        break;
                                    }
                                }
                            }
                        });
        return res;
    }

    template <typename TheSymbolExprType = symbols_expression_empty_t>
    bool hasSymbolDependency( std::string const& symbolStr, TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
    {
        return this->hasSymbolDependency( std::set<std::string>( { symbolStr } ), se );
    }

    template <int Dim, typename TheSymbolExprType = symbols_expression_empty_t>
    bool hasSymbolDependencyOnCoordinatesInSpace( TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
    {
        std::set<std::string> coords( { "x","y","z" } );
        if ( Dim >=1 )
            coords.insert( "x" );
        if ( Dim >=2 )
            coords.insert( "y" );
        if ( Dim ==3 )
            coords.insert( "z" );
        return this->hasSymbolDependency( coords, se );
    }

private :
    boost::optional<expr_scalar_type> M_exprScalar;
    boost::optional<expr_vectorial2_type> M_exprVectorial2;
    boost::optional<expr_vectorial3_type> M_exprVectorial3;
    boost::optional<expr_matrix22_type> M_exprMatrix22;
    boost::optional<expr_matrix33_type> M_exprMatrix33;
};

class FEELPP_EXPORT ModelExpressionScalar : private ModelExpression
{
    typedef ModelExpression super_type;
public :
    static const uint16_type expr_order = super_type::expr_order;
    typedef typename super_type::expr_scalar_type expr_scalar_type;

    ModelExpressionScalar() = default;
    ModelExpressionScalar( ModelExpressionScalar const& ) = default;
    ModelExpressionScalar( ModelExpressionScalar && ) = default;
    ModelExpressionScalar& operator=( ModelExpressionScalar const& ) = default;
    ModelExpressionScalar& operator=( ModelExpressionScalar && ) = default;

    bool isConstant() const { return super_type::isConstant(); }
    double value() const { return super_type::value(); }
    expr_scalar_type const& expr() const { return this->exprScalar(); }
    bool hasExpr() const { return this->hasExprScalar(); }

    void setExpr( expr_scalar_type const& expr ) { this->setExprScalar( expr ); }

    void setParameterValues( std::map<std::string,double> const& mp ) { super_type::setParameterValues( mp ); }

    template<typename ExprT>
    constexpr auto //Expr< GinacExVF<ExprT,expr_order> >
    expr( std::string const& symb, ExprT const& e ) const
    {
        return this->exprScalar( symb,e );
    }

};

#if 0
template <int M,int N>
class ModelExpressionMatrix : private ModelExpression
{
    typedef ModelExpression super_type;
public :
    static const uint16_type expr_order = super_type::expr_order;
    typedef typename mpl::if_<mpl::bool_< M==2 && N==2 >,
                              typename super_type::expr_matrix22_type,
                              typename super_type::expr_matrix33_type>::type expr_matrix_type;

    ModelExpressionMatrix() = default;
    ModelExpressionMatrix( ModelExpressionMatrix const& ) = default;
    ModelExpressionMatrix( ModelExpressionMatrix && ) = default;
    ModelExpressionMatrix& operator=( ModelExpressionMatrix const& ) = default;
    ModelExpressionMatrix& operator=( ModelExpressionMatrix && ) = default;

    //bool isConstant() const { return !this->hasExprScalar(); }
    //double value() const { return super_type::value(); }
    expr_matrix_type const& expr() const { return this->exprMatrix<M,N>(); }

    void setValue( double d ) { super_type::setValue( d ); }
    void setExpr( expr_matrix_type const& expr ) { this->setExprScalar( expr ); }

    void setParameterValues( std::map<std::string,double> const& mp ) { super_type::setParameterValues( mp ); }
#if 0
    template<typename ExprT>
    Expr< GinacExVF<ExprT,expr_order> >
    expr( std::string const& symb, ExprT const& e ) const
    {
        return this->exprScalar( symb,e );
    }
#endif
};
#endif
} // namespace Feel

#endif
