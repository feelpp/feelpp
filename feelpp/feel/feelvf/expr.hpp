/* -*- Mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file expr.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#ifndef FEELPP_EXPR_HPP
#define FEELPP_EXPR_HPP 1

#undef max
#include <algorithm>

//#include <boost/version.hpp>
#include <boost/none.hpp>
#include <boost/static_assert.hpp>
#include <boost/foreach.hpp>

#include <boost/multi_array.hpp>

#include <Eigen/Core>

#include <feel/feelcore/environment.hpp>
#include <feel/feelpoly/policy.hpp>
//#include <feel/feelpoly/context.hpp>

#include <feel/feelvf/exprbase.hpp>
#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/shape.hpp>
#include <feel/feelvf/lambda.hpp>
#include <feel/feelvf/symbolsexpr.hpp>

namespace Feel
{
namespace vf
{
class GiNaCBase;

/// \cond detail
typedef node<double>::type node_type;

enum
{
    CONTEXT_1 = ( 1<<0 ), /**< identifier 1 for the context */
    CONTEXT_2 = ( 1<<1 )  /**< identifier 2 for the context */
};

template<typename ExprT>
class ComponentsExpr
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    template<typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;



    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;
    typedef ComponentsExpr<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    ComponentsExpr()
        :
        M_expr()
    {}

    explicit ComponentsExpr( expression_type const & __expr, int c1, int c2 )
        :
        M_expr( __expr ),
        M_c1( c1 ),
        M_c2( c2 )
    {}
    ~ComponentsExpr()
    {}

    //@}

    //! dynamic context
    size_type dynamicContext() const { return Feel::vf::dynamicContext( M_expr ); }

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    evaluate_type
    evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            return evaluate_type::Constant( M_expr.evaluate( p, worldcomm )(M_c1,M_c2) );
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            M_expr.setParameterValues( mp );
        }
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            M_expr.updateParameterValues( pv );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newExpr =  M_expr.applySymbolsExpr( se );
            using new_expr_type = std::decay_t<decltype(newExpr)>;
            return ComponentsExpr<new_expr_type>( newExpr, M_c1, M_c2 );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr.hasSymbolDependency( symb, se );
        }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_expr.dependentSymbols( symb, res, se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            auto theDiffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
            using new_expr_type = std::decay_t<decltype(theDiffExpr)>;
            return ComponentsExpr<new_expr_type>( theDiffExpr,M_c1,M_c2 );
        }

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {

        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::NDim, Scalar, false> shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu ),
            M_c1( expr.M_c1 ),
            M_c2( expr.M_c2 )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev ),
            M_c1( expr.M_c1 ),
            M_c2( expr.M_c2 )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom ),
            M_c1( expr.M_c1 ),
            M_c2( expr.M_c2 )
        {
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... ),
            M_c1( exprExpanded.M_c1 ),
            M_c2( exprExpanded.M_c2 )
            {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
            }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evalijq( i, j, M_c1, M_c2, q );
        }

        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr.evalijq( i, j, M_c1, M_c2, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evaliq( i, M_c1, M_c2, q );
        }

        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_tensor_expr.evalq( M_c1, M_c2, q );
        }

        tensor_expr_type M_tensor_expr;
        const int M_c1, M_c2;
    };
    expression_type M_expr;
    int M_c1, M_c2;
};

class IntegratorBase {};


// type T can be used with vf expr
template <typename T, typename = void>
struct is_vf_expr : std::false_type {};
template <typename T>
struct is_vf_expr<T, std::void_t<decltype(std::declval<T>().context,
                                          std::declval<T>().is_terminal) >> : std::true_type {};
template <typename T>
constexpr bool is_vf_expr_v = is_vf_expr<T>::value;


template <typename T, typename = void>
struct has_evaluate_without_context : std::false_type {};
template <typename T>
struct has_evaluate_without_context<T, std::void_t<decltype(std::declval<T>().evaluate( true )) >>
    : std::true_type {};
template <typename T>
constexpr bool has_evaluate_without_context_v = has_evaluate_without_context<T>::value;

template <typename T, typename = void>
struct evaluate_expression_type
{
    using type = Eigen::Matrix<typename T::value_type,Eigen::Dynamic,Eigen::Dynamic>;
};
template <typename T>
struct evaluate_expression_type<T, std::void_t<typename T::evaluate_type>>
{
    using type = typename T::evaluate_type;
};
template <typename T>
using evaluate_expression_t = typename evaluate_expression_type<T>::type;

template <typename T, typename = void>
struct has_symbolic_parameter_values_type : std::false_type {};
template <typename T>
struct has_symbolic_parameter_values_type <T, std::void_t<decltype(std::declval<T>().setParameterValues( std::map<std::string,double/*typename T::value_type*/>{} )) >>
    : std::true_type {};
template <typename T>
constexpr bool has_symbolic_parameter_values_v = has_symbolic_parameter_values_type<T>::value;

template <typename T, typename TheSymbolExprType, typename = void>
struct has_symbol_dependency_type : std::false_type {};
template <typename T,typename TheSymbolExprType>
struct has_symbol_dependency_type <T,TheSymbolExprType,std::void_t<decltype(std::declval<T>().hasSymbolDependency( "", TheSymbolExprType{} ) ) >>
    : std::true_type {};
template <typename T,typename TheSymbolExprType>
constexpr bool has_symbol_dependency_v = has_symbol_dependency_type<T,TheSymbolExprType>::value;

template <typename T, typename TheSymbolExprType, typename = void>
struct has_dependent_symbols_type : std::false_type {};
template <typename T,typename TheSymbolExprType>
struct has_dependent_symbols_type <T,TheSymbolExprType,std::void_t<decltype(std::declval<T>().dependentSymbols( "", std::declval< std::map<std::string,std::set<std::string>> &>(), TheSymbolExprType{} ) ) >>
    : std::true_type {};
template <typename T,typename TheSymbolExprType>
constexpr bool has_dependent_symbols_v = has_dependent_symbols_type<T,TheSymbolExprType>::value;

#if 0
template <typename T, int diffOrder, typename TheSymbolExprType, typename = void>
struct has_symbolic_diff_type : std::false_type {};
template <typename T, int diffOrder, typename TheSymbolExprType>
struct has_symbolic_diff_type <T, diffOrder, TheSymbolExprType, std::void_t<decltype(std::declval<T>().template diff<diffOrder>( "", Feel::worldcomm_t{},"", TheSymbolExprType{} )) >>
    : std::true_type {};
template <typename T, int diffOrder, typename TheSymbolExprType>
constexpr bool has_symbolic_diff_v = has_symbolic_diff_type<T,diffOrder,TheSymbolExprType>::value;
#endif

// forward declarations
template<typename ExprT>
class Expr;

template <typename ExprT>
Expr<ExprT>
expr( ExprT const& exprt, typename std::enable_if_t<is_vf_expr_v<ExprT> >* = nullptr );

template <typename ExprT>
Expr<ExprT>
expr( ExprT && exprt, typename std::enable_if_t<is_vf_expr_v<ExprT> >* = nullptr );

/*!
  \class Expr
  \brief Variational Formulation Expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class Expr : public ExprBase, public ExprDynamicBase //: public std::enable_shared_from_this<Expr<ExprT> >
{
public:

    inline static const size_type context = ExprT::context;
    static const bool is_terminal = ExprT::is_terminal;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    using evaluate_type = evaluate_expression_t<expression_type>;
    typedef Expr<ExprT> this_type;
    typedef std::shared_ptr<this_type> this_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Expr()
        :
        M_expr()
    {}

    explicit Expr( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    explicit Expr( expression_type && __expr )
        :
        M_expr( __expr )
    {}
    ~Expr() override
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename... TheExprs>
    struct Lambda
    {
        typedef typename ExprT::template Lambda<TheExprs...>::type expr_type;
        typedef Expr<expr_type> type;
        //typedef expr_type type;
    };


#if 0
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr...  e  )
        {
            return expr( M_expr( e... ) );
        }
#else
    template<typename TheExpr1>
    typename Lambda<TheExpr1>::type
    operator()( TheExpr1 const&  e   )
        {
            return expr( M_expr( e ) );
        }
    template<typename TheExpr1,typename TheExpr2>
    typename Lambda<Expr<TheExpr1>,Expr<TheExpr2>>::type
    operator()( Expr<TheExpr1> const&  e1, Expr<TheExpr2> const& e2 )
        {
            return expr( M_expr( e1, e2 ) );
        }

    template<typename TheExpr1,typename TheExpr2,typename TheExpr3>
    typename Lambda<TheExpr1,TheExpr2,TheExpr3>::type
    operator()( TheExpr1 const&  e1, TheExpr2 const& e2, TheExpr3 const& e3  )
        {
            return expr( M_expr( e1, e2, e3 ) );
        }
#endif

#if 0
    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e  ) const { return expr(M_expr(e...)); }
#else
    template<typename TheExpr1>
    typename Lambda<TheExpr1>::type
    operator()( TheExpr1 const&  e   ) const
        {
            return expr( M_expr( e ) );
        }
    template<typename TheExpr1,typename TheExpr2>
    typename Lambda<Expr<TheExpr1>,Expr<TheExpr2>>::type
        operator()( Expr<TheExpr1> const&  e1, Expr<TheExpr2> const& e2 ) const
        {
            return expr( M_expr( e1, e2 ) );
        }
    template<typename TheExpr1,typename TheExpr2,typename TheExpr3>
    typename Lambda<TheExpr1,TheExpr2,TheExpr3>::type
    operator()( TheExpr1 const&  e1, TheExpr2 const& e2, TheExpr3 const& e3  ) const
        {
            return expr( M_expr( e1, e2, e3 ) );
        }

#endif

    Expr<ComponentsExpr<Expr<ExprT> > >
    operator()( int c1, int c2 )
    {
        auto ex = ComponentsExpr<Expr<ExprT> >( Expr<ExprT>( M_expr ), c1, c2 );
        return Expr<ComponentsExpr<Expr<ExprT> > >( ex );
    }

    Expr<ComponentsExpr<Expr<ExprT> > >
    operator()( int c1, int c2 ) const
    {
        auto ex = ComponentsExpr<Expr<ExprT> >( Expr<ExprT>( M_expr ), c1, c2 );
        return Expr<ComponentsExpr<Expr<ExprT> > >( ex );
    }

    void setParameterValues( std::pair<std::string,double/*value_type*/> const& mp )
        {
            this->setParameterValues( { { mp.first, mp.second } } );
        }
    void setParameterValues( std::map<std::string,double/*value_type*/> const& mp )
        {
            if constexpr ( has_symbolic_parameter_values_v<expression_type> )
                 M_expr.setParameterValues( mp );
        }
#if 0
    void setParameterValues( std::map<std::string,value_type> const& mp, mpl::bool_<true> )
        {
            M_expr.setParameterValues( mp );
        }
    void setParameterValues( std::map<std::string,value_type> const& mp, mpl::bool_<false> )
        {
        }
#endif
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            if constexpr ( has_symbolic_parameter_values_v<expression_type> )
                  M_expr.updateParameterValues( pv );
        }

#if 0
    template<typename ExprTT>
    explicit Expr( ExprTT const& )
        {
            
        }
    template<typename ExprTT>
    ExprT operator=( ExprTT const& e )
        {
            
        }
#endif
    //! @return the dynamic context of the expression
    size_type dynamicContext() const
        {
            //std::cout << "dynctx:" << Feel::vf::dynamicContext( M_expr ) << " hasp:" << vm::hasPOINT(Feel::vf::dynamicContext( M_expr )) << std::endl;
            return Feel::vf::dynamicContext( M_expr );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto theNewExpr = M_expr.applySymbolsExpr( se );
            if constexpr( std::is_base_of_v<ExprBase, std::decay_t<decltype(theNewExpr)> > )
                return theNewExpr;
            else
                return Feel::vf::expr( std::move( theNewExpr ) );
            //return Feel::vf::expr( M_expr.applySymbolsExpr( se ) );
        }

    //! return true if the symbol \symb is used in the current expression
    //! note: \se is generally not given if the expressionn has already symbols expr, otherwise it allows to link to depencies of expr
    template <typename TheSymbolExprType = symbols_expression_empty_t>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
        {
            if constexpr ( has_symbol_dependency_v<expression_type,TheSymbolExprType> )
                return M_expr.hasSymbolDependency( symb, se );
            else
                return false;
        }

    //! return true if one symbol of the set \symbs is used in the current expression
    //! note: \se is generally not given if the expressionn has already symbols expr, otherwise it allows to link to depencies of expr
    template <typename TheSymbolExprType = symbols_expression_empty_t>
    bool hasSymbolDependency( std::set<std::string> const& symbs, TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
        {
            for ( std::string const& symb : symbs )
            {
                if ( this->hasSymbolDependency( symb, se ) )
                    return true;
            }
            return false;
        }

    //! return true if the expr depends on x,y,z
    //! note: \se is generally not given if the expressionn has already symbols expr, otherwise it allows to link to depencies of expr
    template <int Dim, typename TheSymbolExprType = symbols_expression_empty_t>
    bool hasSymbolDependencyOnCoordinatesInSpace( TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
        {
            std::vector<std::string> coords( { "x","y","z" } );
            coords.resize(Dim);
            for ( std::string const& c : coords )
                if ( this->hasSymbolDependency( c, se ) )
                    return true;
            return false;
        }


    //! update the list of symbol used in the current expression that have a dependency with symbol \symb
    //! note: \se is generally not given if the expressionn has already symbols expr, otherwise it allows to link to depencies of expr
    template <typename TheSymbolExprType = symbols_expression_empty_t>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
        {
            if constexpr ( has_dependent_symbols_v<expression_type,TheSymbolExprType> )
                 M_expr.dependentSymbols( symb, res, se );
        }

    //! symbolic differiantiation
    //! note: \se is generally not given, it's an internal use linked to depencies of expr
    template <int diffOrder,typename TheSymbolExprType = symbols_expression_empty_t>
    auto diff( std::string const& diffSymbol,
               WorldComm const& world = Environment::worldComm(), std::string const& dirLibExpr = "",
               TheSymbolExprType const& se = symbols_expression_empty_t{} ) const
        {
            auto theDiffExpr = M_expr.template diff<diffOrder>( diffSymbol, world, dirLibExpr, se );
            if constexpr( std::is_base_of_v<ExprBase, std::decay_t<decltype(theDiffExpr)> > )
                return theDiffExpr;
            else
                return Feel::vf::expr( std::move( theDiffExpr ) );
        }

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > >,fusion::pair<vf::detail::gmc<1>,std::shared_ptr<vf::detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {

        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        using expr_type = typename this_type::expression_type;
        using key_type = key_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        using gmc_ptrtype = gmc_ptr_t<Geo_t>;
        using shape = typename tensor_expr_type::shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom )
        {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
        {}

        gmc_ptrtype geom() const { return M_geo; }

        int nPoints() const { return M_geo->nPoints(); }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) noexcept
        {
            M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev ) noexcept 
        {
            M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom ) noexcept 
        {
            M_tensor_expr.update( geom );
        }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx ) noexcept 
        {
            M_tensor_expr.updateContext( ctx... );
        }


        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
            }

        value_type
        evalij( uint16_type i, uint16_type j ) const noexcept
        {
            return M_tensor_expr.evalij( i, j );
        }

        Eigen::Map<const Eigen::Matrix<value_type, shape::M,shape::N>>
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const noexcept 
        {
            return M_tensor_expr.evalijq( i, j, q );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const noexcept 
        {
            return M_tensor_expr.evalijq( i, j, c1, c2, q );
        }

        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const noexcept
        {
            return M_tensor_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const noexcept
        {
            return M_tensor_expr.evaliq( i, c1, c2, q );
        }
        Eigen::Map<const Eigen::Matrix<value_type, shape::M,shape::N>>
        evaliq( uint16_type i, uint16_type q ) const noexcept
        {
            return M_tensor_expr.evaliq( i, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const noexcept
        {
            value_type e = M_tensor_expr.evalq( c1, c2, q );
            return e;
        }
        Eigen::Map<const Eigen::Matrix<value_type, shape::M,shape::N>>
        evalq( uint16_type q ) const noexcept
        {
            return M_tensor_expr.evalq( q );
        }

        gmc_ptrtype M_geo;
        //Geo_t const& M_geo;
        tensor_expr_type M_tensor_expr;
    };

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptr<vf::detail::gmc<0> > >,fusion::pair<vf::detail::gmc<1>,std::shared_ptr<vf::detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensorPermutation
    {

        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        using expr_type = typename this_type::expression_type;
        using key_type = key_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        using gmc_ptrtype = gmc_ptr_t<Geo_t>;
        using shape = typename tensor_expr_type::shape;

        tensorPermutation( this_type const& expr, Geo_t const& geom )
            :
            M_geo( fusion::at_key<key_type>( geom ).get() ),
            M_tensor_expr( expr.expression(), geom )
            {
            }

        gmc_ptrtype geom() const { return M_geo; }

        int nPoints() const { return M_geo->nPoints(); }

        void setPermutation( std::vector<uint16_type> const& perm ) { M_permutation = perm; }

        void update( Geo_t const& geom )
            {
                M_tensor_expr.update( geom );
            }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return M_tensor_expr.evalq( c1, c2, M_permutation[q] );
            }
    private :
        gmc_ptrtype M_geo;
        tensor_expr_type M_tensor_expr;
        std::vector<uint16_type> M_permutation;
    };

    template<typename Geo_t,typename ... Basis_t>
    tensor<Geo_t,Basis_t...> evaluator( Geo_t const& geo,const Basis_t&... basisArg  ) const { return tensor<Geo_t,Basis_t...>( *this, geo, basisArg... ); }

    template<typename Geo_t>
    tensorPermutation<Geo_t> evaluatorWithPermutation( Geo_t geo ) const { return tensorPermutation<Geo_t>( *this, geo ); }


private :
    template<typename ElementType>
    struct EvaluatorTraits
    {
    private :
        typedef ElementType element_type;
        typedef typename element_type::gm_type gm_type;
        typedef typename gm_type::template Context<element_type> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename this_type::template tensor<map_gmc_type> eval_expr_type;
        //typedef typename eval_expr_type::shape shape;
    public :
        typedef eval_expr_type type;
    };

public :
    template<typename ElementType>
    using evaluator_t = typename EvaluatorTraits<ElementType>::type;

#if 0
    class Diff
    {
    public:

        value_type value() const
        {
            return __expression.value();
        }

        value_type grad( int __ith ) const
        {
            return __expression.grad( __ith );
        }

        value_type hessian( int __i, int __j ) const
        {
            return __expression.hessian( __i, __j );
        }

    };
#endif /**/
    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    expression_type& expression()
    {
        return M_expr;
    }

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const
    {
        DVLOG(2) << "calling assemble(u,v)\n";
        M_expr.assemble( __u, __v, __f );
        DVLOG(2) << "calling assemble(u,v) done\n";
    }

    template<typename Elem1, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& __v,
                   FormType& __f ) const
    {
        DVLOG(2) << "calling assemble(v)\n";
        M_expr.assemble( __v, __f );
        DVLOG(2) << "calling assemble(v) done\n";
    }

    template<typename P0hType>
    typename P0hType::element_type
    broken( std::shared_ptr<P0hType>& P0h ) const
    {
        return M_expr.broken( P0h );
    }
    //__typeof__( M_expr.evaluate() )
    //ublas::matrix<typename expression_type::value_type>

    evaluate_type
    evaluate( std::pair<std::string,value_type> const& mp  )
    {
        return M_expr.evaluate( { { mp.first, mp.second } } );
    }
    evaluate_type
    evaluate( std::map<std::string,value_type> const& mp  )
    {
        return M_expr.evaluate( mp );
    }
    evaluate_type
    evaluate( bool parallel = true ) const
    {
        if constexpr ( has_evaluate_without_context_v<expression_type> )
        {
            return M_expr.evaluate( parallel );
        }
        else
        {
            CHECK( false ) << "expression can not be evaluated without context";
            return evaluate_type{};
        }
    }
    template<typename T, int M, int N=1>
    decltype(auto)
    evaluate( std::vector<Eigen::Matrix<T,M,N>> const& v, bool parallel = true ) const
    {
        return M_expr.evaluate( v, true );
    }
    typename expression_type::value_type
    evaluateAndSum() const
    {
        return M_expr.evaluateAndSum();
    }
    std::string expressionStr() const
    {
        return std::string();
        //return M_expr.expressionStr();
    }


    //@}

protected:

private:

    mutable expression_type  M_expr;
};


template <typename ExprT>
Expr<ExprT>
expr( ExprT const& exprt, typename std::enable_if_t<is_vf_expr_v<ExprT> >* /*= nullptr*/ )
{
    return Expr<ExprT>( exprt );
}

template <typename ExprT>
Expr<ExprT>
expr( ExprT && exprt, typename std::enable_if_t<is_vf_expr_v<ExprT> >* /*= nullptr*/ )
{
    return Expr<ExprT>( std::forward<ExprT>( exprt ) );
}

template <typename ExprT>
std::shared_ptr<Expr<ExprT> >
exprPtr( ExprT const& exprt )
{
    return std::shared_ptr<Expr<ExprT> >( new Expr<ExprT>( exprt ) );
}

template <typename ExprT>
std::ostream&
operator<<( std::ostream& os, Expr<ExprT> const& exprt )
{
    os << exprt.expression();
    return os;
}

template <typename ExprT>
std::string
str( Expr<ExprT> && exprt )
{
    return str(std::forward<Expr<ExprT>>(exprt).expression());
}


extern Expr<LambdaExpr1> _e1;
extern Expr<LambdaExpr2> _e2;
extern Expr<LambdaExpr3> _e3;
extern Expr<LambdaExpr1V> _e1v;
extern Expr<LambdaExpr2V> _e2v;
extern Expr<LambdaExpr3V> _e3v;

/**
 * \class ExpressionOrder
 *
 * Class that compute the expression polynomial order of \p ExprT
 *
 * \tparam ExprT expression whose approximate order must be computed
 *
 * Note that if the expression is polynomial then the order is exact, however if
 * analytic functions such as exp, cos, sin ... then the order is only an
 * approximation.
 */
template<typename IntElts,typename ExprT>
struct ExpressionOrder
{

    typedef typename boost::tuples::template element<1, IntElts>::type element_iterator_type;
    typedef typename boost::remove_reference<typename element_iterator_type::reference>::type const_t;
    typedef typename boost::unwrap_reference<typename boost::remove_const<const_t>::type>::type the_face_element_type;
    typedef typename the_face_element_type::super2::template Element<the_face_element_type>::type the_element_type;

    static inline const uint16_type nOrderGeo = the_element_type::nOrder;
#if 0
    static const bool is_polynomial = ExprT::imIsPoly;
#if 0
    static const int value = boost::mpl::if_< boost::mpl::bool_< ExprT::imIsPoly > ,
                     typename boost::mpl::if_< boost::mpl::greater< boost::mpl::int_<ExprT::imorder>,
                     boost::mpl::int_<19> > ,
                     boost::mpl::int_<19>,
                     boost::mpl::int_<ExprT::imorder> >::type,
                     boost::mpl::int_<10> >::type::value;
#else
    // this is a very rough approximation
    static inline const uint16_type value = ( ExprT::imorder )?( ExprT::imorder*nOrderGeo ):( nOrderGeo );
    static inline const uint16_type value_1 = ExprT::imorder+(the_element_type::is_hypercube?nOrderGeo:0);
#endif
#else
    static bool isPolynomial( ExprT const& expr ) { return expr.isPolynomial(); }
    static quad_order_type value( ExprT const& expr ) { return ( expr.polynomialOrder() )?( expr.polynomialOrder()*nOrderGeo ):( nOrderGeo ); }
    static quad_order_type value_1( ExprT const& expr ) { return expr.polynomialOrder()+(the_element_type::is_hypercube?nOrderGeo:0); }
#endif
    ExpressionOrder() = default;
};




#if 0
template < class Element, int Type>
class GElem
{
public:

    static const size_type context = vm::JACOBIAN |vm::POINT;
    static const bool is_terminal = false;

    typedef Element element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef GElem<element_type, Type> this_type;
    typedef this_type self_type;

    typedef typename element_type::functionspace_type functionspace_type;
    typedef typename functionspace_type::reference_element_type* fe_ptrtype;
    typedef typename functionspace_type::reference_element_type fe_type;
    typedef typename functionspace_type::geoelement_type geoelement_type;
    typedef typename functionspace_type::gm_type gm_type;
    typedef typename functionspace_type::value_type value_type;
    static inline const uint16_type rank = fe_type::rank;
    static inline const uint16_type nComponents1 = fe_type::nComponents1;
    static inline const uint16_type nComponents2 = fe_type::nComponents2;
    typedef std::map<size_type,std::vector<element_ptrtype> > basis_type;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ( Type==0 );
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ( Type==1 );
    };


    GElem ( std::map<size_type,std::vector<element_ptrtype> > const& v )
        :
        M_basis ( v )
    {
        typename basis_type::iterator it = M_basis.begin();
        typename basis_type::iterator en = M_basis.end();

        for ( ; it != en; ++it )
            for ( uint16_type i = 0; i < it->second.size(); ++i )
                it->second[i]->updateGlobalValues();
    }
    GElem( GElem const& op )
        :
        M_basis ( op.M_basis )
    {

    }

    basis_type const&  basis() const
    {
        return M_basis;
    }

    //! polynomial order
    constexpr uint16_type polynomialOrder() const { return element_type::functionspace_type::basis_type::nOrder; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        typedef typename mpl::if_<mpl::equal_to<mpl::int_<rank>,
                mpl::int_<0> >,
                mpl::identity<Shape<gmc_type::NDim, Scalar, false> >,
                typename mpl::if_<mpl::equal_to<mpl::int_<rank>,
                mpl::int_<1> >,
                mpl::identity<Shape<gmc_type::NDim, Vectorial, false> >,
                mpl::identity<Shape<gmc_type::NDim, Tensor2, false> > >::type>::type::type shape;
        typedef typename fe_type::PreCompute pc_type;
        typedef std::shared_ptr<pc_type> pc_ptrtype;
        typedef typename fe_type::template Context<context, fe_type, gm_type,geoelement_type,gmc_type::context> ctx_type;
        typedef std::shared_ptr<ctx_type> ctx_ptrtype;

        typedef typename expression_type::value_type value_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_expr( expr ),
            M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            M_ctx( new ctx_type( expr.basis().begin()->second[0]->functionSpace()->fe(),
                                 fusion::at_key<key_type>( geom ), ( pc_ptrtype const& )M_pc ) ),
            M_loc( expr.basis().begin()->second.size() )
        {
            for ( uint16_type i = 0; i < M_loc.size(); ++i )
                M_loc[i].resize( boost::extents[M_pc.nPoints()][nComponents1][nComponents2] );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            M_expr( expr ),
            M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            M_ctx( new ctx_type( expr.basis().begin()->second[0]->functionSpace()->fe(),
                                 fusion::at_key<key_type>( geom ), ( pc_ptrtype const& )M_pc ) ),
            M_loc( expr.basis().begin()->second.size() )
        {
            for ( uint16_type i = 0; i < M_loc.size(); ++i )
                M_loc[i].resize( boost::extents[M_pc.nPoints()] );
        }
        tensor( expression_type const& expr,
                Geo_t const& geom )
            :
            M_expr( expr ),
            M_pc( expr.basis().begin()->second[0]->functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ),
            M_ctx( new ctx_type( expr.basis().begin()->second[0]->functionSpace()->fe(),
                                 fusion::at_key<key_type>( geom ), ( pc_ptrtype const& )M_pc ) ),
            M_loc( expr.basis().begin()->second.size() )
        {
            for ( uint16_type i = 0; i < M_loc.size(); ++i )
                M_loc[i].resize( boost::extents[M_pc.nPoints()] );
        }
        template<typename IM>
        void init( IM const& /*im*/ )
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
            //VLOG(1) << "[GElem] updating element " << fusion::at_key<key_type>( geom )->id() << "\n";
            typename basis_type::iterator it = const_cast<basis_type&>( M_expr.basis() ).find( fusion::at_key<key_type>( geom )->id() );
            typename basis_type::iterator en = const_cast<basis_type&>( M_expr.basis() ).end();

            FEELPP_ASSERT( it != en )( fusion::at_key<key_type>( geom )->id() ).error ( "invalid basis function to integrate" );

            for ( uint16_type i = 0; i < M_loc.size(); ++i )
            {
                //M_loc[i] = it->second[i]->id( *M_ctx, M_pc, M_loc[i] );
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( Type == 0 )
                return M_loc[i]( c1,c2,q );

            return M_loc[j]( c1,c2,q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext> ) const
        {
            if ( Type == 0 )
                return M_loc[i]( c1,c2,q );

            return M_loc[j]( c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_loc[i]( c1,c2,q );
        }
        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {

        }
    private:
        this_type const& M_expr;
        pc_type M_pc;
        ctx_ptrtype M_ctx;
        std::vector<typename element_type::id_type> M_loc;
    };
private:
    basis_type M_basis;

};

template<typename Elem>
inline
Expr< GElem<Elem,1> >
basist( std::map<size_type,std::vector<std::shared_ptr<Elem> > > const& v )
{
    typedef GElem<Elem,1> expr_t;
    return Expr< expr_t >(  expr_t( v ) );
}
template<typename Elem>
inline
Expr< GElem<Elem,0> >
basis( std::map<size_type,std::vector<std::shared_ptr<Elem> > > const& v )
{
    typedef GElem<Elem,0> expr_t;
    return Expr< expr_t >(  expr_t( v ) );
}
#endif

/// \endcond
} // vf


using namespace vf;

} // feel
#endif /* FEELPP_EXPR_HPP */
