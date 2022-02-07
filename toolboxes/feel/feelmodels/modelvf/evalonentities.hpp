/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_MODELS_VF_EVALONFENTITIES_H
#define FEELPP_MODELS_VF_EVALONFENTITIES_H

#include <feel/feelvf/expr.hpp>


namespace Feel
{
namespace vf
{

template<typename GmcT>
fusion::map<fusion::pair<vf::detail::gmc<0>, std::shared_ptr<GmcT>>>
mapgmcFix( std::shared_ptr<GmcT> const& ctx )
{
    return { fusion::make_pair<vf::detail::gmc<0> >( ctx ) };
}




template <typename ExprType>
class EvalOnFaces
{
public :
    using this_type = EvalOnFaces<ExprType>;
    using expression_type = ExprType;

    enum class InternalFacesEvalType { Mean=0,Max,Sum,One_Side };

    static const size_type context = expression_type::context;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = expression_type::template HasTestFunction<Func>::result;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = expression_type::template HasTrialFunction<Func>::result;
    };
    template<typename Funct>
    static const bool has_test_basis = expression_type::template has_test_basis<Funct>;
    template<typename Funct>
    static const bool has_trial_basis = expression_type::template has_trial_basis<Funct>;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;

    using value_type = typename expression_type::value_type;
    using evaluate_type = typename expression_type::evaluate_type;

#if 0
    template <typename TheExprType>
    EvalOnFaces( TheExprType && e )
        :
        M_expr( std::forward<TheExprType>(e) )
        {}
#else
    EvalOnFaces( expression_type const& e, std::string const& type, std::set<std::string> const& markers )
        :
        M_expr( e ),
        M_type( InternalFacesEvalType::One_Side )
        {
            if ( type == "mean" )
                M_type = InternalFacesEvalType::Mean;
            else if ( type == "max" )
                M_type = InternalFacesEvalType::Max;
            else if ( type == "sum" )
                M_type = InternalFacesEvalType::Sum;
            else if ( type == "one_side" )
                M_type = InternalFacesEvalType::One_Side;
 
            // if ( M_type == EvalType::Markers_Connection )
            // {
            //     M_markersConnection = markers;
            //     CHECK( !M_markersConnection.empty() ) << "no markers";
            // }
        }
#endif
    EvalOnFaces( EvalOnFaces const& ) = default;
    EvalOnFaces( EvalOnFaces && ) = default;


    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    size_type dynamicContext() const { return M_expr.dynamicContext(); }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newexpr = M_expr.applySymbolsExpr( se );
            using newexpr_type = std::decay_t<decltype(newexpr)>;
            return EvalOnFaces<newexpr_type>( newexpr );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const { return M_expr.hasSymbolDependency( symb,se ); }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_expr.dependentSymbols( symb, res, se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            CHECK( false ) << "TODO";
            return *this;
        }


    expression_type const& expr() const { return M_expr; }
    InternalFacesEvalType type() const { return M_type; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
        static constexpr bool has_two_side = map_size::value == 2;
        using gmcKey0 = vf::detail::gmc<0>;
        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
               vf::detail::gmc<1>,
               vf::detail::gmc<0> >::type gmcKey1;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmcKey0>::type left_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmcKey0>::type::element_type left_gmc_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmcKey1>::type right_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmcKey1>::type::element_type right_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, left_gmc_ptrtype> > map_left_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, right_gmc_ptrtype> > map_right_gmc_type;
        typedef typename expression_type::template tensor<map_left_gmc_type, Basis_i_t, Basis_j_t> left_tensor_expr_type;
        typedef typename expression_type::template tensor<map_right_gmc_type, Basis_i_t, Basis_j_t> right_tensor_expr_type;

        using value_type = typename left_tensor_expr_type::value_type;
        using expr_shape = typename left_tensor_expr_type::shape;
        using shape = expr_shape;

        struct is_zero
        {
            static const bool value = false;
        };

        template <typename ... TheArgsType>
        tensor( this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_expr( expr ),
            M_useLeft( false ), M_useRight( false ),
            M_type( expr.type() )
            {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_expr( expr ),
            M_useLeft( false ), M_useRight( false ),
            M_type( expr.type() )
            {
                CHECK( false ) << "TODO";
            }

        template<typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        // template <typename ... TheArgsType>
        // void update( Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        // {
        // }
#if 1
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            CHECK( false ) << "TODO";
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            CHECK( false ) << "TODO";
        }
        void update( Geo_t const& geom )
        {
            M_useLeft = true;
            M_useRight = has_two_side;
#if 0
            if ( M_type == EvalType::Markers_Connection )
                upMarkerId; // only if mesh changed or not init
#endif
            if constexpr( has_two_side )
            {
                if ( M_type == InternalFacesEvalType::One_Side )
                    M_useRight = false;
#if 0
                else if ( M_type == EvalType::Markers_Connection )
                {
                    auto leftGmc = fusion::at_key<gmcKey0>( geom );
                    auto const& leftElement = leftGmc->element();

                    //leftElement.marker()
                    leftElement.mesh()->markerName( leftElement.marker() );
                }
#endif
            }

            if ( M_useLeft )
            {
                auto leftGeom = Feel::vf::mapgmcFix( fusion::at_key<gmcKey0>( geom ) );
                if ( !M_leftTensor )
                    M_leftTensor.emplace( M_expr.expr(), leftGeom );
                M_leftTensor->update( leftGeom );
            }

            if constexpr( has_two_side )
            {
                if ( M_useRight )
                {
                    auto rightGeom = Feel::vf::mapgmcFix( fusion::at_key<gmcKey1>( geom ) );
                    if ( !M_rightTensor )
                        M_rightTensor.emplace( M_expr.expr(), rightGeom );
                    M_rightTensor->update( rightGeom );
                }
            }
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            this->update( geom );
        }
#endif
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            CHECK( false ) << "TODO";
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            CHECK( false ) << "TODO";
            value_type res(0);
            return res;
        }
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                CHECK( false ) << "TODO";
                value_type res(0);
                return res;
            }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                value_type res(0);

                value_type evalLeft = M_useLeft ? M_leftTensor->evalq( c1,c2,q ) : 0;
                value_type evalRight = M_useRight ? M_rightTensor->evalq( c1,c2,q ) : 0;
                switch ( M_type )
                {
                case InternalFacesEvalType::Mean:
                    res = 0.5*(evalLeft+evalRight);
                    break;
                case InternalFacesEvalType::Max:
                    res = std::max(evalLeft,evalRight);
                    break;
                default: // all others cases can be take into account here
                    res = evalLeft+evalRight;
                    break;
                }
                return res;
            }
    private:
        this_type const& M_expr;
        InternalFacesEvalType M_type;
        std::optional<left_tensor_expr_type> M_leftTensor;
        std::optional<right_tensor_expr_type> M_rightTensor;
        bool M_useLeft, M_useRight;
        //map_left_gmc_type M_leftGeom;
        //map_right_gmc_type M_rightGeom;
        //std::set<flag_type> M_markersId;

    };

private :
    expression_type M_expr;
    InternalFacesEvalType M_type;
    std::set<std::string> M_markersConnection;
};


template <typename ExprType>
auto evalOnFaces( ExprType && e, std::set<std::string> const& requiresMarkersConnection, std::string const& internalFacesEvalutationType = "" )
{
    return Feel::vf::expr( EvalOnFaces<std::decay_t<ExprType>>( std::forward<ExprType>(e), internalFacesEvalutationType, requiresMarkersConnection ) );
}


template <typename RangeType,typename ExprType>
decltype(auto) evalOnEntities( RangeType const& range, ExprType && e, std::set<std::string> const& requiresMarkersConnection, std::string const& internalFacesEvalutationType = "" )
{
    if constexpr( filter_enum_t<RangeType>::value == MESH_FACES )
        return evalOnFaces( std::forward<ExprType>( e ), requiresMarkersConnection, internalFacesEvalutationType );
    else
        return std::forward<ExprType>( e );
}

} // namespace vf
} // namespace Feel

#endif
