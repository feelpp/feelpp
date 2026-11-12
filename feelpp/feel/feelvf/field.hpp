/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELVF_FIELD_HPP
#define FEELPP_FEELVF_FIELD_HPP 1

#include <type_traits>
#include <utility>

#include <boost/ref.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/operators.hpp>

namespace Feel::vf
{

template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
class FieldValue;

/**
 * @brief Non-owning value terminal for a finite-element field.
 *
 * This terminal gives a function-space element first-class tensor support in
 * the variational DSEL without changing the element class itself. The element
 * remains the coefficient/state object, while this view provides the
 * `__VALUE` expression semantics used by `idv`, `gradv`, and integration.
 */
template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
class ElementValueTerminal : public OpId<ElementType, __VALUE>
{
private:
    using op_base_type = OpId<ElementType, __VALUE>;

public:
    using element_type = ElementType;
    using value_type = typename op_base_type::value_type;

    /**
     * @brief Build a value terminal referencing an existing element.
     */
    explicit ElementValueTerminal( element_type const& element, bool useInterpWithConfLoc = false )
        : op_base_type( element, useInterpWithConfLoc )
    {}

    /**
     * @brief Referenced finite-element field.
     */
    element_type const& element() const noexcept { return this->e(); }

    /**
     * @brief Tensor adaptor for the finite-element value terminal.
     */
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
        : public op_base_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>
    {
        using super_type = typename op_base_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        tensor( ElementValueTerminal const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            : super_type( static_cast<op_base_type const&>( expr ), geom, fev, feu )
        {}

        tensor( ElementValueTerminal const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            : super_type( static_cast<op_base_type const&>( expr ), geom, fev )
        {}

        tensor( ElementValueTerminal const& expr,
                Geo_t const& geom )
            : super_type( static_cast<op_base_type const&>( expr ), geom )
        {}

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType,
                  typename... TheArgsType>
        tensor( std::true_type,
                TheExprExpandedType const& exprExpanded,
                TupleTensorSymbolsExprType& ttse,
                ElementValueTerminal const& expr,
                Geo_t const& geom,
                TheArgsType const&... theInitArgs )
            : super_type( std::true_type{}, exprExpanded, ttse,
                          static_cast<op_base_type const&>( expr ), geom,
                          theInitArgs... )
        {}
    };
};

/**
 * @brief First-class value view of a finite-element field.
 *
 * This is the compatibility-safe spelling for value gradients:
 * `grad(value(uh))` is equivalent to `gradv(uh)`. Raw `grad(uh)` keeps its
 * existing Feel++ meaning.
 */
template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
class FieldValue : public Expr<ElementValueTerminal<ElementType>>
{
private:
    using value_expression_type = ElementValueTerminal<ElementType>;
    using expr_base_type = Expr<value_expression_type>;

public:
    using element_type = ElementType;
    using value_type = typename expr_base_type::value_type;

    /**
     * @brief Build a value expression view of a finite-element field.
     */
    explicit FieldValue( element_type const& element, bool useInterpWithConfLoc = false )
        : expr_base_type( value_expression_type( element, useInterpWithConfLoc ) ),
          M_element( boost::cref( element ) )
    {}

    /**
     * @brief Field element referenced by this value expression.
     */
    element_type const& element() const noexcept { return M_element.get(); }

    /**
     * @brief Access this object as the underlying Feel++ expression.
     */
    expr_base_type const& expressionBase() const noexcept { return *this; }

private:
    boost::reference_wrapper<element_type const> M_element;
};

namespace detail
{
template <typename T>
struct is_field_value : std::false_type {};

template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
struct is_field_value<FieldValue<ElementType>> : std::true_type {};
} // namespace detail

template <typename T>
inline constexpr bool is_field_value_v =
    detail::is_field_value<std::remove_cvref_t<T>>::value;

/**
 * @brief Create a first-class value expression from a finite-element field.
 */
template <typename ElementType>
    requires ( Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>> &&
               !is_field_value_v<ElementType> )
[[nodiscard]] inline FieldValue<ElementType>
value( ElementType const& element, bool useInterpWithConfLoc = false )
{
    return FieldValue<ElementType>( element, useInterpWithConfLoc );
}

/**
 * @brief Return the expression associated with a value-field object.
 */
template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
[[nodiscard]] inline auto const&
idv( FieldValue<ElementType> const& field )
{
    return field.expressionBase();
}

/**
 * @brief Gradient of a first-class field value.
 */
template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
[[nodiscard]] inline auto
grad( FieldValue<ElementType> const& field )
{
    return gradv( field.element() );
}

/**
 * @brief Divergence of a first-class field value.
 */
template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
[[nodiscard]] inline auto
div( FieldValue<ElementType> const& field )
{
    return divv( field.element() );
}

/**
 * @brief Curl of a first-class field value.
 */
template <typename ElementType>
    requires Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>>
[[nodiscard]] inline auto
curl( FieldValue<ElementType> const& field )
{
    return curlv( field.element() );
}

/**
 * @brief Promote finite-element fields to value expressions where unambiguous.
 *
 * DSEL entry points use this helper to accept `_expr = uh` as `_expr =
 * value(uh)`. Differential operators stay explicit: use `grad(value(uh))` or
 * the supported `gradv(uh)`.
 */
template <typename ExprType>
[[nodiscard]] inline decltype( auto )
toValueExpression( ExprType&& expr )
{
    if constexpr ( Feel::is_functionspace_element_v<std::remove_cvref_t<ExprType>> &&
                   !is_field_value_v<ExprType> )
        return value( expr );
    else
        return std::forward<ExprType>( expr );
}

#define FEELPP_VF_FIELD_BINARY_OPERATOR(OP)                                            \
    template <typename ElementType, typename RhsType>                                  \
        requires ( Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>> && \
                   !is_field_value_v<ElementType> )                                    \
    [[nodiscard]] inline auto operator OP( ElementType const& element, RhsType const& rhs ) \
    {                                                                                  \
        return value( element ) OP rhs;                                                \
    }                                                                                  \
    template <typename LhsType, typename ElementType>                                  \
        requires ( Feel::is_functionspace_element_v<std::remove_cvref_t<ElementType>> && \
                   !is_field_value_v<ElementType> &&                                   \
                   !( Feel::is_functionspace_element_v<std::remove_cvref_t<LhsType>> && \
                      !is_field_value_v<LhsType> ) )                                   \
    [[nodiscard]] inline auto operator OP( LhsType const& lhs, ElementType const& element ) \
    {                                                                                  \
        return lhs OP value( element );                                                \
    }                                                                                  \
    /**/

FEELPP_VF_FIELD_BINARY_OPERATOR( + )
FEELPP_VF_FIELD_BINARY_OPERATOR( - )
FEELPP_VF_FIELD_BINARY_OPERATOR( * )
FEELPP_VF_FIELD_BINARY_OPERATOR( / )

#undef FEELPP_VF_FIELD_BINARY_OPERATOR

} // namespace Feel::vf

#endif /* FEELPP_FEELVF_FIELD_HPP */
