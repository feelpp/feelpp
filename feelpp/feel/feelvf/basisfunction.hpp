/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELVF_BASISFUNCTION_HPP
#define FEELPP_FEELVF_BASISFUNCTION_HPP 1

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/operators.hpp>

namespace Feel::vf
{

/**
 * @brief Role carried by a first-class variational basis object.
 */
enum class basis_role
{
    test,
    trial
};

template <typename SpaceType, basis_role Role>
class BasisFunction;

/**
 * @brief Owning variant of a supported Feel++ basis operator terminal.
 *
 * The existing `id`, `idt`, `grad`, `gradt`, and related operators store a
 * reference to a finite-element object. That is correct when the element is
 * owned by user code. First-class `test(Xh)`/`trial(Xh)` create their own
 * internal element, so this wrapper keeps that element alive inside copied
 * expression trees and provides the tensor terminal used directly by the
 * variational DSEL.
 */
template <template <typename, OperatorType> class OperatorTemplate,
          typename ElementType, OperatorType Type>
class OpOwned;

namespace detail
{
template <typename ElementType>
class OpOwnedStorage
{
public:
    using element_type = ElementType;
    using element_ptrtype = std::shared_ptr<element_type>;

    explicit OpOwnedStorage( element_ptrtype element )
        : M_element( std::move( element ) )
    {
        if ( !M_element )
            throw std::invalid_argument( "Feel::vf::OpOwned requires a non-null element" );
    }

    element_type const& element() const noexcept { return *M_element; }
    element_ptrtype const& elementPtr() const noexcept { return M_element; }

private:
    element_ptrtype M_element;
};

template <typename T>
struct is_basis_function : std::false_type {};

template <typename SpaceType, basis_role Role>
struct is_basis_function<BasisFunction<SpaceType, Role>> : std::true_type {};

template <typename T, basis_role Role>
struct is_basis_role : std::false_type {};

template <typename SpaceType, basis_role BasisRole, basis_role Role>
struct is_basis_role<BasisFunction<SpaceType, BasisRole>, Role>
    : std::bool_constant<BasisRole == Role> {};

template <typename SpaceType, basis_role Role>
class BasisFunctionStorage
{
public:
    using space_type = SpaceType;
    using space_ptrtype = std::shared_ptr<space_type>;
    using element_type = typename space_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using basis_type = typename space_type::basis_type;

    /**
     * @brief Build the storage from a function space.
     */
    explicit BasisFunctionStorage( space_ptrtype space,
                                   std::string name = {}, std::string desc = {} )
        : M_space( std::move( space ) )
    {
        if ( !M_space )
            throw std::invalid_argument( "Feel::vf::BasisFunction requires a non-null function space" );
        M_element = std::make_shared<element_type>( M_space->element( std::move( name ), std::move( desc ) ) );
    }

    /**
     * @brief Function space associated with the basis object.
     */
    space_ptrtype const& space() const noexcept { return M_space; }

    /**
     * @brief Internal element used by the existing Feel++ expression engine.
     */
    element_type const& element() const noexcept { return *M_element; }

    /**
     * @brief Internal element handle kept alive by the basis object.
     */
    element_ptrtype const& elementPtr() const noexcept { return M_element; }

    /**
     * @brief Polynomial basis carried by the function space.
     */
    basis_type const& basis() const { return M_space->fe(); }

    /**
     * @brief Mesh support of the function space.
     */
    auto support() const { return Feel::support( M_space ); }

private:
    space_ptrtype M_space;
    element_ptrtype M_element;
};

template <typename SpaceType, basis_role Role>
struct basis_expr
{
    using element_type = typename SpaceType::element_type;
    using type = std::conditional_t<Role == basis_role::trial,
                                    OpOwned<OpId, element_type, __TRIAL>,
                                    OpOwned<OpId, element_type, __TEST>>;
};
} // namespace detail

template <template <typename, OperatorType> class OperatorTemplate,
          typename ElementType, OperatorType Type>
class OpOwned
    : private detail::OpOwnedStorage<ElementType>,
      public OperatorTemplate<ElementType, Type>
{
private:
    using storage_type = detail::OpOwnedStorage<ElementType>;
    using op_base_type = OperatorTemplate<ElementType, Type>;

public:
    using element_type = ElementType;
    using element_ptrtype = typename storage_type::element_ptrtype;

    /**
     * @brief Build an owning terminal from a shared finite-element object.
     */
    explicit OpOwned( element_ptrtype element, bool useInterpWithConfLoc = false )
        : storage_type( std::move( element ) ),
          op_base_type( storage_type::element(), useInterpWithConfLoc )
    {}

    /**
     * @brief Copy while rebinding the supported operator to the owned element.
     */
    OpOwned( OpOwned const& other )
        : storage_type( other.elementPtr() ),
          op_base_type( storage_type::element(), other.useInterpWithConfLoc() )
    {}

    /**
     * @brief Move while rebinding the supported operator to the moved element.
     */
    OpOwned( OpOwned&& other ) noexcept
        : storage_type( std::move( static_cast<storage_type&>( other ) ) ),
          op_base_type( storage_type::element(), other.useInterpWithConfLoc() )
    {}

    /**
     * @brief Assign and rebind the supported operator to the owned element.
     */
    OpOwned& operator=( OpOwned const& other )
    {
        static_cast<storage_type&>( *this ) = static_cast<storage_type const&>( other );
        static_cast<op_base_type&>( *this ) = op_base_type( storage_type::element(), other.useInterpWithConfLoc() );
        return *this;
    }

    /**
     * @brief Move-assign and rebind the supported operator to the owned element.
     */
    OpOwned& operator=( OpOwned&& other ) noexcept
    {
        bool const useInterpWithConfLoc = other.useInterpWithConfLoc();
        static_cast<storage_type&>( *this ) = std::move( static_cast<storage_type&>( other ) );
        static_cast<op_base_type&>( *this ) = op_base_type( storage_type::element(), useInterpWithConfLoc );
        return *this;
    }

    /**
     * @brief Internal element handle owned by this terminal.
     */
    element_ptrtype const& elementPtr() const noexcept { return storage_type::elementPtr(); }

    /**
     * @brief Tensor adaptor for the owned basis operator terminal.
     *
     * This is the point where `test(Xh)`, `trial(Xh)`, and their derivatives
     * become real Feel++ expression terminals. The implementation reuses the
     * existing operator tensor evaluator but binds it to the owned element
     * stored by this expression node.
     */
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
        : public op_base_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>
    {
        using super_type = typename op_base_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        tensor( OpOwned const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            : super_type( static_cast<op_base_type const&>( expr ), geom, fev, feu )
        {}

        tensor( OpOwned const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            : super_type( static_cast<op_base_type const&>( expr ), geom, fev )
        {}

        tensor( OpOwned const& expr,
                Geo_t const& geom )
            : super_type( static_cast<op_base_type const&>( expr ), geom )
        {}

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType,
                  typename... TheArgsType>
        tensor( std::true_type,
                TheExprExpandedType const& exprExpanded,
                TupleTensorSymbolsExprType& ttse,
                OpOwned const& expr,
                Geo_t const& geom,
                TheArgsType const&... theInitArgs )
            : super_type( std::true_type{}, exprExpanded, ttse,
                          static_cast<op_base_type const&>( expr ), geom,
                          theInitArgs... )
        {}
    };
};

template <typename ElementType, OperatorType Type>
using OpOwnedId = OpOwned<OpId, ElementType, Type>;

namespace detail
{
template <typename BasisType>
inline constexpr OperatorType basis_operator_type_v =
    std::remove_cvref_t<BasisType>::role_v == basis_role::trial ? __TRIAL : __TEST;

template <template <typename, OperatorType> class OperatorTemplate, typename BasisType>
[[nodiscard]] inline auto
basisOperatorExpression( BasisType const& basis )
{
    using basis_type = std::remove_cvref_t<BasisType>;
    static_assert( is_basis_function<basis_type>::value,
                   "basisOperatorExpression requires a Feel::vf::BasisFunction" );

    using expression_type = OpOwned<OperatorTemplate,
                                    typename basis_type::element_type,
                                    basis_operator_type_v<basis_type>>;
    return Expr<expression_type>( expression_type( basis.elementPtr() ) );
}
} // namespace detail

/**
 * @brief First-class test or trial basis function for a Feel++ function space.
 *
 * `BasisFunction` is a real DSEL expression terminal backed by an owning
 * operator tensor. It carries the function space, polynomial basis, support,
 * and mathematical role explicitly while interoperating with the supported
 * `id`, `idt`, `grad`, `gradt`, `idv`, and `gradv` API.
 */
template <typename SpaceType, basis_role Role>
class BasisFunction
    : private detail::BasisFunctionStorage<SpaceType, Role>,
      public Expr<typename detail::basis_expr<SpaceType, Role>::type>
{
private:
    using storage_type = detail::BasisFunctionStorage<SpaceType, Role>;
    using basis_expression_type = typename detail::basis_expr<SpaceType, Role>::type;
    using expr_base_type = Expr<basis_expression_type>;

public:
    using space_type = typename storage_type::space_type;
    using space_ptrtype = typename storage_type::space_ptrtype;
    using element_type = typename storage_type::element_type;
    using element_ptrtype = typename storage_type::element_ptrtype;
    using basis_type = typename storage_type::basis_type;
    using value_type = typename expr_base_type::value_type;
    static constexpr basis_role role_v = Role;

    /**
     * @brief Build a basis object from a function space.
     */
    explicit BasisFunction( space_ptrtype const& space,
                            std::string name = {}, std::string desc = {} )
        : storage_type( space, std::move( name ), std::move( desc ) ),
          expr_base_type( basis_expression_type( storage_type::elementPtr() ) )
    {}

    /**
     * @brief Copy constructor rebinding the expression to the copied storage.
     */
    BasisFunction( BasisFunction const& other )
        : storage_type( other ),
          expr_base_type( basis_expression_type( storage_type::elementPtr() ) )
    {}

    /**
     * @brief Move constructor rebinding the expression to the moved storage.
     */
    BasisFunction( BasisFunction&& other ) noexcept
        : storage_type( std::move( other ) ),
          expr_base_type( basis_expression_type( storage_type::elementPtr() ) )
    {}

    /**
     * @brief Copy assignment rebinding the expression to the assigned storage.
     */
    BasisFunction& operator=( BasisFunction const& other )
    {
        static_cast<storage_type&>( *this ) = static_cast<storage_type const&>( other );
        this->expression() = basis_expression_type( storage_type::elementPtr() );
        return *this;
    }

    /**
     * @brief Move assignment rebinding the expression to the assigned storage.
     */
    BasisFunction& operator=( BasisFunction&& other ) noexcept
    {
        static_cast<storage_type&>( *this ) = std::move( static_cast<storage_type&>( other ) );
        this->expression() = basis_expression_type( storage_type::elementPtr() );
        return *this;
    }

    /**
     * @brief Function space associated with the basis object.
     */
    space_ptrtype const& space() const noexcept { return storage_type::space(); }

    /**
     * @brief Internal element used by the existing expression engine.
     */
    element_type const& element() const noexcept { return storage_type::element(); }

    /**
     * @brief Internal element handle kept alive by this basis object.
     */
    element_ptrtype const& elementPtr() const noexcept { return storage_type::elementPtr(); }

    /**
     * @brief Polynomial basis carried by the function space.
     */
    basis_type const& basis() const { return storage_type::basis(); }

    /**
     * @brief Mesh support of the function space.
     */
    auto support() const { return storage_type::support(); }

    /**
     * @brief Mathematical role of this basis object.
     */
    constexpr basis_role role() const noexcept { return Role; }

    /**
     * @brief Human-readable role name.
     */
    constexpr char const* roleName() const noexcept
    {
        if constexpr ( Role == basis_role::trial )
            return "trial";
        else
            return "test";
    }

    /**
     * @brief Access this object as the underlying Feel++ expression.
     */
    expr_base_type const& expressionBase() const noexcept { return *this; }
};

template <typename SpaceType>
using TestBasis = BasisFunction<SpaceType, basis_role::test>;

template <typename SpaceType>
using TrialBasis = BasisFunction<SpaceType, basis_role::trial>;

/**
 * @brief Create a first-class test basis function for `Xh`.
 */
template <typename SpaceType>
[[nodiscard]] inline TestBasis<SpaceType>
test( std::shared_ptr<SpaceType> const& Xh )
{
    return TestBasis<SpaceType>( Xh, "v", "v" );
}

/**
 * @brief Create a first-class trial basis function for `Xh`.
 */
template <typename SpaceType>
[[nodiscard]] inline TrialBasis<SpaceType>
trial( std::shared_ptr<SpaceType> const& Xh )
{
    return TrialBasis<SpaceType>( Xh, "u", "u" );
}

/** @brief Create an owning basis while preserving the named-space API. */
template <typename SpaceType>
[[nodiscard]] inline TestBasis<SpaceType>
test( std::shared_ptr<SpaceType> const& Xh, std::string name, std::string desc = "v" )
{
    return TestBasis<SpaceType>( Xh, std::move( name ), std::move( desc ) );
}

/** @brief Create an owning basis while preserving the named-space API. */
template <typename SpaceType>
[[nodiscard]] inline TrialBasis<SpaceType>
trial( std::shared_ptr<SpaceType> const& Xh, std::string name, std::string desc = "u" )
{
    return TrialBasis<SpaceType>( Xh, std::move( name ), std::move( desc ) );
}

/**
 * @brief Return the expression associated with a test basis object.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
id( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpId>( basis );
}

/**
 * @brief Return the expression associated with a trial basis object.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
idt( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpId>( basis );
}

/** @brief Preserve role-aware id() for trial bases with owning storage. */
template <typename SpaceType>
[[nodiscard]] inline auto
id( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpId>( basis );
}

/**
 * @brief Boundary normal component of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
normal( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpN>( basis );
}

/**
 * @brief Boundary normal component of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
normal( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpN>( basis );
}

/**
 * @brief x-derivative of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dx( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDx>( basis );
}

/**
 * @brief x-derivative of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dx( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDx>( basis );
}

/**
 * @brief y-derivative of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dy( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDy>( basis );
}

/**
 * @brief y-derivative of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dy( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDy>( basis );
}

/**
 * @brief z-derivative of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dz( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDz>( basis );
}

/**
 * @brief z-derivative of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dz( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDz>( basis );
}

/**
 * @brief Normal derivative of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dn( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDn>( basis );
}

/**
 * @brief Normal derivative of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
dn( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDn>( basis );
}

/**
 * @brief Gradient of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
grad( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpGrad>( basis );
}

/**
 * @brief Gradient of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
grad( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpGrad>( basis );
}

/**
 * @brief Symmetric gradient of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
symm_grad( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpSymmGrad>( basis );
}

/**
 * @brief Symmetric gradient of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
symm_grad( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpSymmGrad>( basis );
}

/**
 * @brief Divergence of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
div( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDiv>( basis );
}

/**
 * @brief Divergence of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
div( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpDiv>( basis );
}

/**
 * @brief Curl of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curl( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurl>( basis );
}

/**
 * @brief Curl of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curl( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurl>( basis );
}

/**
 * @brief x-component of the curl of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curlx( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurlX>( basis );
}

/**
 * @brief x-component of the curl of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curlx( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurlX>( basis );
}

/**
 * @brief y-component of the curl of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curly( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurlY>( basis );
}

/**
 * @brief y-component of the curl of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curly( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurlY>( basis );
}

/**
 * @brief z-component of the curl of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curlz( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurlZ>( basis );
}

/**
 * @brief z-component of the curl of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
curlz( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpCurlZ>( basis );
}

/**
 * @brief Hessian of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
hess( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpHess>( basis );
}

/**
 * @brief Hessian of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
hess( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpHess>( basis );
}

/**
 * @brief Laplacian of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
laplacian( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpLap>( basis );
}

/**
 * @brief Laplacian of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
laplacian( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpLap>( basis );
}

/**
 * @brief Trace of a first-class test basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
trace( TestBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpTrace>( basis );
}

/**
 * @brief Trace of a first-class trial basis function.
 */
template <typename SpaceType>
[[nodiscard]] inline auto
trace( TrialBasis<SpaceType> const& basis )
{
    return detail::basisOperatorExpression<OpTrace>( basis );
}

template <typename T>
inline constexpr bool is_basis_function_v =
    detail::is_basis_function<std::remove_cvref_t<T>>::value;

template <typename T>
inline constexpr bool is_test_basis_v =
    detail::is_basis_role<std::remove_cvref_t<T>, basis_role::test>::value;

template <typename T>
inline constexpr bool is_trial_basis_v =
    detail::is_basis_role<std::remove_cvref_t<T>, basis_role::trial>::value;

} // namespace Feel::vf

#endif /* FEELPP_FEELVF_BASISFUNCTION_HPP */
