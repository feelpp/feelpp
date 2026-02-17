/**
 * @file modernization_examples.hpp
 * @brief Before/After examples of C++20 concept modernization
 * 
 * This file demonstrates concrete improvements from replacing SFINAE
 * with C++20 concepts in Feel++ code.
 */

#ifndef FEELPP_MODERNIZATION_EXAMPLES_HPP
#define FEELPP_MODERNIZATION_EXAMPLES_HPP

#include <type_traits>
#include <memory>

namespace Feel
{
namespace Examples
{

//=============================================================================
// Example 1: Simple Type Constraint
//=============================================================================

// BEFORE (C++17 SFINAE - verbose, unclear)
//-----------------------------------------------------------------------------
template <typename T, typename = void>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),
                                   decltype(std::declval<T>().end())>>
    : std::true_type {};

template <typename T>
inline constexpr bool is_iterable_v = is_iterable<T>::value;

// Usage with SFINAE
template <typename Container>
typename std::enable_if_t<is_iterable_v<Container>, void>
process_old(Container const& c)
{
    // Error here is 50+ lines of template instantiation stack
    for (auto const& item : c) { }
}

// AFTER (C++20 Concepts - clear, concise)
//-----------------------------------------------------------------------------
#if __cplusplus >= 202002L
template <typename T>
concept Iterable = requires(T t) {
    { t.begin() } -> std::input_or_output_iterator;
    { t.end() } -> std::sentinel_for<decltype(t.begin())>;
};

// Usage with concept
template <Iterable Container>
void process_new(Container const& c)
{
    // Error here: "constraint not satisfied: Container does not satisfy Iterable"
    // Clear, points to exact problem
    for (auto const& item : c) { }
}
#endif

//=============================================================================
// Example 2: Multiple Constraints (Common in Feel++)
//=============================================================================

// BEFORE - Unreadable mess
//-----------------------------------------------------------------------------
template <typename T, 
          typename = std::enable_if_t<std::is_base_of_v<ExprBase, T>>,
          typename = std::enable_if_t<std::is_arithmetic_v<typename T::value_type>>,
          typename = std::enable_if_t<T::rank == 0>>
auto compute_norm_old(T const& expr)
{
    return /* ... */;
}

// AFTER - Crystal clear
//-----------------------------------------------------------------------------
#if __cplusplus >= 202002L
template <typename T>
concept ScalarExpr = requires {
    requires std::derived_from<T, ExprBase>;
    typename T::value_type;
    requires std::is_arithmetic_v<typename T::value_type>;
    requires T::rank == 0;
};

template <ScalarExpr T>
auto compute_norm_new(T const& expr)
{
    return /* ... */;
}
#endif

//=============================================================================
// Example 3: Overload Resolution (blockforms.hpp pattern)
//=============================================================================

class ProductSpacesBase {};
class ProductSpaceBase {};

// BEFORE - 40+ occurrences in blockforms.hpp
//-----------------------------------------------------------------------------
template <typename T>
class BlockBilinearFormOld
{
public:
    // Constructor 1: For ProductSpaces
    template <typename PS>
    BlockBilinearFormOld(PS&& ps, 
                        std::enable_if_t<std::is_base_of<ProductSpacesBase, 
                                                         std::decay_t<PS>>::value>* = nullptr)
    { /* plural */ }

    // Constructor 2: For ProductSpace (singular)
    template <typename PS>
    BlockBilinearFormOld(PS&& ps,
                        std::enable_if_t<std::is_base_of<ProductSpaceBase,
                                                         std::decay_t<PS>>::value>* = nullptr)
    { /* singular */ }
    
    // Error: "no matching function for call to BlockBilinearForm"
    // (followed by 200 lines of template candidates)
};

// AFTER - Elegant, self-documenting
//-----------------------------------------------------------------------------
#if __cplusplus >= 202002L
template <typename T>
concept ProductSpaces = std::derived_from<std::decay_t<T>, ProductSpacesBase>;

template <typename T>
concept ProductSpace = std::derived_from<std::decay_t<T>, ProductSpaceBase>;

template <typename T>
class BlockBilinearFormNew
{
public:
    // Constructor 1: For ProductSpaces (clear from signature!)
    template <ProductSpaces PS>
    BlockBilinearFormNew(PS&& ps)
    { /* plural */ }

    // Constructor 2: For ProductSpace
    template <ProductSpace PS>
    BlockBilinearFormNew(PS&& ps)
    { /* singular */ }
    
    // Error: "constraint not satisfied: T does not satisfy ProductSpaces or ProductSpace"
    // (concise, helpful)
};
#endif

//=============================================================================
// Example 4: Dimension-based Dispatch (integrator.hpp pattern)
//=============================================================================

constexpr int MESH_ELEMENTS = 0;
constexpr int MESH_FACES = 1;
constexpr int MESH_POINTS = 2;

// BEFORE - Ugly template tricks
//-----------------------------------------------------------------------------
template <int iDim>
class IntegratorOld
{
public:
    // Elements version
    template <int iDimDummy = iDim, 
              std::enable_if_t<iDimDummy == MESH_ELEMENTS, bool> = true>
    void integrate_impl()
    { /* element integration */ }

    // Faces version
    template <int iDimDummy = iDim,
              std::enable_if_t<iDimDummy == MESH_FACES, bool> = true>
    void integrate_impl()
    { /* face integration */ }

    // Error: ambiguous overload, incomprehensible template error dump
};

// AFTER - Readable and maintainable
//-----------------------------------------------------------------------------
#if __cplusplus >= 202002L
template <int Dim>
concept MeshElements = (Dim == MESH_ELEMENTS);

template <int Dim>
concept MeshFaces = (Dim == MESH_FACES);

template <int Dim>
concept MeshPoints = (Dim == MESH_POINTS);

template <int iDim>
class IntegratorNew
{
public:
    // Elements version - clear from signature
    void integrate_impl() requires MeshElements<iDim>
    { /* element integration */ }

    // Faces version - clear from signature
    void integrate_impl() requires MeshFaces<iDim>
    { /* face integration */ }

    // Error: "no matching overload for integrate_impl"
    //        "constraint not satisfied: iDim=3 does not satisfy any variant"
    // (helpful, suggests valid values)
};
#endif

//=============================================================================
// Example 5: VF Expression Constraints (200+ uses in feelvf/)
//=============================================================================

class ExprBase {};

// BEFORE - Pointer parameter trick for SFINAE
//-----------------------------------------------------------------------------
template <typename ExprT>
auto expr_old(ExprT const& e, 
              std::enable_if_t<std::is_base_of_v<ExprBase, ExprT>>* = nullptr)
{
    return /* ... */;
}

template <typename ExprT>
auto trace_old(ExprT const& v,
               std::enable_if_t<std::is_base_of_v<ExprBase, ExprT>>* = nullptr)
{
    return /* ... */;
}

template <typename ExprT>
auto vonmises_old(ExprT const& v,
                  std::enable_if_t<std::is_base_of_v<ExprBase, ExprT>>* = nullptr)
{
    return /* ... */;
}

// AFTER - Consistent, clear pattern
//-----------------------------------------------------------------------------
#if __cplusplus >= 202002L
template <typename T>
concept VfExpr = std::derived_from<T, ExprBase>;

template <VfExpr ExprT>
auto expr_new(ExprT const& e)
{
    return /* ... */;
}

template <VfExpr ExprT>
auto trace_new(ExprT const& v)
{
    return /* ... */;
}

template <VfExpr ExprT>
auto vonmises_new(ExprT const& v)
{
    return /* ... */;
}
#endif

//=============================================================================
// Example 6: Quadrature Overload Resolution
//=============================================================================

// BEFORE - Complex SFINAE chains
//-----------------------------------------------------------------------------
template <typename QuadType, typename Quad1Type, typename ExprType>
auto integrate_old(
    QuadType const& quad,
    Quad1Type const& quad1,
    ExprType const& expr,
    std::enable_if_t<std::is_integral<QuadType>::value && 
                     std::is_integral<Quad1Type>::value>* = nullptr)
{
    return /* both are integral orders */;
}

template <typename QuadType, typename Quad1Type, typename ExprType>
auto integrate_old(
    QuadType const& quad,
    Quad1Type const& quad1,
    ExprType const& expr,
    std::enable_if_t<std::is_integral<QuadType>::value && 
                     !std::is_integral<Quad1Type>::value>* = nullptr)
{
    return /* quad is order, quad1 is quadrature object */;
}

// ... 2 more overloads with different combinations

// AFTER - Self-documenting
//-----------------------------------------------------------------------------
#if __cplusplus >= 202002L
template <typename QuadType, typename Quad1Type, typename ExprType>
requires std::integral<QuadType> && std::integral<Quad1Type>
auto integrate_new(QuadType quad, Quad1Type quad1, ExprType const& expr)
{
    return /* both are integral orders */;
}

template <typename QuadType, typename Quad1Type, typename ExprType>
requires std::integral<QuadType> && (!std::integral<Quad1Type>)
auto integrate_new(QuadType quad, Quad1Type const& quad1, ExprType const& expr)
{
    return /* quad is order, quad1 is quadrature object */;
}
#endif

//=============================================================================
// Example 7: Ranges and C++20 Views
//=============================================================================

#if __cplusplus >= 202002L
#include <ranges>

// Modern iteration with ranges
template <typename Mesh>
void process_boundary_with_ranges(std::shared_ptr<Mesh> mesh)
{
    // Old style
    for (auto it = mesh->beginFace(); it != mesh->endFace(); ++it)
    {
        if (it->isOnBoundary())
        {
            // process
        }
    }

    // Modern C++20 style
    for (auto& face : mesh->faces() | std::views::filter([](auto& f) { 
        return f.isOnBoundary(); 
    }))
    {
        // process - clearer intent
    }

    // Composable transformations
    auto boundary_areas = mesh->faces()
        | std::views::filter([](auto& f) { return f.isOnBoundary(); })
        | std::views::transform([](auto& f) { return f.measure(); });
}
#endif

//=============================================================================
// Comparison: Error Message Quality
//=============================================================================

/*
ERROR with SFINAE (C++17):
--------------------------
In file included from test.cpp:1:
In instantiation of 'typename std::enable_if_t<is_vf_expr_v<ExprT>, void> 
    process_old(const Container&) [with Container = int; 
    typename std::enable_if_t<is_vf_expr_v<ExprT>, void> = <type error>]':
test.cpp:42:23:   required from here
error: no type named 'type' in 'struct std::enable_if<false, void>'
note: candidate: 'template<class Container> typename std::enable_if_t<is_vf_expr_v<ExprT>, void> process_old(const Container&)'
note:   template argument deduction/substitution failed:
note:   couldn't deduce template parameter 'Container'
[... 200 more lines of template instantiation stack ...]

ERROR with Concepts (C++20):
----------------------------
error: cannot call process_new(int)
note: constraints not satisfied
note: the required expression 't.begin()' would be ill-formed
note: the required expression 't.end()' would be ill-formed
note: concept 'Iterable<int>' was not satisfied

MUCH CLEARER! ✅
*/

} // namespace Examples
} // namespace Feel

#endif /* FEELPP_MODERNIZATION_EXAMPLES_HPP */
