/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-05-31

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2026 Feel++ Consortium

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
   \file projectors.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-05-31

   Modernized to C++23 with concepts (2026-01-04)
 */
#ifndef FEELPP_VF_PROJECTORS_H
#define FEELPP_VF_PROJECTORS_H

#include <memory>
#include <type_traits>

#include <feel/feelcore/parameter.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/concepts.hpp>
#include <feel/feelvf/concepts.hpp>

#include <feel/feelvf/detail/clean.hpp>
#include <feel/feelvf/expr.hpp>

namespace Feel
{
namespace vf
{
namespace details
{

/**
 * @brief Base class for projectors using C++20/23 features
 *
 * @tparam iDim Projector type (NODAL, L2, etc.)
 * @tparam FunctionSpaceType The function space type (must satisfy FunctionSpaceConcept)
 * @tparam IteratorRange The range type for iteration (must satisfy RangeConcept)
 * @tparam ExprT The expression type to project
 *
 * @note Uses if constexpr for compile-time dispatch instead of tag dispatching.
 * @note Uses C++20 concepts for improved error messages and compile-time checking.
 */
template <ProjectorType iDim,
          FunctionSpaceConcept FunctionSpaceType,
          RangeConcept IteratorRange,
          typename ExprT>
class Projector
{
public:
    //
    // Type aliases (modern style)
    //
    static constexpr size_type context = ExprT::context | vm::POINT;

    using functionspace_type = FunctionSpaceType;
    using functionspace_ptrtype = std::shared_ptr<functionspace_type>;
    using element_type = typename functionspace_type::element_type;
    using basis_type = typename functionspace_type::basis_type;
    using expression_type = ExprT;
    using value_type = typename expression_type::value_type;

    using range_iterator = IteratorRange;
    using idim_type = typename IteratorRange::idim_t;
    using iterator_type = typename IteratorRange::iterator_t;

    //
    // Constructors
    //

    /**
     * @brief Construct a projector from a function space, range, and expression
     */
    constexpr Projector( functionspace_ptrtype const& functionspace,
                         IteratorRange const& range,
                         expression_type const& expr,
                         GeomapStrategyType geomap_strategy )
        : M_functionspace( functionspace )
        , M_range( range )
        , M_expr( expr )
        , M_geomap_strategy( geomap_strategy )
    {
        DVLOG( 2 ) << "Projector constructor from expression\n";
    }

    /**
     * @brief Copy constructor
     */
    constexpr Projector( Projector const& other ) = default;

    /**
     * @brief Move constructor
     */
    constexpr Projector( Projector&& other ) noexcept = default;

    /**
     * @brief Destructor
     */
    ~Projector() = default;

    //
    // Operators
    //

    /**
     * @brief Execute the projection
     * @param sum If true, accumulate contributions instead of replacing
     * @return The projected element
     */
    [[nodiscard]] element_type operator()( bool sum = false ) const
    {
        return applyProjection( sum );
    }

    //
    // Accessors
    //

    /**
     * @brief Get the variational expression
     */
    [[nodiscard]] constexpr expression_type const& expression() const noexcept
    {
        return M_expr;
    }

    /**
     * @brief Get the function space
     */
    [[nodiscard]] constexpr functionspace_ptrtype const& functionSpace() const noexcept
    {
        return M_functionspace;
    }

    /**
     * @brief Get the range
     */
    [[nodiscard]] constexpr range_iterator const& range() const noexcept
    {
        return M_range;
    }

    //
    // Compile-time properties
    //

    /**
     * @brief Polynomial order of the projection
     */
    [[nodiscard]] static constexpr uint16_type polynomialOrder() noexcept
    {
        return functionspace_type::basis_type::nOrder;
    }

    /**
     * @brief Expression is polynomial (always true for projectors)
     */
    [[nodiscard]] static constexpr bool isPolynomial() noexcept
    {
        return true;
    }

private:
    /**
     * @brief Apply projection using if constexpr for compile-time dispatch
     *
     * Replaces the old mpl::size_t<> tag dispatching pattern with modern
     * C++17/20 if constexpr, which is cleaner and generates better error messages.
     */
    [[nodiscard]] element_type applyProjection( bool sum ) const
    {
        constexpr auto dim_value = idim_type::value;

        if constexpr ( dim_value == MESH_ELEMENTS )
        {
            return projectOnElements( sum );
        }
        else if constexpr ( dim_value == MESH_FACES )
        {
            return projectOnFaces( sum );
        }
        else if constexpr ( dim_value == MESH_EDGES )
        {
            return projectOnEdges( sum );
        }
        else if constexpr ( dim_value == MESH_POINTS )
        {
            return projectOnPoints( sum );
        }
        else
        {
            static_assert( dim_value == MESH_ELEMENTS || dim_value == MESH_FACES ||
                               dim_value == MESH_EDGES || dim_value == MESH_POINTS,
                           "Unsupported mesh entity type for projection" );
        }
    }

    /**
     * @brief Project onto mesh elements
     */
    [[nodiscard]] element_type projectOnElements( bool sum ) const
    {
        element_type v( M_functionspace );
        FEELPP_ASSERT( v.size() == M_functionspace->dof()->nDof() )
            ( v.size() )( M_functionspace->dof()->nDof() )
                .warn( "invalid size" );
        v.setZero();
        v.on( _range = M_range, _expr = M_expr, _geomap = M_geomap_strategy, _accumulate = sum );
        return v;
    }

    /**
     * @brief Project onto mesh faces
     */
    [[nodiscard]] element_type projectOnFaces( bool sum ) const
    {
        element_type v( M_functionspace );
        v.setZero();
        v.on( _range = M_range, _expr = M_expr, _geomap = M_geomap_strategy, _accumulate = sum );
        return v;
    }

    /**
     * @brief Project onto mesh edges
     */
    [[nodiscard]] element_type projectOnEdges( bool sum ) const
    {
        element_type v( M_functionspace );
        v.setZero();
        v.on( _range = M_range, _expr = M_expr, _geomap = M_geomap_strategy, _accumulate = sum );
        return v;
    }

    /**
     * @brief Project onto mesh points
     */
    [[nodiscard]] element_type projectOnPoints( [[maybe_unused]] bool sum ) const
    {
        element_type v( M_functionspace );
        v.setZero();

        auto [ignore, pt_it, pt_en] = M_range;  // C++17 structured bindings

        if ( pt_it == pt_en )
            return v;

        // TODO: Implement point projection
        // For each DOF at the marked points, evaluate the expression
        // and set the corresponding DOF value

        return v;
    }

    //
    // Member data
    //
    functionspace_ptrtype const& M_functionspace;
    range_iterator M_range;
    expression_type const& M_expr;
    GeomapStrategyType M_geomap_strategy;
};

} // namespace details

//
// Free function interface with concepts
//

/**
 * @brief Nodal projection of an expression onto a function space subrange
 *
 * @tparam FunctionSpaceType Function space type (must satisfy FunctionSpaceConcept)
 * @tparam IteratorRange Range type (must satisfy RangeConcept)
 * @tparam ExprT Expression type
 *
 * @param functionspace The target function space
 * @param range_it The mesh range to project onto
 * @param expr The expression to project
 * @param geomap Geometric mapping strategy
 *
 * @return Element of the function space containing the projected expression
 */
template <FunctionSpaceConcept FunctionSpaceType,
          RangeConcept IteratorRange,
          typename ExprT>
[[nodiscard]] auto
project( std::shared_ptr<FunctionSpaceType> const& functionspace,
         IteratorRange const& range_it,
         Expr<ExprT> const& expr,
         GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
    -> typename FunctionSpaceType::element_type
{
    using projector_type = details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT>>;
    projector_type p( functionspace, range_it, expr, geomap );
    return p();
}

/**
 * @brief Implementation detail for project with named parameters
 */
template <FunctionSpaceConcept FunctionSpaceType,
          RangeConcept IteratorRange,
          typename ExprT>
[[nodiscard]] auto
project_impl( std::shared_ptr<FunctionSpaceType> const& functionspace,
              IteratorRange const& range_it,
              Expr<ExprT> const& expr,
              GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
    -> typename FunctionSpaceType::element_type
{
    using projector_type = details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT>>;
    projector_type p( functionspace, range_it, expr, geomap );
    return p();
}

/**
 * @brief Nodal projection of an expression onto an entire function space
 *
 * @tparam FunctionSpaceType Function space type (must satisfy FunctionSpaceConcept)
 * @tparam ExprT Expression type
 *
 * @param functionspace The target function space
 * @param expr The expression to project
 * @param geomap Geometric mapping strategy
 *
 * @return Element of the function space containing the projected expression
 */
template <FunctionSpaceConcept FunctionSpaceType, typename ExprT>
[[nodiscard]] auto
project( std::shared_ptr<FunctionSpaceType> const& functionspace,
         Expr<ExprT> const& expr,
         GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
    -> typename FunctionSpaceType::element_type
{
    return project( functionspace, elements( functionspace->mesh() ), expr, geomap );
}

/**
 * @brief Sum projection (accumulate nodal contributions)
 *
 * @tparam FunctionSpaceType Function space type (must satisfy FunctionSpaceConcept)
 * @tparam IteratorRange Range type (must satisfy RangeConcept)
 * @tparam ExprT Expression type
 *
 * @param functionspace The target function space
 * @param range_it The mesh range
 * @param expr The expression to project
 * @param geomap Geometric mapping strategy
 * @param parallelSync Whether to synchronize in parallel
 *
 * @return Element with accumulated contributions
 */
template <FunctionSpaceConcept FunctionSpaceType,
          RangeConcept IteratorRange,
          typename ExprT>
[[nodiscard]] auto
sum( std::shared_ptr<FunctionSpaceType> const& functionspace,
     IteratorRange const& range_it,
     Expr<ExprT> const& expr,
     GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_OPT,
     bool parallelSync = true )
    -> typename FunctionSpaceType::element_type
{
    using projector_type = details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT>>;
    projector_type p( functionspace, range_it, expr, geomap );
    auto res = p( true );
    if ( parallelSync )
        sync( res, "+" );
    return res;
}

/**
 * @brief Sum projection onto entire function space
 *
 * @tparam FunctionSpaceType Function space type (must satisfy FunctionSpaceConcept)
 * @tparam ExprT Expression type
 */
template <FunctionSpaceConcept FunctionSpaceType, typename ExprT>
[[nodiscard]] auto
sum( std::shared_ptr<FunctionSpaceType> const& functionspace,
     Expr<ExprT> const& expr,
     GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_OPT,
     bool parallelSync = true )
    -> typename FunctionSpaceType::element_type
{
    return sum( functionspace, elements( functionspace->mesh() ), expr, geomap, parallelSync );
}

/// \cond DETAIL
namespace detail
{

template <typename S>
struct space_ptr
{
    using type = typename S::element_type;
};

template <typename S>
struct space_value
{
    using type = S;
};

} // namespace detail
/// \endcond

/**
 * @brief Named parameter interface for projection
 *
 * @param _space The function space to project onto
 * @param _range The range of mesh elements (default: all elements)
 * @param _expr The expression to project
 * @param _geomap Geometric mapping strategy
 * @param _accumulate Whether to accumulate contributions
 *
 * @example
 * @code
 * auto u = project(_space=Vh, _expr=sin(Px()*Py()));
 * auto u = project(_space=Vh, _range=markedfaces(mesh, "inlet"), _expr=cst(1.0));
 * @endcode
 */
template <typename... Ts>
[[nodiscard]] auto project( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& space = args.get( _space );
    auto&& expr = args.get( _expr );
    auto&& range = args.get_else_invocable( _range, [&space]() { return elements( support( space ) ); } );
    GeomapStrategyType geomap = args.get_else( _geomap, GeomapStrategyType::GEOMAP_OPT );
    [[maybe_unused]] bool accumulate = args.get_else( _accumulate, false );

    return project_impl( space, range, expr, Feel::detail::geomapStrategy( range, geomap ) );
}

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_PROJECTORS_H */
