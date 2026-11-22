/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELPDE_BOUNDARYCONDITIONDESCRIPTORS_HPP
#define FEELPP_FEELPDE_BOUNDARYCONDITIONDESCRIPTORS_HPP 1

#include <concepts>
#include <cstddef>
#include <functional>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <feel/feelcore/parameter.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel::pde
{

/**
 * @brief Boundary condition descriptor bound to a mesh.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 */
template <typename MeshType>
class BoundaryConditionDescriptor
{
public:
    using mesh_type = std::decay_t<MeshType>;
    using id_type = ModelBoundaryId;
    using condition_type = ModelBoundaryCondition;

    /** @brief Build a descriptor from a model id, condition, and mesh. */
    BoundaryConditionDescriptor( id_type id, condition_type condition, mesh_type mesh )
        : M_id( std::move( id ) ),
          M_condition( std::move( condition ) ),
          M_mesh( std::move( mesh ) )
    {}

    /** @brief Return the boundary condition id. */
    id_type const& id() const noexcept { return M_id; }

    /** @brief Return the field name targeted by this condition. */
    std::string const& field() const noexcept { return M_id.field(); }

    /** @brief Return the boundary condition type. */
    std::string const& type() const noexcept { return M_id.type(); }

    /** @brief Return the boundary condition type. */
    std::string const& kind() const noexcept { return this->type(); }

    /** @brief Return the user-defined boundary condition name. */
    std::string const& name() const noexcept { return M_id.name(); }

    /** @brief Return the underlying model boundary condition. */
    condition_type const& condition() const noexcept { return M_condition; }

    /** @brief Return the first marker name when present. */
    std::string const& marker() const
    {
        static std::string const empty;
        auto const& markerSet = this->markers();
        return markerSet.empty() ? empty : *markerSet.begin();
    }

    /** @brief Return the marker set associated with the condition. */
    std::set<std::string> const& markers() const noexcept { return M_condition.markers(); }

    /** @brief Return the material name associated with the condition. */
    std::string const& material() const noexcept { return M_condition.material(); }

    /** @brief Return true when the primary expression is present. */
    bool hasExpression() const { return M_condition.hasExpression(); }

    /** @brief Return true when expression 1 is present. */
    bool hasExpression1() const { return M_condition.hasExpression1(); }

    /** @brief Return true when expression 2 is present. */
    bool hasExpression2() const { return M_condition.hasExpression2(); }

    /** @brief Return the primary model expression. */
    template <int M = 1, int N = 1>
    auto expr() const { return M_condition.template expr<M, N>(); }

    /** @brief Return expression 1. */
    template <int M = 1, int N = 1>
    auto expr1() const { return M_condition.template expr1<M, N>(); }

    /** @brief Return expression 2. */
    template <int M = 1, int N = 1>
    auto expr2() const { return M_condition.template expr2<M, N>(); }

    /** @brief Return the mesh handle bound to this descriptor. */
    mesh_type const& mesh() const noexcept { return M_mesh; }

    /** @brief Return marked faces matching the descriptor markers. */
    auto faces() const { return markedfaces( M_mesh, this->markers() ); }

    /** @brief Return marked elements matching the descriptor markers. */
    auto elements() const { return markedelements( M_mesh, this->markers() ); }

private:
    id_type M_id;
    condition_type M_condition;
    mesh_type M_mesh;
};

/**
 * @brief Range-like view over boundary condition descriptors.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 */
template <typename MeshType>
class BoundaryConditionView
{
public:
    using mesh_type = std::decay_t<MeshType>;
    using descriptor_type = BoundaryConditionDescriptor<mesh_type>;
    using container_type = std::vector<descriptor_type>;
    using const_iterator = typename container_type::const_iterator;

    /** @brief Build an empty boundary condition view. */
    BoundaryConditionView() = default;

    /** @brief Build a boundary condition view from descriptors. */
    explicit BoundaryConditionView( container_type items )
        : M_items( std::move( items ) )
    {}

    /** @brief Return an iterator to the first descriptor. */
    const_iterator begin() const noexcept { return M_items.begin(); }

    /** @brief Return an iterator past the last descriptor. */
    const_iterator end() const noexcept { return M_items.end(); }

    /** @brief Return an iterator to the first descriptor. */
    const_iterator cbegin() const noexcept { return M_items.cbegin(); }

    /** @brief Return an iterator past the last descriptor. */
    const_iterator cend() const noexcept { return M_items.cend(); }

    /** @brief Return the number of descriptors in the view. */
    std::size_t size() const noexcept { return M_items.size(); }

    /** @brief Return true when the view contains no descriptor. */
    bool empty() const noexcept { return M_items.empty(); }

    /** @brief Access a descriptor by index. */
    descriptor_type const& operator[]( std::size_t index ) const { return M_items[index]; }

    /** @brief Return a filtered view matching a predicate. */
    template <typename Predicate>
        requires std::predicate<Predicate&, descriptor_type const&>
    BoundaryConditionView where( Predicate predicate ) const
    {
        container_type items;
        for ( auto const& item : M_items )
            if ( std::invoke( predicate, item ) )
                items.push_back( item );
        return BoundaryConditionView{ std::move( items ) };
    }

    /** @brief Return a view containing only conditions for one field. */
    BoundaryConditionView field( std::string const& name ) const
    {
        return this->where( [&name]( descriptor_type const& item ) { return item.field() == name; } );
    }

    /** @brief Return a view containing only conditions of one type. */
    BoundaryConditionView type( std::string const& name ) const
    {
        return this->where( [&name]( descriptor_type const& item ) { return item.type() == name; } );
    }

    /** @brief Return a view containing only conditions of one type. */
    BoundaryConditionView kind( std::string const& name ) const { return this->type( name ); }

    /** @brief Return Dirichlet boundary conditions. */
    BoundaryConditionView dirichlet() const { return this->type( "Dirichlet" ); }

    /** @brief Return Neumann boundary conditions. */
    BoundaryConditionView neumann() const { return this->type( "Neumann" ); }

    /** @brief Return Robin boundary conditions. */
    BoundaryConditionView robin() const { return this->type( "Robin" ); }

    /** @brief Return Fourier boundary conditions. */
    BoundaryConditionView fourier() const { return this->type( "Fourier" ); }

    /** @brief Return volumic load conditions. */
    BoundaryConditionView load() const { return this->type( "Load" ); }

private:
    container_type M_items;
};

/**
 * @brief Build a boundary condition view from model properties.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 * @param model Model properties containing BoundaryConditions2 data.
 * @param field Field name to keep; an empty string keeps all fields.
 * @param mesh Mesh handle bound to each descriptor.
 * @return Boundary condition descriptor view.
 */
template <typename MeshType>
auto
boundaryConditions( ModelProperties const& model, std::string const& field, MeshType mesh )
{
    using mesh_type = std::decay_t<MeshType>;
    using view_type = BoundaryConditionView<mesh_type>;
    using descriptor_type = typename view_type::descriptor_type;

    std::vector<descriptor_type> items;
    for ( auto const& [id, condition] : model.boundaryConditions2().flatten() )
    {
        if ( !field.empty() && id.field() != field )
            continue;
        items.emplace_back( id, condition, mesh );
    }
    return view_type{ std::move( items ) };
}

/**
 * @brief Build a boundary condition view from model properties.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 * @param model Model properties containing BoundaryConditions2 data.
 * @param mesh Mesh handle bound to each descriptor.
 * @return Boundary condition descriptor view.
 */
template <typename MeshType>
auto
boundaryConditions( ModelProperties const& model, MeshType mesh )
{
    return boundaryConditions( model, std::string{}, std::move( mesh ) );
}

/**
 * @brief Build a boundary condition view from named arguments.
 *
 * Supported arguments are `_model`, `_mesh`, and optional `_name` for the
 * selected field.
 *
 * @tparam Args Named argument types.
 * @param args Named arguments.
 * @return Boundary condition descriptor view.
 */
template <typename... Args>
    requires ( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
auto
boundaryConditions( Args&&... args )
{
    auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
    auto const& model = namedArgs.get( _model );
    auto mesh = namedArgs.get( _mesh );
    std::string field = namedArgs.get_else( _name, std::string{} );
    return boundaryConditions( model, field, std::move( mesh ) );
}

} // namespace Feel::pde

#endif
