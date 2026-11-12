/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELPDE_MATERIALDESCRIPTORS_HPP
#define FEELPP_FEELPDE_MATERIALDESCRIPTORS_HPP 1

#include <concepts>
#include <cstddef>
#include <functional>
#include <initializer_list>
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
 * @brief Material descriptor bound to a mesh.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 */
template <typename MeshType>
class MaterialDescriptor
{
public:
    using mesh_type = std::decay_t<MeshType>;
    using material_type = ModelMaterial;
    using property_type = ModelMaterialProperty;

    /** @brief Build a material descriptor. */
    MaterialDescriptor( std::string key, material_type material, mesh_type mesh )
        : M_key( std::move( key ) ),
          M_material( std::move( material ) ),
          M_mesh( std::move( mesh ) )
    {}

    /** @brief Return the key used in the material map. */
    std::string const& key() const noexcept { return M_key; }

    /** @brief Return the material name. */
    std::string const& name() const noexcept { return M_material.name(); }

    /** @brief Return the underlying model material. */
    material_type const& material() const noexcept { return M_material; }

    /** @brief Return the first mesh marker when present. */
    std::string const& marker() const
    {
        static std::string const empty;
        auto const& markerSet = this->markers();
        return markerSet.empty() ? empty : *markerSet.begin();
    }

    /** @brief Return mesh markers associated with the material. */
    std::set<std::string> const& markers() const noexcept { return M_material.meshMarkers(); }

    /** @brief Return physics attached to the material. */
    std::set<std::string> const& physics() const noexcept { return M_material.physics(); }

    /** @brief Return the first physic attached to the material. */
    std::string physic() const { return M_material.physic(); }

    /** @brief Return true when the material has any physic. */
    bool hasPhysics() const { return M_material.hasPhysics(); }

    /** @brief Return true when the material has a given physic. */
    bool hasPhysics( std::string const& physic ) const { return M_material.hasPhysics( physic ); }

    /** @brief Return true when the material has one of the requested physics. */
    bool hasPhysics( std::set<std::string> const& physics ) const { return M_material.hasPhysics( physics ); }

    /** @brief Return true when a material property exists. */
    bool hasProperty( std::string const& name ) const { return M_material.hasProperty( name ); }

    /** @brief Return a material property by name. */
    property_type const& property( std::string const& name ) const { return M_material.property( name ); }

    /** @brief Return the mesh handle bound to this descriptor. */
    mesh_type const& mesh() const noexcept { return M_mesh; }

    /** @brief Return marked elements matching the material markers. */
    auto elements() const { return markedelements( M_mesh, this->markers() ); }

    /** @brief Return marked faces matching the material markers. */
    auto faces() const { return markedfaces( M_mesh, this->markers() ); }

private:
    std::string M_key;
    material_type M_material;
    mesh_type M_mesh;
};

/**
 * @brief Range-like view over material descriptors.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 */
template <typename MeshType>
class MaterialView
{
public:
    using mesh_type = std::decay_t<MeshType>;
    using descriptor_type = MaterialDescriptor<mesh_type>;
    using container_type = std::vector<descriptor_type>;
    using const_iterator = typename container_type::const_iterator;

    /** @brief Build an empty material view. */
    MaterialView() = default;

    /** @brief Build a material view from descriptors. */
    explicit MaterialView( container_type items )
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
    MaterialView where( Predicate predicate ) const
    {
        container_type items;
        for ( auto const& item : M_items )
            if ( std::invoke( predicate, item ) )
                items.push_back( item );
        return MaterialView{ std::move( items ) };
    }

    /** @brief Return materials that define a property. */
    MaterialView hasProperty( std::string const& name ) const
    {
        return this->where( [&name]( descriptor_type const& item ) { return item.hasProperty( name ); } );
    }

    /** @brief Return materials with a given physic. */
    MaterialView withPhysic( std::string const& physic ) const
    {
        return this->where( [&physic]( descriptor_type const& item ) { return item.hasPhysics( physic ); } );
    }

    /** @brief Return materials with one of the requested physics. */
    MaterialView withPhysics( std::set<std::string> const& physics ) const
    {
        return this->where( [&physics]( descriptor_type const& item ) { return item.hasPhysics( physics ); } );
    }

    /** @brief Return materials with one of the requested physics. */
    MaterialView withPhysics( std::initializer_list<std::string> physics ) const
    {
        return this->withPhysics( std::set<std::string>( physics ) );
    }

private:
    container_type M_items;
};

/**
 * @brief Build a material view from model properties.
 *
 * @tparam MeshType Mesh handle type used by Feel++ filters.
 * @param model Model properties containing material data.
 * @param mesh Mesh handle bound to each descriptor.
 * @return Material descriptor view.
 */
template <typename MeshType>
auto
materials( ModelProperties const& model, MeshType mesh )
{
    using mesh_type = std::decay_t<MeshType>;
    using view_type = MaterialView<mesh_type>;
    using descriptor_type = typename view_type::descriptor_type;

    std::vector<descriptor_type> items;
    for ( auto const& [key, material] : model.materials() )
        items.emplace_back( key, material, mesh );
    return view_type{ std::move( items ) };
}

/**
 * @brief Build a material view from named arguments.
 *
 * Supported arguments are `_model` and `_mesh`.
 *
 * @tparam Args Named argument types.
 * @param args Named arguments.
 * @return Material descriptor view.
 */
template <typename... Args>
    requires ( sizeof...( Args ) > 0 && ( NA::is_named_argument_v<std::remove_cvref_t<Args>> && ... ) )
auto
materials( Args&&... args )
{
    auto namedArgs = NA::make_arguments( std::forward<Args>( args )... );
    auto const& model = namedArgs.get( _model );
    auto mesh = namedArgs.get( _mesh );
    return materials( model, std::move( mesh ) );
}

} // namespace Feel::pde

#endif
