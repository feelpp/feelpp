/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_HDGSPACES_HPP
#define FEELPP_FEELSYSTEM_HDGSPACES_HPP 1

#include <concepts>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <boost/hana/integral_constant.hpp>

#include <feel/feelalg/productspaceconcepts.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelsystem/roles.hpp>
#include <feel/feelvf/basis.hpp>

namespace Feel::systems
{

namespace detail
{
template <typename Role, typename... RoleTypes>
struct role_index_impl;

template <typename Role>
struct role_index_impl<Role>
{
    inline static constexpr int value = -1;
};

template <typename Role, typename FirstRole, typename... RestRoles>
struct role_index_impl<Role, FirstRole, RestRoles...>
{
    inline static constexpr int next = role_index_impl<Role, RestRoles...>::value;
    inline static constexpr int value =
        std::same_as<std::remove_cvref_t<Role>, std::remove_cvref_t<FirstRole>> ? 0 :
        ( next < 0 ? -1 : 1 + next );
};

template <typename Role, typename... RoleTypes>
inline constexpr int role_index_v = role_index_impl<std::remove_cvref_t<Role>, RoleTypes...>::value;

template <int Index, typename ProductSpaceType>
using product_space_component_t =
    std::remove_cvref_t<decltype( std::declval<ProductSpaceType const&>()[boost::hana::int_c<Index>] )>;
} // namespace detail

/**
 * @brief View of one named role space in an HDG/product-space layout.
 *
 * The view is intentionally lightweight. It owns only the function-space
 * pointer copied from the product space and delegates basis construction to the
 * first-class `Feel::vf::trial(Xh)` and `Feel::vf::test(Xh)` helpers.
 */
template <BlockRole Role, int Index, typename SpacePtrType>
class HdgRoleSpace
{
public:
    using role_type = std::remove_cvref_t<Role>;
    using space_ptrtype = std::remove_cvref_t<SpacePtrType>;
    using space_type = typename space_ptrtype::element_type;

    inline static constexpr int index_value = Index;

    /**
     * @brief Construct a role view from its function space.
     */
    explicit HdgRoleSpace( space_ptrtype space )
        : M_space( std::move( space ) )
    {
        if ( !M_space )
            throw std::invalid_argument( "HdgRoleSpace requires a non-null function space" );
    }

    /**
     * @brief Function space associated with this role.
     */
    space_ptrtype const& space() const noexcept { return M_space; }

    /**
     * @brief Compile-time block index for this role.
     */
    constexpr auto index() const noexcept { return boost::hana::int_c<Index>; }

    /**
     * @brief Human-readable role name.
     */
    constexpr std::string_view name() const noexcept { return blockRoleName( role_type{} ); }

    /**
     * @brief Create a first-class trial basis on this role space.
     */
    [[nodiscard]] auto trial() const { return Feel::vf::trial( M_space ); }

    /**
     * @brief Create a first-class test basis on this role space.
     */
    [[nodiscard]] auto test() const { return Feel::vf::test( M_space ); }

    /**
     * @brief Create an unnamed finite-element field in this role space.
     */
    [[nodiscard]] auto element() const { return M_space->element(); }

    /**
     * @brief Create a named finite-element field in this role space.
     */
    [[nodiscard]] auto element( std::string const& name ) const { return M_space->element( name ); }

private:
    /** @brief Function space associated with the semantic role. */
    space_ptrtype M_space;
};

/**
 * @brief Static HDG product-space descriptor with semantic block roles.
 *
 * `HdgSpaces` stores the volume mesh, skeleton/facet mesh, product space, and a
 * compile-time role list. It does not assemble HDG equations; it only makes the
 * existing product-space layout explicit and testable.
 */
template <typename MeshPtrType, typename FacetMeshPtrType,
          StaticProductSpacesType ProductSpaceType, BlockRole... RoleTypes>
class HdgSpaces
{
public:
    using mesh_ptrtype = std::remove_cvref_t<MeshPtrType>;
    using facet_mesh_ptrtype = std::remove_cvref_t<FacetMeshPtrType>;
    using product_space_type = std::remove_cvref_t<ProductSpaceType>;
    using element_type = typename product_space_type::element_type;
    using value_type = typename product_space_type::value_type;

    inline static constexpr int nroles = sizeof...( RoleTypes );

    static_assert( nroles > 0, "HdgSpaces requires at least one semantic role" );
    static_assert( product_space_type::element_type::nspaces == nroles,
                   "HdgSpaces role count must match the static product-space block count" );

    /**
     * @brief Construct HDG spaces from the volume mesh, facet mesh, and product.
     */
    HdgSpaces( mesh_ptrtype mesh, facet_mesh_ptrtype facetMesh, product_space_type product )
        : M_mesh( std::move( mesh ) ),
          M_facetMesh( std::move( facetMesh ) ),
          M_product( std::move( product ) )
    {
        if ( !M_mesh )
            throw std::invalid_argument( "HdgSpaces requires a non-null mesh" );
        if ( !M_facetMesh )
            throw std::invalid_argument( "HdgSpaces requires a non-null facet mesh" );
    }

    /**
     * @brief Volume mesh.
     */
    mesh_ptrtype const& mesh() const noexcept { return M_mesh; }

    /**
     * @brief Facet mesh used by trace spaces.
     */
    facet_mesh_ptrtype const& facetMesh() const noexcept { return M_facetMesh; }

    /**
     * @brief Product space consumed by block forms.
     */
    product_space_type const& product() const noexcept { return M_product; }

    /**
     * @brief Mutable product space consumed by block forms.
     */
    product_space_type& product() noexcept { return M_product; }

    /**
     * @brief Alias for product-space consumers.
     */
    product_space_type const& space() const noexcept { return M_product; }

    /**
     * @brief Mutable product-space alias.
     */
    product_space_type& space() noexcept { return M_product; }

    /**
     * @brief Number of semantic roles/product blocks.
     */
    constexpr int numberOfRoles() const noexcept { return nroles; }

    /**
     * @brief Number of product-space blocks.
     */
    int numberOfSpaces() const noexcept { return M_product.numberOfSpaces(); }

    /**
     * @brief Create a product-space element.
     */
    [[nodiscard]] element_type element() { return M_product.element(); }

    /**
     * @brief True when this layout contains a role.
     */
    template <typename Role>
    inline static constexpr bool has_role = detail::role_index_v<Role, RoleTypes...> >= 0;

    /**
     * @brief True when this layout contains a role.
     */
    template <BlockRole Role>
    constexpr bool hasRole( Role = {} ) const noexcept
    {
        return has_role<std::remove_cvref_t<Role>>;
    }

    /**
     * @brief Numeric block index associated with a role.
     */
    template <BlockRole Role>
        requires( has_role<std::remove_cvref_t<Role>> )
    static consteval int indexValue( Role = {} ) noexcept
    {
        return detail::role_index_v<Role, RoleTypes...>;
    }

    /**
     * @brief Hana block index associated with a role.
     */
    template <BlockRole Role>
        requires( has_role<std::remove_cvref_t<Role>> )
    constexpr auto index( Role = {} ) const noexcept
    {
        return boost::hana::int_c<indexValue( Role{} )>;
    }

    /**
     * @brief Role-space view for a semantic role.
     */
    template <BlockRole Role>
        requires( has_role<std::remove_cvref_t<Role>> )
    [[nodiscard]] auto component( Role = {} ) const
    {
        constexpr int roleIndex = indexValue( Role{} );
        using space_ptrtype = detail::product_space_component_t<roleIndex, product_space_type>;
        return HdgRoleSpace<std::remove_cvref_t<Role>, roleIndex, space_ptrtype>(
            M_product[boost::hana::int_c<roleIndex>] );
    }

    /**
     * @brief Flux role-space view.
     */
    [[nodiscard]] auto flux() const
        requires( has_role<roles::flux_t> )
    {
        return this->component( roles::flux );
    }

    /**
     * @brief Potential role-space view.
     */
    [[nodiscard]] auto potential() const
        requires( has_role<roles::potential_t> )
    {
        return this->component( roles::potential );
    }

    /**
     * @brief Trace role-space view.
     */
    [[nodiscard]] auto trace() const
        requires( has_role<roles::trace_t> )
    {
        return this->component( roles::trace );
    }

    /**
     * @brief Stress role-space view.
     */
    [[nodiscard]] auto stress() const
        requires( has_role<roles::stress_t> )
    {
        return this->component( roles::stress );
    }

    /**
     * @brief Velocity role-space view.
     */
    [[nodiscard]] auto velocity() const
        requires( has_role<roles::velocity_t> )
    {
        return this->component( roles::velocity );
    }

    /**
     * @brief Pressure role-space view.
     */
    [[nodiscard]] auto pressure() const
        requires( has_role<roles::pressure_t> )
    {
        return this->component( roles::pressure );
    }

    /**
     * @brief Displacement role-space view.
     */
    [[nodiscard]] auto displacement() const
        requires( has_role<roles::displacement_t> )
    {
        return this->component( roles::displacement );
    }

    /**
     * @brief Constraint role-space view.
     */
    [[nodiscard]] auto constraint() const
        requires( has_role<roles::constraint_t> )
    {
        return this->component( roles::constraint );
    }

private:
    /** @brief Volume mesh used by volume spaces. */
    mesh_ptrtype M_mesh;
    /** @brief Facet mesh used by trace/skeleton spaces. */
    facet_mesh_ptrtype M_facetMesh;
    /** @brief Product space with components ordered like the role list. */
    product_space_type M_product;
};

/**
 * @brief Create an `HdgSpaces` descriptor from a product and explicit roles.
 */
template <BlockRole... RoleTypes, typename MeshPtrType, typename FacetMeshPtrType,
          StaticProductSpacesType ProductSpaceType>
[[nodiscard]] auto
makeHdgSpaces( MeshPtrType mesh, FacetMeshPtrType facetMesh, ProductSpaceType product )
{
    return HdgSpaces<MeshPtrType, FacetMeshPtrType, ProductSpaceType, RoleTypes...>(
        std::move( mesh ), std::move( facetMesh ), std::move( product ) );
}

} // namespace Feel::systems

namespace Feel::hdg
{

/**
 * @brief Construct the current three-block Poisson/Laplacian HDG space layout.
 *
 * The returned layout is semantically equivalent to the existing quickstart
 * sequence `Vh = Pdhv<OrderP>(mesh)`, `Wh = Pdh<OrderP>(mesh)`,
 * `Mh = Pdh<OrderP>(createSubmesh(faces(mesh)))`, and
 * `product(Vh, Wh, Mh)`.
 */
template <int OrderP, typename MeshType>
[[nodiscard]] auto
poissonSpaces( std::shared_ptr<MeshType> const& mesh )
{
    auto Vh = Pdhv<OrderP>( mesh );
    auto Wh = Pdh<OrderP>( mesh );
    auto facetMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    auto Mh = Pdh<OrderP>( facetMesh );
    auto ps = product( Vh, Wh, Mh );
    return systems::makeHdgSpaces<systems::roles::flux_t,
                                  systems::roles::potential_t,
                                  systems::roles::trace_t>(
        mesh, facetMesh, std::move( ps ) );
}

} // namespace Feel::hdg

#endif // FEELPP_FEELSYSTEM_HDGSPACES_HPP
