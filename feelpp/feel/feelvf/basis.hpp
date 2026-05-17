/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-25

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
   \file basis.hpp
   \author Christophe Prud'homme
   \date 2026-03-25
 */
#ifndef FEELPP_VF_BASIS_HPP
#define FEELPP_VF_BASIS_HPP 1

#include <concepts>
#include <string>
#include <type_traits>
#include <utility>

#include <feel/feelvf/operators.hpp>

namespace Feel
{
namespace vf
{

namespace detail
{

template <OperatorType Role>
concept BasisRole = ( Role == __TEST ) || ( Role == __TRIAL );

template <typename T>
concept FunctionSpaceElement =
    std::is_base_of_v<FunctionSpaceBase::ElementBase, std::remove_cvref_t<T>>;

template <typename T>
concept FunctionSpaceHandle = requires( T const& t, std::string const& name, std::string const& desc )
{
    typename std::remove_reference_t<decltype( *t )>::element_type;
    { t->element( name, desc ) } -> std::same_as<typename std::remove_reference_t<decltype( *t )>::element_type>;
};

template <FunctionSpaceElement ElementType>
class BasisStorage
{
protected:
    explicit BasisStorage( ElementType element )
        :
        M_element( std::move( element ) )
    {}

    ElementType const& storedElement() const noexcept
    {
        return M_element;
    }

    ElementType& storedElement() noexcept
    {
        return M_element;
    }

private:
    ElementType M_element;
};

template <FunctionSpaceElement ElementType, OperatorType Role>
    requires BasisRole<Role>
class BasisProxy
    : private BasisStorage<ElementType>,
      public Expr<OpId<ElementType, Role>>
{
    using storage_type = BasisStorage<ElementType>;
    using terminal_type = OpId<ElementType, Role>;
    using expr_base_type = Expr<terminal_type>;

public:
    using element_type = ElementType;
    static constexpr OperatorType role = Role;

    explicit BasisProxy( element_type element )
        :
        storage_type( std::move( element ) ),
        expr_base_type( terminal_type( this->element() ) )
    {}

    // Rebind the terminal reference on copies because OpId stores a reference_wrapper.
    BasisProxy( BasisProxy const& other )
        :
        storage_type( other.element() ),
        expr_base_type( terminal_type( this->element() ) )
    {}

    BasisProxy( BasisProxy&& other ) noexcept( std::is_nothrow_move_constructible_v<element_type> )
        :
        storage_type( std::move( other.storage_type::storedElement() ) ),
        expr_base_type( terminal_type( this->element() ) )
    {}

    BasisProxy& operator=( BasisProxy const& other )
    {
        if ( this != &other )
        {
            this->storage_type::storedElement() = other.element();
            this->rebindExpr();
        }
        return *this;
    }

    BasisProxy& operator=( BasisProxy&& other ) noexcept( std::is_nothrow_move_assignable_v<element_type> )
    {
        if ( this != &other )
        {
            this->storage_type::storedElement() = std::move( other.storage_type::storedElement() );
            this->rebindExpr();
        }
        return *this;
    }

    [[nodiscard]] element_type const& element() const noexcept
    {
        return this->storage_type::storedElement();
    }

private:
    void rebindExpr()
    {
        expr_base_type::operator=( expr_base_type( terminal_type( this->element() ) ) );
    }
};

template <typename T>
concept BasisProxyType = requires( T const& t )
{
    typename std::remove_cvref_t<T>::element_type;
    { std::remove_cvref_t<T>::role } -> std::convertible_to<OperatorType>;
    { t.element() } -> std::same_as<typename std::remove_cvref_t<T>::element_type const&>;
};

template <BasisProxyType ProxyType>
using basis_proxy_type_t = std::remove_cvref_t<ProxyType>;

} // namespace detail

template <detail::FunctionSpaceElement ElementType>
[[nodiscard]] inline auto
trial( ElementType element )
{
    return detail::BasisProxy<ElementType, __TRIAL>( std::move( element ) );
}

template <detail::FunctionSpaceHandle SpaceHandle>
[[nodiscard]] inline auto
trial( SpaceHandle const& space,
       std::string name = "u",
       std::string desc = "u" )
{
    using space_type = std::remove_reference_t<decltype( *space )>;
    using element_type = typename space_type::element_type;
    return detail::BasisProxy<element_type, __TRIAL>( space->element( std::move( name ), std::move( desc ) ) );
}

template <detail::FunctionSpaceElement ElementType>
[[nodiscard]] inline auto
test( ElementType element )
{
    return detail::BasisProxy<ElementType, __TEST>( std::move( element ) );
}

template <detail::FunctionSpaceHandle SpaceHandle>
[[nodiscard]] inline auto
test( SpaceHandle const& space,
      std::string name = "v",
      std::string desc = "v" )
{
    using space_type = std::remove_reference_t<decltype( *space )>;
    using element_type = typename space_type::element_type;
    return detail::BasisProxy<element_type, __TEST>( space->element( std::move( name ), std::move( desc ) ) );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
id( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return idt( proxy.element() );
    else
        return id( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
normal( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return normalt( proxy.element() );
    else
        return normal( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
dx( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return dxt( proxy.element() );
    else
        return dx( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
dy( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return dyt( proxy.element() );
    else
        return dy( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
dz( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return dzt( proxy.element() );
    else
        return dz( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
dn( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return dnt( proxy.element() );
    else
        return dn( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
grad( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return gradt( proxy.element() );
    else
        return grad( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
symm_grad( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return symm_gradt( proxy.element() );
    else
        return symm_grad( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
div( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return divt( proxy.element() );
    else
        return div( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
curl( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return curlt( proxy.element() );
    else
        return curl( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
curlx( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return curlxt( proxy.element() );
    else
        return curlx( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
curly( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return curlyt( proxy.element() );
    else
        return curly( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
curlz( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return curlzt( proxy.element() );
    else
        return curlz( proxy.element() );
}

template <detail::FunctionSpaceElement ElementType>
[[nodiscard]] inline auto
omega( ElementType const& element )
{
    return 0.5*curl( element );
}

template <detail::FunctionSpaceElement ElementType>
[[nodiscard]] inline auto
omegat( ElementType const& element )
{
    return 0.5*curlt( element );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
omega( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return omegat( proxy.element() );
    else
        return omega( proxy.element() );
}

template <detail::FunctionSpaceElement ElementType>
[[nodiscard]] inline auto
omegaz( ElementType const& element )
{
    return 0.5*curlz( element );
}

template <detail::FunctionSpaceElement ElementType>
[[nodiscard]] inline auto
omegazt( ElementType const& element )
{
    return 0.5*curlzt( element );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
omegaz( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return omegazt( proxy.element() );
    else
        return omegaz( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
hess( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return hesst( proxy.element() );
    else
        return hess( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
laplacian( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return laplaciant( proxy.element() );
    else
        return laplacian( proxy.element() );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
trace( ProxyType const& proxy )
{
    if constexpr ( detail::basis_proxy_type_t<ProxyType>::role == __TRIAL )
        return tracet( proxy.element() );
    else
        return trace( proxy.element() );
}

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_BASIS_HPP */
