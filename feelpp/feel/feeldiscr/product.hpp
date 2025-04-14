/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-23

  Copyright (C) 2014-2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file product.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-03-23
 */
#if !defined( FEELPP_FEELDISCR_PRODUCT_HPP )
#define FEELPP_FEELDISCR_PRODUCT_HPP 1

#include <boost/hana/ext/std/vector.hpp>
#include <feel/feelalg/vectorblock.hpp>

#include <boost/fusion/view/joint_view.hpp>

namespace Feel
{

constexpr auto is_void = hana::integral( hana::metafunction<std::is_void> );

struct product_space_tag
{
};

template <typename T, bool same_mesh = false>
class ProductSpace : public std::vector<std::shared_ptr<decay_type<T>>>, ProductSpaceBase
{
  public:
    using feel_tag = product_space_tag;

    using super = std::vector<std::shared_ptr<decay_type<T>>>;
    // using value_type = typename decay_type<decltype(super[0_c])>::value_type;
    using functionspace_type = ProductSpace<T, same_mesh>;
    using underlying_functionspace_type = decay_type<T>;
    using underlying_functionspace_ptrtype = std::shared_ptr<underlying_functionspace_type>;
    using mesh_type = typename underlying_functionspace_type::mesh_type;
    using mesh_ptrtype = typename underlying_functionspace_type::mesh_ptrtype;
    using value_type = typename underlying_functionspace_type::value_type;
    using worldcomm_type = WorldComm;

    /**
     * construct a product of n identical spaces from mesh \p m
     */
    ProductSpace( int n, mesh_ptrtype m )
        : super( 1, underlying_functionspace_type::New( m ) ),
          M_nspaces( n ),
          M_props( M_nspaces, "Ibc" )
    {
    }

    /**
     * construct a product of n spaces from a vector of  n different meshes
     */
    ProductSpace( std::vector<mesh_ptrtype> const& m )
        : super( m.size() ),
          M_nspaces( m.size() ),
          M_props( M_nspaces, "Ibc" )
    {
        std::transform( this->begin(), this->end(), m.begin(),
                        []( auto const& e )
                        { return underlying_functionspace_type::New( e ); } );
    }

    /**
     * construct a product of n spaces from a vector of n different spaces of
     * the same type
     */
    ProductSpace( std::vector<underlying_functionspace_ptrtype> const& m )
        : super( m.begin(), m.end() ),
          M_nspaces( m.size() ),
          M_props( M_nspaces, "Ibc" )
    {
    }

    /**
     * construct a product of n identical spaces
     */
    ProductSpace( int n, underlying_functionspace_ptrtype X )
        : super( 1, X ),
          M_nspaces( n ),
          M_props( M_nspaces, "Ibc" )
    {
    }

    /**
     * @return the number of function spaces in the product
     */
    int numberOfSpaces() const { return M_nspaces; }

    //! \return the total number of degrees of freedom
    size_type nDof() const
    {
        if ( same_mesh )
            return M_nspaces * this->front()->nDof();
        else
            return std::accumulate( this->begin(), this->end(), size_type( 0 ), []( auto i, auto const& e )
                                    { return i + e.nDof(); } );
    }

    //! \return the number of degrees of freedom owned by the process
    size_type nLocalDof() const
    {
        if ( same_mesh )
            return M_nspaces * this->front()->nLocalDof();
        else
            return std::accumulate( this->begin(), this->end(), size_type( 0 ), []( auto i, auto const& e )
                                    { return i + e.nLocalDof(); } );
    }

    worldcomm_type const& worldComm() { return this->front()->worldComm(); }
    worldcomm_type const& worldComm() const { return this->front()->worldComm(); }
    mesh_ptrtype mesh() const { return this->front()->mesh(); }

    underlying_functionspace_ptrtype& operator[]( int i ) { return same_mesh ? this->front() : this->at( i ); }
    underlying_functionspace_ptrtype const& operator[]( int i ) const { return same_mesh ? this->front() : this->at( i ); }

    void setProperties( std::initializer_list<std::string> s )
    {
        M_props = s;
    }
    void setProperties( std::vector<std::string> const& s )
    {
        M_props = s;
    }
    bool isIbcSpace( int i ) const { return M_props[i] == "Ibc"; }
    bool hasIbcSpace() const
    {
        auto f = std::find( M_props.begin(), M_props.end(), "Ibc" );
        return f != std::end( M_props );
    }
    bool hasOtherThanIbcSpace() const
    {
        auto f = std::find_if( M_props.begin(), M_props.end(), []( std::string const& s )
                               { return s != "Ibc"; } );
        return f != std::end( M_props );
    }
    class Element : public BlocksBaseVector<double>, FunctionSpaceBase::ElementBase
    {
      public:
        using super = BlocksBaseVector<double>;
        using value_type = typename underlying_functionspace_type::value_type;
        using underlying_element_type = typename underlying_functionspace_type::element_type;
        using underlying_element_ptrtype = typename underlying_functionspace_type::element_ptrtype;
        using local_interpolant_type = typename underlying_element_type::local_interpolant_type;
        using local_interpolants_type = typename underlying_element_type::local_interpolants_type;

        Element() = default;
        Element( Element const& ) = default;
        Element( Element&& ) = default;

        Element& operator=( Element const& ) = default;
        Element& operator=( Element&& ) = default;

        Element( functionspace_type& X )
            : super( X ), M_fspace( X ) {}

        underlying_element_type&
        operator[]( int n1 )
        {
            return dynamic_cast<underlying_element_type&>( *( ( *this )( n1, 0 ) ) );
        }
        underlying_element_type const&
        operator[]( int n1 ) const
        {
            return dynamic_cast<underlying_element_type const&>( *( ( *this )( n1, 0 ) ) );
        }
        functionspace_type& functionSpace() { return M_fspace; }
        functionspace_type const& functionSpace() const { return M_fspace; }
        // void zero() { super::zero(); }
        functionspace_type M_fspace;
    };
    using element_type = Element;
    using element_ptrtype = std::shared_ptr<element_type>;
    Element element()
    {
        Element u( *this );
        return u;
    }
    std::shared_ptr<Element> elementPtr() { return std::make_shared<Element>( *this ); }
    size_type M_nspaces;
    std::vector<std::string> M_props;
};

template <typename SpaceT, bool same_mesh = true>
using dyn_product_space_t = ProductSpace<SpaceT, same_mesh>;
template <typename SpaceT, bool same_mesh = true>
using dyn_product_space_ptr_t = std::shared_ptr<ProductSpace<SpaceT, same_mesh>>;
template <typename SpaceT, bool same_mesh = true>
using dyn_product_space_element_t = typename ProductSpace<SpaceT, same_mesh>::element_type;
template <typename SpaceT, bool same_mesh = true>
using dyn_product_space_element_ptr_t = typename ProductSpace<SpaceT, same_mesh>::element_ptrtype;

template <typename SpaceT, bool same_mesh = true>
ProductSpace<SpaceT, same_mesh>
dynProduct( int n, typename SpaceT::mesh_ptrtype mesh )
{
    return ProductSpace<SpaceT, same_mesh>( n, mesh );
}
template <typename SpaceT, bool same_mesh = true>
std::shared_ptr<ProductSpace<SpaceT, same_mesh>>
dynProductPtr( int n, typename SpaceT::mesh_ptrtype mesh )
{
    return std::make_shared<ProductSpace<SpaceT, same_mesh>>( n, mesh );
}

template <typename SpaceT, bool same_mesh = true>
std::shared_ptr<ProductSpace<SpaceT, same_mesh>>
dynProductPtr( int n, std::shared_ptr<SpaceT> const& s )
{
    return std::make_shared<ProductSpace<SpaceT, same_mesh>>( n, s );
}

template <typename X, bool SM>
bool operator==( ProductSpace<X, SM> const& x, ProductSpace<X, SM> const& y )
{
    return x == y;
}

template <typename X, bool SM>
bool operator!=( ProductSpace<X, SM> const& x, ProductSpace<X, SM> const& y )
{
    return !( x == y );
}

template <typename... SpaceList>
class ProductSpaces : public ProductSpacesBase
{
  public:
    using tuple_spaces_type = hana::tuple<SpaceList...>;

    // using value_type = typename decay_type<decltype(super[0_c])>::value_type;
    using value_type = double;
    using functionspace_type = ProductSpaces<SpaceList...>;

    using datamap_ptrtype = std::shared_ptr<DataMap<>>;

    static inline constexpr int nSpaces = sizeof...( SpaceList );

    ProductSpaces( SpaceList... l )
        : M_tupleSpaces( l... ),
          M_dm()
    {
        std::vector<datamap_ptrtype> subdm;
        hana::for_each( M_tupleSpaces, [&]( auto& x )
                        { subdm.push_back( x->mapPtr() ); } );
        M_dm = std::make_shared<DataMap<>>( subdm, this->worldCommPtr() );
    }
    constexpr int numberOfSpaces() const { return hana::size( M_tupleSpaces ); }

    tuple_spaces_type const& tupleSpaces() const { return M_tupleSpaces; }
    tuple_spaces_type& tupleSpaces() { return M_tupleSpaces; }

    //! \return the total number of degrees of freedom
    size_type nDof() const
    {
        return hana::fold_left( M_tupleSpaces, 0, [&]( size_type s, auto& e )
                                { return s + e->nDof(); } );
    }
    //! \return the number of degrees of freedom owned by the process
    size_type nLocalDof() const
    {
        return hana::fold_left( M_tupleSpaces, 0, [&]( size_type s, auto& e )
                                { return s + e->nLocalDof(); } );
    }

    template <typename N>
    decltype( auto )
    operator[]( N const& n1 ) const
    {
        return M_tupleSpaces[n1];
    }
    template <typename N>
    decltype( auto )
    operator[]( N const& n1 )
    {
        return M_tupleSpaces[n1];
    }

    /**
     * @brief get the WorldComm
     * @return worldComm
     */
    auto worldComm() const { return mesh()->worldComm(); }

    /**
     * @brief get the shared pointer to the WorldComm
     *
     * @return shared_ptr to worldComm
     */
    auto worldCommPtr() const { return mesh()->worldCommPtr(); }

    /**
     * @brief get the mesh of
     *
     * @return mesh shared_ptr
     */
    auto mesh() const { return M_tupleSpaces[0_c]->mesh(); }

    class Element : public BlocksBaseVector<double>, FunctionSpaceBase::ElementBase
    {
      public:
        using super = BlocksBaseVector<double>;
        static const int nspaces = sizeof...( SpaceList );
        Element() = default;
        Element( Element const& ) = default;
        Element( Element&& ) = default;

        Element& operator=( Element const& ) = default;
        Element& operator=( Element&& ) = default;

        Element( functionspace_type& X )
            : super( X ), M_fspace( X ) { this->buildVector(); }

        template <typename N>
        decltype( auto )
        operator[]( N const& n1 ) const
        {
            return dynamic_cast<decltype( M_fspace[n1]->element() ) const&>( *( ( *this )( n1, 0 ) ) );
        }
        template <typename N>
        decltype( auto )
        operator[]( N const& n1 )
        {
            return dynamic_cast<decltype( M_fspace[n1]->element() )&>( *( ( *this )( n1, 0 ) ) );
        }

        template <typename N>
        decltype( auto )
        operator()( N const& n1 ) const
        {
            if constexpr ( is_shared_ptr_v<std::decay_t<decltype( M_fspace[n1] )>> )
                return dynamic_cast<decltype( M_fspace[n1]->element() ) const&>( *( super::operator()( int( n1 ), 0 ) ) );
            else
                return dynamic_cast<decltype( M_fspace[n1].element() ) const&>( *( super::operator()( int( n1 ), 0 ) ) );
        }

        template <typename N>
        decltype( auto )
        operator()( N const& n1 )
        {
            if constexpr ( is_shared_ptr_v<std::decay_t<decltype( M_fspace[n1] )>> )
                return dynamic_cast<decltype( M_fspace[n1]->element() )&>( *( super::operator()( int( n1 ), 0 ) ) );
            else
                return dynamic_cast<decltype( M_fspace[n1].element() )&>( *( super::operator()( int( n1 ), 0 ) ) );
        }

        functionspace_type& functionSpace() { return M_fspace; }
        functionspace_type const& functionSpace() const { return M_fspace; }
        // void zero() { super::zero(); }
        functionspace_type& M_fspace;
    };

    Element element( std::string name = "U" )
    {
        Element u( *this );
        return u;
    }
    using element_type = Element;
    using element_ptrtype = std::shared_ptr<element_type>;

  private:
    tuple_spaces_type M_tupleSpaces;
    datamap_ptrtype M_dm;
};

template <typename... SList1, typename... SList2>
bool operator==( ProductSpaces<SList1...> const& x, ProductSpaces<SList2...> const& y )
{
    return x.tupleSpaces() == y.tupleSpaces();
}
template <typename... SList1, typename... SList2>
bool operator!=( ProductSpaces<SList1...> const& x, ProductSpaces<SList2...> const& y )
{
    return !( x == y );
}

//!
//! class mixing dynamic and compile-time space product
//!
template <typename T, typename... SpaceList>
class ProductSpaces2 : public ProductSpacesBase
{
  public:
    using tuple_spaces_type = typename mpl::if_<mpl::is_void_<T>,
                                                mpl::identity<hana::tuple<SpaceList...>>,
                                                mpl::identity<hana::tuple<SpaceList..., std::shared_ptr<ProductSpace<T, true>>>>>::type::type; // hana::if_(is_void(T),,>);

    // using value_type = typename decay_type<decltype(super[0_c])>::value_type;
    using value_type = double;
    using functionspace_type = ProductSpaces2<T, SpaceList...>;
    using mesh_type = typename decay_type<T>::mesh_type;
    using mesh_ptrtype = typename decay_type<T>::mesh_ptrtype;

    ProductSpaces2() = default;
    ProductSpaces2( std::shared_ptr<ProductSpace<T, true>> const& p, SpaceList... l )
        : M_tupleSpaces( l..., p ) {}

    int numberOfSpaces() const { return int( hana::size( M_tupleSpaces ) ) + hana::back( M_tupleSpaces )->numberOfSpaces() - 1; }

    tuple_spaces_type const& tupleSpaces() const { return M_tupleSpaces; }
    tuple_spaces_type& tupleSpaces() { return M_tupleSpaces; }

    //! \return the total number of degrees of freedom
    // size_type nDof() const { return hana::fold_left( *this, 0, [&](size_type s, auto& e ) { return s + e->nDof(); } ); }
    //! \return the number of degrees of freedom owned by the process
    // size_type nLocalDof() const { return hana::fold_left( *this, 0, [&](size_type s, auto& e ) { return s + e->nLocalDof(); } ); }

    template <typename N>
    decltype( auto )
    operator[]( N const& n1 ) const
    {
        return M_tupleSpaces[n1];
    }
    template <typename N>
    decltype( auto )
    operator[]( N const& n1 )
    {
        return M_tupleSpaces[n1];
    }

    class Element : public BlocksBaseVector<double>, FunctionSpaceBase::ElementBase
    {
      public:
        using super = BlocksBaseVector<double>;
        static const int nspaces = sizeof...( SpaceList ) + 1;

        Element() = default;
        Element( Element const& ) = default;
        Element( Element&& ) = default;

        Element& operator=( Element const& ) = default;
        Element& operator=( Element&& ) = default;

        Element( functionspace_type& X )
            : super( X ), M_fspace( X ) {}

        template <typename N>
        decltype( auto )
        operator[]( N const& n1 ) const
        {
            return dynamic_cast<decltype( M_fspace[n1]->element() ) const&>( *( ( *this )( n1, 0 ) ) );
        }
        template <typename N>
        decltype( auto )
        operator[]( N const& n1 )
        {
            return dynamic_cast<decltype( M_fspace[n1]->element() )&>( *( ( *this )( n1, 0 ) ) );
        }
#if 0
        template<typename N>
        //decltype(auto)
        auto const &
            operator()( N const& n1, int i = 0 ) const
            {
                return hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(M_fspace[n1])>>{},
                                 [&] (auto&& x, auto& s) {
                                     return dynamic_cast<decltype((*s)[i]->element()) const&>(x);
                                 },
                                 [&] (auto&& x, auto& s){
                                     return dynamic_cast<decltype(s->element()) const&>(x);
                                 })(*(super::operator()(int(n1)+i,0)), M_fspace[n1] );

            }
#endif
        template <typename N>
        decltype( auto )
        operator()( N const& n1, int i )
        {
            return dynamic_cast<decltype( ( *M_fspace[n1] )[i]->element() )&>( *super::operator()( int( n1 ) + i, 0 ) );
        }
        //!
        //! \return the n-th element
        //!
        template <typename N>
        decltype( auto )
        operator()( N const& n1 )
        {
            return dynamic_cast<decltype( M_fspace[n1]->element() )&>( *super::operator()( int( n1 ), 0 ) );
        }
        //!
        //! \return the i-th element of the n-th dynamic function space product of the product space
        //!
        template <typename N>
        decltype( auto )
        operator()( N const& n1, int i ) const
        {
            return dynamic_cast<decltype( ( *M_fspace[n1] )[i]->element() ) const&>( *super::operator()( int( n1 ) + i, 0 ) );
        }
        //!
        //! \return the n-th element
        //!
        template <typename N>
        decltype( auto )
        operator()( N const& n1 ) const
        {
            return dynamic_cast<decltype( M_fspace[n1]->element() ) const&>( *super::operator()( int( n1 ), 0 ) );
        }

        functionspace_type& functionSpace() { return M_fspace; }
        functionspace_type const& functionSpace() const { return M_fspace; }
        template <typename N>
        decltype( auto ) functionSpace( N const& n ) { return M_fspace[n]; }
        template <typename N>
        decltype( auto ) functionSpace( N const& n ) const { return M_fspace[n]; }
        // void zero() { super::zero(); }
        functionspace_type M_fspace;
    };

    Element element()
    {
        Element u( *this );
        return u;
    }
    std::shared_ptr<Element> elementPtr() { return std::make_shared<Element>( *this ); }
    using element_type = Element;
    using element_ptrtype = std::shared_ptr<element_type>;

  private:
    tuple_spaces_type M_tupleSpaces;
};

template <typename T, typename... SList1, typename... SList2>
bool operator==( ProductSpaces2<T, SList1...> const& x, ProductSpaces2<T, SList2...> const& y )
{
    return x.tupleSpaces() == y.tupleSpaces();
}
template <typename T, typename... SList1, typename... SList2>
bool operator!=( ProductSpaces2<T, SList1...> const& x, ProductSpaces2<T, SList2...> const& y )
{
    return !( x == y );
}

template <typename PS>
constexpr auto is_product_spaces( PS&& ps )
{
    return hana::traits::is_base_of( ps );
}
template <typename... SpaceList>
std::shared_ptr<ProductSpaces<SpaceList...>>
productPtr( SpaceList... spaces )
{
    return std::make_shared<ProductSpaces<SpaceList...>>( spaces... );
}
template <typename... SpaceList>
using product_spaces_t = ProductSpaces<SpaceList...>;
template <typename... SpaceList>
using product_spaces_ptr_t = std::shared_ptr<ProductSpaces<SpaceList...>>;
template <typename... SpaceList>
using product_spaces_element_t = typename ProductSpaces<SpaceList...>::element_type;
template <typename... SpaceList>
using product_spaces_element_ptr_t = typename ProductSpaces<SpaceList...>::element_ptrtype;

/**
 * compile time product of function spaces
 * \param spaces expansion of arbitrary number of spaces
 * \return the compile time product
 */
template <typename... SpaceList>
ProductSpaces<SpaceList...>
product( SpaceList... spaces )
{
    return ProductSpaces<SpaceList...>( spaces... );
}

template <typename SpaceT, typename... SpaceList>
using dyn_product_spaces_t = ProductSpaces2<SpaceT, SpaceList...>;
template <typename SpaceT, typename... SpaceList>
using dyn_product_spaces_ptr_t = std::shared_ptr<ProductSpaces2<SpaceT, SpaceList...>>;
template <typename SpaceT, typename... SpaceList>
using dyn_product_spaces_element_t = typename ProductSpaces2<SpaceT, SpaceList...>::element_type;
template <typename SpaceT, typename... SpaceList>
using dyn_product_spaces_element_ptr_t = typename ProductSpaces2<SpaceT, SpaceList...>::element_ptrtype;

/**
 * mixed compile time and runtime product of function space
 * \param spaces expansion of arbitrary number of function spaces
 * \param ps dynamic product of function spaces of the same type
 * \return the space product
 */
template <typename T, typename... SpaceList>
ProductSpaces2<T, SpaceList...>
product2( std::shared_ptr<ProductSpace<T, true>> const& ps, SpaceList... spaces )
{
    return ProductSpaces2<T, SpaceList...>( ps, spaces... );
}
template <typename T, typename... SpaceList>
std::shared_ptr<ProductSpaces2<T, SpaceList...>>
product2Ptr( std::shared_ptr<ProductSpace<T, true>> const& ps, SpaceList... spaces )
{
    return std::make_shared<ProductSpaces2<T, SpaceList...>>( ps, spaces... );
}
/**
 * mixed compile time and runtime product of function space
 * \param spaces expansion of arbitrary number of function spaces
 * \param ps dynamic product of function spaces of the same type
 * \return the space product
 */
template <typename T, typename... SpaceList>
ProductSpaces2<T, SpaceList...>
product2( SpaceList... spaces, std::shared_ptr<ProductSpace<T, true>> const& ps )
{
    return product2( ps, spaces... );
}

#if 0
namespace detail {

namespace product {

template<typename T>
struct clean
{
    typedef typename boost::remove_reference<typename boost::remove_const<typename boost::remove_pointer<T>::type>::type >::type type;
};
template<typename Space>
struct GetMesh
{
    typedef typename clean<typename Space::element_type>::type::mesh_type type;
};
template<typename Space>
struct GetBasis
{
    typedef typename clean<typename Space::element_type>::type::basis_type type;
};
template<typename Space>
struct GetPeriodicity
{
    typedef typename clean<typename Space::element_type>::type::periodicity_type type;
};
template<typename Space>
struct GetMortar
{
    typedef typename clean<typename Space::element_type>::type::mortar_0_type type;
};

template<typename... SpaceList>
struct Product
{
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetMesh<mpl::_1>, mpl::back_inserter<meshes<> > >::type mesh_type;
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetBasis<mpl::_1>, mpl::back_inserter<bases<> > >::type basis_type;
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetPeriodicity<mpl::_1>, mpl::back_inserter<Periodicity<> > >::type periodicity_type;
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetMortar<mpl::_1>, mpl::back_inserter<mortars<> > >::type mortar_type;

    typedef FunctionSpace<typename mpl::front<mesh_type>::type,basis_type,double,periodicity_type,mortar_type> type;
    typedef std::shared_ptr<type> ptrtype;
};
} }

template<typename... SpaceList>
typename Feel::detail::product::Product<SpaceList...>::ptrtype
product( SpaceList... spaces )
{
    typedef typename Feel::detail::product::Product<SpaceList...> product_type;
    typedef typename product_type::type space_type;
    typedef typename product_type::ptrtype space_ptrtype;
    std::cout << "sizeof " << sizeof...(SpaceList) << "\n";
    space_ptrtype Rh( space_type::NewFromList( spaces... ) );
    return Rh;
}

#endif

} // namespace Feel
#endif
