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
#if !defined(FEELPP_PRODUCT_HPP)
#define FEELPP_PRODUCT_HPP 1

#include <feel/feelalg/vectorblock.hpp>

#include <boost/fusion/view/joint_view.hpp>


namespace Feel {

constexpr auto is_void = hana::integral(hana::metafunction<std::is_void>);

template<typename T, bool same_mesh = false>
class ProductSpace : public std::vector<boost::shared_ptr<decay_type<T>>>, ProductSpaceBase
{
public:
    using super =  std::vector<boost::shared_ptr<decay_type<T>>>;
    //using value_type = typename decay_type<decltype(super[0_c])>::value_type;
    using value_type = double;
    using functionspace_type = ProductSpace<T,same_mesh>;
    using underlying_functionspace_type = decay_type<T>;
    using underlying_functionspace_ptrtype = boost::shared_ptr<underlying_functionspace_type>;
    using mesh_type = typename underlying_functionspace_type::mesh_type;
    using mesh_ptrtype = typename underlying_functionspace_type::mesh_ptrtype;

    ProductSpace( int n, mesh_ptrtype m  )
        :
        super( 1, underlying_functionspace_type::New( m ) ),
        M_nspaces(n)
        {

        }

    int numberOfSpaces() const { return M_nspaces; }

    //! \return the total number of degrees of freedom
    size_type nDof() const
        {
            if ( same_mesh )
                return M_nspaces*this->front()->nDof();
            else
                return std::accumulate( this->begin(), this->end(), size_type(0), []( auto i, auto const& e ) { return i+e.nDof(); } );

        }

    //! \return the number of degrees of freedom owned by the process
    size_type nLocalDof() const
        {
            if ( same_mesh )
                return M_nspaces*this->front()->nLocalDof();
            else
                return std::accumulate( this->begin(), this->end(), size_type(0), []( auto i, auto const& e ) { return i+e.nLocalDof(); } );

        }

    underlying_functionspace_ptrtype& operator[]( int i ) { return same_mesh?this->front():this->at(i); }
    underlying_functionspace_ptrtype const& operator[]( int i ) const { return same_mesh?this->front():this->at(i); }

    class Element : public BlocksBaseVector<double>, FunctionSpaceBase::ElementBase
    {
    public:
        using super = BlocksBaseVector<double>;
        using underlying_element_type = typename underlying_functionspace_type::element_type;
        using underlying_element_ptrtype = typename underlying_functionspace_type::element_ptrtype;
        Element() = default;
        Element( Element const& ) = default;
        Element( Element && ) = default;

        Element& operator=( Element const& ) = default;
        Element& operator=( Element && ) = default;

        Element( functionspace_type X ): super(X),  M_fspace(X) {}

        underlying_element_type&
        operator[]( int n1 )
            {
                return dynamic_cast<underlying_element_type &>(*((*this)(n1,0)));
            }
        underlying_element_type const&
        operator[]( int n1 ) const
            {
                return dynamic_cast<underlying_element_type const&>(*((*this)(n1,0)));
            }
        functionspace_type functionSpace() { return M_fspace; }
        //void zero() { super::zero(); }
        functionspace_type M_fspace;
    };
    using element_type = Element;
    using element_ptrtype = boost::shared_ptr<element_type>;
    Element element() { Element u( *this ); return u; }
    boost::shared_ptr<Element> elementPtr() { return boost::make_shared<Element>( *this ); }
    size_type M_nspaces;
};

template<typename... SpaceList>
class ProductSpaces : public  hana::tuple<SpaceList...>, ProductSpacesBase
{
public:
#if 0
    using super = typename mpl::if_<mpl::is_void_<T>,
                                    mpl::identity<hana::tuple<SpaceList...>>,
                                    mpl::identity<hana::tuple<SpaceList...,ProductSpace<T,true>>>
                                    >::type::type ;//hana::if_(is_void(T),,>);
#else
    using super = hana::tuple<SpaceList...>;
#endif
    //using value_type = typename decay_type<decltype(super[0_c])>::value_type;
    using value_type = double;
    using functionspace_type = ProductSpaces<SpaceList...>;

    ProductSpaces( SpaceList... l ) : super( l... ){}
    int numberOfSpaces() const { return hana::size( *this ); }

    //! \return the total number of degrees of freedom
    size_type nDof() const { return hana::fold_left( *this, 0, [&](size_type s, auto& e ) { return s + e->nDof(); } ); }
    //! \return the number of degrees of freedom owned by the process
    size_type nLocalDof() const { return hana::fold_left( *this, 0, [&](size_type s, auto& e ) { return s + e->nLocalDof(); } ); }


    class Element : public BlocksBaseVector<double>, FunctionSpaceBase::ElementBase
    {
    public:
        using super = BlocksBaseVector<double>;
        Element() = default;
        Element( Element const& ) = default;
        Element( Element && ) = default;

        Element& operator=( Element const& ) = default;
        Element& operator=( Element && ) = default;

        Element( functionspace_type X ): super(X),  M_fspace(X) {}

        template<typename N>
        decltype(auto)
        operator[]( N const& n1 ) const
            {
                return dynamic_cast<decltype(M_fspace[n1]->element()) const&>(*((*this)(n1,0)));
            }
        functionspace_type functionSpace() { return M_fspace; }
        //void zero() { super::zero(); }
        functionspace_type M_fspace;
    };

    Element element() { Element u( *this ); return u; }
    using element_type = Element;
    using element_ptrtype = boost::shared_ptr<element_type>;
};

template<typename T,typename... SpaceList>
class ProductSpaces2 : public  hana::tuple<SpaceList...,boost::shared_ptr<ProductSpace<T,true>>>, ProductSpacesBase
{
public:

    using super = typename mpl::if_<mpl::is_void_<T>,
                                    mpl::identity<hana::tuple<SpaceList...>>,
                                    mpl::identity<hana::tuple<SpaceList...,boost::shared_ptr<ProductSpace<T,true>>>>
                                    >::type::type ;//hana::if_(is_void(T),,>);

    //using value_type = typename decay_type<decltype(super[0_c])>::value_type;
    using value_type = double;
    using functionspace_type = ProductSpaces2<T,SpaceList...>;
    using mesh_type = typename decay_type<T>::mesh_type;
    using mesh_ptrtype = typename decay_type<T>::mesh_ptrtype;

    ProductSpaces2( boost::shared_ptr<ProductSpace<T,true>> const& p, SpaceList... l ) : super( l..., p) {}
    int numberOfSpaces() const { return int(hana::size( *this ))+hana::back(*this)->numberOfSpaces()-1; }

    //! \return the total number of degrees of freedom
    //size_type nDof() const { return hana::fold_left( *this, 0, [&](size_type s, auto& e ) { return s + e->nDof(); } ); }
    //! \return the number of degrees of freedom owned by the process
    //size_type nLocalDof() const { return hana::fold_left( *this, 0, [&](size_type s, auto& e ) { return s + e->nLocalDof(); } ); }

    class Element : public BlocksBaseVector<double>, FunctionSpaceBase::ElementBase
    {
    public:
        using super = BlocksBaseVector<double>;
        Element() = default;
        Element( Element const& ) = default;
        Element( Element && ) = default;

        Element& operator=( Element const& ) = default;
        Element& operator=( Element && ) = default;

        Element( functionspace_type& X ): super(X),  M_fspace(X) {}

        template<typename N>
        decltype(auto)
            operator[]( N const& n1 ) const
            {
                return dynamic_cast<decltype(M_fspace[n1]->element()) const&>(*((*this)(n1,0)));
            }
        template<typename N>
        decltype(auto)
            operator[]( N const& n1 )
            {
                return dynamic_cast<decltype(M_fspace[n1]->element()) const&>(*((*this)(n1,0)));
            }
        template<typename N>
        decltype(auto)
            operator()( N const& n1, int i = 0 ) const
            {
                return hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(M_fspace[n1])>>{},
                                 [&] (auto&& x) {
                                     return dynamic_cast<decltype((*M_fspace[n1])[i]->element()) const&>(x);
                                 },
                                 [&] (auto&& x){
                                     return dynamic_cast<decltype(M_fspace[n1]->element()) const&>(x);
                                 })(*(super::operator()(int(n1)+i,0)));

            }

        template<typename N>
        decltype(auto)
            operator()( N const& n1, int i = 0 )
            {
                return hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(M_fspace[n1])>>{},
                                 [&] (auto&& x, auto & s) {
                                     return dynamic_cast<decltype((*s)[i]->element())&>(std::forward<decltype(x)>(x));
                                 },
                                 [&] (auto&& x, auto &s){
                                     return dynamic_cast<decltype(s->element()) &>(std::forward<decltype(x)>(x));
                                 })(*(super::operator()(int(n1)+i,0)), M_fspace[n1] );

            }
        functionspace_type functionSpace() { return M_fspace; }
        //void zero() { super::zero(); }
        functionspace_type M_fspace;
    };

    Element element() { Element u( *this ); return u; }
    boost::shared_ptr<Element> elementPtr() { return boost::make_shared<Element>( *this );  }
    using element_type = Element;
    using element_ptrtype = boost::shared_ptr<element_type>;
};

template<typename PS>
constexpr auto is_product_spaces( PS&& ps )
{
    return hana::traits::is_base_of(ps);
}
template<typename... SpaceList>
std::unique_ptr<ProductSpaces<SpaceList...>>
productPtr( SpaceList... spaces )
{
    return std::make_unique<ProductSpaces<SpaceList...>>( spaces... );
}
template<typename... SpaceList>
ProductSpaces<SpaceList...>
product( SpaceList... spaces )
{
    return ProductSpaces<SpaceList...>( spaces... );
}

template<typename T, typename... SpaceList>
ProductSpaces2<T,SpaceList...>
product2( boost::shared_ptr<ProductSpace<T,true>> const& ps, SpaceList... spaces )
{
    return ProductSpaces2<T,SpaceList...>( ps, spaces... );
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
    typedef boost::shared_ptr<type> ptrtype;
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

} // Feel++
#endif
