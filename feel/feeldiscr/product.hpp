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



template<typename... SpaceList>
class ProductSpaces : public  hana::tuple<SpaceList...>, ProductSpacesBase
{
public:
    using super =  hana::tuple<SpaceList...>;
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
