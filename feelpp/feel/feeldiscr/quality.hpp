/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 09 Mar 2020

 Copyright (C) 2020 Feel++ Consortium

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
#pragma once


#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/measures.hpp>

namespace Feel
{

/**
 * Return default triangle/tetrahedra quality (eta) of a given cell.
 
 The quality measure relates the area of the triangle (a)
 to its edge lengths (l1, l2, l3).
 \f$
 \eta = \frac{4\sqrt{3}a}{l_1^2 + l_2^2 + l_3^2}
 \f$
 */
template<typename MeshType>
std::shared_ptr<Pdh_element_t<MeshType,0>>
etaQ( Pdh_ptrtype<MeshType,0> const& Xh, std::shared_ptr<MeshType> const& m )
{
    using value_type=value_t<MeshType>;
    
    auto q = Xh->elementPtr();

    local_interpolant_t<decltype(Xh)> Ihloc( 1 );
    for ( auto const& eltWrap : elements( m ) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if constexpr ( dimension_v<MeshType> == 2 )
            // use formula 7 in D Field paper
            Ihloc( 0 ) = 4*math::sqrt(3.)*elt.measure()/elt.faceMeasures().array().square().sum();
        else if constexpr ( dimension_v<MeshType> == 3 )
            // use formula 33 in D Field 
            Ihloc( 0 ) = math::pow(2.,3)*math::pow(3.,8)*math::sqrt(3.)*math::pow(elt.measure(),4)/math::pow(elt.faceMeasures().array().square().sum(),3);
        q->assign( elt, Ihloc );
    }
    return q;
    
}
template<typename MeshType>
std::shared_ptr<Pdh_element_t<MeshType,0>>
etaQ( std::shared_ptr<MeshType> const& m )
{
    auto Xh = Pdh<0>( m );
    return etaQ( Xh, m );
}

/**
 * Return the normalized shape ratio (NSR) for a given element.
 * it is referred as the radius ratio, as it is described by
 * the ratio between the inradius (r) and the circumradius (R).
 * \f$\rho = \frac{dr}{R}\f$ where d is the dimension
 *
 */
template<typename MeshType>
std::shared_ptr<Pdh_element_t<MeshType,0>>
nsrQ( Pdh_ptrtype<MeshType,0> const& Xh, std::shared_ptr<MeshType> const& m )
{
    using value_type=value_t<MeshType>;
    
    auto q = Xh->elementPtr();

    local_interpolant_t<decltype(Xh)> Ihloc( 1 );
    for ( auto const& eltWrap : elements( m ) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if constexpr ( dimension_v<MeshType> == 2 )
        {
            // use formula 3 in D Field paper
            auto l = elt.faceMeasures();
            auto r = 2 * elt.measure() / l.sum(); // inradius
            auto R = 0.25 * l.prod() / elt.measure();  // circumradius
            Ihloc( 0 ) = 2*r/R;
        }
        else if constexpr ( dimension_v<MeshType> == 3 )
        {
            Ihloc( 0 ) = 0;
        }
        q->assign( elt, Ihloc );
    }
    return q;
    
}
template<typename MeshType>
std::shared_ptr<Pdh_element_t<MeshType,0>>
nsrQ( std::shared_ptr<MeshType> const& m )
{
    auto Xh = Pdh<0>( m );
    return nsrQ( Xh, m );
}

}
