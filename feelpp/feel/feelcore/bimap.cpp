/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2020-05-04

   Copyright (C) 2020 Université de Strasbourg

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
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wredeclared-class-member"
#endif
#include <boost/bimap.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/bimap.hpp>

namespace Feel
{
template<typename LeftT, typename RightT>
class bimap<LeftT,RightT>::Pimpl
{
public:
    using left_t = LeftT;
    using right_t = RightT;
    using bm_type = boost::bimap< left_t, right_t >;
    
    bm_type bm;
};
template<typename LeftT, typename RightT>
bimap<LeftT, RightT>::bimap() 
  : pimpl_( new Pimpl ) 
{
}
template<typename LeftT, typename RightT>
bimap<LeftT, RightT>::~bimap() 
{
}
template<typename LeftT, typename RightT>
bimap<LeftT, RightT>&
bimap<LeftT, RightT>::operator=( std::map<left_t, right_t> const& m ) 
{
    for( auto i : m )   
        this->insert( i );
    return *this;
}
template<typename LeftT, typename RightT>
int
bimap<LeftT, RightT>::size() const
{
    return pimpl_->bm.size(); 
}
template<typename LeftT, typename RightT>
void
bimap<LeftT, RightT>::insert( std::pair<left_t,right_t>&& p ) 
{
    pimpl_->bm.insert( typename Pimpl::bm_type::value_type( p.first, p.second ) ); 
}
template<typename LeftT, typename RightT>
std::optional<typename bimap<LeftT, RightT>::right_t>
bimap<LeftT, RightT>::leftFind( left_t i ) const
{
    if ( auto l = pimpl_->bm.left.find( i ); l != pimpl_->bm.left.end() )
        return l->second;
    else  
        return std::nullopt;
}
template<typename LeftT, typename RightT>
std::optional<typename bimap<LeftT, RightT>::left_t>
bimap<LeftT, RightT>::rightFind( right_t i ) const
{
    if ( auto l = pimpl_->bm.right.find( i ); l != pimpl_->bm.right.end() )
        return l->second;
    else  
        return std::nullopt;
}
template<typename LeftT, typename RightT>
typename bimap<LeftT, RightT>::right_t
bimap<LeftT, RightT>::leftAt( left_t i ) const
{
    return pimpl_->bm.left.at( i );
}
template<typename LeftT, typename RightT>
typename bimap<LeftT, RightT>::left_t
bimap<LeftT, RightT>::rightAt( right_t i ) const
{
    return pimpl_->bm.left.at( i );
} 
      
template class bimap<uint32_type, uint32_type>;
template class bimap<int32_type, int32_type>;
}