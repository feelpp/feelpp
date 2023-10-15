/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 20 April 2022

 Copyright (C) 2022 University of Strasbourg 
 Copyright (C) 2022 Feel++ Consortium
 
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
#include <tuple>
#include <type_traits>
#include <boost/mp11/utility.hpp>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/meshbase.hpp>

namespace Feel 
{
// clang-format off
template<typename MeshType, int MESH_ENTITIES,std::enable_if_t<std::is_base_of_v<MeshBase<>,decay_type<std::remove_pointer_t<MeshType>>>,int> = 0>
using entities_reference_wrapper_t = boost::mp11::mp_if_c<MESH_ENTITIES==MESH_ELEMENTS,
                                                          elements_reference_wrapper_t<MeshType>,
                                                          boost::mp11::mp_if_c<MESH_ENTITIES==MESH_FACES,
                                                                        faces_reference_wrapper_t<MeshType>,
                                                                        boost::mp11::mp_if_c<MESH_ENTITIES==MESH_EDGES,
                                                                                      edges_reference_wrapper_t<MeshType>,
                                                                                      points_reference_wrapper_t<MeshType>
                                                                                     > 
                                                                        > 
                                                          >;
// clang-format on

/**
 * @brief RangeBase class
 * 
 */
template<typename IndexT = uint32_type>
class FEELPP_EXPORT RangeBase : public CommObject 
{
public:
    using index_t = IndexT;
    RangeBase() = default;
    RangeBase( RangeBase const& ) = default;
    RangeBase( RangeBase && ) = default;
    RangeBase& operator=( RangeBase const& ) = default;
    RangeBase& operator=( RangeBase && ) = default;
    virtual ~RangeBase() = default;
    RangeBase( std::shared_ptr<MeshBase<index_t>> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b ) {}

    std::shared_ptr<MeshBase<index_t>> meshBase() { return M_mesh_base; }
    std::shared_ptr<MeshBase<index_t>> meshBase() const { return M_mesh_base; }

protected:
    std::shared_ptr<MeshBase<index_t>> M_mesh_base;

};

template <typename MeshType, int MESH_ENTITIES, std::enable_if_t<std::is_base_of_v<MeshBase<>, decay_type<std::remove_pointer_t<MeshType>>>, int> = 0>
class FEELPP_EXPORT Range
    : public entities_reference_wrapper_t<MeshType, MESH_ENTITIES>, public RangeBase<typename decay_type<std::remove_pointer_t<MeshType>>::index_t>
{
    using super = entities_reference_wrapper_t<MeshType,MESH_ENTITIES>;
    using super_range = RangeBase<typename decay_type<std::remove_pointer_t<MeshType>>::index_t>;
  public:
    
    using mesh_t = decay_type<std::remove_pointer_t<MeshType>>;
    using mesh_ptr_t = std::shared_ptr<const mesh_t>;
    using index_type = typename mesh_t::index_type;
    using index_t = index_type;

    static constexpr int mesh_entities = MESH_ENTITIES;
    static constexpr int nDim = mesh_t::nDim;
    static constexpr int nRealDim = mesh_t::nRealDim;


    Range() = default;
    Range( Range const& e ) = default;
    Range( Range && e ) = default;
    Range( super const& e ) : super( e ) {}
    Range( super && e ) : super( std::move( e ) ) {}
    ~Range() = default;

    static constexpr int entities() { return MESH_ENTITIES; }
    static constexpr bool isOnElements() { return mesh_entities == MESH_ELEMENTS;  }
    static constexpr bool isOnFaces() { return mesh_entities == MESH_FACES;  }
    static constexpr bool isOnEdges() { return mesh_entities == MESH_EDGES;  }
    static constexpr bool isOnFacets() { return mesh_entities == MESH_FACES; }
    static constexpr bool isOnPoints() { return mesh_entities == MESH_POINTS;  }

    Range( super const& er, MeshType const& m, int marker, int pid )
        : super( er ), M_mesh( unwrap_ptr(m).shared_from_this() ), M_marker( marker ), M_pid( pid )
    {
    }
    Range( super const& er, std::shared_ptr<const MeshType> const& m, int marker, int pid )
        : super( er ), M_mesh( m ), M_marker( marker ), M_pid( pid )
    {
    }
    Range( super const& er, std::shared_ptr<std::remove_const_t<MeshType>> const& m, int marker, int pid ): super( er ), M_mesh(m), M_marker(marker), M_pid(pid)
    {
    }
    Range( super && er, MeshType const& m, int marker, int pid )
        : super( std::move( er ) ), M_mesh( unwrap_ptr(m).shared_from_this() ), M_marker( marker ), M_pid( pid )
    {
    }
    Range( super && er, std::shared_ptr<std::remove_const_t<MeshType>> const& m, int marker, int pid )
        : super( std::move( er ) ), M_mesh( m ), M_marker( marker ), M_pid( pid )
    {
    }
    Range( super && er, std::shared_ptr<const MeshType> const& m, int marker, int pid )
        : super( std::move( er ) ), M_mesh( m ), M_marker( marker ), M_pid( pid )
    {
    }

    Range& operator=( Range const& ) = default;
    Range& operator=( Range && ) = default;

    auto begin() { return this->template get<1>(); }
    auto end() { return this->template get<2>(); }
    auto const& begin() const { return this->template get<1>(); }
    auto const& end() const { return this->template get<2>(); }

    int marker() const { return M_marker; }
    void setMarker( int m ) { M_marker = m; } 
    bool hasMarker() const { return M_marker != -1; }

    void setProcessId( int id ) { M_pid = id; } 
    int pid() const { return M_pid; }
    bool hasPid() const { return M_pid != -1; } 

private:
    mesh_ptr_t M_mesh;
    int M_marker;
    int M_pid;
};

template <typename... Tv>
auto
range( Tv&&... v )
{
    auto args = NA::make_arguments( std::forward<Tv>( v )... );
    using mesh_t = decay_type<decltype( args.get( _mesh ) )>;
    using entities_t = typename boost::tuples::element<0,decay_type<decltype( args.get( _range ) )>>::type;
    return Range<mesh_t,entities_t::value>{ args.get( _range ),
                                            args.get( _mesh ),                  
                                            args.get_else( _marker1, -1 ),
                                            args.get_else( _pid, -1 ) };
}

/**
 * @brief get the begin iterator of a range
 * 
 * @tparam MeshType the mesh type
 * @tparam MESH_ENTITIES the mesh entities
 * @param range the range of elements
 * @return auto 
 */
template<typename MeshType, int MESH_ENTITIES>
auto begin( Range<MeshType,MESH_ENTITIES> &range )
{
    return range.begin();
}

template<typename MeshType, int MESH_ENTITIES>
auto end( Range<MeshType,MESH_ENTITIES> &range )
{
    return range.end();
}
template<typename MeshType, int MESH_ENTITIES>
auto begin( Range<MeshType,MESH_ENTITIES> const&range )
{
    return range.begin();
}

template<typename MeshType, int MESH_ENTITIES>
auto end( Range<MeshType,MESH_ENTITIES> const&range )
{
    return range.end();
}

}