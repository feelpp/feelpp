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
#include <feel/feelcore/enums.hpp>
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
template<typename MeshType, int MESH_ENTITIES,std::enable_if_t<std::is_base_of_v<MeshBase<>,decay_type<std::remove_pointer_t<MeshType>>>,int> = 0>
using entities_t = boost::mp11::mp_if_c<MESH_ENTITIES==MESH_ELEMENTS,
                                        typename MeshTraits<MeshType>::element_type,
                                        boost::mp11::mp_if_c<MESH_ENTITIES==MESH_FACES,
                                                            typename MeshTraits<MeshType>::face_type,
                                                            boost::mp11::mp_if_c<MESH_ENTITIES==MESH_EDGES,
                                                                                typename MeshTraits<MeshType>::edge_type,
                                                                                typename MeshTraits<MeshType>::point_type
                                                                               >
                                                            >
                                        >;
template <typename GeoShape, typename T, int Tag, typename IndexT>
class Mesh;
template <typename GeoShape, typename T, typename IndexT>
class MeshStructured;
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

    template <typename GeoShape, typename T = double, int Tag = 0>
    RangeBase( std::shared_ptr<Mesh<GeoShape,T,Tag,IndexT>> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    template <typename GeoShape, typename T = double, int Tag = 0>
    RangeBase( std::shared_ptr<Mesh<GeoShape,T,Tag,IndexT> const> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    template <typename GeoShape, typename T = double>
    RangeBase( std::shared_ptr<MeshStructured<GeoShape,T,IndexT> const> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    template <typename GeoShape, typename T = double>
    RangeBase( std::shared_ptr<MeshStructured<GeoShape,T,IndexT>> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    RangeBase( MeshBase<index_t> const* b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b ) {}

    MeshBase<index_t> const* meshBase() const { return M_mesh_base; }

protected:
    //std::weak_ptr<MeshBase<index_t>> M_mesh_base;
    //std::weak_ptr<MeshBase<index_t> const> M_mesh_base;
    MeshBase<index_t> const* M_mesh_base;

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
    using idim_t = typename boost::tuples::template element<0, super>::type;
    using iterator_t = typename boost::tuples::template element<1, super>::type;
    using element_t = entities_t<mesh_t, MESH_ENTITIES>;

    static constexpr int mesh_entities = MESH_ENTITIES; // idim_t
    static constexpr int nDim = mesh_t::nDim;
    static constexpr int nRealDim = mesh_t::nRealDim;


    Range() = default;
    Range( Range const& e ) = default;
    Range( Range && e ) =  default;
    ~Range() = default;
    template <typename OtherMeshType>
    Range( Range<OtherMeshType, MESH_ENTITIES> const& other ) : super(other), super_range(other) {}

    static constexpr int entities() { return MESH_ENTITIES; }
    static constexpr bool isOnElements() { return mesh_entities == MESH_ELEMENTS;  }
    static constexpr bool isOnFaces() { return mesh_entities == MESH_FACES;  }
    static constexpr bool isOnEdges() { return mesh_entities == MESH_EDGES;  }
    static constexpr bool isOnFacets() { return mesh_entities == MESH_FACES; }
    static constexpr bool isOnPoints() { return mesh_entities == MESH_POINTS;  }

    // Constructor for Mesh
    Range( super const& er, MeshType const& m, int marker, int pid )
        : super( er ), super_range( unwrap_ptr(m).shared_from_this() ), M_marker( marker ), M_pid( pid )
    {}

    // Constructor for Mesh*
    Range( super const& er, MeshType* m, int marker, int pid )
        : super( er ), super_range( unwrap_ptr(m).shared_from_this() ), M_marker( marker ), M_pid( pid )
    {}

    // Constructor for std::shared_ptr<Mesh>
    Range( super const& er, std::shared_ptr<std::remove_const_t<MeshType>> const& m, int marker, int pid )
        : super( er ), super_range( m ), M_marker( marker ), M_pid( pid )
    {}

    // Constructor for std::shared_ptr<const Mesh>
    Range( super const& er, std::shared_ptr<const MeshType> const& m, int marker, int pid )
        : super( er ), super_range( m ), M_marker( marker ), M_pid( pid )
    {}

    Range& operator=( Range const& ) = default;
    Range& operator=( Range && ) = default;

    bool isEmpty() const { return begin() == end(); }

    auto begin() { return this->template get<1>(); }
    auto end() { return this->template get<2>(); }
    auto const& begin() const { return this->template get<1>(); }
    auto const& end() const { return this->template get<2>(); }

    element_t const& front() const { return *begin(); }
    element_t const& back() const { return *std::prev( end() ); }


    int marker() const { return M_marker; }
    void setMarker( int m ) { M_marker = m; }
    bool hasMarker() const { return M_marker != -1; }

    void setProcessId( int id ) { M_pid = id; }
    int pid() const { return M_pid; }
    bool hasPid() const { return M_pid != -1; }

    mesh_t const& mesh() const { return  dynamic_cast<mesh_t const&>( *this->meshBase() ); }
private:
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

template <typename RangeT>
inline constexpr bool is_range_v = std::is_base_of_v<RangeBase<>,decay_type<RangeT>>;
//inline constexpr bool is_range_v = is_tuple_v<decay_type<RangeT>> || std::is_base_of_v<RangeBase<>,decay_type<RangeT>>;//std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<range_iterators_t<RangeT>>::iterator_category>;


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

/**
 * @brief Utility function for MPI global reduction.
 *
 * @param local The local data to be reduced.
 * @param global Whether to perform global reduction.
 * @param r The range providing the WorldComm.
 *
 * @return Reduced result.
 */
template<typename T, typename RangeT, std::enable_if_t<is_range_v<RangeT>,int> = 0>
T
globalReduce(T local, bool global, RangeT const& r)
{
    T result = local;
    if (global)
        mpi::all_reduce(r.worldComm(),
                        local,
                        result,
                        std::plus<T>());
    return result;
}

/**
 * @brief Lambda function factory to check if a cell is a ghost.
 *
 * @tparam RangeT Type of the range.
 *
 * @return Lambda function.
 */
template <typename RangeT, std::enable_if_t<is_range_v<RangeT>,int> = 0>
auto isGhostCell()
{
    return [](auto const& cell)
    {
        if constexpr (RangeT::mesh_entities == ElementsType::MESH_FACES)
            return cell.isGhostFace();
        else if constexpr (RangeT::mesh_entities == ElementsType::MESH_ELEMENTS)
            return cell.isGhostCell();
        return false;
    };
}

/**
 * @brief Lambda function factory to count cells excluding ghosts.
 *
 * @param is_ghost The function to check if a cell is ghost.
 *
 * @return Lambda function.
 */
template <typename Predicate>
auto countWithoutGhost(Predicate is_ghost)
{
    return [is_ghost](size_type x, auto const& c)
    {
        auto const& cell = unwrap_ref(c);
        return !is_ghost(cell) ? x + 1 : x;
    };
}

/**
 * @brief Compute the number of elements in a range.
 *
 * @tparam RangeT Type of the range.
 * @param r The range.
 * @param global Whether to compute the number globally.
 *
 * @return Number of elements.
 */
template <typename RangeT, std::enable_if_t<is_filter_v<RangeT>,int> = 0>
size_type
nelements(RangeT const& r, bool global = false)
{
    size_type d = std::accumulate(r.begin(), r.end(), 0, countWithoutGhost(isGhostCell<RangeT>()));
    return globalReduce(d, global, r);
}

/**
 * @brief [Deprecated] Compute the number of elements in a list of ranges.
 *
 * @tparam MeshType Type of the mesh.
 * @tparam Entities The entity type.
 * @param r The list of ranges.
 * @param global Whether to compute the number globally.
 *
 * @return Number of elements.
 */
template<typename MeshType, int Entities>
FEELPP_DEPRECATED size_type
nelements(std::list<Range<MeshType,Entities>> const& r, bool global)
{
    return nelements(r, global);
}

/**
 * @brief Compute the number of elements in a collection of ranges.
 *
 * @tparam CollectionOfRangeT Collection type (e.g., list, vector).
 * @param its The collection of ranges.
 * @param global Whether to compute the number globally.
 *
 * @return Number of elements.
 */
template <typename CollectionOfRangeT, std::enable_if_t<is_range_v<typename CollectionOfRangeT::value_type>,int> = 0>
size_type
nelements(CollectionOfRangeT const& its, bool global = false)
{
    size_type d = 0;
    std::for_each(its.begin(), its.end(),
                  [&d](auto const& t)
                  {
                      d += nelements(t, false);
                  });
    return globalReduce(d, global, its.front());  // Assuming all ranges in the collection have the same WorldComm
}

/**
 * @brief Compute the number of elements in a range for a given zone.
 *
 * @tparam RangeT Type of the range.
 * @param its The range.
 * @param z The zone.
 *
 * @return Number of elements.
 */
template <typename RangeT, std::enable_if_t<is_range_v<RangeT>,int> = 0>
size_type
nelements(RangeT const& its, Zone const& z)
{
    return nelements(its, (z == Zone::GLOBAL));
}


} // namespace Feel