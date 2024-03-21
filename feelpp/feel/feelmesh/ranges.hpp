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
#include <iostream>
#include <fmt/ostream.h>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelcore/enums.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feeldiscr/mesh_fwd.hpp>

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
using container_reference_wrapper_t = boost::mp11::mp_if_c<MESH_ENTITIES==MESH_ELEMENTS,
                                                          typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype,
                                                          boost::mp11::mp_if_c<MESH_ENTITIES==MESH_FACES,
                                                                        typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype,
                                                                        boost::mp11::mp_if_c<MESH_ENTITIES==MESH_EDGES,
                                                                                      typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype,
                                                                                      typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype
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

#if 0
    template <typename GeoShape, typename T = double, int Tag = 0>
    RangeBase( std::shared_ptr<Mesh<GeoShape,T,Tag,IndexT>> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    template <typename GeoShape, typename T = double, int Tag = 0>
    RangeBase( std::shared_ptr<Mesh<GeoShape,T,Tag,IndexT> const> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    template <typename GeoShape, typename T = double>
    RangeBase( std::shared_ptr<MeshStructured<GeoShape,T,IndexT> const> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}

    template <typename GeoShape, typename T = double>
    RangeBase( std::shared_ptr<MeshStructured<GeoShape,T,IndexT>> const& b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b.get() ) {}
#else

    template <typename MeshType>
    RangeBase(std::shared_ptr<MeshType> const& b)
        : CommObject(b->worldCommPtr()),
          M_mesh_base(b.get())
    {}

    template <typename MeshType>
    RangeBase(std::shared_ptr<const MeshType> const& b)
        : CommObject(b->worldCommPtr()),
          M_mesh_base(b.get())
    {}

    template <typename MeshType, std::enable_if_t<!std::is_pointer_v<MeshType>, int> = 0>
    RangeBase(MeshType const& b)
        : CommObject(b.worldCommPtr()),
          M_mesh_base(&b)
    {}

    template <typename MeshType, std::enable_if_t<std::is_pointer_v<MeshType>, int> = 0>
    RangeBase(MeshType const b)
        : CommObject(b->worldCommPtr()),
          M_mesh_base(b)
    {}

#endif
    RangeBase( MeshBase<index_t> const* b ) : CommObject( b->worldCommPtr() ), M_mesh_base( b ) {}

    MeshBase<index_t> const* meshBase() const { return M_mesh_base; }

    template<typename MeshType>
    void setMesh( std::shared_ptr<MeshType> const& m )
    {
        this->setWorldCommPtr( m->worldCommPtr() );
        M_mesh_base = m.get();
    }


protected:
    //std::weak_ptr<MeshBase<index_t>> M_mesh_base;
    //std::weak_ptr<MeshBase<index_t> const> M_mesh_base;
    MeshBase<index_t> const* M_mesh_base;

};

template<typename MeshType>
std::shared_ptr<MeshType> shared_from_this( std::shared_ptr<MeshType> const& m )
{
    return m;
}
template<typename MeshType, typename = std::enable_if_t<std::is_base_of_v<std::enable_shared_from_this<MeshType>,MeshType>>>
std::shared_ptr<MeshType const> shared_from_this( MeshType const& m )
{
    return m.shared_from_this();
}
template<typename MeshType, typename = std::enable_if_t<!std::is_base_of_v<std::enable_shared_from_this<MeshType>,MeshType>>>
MeshType const* shared_from_this( MeshType const& m )
{
    return &m;
}
template <typename MeshType, int MESH_ENTITIES, std::enable_if_t<std::is_base_of_v<MeshBase<>, decay_type<std::remove_pointer_t<MeshType>>>, int> = 0>
class FEELPP_EXPORT Range
    : public RangeBase<typename decay_type<std::remove_pointer_t<MeshType>>::index_t>
    //public entities_reference_wrapper_t<MeshType, MESH_ENTITIES>, public RangeBase<typename decay_type<std::remove_pointer_t<MeshType>>::index_t>
{
    using super = entities_reference_wrapper_t<MeshType,MESH_ENTITIES>;
    using super_range = RangeBase<typename decay_type<std::remove_pointer_t<MeshType>>::index_t>;
  public:

    using mesh_non_const_t = std::remove_const_t<decay_type<std::remove_pointer_t<MeshType>>>;
    using mesh_t = mesh_non_const_t;
    using mesh_const_t = std::add_const_t<mesh_t>;
    using mesh_ptr_non_const_t = std::shared_ptr<mesh_non_const_t>;
    using mesh_ptr_const_t = std::shared_ptr<const mesh_non_const_t>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    using index_type = typename mesh_t::index_type;
    using index_t = index_type;
    using container_t = typename boost::tuples::element<3,super>::type::element_type;
    using container_ptr_t = typename boost::tuples::element<3,super>::type;
    using idim_t = typename boost::tuples::template element<0, super>::type;
    using iterator_t = typename boost::tuples::template element<1, super>::type;
    using element_t = typename  boost::unwrap_reference<typename iterator_t::value_type>::type;
    using value_t = typename mesh_t::value_type;
    using mesh_support_ptr_t = std::shared_ptr < MeshSupport<decay_type<MeshType> >>;
    using mesh_support_t = MeshSupport<decay_type<MeshType>>;

    static constexpr int mesh_entities = MESH_ENTITIES; // idim_t
    static constexpr int nDim = mesh_t::nDim;
    static constexpr int nRealDim = mesh_t::nRealDim;


    Range()  : super_range(), mesh_(nullptr), mesh_support_(), cont_( std::make_shared<container_t>() )
    {
    }
    Range( Range const& e ) : super_range(e.mesh_), mesh_( e.mesh_ ), mesh_support_( e.mesh_support_ ), cont_( e.cont_ ) {}
    Range( Range && e ) = default;
    ~Range() = default;

    // Constructor for when OtherMeshType is a std::shared_ptr<MeshType>
    template <typename OtherMeshType>
    Range( const Range<OtherMeshType, MESH_ENTITIES>& other,
           std::enable_if_t<std::is_same_v<std::shared_ptr<MeshType>, std::decay_t<std::remove_pointer_t<OtherMeshType>>>, int> = 0 )
        : super_range( other.mesh() ), mesh_( other.mesh() ), mesh_support_( other.meshSupport() ), cont_( other.container() ) 
        {
            //std::cout << fmt::format( "Range<OtherMeshType> shared constructor") << std::endl;
        }

    // Constructor for when OtherMeshType is not a std::shared_ptr<MeshType>
    template <typename OtherMeshType>
    Range( const Range<OtherMeshType, MESH_ENTITIES>& other,
           std::enable_if_t <is_simplex_v<MeshType> && !std::is_same_v<std::shared_ptr<MeshType>, std::decay_t<std::remove_pointer_t<OtherMeshType>>>, int> = 0 )
            : super_range( other.mesh() ), mesh_( other.mesh() ), mesh_support_( other.meshSupport() ), cont_( other.container() )
    {
            // std::cout << fmt::format("Range<OtherMeshType> not shared constructor") << std::endl;
    }

    template <typename OtherMeshType>
    Range( const Range<OtherMeshType, MESH_ENTITIES>& other,
           std::enable_if_t<!is_simplex_v<MeshType> && !std::is_same_v<std::shared_ptr<MeshType>, std::decay_t<std::remove_pointer_t<OtherMeshType>>>, int> = 0 )
        : super_range( other.mesh() ), mesh_( nullptr ), mesh_support_( nullptr ), cont_( other.container() )
    {
        // std::cout << fmt::format("Range<OtherMeshType> not shared constructor") << std::endl;
    }

    static constexpr idim_t idim() { return idim_t(); }
    static constexpr int iDim() { return idim_t::value; }
    static constexpr int entities() { return MESH_ENTITIES; }
    static constexpr bool isOnElements() { return mesh_entities == MESH_ELEMENTS;  }
    static constexpr bool isOnFaces() { return mesh_entities == MESH_FACES;  }
    static constexpr bool isOnEdges() { return mesh_entities == MESH_EDGES;  }
    static constexpr bool isOnFacets() { return mesh_entities == MESH_FACES; }
    static constexpr bool isOnPoints() { return mesh_entities == MESH_POINTS;  }

    // Constructors handling both MeshType and MeshType const
    Range(mesh_non_const_t const& m) : Range(shared_from_this(m)), mesh_(&m)
    {}

    Range(mesh_ptr_non_const_t const& m) : super_range(m), mesh_(m.get()), cont_(std::make_shared<container_t>())
    {}

    Range(mesh_ptr_const_t const& m) : super_range(m), mesh_(m.get()), cont_(std::make_shared<container_t>())
    {}

    Range(mesh_non_const_t const* m) : super_range(m), mesh_(m), cont_(std::make_shared<container_t>())
    {}

    // Constructor for Mesh
    Range( super const& er, mesh_non_const_t const& m )
        : super_range( unwrap_ptr(m).shared_from_this() ), mesh_(&m), cont_( er.template get<3>() )
    {}

    // Constructor for Mesh*
    Range( super const& er, mesh_non_const_t const* m)
        : super_range( unwrap_ptr(m).shared_from_this() ), mesh_(m), cont_( er.template get<3>() )
    {}

    // Constructor for std::shared_ptr<Mesh>
    Range( super const& er, std::shared_ptr<std::remove_const_t<MeshType>> const& m )
        : super_range( m ), mesh_(m.get()), cont_( er.template get<3>() )
    {}

    // Constructor for std::shared_ptr<const Mesh>
    Range( super const& er, std::shared_ptr<const MeshType> const& m)
        : super_range( m ), mesh_(m.get()), cont_( er.template get<3>() )
    {}
    Range( mesh_non_const_t const& m, iterator_t beg, iterator_t end )
        : Range(shared_from_this(m)), mesh_(&m), cont_(std::make_shared<container_t>(beg,end))
    {}
    Range(mesh_ptr_non_const_t const& m, iterator_t beg, iterator_t end) : super_range(m), mesh_(m.get()), cont_(std::make_shared<container_t>(beg,end))
    {}

    Range(mesh_ptr_const_t const& m, iterator_t beg, iterator_t end) : super_range(m), mesh_(m.get()), cont_(std::make_shared<container_t>(beg,end))
    {}

    Range(mesh_non_const_t const* m, iterator_t beg, iterator_t end) : super_range(m), mesh_(m), cont_(std::make_shared<container_t>(beg,end))
    {}


    Range& operator=( Range const& ) = default;
    Range& operator=( Range && ) = default;

    bool isEmpty() const { return cont_->empty(); }

    template<int N>
    auto get() 
    { 
        if constexpr ( N==1 ) 
            return cont_->begin();
        else 
            return cont_->end();
    }
    auto begin() { return cont_->begin(); }
    auto end() { return cont_->end(); }
    auto begin() const { return cont_->begin(); }
    auto end() const { return cont_->end(); }

    element_t const& front() const { return boost::unwrap_ref(cont_->front()); }
    element_t const& back() const { return boost::unwrap_ref(cont_->back()); }

    int size() const { return cont_->size(); }
    
    container_ptr_t const& container() const { return cont_; }
    container_ptr_t container() { return cont_; }
    void clear()
    {
        this->container()->clear();
    }
    void push_back( element_t const& e )
    {
        this->container()->push_back( boost::cref( e ) );
    }
    void shrink_to_fit() { this->container()->shrink_to_fit(); }

    auto mesh() const { return mesh_; }

    /**
     * @brief Check if the range has mesh support
     */
    bool hasMeshSupport() const { return mesh_support_ != nullptr; }

    /**
     * @brief Get the mesh support
     */
    mesh_support_ptr_t meshSupport() const { return mesh_support_; }

    void setMeshSupport( mesh_support_ptr_t const& ms ) { mesh_support_ = ms; }
private:
    //std::weak_ptr<const mesh_t> mesh_;
    mesh_non_const_t const* mesh_;
    mesh_support_ptr_t mesh_support_;
    container_ptr_t cont_;
};




template <typename... Tv>
auto
range( Tv&&... v )
{
    auto args = NA::make_arguments( std::forward<Tv>( v )... );
    using mesh_t = decay_type<decltype( args.get( _mesh ) )>;
    using entities_t = typename boost::tuples::element<0,decay_type<decltype( args.get( _range ) )>>::type;
    return Range<mesh_t,entities_t::value>{ args.get( _range ),
                                            args.get( _mesh ) };
}


template <typename RangeT>
inline constexpr bool is_range_v = std::is_base_of_v<RangeBase<>,decay_type<RangeT>>;
//inline constexpr bool is_range_v = is_tuple_v<decay_type<RangeT>> || std::is_base_of_v<RangeBase<>,decay_type<RangeT>>;//std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<range_iterators_t<RangeT>>::iterator_category>;

template <typename Type>
struct value_type_trait<Type, std::enable_if_t<is_range_v<Type>>>
{
    using type = typename decay_type<Type>::value_t;
};

/**
 * @brief Specialization for Range
 * 
 * @tparam RangeType 
 */
template <typename RangeType>
struct element_type_helper<RangeType, std::enable_if_t<is_range_v<RangeType>>> {
    using type = typename decay_type<RangeType>::element_t;
    using ptrtype = std::shared_ptr<type>;
};

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

template <typename MeshType, int MESH_ENTITIES, typename = std::enable_if_t<std::is_base_of_v<MeshBase<>, decay_type<std::remove_pointer_t<MeshType>>>, int>>
std::ostream& operator<<(std::ostream& os, const Range<MeshType, MESH_ENTITIES>& range)
{
    os << "Range of " << MESH_ENTITIES << " entities with " << nelements( range) << " elements.";
    return os;
}

/**
 * @brief 
 */
template <typename RangeType, std::enable_if_t<is_range_v<RangeType>,int> = 0>
std::vector<decay_type<RangeType>>
partitionRange(RangeType&& range, int nParts) 
{
    using range_t = decay_type<RangeType>;
    std::vector<range_t> partitions;
    partitions.reserve(nParts);

    auto partSize = std::forward<RangeType>(range).size() / nParts;
    auto partBegin = std::forward<RangeType>(range).begin();

    for (int i = 0; i < nParts; ++i) 
    {
        auto start = partBegin;
        std::advance(partBegin, partSize); // Advance partBegin for the next iteration
        auto end = (i == nParts - 1) ? std::forward<RangeType>(range).end() : partBegin;
        
        partitions.emplace_back( std::forward<RangeType>(range).mesh(), start, end );
    }

    return partitions;
}


} // namespace Feel