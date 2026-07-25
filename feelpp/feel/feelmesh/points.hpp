/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2005-2006 EPFL
    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
   @file points.hpp
   @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   @date 2005-09-03
   @brief Points container class (C++20/23 modernized)
 */
#ifndef FEELPP_MESH_POINTS_HPP
#define FEELPP_MESH_POINTS_HPP

#include <concepts>
#include <functional>
#include <type_traits>
#include <unordered_map>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmesh/geocoordinatestorage.hpp>
#include <feel/feelmesh/geoelement.hpp>

namespace Feel
{
/// \cond detail
/*!
  \class Points
  \brief Points container class

  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @see
*/
template <uint16_type nDim, typename T = double, typename IndexT = uint32_type, typename SubFace = SubFaceOfNone<0> >
class Points
{
  public:
    /** @name Typedefs
     */
    //@{
    using index_type = IndexT;
    using size_type = index_type;
    using point_type = GeoElement0D<nDim, SubFace, T, IndexT>;
    /** Mesh-owned sparse metadata storage for points. */
    using sparse_metadata_storage_type =
        typename point_type::sparse_metadata_storage_type;
    /** Aligned mesh-owned coordinate storage. */
    using coordinate_storage_type =
        GeoCoordinateStorage<nDim, T, index_type>;
    /** Read-only zero-copy Eigen coordinate map. */
    using coordinate_const_map_type =
        typename coordinate_storage_type::const_map_type;
    /** Mutable zero-copy Eigen coordinate map. */
    using coordinate_map_type =
        typename coordinate_storage_type::map_type;

    using points_type = std::unordered_map<size_type, point_type>;

    using point_iterator = typename points_type::iterator;
    using point_const_iterator = typename points_type::const_iterator;

    using points_reference_wrapper_type = std::vector<boost::reference_wrapper<point_type const>>;
    using points_reference_wrapper_ptrtype = std::shared_ptr<points_reference_wrapper_type>;
    using point_reference_wrapper_iterator = typename points_reference_wrapper_type::iterator;
    using point_reference_wrapper_const_iterator = typename points_reference_wrapper_type::const_iterator;

    using ordered_points_reference_wrapper_type = std::vector<boost::reference_wrapper<point_type>>;
    using ordered_point_reference_wrapper_iterator = typename ordered_points_reference_wrapper_type::iterator;
    using ordered_point_reference_wrapper_const_iterator = typename ordered_points_reference_wrapper_type::const_iterator;

    using point_interprocess_map_type = std::unordered_map<index_type, std::set<rank_type>>;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Points( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        : M_worldComm( worldComm ),
          M_points(),
          M_sparseMetadata(),
          M_coordinateStorage(),
          M_needToOrderPoints( false )
    {
    }

    Points( Points const& f )
        : M_worldComm( f.M_worldComm ),
          M_points( f.M_points ),
          M_sparseMetadata(),
          M_coordinateStorage(),
          M_needToOrderPoints( false )
    {
        this->bindSparseMetadata();
        this->buildOrderedPoints();
        this->synchronizePointCoordinates();
    }

    virtual ~Points() {}

    void clear()
    {
        DVLOG( 1 ) << "deleting points...\n";
        M_orderedPoints.clear();
        M_points.clear();
        M_sparseMetadata.clear();
        M_coordinateStorage.clear();
        M_needToOrderPoints = false;
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    Points& operator=( Points const& e )
    {
        if ( this != &e )
        {
            M_worldComm = e.M_worldComm;
            M_sparseMetadata.clear();
            M_points = e.M_points;
            this->bindSparseMetadata();
            this->buildOrderedPoints();
            this->synchronizePointCoordinates();
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @brief Return the mutable compatibility point container.
     *
     * Direct coordinate edits require `synchronizePointCoordinates()` before
     * compact-coordinate consumers run. Prefer `editPointCoordinates()` for
     * transactional zero-copy coordinate updates.
     *
     * @return mutable point container
     */
    points_type& points()
    {
        return M_points;
    }

    /**
     * \return the points container
     */
    points_type const& points() const
    {
        return M_points;
    }

    /**
     * @brief Return populated point sparse-sidecar counts.
     */
    [[nodiscard]] GeoEntitySparseMetadataStatistics
    pointSparseMetadataStatistics() const noexcept
    {
        return M_sparseMetadata.statistics();
    }

    /** @brief Return compact point-connectivity occupancy counters. */
    [[nodiscard]] GeoEntityConnectivityStatistics
    pointConnectivityStatistics() const noexcept
    {
        return M_sparseMetadata.connectivityStatistics();
    }

    /**
     * @brief Return an immutable zero-copy Eigen map of point coordinates.
     *
     * The map is invalidated by insertion, erase, clear, or explicit
     * coordinate synchronization.
     */
    [[nodiscard]] coordinate_const_map_type
    pointCoordinates() const noexcept
    {
        return M_coordinateStorage.coordinates();
    }

    /**
     * @brief Return a mutable zero-copy Eigen map of point coordinates.
     *
     * Mutating this map updates the compact coordinate storage only. Prefer
     * `editPointCoordinates()` when compatibility point objects must observe
     * the edit immediately, or call `synchronizeLegacyPointCoordinates()`.
     */
    [[nodiscard]] coordinate_map_type
    pointCoordinates() noexcept
    {
        return M_coordinateStorage.coordinates();
    }

    /**
     * @brief Return the dense local index of an external point identifier.
     * @param id external point identifier
     */
    [[nodiscard]] uint32_type pointLocalIndex( index_type id ) const
    {
        return M_coordinateStorage.localIndex( id );
    }

    /** @brief Return coordinate-storage occupancy and lookup counters. */
    [[nodiscard]] GeoCoordinateStorageStatistics
    pointCoordinateStorageStatistics() const noexcept
    {
        return M_coordinateStorage.statistics();
    }

    /**
     * @brief Rebuild aligned coordinates from the legacy point facade.
     *
     * This compatibility boundary is required after code mutates legacy
     * point nodes directly. It can be removed once points are non-owning
     * views of the compact coordinate storage.
     */
    void synchronizePointCoordinates()
    {
        M_coordinateStorage.clear();
        M_coordinateStorage.reserve( M_orderedPoints.size() );
        for ( auto const& pointReference : M_orderedPoints )
        {
            auto const& point = unwrap_ref( pointReference );
            auto const& node = point.node();
            (void)M_coordinateStorage.append(
                point.id(),
                std::span<T const>(
                    node.data().begin(), node.size() ) );
        }
    }

    /**
     * @brief Copy compact coordinate edits back to legacy point objects.
     *
     * This compatibility operation is temporary while point objects still
     * own coordinate values independently of the compact storage.
     */
    void synchronizeLegacyPointCoordinates()
    {
        auto coordinates = M_coordinateStorage.coordinates();
        for ( uint32_type localIndex = 0;
              localIndex < M_coordinateStorage.size();
              ++localIndex )
        {
            auto& point = this->point(
                M_coordinateStorage.externalId( localIndex ) );
            for ( uint16_type component = 0;
                  component < nDim; ++component )
                point( component ) =
                    coordinates( component, localIndex );
        }
    }

    /**
     * @brief Apply an edit to compact coordinates and update legacy points.
     *
     * The editor receives a mutable zero-copy Eigen map. Synchronization is
     * performed on normal return and during exception unwinding.
     *
     * @tparam Editor callable accepting `coordinate_map_type&`
     * @param editor coordinate edit operation
     */
    template <typename Editor>
        requires std::invocable<Editor&, coordinate_map_type&>
    void editPointCoordinates( Editor&& editor )
    {
        auto coordinates = M_coordinateStorage.coordinates();
        try
        {
            std::invoke(
                std::forward<Editor>( editor ), coordinates );
        }
        catch ( ... )
        {
            this->synchronizeLegacyPointCoordinates();
            throw;
        }
        this->synchronizeLegacyPointCoordinates();
    }

    virtual bool isEmpty() const
    {
        return M_points.empty();
    }
    bool isBoundaryPoint( point_type const& e ) const
    {
        return e.isOnBoundary();
    }
    bool isBoundaryPoint( size_type const& id ) const
    {
        auto itFindPt = M_points.find( id );
        if ( itFindPt == M_points.end() )
            return false;
        return itFindPt->isOnBoundary();
    }

    point_type& point( size_type i )
    {
        auto itFindPt = M_points.find( i );
        CHECK( itFindPt != M_points.end() ) << " point " << i << "does not found";
        return itFindPt->second;
    }
    point_type const& point( size_type i ) const
    {
        auto itFindPt = M_points.find( i );
        CHECK( itFindPt != M_points.end() ) << " point " << i << "does not found";
        return itFindPt->second;
    }

    point_const_iterator pointIterator( size_type i ) const
    {
        return M_points.find( i );
    }
    point_iterator pointIterator( size_type i )
    {
        return M_points.find( i );
    }

    bool hasPoint( size_type i ) const
    {
        return M_points.find( i ) != M_points.end();
    }

    point_iterator beginPoint()
    {
        return M_points.begin();
    }
    point_const_iterator beginPoint() const
    {
        return M_points.begin();
    }
    point_iterator endPoint()
    {
        return M_points.end();
    }
    point_const_iterator endPoint() const
    {
        return M_points.end();
    }

    ordered_points_reference_wrapper_type const& orderedPoints() const noexcept { return M_orderedPoints; }

    ordered_point_reference_wrapper_iterator beginOrderedPoint()
    {
        return M_orderedPoints.begin();
    }
    ordered_point_reference_wrapper_const_iterator beginOrderedPoint() const
    {
        return M_orderedPoints.begin();
    }
    ordered_point_reference_wrapper_iterator endOrderedPoint()
    {
        return M_orderedPoints.end();
    }
    ordered_point_reference_wrapper_const_iterator endOrderedPoint() const
    {
        return M_orderedPoints.end();
    }

    /**
     * \return the range of iterator \c (begin,end) over the points
     * with marker \p m on processor \p p
     */
    std::pair<point_iterator, point_iterator>
    pointsRange()
    {
        return std::make_pair( M_points.begin(), M_points.end() );
    }

    /**
     * \return the range of iterator \c (begin,end) over the points
     * with marker \p m on processor \p p
     */
    std::pair<point_const_iterator, point_const_iterator>
    pointsRange() const
    {
        return std::make_pair( M_points.begin(), M_points.end() );
    }

    /**
     * \return the range of iterator \c (begin,end) over the points
     * with \c Id \p m
     */
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithId( size_type m ) const
    {
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        if ( this->hasPoint( m ) )
            mypoints->push_back( boost::cref( this->point( m ) ) );
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    /**
     * \return iterator over marked points
     */
    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithAllMarkerByType( uint16_type markerType, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( !point.hasMarkerType( markerType ) )
                continue;
            if ( point.marker().isOff() )
                continue;
            if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
            {
                if ( !Feel::detail::checkPartitionPredicate<EPT>( point, part, M_interprocessPoints ) )
                    continue;
            }
            mypoints->push_back( boost::cref( point ) );
        }
        mypoints->shrink_to_fit();
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    /**
     * \return iterator over marked points
     */
    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( !point.hasMarkerType( markerType ) )
                continue;
            if ( !point.marker( markerType ).hasOneOf( markerFlags ) )
                continue;
            if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
            {
                if ( !Feel::detail::checkPartitionPredicate<EPT>( point, part, M_interprocessPoints ) )
                    continue;
            }
            mypoints->push_back( boost::cref( point ) );
        }
        mypoints->shrink_to_fit();
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    /**
     * \return iterator over marked points
     */
    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
    {
        if ( m == invalid_flag_type_value )
            return this->pointsWithAllMarkerByType<EPT>( markerType, p );
        else
            return this->pointsWithMarkerByType<EPT>( markerType, std::set<flag_type>( {m} ), p );
    }

    /**
     * \return iterator over marked points
     */
    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithMarker( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
    {
        return this->pointsWithMarkerByType<EPT>( 1, m, p );
    }

    /**
     * get iterator over internal points
     *
     *
     * @return iterator over internal points
     */
    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    internalPoints( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( !point.isInternal() )
                continue;
            if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
            {
                if ( !Feel::detail::checkPartitionPredicate<EPT>( point, part, M_interprocessPoints ) )
                    continue;
            }
            mypoints->push_back( boost::cref( point ) );
        }
        mypoints->shrink_to_fit();
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }
    /**
     * get iterator over boundary points
     *
     *
     * @return iterator over boundary points
     */
    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    boundaryPoints( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if ( !point.isOnBoundary() )
                continue;
            if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
            {
                if ( !Feel::detail::checkPartitionPredicate<EPT>( point, part, M_interprocessPoints ) )
                    continue;
            }
            mypoints->push_back( boost::cref( point ) );
        }
        mypoints->shrink_to_fit();
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }

    template <entity_process_t EPT = entity_process_t::LOCAL_ONLY>
    std::tuple<point_reference_wrapper_const_iterator, point_reference_wrapper_const_iterator, points_reference_wrapper_ptrtype>
    pointsWithProcessId( rank_type p = invalid_rank_type_value ) const
    {
        const rank_type part = ( p == invalid_rank_type_value ) ? this->worldCommPoints().localRank() : p;
        points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
        auto it = this->beginOrderedPoint();
        auto en = this->endOrderedPoint();
        for ( ; it != en; ++it )
        {
            auto const& point = unwrap_ref( *it );
            if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
            {
                if ( !Feel::detail::checkPartitionPredicate<EPT>( point, part, M_interprocessPoints ) )
                    continue;
            }
            mypoints->push_back( boost::cref( point ) );
        }
        mypoints->shrink_to_fit();
        return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
    }
#if 0
    template <typename SF = SubFace>
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    interProcessPoints( rank_type neighbor_pid = invalid_rank_type_value,
                        std::enable_if_t< std::is_base_of<SubFaceOfBase, SF>::value >* = nullptr ) const
        {
            bool allNeighbor = ( neighbor_pid == invalid_rank_type_value );
            const rank_type part = this->worldCommPoints().localRank();
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginOrderedPoint();
            auto en = this->endOrderedPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = unwrap_ref( *it );
                if ( !point.isInterProcessDomain() )
                    continue;
                if ( point.partition1() != part )
                    continue;
                if ( !allNeighbor && point.partition2() != neighbor_pid )
                    continue;
                mypoints->push_back( boost::cref( point ) );
            }
            mypoints->shrink_to_fit();
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }

#else

    template <entity_process_t EPT = entity_process_t::ALL>
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    interProcessPoints( rank_type neighbor_pid = invalid_rank_type_value ) const
        {
            bool allNeighbor = ( neighbor_pid == invalid_rank_type_value );
            const rank_type part = this->worldCommPoints().localRank();
            points_reference_wrapper_ptrtype mypoints( new points_reference_wrapper_type );
            auto it = this->beginOrderedPoint();
            auto en = this->endOrderedPoint();
            for ( ; it!=en;++it )
            {
                auto const& point = unwrap_ref( *it );
                auto itFindIP = M_interprocessPoints.find( point.id() );
                if ( itFindIP == M_interprocessPoints.end() )
                    continue;
                if ( !allNeighbor && itFindIP->second.find( neighbor_pid ) == itFindIP->second.end() )
                    continue;
                if constexpr ( EPT == entity_process_t::LOCAL_ONLY || EPT == entity_process_t::GHOST_ONLY || EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
                {
                    if ( !Feel::detail::checkPartitionPredicate<EPT>( point, part, M_interprocessPoints ) )
                        continue;
                }
                mypoints->push_back( boost::cref( point ) );
            }
            mypoints->shrink_to_fit();
            return std::make_tuple( mypoints->begin(), mypoints->end(), mypoints );
        }
#endif

    template <entity_filter_t FF, entity_process_t EPT, typename ... Ts>
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    pointsFilterImpl( Ts&&... ts ) const
        {
            if constexpr ( FF == entity_filter_t::PROCESS_ID )
                return this->pointsWithProcessId<EPT>( std::forward<Ts>( ts )... );
            else if constexpr ( FF == entity_filter_t::MARKER )
                return this->pointsWithMarkerByType<EPT>( std::forward<Ts>( ts )... );
            else if constexpr ( FF == entity_filter_t::ON_BOUNDARY )
                return this->boundaryPoints<EPT>( std::forward<Ts>( ts )... );
            else if constexpr ( FF == entity_filter_t::INTERNAL )
                return this->internalPoints<EPT>( std::forward<Ts>( ts )... );
            return {};
        }
    template <entity_filter_t FF, typename ... Ts>
    std::tuple<point_reference_wrapper_const_iterator,point_reference_wrapper_const_iterator,points_reference_wrapper_ptrtype>
    pointsFilter( entity_process_t ept, Ts&&... ts ) const
        {
            return std::invoke(
                [this,&ept](auto&& ... args)
                    {
                        switch ( ept )
                        {
                        default:
                        case entity_process_t::LOCAL_ONLY:
                            return this->pointsFilterImpl<FF, entity_process_t::LOCAL_ONLY>( std::forward<decltype(args)>(args) ... );
                        case entity_process_t::LOCAL_AND_INTERPROCESS_ONLY:
                            return this->pointsFilterImpl<FF, entity_process_t::LOCAL_AND_INTERPROCESS_ONLY>( std::forward<decltype(args)>(args) ... );
                        case entity_process_t::GHOST_ONLY:
                            return this->pointsFilterImpl<FF, entity_process_t::GHOST_ONLY>( std::forward<decltype(args)>(args) ... );
                        case entity_process_t::ALL:
                            return this->pointsFilterImpl<FF, entity_process_t::ALL>( std::forward<decltype(args)>(args) ... );
                        }
                    },
                std::forward<Ts>( ts )... );
        }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //!
    //! reserve size of container
    //! @param nPoint : the size reserved
    //!
    void reserveNumberOfPoint( size_type nPoint )
    {
        M_points.reserve( nPoint );
        M_orderedPoints.reserve( nPoint );
    }

    /**
     * add a new point in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    std::pair<point_iterator,bool> addPoint( point_type const& f )
    {
        //return M_points.insert( std::make_pair( f.id(), f ) ).first->second;
        auto ret = M_points.emplace( std::make_pair( f.id(), f ) );

        auto& newPoint = ret.first->second;
        if ( ret.second )
        {
            newPoint.attachSparseMetadataStorage( M_sparseMetadata );
            auto const& node = newPoint.node();
            (void)M_coordinateStorage.append(
                newPoint.id(),
                std::span<T const>(
                    node.data().begin(), node.size() ) );
            if ( !M_needToOrderPoints && !M_orderedPoints.empty() && unwrap_ref( M_orderedPoints.back() ).id() > newPoint.id() )
                M_needToOrderPoints = true;
            M_orderedPoints.push_back( boost::ref( newPoint ) );
        }
        return ret;
    }

    /**
     * add a new point in the mesh
     * @param f a new point
     * @return the new point from the list
     */
    std::pair<point_iterator,bool> addPoint( point_type&& f )
    {
        //return M_points.insert( std::make_pair( f.id(), f ) ).first->second;
        auto ret = M_points.emplace(
            std::make_pair( f.id(), std::move( f ) ) );

        auto& newPoint = ret.first->second;
        if ( ret.second )
        {
            newPoint.attachSparseMetadataStorage( M_sparseMetadata );
            auto const& node = newPoint.node();
            (void)M_coordinateStorage.append(
                newPoint.id(),
                std::span<T const>(
                    node.data().begin(), node.size() ) );
            if ( !M_needToOrderPoints && !M_orderedPoints.empty() && unwrap_ref( M_orderedPoints.back() ).id() > newPoint.id() )
                M_needToOrderPoints = true;
            M_orderedPoints.push_back( boost::ref( newPoint ) );
        }
        return ret;
    }

    /**
     * erase point at position \p position
     *
     * @param position \p position is a valid dereferenceable iterator of the index.
     *
     * @return An iterator pointing to the point immediately
     * following the one that was deleted, or \c end() if no such point
     * exists.
     */
    point_iterator erasePoint( point_iterator it )
    {
        size_type erasedId = it->first;
        auto itOrdered = std::find_if(
            M_orderedPoints.begin(), M_orderedPoints.end(),
            [&erasedId]( auto& pointWrap )
            {
                return unwrap_ref( pointWrap ).id() == erasedId;
            } );
        it->second.clearEntityStorage();
        auto itret = M_points.erase( it );
        M_orderedPoints.erase( itOrdered );
        this->synchronizePointCoordinates();
        return itret;
    }

    WorldComm const& worldCommPoints() const
    {
        return *M_worldComm;
    }

    void setWorldCommPoints( worldcomm_ptr_t const& _worldComm )
    {
        M_worldComm = _worldComm;
    }

    //@}

    void updateOrderedPoints()
    {
        if ( !M_needToOrderPoints )
            return;
        std::sort( M_orderedPoints.begin(), M_orderedPoints.end(),
                   []( auto const& a, auto const& b ) -> bool {
                       return unwrap_ref( a ).id() < unwrap_ref( b ).id();
                   } );
        M_needToOrderPoints = false;
    }

    //! return true if the point id is interprocess of current partition
    bool isInterprocessPoints( index_type pointId ) const
        {
            return this->findInterprocessPoints( pointId ).first;
        }
    //! try to find data of interprocess of current partition point id and return pair(bool,iterator)
    std::pair<bool,typename point_interprocess_map_type::const_iterator> findInterprocessPoints( index_type pointId ) const
        {
            auto itFind = M_interprocessPoints.find( pointId );
            return std::make_pair( itFind != M_interprocessPoints.end(), itFind );
        }
protected:
    //! update interprocess points from mapping ( pt id -> ( isOnActiveElt, isOnGhostEltRanks ) )
    void updateInterprocessPoints( std::unordered_map<index_type,std::tuple<bool,std::set<rank_type>>> const& pointsInterprocessDetection )
        {
            M_interprocessPoints.clear();
            for ( auto const& [pointId,ipData] : pointsInterprocessDetection )
            {
                if ( !std::get<0>( ipData ) ) // not on current process
                    continue;
                if ( std::get<1>( ipData ).empty() ) // not on neighbor process
                    continue;
                //M_interprocessPoints.try_emplace( pointId, std::move( std::get<1>( ipData ) ) );
                M_interprocessPoints.try_emplace( pointId, std::get<1>( ipData ) );
            }
        }

  private:
    void bindSparseMetadata()
    {
        for ( auto& [id, point] : M_points )
        {
            (void)id;
            point.attachSparseMetadataStorage( M_sparseMetadata );
        }
    }

    void buildOrderedPoints()
    {
        M_orderedPoints.clear();
        auto it = beginPoint(), en = endPoint();
        size_type nPoint = std::distance( it, en );
        M_orderedPoints.reserve( nPoint );
        for ( ; it != en; ++it )
            M_orderedPoints.push_back( boost::ref( it->second ) );
        M_needToOrderPoints = true;
        this->updateOrderedPoints();
    }

    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        if constexpr ( Archive::is_loading::value )
        {
            M_points.clear();
            M_sparseMetadata.clear();
            M_coordinateStorage.clear();
            M_orderedPoints.clear();
            M_needToOrderPoints = false;
            size_type nPoints = 0;
            ar& BOOST_SERIALIZATION_NVP( nPoints );
            point_type newPoint;
            for ( size_type k = 0; k < nPoints; ++k )
            {
                ar& boost::serialization::make_nvp( "point", newPoint );
                this->addPoint( std::move( newPoint ) );
            }
        }
        else
        {
            auto it = beginOrderedPoint(), en = endOrderedPoint();
            size_type nPoints = std::distance( it, en );
            ar& BOOST_SERIALIZATION_NVP( nPoints );
            for ( ; it != en; ++it )
                ar& boost::serialization::make_nvp( "point", unwrap_ref( *it ) );
        }
    }

  private:
    worldcomm_ptr_t M_worldComm;

    points_type M_points;
    sparse_metadata_storage_type M_sparseMetadata;
    coordinate_storage_type M_coordinateStorage;
    ordered_points_reference_wrapper_type M_orderedPoints;
    bool M_needToOrderPoints;
    point_interprocess_map_type M_interprocessPoints;
};
/// \endcond
} // namespace Feel
#endif /* __points_H */
