/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file localization.hpp
 * @brief Physical-to-reference point localization on finite-element meshes.
 */
#ifndef FEELPP_LOCALIZATION_HPP
#define FEELPP_LOCALIZATION_HPP 1

#include <type_traits>

namespace Feel
{

/**
 * @enum LocalizationPolicy
 * @brief Compile-time acceptance strategy for localized points.
 */
enum class LocalizationPolicy
{
    /** Accept every point contained by the candidate element. */
    Unconstrained,
    /** Require supplied geometric points to belong to the candidate. */
    RequireContainedPoints
};

template<typename MeshType>
class Localization
{
public:
    using localization_type = Localization;
    using localization_ptrtype = std::shared_ptr<localization_type>;

    using mesh_type = MeshType;
    using mesh_ptrtype = std::weak_ptr<mesh_type>;
    using range_elements_type = Range<mesh_type,MESH_ELEMENTS>;

    using index_type = typename MeshType::index_type;
    using size_type = typename MeshType::size_type;
    using node_type = typename MeshType::node_type;
    using matrix_node_type = typename matrix_node<typename node_type::value_type>::type;
    
    static inline const uint16_type nDim = mesh_type::nDim;
    static inline const uint16_type nRealDim = mesh_type::nRealDim;
    static inline const uint16_type nOrder = mesh_type::nOrder;
    using value_type = typename mesh_type::value_type;
    using kdtree_type = KDTree;
    using kdtree_ptrtype = std::shared_ptr<KDTree>;

    //!  a node x => a list of id elt which contain the node x
    using node_elem_type = boost::tuple<node_type, std::list<size_type>>;
    
    using container_output_type = std::list<boost::tuple<size_type,node_type>>;
    using container_output_iterator_type = typename container_output_type::iterator;
    //! map between element id and list of node described in the reference elt
    //! typedef std::map<size_type, std::list<boost::tuple<size_type,node_type> > > container_search_type
    using container_search_type = std::map<size_type, container_output_type>;
    using container_search_const_iterator_type = typename container_search_type::const_iterator;
    using container_search_iterator_type = typename container_search_type::iterator;
    
    //!  geomap inverse
    using gm_inverse_type = std::conditional_t<
        mesh_type::element_type::is_simplex,
        GeoMapInverse<nDim,nOrder,nRealDim,value_type,Simplex>,
        GeoMapInverse<nDim,nOrder,nRealDim,value_type,Hypercube>>;
    using gmc_inverse_type = typename gm_inverse_type::gic_type;
    
    using gm1_inverse_type = std::conditional_t<
        mesh_type::element_type::is_simplex,
        GeoMapInverse<nDim,1,nRealDim,value_type,Simplex>,
        GeoMapInverse<nDim,1,nRealDim,value_type,Hypercube>>;
    using gmc1_inverse_type = typename gm1_inverse_type::gic_type;
    
    //!  reference convex
    using ref_convex_type = typename mesh_type::gm_type::reference_convex_type;
    using ref_convex1_type = typename mesh_type::gm1_type::reference_convex_type;

    /**
     * @struct PointLocation
     * @brief Result of mapping one physical point into a specified element.
     */
    struct PointLocation
    {
        /** Whether the point satisfies the selected localization policy. */
        bool found = false;
        /** Reference coordinates produced by the inverse geometric mapping. */
        node_type referencePoint;
        /** Distance returned by the reference-convex containment check. */
        value_type distance = value_type{};
    };

    /**
     * @struct ElementLocation
     * @brief Result of searching the mesh for one physical point.
     */
    struct ElementLocation
    {
        /** Whether a containing or extrapolation element was selected. */
        bool found = false;
        /** Identifier of the selected element. */
        size_type elementId = invalid_v<size_type>;
        /** Reference coordinates in the selected element. */
        node_type referencePoint;
    };

    /**
     * @struct BatchLocation
     * @brief Summary of a batch localization operation.
     *
     * Detailed per-element reference coordinates remain available through
     * `resultAnalysis()` while the compatibility cache is retained.
     */
    struct BatchLocation
    {
        /** Per-input-point localization status. */
        std::vector<bool> found;
        /** Last element selected by the batch-localization walk. */
        size_type lastElementId = invalid_v<size_type>;
    };
    
    //! --------------------------------------------------------------
    //!  Constructors
    //!
    Localization()
        :
        M_mesh (),
        M_kd_tree( new kdtree_type() ),
        M_isInit( false ),
        M_isInitBoundaryFaces( false ),
        M_doExtrapolation( boption( _name=(boost::format("mesh%1%d.localisation.use-extrapolation") % nDim).str() ) ),
        M_barycenter(),
        M_barycentersWorld()
        {
            DVLOG(2) << "[Mesh::Localization] create Localization tool\n";
            int optNbNeighbor = ioption( _name=(boost::format("mesh%1%d.localisation.nelt-in-leaf-kdtree") % nDim).str() );
            int usedNbNeighbor = ( optNbNeighbor < 0 )? 2*mesh_type::element_type::numPoints : optNbNeighbor;
            M_kd_tree->nbNearNeighbor( usedNbNeighbor );

            M_resultAnalysis.clear();
            DVLOG(2) << "[Mesh::Localization] create Localization tool done\n";
        }

    Localization( Localization const & L ) :
            M_mesh( L.M_mesh ),
            M_kd_tree( new kdtree_type( *( L.M_kd_tree ) ) ),
            M_geoGlob_Elts( L.M_geoGlob_Elts ),
            M_isInit( L.M_isInit ),
            M_isInitBoundaryFaces( L.M_isInitBoundaryFaces ),
            M_resultAnalysis( L.M_resultAnalysis ),
            M_doExtrapolation( L.M_doExtrapolation ),
            M_gic( L.M_gic ), M_gic1( L.M_gic1 ),
            M_barycenter( L.M_barycenter ),
            M_barycentersWorld( L.M_barycentersWorld )
        {}

        //! --------------------------------------------------------------
         //!  Define the mesh whith or not init
         //!
        void
        setMesh( std::shared_ptr<mesh_type> m,bool b=true )
        {
            M_mesh = m;
            if ( b )
                this->init( elements(m) );
            else M_isInit=b;

            M_resultAnalysis.clear();
        }

        void
        setMesh( std::shared_ptr<mesh_type> m, range_elements_type const& range, bool b=true )
        {
            M_mesh = m;
            M_rangeElements = range;
            if ( b )
                this->init( elements(m) );
            else
                M_isInit=b;

            M_resultAnalysis.clear();
        }

        //! --------------------------------------------------------------
         //!  Define if necessary to use extrapolation
         //!
        void
        setExtrapolation( bool b )
        {
            M_doExtrapolation = b;
        }

        //! --------------------------------------------------------------
         //!  Run the init function if necessary
         //!
        void updateForUse()
        {
            if ( this->isInit() )
                return;
            auto mesh = M_mesh.lock();
            if ( !mesh )
                return;
            auto rangeElt = M_rangeElements? *M_rangeElements : elements(mesh);
            this->init( rangeElt );
        }

        //! --------------------------------------------------------------
         //!  Run the init function if necessary
         //!
        void updateForUseBoundaryFaces()
        {
            if ( !this->isInitBoundaryFaces() )
                this->initBoundaryFaces();
        }

        //! --------------------------------------------------------------
         //!  Access
         //!
        bool isInit() const
        {
            return M_isInit;
        }

        bool isInitBoundaryFaces() const
        {
            return M_isInitBoundaryFaces;
        }

        bool doExtrapolation() const
        {
            return M_doExtrapolation;
        }

        mesh_ptrtype mesh()
        {
            return M_mesh;
        }
        mesh_ptrtype mesh() const
        {
            return M_mesh;
        }

        kdtree_ptrtype kdtree()
        {
            return M_kd_tree;
        }

        kdtree_ptrtype const& kdtree() const
        {
            return M_kd_tree;
        }

        node_type const& barycenter() const
        {
            CHECK( this->isInit() || this->isInitBoundaryFaces()  ) << " localization tool not init \n";
            return M_barycenter;
        }

        void computeBarycenter();

        bool hasComputedBarycentersWorld()
        {
#if BOOST_VERSION >= 105600
            return M_barycentersWorld != boost::none;
#else
            return M_barycentersWorld;
#endif
        }

        std::vector<boost::tuple<bool,node_type> > const& barycentersWorld() const
        {
            CHECK( M_barycentersWorld ) << " you must call computeBarycentersWorld() before barycentersWorld() \n";
            return M_barycentersWorld.get();
        }

        void computeBarycentersWorld();

        /**
         * @brief Return reference points grouped by containing element.
         *
         * @return compatibility result cache populated by `localize()`
         */
        [[nodiscard]] container_search_type const& resultAnalysis() const
        {
            return M_resultAnalysis;
        }

        /** @brief Compatibility spelling for `resultAnalysis()`. */
        [[nodiscard]] container_search_type const& result_analysis() const
        {
            return this->resultAnalysis();
        }

        /** @brief Return the beginning of the compatibility result cache. */
        container_search_iterator_type resultAnalysisBegin()
        {
            return M_resultAnalysis.begin();
        }

        /** @brief Return the end of the compatibility result cache. */
        container_search_iterator_type resultAnalysisEnd()
        {
            return M_resultAnalysis.end();
        }

        /** @brief Compatibility spelling for `resultAnalysisBegin()`. */
        container_search_iterator_type result_analysis_begin()
        {
            return this->resultAnalysisBegin();
        }

        /** @brief Compatibility spelling for `resultAnalysisEnd()`. */
        container_search_iterator_type result_analysis_end()
        {
            return this->resultAnalysisEnd();
        }

        /**
         * @brief Locate a physical point in a specified mesh element.
         * @param elementId candidate element identifier
         * @param point physical point
         * @return point localization result
         */
        [[nodiscard]] PointLocation
        locateInElement( size_type elementId, node_type const& point ) const
        {
            return this->locateInElementImpl( elementId, point );
        }

        /**
         * @brief Locate a point using a compile-time conformity strategy.
         * @tparam Policy point-set acceptance policy
         * @param elementId candidate element identifier
         * @param point physical point
         * @param requiredPoints points that must belong to a conforming candidate
         * @return point localization result
         */
        template <LocalizationPolicy Policy>
        [[nodiscard]] PointLocation
        locateInElement( size_type elementId,
                         node_type const& point,
                         matrix_node_type const& requiredPoints ) const
        {
            if constexpr ( Policy == LocalizationPolicy::RequireContainedPoints )
                return this->locateInElementConformingImpl(
                    elementId, point, requiredPoints );
            else
                return this->locateInElementImpl( elementId, point );
        }

        /**
         * @brief Search the mesh for an element containing a physical point.
         * @param point physical point
         * @return element localization result
         */
        [[nodiscard]] ElementLocation locateElement( node_type const& point )
        {
            return this->locateElementImpl( point );
        }

        /**
         * @brief Search using a compile-time conformity strategy.
         * @tparam Policy point-set acceptance policy
         * @param point physical point
         * @param requiredPoints points that must belong to a conforming candidate
         * @return element localization result
         */
        template <LocalizationPolicy Policy>
        [[nodiscard]] ElementLocation
        locateElement( node_type const& point,
                       matrix_node_type const& requiredPoints )
        {
            if constexpr ( Policy == LocalizationPolicy::RequireContainedPoints )
                return this->locateElementConformingImpl( point, requiredPoints );
            else
                return this->locateElementImpl( point );
        }

        /**
         * @brief Localize a batch of physical points.
         * @param points physical points stored by column
         * @param hypotheticalElement first candidate element, or `invalid_v`
         * @return batch-localization summary
         */
        BatchLocation
        localize( matrix_node_type const& points,
                  size_type hypotheticalElement )
        {
            return this->localizeImpl( points, hypotheticalElement );
        }

        /**
         * @brief Localize a batch using a compile-time conformity strategy.
         * @tparam Policy point-set acceptance policy
         * @param points physical points stored by column
         * @param hypotheticalElement first candidate element, or `invalid_v`
         * @param requiredPoints points that must belong to a conforming candidate
         * @return batch-localization summary
         */
        template <LocalizationPolicy Policy>
        BatchLocation
        localize( matrix_node_type const& points,
                  size_type hypotheticalElement,
                  matrix_node_type const& requiredPoints )
        {
            if constexpr ( Policy == LocalizationPolicy::RequireContainedPoints )
                return this->localizeConformingImpl(
                    points, hypotheticalElement, requiredPoints );
            else
                return this->localizeImpl( points, hypotheticalElement );
        }

        /**
         * @brief Localize using a standard compile-time Boolean policy.
         * @tparam Conforming whether required points constrain candidates
         * @param points physical points stored by column
         * @param hypotheticalElement first candidate element, or `invalid_v`
         * @param requiredPoints points required by conforming localization
         * @return batch-localization summary
         */
        template <bool Conforming>
        BatchLocation
        localize( matrix_node_type const& points,
                  size_type hypotheticalElement,
                  matrix_node_type const& requiredPoints,
                  std::bool_constant<Conforming> )
        {
            constexpr auto policy =
                Conforming
                    ? LocalizationPolicy::RequireContainedPoints
                    : LocalizationPolicy::Unconstrained;
            return this->template localize<policy>(
                points, hypotheticalElement, requiredPoints );
        }

        /** @brief Legacy tuple adapter for `locateInElement()`. */
        boost::tuple<bool,node_type,double>
        isIn( size_type elementId, node_type const& point ) const
        {
            auto location = this->locateInElement( elementId, point );
            return boost::make_tuple(
                location.found, location.referencePoint, location.distance );
        }

        /** @brief Legacy MPL adapter for conforming `locateInElement()`. */
        boost::tuple<bool,node_type,double>
        isIn( size_type elementId,
              node_type const& point,
              matrix_node_type const& requiredPoints,
              mpl::int_<1> ) const
        {
            auto location = this->template locateInElement<
                LocalizationPolicy::RequireContainedPoints>(
                    elementId, point, requiredPoints );
            return boost::make_tuple(
                location.found, location.referencePoint, location.distance );
        }

        boost::tuple<uint16_type,std::vector<bool> > isIn( std::vector<size_type> _ids, const node_type & _pt );

        /** @brief Legacy tuple adapter for `locateElement()`. */
        boost::tuple<bool, size_type,node_type>
        searchElement( node_type const& point )
        {
            auto location = this->locateElement( point );
            return boost::make_tuple(
                location.found, location.elementId, location.referencePoint );
        }

        /** @brief Legacy MPL adapter for unconstrained `locateElement()`. */
        boost::tuple<bool, size_type,node_type>
        searchElement( node_type const& point,
                       matrix_node_type const&,
                       mpl::int_<0> )
        {
            return this->searchElement( point );
        }

        /** @brief Legacy MPL adapter for conforming `locateElement()`. */
        boost::tuple<bool, size_type,node_type>
        searchElement( node_type const& point,
                       matrix_node_type const& requiredPoints,
                       mpl::int_<1> )
        {
            auto location = this->template locateElement<
                LocalizationPolicy::RequireContainedPoints>(
                    point, requiredPoints );
            return boost::make_tuple(
                location.found, location.elementId, location.referencePoint );
        }

        //! ---------------------------------------------------------------
         //!  Research all elements which contains the node p
         //!
        boost::tuple<bool, std::list<boost::tuple<size_type,node_type> > > searchElements( const node_type & p );

        //! ---------------------------------------------------------------
         //!  Research the element which contains the node p, forall p in the
         //!  matrix_node_type m. The result is save by this object
         //!
        boost::tuple<std::vector<bool>, size_type>
        run_analysis( matrix_node_type const& points,
                      size_type const& hypotheticalElement )
        {
            auto result = this->localize( points, hypotheticalElement );
            return boost::make_tuple(
                std::move( result.found ), result.lastElementId );
        }

        //! ---------------------------------------------------------------
         //!  Research the element which contains the node p, forall p in the
         //!  matrix_node_type m. The result is save by this object
         //!
        boost::tuple<std::vector<bool>, size_type>
        run_analysis( matrix_node_type const& points,
                      size_type const& hypotheticalElement,
                      matrix_node_type const&,
                      mpl::int_<0> )
        {
            return this->run_analysis( points, hypotheticalElement );
        }

        //! ---------------------------------------------------------------
         //!  Research the element which contains the node p, forall p in the
         //!  matrix_node_type m. The result is save by this object
         //!
        boost::tuple<std::vector<bool>, size_type>
        run_analysis( matrix_node_type const& points,
                      size_type const& hypotheticalElement,
                      matrix_node_type const& requiredPoints,
                      mpl::int_<1> )
        {
            auto result = this->template localize<
                LocalizationPolicy::RequireContainedPoints>(
                    points, hypotheticalElement, requiredPoints );
            return boost::make_tuple(
                std::move( result.found ), result.lastElementId );
        }

        /**
         * @brief Localize points with a standard compile-time conformity tag.
         *
         * @tparam Conforming whether candidate vertices must belong to the
         * localized element
         *
         * This C++20-facing overload preserves the established localization
         * implementation while allowing clients to avoid Boost.MPL dispatch
         * tags.
         */
        template <bool Conforming>
        boost::tuple<std::vector<bool>, size_type>
        run_analysis( matrix_node_type const& points,
                      size_type const& hypotheticalElement,
                      matrix_node_type const& requiredPoints,
                      std::bool_constant<Conforming> )
        {
            auto result = this->localize(
                points, hypotheticalElement, requiredPoints,
                std::bool_constant<Conforming>{} );
            return boost::make_tuple(
                std::move( result.found ), result.lastElementId );
        }

        //! ---------------------------------------------------------------
         //!  Reset all data
         //!
        void reset()
        {
            M_isInit=false;
            M_isInitBoundaryFaces=false;

            //clear data
            M_geoGlob_Elts.clear();
            if ( M_kd_tree )
                M_kd_tree->clear();
            //this->updateForUse();
        }

        void resetBoundaryFaces()
        {
            M_isInit=false;
            M_isInitBoundaryFaces=false;
            this->initBoundaryFaces();
        }

    private :

        //! ---------------------------------------------------------------
        //! initializes the kd tree and the map between node and list elements(all elements)
        //!
        FEELPP_NO_EXPORT void init( range_elements_type const& range );

        //! ---------------------------------------------------------------
        //! initializes the kd tree and the map between node and list elements(only on boundary)
        //!
        FEELPP_NO_EXPORT void initBoundaryFaces();

        //! ---------------------------------------------------------------
        //! search near elt in kd tree and get a sorted list
        //!
        FEELPP_NO_EXPORT void searchInKdTree( const node_type & p,
                                              std::list< std::pair<size_type, uint> > & listTri );

        /** @brief Implement unconstrained localization in one element. */
        FEELPP_NO_EXPORT PointLocation
        locateInElementImpl( size_type elementId,
                             node_type const& point ) const;

        /** @brief Implement conforming localization in one element. */
        FEELPP_NO_EXPORT PointLocation
        locateInElementConformingImpl(
            size_type elementId,
            node_type const& point,
            matrix_node_type const& requiredPoints ) const;

        /** @brief Implement unconstrained single-element search. */
        FEELPP_NO_EXPORT ElementLocation
        locateElementImpl( node_type const& point );

        /** @brief Implement conforming single-element search. */
        FEELPP_NO_EXPORT ElementLocation
        locateElementConformingImpl(
            node_type const& point,
            matrix_node_type const& requiredPoints );

        /** @brief Implement unconstrained batch localization. */
        FEELPP_NO_EXPORT BatchLocation
        localizeImpl( matrix_node_type const& points,
                      size_type hypotheticalElement );

        /** @brief Implement conforming batch localization. */
        FEELPP_NO_EXPORT BatchLocation
        localizeConformingImpl(
            matrix_node_type const& points,
            size_type hypotheticalElement,
            matrix_node_type const& requiredPoints );

        /** @brief Compute the KD-tree barycenter in any physical dimension. */
        FEELPP_NO_EXPORT node_type computeBarycenterImpl() const;

    private:

        mesh_ptrtype M_mesh;
        std::optional<range_elements_type> M_rangeElements;
        kdtree_ptrtype M_kd_tree;
        //! map between node and list elements
        std::map<size_type, node_elem_type > M_geoGlob_Elts;
        bool M_isInit,M_isInitBoundaryFaces;
        container_search_type M_resultAnalysis;
        bool M_doExtrapolation;

        ref_convex_type M_refelem;
        ref_convex1_type M_refelem1;
        mutable std::shared_ptr<gmc_inverse_type> M_gic;
        mutable std::shared_ptr<gmc1_inverse_type> M_gic1;

        node_type M_barycenter;
        boost::optional<std::vector<boost::tuple<bool,node_type> > > M_barycentersWorld;
        
    };


template<typename MeshType>
void
Localization<MeshType>::init( range_elements_type const& range )
{
    this->reset();

    auto mesh = M_mesh.lock();
    if ( !mesh ) return;

    VLOG(1) << "Localization apply init";

    for ( auto const& eltWrap : range )
    {
        auto const& elt = unwrap_ref( eltWrap );
        M_gic.reset( new gmc_inverse_type( mesh->gm(), elt, mesh->worldComm().subWorldCommSeqPtr() ) );
        M_gic1.reset( new gmc1_inverse_type(
            mesh->gm1(), elt, geometry::AffineGeometryOrder{},
            mesh->worldComm().subWorldCommSeqPtr() ) );
        break;
    }
    for ( auto const& eltWrap : range )
    {
        auto const& elt = unwrap_ref( eltWrap );
        index_type eltId = elt.id();
        for ( int i=0; i<elt.nPoints(); ++i )
        {
            auto const& thePoint = elt.point( i );
            index_type thePointId = thePoint.id();

            auto itFindPoint = M_geoGlob_Elts.find( thePointId );
            if ( itFindPoint == M_geoGlob_Elts.end() )
            {
                node_elem_type theData = boost::make_tuple( thePoint.node(), std::list<size_type>( {eltId} ) );
                M_geoGlob_Elts.emplace( std::make_pair( thePointId, std::move( theData ) ) );
                M_kd_tree->addPoint( thePoint.node(), thePointId );
            }
            else
                boost::get<1>( itFindPoint->second ).push_back( eltId );
        }
    }

    this->computeBarycenter();

    M_isInit=true;
}

template<typename MeshType>
void
Localization<MeshType>::initBoundaryFaces()
{
    auto mesh = M_mesh.lock();
    if ( !mesh ) return;

    DLOG_IF( WARNING, this->isInitBoundaryFaces() == false ) << "You have already initialized the tool of localization\n";

    //clear data
    M_geoGlob_Elts.clear();
    M_kd_tree->clear()
;
    // typename mesh_type::location_face_iterator face_it;
    // typename mesh_type::location_face_iterator face_en;
    // boost::tie( boost::tuples::ignore, face_it, face_en ) = Feel::boundaryfaces( mesh );
    auto rangeBoundaryFaces = Feel::boundaryfaces( mesh );
    auto face_it = rangeBoundaryFaces.begin();
    auto face_en = rangeBoundaryFaces.end();
    bool hasInitGic=false;
    for ( ; face_it != face_en; ++face_it )
    {
        auto const& face = boost::unwrap_ref( *face_it );
        for ( int i=0; i<face.nPoints(); ++i )
        {
            if ( face.isConnectedTo0() )
                {
                    if ( boost::get<1>( M_geoGlob_Elts[face.point( i ).id()] ).size()==0 )
                        {
                            boost::get<0>( M_geoGlob_Elts[face.point( i ).id()] ) = face.point( i ).node();
                            M_kd_tree->addPoint( face.point( i ).node(),face.point( i ).id() );
                        }
                    boost::get<1>( M_geoGlob_Elts[face.point( i ).id()] ).push_back( face.element( 0 ).id() );

                    if ( !hasInitGic )
                    {
                        M_gic.reset( new gmc_inverse_type( mesh->gm(), face.element( 0 ), mesh->worldComm().subWorldCommSeqPtr() ) );
                        M_gic1.reset( new gmc1_inverse_type(
                            mesh->gm1(), face.element( 0 ),
                            geometry::AffineGeometryOrder{},
                            mesh->worldComm().subWorldCommSeqPtr() ) );
                        hasInitGic=true;
                    }
                }
        }
    }

    this->computeBarycenter();

    M_isInitBoundaryFaces=true;
    M_isInit=false;

}



template<typename MeshType>
typename Localization<MeshType>::PointLocation
Localization<MeshType>::locateInElementImpl(
    size_type elementId, node_type const& point ) const
{
    bool isin=false;
    double dmin=0;
    node_type x_ref;

    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt = mesh->element( elementId );
    auto const eltView = mesh->elementView( elementId );

    if ( elt.isOnBoundary() || nDim != nRealDim )
        {
#if 0
            // get inverse geometric transformation
            gmc_inverse_type gic( mesh->gm(), elt, mesh->worldComm().subWorldCommSeqPtr() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( point);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem.isIn( gic.xRef() );
#else
            M_gic->update( eltView );
            M_gic->setXReal( point);
            x_ref=M_gic->xRef();
            if ( nDim == nRealDim )
                isin = M_gic->isIn();
            else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                boost::tie( isin, dmin ) = M_refelem.isIn( M_gic->xRef() );
#endif
        }
    else
        {
#if 0
            // get inverse geometric transformation
            gmc1_inverse_type gic(
                mesh->gm1(), elt, geometry::AffineGeometryOrder{},
                mesh->worldComm().subWorldCommSeqPtr() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( point);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem1.isIn( gic.xRef() );
#else
            M_gic1->update( eltView, geometry::AffineGeometryOrder{} );
            M_gic1->setXReal( point);
            x_ref=M_gic1->xRef();
            if ( nDim == nRealDim )
                isin = M_gic1->isIn();
            else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                boost::tie( isin, dmin ) = M_refelem1.isIn( M_gic1->xRef() );
#endif
        }

    return PointLocation{ isin, std::move( x_ref ), dmin };
}

template<typename MeshType>
boost::tuple<uint16_type,std::vector<bool> >
Localization<MeshType>::isIn( std::vector<size_type> _ids, const node_type & _pt )
{
    typedef typename mesh_type::gm_type::reference_convex_type ref_convex_type;
    typedef typename mesh_type::gm1_type::reference_convex_type ref_convex1_type;

    uint16_type nbId = _ids.size();
    std::vector<bool> isin( _ids.size(),false );
    bool isin2=false;
    double dmin;
    node_type __x_ref;

    auto mesh = M_mesh.lock();
    uint16_type nbIsIn=0;

    for ( uint16_type i = 0; i< nbId ; ++i )
    {
        //get element with the id
        auto const& elt = mesh->element( _ids[i] );
        auto const eltView = mesh->elementView( _ids[i] );

        if ( elt.isOnBoundary() || nDim != nRealDim )
            {
#if 0
                // get inverse geometric transformation
                gmc_inverse_type gic( mesh->gm(), elt );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem.isIn( gic.xRef() );
#else
                M_gic->update( eltView );
                M_gic->setXReal( _pt);
                __x_ref=M_gic->xRef();
                if ( nDim == nRealDim )
                    isin2 = M_gic->isIn();
                else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                    boost::tie( isin2, dmin ) = M_refelem.isIn( M_gic->xRef() );
#endif
                isin[i] = isin2;
            }
        else
            {
#if 0
                // get inverse geometric transformation
                gmc1_inverse_type gic(
                    mesh->gm1(), elt, geometry::AffineGeometryOrder{} );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem1.isIn( gic.xRef() );
#else
                M_gic1->update( eltView, geometry::AffineGeometryOrder{} );
                M_gic1->setXReal( _pt);
                __x_ref=M_gic1->xRef();
                if ( nDim == nRealDim )
                    isin2 = M_gic1->isIn();
                else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                    boost::tie( isin2, dmin ) = M_refelem1.isIn( M_gic1->xRef() );
#endif
                isin[i] = isin2;
            }
        if (isin[i]) ++nbIsIn;
    }

    return boost::make_tuple( nbIsIn,isin );
}


template<typename MeshType>
typename Localization<MeshType>::ElementLocation
Localization<MeshType>::locateElementImpl( node_type const& point )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    auto mesh = M_mesh.lock();
    bool isin=false;double dmin=0;
    node_type x_ref;
    size_type idEltFound = invalid_v<size_type>;

    std::list< std::pair<size_type, uint> > ListTri;
    this->searchInKdTree(point,ListTri);

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto const itLT_end=ListTri.end();

    DCHECK( std::distance( itLT,itLT_end )>0 ) << "problem in localization : listTri is empty\n";

    while ( itLT != itLT_end && !isin  )
    {
        //get element with the id
        //elt = M_mesh->element( itLT->first );

        // search point in this elt
        auto location = this->locateInElement( itLT->first, point );
        isin = location.found;
        x_ref = std::move( location.referencePoint );
        dmin = location.distance;
        // if not inside, continue the research with an other element
        if (!isin) ++itLT;
        else idEltFound=itLT->first;
    }


    if( !isin && this->doExtrapolation() )
        {
            // first elt
            idEltFound = ListTri.begin()->first;
            auto const& eltUsedForExtrapolation = mesh->element(idEltFound);
            DVLOG(1) << "localisation tool use extrapolation for the point" << point
                     << " with elt.id() " << eltUsedForExtrapolation.id()
                     << " and elt.G() " << eltUsedForExtrapolation.G()
                     << "\n";
            gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeqPtr() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( point);
            isin = true;
            idEltFound = eltUsedForExtrapolation.id();
            x_ref = gic.xRef();
        }

    return ElementLocation{ isin, idEltFound, std::move( x_ref ) };

}

template<typename MeshType>
typename Localization<MeshType>::BatchLocation
Localization<MeshType>::localizeImpl(
    matrix_node_type const& points,
    size_type hypotheticalElement )
{
    // if no init then init with all geo point
    if ( !( this->isInit() || this->isInitBoundaryFaces() ) )
        this->updateForUse();

    bool find_x;
    size_type cv_id=hypotheticalElement;
    node_type x_ref;double dmin=0;
    std::vector<bool> hasFindPts(points.size2(),false);

    M_resultAnalysis.clear();

    bool doExtrapolationAtStart = this->doExtrapolation();
    auto nPtMaxNearNeighborAtStart = this->kdtree()->nPtMaxNearNeighbor();


    // first step : no extrapolation
    if ( doExtrapolationAtStart ) this->setExtrapolation( false );

    // first currentEltHypothetical
    auto currentEltHypothetical = hypotheticalElement;
    for ( size_type i=0; i< points.size2(); ++i )
        {
            bool testHypothetical_find = false;

            if ( hypotheticalElement!=invalid_v<size_type> )
                {
                    auto location = this->locateInElement(
                        currentEltHypothetical, ublas::column( points, i ) );
                    testHypothetical_find = location.found;
                    x_ref = std::move( location.referencePoint );
                    dmin = location.distance;
                }
            if ( testHypothetical_find )
                {
                    cv_id = currentEltHypothetical;
                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                    hasFindPts[i]=true;
                    DVLOG(1) << "localisation tool (hypthetical elt) : has found pt " << ublas::column( points, i )
                             << "in elt id " << cv_id
                             << "with G() " << M_mesh.lock()->element(cv_id).G()
                             << "\n";
                }
            else // search kdtree
                {
                    // kdtree call
                    auto location = this->locateElement(
                        ublas::column( points, i ) );
                    find_x = location.found;
                    cv_id = location.elementId;
                    x_ref = std::move( location.referencePoint );
                    // traitement
                    if (find_x) // if find : OK
                        {
                            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                            currentEltHypothetical = cv_id;
                            hasFindPts[i]=true;
                            DVLOG(1)  << "localisation tool (first pass) : has found pt " << ublas::column( points, i )
                                      << " in elt id " << cv_id
                                      << " with G() " << M_mesh.lock()->element(cv_id).G()
                                      << " and ref point " << x_ref
                                      << "\n";
                        }
                    else// if (false) // try an other method (no efficient but maybe a solution)
                        {
                            // search in all element
                            this->kdtree()->nbNearNeighbor(2*nPtMaxNearNeighborAtStart /*15*/ /*this->mesh()->numElements()*/);

                            location = this->locateElement(
                                ublas::column( points, i ) );
                            find_x = location.found;
                            cv_id = location.elementId;
                            x_ref = std::move( location.referencePoint );
                            //revert parameter
                            this->kdtree()->nbNearNeighbor(nPtMaxNearNeighborAtStart);
                            // if find : OK (but strange!)
                            if (find_x)
                                {
                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;
                                    DVLOG(1)  << "localisation tool (second pass) : has found pt " << ublas::column( points, i )
                                              << " in elt id " << cv_id
                                              << " with G() " << M_mesh.lock()->element(cv_id).G()
                                              << " and ref point " << x_ref
                                              << "\n";
                                }
                            else if (doExtrapolationAtStart)
                                {
                                    this->setExtrapolation(true);

                                    location = this->locateElement(
                                        ublas::column( points, i ) );
                                    find_x = location.found;
                                    cv_id = location.elementId;
                                    x_ref = std::move( location.referencePoint );

                                    CHECK( find_x ) << "localisation tool : invalid extrapolation \n";

                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;

                                    this->setExtrapolation(false);
                                }
                        }
                } // search kdtree

            DLOG_IF(WARNING, !hasFindPts[i]) << "the localisation tool fails to find the point " << ublas::column( points, i ) << "\n";

        } // for (size_type i=0;i< m.size2();++i)

    //revert parameter
    this->setExtrapolation(doExtrapolationAtStart);

    return BatchLocation{ std::move( hasFindPts ), cv_id };

} // localizeImpl



template<typename MeshType>
boost::tuple<bool, std::list<boost::tuple<typename Localization<MeshType>::size_type, typename MeshType::node_type> > >
Localization<MeshType>::searchElements( const node_type & p )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    //this->kdtree()->nbNearNeighbor(this->mesh()->numElements());

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( p,ListTri );

    typename mesh_type::element_type elt;
    typename mesh_type::gm_type::reference_convex_type refelem;
    typename mesh_type::gm1_type::reference_convex_type refelem1;

    bool isin=false;
    double dmin;
    node_type x_ref;

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    CHECK( std::distance( itLT,itLT_end ) > 0 )
        << "problem in list localization: empty candidate list";
#endif

    std::list<boost::tuple<size_type,node_type> > newlistelts;
    newlistelts.clear();
    bool find = false;
    bool finishSearch=false;
    while ( itLT != itLT_end && !finishSearch /*&& !isin*/  )
    {
#if 0
        //get element with the id
        elt= M_mesh->element( itLT->first );

        if ( elt.isOnBoundary() )
        {
            // get inverse geometric transformation
            typename mesh_type::Inverse::gic_type gic( M_mesh->gm(), elt );

            //apply the inverse geometric transformation for the point p
            gic.setXReal( p );
            __x_ref=gic.xRef();

            // the point is in the reference element ?
            boost::tie( isin, dmin ) = refelem.isIn( gic.xRef() );
        }

        else
        {
            // get inverse geometric transformation
            typename mesh_type::Inverse::gic1_type gic(
                M_mesh->gm1(), elt, geometry::AffineGeometryOrder{} );

            //apply the inverse geometric transformation for the point p
            gic.setXReal( p );
            __x_ref=gic.xRef();

            // the point is in the reference element ?
            boost::tie( isin, dmin ) = refelem1.isIn( gic.xRef() );
            //std::cout << "gic.xRef()" << gic.xRef() << std::endl;
        }
#else
        boost::tie(isin,x_ref, dmin) = this->isIn(itLT->first,p);
#endif
        if ( isin )
        {
            newlistelts.push_back( boost::make_tuple( itLT->first,x_ref ) );
            find = true;
            // if not on boundary -> finish for this point
            if (dmin>1e-7) finishSearch=true;
        }

        //if (find) std::cout << elt.G() << std::endl;

        //if not inside, continue the research with an other element
        //if (!isin) ++itLT;
        ++itLT;
    }

    if ( !find ) std::cout << "\n WARNING EXTRAPOLATION IN SEARCHELEMENTS!!!"<<std::endl;

    if ( find )
        return boost::make_tuple( true,newlistelts );

    else if ( !find && !M_doExtrapolation )
        return boost::make_tuple( false,newlistelts );

    else
    {
        auto mesh = M_mesh.lock();
        //std::cout << "\n WARNING EXTRAPOLATION \n";
        itLT=ListTri.begin();
        elt= mesh->element( itLT->first );
        typename mesh_type::Inverse::gic_type gic( mesh->gm(), elt );
        //apply the inverse geometric transformation for the point p
        //gic.setXReal(boost::get<0>(*ptsNN.begin()));
        gic.setXReal( p );
        x_ref=gic.xRef();
        //return boost::make_tuple( true, itLT->first, __x_ref);
        newlistelts.push_back( boost::make_tuple( itLT->first,x_ref ) );
        find = true;
        return boost::make_tuple( true,newlistelts );
    }
} // searchElements


template<typename MeshType>
void
Localization<MeshType>::searchInKdTree( const node_type & p,
        std::list< std::pair<size_type, uint> > & ListTri )
{
    //search for nearest points
    M_kd_tree->search( p );

    //get the results of research
    typename KDTree::points_search_type ptsNN = M_kd_tree->pointsNearNeighbor();

    typename KDTree::points_search_const_iterator itNN = ptsNN.begin();
    typename KDTree::points_search_const_iterator itNN_end = ptsNN.end();

#if !defined( NDEBUG )
    CHECK( std::distance( itNN, itNN_end ) > 0 )
        << "no near-neighbor point found";
#endif

    //iterator on a l(ist index element
    typename std::list<size_type>::iterator itL;
    typename std::list<size_type>::iterator itL_end;

    //ListTri will contain the indices of elements (size_type)
    //and the number of occurrence(uint)
    //std::list< std::pair<size_type, uint> > ListTri;
    typename std::list< std::pair<size_type, uint> >::iterator itLT;
    typename std::list< std::pair<size_type, uint> >:: iterator itLT_end;

    //create of ListTri : sort largest to smallest occurrences
    //In case of equality : if the point is closer than another then it will be before
    //                      if it is in the same point then  the lowest index will be before
    for ( ; itNN != itNN_end; ++itNN )
    {
        itL= boost::get<1>( M_geoGlob_Elts[boost::get<3>( *itNN )] ).begin();
        itL_end= boost::get<1>( M_geoGlob_Elts[boost::get<3>( *itNN )] ).end();

        for ( ; itL != itL_end; ++itL )
        {
            itLT=ListTri.begin();
            itLT_end=ListTri.end();
            bool find=false;

            while ( itLT != itLT_end && !find )
            {
                if ( itLT->first == *itL ) find=true;

                else ++itLT;
            }

            if ( find )
            {
                uint nb=itLT->second+1;
                size_type numEl=itLT->first;
                ListTri.remove( *itLT );
                itLT=ListTri.begin();
                itLT_end=ListTri.end();
                bool find=false;

                while ( itLT != itLT_end && !find )
                {
                    if ( itLT->second < nb ) find=true;

                    else ++itLT;
                }

                ListTri.insert( itLT,std::make_pair( numEl,nb ) );
            }

            else ListTri.push_back( std::make_pair( *itL,1 ) );
        }
    }

}

template<typename MeshType>
typename Localization<MeshType>::PointLocation
Localization<MeshType>::locateInElementConformingImpl(
    size_type elementId,
    node_type const& point,
    matrix_node_type const& requiredPoints ) const
{
    constexpr value_type pointTolerance = static_cast<value_type>( 1e-8 );
    auto mesh = M_mesh.lock();
    auto const& elt= mesh->element( elementId );
    auto const& eltG = elt.G();

    for ( size_type i = 0; i < requiredPoints.size2(); ++i )
    {
        auto const requiredPoint = ublas::column( requiredPoints, i );
        bool pointFound = false;
        for ( size_type j = 0; j < eltG.size2() && !pointFound; ++j )
        {
            auto const candidatePoint = ublas::column( eltG, j );
            if ( requiredPoint.size() != candidatePoint.size() )
                continue;

            pointFound = true;
            for ( size_type component = 0;
                  component < requiredPoint.size();
                  ++component )
            {
                if ( std::abs( requiredPoint( component ) -
                              candidatePoint( component ) ) >= pointTolerance )
                {
                    pointFound = false;
                    break;
                }
            }
        }

        if ( !pointFound )
            return PointLocation{};
    }

    auto location = this->locateInElement( elementId, point );
    LOG_IF( ERROR, !location.found )
        << "Localization conforming point-set check succeeded but inverse "
           "geometric localization failed";
    return location;
}



template<typename MeshType>
typename Localization<MeshType>::ElementLocation
Localization<MeshType>::locateElementConformingImpl(
    node_type const& point,
    matrix_node_type const& requiredPoints )
{
    bool isin=false;
    node_type x_ref;
    auto mesh=this->mesh().lock();
    size_type idEltFound = invalid_v<size_type>;

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( point,ListTri );

    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    CHECK( std::distance( itLT,itLT_end ) > 0 )
        << "problem in list localization: empty candidate list";
#endif


    //research the element which contains the point p
    while ( itLT != itLT_end && !isin  )
    {
        auto location = this->template locateInElement<
            LocalizationPolicy::RequireContainedPoints>(
                itLT->first, point, requiredPoints );
        isin = location.found;
        x_ref = std::move( location.referencePoint );

        //if not inside, continue the research with an other element
        if ( !isin ) ++itLT;
        else idEltFound=itLT->first;
    } //while ( itLT != itLT_end && !isin  )

    if (!isin)
        {
            if( this->doExtrapolation() )
                {
                    LOG(WARNING) << "WARNING EXTRAPOLATION for the point" << point;
                    //std::cout << "W";
                    auto const& eltUsedForExtrapolation = mesh->element(ListTri.begin()->first);
                    gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeqPtr() );
                    //apply the inverse geometric transformation for the point p
                    gic.setXReal( point);
                    isin = true;
                    idEltFound = eltUsedForExtrapolation.id();
                    x_ref = gic.xRef();
                }
            else
                {
                    // idEltFound = invalid_v<size_type>;//mesh->beginElementWithId(mesh->worldComm().localRank())->id();
                    isin = false;
                    //x_ref=?
                }
        }

    return ElementLocation{ isin, idEltFound, std::move( x_ref ) };

} // locateElementConformingImpl


template<typename MeshType>
typename Localization<MeshType>::BatchLocation
Localization<MeshType>::localizeConformingImpl(
    matrix_node_type const& points,
    size_type hypotheticalElement,
    matrix_node_type const& requiredPoints )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    bool find_x=false;
    size_type cv_id=hypotheticalElement;
    node_type x_ref;
    std::vector<bool> hasFindPts(points.size2(),false);

    M_resultAnalysis.clear();
    auto currentEltHypothetical = hypotheticalElement;
    for ( size_type i=0; i< points.size2(); ++i )
    {

        bool testHypothetical_find = false;
        if ( hypotheticalElement!=invalid_v<size_type> )
        {
            auto location = this->template locateInElement<
                LocalizationPolicy::RequireContainedPoints>(
                    currentEltHypothetical,
                    ublas::column( points, i ),
                    requiredPoints );
            testHypothetical_find = location.found;
            x_ref = std::move( location.referencePoint );
        }
        if ( testHypothetical_find )
        {
            cv_id = currentEltHypothetical;
            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
            hasFindPts[i]=true;
        }
        else
        {
            auto location = this->template locateElement<
                LocalizationPolicy::RequireContainedPoints>(
                    ublas::column( points, i ), requiredPoints );
            find_x = location.found;
            cv_id = location.elementId;
            x_ref = std::move( location.referencePoint );

            if ( find_x )
            {
                M_resultAnalysis[cv_id].push_back( boost::make_tuple( i,x_ref ) );
                hasFindPts[i]=true;
                currentEltHypothetical = cv_id;
            }
        }
        //else std::cout<<"\nNew Probleme Localization\n" << std::endl;
    }

    return BatchLocation{ std::move( hasFindPts ), cv_id };

} // localizeConformingImpl


template<typename MeshType>
void
Localization<MeshType>::computeBarycenter()
{
    M_barycenter = this->computeBarycenterImpl();
}

template<typename MeshType>
typename MeshType::node_type
Localization<MeshType>::computeBarycenterImpl() const
{
    node_type result( nRealDim );
    result.clear();

    auto const pointCount = this->kdtree()->nPoints();
    if ( pointCount == 0 )
        return result;

    for ( size_type pointIndex = 0; pointIndex < pointCount; ++pointIndex )
    {
        auto const& point =
            this->kdtree()->points()[pointIndex].template get<0>();
        for ( uint16_type component = 0; component < nRealDim; ++component )
            result( component ) += point( component );
    }

    for ( uint16_type component = 0; component < nRealDim; ++component )
        result( component ) /= pointCount;
    return result;
}


template<typename MeshType>
void
Localization<MeshType>::computeBarycentersWorld()
{
    LOG(INFO) << "computeBarycentersWorld : run mpi::all_gather\n";
    auto mesh = M_mesh.lock();
    M_barycentersWorld = std::vector<boost::tuple<bool,node_type> >( mesh->worldComm().localSize() );
    mpi::all_gather( mesh->worldComm().localComm(),
                     boost::make_tuple( this->kdtree()->nPoints() > 0 , this->barycenter() ),
                     *M_barycentersWorld );
}


} // Feel
#endif
