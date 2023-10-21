//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 22 Oct 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_LOCALIZATION_HPP
#define FEELPP_LOCALIZATION_HPP 1



namespace Feel {

template<typename MeshType>
class Localization
{
public:    
    typedef Localization localization_type;
    typedef std::shared_ptr<localization_type> localization_ptrtype;

    
        

    using mesh_type = MeshType;
    typedef std::weak_ptr<mesh_type> mesh_ptrtype;
    using range_elements_type = Range<mesh_type,MESH_ELEMENTS>;

    using index_type = typename MeshType::index_type;
    using size_type = typename MeshType::size_type;
    using node_type = typename MeshType::node_type;
    typedef typename matrix_node<typename node_type::value_type>::type matrix_node_type;
    
    static inline const uint16_type nDim = mesh_type::nDim;
    static inline const uint16_type nRealDim = mesh_type::nRealDim;
    static inline const uint16_type nOrder = mesh_type::nOrder;
    using value_type = typename mesh_type::value_type;
    typedef KDTree kdtree_type;
    typedef typename std::shared_ptr<KDTree> kdtree_ptrtype;

    //!  a node x => a list of id elt which contain the node x
    typedef boost::tuple<node_type, std::list<size_type> > node_elem_type;
    
    typedef typename std::list<boost::tuple<size_type,node_type> > container_output_type;
    typedef typename std::list<boost::tuple<size_type,node_type> >::iterator container_output_iterator_type;
    //! map between element id and list of node described in the reference elt
    //! typedef std::map<size_type, std::list<boost::tuple<size_type,node_type> > > container_search_type
    typedef std::map<size_type, container_output_type > container_search_type;
    typedef typename container_search_type::const_iterator container_search_const_iterator_type;
    typedef typename container_search_type::iterator container_search_iterator_type;
    
    //!  geomap inverse
    typedef typename mpl::if_<mpl::bool_<mesh_type::element_type::is_simplex>,
                              mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,value_type,Simplex> >,
                              mpl::identity<GeoMapInverse<nDim,nOrder,nRealDim,value_type,Hypercube> > >::type::type gm_inverse_type;
    typedef typename gm_inverse_type::gic_type gmc_inverse_type;
    
    typedef typename mpl::if_<mpl::bool_<mesh_type::element_type::is_simplex>,
                              mpl::identity<GeoMapInverse<nDim,1,nRealDim,value_type,Simplex> >,
                              mpl::identity<GeoMapInverse<nDim,1,nRealDim,value_type,Hypercube> > >::type::type gm1_inverse_type;
    typedef typename gm1_inverse_type::gic_type gmc1_inverse_type;
    
    //!  reference convex
    typedef typename mesh_type::gm_type::reference_convex_type ref_convex_type;
    typedef typename mesh_type::gm1_type::reference_convex_type ref_convex1_type;
    
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

        container_search_type const & result_analysis() const { return M_resultAnalysis;}

        container_search_iterator_type result_analysis_begin()
        {
            return M_resultAnalysis.begin();
        }
        container_search_iterator_type result_analysis_end()
        {
            return M_resultAnalysis.end();
        }

        //! ---------------------------------------------------------------
         //!  True if the node p is in mesh->element(id)
         //!
        boost::tuple<bool,node_type,double> isIn( size_type _id, const node_type & _pt ) const;
        boost::tuple<bool,node_type,double> isIn( size_type _id, const node_type & _pt, const matrix_node_type & setPoints, mpl::int_<1> /**/ ) const;
        boost::tuple<uint16_type,std::vector<bool> > isIn( std::vector<size_type> _ids, const node_type & _pt );

        //! ---------------------------------------------------------------
         //!  Research only one element which contains the node p
         //!
        boost::tuple<bool, size_type,node_type> searchElement(const node_type & p);
        //! ---------------------------------------------------------------
         //!  Research only one element which contains the node p
         //!
        boost::tuple<bool, size_type,node_type> searchElement(const node_type & p,
                                                              const matrix_node_type & setPoints,
                                                              mpl::int_<0> /**/ )
        {
            return searchElement( p );
        }

        //! ---------------------------------------------------------------
         //!  Research only one element which contains the node p and which this elt have as geometric point contain setPoints
         //!
        boost::tuple<bool, size_type,node_type> searchElement( const node_type & p,
                                                               const matrix_node_type & setPoints,
                                                               mpl::int_<1> /**/ );

        //! ---------------------------------------------------------------
         //!  Research all elements which contains the node p
         //!
        boost::tuple<bool, std::list<boost::tuple<size_type,node_type> > > searchElements( const node_type & p );

        //! ---------------------------------------------------------------
         //!  Research the element which contains the node p, forall p in the
         //!  matrix_node_type m. The result is save by this object
         //!
        boost::tuple<std::vector<bool>, size_type> run_analysis(const matrix_node_type & m,
                                                                const size_type & eltHypothetical);

        //! ---------------------------------------------------------------
         //!  Research the element which contains the node p, forall p in the
         //!  matrix_node_type m. The result is save by this object
         //!
        boost::tuple<std::vector<bool>, size_type>  run_analysis(const matrix_node_type & m,
                                                                 const size_type & eltHypothetical,
                                                                 const matrix_node_type & setPoints,
                                                                 mpl::int_<0> /**/)
        {
            return run_analysis( m,eltHypothetical );
        }

        //! ---------------------------------------------------------------
         //!  Research the element which contains the node p, forall p in the
         //!  matrix_node_type m. The result is save by this object
         //!
        boost::tuple<std::vector<bool>, size_type> run_analysis(const matrix_node_type & m,
                                                                const size_type & eltHypothetical,
                                                                const matrix_node_type & setPoints,
                                                                mpl::int_<1> /**/);

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

        //! ---------------------------------------------------------------
        //!  computed barycenter
        //!
        FEELPP_NO_EXPORT node_type computeBarycenter(mpl::int_<1> /**/) const;
        FEELPP_NO_EXPORT node_type computeBarycenter(mpl::int_<2> /**/) const;
        FEELPP_NO_EXPORT node_type computeBarycenter(mpl::int_<3> /**/) const;

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
        M_gic1.reset( new gmc1_inverse_type( mesh->gm1(), elt, mpl::int_<1>(), mesh->worldComm().subWorldCommSeqPtr() ) );
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
    auto face_it = boost::get<1>( rangeBoundaryFaces );
    auto face_en = boost::get<2>( rangeBoundaryFaces );
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
                        M_gic1.reset( new gmc1_inverse_type( mesh->gm1(), face.element( 0 ), mpl::int_<1>(), mesh->worldComm().subWorldCommSeqPtr() ) );
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
boost::tuple<bool,typename MeshType::node_type,double>
Localization<MeshType>::isIn( size_type _id, const node_type & _pt ) const
{
    bool isin=false;
    double dmin=0;
    node_type x_ref;

    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt = mesh->element( _id );

    if ( elt.isOnBoundary() || nDim != nRealDim )
        {
#if 0
            // get inverse geometric transformation
            gmc_inverse_type gic( mesh->gm(), elt, mesh->worldComm().subWorldCommSeqPtr() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem.isIn( gic.xRef() );
#else
            M_gic->update( elt );
            M_gic->setXReal( _pt);
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
            gmc1_inverse_type gic( mesh->gm1(), elt, mpl::int_<1>(), mesh->worldComm().subWorldCommSeqPtr() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( _pt);
            x_ref=gic.xRef();
            // the point is in the reference element ?
            boost::tie( isin, dmin ) = M_refelem1.isIn( gic.xRef() );
#else
            M_gic1->update( elt, mpl::int_<1>() );
            M_gic1->setXReal( _pt);
            x_ref=M_gic1->xRef();
            if ( nDim == nRealDim )
                isin = M_gic1->isIn();
            else // in this case, result given with gic->isIn() seems not work (see geomap.hpp)
                boost::tie( isin, dmin ) = M_refelem1.isIn( M_gic1->xRef() );
#endif
        }

    return boost::make_tuple(isin,x_ref,dmin);
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
                M_gic->update( elt );
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
                gmc1_inverse_type gic( mesh->gm1(), elt, mpl::int_<1>() );
                //apply the inverse geometric transformation for the point p
                gic.setXReal( _pt);
                __x_ref=gic.xRef();
                // the point is in the reference element ?
                boost::tie( isin2, dmin ) = M_refelem1.isIn( gic.xRef() );
#else
                M_gic1->update( elt, mpl::int_<1>() );
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
boost::tuple<bool, typename Localization<MeshType>::size_type, typename MeshType::node_type>
Localization<MeshType>::searchElement( const node_type & p )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    auto mesh = M_mesh.lock();
    bool isin=false;double dmin=0;
    node_type x_ref;
    size_type idEltFound = 0;//mesh->beginElementWithId(mesh->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    this->searchInKdTree(p,ListTri);

    //research the element which contains the point p
    auto itLT=ListTri.begin();
    auto const itLT_end=ListTri.end();

    DCHECK( std::distance( itLT,itLT_end )>0 ) << "problem in localization : listTri is empty\n";

    while ( itLT != itLT_end && !isin  )
    {
        //get element with the id
        //elt = M_mesh->element( itLT->first );

        // search point in this elt
        boost::tie(isin,x_ref,dmin) = this->isIn(itLT->first,p);
        // if not inside, continue the research with an other element
        if (!isin) ++itLT;
        else idEltFound=itLT->first;
    }


    if( !isin && this->doExtrapolation() )
        {
            // first elt
            idEltFound = ListTri.begin()->first;
            auto const& eltUsedForExtrapolation = mesh->element(idEltFound);
            DVLOG(1) << "localisation tool use extrapolation for the point" << p
                     << " with elt.id() " << eltUsedForExtrapolation.id()
                     << " and elt.G() " << eltUsedForExtrapolation.G()
                     << "\n";
            gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeqPtr() );
            //apply the inverse geometric transformation for the point p
            gic.setXReal( p);
            boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

}

template<typename MeshType>
boost::tuple<std::vector<bool>, typename Localization<MeshType>::size_type>
Localization<MeshType>::run_analysis( const matrix_node_type & m,
                                      const size_type & eltHypothetical )
{
    // if no init then init with all geo point
    if ( !( this->isInit() || this->isInitBoundaryFaces() ) )
        this->updateForUse();

    bool find_x;
    size_type cv_id=eltHypothetical;
    node_type x_ref;double dmin=0;
    std::vector<bool> hasFindPts(m.size2(),false);

    M_resultAnalysis.clear();

    bool doExtrapolationAtStart = this->doExtrapolation();
    auto nPtMaxNearNeighborAtStart = this->kdtree()->nPtMaxNearNeighbor();


    // first step : no extrapolation
    if ( doExtrapolationAtStart ) this->setExtrapolation( false );

    // first currentEltHypothetical
    auto currentEltHypothetical = eltHypothetical;//this->mesh()->beginElement()->id();//eltHypothetical;
    for ( size_type i=0; i< m.size2(); ++i )
        {
            bool testHypothetical_find = false;

            if ( eltHypothetical!=invalid_v<size_type> )
                {
                    boost::tie( testHypothetical_find,x_ref,dmin ) = this->isIn( currentEltHypothetical,ublas::column( m, i ) );
                }
            if ( testHypothetical_find )
                {
                    cv_id = currentEltHypothetical;
                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                    hasFindPts[i]=true;
                    DVLOG(1) << "localisation tool (hypthetical elt) : has found pt " << ublas::column( m, i )
                             << "in elt id " << cv_id
                             << "with G() " << M_mesh.lock()->element(cv_id).G()
                             << "\n";
                }
            else // search kdtree
                {
                    // kdtree call
                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                    // traitement
                    if (find_x) // if find : OK
                        {
                            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                            currentEltHypothetical = cv_id;
                            hasFindPts[i]=true;
                            DVLOG(1)  << "localisation tool (first pass) : has found pt " << ublas::column( m, i )
                                      << " in elt id " << cv_id
                                      << " with G() " << M_mesh.lock()->element(cv_id).G()
                                      << " and ref point " << x_ref
                                      << "\n";
                        }
                    else// if (false) // try an other method (no efficient but maybe a solution)
                        {
                            // search in all element
                            this->kdtree()->nbNearNeighbor(2*nPtMaxNearNeighborAtStart /*15*/ /*this->mesh()->numElements()*/);

                            boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));
                            //revert parameter
                            this->kdtree()->nbNearNeighbor(nPtMaxNearNeighborAtStart);
                            // if find : OK (but strange!)
                            if (find_x)
                                {
                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;
                                    DVLOG(1)  << "localisation tool (second pass) : has found pt " << ublas::column( m, i )
                                              << " in elt id " << cv_id
                                              << " with G() " << M_mesh.lock()->element(cv_id).G()
                                              << " and ref point " << x_ref
                                              << "\n";
                                }
                            else if (doExtrapolationAtStart)
                                {
                                    this->setExtrapolation(true);

                                    boost::tie( find_x, cv_id, x_ref ) = this->searchElement(ublas::column( m, i ));

                                    CHECK( find_x ) << "localisation tool : invalid extrapolation \n";

                                    M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
                                    currentEltHypothetical = cv_id;
                                    hasFindPts[i]=true;

                                    this->setExtrapolation(false);
                                }
                        }
                } // search kdtree

            DLOG_IF(WARNING, !hasFindPts[i]) << "the localisation tool fails to find the point " << ublas::column( m, i ) << "\n";

        } // for (size_type i=0;i< m.size2();++i)

    //revert parameter
    this->setExtrapolation(doExtrapolationAtStart);

    return boost::make_tuple(hasFindPts,cv_id);

} //run_analysis



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
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
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
            typename mesh_type::Inverse::gic1_type gic( M_mesh->gm1(), elt,mpl::int_<1>() );

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
    FEELPP_ASSERT( std::distance( itNN,itNN_end )>0 ).error( "none Near Neighbor Points are find" );
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
boost::tuple<bool,typename MeshType::node_type,double>
Localization<MeshType>::isIn( size_type _id,
                                         const node_type & _pt,
                                         const matrix_node_type & setPoints,
                                         mpl::int_<1> /**/ ) const
{
    bool isin=true; // warning : start with true
    double dmin=0.;
    node_type x_ref;
    auto mesh = M_mesh.lock();
    //get element with the id
    auto const& elt= mesh->element( _id );
    auto const& eltG = elt.G();

    // check conformity between setPoints (given) and eltG (localize)
    std::vector<bool> find( setPoints.size2(),false );
    for ( size_type i=0; i< setPoints.size2() && isin; ++i )
    {
        auto const& thePt = ublas::column( setPoints,i );

        for ( size_type j=0; j<eltG.size2(); ++j )
        {
            auto ptjeltG = ublas::column( eltG,j );

            if ( ptjeltG.size()==1 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 )
                    find[i]=true;
            }
            else if ( ptjeltG.size()==2 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 &&
                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-8 )
                    find[i]=true;
            }
            else if ( ptjeltG.size()==3 )
            {
                if ( std::abs( thePt( 0 )-ptjeltG( 0 ) )<1e-8 &&
                     std::abs( thePt( 1 )-ptjeltG( 1 ) )<1e-8 &&
                     std::abs( thePt( 2 )-ptjeltG( 2 ) )<1e-8 )
                    find[i]=true;
            }
        }
        // up checking
        isin &= find[i];
    }

    // if find -> get ref point and check
    if ( isin )
    {
        bool isin2=false;
        boost::tie(isin2,x_ref,dmin) = this->isIn(_id,_pt);
        LOG_IF(ERROR, !isin2) << "Mesh::Localization::isIn<Conformal> : check fail -> maybe x_ref is not correct";
    }

    return boost::make_tuple(isin,x_ref,dmin);
}



template<typename MeshType>
boost::tuple<bool, typename Localization<MeshType>::size_type, typename MeshType::node_type>
Localization<MeshType>::searchElement( const node_type & p,
        const matrix_node_type & setPoints,
        mpl::int_<1> /**/ )
{
    typename mesh_type::element_type elt;
    typename mesh_type::gm_type::reference_convex_type refelem;
    typename mesh_type::gm1_type::reference_convex_type refelem1;

    bool isin=false,isin2=false;
    double dmin=0.;
    node_type x_ref;
    auto mesh=this->mesh().lock();
    size_type idEltFound = invalid_v<size_type>;//mesh->beginElementWithId(mesh->worldComm().localRank())->id();

    std::list< std::pair<size_type, uint> > ListTri;
    searchInKdTree( p,ListTri );

    auto itLT=ListTri.begin();
    auto itLT_end=ListTri.end();

#if !defined( NDEBUG )
    //if(std::distance(itLT,itLT_end)==0) std::cout<<"\nListTri vide\n";
    FEELPP_ASSERT( std::distance( itLT,itLT_end )>0 ).error( " problem in list localization : is empty" );
#endif


    //research the element which contains the point p
    while ( itLT != itLT_end && !isin  )
    {
        boost::tie(isin,x_ref,dmin) = this->isIn(itLT->first,p,setPoints,mpl::int_<1>());

        //if not inside, continue the research with an other element
        if ( !isin ) ++itLT;
        else idEltFound=itLT->first;
    } //while ( itLT != itLT_end && !isin  )

    if (!isin)
        {
            if( this->doExtrapolation() )
                {
                    LOG(WARNING) << "WARNING EXTRAPOLATION for the point" << p;
                    //std::cout << "W";
                    auto const& eltUsedForExtrapolation = mesh->element(ListTri.begin()->first);
                    gmc_inverse_type gic( mesh->gm(), eltUsedForExtrapolation, mesh->worldComm().subWorldCommSeqPtr() );
                    //apply the inverse geometric transformation for the point p
                    gic.setXReal( p);
                    boost::tie(isin,idEltFound,x_ref) = boost::make_tuple(true,eltUsedForExtrapolation.id(),gic.xRef());
                }
            else
                {
                    // idEltFound = invalid_v<size_type>;//mesh->beginElementWithId(mesh->worldComm().localRank())->id();
                    isin = false;
                    //x_ref=?
                }
        }

    return boost::make_tuple( isin, idEltFound, x_ref);

} //searchElement


template<typename MeshType>
boost::tuple<std::vector<bool>, typename Localization<MeshType>::size_type>
Localization<MeshType>::run_analysis( const matrix_node_type & m,
        const size_type & eltHypothetical,
        const matrix_node_type & setPoints,
        mpl::int_<1> /**/ )
{

    DCHECK( this->isInit() || this->isInitBoundaryFaces() ) << "You don't have initialized the tool of localization\n";

    bool find_x=false;
    size_type cv_id=eltHypothetical;
    node_type x_ref;
    double dmin;
    std::vector<bool> hasFindPts(m.size2(),false);

    M_resultAnalysis.clear();
    auto currentEltHypothetical = eltHypothetical;
    for ( size_type i=0; i< m.size2(); ++i )
    {

        bool testHypothetical_find = false;
        if ( eltHypothetical!=invalid_v<size_type> )
        {
            boost::tie( testHypothetical_find,x_ref,dmin ) = this->isIn( currentEltHypothetical,ublas::column( m, i ),setPoints,mpl::int_<1>() );
        }
        if ( testHypothetical_find )
        {
            cv_id = currentEltHypothetical;
            M_resultAnalysis[cv_id].push_back( boost::make_tuple(i,x_ref) );
            hasFindPts[i]=true;
        }
        else
        {
            boost::tie( find_x, cv_id, x_ref ) = this->searchElement( ublas::column( m, i ),setPoints,mpl::int_<1>() );

            if ( find_x )
            {
                M_resultAnalysis[cv_id].push_back( boost::make_tuple( i,x_ref ) );
                hasFindPts[i]=true;
                currentEltHypothetical = cv_id;
            }
        }
        //else std::cout<<"\nNew Probleme Localization\n" << std::endl;
    }

    return boost::make_tuple(hasFindPts,cv_id);

} // run_analysis


template<typename MeshType>
void
Localization<MeshType>::computeBarycenter()
{
    M_barycenter = this->computeBarycenter(mpl::int_<nRealDim>());
}

template<typename MeshType>
typename MeshType::node_type
Localization<MeshType>::computeBarycenter(mpl::int_<1> /**/) const
{
    node_type res(1);
    res(0)=0;
    if ( this->kdtree()->nPoints()==0 ) return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);
        }
    res(0)/=this->kdtree()->nPoints();
    return res;
}
template<typename MeshType>
typename MeshType::node_type
Localization<MeshType>::computeBarycenter(mpl::int_<2> /**/) const
{
    node_type res(2);
    res(0)=0;res(1)=0;
    if ( this->kdtree()->nPoints()==0 )
        return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();
    return res;
}
template<typename MeshType>
typename MeshType::node_type
Localization<MeshType>::computeBarycenter(mpl::int_<3> /**/) const
{
    node_type res(3);
    res(0)=0;res(1)=0;res(2)=0;
    if ( this->kdtree()->nPoints()==0 )
        return res;
    for (size_type i = 0 ; i<this->kdtree()->nPoints() ; ++i)
        {
            auto const& pt = this->kdtree()->points()[i].template get<0>();
            res(0)+=pt(0);res(1)+=pt(1);res(2)+=pt(2);
        }
    res(0)/=this->kdtree()->nPoints();res(1)/=this->kdtree()->nPoints();res(2)/=this->kdtree()->nPoints();
    return res;
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
