/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-31

  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

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
/**
   \file lagrangep1.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-01-31
 */
#ifndef FEELPP_DISCR_OPERATORLAGRANGEP1_H
#define FEELPP_DISCR_OPERATORLAGRANGEP1_H 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/pointsettomesh.hpp>
//#include <feel/feelfilters/exporterquick.hpp>
#ifdef FEELPP_HAS_GMSH
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#endif
//#include <feel/feelfilters/savegmshmesh.hpp>

namespace Feel
{
namespace detailOpLagP1
{
/**
 * \internal
 * class that allows to compute the basis function type
 */
template<typename SpaceType>
struct SpaceToLagrangeP1Space
{
    typedef SpaceType domain_space_type;
    typedef typename domain_space_type::fe_type::convex_type convex_type;
    typedef typename domain_space_type::mesh_type domain_mesh_type;
    //typedef typename domain_mesh_type::shape_type convex_type;

    typedef typename mpl::if_<mpl::and_<mpl::equal_to<mpl::bool_<convex_type::is_simplex>, mpl::bool_<true> >,
                                       mpl::equal_to<mpl::int_<convex_type::nDim>, mpl::int_<1> > >,
                              mpl::identity<Hypercube<convex_type::nDim,convex_type::nOrder,convex_type::nRealDim > >,
                              mpl::identity<convex_type> >::type::type domain_reference_convex_type;

    template< template<uint16_type Dim> class PsetType>
    struct SelectConvex
    {
        //typedef typename PsetType<convex_type::nDim>::value_type value_type;
        typedef typename SpaceType::value_type value_type;

        typedef Lagrange<1,PsetType> type;
    };
    typedef typename convex_type::template shape<domain_mesh_type::nDim, 1, domain_mesh_type::nRealDim> image_convex_type;
    typedef Mesh<image_convex_type> image_mesh_type;

    typedef typename mpl::if_<mpl::bool_<domain_space_type::is_scalar>,
            mpl::identity<typename SelectConvex<Scalar>::type>,
            typename mpl::if_<mpl::bool_<domain_space_type::is_vectorial>,
            mpl::identity<typename SelectConvex<Vectorial>::type>,
            mpl::identity<typename SelectConvex<Tensor2>::type> >::type>::type::type basis_type;


    typedef FunctionSpace<image_mesh_type, bases<basis_type> > image_space_type;
};

}
/**
 * \class OperatorLagrangeP1
 * \ingroup SpaceTime
 * @author build a P1 Lagrange dof-equivalent space from a \c FunctionSpace
 *
 * \see FunctionSpace
 */
template<typename SpaceType>
class OperatorLagrangeP1
    :
    //public OperatorInterpolation<SpaceType, typename detail::SpaceToLagrangeP1Space<SpaceType>::image_space_type>
public OperatorLinear<SpaceType, typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::image_space_type>
{
    //typedef OperatorInterpolation<SpaceType, typename detail::SpaceToLagrangeP1Space<SpaceType>::image_space_type> super;
    typedef OperatorLinear<SpaceType, typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::image_space_type> super;

public:


    /** @name Typedefs
     */
    //@{
    /**
     * domain space
     */
    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename domain_space_type::mesh_type domain_mesh_type;
    typedef typename domain_space_type::mesh_ptrtype domain_mesh_ptrtype;
    typedef typename domain_space_type::value_type value_type;
    typedef typename domain_space_type::fe_type domain_fe_type;


    typedef typename super::backend_ptrtype backend_ptrtype;
    /**
     * image space
     */
    typedef typename super::dual_image_space_type dual_image_space_type;
    typedef typename super::dual_image_space_ptrtype dual_image_space_ptrtype;
    typedef typename dual_image_space_type::mesh_type image_mesh_type;
    typedef typename dual_image_space_type::mesh_ptrtype image_mesh_ptrtype;
    typedef typename dual_image_space_type::fe_type image_fe_type;

    /**
     * type which instances will hold the correspondence between the
     * mesh elements associated to the domain space and the ones
     * associated to the image space
     */
    typedef std::vector<std::list<size_type> > el2el_type;

    typedef typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::convex_type domain_convex_type;
    typedef typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::domain_reference_convex_type domain_reference_convex_type;
    typedef typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::image_convex_type image_convex_type;
    //typedef PointSetWarpBlend<domain_convex_type, domain_space_type::basis_type::nOrder, value_type> pset_type;
    //typedef PointSetFekete<domain_reference_convex_type/*domain_convex_type*/, domain_space_type::basis_type::nOrder, value_type> pset_type;
    typedef PointSetToMesh<domain_reference_convex_type/*domain_convex_type*/, value_type> p2m_type;
    typedef typename p2m_type::mesh_ptrtype domain_reference_mesh_ptrtype;

    typedef typename domain_mesh_type::gm_type gm_type;
    typedef typename domain_mesh_type::element_type element_type;
    static const size_type gmc_context_v = vm::POINT;
    typedef typename gm_type::template Context<element_type> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gm_type::precompute_ptrtype gmpc_ptrtype;

    //typedef typename mpl::at_c<functionspace_vector_type,SpaceIndex>::type functionspace_ptrtype;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Construct a P1 Lagrange space (ie spanned by a P1 Lagrange
     * basis) out of the space \p space
     */
    OperatorLagrangeP1( domain_space_ptrtype const& space,
                        backend_ptrtype const& backend,
                        std::string pathMeshLagP1=".",
                        std::string prefix="",
                        bool rebuild=true,
                        bool parallelBuild=true,
                        size_type meshUpdate=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK
            );

    /**
     * destructor. nothing really to be done here
     */
    ~OperatorLagrangeP1() override
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * copy operator
     */
    OperatorLagrangeP1 const&  operator=( OperatorLagrangeP1 const& lp1 )
    {
        if ( this != &lp1 )
        {
            M_el2el = lp1.M_el2el;
            M_el2pt = lp1.M_el2pt;
            M_gmpc = lp1.M_gmpc;
            M_p2m = lp1.M_p2m;
        }

        return *this;
    }


    //@}

    /** @name Accessors
     */
    //@{

    uint16_type
    localDof( uint16_type localptid ) const
    {
        CHECK( M_relationRefMeshPointIdToRefFePointId.find(localptid) != M_relationRefMeshPointIdToRefFePointId.end() ) << "point id " << localptid << " not found in relation map";
        return  M_relationRefMeshPointIdToRefFePointId.find(localptid)->second;
    }

    uint16_type
    localDof( typename domain_mesh_type::element_type const& elt,
              uint16_type localptid ) const
    {
        uint16_type localIdPtSet = this->localDof( localptid );

#if !defined(NDEBUG)
        node_type ptRealMesh( domain_fe_type::nRealDim );
        ptRealMesh = M_gmc->xReal( localIdPtSet );

        double dofPtCompareTol = std::max(1e-15,elt.hMin()*1e-5);
        uint16_type localptid_dof = invalid_uint16_type_value;
        bool findPtRelation=false;
        for ( uint16_type i = 0; i < domain_fe_type::nLocalDof && !findPtRelation ; ++i )
        {
            size_type searchdofid = this->domainSpace()->dof()->localToGlobal( elt.id(),i, 0 ).index();
            auto const& ptDof = boost::get<0>( this->domainSpace()->dof()->dofPoint( searchdofid ) );

            bool isSamePoints=true;
            for ( uint16_type d=0;d<domain_fe_type::nRealDim && isSamePoints;++d )
                isSamePoints = isSamePoints && (std::abs( ptDof[d]-ptRealMesh[d] )<dofPtCompareTol);

            if ( isSamePoints )
            {
                localptid_dof =i;
                findPtRelation = true;
            }
        }
        CHECK( localptid_dof != invalid_uint16_type_value ) << "point relation (real dof / real mesh) does not found\n";
        //std::cout << "localptid_dof " << localptid_dof << " vs " << localIdPtSet << "\n";
        CHECK( localptid_dof == localIdPtSet ) << "fast access wrong\n";
        //return localptid_dof;
#endif
        return localIdPtSet;
    }

    image_mesh_ptrtype
    mesh()
    {
        return M_mesh;
    }

    domain_reference_mesh_ptrtype //const&
    referenceMesh() const
    {
        return M_p2m.mesh();
    }
    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * check that the P1 mesh is properly constructed
     */
    void check() const;


    //@}

private:
    OperatorLagrangeP1();
    OperatorLagrangeP1( OperatorLagrangeP1 const & );

    void buildReferenceMesh( bool rebuild, std::string pathMeshLagP1, std::string prefix );
    void buildLagrangeP1Mesh( bool parallelBuild, size_type meshUpdate );
    void buildOperator();

private:
    image_mesh_ptrtype M_mesh;
    el2el_type M_el2el;
    std::vector<std::vector<size_type> > M_el2pt;

    gmpc_ptrtype M_gmpc;
    gmc_ptrtype M_gmc;

    mutable p2m_type M_p2m;

    std::map<uint16_type,uint16_type> M_relationRefMeshPointIdToRefFePointId;
};


//
// OperatorLagrangeP1
//
template<typename space_type>
OperatorLagrangeP1<space_type>::OperatorLagrangeP1( domain_space_ptrtype const& space,
                                                    backend_ptrtype const& backend,
                                                    std::string pathMeshLagP1,
                                                    std::string prefix,
                                                    bool rebuild,
                                                    bool parallelBuild,
                                                    size_type meshUpdate )
    :
    super( space,
           dual_image_space_ptrtype( dual_image_space_type::New( _mesh=image_mesh_ptrtype( new image_mesh_type ) ) ),
           backend,
           false ),
    M_mesh( new image_mesh_type(space->worldCommPtr()) ),
    M_el2el(),
    M_el2pt(),
    M_gmpc( new typename gm_type::precompute_type( this->domainSpace()->mesh()->gm(),
                                                   space->fe()->points()/*M_pset.points()*/ ) ),
    M_p2m()
{

    this->buildReferenceMesh( rebuild, pathMeshLagP1, prefix );

    this->buildLagrangeP1Mesh( parallelBuild, meshUpdate );

    // not work at this time!
    if (false)
        this->buildOperator();

    this->check();

    //
    // Update generates the matrix associated with the interpolation operator
    // comment out for now
#if 0
    this->update();
#endif // 0

}


template<typename space_type>
void
OperatorLagrangeP1<space_type>::buildReferenceMesh( bool rebuild, std::string pathMeshLagP1, std::string prefix )
{
    std::string tagMeshLagP1base = this->domainSpace()->basisName() + (boost::format("-order%1%")%domain_fe_type::nOrder).str();
    if constexpr ( is_lagrange_polynomialset_v<domain_fe_type> )
        {
            if ( domain_convex_type::is_simplex && domain_fe_type::nOrder > 2 )
            {
                if ( domain_fe_type::dual_space_type::template is_pointset_v<PointSetEquiSpaced> )
                    tagMeshLagP1base += "-equispaced";
                else if ( domain_fe_type::dual_space_type::template is_pointset_v<PointSetWarpBlend> )
                    tagMeshLagP1base += "-warpblend";
                else if ( domain_fe_type::dual_space_type::template is_pointset_v<PointSetFeketeSimplex> ) //PointSetFeketeSimplex // PointSetFekete
                    tagMeshLagP1base += "-fekete";
            }
        }

    std::string fileNameMeshDataBase = (boost::format("%1%-%2%d-%3%.msh")%domain_reference_convex_type::type() %domain_reference_convex_type::topological_dimension %tagMeshLagP1base).str();
    // WARNING, if we use a feelpp installed lib, datadir is wrong and can be not exist
    std::string pathMeshDataBase = Environment::expand("$datadir/operatorlagrangep1/"+fileNameMeshDataBase);

    VLOG(1) << "search reference mesh in : "<< pathMeshDataBase;
    bool useMeshInDataBase = false;

    if ( fs::exists( pathMeshDataBase ) )
    {
        VLOG(1) << "[P1 Lagrange] load reference mesh in database";
        useMeshInDataBase = true;
        loadGMSHMesh(_mesh=M_p2m.mesh(),
                     _filename=pathMeshDataBase,
                     _straighten=0,_refine=0,
                     _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES,
                     _worldcomm=M_p2m.mesh()->worldCommPtr(),
                     _respect_partition=0,_rebuild_partitions=0 );
        //if (M_mesh->worldComm().isMasterRank() )
        //    saveGMSHMesh(_mesh=M_p2m.mesh(),_filename="toto.msh");
    }
    else
    {
        VLOG(1) << "[P1 Lagrange] use VTK triangulation load";

        std::string nameMeshLagP1base = "meshLagP1-" + tagMeshLagP1base;
        std::string nameMeshLagP1 = prefixvm(prefix,nameMeshLagP1base);
        fs::path fullPathMesh = fs::path(pathMeshLagP1)/fs::path(nameMeshLagP1+"-1.0.fdb");

        bool masterRankMustLoadMesh = false;
        if ( M_mesh->worldComm().isMasterRank() )
        {
            if ( rebuild || !fs::exists( fullPathMesh ) )
            {
                VLOG(1) << "build lagP1\n";
                // create a triangulation of the current convex
                // using equispaced points defined on the reference
                // element
                //M_p2m.addBoundaryPoints( M_pset.points() );
                PointSet<domain_reference_convex_type,double> thepset( this->domainSpace()->fe()->points() );
                M_p2m.visit( &thepset );

                // #if !defined( NDEBUG )
                //     ExporterQuick<image_mesh_type> exp( "vtk", "ensight" );
                //     exp.save( M_p2m.mesh() );
                // #endif

                // do not renumber the mesh entities
                M_p2m.mesh()->components().clear ( MESH_RENUMBER );
                M_p2m.mesh()->components().clear ( MESH_CHECK );
                M_p2m.mesh()->updateForUse();

                M_p2m.mesh()->save( _name=nameMeshLagP1,_path=pathMeshLagP1 );
            }
            else
                masterRankMustLoadMesh = true;
        }
        // waiting master process write the meh on disk
        M_mesh->worldComm().barrier();

        // load mesh from the disk file
        if ( !M_mesh->worldComm().isMasterRank() || masterRankMustLoadMesh )
        {
            VLOG(1) << "load lagP1\n";
            M_p2m.mesh()->load( _name=nameMeshLagP1,_path=pathMeshLagP1,
                                _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
        }
        //if (M_mesh->worldComm().isMasterRank() )
        //    saveGMSHMesh(_mesh=M_p2m.mesh(),_filename=(fs::path(pathMeshLagP1)/fs::path(nameMeshLagP1+".msh")).string());
    }

    VLOG(2) << "[P1 Lagrange] load mesh done Pointset";

    double tolPtsDiff = (useMeshInDataBase)?1e-9:1e-4;
    // compute relation (ref mesh / ref dof)
    auto itl = M_p2m.mesh()->beginPoint();
    auto const enl = M_p2m.mesh()->endPoint();
    VLOG(2) << "rank "<< Environment::worldComm().rank() << " nPoint " << M_p2m.mesh()->numPoints() << " nElts " << M_p2m.mesh()->numGlobalElements() << "\n";

    for ( ; itl!=enl ; ++itl )
    {
        auto const& ptRefMesh = itl->second;
        size_type __npts = M_gmpc->nPoints();
        uint16_type localIdGeoPc = invalid_uint16_type_value;

        for ( uint16_type __i = 0; __i < __npts; ++__i )
        {
            node_type ptRefFe(domain_fe_type::nRealDim);
            ptRefFe = M_gmpc->node( __i );
            bool isSamePoints=true;
            for (uint16_type d=0;d<domain_fe_type::nRealDim && isSamePoints;++d)
                isSamePoints = isSamePoints && (std::abs( ptRefFe[d]-ptRefMesh[d] )<tolPtsDiff);
            if ( isSamePoints )
            {
                localIdGeoPc = __i;
                break;
            }
        }
        CHECK( localIdGeoPc != invalid_uint16_type_value ) << "point relation (ref mesh / ref dof) does not found\n";
        M_relationRefMeshPointIdToRefFePointId[(uint16_type)ptRefMesh.id()] = localIdGeoPc;
    }
}

template<typename space_type>
void
OperatorLagrangeP1<space_type>::buildLagrangeP1Mesh( bool parallelBuild, size_type meshUpdate )
{
    typedef typename image_mesh_type::point_type point_type;
    typedef typename image_mesh_type::element_type element_type;
    typedef typename image_mesh_type::face_type face_type;

    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, this->domainSpace()->mesh()->markerNames() )
    {
        M_mesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    const rank_type nProc = M_mesh->worldComm().localSize();
    const rank_type procId =  M_mesh->worldComm().localRank();
    bool doParallelBuild = ( nProc==1 )? false : parallelBuild;

    auto meshDomain = this->domainSpace()->mesh();
    auto dofDomain = this->domainSpace()->dof();

    // iterate over active element on process
    auto rangeElements = meshDomain->elementsWithProcessId();
    auto it = std::get<0>( rangeElements );
    auto en = std::get<1>( rangeElements );

    // memory nodes
    std::vector<bool> dofDone( this->domainSpace()->nLocalDof(), false );
    size_type nNewElem  = 0, nNewFaces = 0;


    // data used for parallelism
    // procId -> ( eltId, parentEltId, (pointId_1, pointIdInRef_1) ... )
    std::map<rank_type, std::vector<boost::tuple<size_type,size_type,std::vector<boost::tuple<size_type,uint16_type>>> > > dataToSend;
    std::map<rank_type, std::vector<boost::tuple<size_type,size_type,std::vector<boost::tuple<size_type,uint16_type>>> > > dataToRecv;

    boost::tuple<size_type,size_type,std::vector<boost::tuple<size_type,uint16_type>>> tmpGhostData;
    boost::get<2>( tmpGhostData ).resize( image_mesh_type::element_type::numVertices );

    for ( ; it != en ; ++it )
    {
        auto const& curelt = unwrap_ref( *it );
        DVLOG(2) << "=========================================\n";
        DVLOG(2) << "global element " << curelt.id() << " oriented ok ? : " << curelt.isAnticlockwiseOriented() << "\n";
        DVLOG(2) << "global element G=" << curelt.G() << "\n";

        // update the gemap context
        if ( !M_gmc )
            M_gmc = meshDomain->gm()->template context<gmc_context_v>( curelt, M_gmpc );
        else
            M_gmc->template update<gmc_context_v>( curelt );

        // iterate on each new element (from ref elt)
        auto itl = M_p2m.mesh()->beginElement();
        auto const enl = M_p2m.mesh()->endElement();
        for ( ; itl != enl; ++itl )
        {
            auto const& eltRef = itl->second;
            DVLOG(2) << "************************************\n";
            DVLOG(2) << "local element " << eltRef.id() << " oriented ok ? : " << eltRef.isAnticlockwiseOriented() << "\n";
            DVLOG(2) << "local element G=" << eltRef.G() << "\n";

            element_type elt;
            elt.setId ( nNewElem++ );
            elt.setMarkers( curelt.markers() );
            elt.setProcessIdInPartition( curelt.pidInPartition() );
            elt.setProcessId( curelt.processId() );
            if ( doParallelBuild )
                elt.setNeighborPartitionIds( curelt.neighborPartitionIds() );

            // accumulate the points
            for ( int p = 0; p < image_mesh_type::element_type::numVertices; ++p )
            {
                DVLOG(2) << "local In original element, vertex number " << eltRef.point( p ).id() << "\n";
                uint16_type localptid = eltRef.point( p ).id();
                uint16_type localptid_dof = this->localDof( curelt, localptid );
                size_type gpdof = dofDomain->localToGlobal( curelt.id(),localptid_dof, 0 ).index();
                size_type newPtId = dofDomain->mapGlobalProcessToGlobalCluster( gpdof );

                if ( doParallelBuild )
                    boost::get<2>( tmpGhostData )[p] = boost::make_tuple( newPtId,localptid );

                // add new node if not inserted before
                if ( !dofDone[gpdof] )
                    {
                        dofDone[gpdof] = true;
#if 0
                        point_type __pt( newPtId, boost::get<0>( dofDomain->dofPoint( gpdof ) )  );
#else
                        point_type __pt( newPtId, M_gmc->xReal( localptid_dof ) );
#endif
                        __pt.setProcessId( curelt.processId() );
                        __pt.setProcessIdInPartition( curelt.pidInPartition() );

                        DVLOG(2) << "[OperatorLagrangeP1] element id "
                                 << elt.id() << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] local point id "
                                 << localptid << " coord " << eltRef.point( p ).node() << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] point local id "
                                 << localptid << " global id " << newPtId << "\n"
                                 << " localptid_dof = " << localptid_dof << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] adding point "
                                 << __pt.id() << " : " << __pt.node() << "\n";

                        M_mesh->addPoint( __pt );
                    }

                elt.setPoint( p, M_mesh->point( newPtId ) );
            }

            // add element in mesh
            auto [eit,inserted] = M_mesh->addElement ( elt );
            auto const& [eid,theNewElt] = *eit;

            // store // infos
            if ( doParallelBuild )
            {
                boost::get<0>( tmpGhostData ) = theNewElt.id();
                for ( auto const& [ gpid, idiop ] : curelt.idInOthersPartitions() )
                {
                    boost::get<1>( tmpGhostData ) = idiop;
                    dataToSend[gpid].push_back( tmpGhostData );
                }
            }

        } //  for (  ; itl != enl; ++itl )

    } // for ( size_type elid = 0, pt_image_id = 0; it != en; ++it )


    // add marked faces
    flag_type markerType = 1;
    auto rangeMarkedFaces = meshDomain->facesWithMarkerByType( markerType );
    auto itMarkedFaces = std::get<0>( rangeMarkedFaces );
    auto enMarkedFaces = std::get<1>( rangeMarkedFaces );
    for ( ; itMarkedFaces != enMarkedFaces ; ++itMarkedFaces )
    {
        auto const& curFace = unwrap_ref( *itMarkedFaces );

        std::vector<std::pair<uint16_type,size_type>> dofIdsInElementFromFace( dofDomain->nLocalDofOnFace( true ) );
        for ( int lFaceDofId=0 ; lFaceDofId < dofDomain->nLocalDofOnFace( true ); ++lFaceDofId )
        {
            auto localFaceDof = this->domainSpace()->dof()->localToGlobal( curFace,lFaceDofId,0 );
            size_type gpdof = localFaceDof.index();
            uint16_type ldofInElement = localFaceDof.localDofInElement();
            dofIdsInElementFromFace[lFaceDofId] = std::make_pair(ldofInElement, gpdof);
        }

        auto itfRef = M_p2m.mesh()->beginFace();
        auto const enfRef = M_p2m.mesh()->endFace();
        for ( ; itfRef != enfRef; ++itfRef )
        {
            auto const& faceRef = itfRef->second;
            bool foundFace = true;
            std::vector<size_type> gpdofs( faceRef.nPoints() );
            for ( uint16_type p = 0; p < faceRef.nPoints(); ++p )
            {
                size_type ptIdRef = faceRef.point( p ).id();
                uint16_type ldof = localDof( ptIdRef );
                //ldofs[p] = ldof;
                auto findLocalDof = std::find_if( dofIdsInElementFromFace.begin(), dofIdsInElementFromFace.end(),
                                                  [&ldof]( auto const& i ) { return ldof == i.first; });
                if ( findLocalDof == dofIdsInElementFromFace.end() )
                {
                    foundFace = false;
                    break;
                }
                else
                    gpdofs[p] = findLocalDof->second;
            }
            if ( foundFace )
            {
                face_type newFace;
                newFace.setId( nNewFaces++ );
                //newFace.setProcessId( curelt.processId() );
                newFace.setProcessIdInPartition( curFace.pidInPartition() );
                newFace.setMarker( markerType, curFace.marker( markerType ).value() );
                for ( uint16_type p = 0; p < faceRef.nPoints(); ++p )
                {
                    size_type ptId = dofDomain->mapGlobalProcessToGlobalCluster( gpdofs[ p ] );
                    newFace.setPoint( p, M_mesh->point( ptId ) );
                }
                M_mesh->addFace( newFace );
            }
        }
    }


    if ( doParallelBuild )
    {
        VLOG(2) << "[P1 Lagrange] start build of ghost elements \n";
        int neighborSubdomains = meshDomain->neighborSubdomains().size();
        int nbRequest = 2 * neighborSubdomains;
        mpi::request* reqs = new mpi::request[nbRequest];
        int cptRequest = 0;
        // first send/recv
        for ( rank_type neighborRank : meshDomain->neighborSubdomains() )
        {
            reqs[cptRequest++] = meshDomain->worldComm().localComm().isend( neighborRank, 0, dataToSend[neighborRank] );
            reqs[cptRequest++] = meshDomain->worldComm().localComm().irecv( neighborRank, 0, dataToRecv[neighborRank] );
        }
        // wait all requests
        mpi::wait_all( reqs, reqs + nbRequest );

        // add ghost elements
        for ( auto const& [pid,data] : dataToRecv )
        {
            for ( auto const& dataByElt : data )
            {
                size_type eltId = boost::get<0>( dataByElt );
                size_type parentEltId = boost::get<1>( dataByElt );
                auto const& parentElt = meshDomain->element( parentEltId );

                element_type elt;
                elt.setId( nNewElem++ );
                elt.setMarkers( parentElt.markers() );
                elt.setProcessIdInPartition( procId );
                elt.setProcessId( pid );
                elt.setIdInOtherPartitions( pid, eltId );
                elt.setNeighborPartitionIds( parentElt.neighborPartitionIds() );

                if ( !M_gmc )
                    M_gmc = meshDomain->gm()->template context<gmc_context_v>( parentElt, M_gmpc );
                else
                    M_gmc->template update<gmc_context_v>( parentElt );

                auto const& pointData = boost::get<2>( dataByElt );
                bool isConnectedToActiveAnElement = false;
                for ( uint16_type p = 0; p<pointData.size(); ++p )
                {
                    size_type ptId = boost::get<0>( pointData[p] );
                    uint16_type ptIdInRef = boost::get<1>( pointData[p] );
                    auto ptIterator = M_mesh->pointIterator( ptId );
                    if ( ptIterator == M_mesh->endPoint() )
                    {
                        point_type __pt( ptId, M_gmc->xReal( this->localDof( ptIdInRef ) ) );
                        __pt.setProcessIdInPartition( procId );
                        ptIterator = M_mesh->addPoint( __pt ).first;
                    }
                    else if ( ptIterator->second.processId() != invalid_v<rank_type> )
                        isConnectedToActiveAnElement = true;
                    elt.setPoint( p, ptIterator->second );
                }
                if ( isConnectedToActiveAnElement )
                    M_mesh->addElement( elt );
            }
        }
    }

    DVLOG(2) << "[P1 Lagrange] Number of points in mesh: " << M_mesh->numPoints() << "\n";

    M_mesh->setNumVertices( M_mesh->numPoints() );
    M_mesh->components().reset();
    M_mesh->components().set ( meshUpdate );
    M_mesh->updateForUse();
}


template<typename space_type>
void
OperatorLagrangeP1<space_type>::buildOperator()
{
    // construct the p1 space and set the operator

    //auto Xh_image = dual_image_space_type::New(_mesh=M_mesh,_worldscomm=worldsComm);
    auto Xh_image = dual_image_space_type::New(_mesh=M_mesh,_worldscomm=makeWorldsComm(dual_image_space_type::nSpaces,M_mesh->worldCommPtr() ) );

    this->init( this->domainSpace(),
                Xh_image,
                this->backend(),
                false );
}

template<typename space_type>
void
OperatorLagrangeP1<space_type>::check() const
{
#if !defined(NDEBUG)
    for ( size_type i = 0; i < this->domainSpace()->nLocalDof()/domain_space_type::nComponents; ++i )
    {
        // this test is not good: need to investigate
#if 0
        DLOG_IF( WARNING, (ublas::norm_2( boost::get<0>( this->domainSpace()->dof()->dofPoint( i ) )- M_mesh->point( i ).node() ) > 1e-10 ) )
            << "inconsistent point coordinates at index "   << i
            << "dofpt : " << boost::get<0>( this->domainSpace()->dof()->dofPoint( i ) )
            << "mesh pt : " << M_mesh->point( i ).node() << "\n";
#endif
    }
#endif

}
#if 0

template<typename space_type>
typename OperatorLagrangeP1<space_type>::image_space_type::element_type
OperatorLagrangeP1<space_type>::operator()( element_type const& u ) const
{
    typename image_space_type::element_type res( this->dualImageSpace(), u.name() );
    //res.resize( M_mesh->numPoints() );
    res = ublas::scalar_vector<value_type>( res.size(), value_type( 0 ) );
    // construct the mesh
    typename domain_mesh_type::element_const_iterator it = this->domainSpace()->mesh()->beginElement();
    typename domain_mesh_type::element_const_iterator en = this->domainSpace()->mesh()->endElement();

    for ( ; it != en; ++it )
    {
        this->domainSpace()->mesh()->gm()->template context<gmc_context_v>( unwrap_ref( *it ), M_gmpc );
        typename matrix_node<value_type>::type u_at_pts = u.id( *gmc );

        typename std::list<size_type>::const_iterator ite = M_el2el[it->id()].begin();
        typename std::list<size_type>::const_iterator ene = M_el2el[it->id()].end();

        for ( ; ite != ene; ++ite )
        {
            typename domain_mesh_type::element_type const& elt = M_mesh->element( *ite );

            for ( int c = 0; c < nComponents; ++c )
                for ( int p = 0; p < domain_mesh_type::element_type::numVertices; ++p )
                {
                    size_type ptid = this->dualImageSpace()->dof()->localToGlobal( elt.id(), p, c ).index();
                    res( ptid ) = ublas::column( u_at_pts, M_el2pt[*ite][p] )( nComponents*c+c );
                }
        }
    }

    return res;
}
#endif

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
/**
 * \return the P1 Lagrange adaptor  associated to the space \p Xh
 */
template<typename space_type>
std::shared_ptr<OperatorLagrangeP1<space_type> >
opLagrangeP1_impl( std::shared_ptr<space_type> const& Xh,
                   typename OperatorLagrangeP1<space_type>::backend_ptrtype const& backend,
                   std::string pathMeshLagP1,
                   std::string prefix,
                   bool rebuild,
                   bool parallel,
                   size_type meshUpdate
                   )
{
    return std::shared_ptr<OperatorLagrangeP1<space_type> >( new OperatorLagrangeP1<space_type>( Xh,backend,pathMeshLagP1,prefix,rebuild,parallel,meshUpdate ) );
}

template <typename ... Ts>
auto lagrangeP1( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && space = args.get(_space);
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );
    std::string const& path = args.get_else(_path,".");
    std::string const& prefix = args.get_else(_prefix,"");
    bool rebuild = args.get_else(_rebuild, true );
    bool parallel = args.get_else(_parallel, true );
    size_type update = args.get_else(_update, /*MESH_RENUMBER|*/MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    return opLagrangeP1_impl(space,backend,path,prefix,rebuild,parallel,update);
}


} // Feel
#endif /* __OperatorLagrangeP1_H */
