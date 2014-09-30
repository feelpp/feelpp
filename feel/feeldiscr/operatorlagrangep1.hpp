/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
#ifndef __OperatorLagrangeP1_H
#define __OperatorLagrangeP1_H 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/pointsettomesh.hpp>
//#include <feel/feelfilters/exporterquick.hpp>

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
    typedef typename convex_type::template shape<domain_mesh_type::nDim, 1, domain_mesh_type::nRealDim>::type image_convex_type;
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
     * type which instances will hold the correspondance between the
     * mesh elements associated to the domain space and the ones
     * associated to the image space
     */
    typedef std::vector<std::list<size_type> > el2el_type;

    typedef typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::convex_type domain_convex_type;
    typedef typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::domain_reference_convex_type domain_reference_convex_type;
    typedef typename detailOpLagP1::SpaceToLagrangeP1Space<SpaceType>::image_convex_type image_convex_type;
    //typedef PointSetWarpBlend<domain_convex_type, domain_space_type::basis_type::nOrder, value_type> pset_type;
    typedef PointSetFekete<domain_reference_convex_type/*domain_convex_type*/, domain_space_type::basis_type::nOrder, value_type> pset_type;
    typedef PointSetToMesh<domain_reference_convex_type/*domain_convex_type*/, value_type> p2m_type;
    typedef typename p2m_type::mesh_ptrtype domain_reference_mesh_ptrtype;

    typedef typename domain_mesh_type::gm_type gm_type;
    typedef typename domain_mesh_type::element_type element_type;
    typedef typename gm_type::template Context<vm::POINT, element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
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
                        bool parallelBuild=true );

    /**
     * destructor. nothing really to be done here
     */
    ~OperatorLagrangeP1()
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
            M_pset = lp1.M_pset;
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
    localDof( typename domain_mesh_type::element_const_iterator it,
              uint16_type localptid ) const
    {
        //return localDof( it, localptid, mpl::bool_<is_modal>() );
        return localDof( it, localptid, mpl::bool_<false>() );
    }

    uint16_type
    localDof( typename domain_mesh_type::element_const_iterator /*it*/,
              uint16_type localptid,
              mpl::bool_<false> ) const
    {
        return localptid;
    }

    uint16_type
    localDof( typename domain_mesh_type::element_const_iterator it,
              uint16_type localptid,
              mpl::bool_<true> ) const;


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
    void buildLagrangeP1Mesh( bool parallelBuild );
    void buildOperator();

private:
    image_mesh_ptrtype M_mesh;
    el2el_type M_el2el;
    std::vector<std::vector<size_type> > M_el2pt;

    pset_type M_pset;
    gmpc_ptrtype M_gmpc;

    mutable p2m_type M_p2m;

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
                                                    bool parallelBuild )
    :
    super( space,
           dual_image_space_ptrtype( new dual_image_space_type( image_mesh_ptrtype( new image_mesh_type ) ) ),
           backend,
           false ),
    M_mesh( new image_mesh_type(space->worldComm()) ),
    M_el2el(),
    M_el2pt(),
    M_pset( 0 ),
    M_gmpc( new typename gm_type::precompute_type( this->domainSpace()->mesh()->gm(),
             M_pset.points() ) ),
    M_p2m()
{

    this->buildReferenceMesh( rebuild, pathMeshLagP1, prefix );

    this->buildLagrangeP1Mesh( parallelBuild );

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


#if 0

    for ( size_type i = 0; i < this->domainSpace()->nLocalDof()/domain_space_type::nComponents; ++i )
    {
        FEELPP_ASSERT( boost::get<1>( this->domainSpace()->dof()->dofPoint( i ) ) ==
                       boost::get<1>( this->dualImageSpace()->dof()->dofPoint( i ) ) )
        ( boost::get<0>( this->domainSpace()->dof()->dofPoint( i ) ) )
        ( boost::get<0>( this->dualImageSpace()->dof()->dofPoint( i ) ) )
        ( i ).warn( "check inconsistent point id" );
        FEELPP_ASSERT( ublas::norm_2( boost::get<0>( this->domainSpace()->dof()->dofPoint( i ) )-
                                      boost::get<0>( this->dualImageSpace()->dof()->dofPoint( i ) )
                                    ) < 1e-10 )
        ( boost::get<0>( this->domainSpace()->dof()->dofPoint( i ) ) )
        ( boost::get<0>( this->dualImageSpace()->dof()->dofPoint( i ) ) )
        ( i ).warn( "check inconsistent point coordinates" );
    }

#endif

}


template<typename space_type>
void
OperatorLagrangeP1<space_type>::buildReferenceMesh( bool rebuild, std::string pathMeshLagP1, std::string prefix )
{
    std::string nameMeshLagP1base = "meshLagP1-"+this->domainSpace()->basisName()+(boost::format("-order%1%")%domain_space_type::basis_type::nOrder).str();
    std::string nameMeshLagP1 = prefixvm(prefix,nameMeshLagP1base);
    if ( rebuild )
    {
        if (M_mesh->worldComm().isMasterRank() )
        {
            // create a triangulation of the current convex
            // using equispaced points defined on the reference
            // element
            //M_p2m.addBoundaryPoints( M_pset.points() );
            M_p2m.visit( &M_pset );

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

        M_mesh->worldComm().barrier();

        if (M_mesh->worldComm().globalRank() != M_mesh->worldComm().masterRank() )
        {
            M_p2m.mesh()->load( _name=nameMeshLagP1,_path=pathMeshLagP1,
                                 _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
        }
    }
    else
    {
        M_p2m.mesh()->load( _name=nameMeshLagP1,_path=pathMeshLagP1,
                             _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
    }

    VLOG(2) << "[P1 Lagrange] Pointset " << M_pset.points() << "\n";
}

template<typename space_type>
void
OperatorLagrangeP1<space_type>::buildLagrangeP1Mesh( bool parallelBuild )
{
    typedef typename image_mesh_type::point_type point_type;
    typedef typename image_mesh_type::element_type element_type;
    typedef typename image_mesh_type::face_type face_type;

    // inherit the table of markersName
    BOOST_FOREACH( auto itMark, this->domainSpace()->mesh()->markerNames() )
    {
        M_mesh->addMarkerName( itMark.first,itMark.second[0],itMark.second[1] );
    }

    bool doParallelBuild = (M_mesh->worldComm().size()==1)?false:parallelBuild;

    // iterate over active element on process
    auto it = this->domainSpace()->mesh()->beginElementWithProcessId();
    auto en = this->domainSpace()->mesh()->endElementWithProcessId();

    // memory nodes
    std::vector<size_type> new_node_numbers( this->domainSpace()->nLocalDof(),invalid_size_type_value );
    // the number of nodes on the new mesh, will be incremented
    size_type nNewNodes = 0, nNewElem  = 0, nNewFaces = 0;
    // memory for parallelism
    const int nProc = M_mesh->worldComm().localSize();
    std::map<size_type,size_type> mapGhostDofIdClusterToProcess;
    std::map<size_type,std::vector<size_type> > mapActiveEltWhichAreGhostInOtherPartition;
    std::vector<int> nbMsgToRecv( nProc , 0 );

    // use the geometric transformation to transform
    // the local equispace mesh to a submesh of the
    // current element
    gmc_ptrtype gmc;
    if ( it!=en )
        gmc = gmc_ptrtype( new gmc_type( this->domainSpace()->mesh()->gm(), *it, M_gmpc ) );

    for ( size_type elid = 0, pt_image_id = 0; it != en; ++it )
    {
        DVLOG(2) << "=========================================\n";
        DVLOG(2) << "global element " << it->id() << " oriented ok ? : " << it->isAnticlockwiseOriented() << "\n";
        DVLOG(2) << "global element G=" << it->G() << "\n";

        gmc->update( *it );

        // accumulate the local mesh in element *it in the new mesh
        auto itl = M_p2m.mesh()->beginElement();
        auto const enl = M_p2m.mesh()->endElement();

        // init parallel data
        if ( doParallelBuild && it->numberOfNeighborPartitions() > 0 )
        {
            // memory elt wich is send
            mapActiveEltWhichAreGhostInOtherPartition[it->id()] = std::vector<size_type>(std::distance(itl,enl),invalid_size_type_value);
            // counter of mpi msg recv
            auto itneighbor = it->neighborPartitionIds().begin();
            auto const enneighbor = it->neighborPartitionIds().end();
            for ( ; itneighbor!=enneighbor ; ++itneighbor )
                nbMsgToRecv[*itneighbor]++;
        }

        // get dofs in each faces of this ref element
        std::vector< std::set<size_type> > dofsInFace( it->numTopologicalFaces );
        for ( uint16_type f = 0; f < it->numTopologicalFaces ; f++ )
        {
            auto const& theFaceBase = it->face( f );
            for ( size_type localFaceDof = 0 ; localFaceDof < this->domainSpace()->dof()->nLocalDofOnFace() ; ++localFaceDof )
            {
                const size_type globFaceDof = this->domainSpace()->dof()->localToGlobal( theFaceBase,localFaceDof,0 ).template get<0>();
                if ( globFaceDof != invalid_size_type_value )
                    dofsInFace[f].insert( globFaceDof );
            }
        }

        // iterate on each new element (from ref elt)
        for ( int cptEltp2m=0 ; itl != enl; ++itl, ++elid,++cptEltp2m )
        {
            DVLOG(2) << "************************************\n";
            DVLOG(2) << "local elt = " << elid << "\n";
            DVLOG(2) << "local element " << itl->id() << " oriented ok ? : " << itl->isAnticlockwiseOriented() << "\n";
            DVLOG(2) << "local element G=" << itl->G() << "\n";

            element_type elt;
            elt.setId( elid );
            elt.setMarker( it->marker().value() );
            elt.setMarker2( it->marker2().value() );
            elt.setMarker3( it->marker3().value() );
            elt.setProcessIdInPartition( it->pidInPartition() );
            elt.setProcessId(it->processId());

            if ( doParallelBuild )
            {
                elt.setNumberOfPartitions( it->numberOfPartitions() );
                elt.setNeighborPartitionIds( it->neighborPartitionIds() );
            }
            else
            {
                elt.setNumberOfPartitions( 1 );
            }

            // accumulate the points
            for ( int p = 0; p < image_mesh_type::element_type::numVertices; ++p )
            {
                DVLOG(2) << "local In original element, vertex number " << itl->point( p ).id() << "\n";
                uint16_type localptid = itl->point( p ).id();
                uint16_type localptid_dof = localDof( it, localptid );
                size_type ptid = boost::get<0>( this->domainSpace()->dof()->localToGlobal( it->id(),
                                                localptid_dof, 0 ) );
                FEELPP_ASSERT( ptid < this->domainSpace()->nLocalDof()/domain_space_type::nComponents )( ptid )( this->domainSpace()->nLocalDof()/domain_space_type::nComponents ).warn( "invalid domain dof index" );

                if (doParallelBuild && !this->domainSpace()->dof()->dofGlobalClusterIsOnProc(this->domainSpace()->dof()->mapGlobalProcessToGlobalCluster(ptid)))
                    mapGhostDofIdClusterToProcess[this->domainSpace()->dof()->mapGlobalProcessToGlobalCluster(ptid)] = ptid;

                // add new node if not inserted before
                if ( new_node_numbers[ptid] == invalid_size_type_value )
                    {
                        new_node_numbers[ptid] = nNewNodes;

                        point_type __pt( nNewNodes, boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) )  );
                        __pt.setProcessId( it->processId() );
                        __pt.setProcessIdInPartition( it->pidInPartition() );


                        DVLOG(2) << "[OperatorLagrangeP1] element id "
                                      << elid << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] local point id "
                                      << localptid << " coord " << itl->point( p ).node() << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] point local id "
                                      << localptid << " global id " << ptid << "\n"
                                      << " localptid_dof = " << localptid_dof << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] adding point "
                                      << __pt.id() << " : " << __pt.node() << "\n";
                        DVLOG(2) << "[OperatorLagrangeP1] domain point gid "
                                      << boost::get<1>( this->domainSpace()->dof()->dofPoint( ptid ) )
                                      << " coords: " << boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) ) << "\n";

                        DLOG_IF( WARNING, ( ublas::norm_2( boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) )-__pt.node() ) > 1e-10 ) )
                            << "inconsistent point coordinates at pt index " << ptid << " in element " << elid << " with "
                            << " dot pt: " << boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) )
                            << " mesh  pt: " << __pt.node()
                            << " itl->id: " << itl->id()
                            << " local pt id : " << localptid
                            << " local pt id dof : " << localptid_dof << "\n";

                        M_mesh->addPoint( __pt );

                        nNewNodes++;
                    }
                //--------------------------------------------------------------------//

                elt.setPoint( p, M_mesh->point( new_node_numbers[ptid] ) );

#if 1
                FEELPP_ASSERT( ublas::norm_2( boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) )-elt.point( p ).node() ) < 1e-10 )
                ( boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) ) )
                ( p )
                ( elt.point( p ).node() )
                ( elid )
                ( itl->id() )
                ( localptid )
                ( localptid_dof )
                ( ptid ).warn( "[after] inconsistent point coordinates" );

                FEELPP_ASSERT( ublas::norm_2( elt.point( p ).node()-ublas::column( elt.G(), p ) ) < 1e-10 )
                ( p )
                ( elt.point( p ).node() )
                ( elid )
                ( elt.G() ).warn( "[after] inconsistent point coordinates" );
#endif
            }

#if 0
            static const int nDim = element_type::nDim;
            ublas::vector<double> D( nDim + 1, 0 );
            typename matrix_node<value_type>::type M( nDim+1,nDim+1 );
            typename matrix_node<value_type>::type P( nDim,nDim+1 );
            typename matrix_node<value_type>::type P0( nDim,nDim+1 );
            std::fill( M.data().begin(), M.data().end(), value_type( 1.0 ) );

            for ( int p = 0; p < element_type::numVertices; ++p )
            {
                ublas::column( P0, p ) = elt.vertex( p );
            }

            ublas::subrange( M, 0, nDim, 0, nDim+1 ) = P0;
            double meas_times = details::det( M, mpl::int_<nDim+1>() );
            FEELPP_ASSERT( meas_times > 0 )( elt.id() )( elt.G() ).warn( "negative area" );
#else
            //FEELPP_ASSERT( elt.isAnticlockwiseOriented() )( elt.id() )( elt.G() ).warn( "negative area" );
            //FEELPP_ASSERT( ublas::norm_inf( elt.G()- it->G() ) < 1e-10 )( it->G() )( elt.G() ).warn( "global: not same element" );
            //FEELPP_ASSERT( ublas::norm_inf( elt.G()- itl->G() ) < 1e-10 )( itl->G() )( elt.G() ).warn( "local: not same element" );
#endif

            // set id of element and increment the element counter
            elt.setId ( nNewElem++ );
            // add element in mesh
            auto const& theNewElt = M_mesh->addElement ( elt );

            // save data if elt is connected to another partition
            if (doParallelBuild && it->numberOfNeighborPartitions() > 0)
                mapActiveEltWhichAreGhostInOtherPartition[it->id()][cptEltp2m] = theNewElt.id();

#if 1
            // Maybe add faces for this element
            for ( uint16_type s=0; s<itl->numTopologicalFaces; s++ )
            {
                if ( !itl->facePtr( s ) ) continue;

                auto const& theFaceBase = itl->face( s );

                face_type newFace;
                newFace.setOnBoundary( true );
                newFace.setProcessId( it->processId() );
                newFace.setProcessIdInPartition( it->pidInPartition() );
                newFace.setProcessId(it->processId());
                if ( doParallelBuild )
                {
                    newFace.setNumberOfPartitions( it->numberOfPartitions() );
                    newFace.setNeighborPartitionIds( it->neighborPartitionIds() );
                }
                else
                {
                    newFace.setNumberOfPartitions( 1 );
                }

                // set points in face and up counter for connecting with ref faces
                std::vector<uint16_type> nPtInRefFace(it->numTopologicalFaces,0);
                for ( uint16_type p = 0; p < theFaceBase.nPoints(); ++p )
                {
                    uint16_type localptidFace = theFaceBase.point( p ).id();
                    uint16_type localptidFace_dof = localDof( it, localptidFace );
                    const size_type theglobdof = this->domainSpace()->dof()->localToGlobal( it->id(),localptidFace_dof,0 ).template get<0>();
                    newFace.setPoint( p, M_mesh->point( new_node_numbers[theglobdof] ) );
                    // update face point connection with reference faces
                    for ( uint16_type fId = 0; fId <  dofsInFace.size() ; ++fId )
                    {
                        if ( dofsInFace[fId].find( theglobdof ) != dofsInFace[fId].end() )
                            nPtInRefFace[fId]++;
                    }
                }

                // update faces from reference faces
                for ( uint16_type fId = 0; fId < nPtInRefFace.size() ; ++fId )
                {
                    // find connection if all points are in reference faces
                    if ( nPtInRefFace[fId] == theFaceBase.nPoints() )
                    {
                        // update marker from ref
                        newFace.setMarker( it->face( fId ).marker().value() );
                        newFace.setMarker2( it->face( fId ).marker2().value() );
                        newFace.setMarker3( it->face( fId ).marker3().value() );
                    }
                }

                // set id of element and increment the face counter
                newFace.setId( nNewFaces++ );
                // add it to the list of faces
                auto addFaceRes = M_mesh->addFace( newFace );
            }
#endif

        } //  for ( int cptEltp2m=0 ; itl != enl; ++itl, ++elid,++cptEltp2m )






    } // for ( size_type elid = 0, pt_image_id = 0; it != en; ++it )



    if ( doParallelBuild )
    {
        VLOG(2) << "[P1 Lagrange] start parallel build \n";

        //auto const theWorldCommSize = M_mesh->worldComm().size();
        std::vector<int> nbMsgToSend( nProc , 0 );
        std::vector< std::map<int,size_type> > mapMsg( nProc );

        auto iv = this->domainSpace()->mesh()->beginGhostElement();
        auto const en = this->domainSpace()->mesh()->endGhostElement();
        for ( ; iv != en; ++iv )
        {
            auto const procGhost = iv->processId();
            auto const idEltOnGhost = iv->idInOthersPartitions(procGhost);
            M_mesh->worldComm().localComm().send(procGhost, nbMsgToSend[procGhost], idEltOnGhost);
            mapMsg[procGhost].insert( std::make_pair( nbMsgToSend[procGhost],iv->id() ) );
            ++nbMsgToSend[procGhost];
        }

#if 0
        // counter of msg received for each process
        std::vector<int> nbMsgToRecv2;
        mpi::all_to_all( M_mesh->worldComm().localComm(),
                         nbMsgToSend,
                         nbMsgToRecv2 );
        for ( int proc=0; proc<nProc; ++proc )
            CHECK( nbMsgToRecv[proc]==nbMsgToRecv2[proc] ) << "paritioning data incorect "
                                                           << "myrank " << M_mesh->worldComm().localRank() << " proc " << proc
                                                           << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
                                                           << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc] << "\n";
#endif

        VLOG(2) << "[P1 Lagrange] parallel build : finish first send \n";

        // counter of request
        int nbRequest=0;
        for ( int proc=0; proc<nProc; ++proc )
            nbRequest+=nbMsgToRecv[proc]+nbMsgToSend[proc];

        mpi::request * reqs = new mpi::request[nbRequest];
        int cptRequest=0;

        // recv dof asked and re-send dof in this proc
        for ( int proc=0; proc<nProc; ++proc )
        {
            for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
            {
                //recv
                size_type idEltRecv;
                M_mesh->worldComm().localComm().recv( proc, cpt, idEltRecv );
                auto const& theeltIt = this->domainSpace()->mesh()->elementIterator(idEltRecv);

                DCHECK(mapActiveEltWhichAreGhostInOtherPartition.find(idEltRecv) != mapActiveEltWhichAreGhostInOtherPartition.end() ) << "invalid idEltRecv " << idEltRecv << "\n";
                auto const& idOfNewElt = mapActiveEltWhichAreGhostInOtherPartition[idEltRecv];

                std::vector<boost::tuple<size_type,ublas::vector<double> > > resultClusterDofsAndNodesToSend(M_p2m.mesh()->numPoints());

                auto itp = M_p2m.mesh()->beginPoint();
                auto const enp = M_p2m.mesh()->endPoint();
                for ( ; itp!=enp ; ++itp )
                {
                    uint16_type localptid = itp->id();
                    uint16_type localptid_dof = localDof( theeltIt, localptid );
                    size_type ptid = boost::get<0>( this->domainSpace()->dof()->localToGlobal( theeltIt->id(),localptid_dof, 0 ) );
                    size_type idInProcAsked = this->domainSpace()->dof()->mapGlobalProcessToGlobalCluster(ptid);
                    ublas::vector<double> thedofnode = boost::get<0>( this->domainSpace()->dof()->dofPoint( ptid ) );
                    resultClusterDofsAndNodesToSend[localptid] = boost::make_tuple( idInProcAsked, thedofnode);
                }

                auto itl = M_p2m.mesh()->beginElement();
                auto const enl = M_p2m.mesh()->endElement();
                // localIdNewElt -> ( localIdPt, globalIdNewElt )
                std::vector< boost::tuple< std::vector<uint16_type>, size_type > > resultToSendBis(std::distance(itl,enl));
                for ( size_type elidp2m = 0 ; itl != enl; ++itl, ++elidp2m )
                {
                    std::vector<uint16_type> vecLocalIdPtToSend( image_mesh_type::element_type::numVertices );
                    // accumulate the points
                    for ( int p = 0; p < image_mesh_type::element_type::numVertices; ++p )
                    {
                        DVLOG(2) << "local In original element, vertex number " << itl->point( p ).id() << "\n";
                        const uint16_type localptid = itl->point( p ).id();
                        vecLocalIdPtToSend[p]=localptid;
                    }

                    const size_type idEltInNewMesh = idOfNewElt[elidp2m];
                    resultToSendBis[elidp2m] = boost::make_tuple( vecLocalIdPtToSend,idEltInNewMesh );
                    //auto eltToUpdate = M_mesh->elementIterator( idEltInNewMesh );
                    //M_mesh->elements().modify( eltToUpdate, Feel::detail::UpdateNeighborPartition( proc ) );
                }
                auto resultToSend = boost::make_tuple(resultToSendBis,resultClusterDofsAndNodesToSend);
#if 0
                M_mesh->worldComm().localComm().send( proc, cpt, resultToSend);
#else
                reqs[cptRequest] = M_mesh->worldComm().localComm().isend( proc, cpt, resultToSend);
                ++cptRequest;
#endif

            } // for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
        } // for ( int proc=0; proc<nProc; ++proc )

        std::map<int,std::vector< boost::tuple< std::vector< boost::tuple< std::vector<uint16_type>, size_type > >,
                                                std::vector<boost::tuple<size_type,ublas::vector<double> > > > > > MAPresultRecvData;


        for ( int proc=0; proc<nProc; ++proc )
        {
            MAPresultRecvData[proc].resize(nbMsgToSend[proc]);
            for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
            {
                reqs[cptRequest] = M_mesh->worldComm().localComm().irecv( proc, cpt, MAPresultRecvData[proc][cpt] );
                ++cptRequest;
            }
        }
        mpi::wait_all(reqs, reqs + nbRequest);

        VLOG(2) << "[P1 Lagrange] parallel build : finish first recv and resend response \n";

        delete [] reqs;

#if 1
        std::map<size_type,size_type> new_ghost_node_numbers;

        // get response to initial request
        for ( int proc=0; proc<nProc; ++proc )
        {
            for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
            {
#if 0
                boost::tuple<  std::vector< boost::tuple< std::vector<uint16_type>, size_type > >,
                               std::vector<boost::tuple<size_type,ublas::vector<double> > >  > resultRecvData;
                M_mesh->worldComm().localComm().recv( proc, cpt, resultRecvData );
#else
                auto resultRecvData = MAPresultRecvData[proc][cpt];
#endif

                auto const& myGhostEltBase = this->domainSpace()->mesh()->element( mapMsg[proc][cpt],proc );

                auto const& mapLocalToGlobalPointId = resultRecvData.template get<1>();
                auto resultRecv = resultRecvData.template get<0>();

                auto itelt = resultRecv.begin();
                auto const enelt = resultRecv.end();
                for ( ;itelt!=enelt;++itelt )
                {
                    element_type elt;
                    elt.setMarker( myGhostEltBase.marker().value() );
                    elt.setMarker2( myGhostEltBase.marker2().value() );
                    elt.setMarker3( myGhostEltBase.marker3().value() );
                    elt.setProcessIdInPartition( myGhostEltBase.pidInPartition() );
#if 0
                    elt.setNumberOfPartitions(2);
                    DCHECK( proc==myGhostEltBase.processId() ) << "invalid process id\n";
                    elt.setProcessId( proc );
                    elt.setNeighborPartitionIds( std::vector<int>({M_mesh->worldComm().localRank()}) );
                    /*std::vector<int> newNeighborPartitionIds(1);
                    newNeighborPartitionIds[0]=M_mesh->worldComm().localRank();
                    elt.setNeighborPartitionIds( newNeighborPartitionIds );*/

#else
                    elt.setNumberOfPartitions( myGhostEltBase.numberOfPartitions() );
                    DCHECK( proc==myGhostEltBase.processId() ) << "invalid process id\n";
                    elt.setProcessId( proc );
                    elt.setNeighborPartitionIds( myGhostEltBase.neighborPartitionIds() );
#endif

                    auto itpt = itelt->template get<0>().begin();
                    auto const enpt = itelt->template get<0>().end();
                    for ( int p = 0 ; itpt!=enpt ; ++itpt,++p )
                    {
                        auto const globclusterdofRecv = mapLocalToGlobalPointId[(int)*itpt].template get<0>();
                        size_type thenewptid=invalid_size_type_value;
                        if ( this->domainSpace()->dof()->dofGlobalClusterIsOnProc(globclusterdofRecv))
                        {
                            const size_type theoldptid = this->domainSpace()->dof()->mapGlobalClusterToGlobalProcess( globclusterdofRecv - this->domainSpace()->dof()->firstDofGlobalCluster() );
                            thenewptid = new_node_numbers[theoldptid];
                            CHECK( thenewptid!=invalid_size_type_value ) << "--1---invalid point id theoldptid="<<theoldptid << " thenewptid="<<thenewptid<< "\n";
                        }
                        else if (mapGhostDofIdClusterToProcess.find(globclusterdofRecv) != mapGhostDofIdClusterToProcess.end() )
                        {
                            const size_type theoldptid = mapGhostDofIdClusterToProcess[globclusterdofRecv];
                            thenewptid = new_node_numbers[theoldptid];
                            CHECK( thenewptid!=invalid_size_type_value ) << "--2---invalid point id\n";
                        }
                        else
                        {
                            if (new_ghost_node_numbers.find(globclusterdofRecv) == new_ghost_node_numbers.end() )
                            {
                                // add a new ghost point
                                new_ghost_node_numbers[globclusterdofRecv] = nNewNodes;
                                thenewptid = nNewNodes;
                                CHECK( thenewptid!=invalid_size_type_value ) << "--3---invalid point id\n";
                                point_type __pt( nNewNodes, mapLocalToGlobalPointId[(int)*itpt].template get<1>() );
                                __pt.setProcessId( invalid_uint16_type_value );
                                __pt.setProcessIdInPartition( myGhostEltBase.pidInPartition() );
                                M_mesh->addPoint( __pt );
                                nNewNodes++;
                            }
                            else
                            {
                                // use present ghost point
                                thenewptid = new_ghost_node_numbers[globclusterdofRecv];
                                CHECK( thenewptid!=invalid_size_type_value ) << "--4---invalid point id\n";
                            }
                            //elt.setPoint( p, M_mesh->point( new_ghost_node_numbers[globclusterdofRecv] ) );
                        }
                        CHECK( thenewptid!=invalid_size_type_value ) << "invalid point id\n";
                        // add point to elt
                        elt.setPoint( p, M_mesh->point( thenewptid ) );

                    } // for ( int p = 0 ; itpt!=enpt ; ++itpt,++p )

                    // set id of element
                    elt.setId ( nNewElem++ );

                    // increment the new element counter
                    nNewElem++;

                    auto const& theNewElt = M_mesh->addElement ( elt );
                    M_mesh->elements().modify( M_mesh->elementIterator( theNewElt.id(),  proc ),
                                                Feel::detail::updateIdInOthersPartitions( proc, itelt->template get<1>() ) );

                } // for ( ;itelt!=enelt;++itelt )
            } // for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        } // for ( int proc=0; proc<nProc; ++proc )

#endif
    } // if ( parallelBuild )

    DVLOG(2) << "dist ghost elt base " << std::distance( this->domainSpace()->mesh()->beginGhostElement(),this->domainSpace()->mesh()->endGhostElement() )
             << " dist ghost elt new " << std::distance( M_mesh->beginGhostElement(),M_mesh->endGhostElement() ) << "\n";

    DVLOG(2) << "[P1 Lagrange] Number of points in mesh: " << M_mesh->numPoints() << "\n";

    M_mesh->setNumVertices( M_mesh->numPoints() );
    M_mesh->components().reset();
    M_mesh->components().set ( MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    M_mesh->updateForUse();
}


template<typename space_type>
void
OperatorLagrangeP1<space_type>::buildOperator()
{
    // construct the p1 space and set the operator

    //auto Xh_image = dual_image_space_type::New(_mesh=M_mesh,_worldscomm=worldsComm);
    auto Xh_image = dual_image_space_type::New(_mesh=M_mesh,_worldscomm=std::vector<WorldComm>(dual_image_space_type::nSpaces,M_mesh->worldComm() ) );

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
        gmc_ptrtype gmc( new gmc_type( this->domainSpace()->mesh()->gm(), *it, M_gmpc ) );
        typename matrix_node<value_type>::type u_at_pts = u.id( *gmc );

        typename std::list<size_type>::const_iterator ite = M_el2el[it->id()].begin();
        typename std::list<size_type>::const_iterator ene = M_el2el[it->id()].end();

        for ( ; ite != ene; ++ite )
        {
            typename domain_mesh_type::element_type const& elt = M_mesh->element( *ite );

            for ( int c = 0; c < nComponents; ++c )
                for ( int p = 0; p < domain_mesh_type::element_type::numVertices; ++p )
                {
                    size_type ptid = boost::get<0>( this->dualImageSpace()->dof()->localToGlobal( elt.id(), p, c ) );
                    res( ptid ) = ublas::column( u_at_pts, M_el2pt[*ite][p] )( nComponents*c+c );
                }
        }
    }

    return res;
}
#endif

template<typename space_type>
uint16_type
OperatorLagrangeP1<space_type>::localDof( typename domain_mesh_type::element_const_iterator it,
        uint16_type localptid,
        mpl::bool_<true> ) const
{
    typedef typename domain_mesh_type::element_type::edge_permutation_type edge_permutation_type;
    uint16_type localptid_dof = localptid;

    if ( domain_mesh_type::nDim == 2 &&
            M_pset.pointToEntity( localptid ).first == 1 &&
            it->edgePermutation( M_pset.pointToEntity( localptid ).second ) == edge_permutation_type::REVERSE_PERMUTATION )
    {
        //size_type edge_id = it->edge( M_pset.pointToEntity( localptid ).second ).id();
        uint16_type edge_id = M_pset.pointToEntity( localptid ).second;

        uint16_type l = localptid - ( domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
                                      edge_id * domain_fe_type::nDofPerEdge );

        FEELPP_ASSERT( l != invalid_uint16_type_value && ( l < domain_fe_type::nDofPerEdge ) )
        ( l )
        ( domain_fe_type::nDofPerEdge )
        ( localptid )
        ( M_pset.pointToEntity( localptid ).second )
        ( domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
          edge_id * domain_fe_type::nDofPerEdge ).error( "invalid value" );

        localptid_dof = ( domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
                          edge_id * domain_fe_type::nDofPerEdge +
                          domain_fe_type::nDofPerEdge - 1 - l );

        FEELPP_ASSERT( localptid_dof != invalid_uint16_type_value &&
                       ( localptid_dof < domain_fe_type::nLocalDof ) )
        ( l )
        ( localptid )
        ( localptid_dof )
        ( domain_fe_type::nLocalDof )
        ( M_pset.pointToEntity( localptid ).second )
        ( domain_mesh_type::element_type::numVertices * domain_fe_type::nDofPerVertex +
          edge_id * domain_fe_type::nDofPerEdge ).error( "invalid value" );
    }

    return localptid_dof;
}

//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------//
/**
 * \return the P1 Lagrange adaptor  associated to the space \p Xh
 */
template<typename space_type>
boost::shared_ptr<OperatorLagrangeP1<space_type> >
opLagrangeP1_impl( boost::shared_ptr<space_type> const& Xh,
                   typename OperatorLagrangeP1<space_type>::backend_ptrtype const& backend,
                   std::string pathMeshLagP1,
                   std::string prefix,
                   bool rebuild,
                   bool parallel
                   )
{
    return boost::shared_ptr<OperatorLagrangeP1<space_type> >( new OperatorLagrangeP1<space_type>( Xh,backend,pathMeshLagP1,prefix,rebuild,parallel ) );
}


template<typename Args>
struct compute_opLagrangeP1_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::space>::type>::type::element_type space_type;
    typedef boost::shared_ptr<OperatorLagrangeP1<space_type> > type;
};


BOOST_PARAMETER_FUNCTION(
    ( typename compute_opLagrangeP1_return<Args>::type ), // 1. return type
    lagrangeP1,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( space,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    ) // required
    ( optional
      ( backend,    *, Backend<typename compute_opLagrangeP1_return<Args>::space_type::value_type>::build( soption( _name="backend" ) ) )
      ( path,       *( boost::is_convertible<mpl::_,std::string> ), std::string(".") )
      ( prefix,     *( boost::is_convertible<mpl::_,std::string> ), std::string("") )
      ( rebuild,    *( boost::is_integral<mpl::_> ), 1 )
      ( parallel,   *( boost::is_integral<mpl::_> ), 1 )
    ) // optionnal
)
{
    Feel::detail::ignore_unused_variable_warning( args );
    return opLagrangeP1_impl(space,backend,path,prefix,rebuild,parallel);
}










} // Feel
#endif /* __OperatorLagrangeP1_H */
