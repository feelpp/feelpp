/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
  This file is part of the Feel library

  Copyright (C) 2012 Université de Grenoble 1

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
 * \file doftablempi.hpp
 * \author Vincent Chabannes
 */
#if !defined(FEELPP_DOFTABLE_MPI_HPP)
#define FEELPP_DOFTABLE_MPI_HPP 1

namespace Feel
{

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostDofMap( mesh_type& mesh )
{
    wc(mesh)->print(fmt::format("[DofTable::buildGhostDofMap rank={}] starts. hasMeshSupport: {}", rank(mesh), this->hasMeshSupport()), Environment::logVerbosityLevel() > 1, Environment::logVerbosityLevel() > 0, Environment::logVerbosityLevel() > 1 );

    this->buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh );

    DVLOG(2) << "[buildGhostDofMap] finish () with rank "<< this->worldComm().rank();
}


//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::generateDofPoints( Range<mesh_type,MESH_ELEMENTS> const& myrange ) const
{
    if ( fe_type::is_modal )
        return;

    DVLOG(2) << "[Dof::generateDofPoints] generating dof coordinates\n";
    typedef typename gm_type::template Context<element_type> gm_context_type;
    typedef std::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    gm_ptrtype gm( new gm_type );
    fe_type fe;

    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, fe.points() ) );

    gm_context_ptrtype __c;

    std::vector<bool> dof_done( this->nLocalDofWithGhost(), false );

    for (auto const& eltWrap : myrange )
    {
        auto const& elt = unwrap_ref( eltWrap );

        if ( __c )
            __c->template update<vm::POINT>( elt );
        else
            __c = gm->template context<vm::POINT>( elt, __geopc );

        for ( auto const& ldof : this->localDof( elt.id() ) )
        {
            size_type thedof = ldof.second.index();
            if ( dof_done[thedof] )
                continue;
            dof_done[thedof] = true;

            uint16_type ldofId = ldof.first.localDof();
            uint16_type ldofParentId = this->fe().dofParent( ldofId );

            if ( ( thedof >= this->firstDof() ) && ( thedof <= this->lastDof() ) )
            {
                DCHECK( thedof < this->nLocalDofWithGhost() )
                    << "invalid local dof index "
                    <<  thedof << ", " << this->nLocalDofWithGhost() << "," << this->firstDof()  << ","
                    <<  this->lastDof() << "," << elt.id() << "," << ldofId << "," << ldofParentId;

                uint16_type comp = this->fe().component( ldofId );
                M_dof_points[thedof] = boost::make_tuple( __c->xReal( ldofParentId ), thedof, comp );
            }
        }
    }

}


template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh_type& mesh )
{
    DVLOG(2) << "[buildGlobalProcessToGlobalClusterDofMapOthersMesh] start\n";

    static const uint16_type nLocalDofUpToVertices = element_type::numVertices*fe_type::nDofPerVertex;
    static const uint16_type nLocalDofUpToEdges = element_type::numVertices*fe_type::nDofPerVertex + element_type::numEdges*fe_type::nDofPerEdge;
    static const uint16_type nLocalDofUpToFaces =
        (nDim==1)? element_type::numVertices*fe_type::nDofPerVertex :
        (nDim==2)? element_type::numVertices*fe_type::nDofPerVertex + element_type::numEdges*fe_type::nDofPerEdge :
        element_type::numVertices*fe_type::nDofPerVertex + element_type::numEdges*fe_type::nDofPerEdge + element_type::numGeometricFaces*fe_type::nDofPerFace;
    static const uint16_type nDofPerVertexForDivision = mpl::if_<boost::is_same<mpl::int_<fe_type::nDofPerVertex>,mpl::int_<0> >,
                                                                 mpl::int_<1>,
                                                                 mpl::int_<fe_type::nDofPerVertex> >::type::value;

    static const uint16_type nDofPerEdgeForDivision = mpl::if_<boost::is_same<mpl::int_<fe_type::nDofPerEdge>,mpl::int_<0> >,
                                                               mpl::int_<1>,
                                                               mpl::int_<fe_type::nDofPerEdge> >::type::value;
    static const uint16_type nLocalDofBeforeDofsOfTopologicalFaces =
        (nDim==1)? 0 :
        (nDim==2)? element_type::numVertices*fe_type::nDofPerVertex :
        element_type::numVertices*fe_type::nDofPerVertex + element_type::numEdges*fe_type::nDofPerEdge;

    static const uint16_type nDofPerTopologicalFace =
        (nDim==1)?fe_type::nDofPerVertex : (nDim==2)?fe_type::nDofPerEdge : fe_type::nDofPerFace;
    static const uint16_type nDofPerTopologicalFaceForDivision = mpl::if_<boost::is_same<mpl::int_<nDofPerTopologicalFace>,mpl::int_<0> >,
                                                                          mpl::int_<1>,
                                                                          mpl::int_<nDofPerTopologicalFace> >::type::value;

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    std::set<uint16_type> localDofUsedForTensor2symm;
    if ( is_tensor2symm )
    {
        std::map<uint16_type,std::vector<uint16_type> > symm2unsymm;
        for ( uint16_type k=0;k<this->nLocalDof();++k )
            symm2unsymm[this->fe().unsymmToSymm(k)].push_back( k );
        for ( auto const& symmdof : symm2unsymm )
        {
            if ( symmdof.second.empty() )
                continue;
            localDofUsedForTensor2symm.insert( symmdof.second.front() );
        }
    }


    // create inverse mapping of vector_permutation (vector_permutation map elt ordering to face ordering)
    // thus, this mapping will map face ordering to elt ordering with respect to permutation
    std::map<face_permutation_type, permutation_vector_type> mapLocalDofFaceToFaceInElement;
    if ( nDim == 3 && nDofPerTopologicalFace > 1 )
    {
        for ( auto const& [perm,dofsMapping] : this->vector_permutation )
        {
            mapLocalDofFaceToFaceInElement[perm].resize(dofsMapping.size());
            for (int k=0;k<dofsMapping.size();++k)
                mapLocalDofFaceToFaceInElement[perm][dofsMapping[k]] = k;
        }
    }


    size_type nLocalDofWithGhost = this->M_n_localWithGhost_df[myRank];
    std::vector<bool> dofdone( nLocalDofWithGhost,false);
    std::vector<bool> dofIsGhost( nLocalDofWithGhost,false);
    size_type nDofNotPresent=0;


    // maybe a vector instead of of map??
    std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > dataToSend;

    //typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myActiveEltsTouchInterProcess( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );

    bool hasMeshSupportPartial = this->hasMeshSupport() && this->meshSupport()->isPartialSupport();

    std::map<rank_type, std::map<size_type,std::vector<size_type> > > dataMemory;
    if ( is_continuous )
    {
        std::map<size_type, std::tuple<rank_type,size_type,uint16_type>> mapPtIdToDofOwnerElt; //( ptId ->( rank, eltId, ptIdInElt  ))
        std::map<size_type, std::tuple<rank_type,size_type,uint16_type,edge_permutation_type>> mapEdgeIdToDofOwnerElt; //( ptId ->( rank, eltId, edgeIdInElt, edge permutation ))
        std::map<size_type, std::tuple<rank_type,size_type,uint16_type,face_permutation_type>> mapFaceIdToDofOwnerElt; //( faceId ->( rank, eltId, faceIdInElt, face permutation))
        std::vector< std::reference_wrapper<const typename mesh_type::element_type> > activeEltTouchInterprocess;
        auto rangeElements = hasMeshSupportPartial? elements( this->meshSupport(), entity_process_t::ALL ) : elements( mesh, entity_process_t::ALL );
        for ( auto const& eltWrap : rangeElements )
        {
            auto const& elt = unwrap_ref( eltWrap );
            rank_type eltPid = elt.processId();
            size_type eltIdOtherPart = elt.idInOthersPartitions( eltPid );
            bool currentEltTouchInterprocess = false;
            for ( uint16_type n=0; n < elt.nVertices(); n++ )
            {
                auto const& point = elt.point(n);
                bool isInterprocessPoint = hasMeshSupportPartial? this->meshSupport()->isInterprocessPoints( point.id() ) : mesh.isInterprocessPoints( point.id() );
                if ( isInterprocessPoint )
                {
                    currentEltTouchInterprocess = true;
                    auto itFind = mapPtIdToDofOwnerElt.find( point.id() );
                    if ( itFind == mapPtIdToDofOwnerElt.end() )
                        mapPtIdToDofOwnerElt.emplace( point.id(), std::make_tuple( eltPid, eltIdOtherPart, n ) );
                    else
                    {
                        auto & currentData = itFind->second;
                        if ( elt.processId() < std::get<0>( currentData ) )
                            currentData = std::make_tuple( eltPid, eltIdOtherPart, n );
                    }
                }
            }
            if constexpr ( nDim == 3 )
            {
                if (  element_type::numEdges*fe_type::nDofPerEdge > 0 && currentEltTouchInterprocess )
                {
                    for ( size_type j = 0; j < elt.nEdges(); j++ )
                    {
                        if ( !elt.edgePtr( j ) )
                            continue;
                        auto & edge = elt.edge( j );

                        bool isInterprocessEdge = hasMeshSupportPartial? this->meshSupport()->isInterprocessEdges( edge.id() ) : mesh.isInterprocessEdges( edge.id() );
                        if ( isInterprocessEdge )
                        {
                            auto itFind = mapEdgeIdToDofOwnerElt.find( edge.id() );
                            if ( itFind == mapEdgeIdToDofOwnerElt.end() )
                                mapEdgeIdToDofOwnerElt.emplace( edge.id(), std::make_tuple( eltPid, eltIdOtherPart, j, elt.edgePermutation( j ).value() ) );
                            else
                            {
                                auto & currentData = itFind->second;
                                if ( elt.processId() < std::get<0>( currentData ) )
                                    currentData = std::make_tuple( eltPid, eltIdOtherPart, j, elt.edgePermutation( j ).value() );
                            }
                        }
                    }
                }
            }

            if ( nDofPerTopologicalFace > 0 && currentEltTouchInterprocess )
            {
                // face
                for ( uint16_type j = 0; j < elt.nTopologicalFaces(); j++ )
                {
                    if ( !elt.facePtr( j ) )
                        continue;
                    auto const& face = elt.face( j );

                    // TODO: maybe ignore face that don't touch interprocess??

                    auto itFind = mapFaceIdToDofOwnerElt.find( face.id() );
                    if ( itFind == mapFaceIdToDofOwnerElt.end() )
                        mapFaceIdToDofOwnerElt.emplace( face.id(), std::make_tuple( eltPid, eltIdOtherPart, j, elt.facePermutation( j ).value() ) );
                    else
                    {
                        auto & currentData = itFind->second;
                        if ( elt.processId() < std::get<0>( currentData ) )
                            currentData = std::make_tuple( eltPid, eltIdOtherPart, j, elt.facePermutation( j ).value() );
                    }
                }
            }

            if ( !elt.isGhostCell() && currentEltTouchInterprocess )
                activeEltTouchInterprocess.push_back( std::cref( elt ) );
        }

        for ( auto const& eltWrap : activeEltTouchInterprocess )
        {
            auto const& elt = eltWrap.get();
            // prepare local dofs to analyse
            std::map<uint16_type,std::map<uint16_type,std::tuple<uint16_type,size_type> > > mapLocalDofByCompToLocalDofAllComp;
            for ( auto const& ldof : this->localDof( elt.id() ) )
            {
                uint16_type ldofId = ldof.first.localDof();
                if ( is_tensor2symm )
                    if ( localDofUsedForTensor2symm.find( ldofId ) == localDofUsedForTensor2symm.end() )
                        continue;
                uint16_type ldofParentId = this->fe().dofParent( ldofId );
                size_type gdofId = ldof.second.index();
                uint16_type comp = this->fe().component(ldofId);
                mapLocalDofByCompToLocalDofAllComp[ldofParentId][comp] = std::make_tuple( ldofId,gdofId );
            }
            // loop over local dof for the detection of ghost dofs
            for (auto const& localDofDatas : mapLocalDofByCompToLocalDofAllComp )
            {
                CHECK( !localDofDatas.second.empty() ) << "no localdof data is empty";
                auto const& localDofDataFirstComponent = *localDofDatas.second.begin();
                const uint16_type locDof = localDofDatas.first;
                //const uint16_type locDof = std::get<0>( localDofDataFirstComponent.second );
                const size_type theglobdoftest = std::get<1>( localDofDataFirstComponent.second );
                DCHECK( theglobdoftest < nLocalDofWithGhost ) << "invalid globdof " << theglobdoftest << "\n";

                if ( dofdone[theglobdoftest] )
                    continue;
                dofdone[theglobdoftest]=true;

                rank_type eltPidOwnerDof = invalid_v<rank_type>;
                size_type eltIdOwnerDof = invalid_v<size_type>;
                uint16_type localDofOwnerDof = invalid_v<uint16_type>;
                if ( locDof < nLocalDofUpToVertices )
                {
                    uint16_type pointIdInElt = locDof / nDofPerVertexForDivision;
                    auto const& point = elt.point( pointIdInElt );
                    auto itFind = mapPtIdToDofOwnerElt.find( point.id() );
                    if ( itFind == mapPtIdToDofOwnerElt.end() )
                        continue;
                    std::tie( eltPidOwnerDof, eltIdOwnerDof, localDofOwnerDof ) = itFind->second;
                }
                else if ( nDim == 3 && locDof < nLocalDofUpToEdges )
                {
                    if constexpr (nDim == 3)
                    {
                        uint16_type locDofInEgdes = locDof - nLocalDofUpToVertices;
                        uint16_type edgeIdInElt = locDofInEgdes / nDofPerEdgeForDivision;
                        uint16_type locDofInEgde = locDofInEgdes % nDofPerEdgeForDivision;

                        auto edgePtr = elt.edgePtr( edgeIdInElt );
                        if ( !edgePtr )
                            continue;
                        auto const& edge = *edgePtr;
                        auto itFind = mapEdgeIdToDofOwnerElt.find( edge.id() );
                        if ( itFind == mapEdgeIdToDofOwnerElt.end() )
                            continue;

                        edge_permutation_type edgePermutationOwnerDof;
                        uint16_type edgeIdInEltOwnerDof = invalid_v<uint16_type>;
                        std::tie( eltPidOwnerDof, eltIdOwnerDof, edgeIdInEltOwnerDof, edgePermutationOwnerDof ) = itFind->second;
                        // do nothing if owner is on curent process id
                        if ( eltPidOwnerDof == myRank )
                            continue;

                        auto edgePermutation = elt.edgePermutation( edgeIdInElt );
                        // shift local dof to relative dof edge numbering
                        localDofOwnerDof = nLocalDofUpToVertices + edgeIdInEltOwnerDof*nDofPerEdge;
                        localDofOwnerDof += ( edgePermutation == edgePermutationOwnerDof )? locDofInEgde : nDofPerEdge-1-locDofInEgde;
                    }
                }
                else if ( locDof < nLocalDofUpToFaces )
                {
                    if ( nDim < 2 && nDim != nRealDim )
                        continue;
                    uint16_type locDofInTopologicalFaces = locDof - nLocalDofBeforeDofsOfTopologicalFaces;
                    uint16_type faceIdInElt = locDofInTopologicalFaces / nDofPerTopologicalFaceForDivision;
                    uint16_type locDofInTopologicalFace = locDofInTopologicalFaces % nDofPerTopologicalFaceForDivision;

                    auto facePtr = elt.facePtr( faceIdInElt );
                    if ( !facePtr )
                        continue;
                    auto const& face = *facePtr;

                    auto itFind = mapFaceIdToDofOwnerElt.find( face.id() );
                    if ( itFind == mapFaceIdToDofOwnerElt.end() )
                        continue;

                    face_permutation_type facePermutationOwnerDof;
                    uint16_type faceIdInEltOwnerDof = invalid_v<uint16_type>;
                    std::tie( eltPidOwnerDof, eltIdOwnerDof, faceIdInEltOwnerDof, facePermutationOwnerDof ) = itFind->second;

                    // do nothing if owner is on curent process id
                    if ( eltPidOwnerDof == myRank )
                        continue;

                    // shift local dof to relative dof face numbering
                    localDofOwnerDof = nLocalDofBeforeDofsOfTopologicalFaces + faceIdInEltOwnerDof*nDofPerTopologicalFace;

                    if constexpr ( nDofPerTopologicalFace > 1 )
                    {
                        auto facePermutation = elt.facePermutation( faceIdInElt );
                        if constexpr ( nDim <= 2)
                        {
                            localDofOwnerDof += ( facePermutation == facePermutationOwnerDof )? locDofInTopologicalFace : nDofPerTopologicalFace-1-locDofInTopologicalFace;
                        }
                        else
                        {
                            if ( facePermutation == facePermutationOwnerDof )
                                localDofOwnerDof += locDofInTopologicalFace;
                            else
                            {
                                uint16_type locFaceDof = invalid_v<uint16_type>;
                                if ( facePermutation.value() == face_permutation_type::IDENTITY )
                                    locFaceDof = locDofInTopologicalFace;
                                else
                                    locFaceDof = this->vector_permutation.at(facePermutation)[locDofInTopologicalFace];

                                if ( facePermutationOwnerDof.value() == face_permutation_type::IDENTITY )
                                    localDofOwnerDof += locFaceDof;
                                else
                                    localDofOwnerDof += mapLocalDofFaceToFaceInElement.at(facePermutationOwnerDof)[locFaceDof];
                            }
                        }

                        // TODO : add debug check that verify dof point are identical!
                    }
                }

                if ( eltPidOwnerDof == invalid_v<rank_type> )
                    continue;

                if ( eltPidOwnerDof != myRank )// invalid_rank_type_value && pidDofActive < myRank )
                {
                    std::vector<std::pair<uint16_type,size_type> > compglobdofs;
                    for ( auto const& localDofAllComp : localDofDatas.second )
                    {
                        uint16_type comp = localDofAllComp.first;
                        const size_type theglobdof = std::get<1>( localDofAllComp.second );
                        const uint16_type thelocdof = std::get<0>( localDofAllComp.second );
                        dofIsGhost[theglobdof] = true;
                        //compglobdofs.push_back( std::make_pair(thelocdof,theglobdof) );
                        ++nDofNotPresent;
#if 1
                        dataToSend[eltPidOwnerDof][eltIdOwnerDof].push_back( this->localDofId( localDofOwnerDof, comp ) );
                        dataMemory[eltPidOwnerDof][eltIdOwnerDof].push_back( theglobdof );
#endif
                    }
                }

            }

        }
    }

    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//

    // update datamap info
    CHECK( this->M_n_localWithGhost_df[myRank] >= nDofNotPresent ) << "invalid data";
    this->M_n_localWithoutGhost_df[myRank] = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;

    // std::vector<boost::tuple<size_type,size_type,size_type> > dataRecvFromGather;
    // auto dataSendToGather = boost::make_tuple(this->M_first_df[myRank],this->M_n_localWithGhost_df[myRank],this->M_n_localWithoutGhost_df[myRank]);
    std::vector<std::tuple<size_type,size_type> > dataRecvFromGather;
    auto dataSendToGather = std::make_tuple(this->M_n_localWithGhost_df[myRank],this->M_n_localWithoutGhost_df[myRank]);
    mpi::all_gather( this->worldComm(),
                     dataSendToGather,
                     dataRecvFromGather );

    for (int p=0;p<this->worldComm().localSize();++p)
    {
        this->M_n_localWithGhost_df[p] = std::get<0>( dataRecvFromGather[p] );
        this->M_n_localWithoutGhost_df[p] = std::get<1>( dataRecvFromGather[p] );
    }
    // update global nDof
    this->M_n_dofs=0;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        this->M_n_dofs+=this->M_n_localWithoutGhost_df[proc];
    }

    this->M_first_df_globalcluster[0]=0;//this->M_first_df[0];
    if ( this->M_n_localWithoutGhost_df[0] > 0 )
        this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0]+this->M_n_localWithoutGhost_df[0]-1;
    else
        this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0];

    for ( int i=1; i<this->worldComm().size(); ++i )
    {
        if ( this->M_n_localWithoutGhost_df[i-1] >0 )
            this->M_first_df_globalcluster[i]=this->M_last_df_globalcluster[i-1]+1;
        else
            this->M_first_df_globalcluster[i]=this->M_last_df_globalcluster[i-1];

        if ( this->M_n_localWithoutGhost_df[i] >0 )
            this->M_last_df_globalcluster[i]=this->M_first_df_globalcluster[i]+this->M_n_localWithoutGhost_df[i]-1;
        else
            this->M_last_df_globalcluster[i]=this->M_first_df_globalcluster[i];
    }
    //------------------------------------------------------------------------------//
    // init map
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );
    //------------------------------------------------------------------------------//
    // add in map the dofs presents
    size_type firstGlobIndex = this->M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( size_type i=0; i< this->M_n_localWithGhost_df[myRank]; ++i )
    {
        if ( !dofIsGhost[i] )
        {
            this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
            ++nextGlobIndex;
        }
    }
   //------------------------------------------------------------------------------//

    // update parallel mapping
#if 0
    this->buildGlobalProcessToGlobalClusterDofMapOthersMeshNonBlockingComm( mesh,listToSend );
#else
    this->buildGlobalProcessToGlobalClusterInterprocessDofs( mesh, dataToSend, dataMemory );
#endif

   //------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------//

    // extended dof table
    if ( this->hasDofTableExtended() )
    {
        //this->buildGhostDofMapExtended( mesh );
        if ( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
            this->buildGhostDofMapExtended( mesh, elements(this->meshSupport(),entity_process_t::GHOST_ONLY ) );
        else
            this->buildGhostDofMapExtended( mesh, elements(mesh,entity_process_t::GHOST_ONLY ) );
    }
}



template<typename MeshType, typename FEType, typename PeriodicityType,typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType,MortarType>::buildGlobalProcessToGlobalClusterInterprocessDofs( mesh_type& mesh,
                                                                                                           std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > & dataToSend,
                                                                                                           std::map<rank_type, std::map<size_type,std::vector<size_type> > > & dataMemory )
{
    std::map<rank_type, std::map<size_type,std::vector<uint16_type> > > dataToRecv;


    int nbMaxRequest = 2*mesh.neighborSubdomains().size();
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;

    // step 1 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank] );
    }
    // step 1 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    //------------------------------------------------------------------------------//
    // step 2 : treat recv and prepare data for response
    std::vector<size_type> tmpIndicesSet;
    std::map<rank_type, std::map<size_type,std::vector<size_type> > > dataToSendStep2, dataToRecvStep2;
    for ( auto const& [rankRecv,eltIdToLocalDofs] : dataToRecv )
    {
        auto & dataToSendStep2AtRank = dataToSendStep2[rankRecv];
        for ( auto const& [eltId,localDofIds] : eltIdToLocalDofs )
        {
            tmpIndicesSet.resize( this->getIndicesSize( eltId ) );
            this->getIndicesSet( eltId, tmpIndicesSet );

            //this->getIndicesSetOnGlobalCluster( eltId, tmpIndicesSetOnGlobalCluster );
            auto & dataToSendStep2AtElt = dataToSendStep2AtRank[eltId];
            dataToSendStep2AtElt.resize( localDofIds.size(), invalid_v<size_type> );
            //for ( size_type localDofIds : localDofIds )
            for (int k=0;k<localDofIds.size();++k)
            {
                size_type dofIdGlobalProcess = tmpIndicesSet.at( localDofIds[k] );
                dataToSendStep2AtElt[k] = this->mapGlobalProcessToGlobalCluster()[ dofIdGlobalProcess ];
#if 0
                auto ababab = this->localDof( eltId );
                std::cout << fmt::format( "SEND M_mapGlobalProcessToGlobalCluster {} = {} with ldof {} and look size:{} eltId:{}",
                                          dofIdGlobalProcess,  dataToSendStep2AtElt[k], localDofIds[k], std::distance(ababab.first,ababab.second),eltId ) << std::endl;
#endif
                // this->M_activeDofSharedOnCluster[dofIdGlobalProcess].insert(rankRecv);
            }
            // dataToSendStep2AtRank.push_back( this->getIndicesOnGlobalCluster( eltId ) );
            // for ( size_type dofIndex : this->getIndices( eltId ) )
            //     if ( !this->dofGlobalProcessIsGhost( dofIndex ) )
            //         this->M_activeDofSharedOnCluster[dofIndex].insert(rankRecv);
        }
    }
    // step 2 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecvStep2[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSendStep2[neighborRank] );
    }
    // step 2 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    //------------------------------------------------------------------------------//
    // step 3 : treat final recv
    for ( auto const& [rankRecv,eltIdToDofGlobalClusters] : dataToRecvStep2 )
    {
        auto & dataMemoryAtRank = dataMemory[rankRecv];
        CHECK( dataMemoryAtRank.size() == eltIdToDofGlobalClusters.size() ) << "sizes should be equal";
        for ( auto const& [eltId,dofGlobalClusterIds] : eltIdToDofGlobalClusters )
        {
            auto & dataMemoryAtElt = dataMemoryAtRank.at(eltId);
            CHECK( dataMemoryAtElt.size() == dofGlobalClusterIds.size() ) << "sizes should be equal";
            for (int k=0;k<dataMemoryAtElt.size();++k)
            {
                size_type dofGlobalProcessId = dataMemoryAtElt[k];
                size_type dofGlobalClusterId = dofGlobalClusterIds[k];
                CHECK( dofGlobalClusterId != invalid_v<size_type> ) << "invalid global cluster dof id";
                //update data map
                //std::cout << fmt::format( "M_mapGlobalProcessToGlobalCluster {} = {}", dofGlobalProcessId, dofGlobalClusterId ) << std::endl;
                this->M_mapGlobalProcessToGlobalCluster[dofGlobalProcessId] = dofGlobalClusterId;
            }
        }
    }
}





//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostDofMapExtended( mesh_type& mesh,
                                                                                   Range<mesh_type,MESH_ELEMENTS> const& ghostEltRange )
{
    DVLOG(2) << "[buildGhostDofMap] call buildGhostDofMapExtended on rank "<<  this->worldComm().rank();

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    size_type start_next_free_dof = this->M_n_localWithGhost_df[myRank];
    //------------------------------------------------------------------------------//
    // build extended dof table
    size_type next_free_dof = start_next_free_dof;
    DofFromElement<self_type,fe_type> dfe( this, *M_fe );
    for ( auto const& ghostEltWrap : ghostEltRange )
    {
        auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
        // elements doftable
        if ( !this->isElementDone( ghostElt.id() ) )
            dfe.add( ghostEltWrap, next_free_dof, myRank );
    }
    //------------------------------------------------------------------------------//
    // update local datamap
    this->M_nGhostDofAddedInExtendedDofTable = next_free_dof-start_next_free_dof;
    std::vector<size_type> dataRecvFromGather;
    mpi::all_gather( this->worldComm().localComm(),
                     this->M_nGhostDofAddedInExtendedDofTable,
                     dataRecvFromGather );
    for (rank_type p=0;p<nProc;++p)
    {
        this->M_n_localWithGhost_df[p] += dataRecvFromGather[p];
    }
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );

    //------------------------------------------------------------------------------//
    // TODO : maybe we can only apply send/recv from active elt (with idInOtherspartitions.size()>0) to all ghost elements

    std::map< rank_type, std::vector<index_type> > dataToSend, dataToRecv, dataMemory;

    // get dofs in extended part
    for ( auto const& ghostEltWrap : ghostEltRange )
    {
        auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
        const rank_type processIdOfGhost = ghostElt.processId();
        const size_type eltIdOfGhost = ghostElt.id();
        const size_type eltIdInOtherPartOfGhost = ghostElt.idInOthersPartitions(processIdOfGhost);

        dataToSend[processIdOfGhost].push_back( eltIdInOtherPartOfGhost );
        dataMemory[processIdOfGhost].push_back( ghostElt.id() );
    }

    //------------------------------------------------------------------------------//

    int nbMaxRequest = 2*mesh.neighborSubdomains().size();
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;
    // send/recv of data size
    std::map<rank_type,std::size_t> sizeRecv, sizeSended;
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        sizeSended[neighborRank] = dataToSend[neighborRank].size();
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, sizeSended[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;
    // step 1 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        std::size_t nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecv[neighborRank].data(), nRecvData );
        std::size_t nSendData = sizeSended[neighborRank];
        if ( nSendData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSend[neighborRank].data(), nSendData );
    }
    // step 1 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;
    //------------------------------------------------------------------------------//

    // step 2 : treat recv and prepare data for response
    std::map< rank_type, std::vector<std::vector<size_type>> > dataToSendStep2, dataToRecvStep2;
    for ( auto const& [rankRecv,eltIds] : dataToRecv )
    {
        auto & dataToSendStep2AtRank = dataToSendStep2[rankRecv];
        dataToSendStep2AtRank.reserve( eltIds.size() );
        for ( index_type eltId : eltIds )
        {
            dataToSendStep2AtRank.push_back( this->getIndicesOnGlobalCluster( eltId ) );
            DCHECK( dataToSendStep2AtRank.back().size() > 0 ) << "no dof found";
            // for ( size_type dofIndex : this->getIndices( eltId ) )
            //     if ( !this->dofGlobalProcessIsGhost( dofIndex ) )
            //         this->M_activeDofSharedOnCluster[dofIndex].insert(rankRecv);
        }
    }

    // step 2 :send/recv of data
    for ( rank_type neighborRank : mesh.neighborSubdomains() )
    {
        std::size_t nRecvData = sizeSended[neighborRank]; // use inverse size
        dataToRecvStep2[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank , 0, dataToRecvStep2[neighborRank].data(), nRecvData );
        std::size_t nSendData = sizeRecv[neighborRank]; // use inverse size
        DCHECK( dataToSendStep2[neighborRank].size() == nSendData ) << fmt::format("incompatible data size {} vs {}",dataToSendStep2[neighborRank].size(), nSendData);
        if ( nSendData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank , 0, dataToSendStep2[neighborRank].data(), nSendData );
    }
    // step 2 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    for ( auto const& [rankRecv,dofIdInElt] : dataToRecvStep2 )
    {
        if ( dofIdInElt.empty() )
            continue;
        auto & dataMemoryAtRank = dataMemory.at( rankRecv );
        DCHECK( dataMemoryAtRank.size() == dofIdInElt.size() )<< fmt::format("incompatible data size {} vs {}",dataMemoryAtRank.size(), dofIdInElt.size() );
        for (int k=0;k<dofIdInElt.size();++k)
        {
            size_type eltId = dataMemoryAtRank[k];
            auto const& dofGlobalClusterIds = dofIdInElt[k];
            auto dofIndices = this->getIndices( eltId );
            DCHECK( dofGlobalClusterIds.size() == dofIndices.size() ) << fmt::format("incompatible data size {} vs {}",dofGlobalClusterIds.size(), dofIndices.size());
            for (int l=0;l<dofGlobalClusterIds.size();++l)
            {
                size_type dofIndex = dofIndices[l];
                size_type dofGlobalClusterIndex = dofGlobalClusterIds[l];
                DCHECK( dofIndex < this->M_mapGlobalProcessToGlobalCluster.size() )<< fmt::format("incompatible dofIndex {} vs {}",dofIndex, this->M_mapGlobalProcessToGlobalCluster.size() );
                if ( this->M_mapGlobalProcessToGlobalCluster[dofIndex] == invalid_v<size_type> )
                    this->M_mapGlobalProcessToGlobalCluster[dofIndex] = dofGlobalClusterIndex;
                else
                    CHECK( dofGlobalClusterIndex == this->M_mapGlobalProcessToGlobalCluster[dofIndex] ) << fmt::format("dof index should be identical {} vs {} ", dofGlobalClusterIndex, this->M_mapGlobalProcessToGlobalCluster[dofIndex] );

                rank_type activeProcId = this->procOnGlobalCluster( dofGlobalClusterIndex );
                if ( activeProcId != myRank )
                    this->addNeighborSubdomain( activeProcId ); // TODO try outside of loop
                else
                    this->M_activeDofSharedOnCluster[dofIndex].insert(rankRecv);
            }
        }
    }
}

} // namespace Feel

#endif /* FEELPP_DOFTABLE_MPI_HPP */
