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
    wc(mesh)->print(fmt::format("[DofTable::buildGhostDofMap rank={}] starts. hasMeshSupport: {}", rank(mesh), this->hasMeshSupport()), FLAGS_v > 1, FLAGS_v > 0, FLAGS_v > 1 );
    // if ( this->hasMeshSupport() )
    //     this->meshSupport()->updateParallelData();

    if ( true )//!mesh.components().test( MESH_UPDATE_FACES ) && !mesh.components().test( MESH_UPDATE_FACES_MINIMAL ) )
    {
        this->buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh );
    }
    else
    {
#if 0
        if (is_continuous)
        {
            DVLOG(2) << "[buildGhostDofMap] call buildGlobalProcessToGlobalClusterDofMapContinuous() with god rank "<<  this->worldComm().godRank() << "\n";
            buildGlobalProcessToGlobalClusterDofMapContinuous( mesh );
        }
        else
        {
            DVLOG(2) << "[buildGhostDofMap] call buildGlobalProcessToGlobalClusterDofMapDiscontinuous() with rank "<<  this->worldComm().rank() << "\n";
            buildGlobalProcessToGlobalClusterDofMapDiscontinuous();
        }

        if ( this->buildDofTableMPIExtended() )
            this->buildGhostDofMapExtended( mesh );
#endif
    }

#if 0
    DVLOG(2) << "[buildGhostDofMap] call localtoglobalOnCluster() with rank "<<  this->worldComm().rank() << "\n";
    auto it_elt = mesh.beginElementWithProcessId( this->comm().rank() );
    auto en_elt = mesh.endElementWithProcessId( this->comm().rank() );

    for ( ; it_elt != en_elt; ++it_elt )
    {
        size_type elid= it_elt->id();

        for ( int i = 0; i < FEType::nLocalDof; ++i )
        {
            int nc1 = ( is_product?nComponents:1 );

            for ( int c1 =0; c1 < nc1; ++c1 )
            {
                int ind = FEType::nLocalDof*c1+i;
                auto const& dof = localToGlobalOnCluster( elid, i, c1 );

                M_locglobOnCluster_indices[elid][ind] = dof.index();
                M_locglobOnCluster_signs[elid][ind] = dof.sign();
            }
        }
    }
#endif
    DVLOG(2) << "[buildGhostDofMap] finish () with rank "<< this->worldComm().rank() << "\n";

}



//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
#if 0
template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGlobalProcessToGlobalClusterDofMapContinuous( mesh_type& mesh )
{
    //------------------------------------------------------------------------------//
    size_type nbFaceDof = invalid_v<size_type>;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );
    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    if ( nbFaceDof == 0 ) return;
    //------------------------------------------------------------------------------//
    // build GlobalProcessToGlobalClusterDofMap for actif dofs and prepare send for ghost dof
    //std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > listToSend(this->worldComm().size());
    std::vector< std::map<size_type,std::set< std::vector<size_type> > > > listToSend(this->worldComm().size());
    std::set<rank_type> procRecvData;
    this->buildGlobalProcessToGlobalClusterDofMapContinuousActifDof(mesh,listToSend,procRecvData);
    //------------------------------------------------------------------------------//
    // update GlobalProcessToGlobalClusterDofMap for ghost dofs
    if ( false )
        this->buildGlobalProcessToGlobalClusterDofMapContinuousGhostDofBlockingComm(mesh,listToSend,procRecvData);
    else
        this->buildGlobalProcessToGlobalClusterDofMapContinuousGhostDofNonBlockingComm(mesh,listToSend,procRecvData);
    //------------------------------------------------------------------------------//
}
#endif
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

namespace detail {

#if 0
template <typename DofTableType>
boost::tuple<rank_type,size_type >
updateDofOnVertices( DofTableType const& doftable, typename DofTableType::mesh_type::face_type const& theface, const rank_type myIdProcess,
                     const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename DofTableType::mesh_type::element_type const& eltOnProc,
                     const uint16_type locDof, std::set<rank_type> & procRecvData )
{
    typedef typename DofTableType::mesh_type MeshType;
    auto mesh = doftable.mesh();
    rank_type procMin = IdProcessOfGhost;
    size_type idFaceMin = idFaceInPartition;

    uint16_type iFaEl = ( theface.processId() == theface.proc_first() )? theface.pos_first():theface.pos_second();
    //local point number (in element)
    uint16_type iPtEl = MeshType::element_type::fToP( iFaEl, locDof );

    auto const& thept = eltOnProc.point(iPtEl);
    auto const theptId = thept.id();

    bool hasMeshSupportPartial = doftable.hasMeshSupport() && doftable.meshSupport()->isPartialSupport();

    auto itprocghost=thept.elementsGhost().begin();
    auto const enprocghost=thept.elementsGhost().end();
    for ( ; itprocghost!=enprocghost ; ++itprocghost)
    {
        if (procMin>itprocghost->first)
        {
            const rank_type theprocGhost=itprocghost->first;
            bool findFace=false;

            size_type ghostEltId = invalid_v<size_type>;// *itprocghost->second.begin();
            if ( hasMeshSupportPartial )
            {
                for ( size_type _ghostEltId : itprocghost->second )
                {
                    if ( doftable.meshSupport()->hasGhostElement( _ghostEltId ) )
                    {
                        ghostEltId = _ghostEltId;
                        break;
                    }
                }
                if ( ghostEltId == invalid_v<size_type> )
                    continue;
            }
            else
            {
                DCHECK(itprocghost->second.size()>0) << "need to have at least one ghost element\n";
                ghostEltId = *itprocghost->second.begin();
                //auto iteltghost = itprocghost->second.begin();
            }
            auto const& eltGhost = mesh->element(ghostEltId);
            for ( uint16_type f = 0; f < MeshType::element_type::numTopologicalFaces && !findFace; ++f )
            {
                if ( !eltGhost.facePtr(f) )
                    continue;
                auto const& faceOnGhost = eltGhost.face(f);
                for ( uint16_type vv = 0; vv < MeshType::face_type::numVertices && !findFace ; ++vv )
                {
                    if ( faceOnGhost.point(vv).id()==theptId )
                    {
                        procMin=theprocGhost;
                        idFaceMin = faceOnGhost.idInOthersPartitions( procMin );
                        findFace=true;
                    }
                }
            }
            CHECK( findFace ) << "PROBLEM with parallel dof table construction : not find a face contained the point on ghost element\n";
        }
    }

    // if current dof is actif then store set of processId which share the dof
    if ( myIdProcess < procMin )
    {
        itprocghost = thept.elementsGhost().begin();
        for ( ; itprocghost!=enprocghost ; ++itprocghost)
        {
            const rank_type procIdGhost = itprocghost->first;
            if ( hasMeshSupportPartial )
            {
                for ( size_type _ghostEltId : itprocghost->second )
                {
                    if ( doftable.meshSupport()->hasGhostElement( _ghostEltId ) )
                    {
                        procRecvData.insert( procIdGhost );
                        break;
                    }
                }
            }
            else
            {
                procRecvData.insert( procIdGhost );
            }
        }
    }

    return boost::make_tuple(procMin,idFaceMin);
}
#endif
template <typename DofTableType>
boost::tuple<rank_type,size_type >
updateDofOnVertices( DofTableType const& doftable, typename DofTableType::mesh_type::element_type const& theelt, const uint16_type ptIdInElt )
{
    auto const& thept = theelt.point( ptIdInElt );
    if ( thept.numberOfProcGhost() == 0 )
        return boost::make_tuple(theelt.processId(),theelt.id());

    auto mesh = doftable.mesh();
    bool hasMeshSupportPartial = doftable.hasMeshSupport() && doftable.meshSupport()->isPartialSupport();
    //const size_type theptId = thept.id();
    rank_type procMin = theelt.processId();
    size_type IdEltMin = theelt.id();
    for ( auto const& eltGhostPair : thept.elementsGhost() )
    {
        const rank_type theprocGhost = eltGhostPair.first;
        if ( theprocGhost < procMin )
        {
            if ( hasMeshSupportPartial )
            {
                for ( size_type _ghostEltId : eltGhostPair.second )
                {
                    if ( doftable.meshSupport()->hasGhostElement( _ghostEltId ) )
                    {
                        auto const& eltGhost = mesh->element( _ghostEltId );
                        procMin = theprocGhost;
                        IdEltMin = eltGhost.idInOthersPartitions( procMin );
                        break;
                    }
                }
            }
            else
            {
                DCHECK( eltGhostPair.second.size()>0 ) << "need to have at least one ghost element\n";
                auto iteltghost = eltGhostPair.second.begin();
                auto const& eltGhost = mesh->element(*iteltghost);
                procMin = theprocGhost;
                IdEltMin = eltGhost.idInOthersPartitions( procMin );
            }
#if 0
            bool findDofVertice = false;
            DCHECK( eltGhostPair.second.size()>0 ) << "need to have at least one ghost element\n";
            auto iteltghost = eltGhostPair.second.begin();
            auto const& eltGhost = mesh.element(*iteltghost);
            for ( uint16_type n=0; n < eltGhost.nVertices(); n++ )
            {
                if ( eltGhost.point(n).id() == theptId )
                {
                    procMin=theprocGhost;
                    IdEltMin = eltGhost.idInOthersPartitions( procMin );
                    findDofVertice=true;
                    break;
                }
            }
            CHECK( findDofVertice ) << "PROBLEM with parallel dof table construction : not find a vertice in ghost element associated to a dof point at interprocess\n";
#endif
        }
    }
    return boost::make_tuple(procMin,IdEltMin);
}
//--------------------------------------------------------------------------------------------------------//

template <typename DofTableType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( DofTableType const& doftable, typename DofTableType::mesh_type::face_type const& theface, const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename DofTableType::mesh_type::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData,
                  typename std::enable_if< mpl::not_<is_3d<typename DofTableType::mesh_type>>::value >::type* = nullptr )
{
    return boost::make_tuple(0,0);
}
template <typename DofTableType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( DofTableType const& doftable, typename DofTableType::mesh_type::face_type const& theface, const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename DofTableType::mesh_type::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData,
                  typename std::enable_if< is_3d<typename DofTableType::mesh_type>::value >::type* = nullptr )
{
    typedef typename DofTableType::mesh_type MeshType;
    auto mesh = doftable.mesh();
    rank_type procMin = IdProcessOfGhost;
    size_type idFaceMin = idFaceInPartition;

    uint16_type iFaEl = ( theface.processId() == theface.proc_first() )? theface.pos_first():theface.pos_second();
    //local edge number (in element)
    uint16_type iEdEl = MeshType::element_type::fToE(  iFaEl, idEdgesInFace );

    auto const& theedge = eltOnProc.edge(iEdEl);
    auto const theedgeId = theedge.id();

    bool hasMeshSupportPartial = doftable.hasMeshSupport() && doftable.meshSupport()->isPartialSupport();

    auto itprocghost=theedge.elementsGhost().begin();
    auto const enprocghost=theedge.elementsGhost().end();
    for ( ; itprocghost!=enprocghost ; ++itprocghost)
    {
        if ( procMin>itprocghost->first )
        {
            const rank_type theprocGhost=itprocghost->first;

            size_type ghostEltId = invalid_v<size_type>;
            if ( hasMeshSupportPartial )
            {
                for ( size_type _ghostEltId : itprocghost->second )
                {
                    if ( doftable.meshSupport()->hasGhostElement( _ghostEltId ) )
                    {
                        ghostEltId = _ghostEltId;
                        break;
                    }
                }
                if ( ghostEltId == invalid_v<size_type> )
                    continue;
            }
            else
            {
                DCHECK(itprocghost->second.size()>0) << "need to have at least one ghost element\n";
                ghostEltId = *itprocghost->second.begin();
            }
            auto const& eltGhost = mesh->element(ghostEltId);

            bool findFace=false;
            for ( uint16_type f = 0; f < MeshType::element_type::numTopologicalFaces && !findFace; ++f )
            {
                if ( !eltGhost.facePtr(f) )
                    continue;
                auto const& faceOnGhost = eltGhost.face(f);
                for ( uint16_type vv = 0; vv < MeshType::face_type::numEdges && !findFace ; ++vv )
                {
                    if (faceOnGhost.edge(vv).id()==theedgeId)
                    {
                        procMin=theprocGhost;
                        findFace=true;
                        idFaceMin = faceOnGhost.idInOthersPartitions(procMin);
                    }
                }
            }
            CHECK( findFace ) << "\nPROBLEM with parallel dof table construction \n";
        }
    }

    // if current dof is actif then store set of processId which share the dof
    if ( myIdProcess < procMin )
    {
        itprocghost=theedge.elementsGhost().begin();
        for ( ; itprocghost!=enprocghost ; ++itprocghost)
        {
            const rank_type procIdGhost = itprocghost->first;
            if ( hasMeshSupportPartial )
            {
                for ( size_type _ghostEltId : itprocghost->second )
                {
                    if ( doftable.meshSupport()->hasGhostElement( _ghostEltId ) )
                    {
                        procRecvData.insert( procIdGhost );
                        break;
                    }
                }
            }
            else
            {
                procRecvData.insert( procIdGhost );
            }
        }

    }

    return boost::make_tuple(procMin,idFaceMin);
}


template <typename DofTableType>
boost::tuple<rank_type,size_type >
updateDofOnEdges( DofTableType const& doftable, typename DofTableType::mesh_type::element_type const& theelt, const uint16_type edgeIdInElt,
                  typename std::enable_if< mpl::not_<is_3d<typename DofTableType::mesh_type>>::value >::type* = nullptr )
{
    return boost::make_tuple( invalid_rank_type_value,invalid_v<size_type> );
}
template <typename DofTableType>
boost::tuple<rank_type,size_type >
updateDofOnEdges( DofTableType const& doftable, typename DofTableType::mesh_type::element_type const& theelt, const uint16_type edgeIdInElt,
                  typename std::enable_if< is_3d<typename DofTableType::mesh_type>::value >::type* = nullptr )
{
    auto const& theedge = theelt.edge( edgeIdInElt );
    if ( theedge.numberOfProcGhost() == 0 )
        return boost::make_tuple(theelt.processId(),theelt.id());

    auto mesh = doftable.mesh();
    bool hasMeshSupportPartial = doftable.hasMeshSupport() && doftable.meshSupport()->isPartialSupport();
    rank_type procMin = theelt.processId();
    size_type IdEltMin = theelt.id();
    for ( auto const& eltGhostPair : theedge.elementsGhost() )
    {
        const rank_type theprocGhost = eltGhostPair.first;
        if ( theprocGhost < procMin )
        {
            if ( hasMeshSupportPartial )
            {
                for ( size_type _ghostEltId : eltGhostPair.second )
                {
                    if ( doftable.meshSupport()->hasGhostElement( _ghostEltId ) )
                    {
                        auto const& eltGhost = mesh->element( _ghostEltId );
                        procMin = theprocGhost;
                        IdEltMin = eltGhost.idInOthersPartitions( procMin );
                        break;
                    }
                }
            }
            else
            {
                DCHECK( eltGhostPair.second.size()>0 ) << "need to have at least one ghost element\n";
                auto iteltghost = eltGhostPair.second.begin();
                auto const& eltGhost = mesh->element(*iteltghost);
                procMin = theprocGhost;
                IdEltMin = eltGhost.idInOthersPartitions( procMin );
            }
        }
    }
    return boost::make_tuple(procMin,IdEltMin);
}

} // namespace detail

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
#if 0
template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType,MortarType>::buildGlobalProcessToGlobalClusterDofMapContinuousActifDof( mesh_type& mesh,
                                                                                                        std::vector< std::map<size_type,std::set<std::vector<size_type> > > > & listToSend,
                                                                                                        std::set<rank_type> & procRecvData )
{
    // goal init container listToSend
    // std::vector< std::map<size_type,std::set<size_type> > > ( proc,( idFace,(globDof,..)), ...   )) )
    const rank_type myRank = this->worldComm().rank();
    //------------------------------------------------------------------------------//
    // get nbFaceDof
    size_type nbFaceDof = invalid_v<size_type>;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );
    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    DVLOG(2) << "[buildGhostInterProcessDofMap] nbFaceDof " << nbFaceDof << "\n";

    if ( nbFaceDof == 0 ) return;

    const uint16_type ncdof = is_product?nComponents:1;
    DVLOG(2) << "[buildGhostInterProcessDofMap] ncdof " << ncdof << "\n";

    //------------------------------------------------------------------------------//

    std::vector<bool> dofdone(this->M_n_localWithGhost_df[myRank],false);
    std::vector<bool> dofIsGhost(this->M_n_localWithGhost_df[myRank],false);
    size_type nDofNotPresent=0;

    // iteration on all interprocessfaces in order to send requests to the near proc
    auto rangeInterProcessFaces = (this->hasMeshSupport())? interprocessfaces(this->meshSupport()) : interprocessfaces(mesh);
    for ( auto const& faceipWrap : rangeInterProcessFaces )
    {
        auto const& faceip = unwrap_ref(faceipWrap);
        DVLOG(2) << "[buildGhostInterProcessDofMap] face id: " << faceip.id() << "\n";
        auto const& elt0 = faceip.element0();
        auto const& elt1 = faceip.element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") eltOnProc id: "  << eltOnProc.id()  << "G(): " << eltOnProc.G() << "\n";
        DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") eltOffProc id: " << eltOffProc.id() << "G(): " << eltOffProc.G() << "\n";

        //------------------------------------------------------------------------------//

        const rank_type IdProcessOfGhostIP = eltOffProc.processId();
        const size_type idFaceInPartitionIP = faceip.idInOthersPartitions( IdProcessOfGhostIP );
        rank_type IdProcessOfGhost = IdProcessOfGhostIP;
        size_type idFaceInPartition = idFaceInPartitionIP;

        //------------------------------------------------------------------------------//
        // for each dof in face
        for ( uint16_type locDof = 0; locDof < nbFaceDof; ++locDof )
        {
            // check only component 0
            DCHECK( M_face_l2g.find( faceip.id() ) != M_face_l2g.end() ) << "not found the face id "<< faceip << "into the mapping faceLocalToGlobal";
            const size_type theglobdoftest = faceLocalToGlobal( faceip.id(),locDof, 0 ).index();
            CHECK( theglobdoftest < this->M_n_localWithGhost_df[myRank] ) << "invalid globdof " << theglobdoftest << "\n";
            if ( dofdone[theglobdoftest] ) continue;

            IdProcessOfGhost = IdProcessOfGhostIP;
            idFaceInPartition = idFaceInPartitionIP;
            if ( locDof < face_type::numVertices*fe_type::nDofPerVertex)
            {
                const int nDofPerVertexTemp = mpl::if_<boost::is_same<mpl::int_<fe_type::nDofPerVertex>,mpl::int_<0> >,
                                                       mpl::int_<1>,
                                                       mpl::int_<fe_type::nDofPerVertex> >::type::value;
                int pointGetLocDof = locDof / nDofPerVertexTemp;

                boost::tie( IdProcessOfGhost, idFaceInPartition ) = Feel::detail::updateDofOnVertices( *this, faceip, myRank, IdProcessOfGhost, idFaceInPartition, eltOnProc, pointGetLocDof,
                                                                                                       procRecvData );
            }
            else if ( nDim == 3 && locDof < (face_type::numVertices*fe_type::nDofPerVertex + face_type::numEdges*fe_type::nDofPerEdge) )
            {
                int locDofInEgde = locDof - face_type::numVertices*fe_type::nDofPerVertex;
                const int nDofPerEdgeTemp = mpl::if_<boost::is_same<mpl::int_<fe_type::nDofPerEdge>,mpl::int_<0> >,
                                                     mpl::int_<1>,
                                                     mpl::int_<fe_type::nDofPerEdge> >::type::value;
                int edgeGetLocDof = locDofInEgde / nDofPerEdgeTemp;

                boost::tie( IdProcessOfGhost, idFaceInPartition ) = Feel::detail::updateDofOnEdges( *this, faceip, myRank, IdProcessOfGhost, idFaceInPartition, eltOnProc, edgeGetLocDof,
                                                                                                    procRecvData );
            }
            else
            {
                if ( myRank < IdProcessOfGhost ) procRecvData.insert( IdProcessOfGhost );
            }

            // if dof is ghost -> prepare send/recv
            if (IdProcessOfGhost<myRank)
            {
                std::vector<size_type > compglobdofs( ncdof );
                for ( uint16_type c = 0; c < ncdof; ++c )
                {
                    // add dof in subcontainer
                    const size_type theglobdof = faceLocalToGlobal( faceip.id(),locDof,c ).index();
                    dofIsGhost[theglobdof] = true;
                    compglobdofs[c]=theglobdof;
                    //listToSend[IdProcessOfGhost][idFaceInPartition].insert(boost::make_tuple(theglobdof,c));
                    ++nDofNotPresent;
                }
                listToSend[IdProcessOfGhost][idFaceInPartition].insert( compglobdofs );
            }

            dofdone[theglobdoftest]=true;

        } // for ( uint16_type locDof = 0; locDof < nbFaceDof; ++locDof )

    } // for ( ; face_it != face_en ; ++face_it )

    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    // update datamap info
    CHECK( this->M_n_localWithGhost_df[myRank] >= nDofNotPresent ) << "invalid data\n" << std::endl;
    //const size_type mynDofWithoutGhost = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;
    this->M_n_localWithoutGhost_df[myRank] = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;
#if 0
    mpi::all_gather( this->worldComm(),
                     mynDofWithoutGhost,
                     this->M_n_localWithoutGhost_df );
#else
    std::vector<boost::tuple<size_type,size_type,size_type> > dataRecvFromGather;
    auto dataSendToGather = boost::make_tuple(this->M_first_df[myRank],this->M_n_localWithGhost_df[myRank],this->M_n_localWithoutGhost_df[myRank]);
    mpi::all_gather( this->worldComm(),
                     dataSendToGather,
                     dataRecvFromGather );

    for (int p=0;p<this->worldComm().localSize();++p)
    {
        this->M_first_df[p] = dataRecvFromGather[p].template get<0>();
        this->M_n_localWithGhost_df[p] = dataRecvFromGather[p].template get<1>();
        this->M_last_df[p] = (this->M_n_localWithGhost_df[p] > 0)? this->M_first_df[p] + this->M_n_localWithGhost_df[p] - 1 : this->M_first_df[p];
        this->M_n_localWithoutGhost_df[p] = dataRecvFromGather[p].template get<2>();
    }
#endif

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
} // buildGlobalProcessToGlobalClusterDofMapContinuousActifDof

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType,MortarType>::
buildGlobalProcessToGlobalClusterDofMapContinuousGhostDofBlockingComm( mesh_type& mesh,
                                                                       std::vector< std::map<size_type,std::set<std::vector<size_type> > > > const& listToSend,
                                                                       std::set<rank_type> const& procRecvData )
{
    const int myRank = this->worldComm().rank();
    //--------------------------------------------------------------------------------------------------------//
    size_type nbFaceDof = invalid_v<size_type>;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    const uint16_type ncdof = is_product?nComponents:1;
    //--------------------------------------------------------------------------------------------------------//

    std::vector<int> nbMsgToSend( this->worldComm().size(), 0 );

    std::vector< std::vector< std::vector<size_type> > > memoryInitialRequest( this->worldComm().size() );

    typedef std::vector< boost::tuple<uint16_type, ublas::vector<double> > > dofs_in_face_subcontainer_type;
    typedef boost::tuple<size_type, dofs_in_face_subcontainer_type > dofs_in_face_container_type;

    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        auto itFaces = listToSend[proc].begin();
        auto const enFaces = listToSend[proc].end();
        const int nFaceToSend = std::distance(itFaces,enFaces);
        memoryInitialRequest[proc].resize(nFaceToSend);
        if ( nFaceToSend>0 )
        {
            this->worldComm().send( proc , 0, nFaceToSend );
            ++nbMsgToSend[proc];
        }

        for ( int cptFaces=0 ; itFaces!=enFaces ; ++itFaces, ++cptFaces)
        {
            auto itDof = itFaces->second.begin();
            auto const enDof = itFaces->second.end();
            const int nDofsInFace = std::distance(itDof,enDof);
            CHECK( nDofsInFace>0 ) << "error in data to send : nDofsInFace=" << nDofsInFace<<" must be > 0 \n";

            dofs_in_face_subcontainer_type dofsInFaceContainer(nDofsInFace);
            memoryInitialRequest[proc][cptFaces].resize(nDofsInFace);
            for (int cptDof=0 ; itDof!=enDof ; ++itDof/*,++cptDof*/)
            {
                for (uint16_type comp=0; comp<ncdof ; ++comp,++cptDof)
                {
                    const size_type theglobdof = itDof->operator[](comp);
                    //auto const theglobdof = itDof->get<0>();
                    //auto const comp = itDof->get<1>();
                    // save the tag of mpi send
                    memoryInitialRequest[proc][cptFaces][cptDof] = theglobdof;
                    //------------------------------------------------------------------------------//
                    // get info to send
                    ublas::vector<double> nodeDofToSend( nRealDim );
                    auto itFindDofPoint = M_dof_points.find( theglobdof );
                    CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                    nodeDofToSend[0]=itFindDofPoint->second.template get<0>()[0];
                    if ( nRealDim>1 )
                        nodeDofToSend[1]=itFindDofPoint->second.template get<0>()[1];
                    if ( nRealDim>2 )
                        nodeDofToSend[2]=itFindDofPoint->second.template get<0>()[2];
                    // up container
                    dofsInFaceContainer[cptDof] = boost::make_tuple(comp,nodeDofToSend);
                    //------------------------------------------------------------------------------//
                }
            }

            this->worldComm().send( proc , nbMsgToSend[proc], boost::make_tuple(itFaces->first,dofsInFaceContainer) );
            ++nbMsgToSend[proc];
        } // for ( int cptFaces=0 ; itFaces!=enFaces ; ++itFaces, ++cptFaces)

    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )

    //--------------------------------------------------------------------------------------------------------//
#if 0
    // counter of msg received for each process
    std::vector<int> nbMsgToRecv;
    mpi::all_to_all( this->worldComm(),
                     nbMsgToSend,
                     nbMsgToRecv );

    for ( int proc=0; proc<this->worldComm().size(); ++proc )
        {
            //CHECK( nbMsgToRecv[proc]==nbMsgToRecv2[proc] )
            if (nbMsgToRecv[proc]!=nbMsgToRecv2[proc] /*|| true*/  ) std::cout
                                                            << "partitioning data incorect "
                                                            << "myrank " << this->worldComm().localRank() << " proc " << proc
                                                            << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
                                                            << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc]
                                                            << " nbMsgToSend[proc] " << nbMsgToSend[proc]
                                                            << "\n";
        }
#endif

    //--------------------------------------------------------------------------------------------------------//
    // recv dof asked and re-send
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        if ( procRecvData.find(proc) == procRecvData.end() ) continue;
        int nbDataRecv=0;
        this->worldComm().recv( proc, 0, nbDataRecv );
        for ( int cpt=1; cpt<nbDataRecv+1; ++cpt )
        {
            dofs_in_face_container_type dataToRecvVec;
            this->worldComm().recv( proc, cpt, dataToRecvVec );

            auto const idFaceInMyPartition = dataToRecvVec.template get<0>();
            DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                    << "idFaceInMyPartition: " << idFaceInMyPartition << "\n";

            auto itDofInFace = dataToRecvVec.template get<1>().begin();
            auto const enDofInFace = dataToRecvVec.template get<1>().end();
            std::vector< size_type > resAskedWithMultiProcess(std::distance(itDofInFace,enDofInFace));
            for ( int cptDofInFace=0 ; itDofInFace != enDofInFace ; ++itDofInFace,++cptDofInFace )
            {
                auto const comp = itDofInFace->template get<0>();
                auto const nodeDofRecv = itDofInFace->template get<1>();

                //------------------------------------------------------------------------------//

                auto const& theface = mesh.face( idFaceInMyPartition );
                auto const& elt0 = theface.element0();
                auto const& elt1 = theface.element1();

                const bool elt0isGhost = elt0.isGhostCell();
                auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
                auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

                //------------------------------------------------------------------------------//
                // search dof on face recv
                int locDof = nbFaceDof;
                bool find=false;

                if ( false && nDim==1 )
                {
                    auto itdofpt = this->dofPointBegin();
                    auto const endofpt = this->dofPointEnd();
                    for ( ; itdofpt!=endofpt && !find ; ++itdofpt )
                    {
                        const auto thedofPt = itdofpt->second.template get<0>();
                        if ( itdofpt->second.template get<2>() != comp ) continue;

                        DVLOG(3) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                                 << "thedofPt: " << thedofPt << "nodeDofRecv: " << nodeDofRecv << "\n";

                        // test equatlity of dofs point
                        bool find2=true;
                        for (uint16_type d=0;d<nRealDim;++d)
                        {
                            find2 = find2 && (std::abs( thedofPt[d]-nodeDofRecv[d] )<1e-9);
                        }
                        // if find else save local dof
                        if (find2) { locDof = itdofpt->second.template get<1>();find=true; }
                    }
                    // check
                    CHECK( find ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";
                    //------------------------------------------------------------------------------//
                    // get global dof
                    const auto dofGlobAsked = locDof;
                    // save response
                    resAskedWithMultiProcess[cptDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                    this->M_activeDofSharedOnCluster[dofGlobAsked].insert(proc);
                }
                else
                {
                for ( uint16_type l = 0; ( l < nbFaceDof && !find ) ; ++l )
                {
                    // dof point in face
                    auto itFindDofPoint = M_dof_points.find( faceLocalToGlobal( idFaceInMyPartition, l, comp ).index() );
                    CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                    auto const& thedofPtInFace = itFindDofPoint->second.template get<0>();
                    DVLOG(3) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                            << "thedofPtInFace: " << thedofPtInFace << "nodeDofRecv: " << nodeDofRecv << "\n";

                    // test equatlity of dofs point
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                    {
                        find2 = find2 && (std::abs( thedofPtInFace[d]-nodeDofRecv[d] )<1e-9);
                    }
                    // if find else save local dof
                    if (find2)
                    {
                        locDof = l;
                        find=true;
                    }
                } // for ( uint16_type l = 0; ( l < nbFaceDof && !find ) ; ++l )
                //------------------------------------------------------------------------------//
                // check
                CHECK( find ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";
                //------------------------------------------------------------------------------//
                // get global dof
                const auto thedof = faceLocalToGlobal( idFaceInMyPartition, locDof, comp );
                const auto dofGlobAsked = thedof.index();
                // save response
                resAskedWithMultiProcess[cptDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                this->M_activeDofSharedOnCluster[dofGlobAsked].insert(proc);
                }
                //------------------------------------------------------------------------------//
            }
            this->worldComm().send( proc, cpt, resAskedWithMultiProcess );

        } // for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )

    //--------------------------------------------------------------------------------------------------------//
    // get response to initial request and update Feel::Mesh::Faces data
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=1; cpt<nbMsgToSend[proc]; ++cpt )
        {
            //------------------------------------------------------------------------------//
            // recv response
            std::vector< size_type > resultRecvWithMultiProcess;
            this->worldComm().recv( proc, cpt, resultRecvWithMultiProcess );
            //------------------------------------------------------------------------------//
            // iterate on dofs
            auto itDofRes = resultRecvWithMultiProcess.begin();
            auto const enDofRes = resultRecvWithMultiProcess.end();
            for ( int cptDofRes=0 ; itDofRes != enDofRes ; ++itDofRes,++cptDofRes )
            {
                auto const myGlobProcessDof = memoryInitialRequest[proc][cpt-1][cptDofRes];
                const auto dofGlobRecv = *itDofRes;
                //update data map
                this->M_mapGlobalProcessToGlobalCluster[myGlobProcessDof]=dofGlobRecv;
            }
        } // for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )

} // buildGlobalProcessToGlobalClusterDofMapContinuousGhostDof

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType,typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType,MortarType>::buildGlobalProcessToGlobalClusterDofMapContinuousGhostDofNonBlockingComm( mesh_type& mesh,
                                                                            std::vector< std::map<size_type,std::set<std::vector<size_type> > > > const& listToSend,
                                                                                                        std::set<rank_type> const& procRecvData )
{
    typedef std::vector< boost::tuple<uint16_type, ublas::vector<double> > > dofs_in_face_subcontainer_type;
    typedef boost::tuple<size_type, dofs_in_face_subcontainer_type > dofs_in_face_container_type;
    typedef std::vector< dofs_in_face_container_type > dofs_container_to_send_type;

    const int myRank = this->worldComm().localRank();
    const int nProc = this->worldComm().localSize();
    //--------------------------------------------------------------------------------------------------------//
    size_type nbFaceDof = invalid_v<size_type>;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    const uint16_type ncdof = is_product?nComponents:1;
    const bool componentsAreSamePoint=true;
    //--------------------------------------------------------------------------------------------------------//
    // compute size of container to send
    std::map< rank_type, int > nDataInVecToSend;
    for ( rank_type proc=0; proc<this->worldComm().size(); ++proc )
    {
        const int nFaceToSend = listToSend[proc].size();
        if ( nFaceToSend == 0 ) continue;
        nDataInVecToSend[proc] = nFaceToSend;
    }
    //--------------------------------------------------------------------------------------------------------//
    // init and resize the container to send
    std::map< rank_type,  dofs_container_to_send_type> dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //--------------------------------------------------------------------------------------------------------//
    // prepare container to send
    std::map< rank_type, std::vector< std::vector<size_type> > > memoryInitialRequest;
    std::map< rank_type, int > nDataInVecToSendBis;
    for ( rank_type proc=0; proc<this->worldComm().size(); ++proc )
    {
        if ( listToSend[proc].size() == 0 ) continue;

        auto itFaces = listToSend[proc].begin();
        auto const enFaces = listToSend[proc].end();
        const int nFaceToSend = std::distance(itFaces,enFaces);
        memoryInitialRequest[proc].resize(nFaceToSend);

        for ( int cptFaces=0 ; itFaces!=enFaces ; ++itFaces, ++cptFaces)
        {
            auto itDof = itFaces->second.begin();
            auto const enDof = itFaces->second.end();
            const int nDofsInFace = std::distance(itDof,enDof)*ncdof;
            const int nDofsInFaceForComm = (componentsAreSamePoint)?std::distance(itDof,enDof) : nDofsInFace;

            CHECK( nDofsInFace>0 ) << "error in data to send : nDofsInFace=" << nDofsInFace<<" must be > 0 \n";

            dofs_in_face_subcontainer_type dofsInFaceContainer(nDofsInFaceForComm);
            memoryInitialRequest[proc][cptFaces].resize(nDofsInFace);
            for (int cptDof=0, cptDof2=0 ; itDof!=enDof ; ++itDof,++cptDof2)
            {
                for (uint16_type comp=0; comp<ncdof ; ++comp,++cptDof)
                {
                    const size_type theglobdof = itDof->operator[](comp);

                    //auto const theglobdof = itDof->get<0>();
                    //auto const comp = itDof->get<1>();
                    // save the tag of mpi send
                    const int indexDof = (componentsAreSamePoint)? comp*nDofsInFaceForComm + cptDof2 : cptDof;
                    memoryInitialRequest[proc][cptFaces][indexDof/*cptDof*/] = theglobdof;
                    //------------------------------------------------------------------------------//
                    if (!componentsAreSamePoint)
                    {
                        // get info to send
                        ublas::vector<double> nodeDofToSend( nRealDim );
                        auto itFindDofPoint = M_dof_points.find( theglobdof );
                        CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                        nodeDofToSend[0]=itFindDofPoint->second.template get<0>()[0];
                        if ( nRealDim>1 )
                            nodeDofToSend[1]=itFindDofPoint->second.template get<0>()[1];
                        if ( nRealDim>2 )
                            nodeDofToSend[2]=itFindDofPoint->second.template get<0>()[2];
                        // up container
                        dofsInFaceContainer[cptDof] = boost::make_tuple(comp,nodeDofToSend);
                    }
                    else if (comp==0)
                    {
                        // get info to send
                        ublas::vector<double> nodeDofToSend( nRealDim );
                        auto itFindDofPoint = M_dof_points.find( theglobdof );
                        CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                        nodeDofToSend[0]=itFindDofPoint->second.template get<0>()[0];
                        if ( nRealDim>1 )
                            nodeDofToSend[1]=itFindDofPoint->second.template get<0>()[1];
                        if ( nRealDim>2 )
                            nodeDofToSend[2]=itFindDofPoint->second.template get<0>()[2];
                        // up container
                        dofsInFaceContainer[cptDof2] = boost::make_tuple(0,nodeDofToSend);
                    }

                    //------------------------------------------------------------------------------//
                }
            }

            if ( nDataInVecToSendBis.find(proc) == nDataInVecToSendBis.end() )
                nDataInVecToSendBis[proc]=0;
            // update container
            dataToSend[proc][nDataInVecToSendBis[proc]] = boost::make_tuple(itFaces->first,dofsInFaceContainer);
            // update counter
            nDataInVecToSendBis[proc]++;
        } // for ( int cptFaces=0 ; itFaces!=enFaces ; ++itFaces, ++cptFaces)
    }

    //--------------------------------------------------------------------------------------------------------//
    // counter of request
    int nbRequest=0;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( dataToSend.find(proc) != dataToSend.end() )
            ++nbRequest;
        if ( procRecvData.find(proc) != procRecvData.end() )
            ++nbRequest;
    }
    if ( nbRequest ==0 ) return;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    //--------------------------------------------------------------------------------------------------------//
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToSend->first , 0, itDataToSend->second );
        ++cptRequest;
    }
    //--------------------------------------------------------------------------------------------------------//
    // first recv
    std::map<rank_type,dofs_container_to_send_type> dataToRecv;
    auto itProcRecvData = procRecvData.begin();
    auto const enProcRecvData = procRecvData.end();
    for ( ; itProcRecvData!=enProcRecvData ; ++itProcRecvData )
    {
        const rank_type proc = *itProcRecvData;
        reqs[cptRequest] = this->worldComm().localComm().irecv( proc , 0, dataToRecv[proc] );
        ++cptRequest;
    }
    //--------------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //--------------------------------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector< std::vector<size_type> > > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        auto itFaceRecv = itDataRecv->second.begin();
        auto const enFaceRecv = itDataRecv->second.end();
        const int nFaceRecv=  itDataRecv->second.size();
        dataToReSend[idProc].resize( nFaceRecv );
        for ( int cptFace=0 ; itFaceRecv!=enFaceRecv ; ++itFaceRecv,++cptFace )
        {
            auto const idFaceInMyPartition = itFaceRecv->template get<0>();
            DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                    << "idFaceInMyPartition: " << idFaceInMyPartition << "\n";
            auto const& theface = mesh.face( idFaceInMyPartition );
            auto const& elt0 = theface.element0();
            auto const& elt1 = theface.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            auto itDofInFace = itFaceRecv->template get<1>().begin();
            auto const enDofInFace = itFaceRecv->template get<1>().end();
            const int nDofInFace = distance(itDofInFace,enDofInFace);
            const int nDofInFaceOpt = (componentsAreSamePoint)? nDofInFace*ncdof : nDofInFace;
            dataToReSend[idProc][cptFace].resize( nDofInFaceOpt,invalid_v<size_type> );
            for ( int cptDofInFace=0 ; itDofInFace != enDofInFace ; ++itDofInFace,++cptDofInFace )
            {
                auto const comp = itDofInFace->template get<0>();
                auto const nodeDofRecv = itDofInFace->template get<1>();
                //------------------------------------------------------------------------------//
                // search dof on face recv
                int locDof = nbFaceDof;
                bool find=false;
                for ( uint16_type l = 0 ; l < nbFaceDof && !find ; ++l )
                {
                    // dof point in face
                    auto itFindDofPoint = M_dof_points.find( faceLocalToGlobal( idFaceInMyPartition, l, comp ).index() );
                    CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                    auto const& thedofPtInFace = itFindDofPoint->second.template get<0>();
                    DVLOG(3) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                            << "thedofPtInFace: " << thedofPtInFace << "nodeDofRecv: " << nodeDofRecv << "\n";
                    // test equatlity of dofs point
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                    {
                        find2 = find2 && (std::abs( thedofPtInFace[d]-nodeDofRecv[d] )<1e-9);
                    }
                    // if find else save local dof
                    if (find2)
                    {
                        locDof = l;
                        find=true;
                    }
                } // for ( uint16_type l = 0; ( l < nbFaceDof && !find ) ; ++l )
                //------------------------------------------------------------------------------//
                // check
                CHECK( find ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";
                //------------------------------------------------------------------------------//
                if (!componentsAreSamePoint)
                {
                    // get global dof
                    const auto thedof = faceLocalToGlobal( idFaceInMyPartition, locDof, comp );
                    const auto dofGlobAsked = thedof.index();
                    // save response
                    dataToReSend[idProc][cptFace][cptDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                    this->M_activeDofSharedOnCluster[dofGlobAsked].insert(idProc);
                }
                else
                {
                    for (uint16_type comp2=0; comp2<ncdof ; ++comp2)
                    {
                        const auto thedof = faceLocalToGlobal( idFaceInMyPartition, locDof, comp2 );
                        const size_type dofGlobAsked = thedof.index();
                        const int indexDofInFace = comp2*nDofInFace + cptDofInFace;
                        // save response
                        dataToReSend[idProc][cptFace][indexDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                        this->M_activeDofSharedOnCluster[dofGlobAsked].insert(idProc);
                    }
                }
                //------------------------------------------------------------------------------//
            }
        } // for ( int cptFace=0 ... )
    } // for ( ; itDataRecv ... )

    //--------------------------------------------------------------------------------------------------------//
    // send respond to the request
    cptRequest=0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend!=enDataToReSend ; ++itDataToReSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToReSend->first , 0, itDataToReSend->second );
        ++cptRequest;
    }
    //--------------------------------------------------------------------------------------------------------//
    // recv the initial request
    std::map<rank_type, std::vector<std::vector<size_type> > > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc, 0, finalDataToRecv[idProc] );
        ++cptRequest;
    }
    //--------------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //--------------------------------------------------------------------------------------------------------//
    // update datamap for ghost dof
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const rank_type idProc = itFinalDataToRecv->first;
        auto itFaceRecv = itFinalDataToRecv->second.begin();
        auto const enFaceRecv = itFinalDataToRecv->second.end();
        for ( int cptFace=0 ; itFaceRecv!=enFaceRecv ; ++itFaceRecv,++cptFace )
        {
            if (componentsAreSamePoint)
            {
                const int nDofsInFace = itFaceRecv->size()/ncdof;
                for ( int cptDof=0 ; cptDof< nDofsInFace ; ++cptDof )
                {
                    for (uint16_type comp2=0; comp2<ncdof ; ++comp2)
                    {
                        const int myindexDof = comp2*nDofsInFace + cptDof;
                        const size_type myGlobProcessDof = memoryInitialRequest[idProc][cptFace][myindexDof];
                        const size_type dofGlobRecv = itFaceRecv->operator[](myindexDof);
                        //update data map
                        this->M_mapGlobalProcessToGlobalCluster[myGlobProcessDof] = dofGlobRecv;
                    }
                }
            }
            else
            {
                auto itDofInFace=itFaceRecv->begin();
                auto const enDofInFace=itFaceRecv->end();
                for ( int cptDof=0 ; itDofInFace!=enDofInFace ; ++itDofInFace,++cptDof )
                {
                    const size_type myGlobProcessDof = memoryInitialRequest[idProc][cptFace][cptDof];
                    const size_type dofGlobRecv = *itDofInFace;
                    //update data map
                    this->M_mapGlobalProcessToGlobalCluster[myGlobProcessDof] = dofGlobRecv;
                }
            }
        }
    }


} // buildGlobalProcessToGlobalClusterDofMapContinuousGhostDofNonBlockingComm

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGlobalProcessToGlobalClusterDofMapDiscontinuous()
{
    const rank_type myRank = this->worldComm().rank();
    //------------------------------------------------------------------------------//
    // update datamap info
    this->M_n_dofs=0;
    for ( rank_type proc=0; proc<this->worldComm().size(); ++proc )
    {
        this->M_n_localWithoutGhost_df[proc] = this->M_n_localWithGhost_df[proc];
        this->M_n_dofs+=this->M_n_localWithoutGhost_df[proc];
    }

    this->M_first_df_globalcluster[0]=0;
    if ( this->M_n_localWithoutGhost_df[0] > 0 )
        this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0]+this->M_n_localWithoutGhost_df[0]-1;
    else
        this->M_last_df_globalcluster[0] = this->M_first_df_globalcluster[0];

    for ( rank_type i=1; i<this->worldComm().size(); ++i )
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
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );
    //------------------------------------------------------------------------------//
    size_type firstGlobIndex = this->M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( size_type i=0; i< this->M_n_localWithGhost_df[myRank]; ++i )
    {
        this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
        ++nextGlobIndex;
    }
    //------------------------------------------------------------------------------//

} // buildGlobalProcessToGlobalClusterDofMapContinuousGhostDofBlockingComm

#endif
























//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
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
    //bool storeRangeActiveEltsTouchInterProcess = this->buildDofTableMPIExtended() && !mesh.components().test( MESH_UPDATE_FACES ) && !mesh.components().test( MESH_UPDATE_FACES_MINIMAL );

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
    if ( this->buildDofTableMPIExtended() )
    {
        this->buildGhostDofMapExtended( mesh );
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
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostDofMapExtended( mesh_type& mesh )
{
    DVLOG(2) << "[buildGhostDofMap] call buildGhostDofMapExtended on rank "<<  this->worldComm().rank() << "\n";
#if 0
    // extract range of elements
    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myActiveEltsTouchInterProcess( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );
    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myGhostEltsExtended( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );
    std::set<size_type> dofdoneActive, dofdoneGhost;
    auto rangeInterProcessFaces = (this->hasMeshSupport())? this->meshSupport()->rangeInterProcessFaces() : interprocessfaces(mesh);
    //for ( ; face_it!=face_en ; ++face_it )
    for ( auto const& faceWrap : rangeInterProcessFaces )
    {
        auto const& faceip = boost::unwrap_ref( faceWrap );//*face_it );
        auto const& elt0 = faceip.element0();
        auto const& elt1 = faceip.element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
        if ( dofdoneActive.find( eltOnProc.id() ) == dofdoneActive.end() )
        {
            dofdoneActive.insert( eltOnProc.id() );
            myActiveEltsTouchInterProcess->push_back(boost::cref(eltOnProc));
        }
        if ( dofdoneGhost.find( eltOffProc.id() ) == dofdoneGhost.end() )
        {
            myGhostEltsExtended->push_back(boost::cref(eltOffProc));
            dofdoneGhost.insert( eltOffProc.id() );
        }
    }
    auto myrangeActive = range(_range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myActiveEltsTouchInterProcess->begin(),myActiveEltsTouchInterProcess->end(),myActiveEltsTouchInterProcess ),
                               _mesh=mesh );
    auto myrangeGhost = range(_range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myGhostEltsExtended->begin(),myGhostEltsExtended->end(),myGhostEltsExtended ),
                              _mesh=mesh );
#endif

    if ( this->hasMeshSupport() && this->meshSupport()->isPartialSupport() )
        this->buildGhostDofMapExtended( mesh, elements(this->meshSupport(),entity_process_t::GHOST_ONLY ) );
    else
        this->buildGhostDofMapExtended( mesh, elements(mesh,entity_process_t::GHOST_ONLY ) );

}

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostDofMapExtended( mesh_type& mesh,
                                                                                   Range<mesh_type,MESH_ELEMENTS> const& ghostEltRange )
{
    DVLOG(2) << "[buildGhostDofMap] call buildGhostDofMapExtended on rank "<<  this->worldComm().rank() << "\n";

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
