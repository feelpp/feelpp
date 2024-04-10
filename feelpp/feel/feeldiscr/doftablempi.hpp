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
    if ( this->hasMeshSupport() )
        this->meshSupport()->updateParallelData();

    if ( true )//!mesh.components().test( MESH_UPDATE_FACES ) && !mesh.components().test( MESH_UPDATE_FACES_MINIMAL ) )
    {
        this->buildGlobalProcessToGlobalClusterDofMapOthersMesh( mesh );
    }
    else
    {
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

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

namespace detail {

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
    auto rangeInterProcessFaces = (this->hasMeshSupport())? this->meshSupport()->rangeInterProcessFaces() : interprocessfaces(mesh);
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

    size_type nLocalDofWithGhost = this->M_n_localWithGhost_df[myRank];
    std::vector<bool> dofdone( nLocalDofWithGhost,false);
    std::vector<bool> dofIsGhost( nLocalDofWithGhost,false);
    size_type nDofNotPresent=0;
    std::vector< std::map<size_type,std::vector< std::vector<std::pair<uint16_type,size_type> > > > > listToSend(this->worldComm().localSize());

    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myActiveEltsTouchInterProcess( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );

    bool hasMeshSupportPartial = this->hasMeshSupport() && this->meshSupport()->isPartialSupport();
    bool storeRangeActiveEltsTouchInterProcess = this->buildDofTableMPIExtended() && !mesh.components().test( MESH_UPDATE_FACES ) && !mesh.components().test( MESH_UPDATE_FACES_MINIMAL );

    if ( is_continuous )
    {
        auto rangeElements = (this->hasMeshSupport())? this->meshSupport()->rangeElements() : elements(mesh);
        for ( auto const& activeEltWrap : rangeElements )
        {
            auto const& activeElt = boost::unwrap_ref( activeEltWrap );
            bool findActiveEltTouchInterProcess = false;
            for ( uint16_type n=0; n < activeElt.nVertices(); n++ )
            {
                if ( hasMeshSupportPartial )
                {
                    for ( auto const& ghostData : activeElt.point(n).elementsGhost() )
                    {
                        if ( findActiveEltTouchInterProcess )
                            break;
                        for ( size_type _ghostEltId  : ghostData.second )
                        {
                            if ( this->meshSupport()->hasGhostElement( _ghostEltId ) )
                            {
                                findActiveEltTouchInterProcess = true;
                                break;
                            }
                        }
                    }
                }
                else
                {
                    if ( activeElt.point(n).numberOfProcGhost() > 0 )
                    {
                        findActiveEltTouchInterProcess = true;
                        break;
                    }
                }
            }

            if ( !findActiveEltTouchInterProcess )
                continue;

            if ( storeRangeActiveEltsTouchInterProcess )
                myActiveEltsTouchInterProcess->push_back(boost::cref(activeElt));

            // prepare local dofs to analyse
            std::map<uint16_type,std::map<uint16_type,std::tuple<uint16_type,size_type> > > mapLocalDofByCompToLocalDofAllComp;
            for ( auto const& ldof : this->localDof( activeElt.id() ) )
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

                rank_type pidDofActive = invalid_rank_type_value;
                size_type idEltInPartition = invalid_v<size_type>;
                if ( locDof < nLocalDofUpToVertices )
                {
                    uint16_type pointIdInElt = locDof / nDofPerVertexForDivision;
                    boost::tie( pidDofActive, idEltInPartition ) = Feel::detail::updateDofOnVertices( *this,activeElt,pointIdInElt );
                }
                else if ( nDim == 3 && locDof < nLocalDofUpToEdges )
                {
                    uint16_type locDofInEgdes = locDof - nLocalDofUpToVertices;
                    uint16_type edgeIdInElt = locDofInEgdes / nDofPerEdgeForDivision;
                    boost::tie( pidDofActive, idEltInPartition ) = Feel::detail::updateDofOnEdges( *this,activeElt,edgeIdInElt );
                }
                else if ( locDof < nLocalDofUpToFaces )
                {
                    if ( nDim < 2 && nDim != nRealDim )
                        continue;
                    uint16_type locDofInTopologicalFaces = locDof - nLocalDofBeforeDofsOfTopologicalFaces;
                    uint16_type faceIdInElt = locDofInTopologicalFaces / nDofPerTopologicalFaceForDivision;
                    auto facePtr = activeElt.facePtr( faceIdInElt );
                    if ( !facePtr )
                        continue;
                    if ( !facePtr->isInterProcessDomain() )
                        continue;
                    auto const& elt0 = facePtr->element0();
                    auto const& elt1 = facePtr->element1();
                    auto const& eltGhost = ( elt0.isGhostCell() )? elt0 : elt1;
                    if ( eltGhost.processId() < myRank )
                    {
                        if ( hasMeshSupportPartial )
                        {
                            if ( !this->meshSupport()->hasGhostElement( eltGhost.id() ) )
                                continue;
                        }
                        pidDofActive = eltGhost.processId();
                        idEltInPartition = eltGhost.idInOthersPartitions( pidDofActive );
                    }
                }
                else
                {
                    // dof inside the element : do nothing
                }

                // if dof is ghost -> prepare send/recv
                if ( pidDofActive != invalid_rank_type_value && pidDofActive < myRank )
                {
                    std::vector<std::pair<uint16_type,size_type> > compglobdofs;
                    for ( auto const& localDofAllComp : localDofDatas.second )
                    {
                        const size_type theglobdof = std::get<1>( localDofAllComp.second );
                        const uint16_type thelocdof = std::get<0>( localDofAllComp.second );
                        dofIsGhost[theglobdof] = true;
                        compglobdofs.push_back( std::make_pair(thelocdof,theglobdof) );
                        ++nDofNotPresent;
                    }

                    listToSend[pidDofActive][idEltInPartition].push_back( compglobdofs );
                }

            } // loop over local dof
        } // for ( auto const& activeElt : elements(mesh) )
    } // is_continuous
    else if ( storeRangeActiveEltsTouchInterProcess ) // discontinuous case maybe
    {
        auto rangeElements = (this->hasMeshSupport())? this->meshSupport()->rangeElements() : elements(mesh);
        for ( auto const& activeEltWrap : rangeElements )
        {
            auto const& activeElt = boost::unwrap_ref( activeEltWrap );
            bool findActiveEltTouchInterProcess = false;
            for ( uint16_type n=0; n < activeElt.nVertices(); n++ )
            {
                if ( hasMeshSupportPartial )
                {
                    for ( auto const& ghostData : activeElt.point(n).elementsGhost() )
                    {
                        if ( findActiveEltTouchInterProcess )
                            break;
                        for ( size_type _ghostEltId  : ghostData.second )
                        {
                            if ( this->meshSupport()->hasGhostElement( _ghostEltId ) )
                            {
                                findActiveEltTouchInterProcess = true;
                                break;
                            }
                        }
                    }
                }
                else
                {
                    if ( activeElt.point(n).numberOfProcGhost() > 0 )
                    {
                        findActiveEltTouchInterProcess = true;
                        break;
                    }
                }
            }

            if ( findActiveEltTouchInterProcess )
                myActiveEltsTouchInterProcess->push_back(boost::cref(activeElt));
        }
    }


    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//

    // update datamap info
    CHECK( this->M_n_localWithGhost_df[myRank] >= nDofNotPresent ) << "invalid data\n" << std::endl;
    this->M_n_localWithoutGhost_df[myRank] = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;

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
    this->buildGlobalProcessToGlobalClusterDofMapOthersMeshNonBlockingComm( mesh,listToSend );

   //------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------//

    // extended dof table
    if ( this->buildDofTableMPIExtended() )
    {
        if ( mesh.components().test( MESH_UPDATE_FACES ) || mesh.components().test( MESH_UPDATE_FACES_MINIMAL ) )
        {
            this->buildGhostDofMapExtended( mesh );
        }
        else
        {
            CHECK( !hasMeshSupportPartial ) << "TODO";
            auto myrangeActive = range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                                    myActiveEltsTouchInterProcess->begin(),myActiveEltsTouchInterProcess->end(),myActiveEltsTouchInterProcess ),
                                        _mesh=mesh );

            typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );
            auto rangeGhostElement = mesh.ghostElements();
            auto itGhostElt = std::get<0>( rangeGhostElement );
            auto enGhostElt = std::get<1>( rangeGhostElement );
            for ( ; itGhostElt != enGhostElt ; ++itGhostElt )
            {
                auto const& ghostelt = boost::unwrap_ref( *itGhostElt );

                for ( uint16_type f =0;f < ghostelt.nTopologicalFaces(); ++f )
                {
                    bool allFaceVerticesAreInActiveMeshPart = true;
                    for ( uint16_type p=0;p<mesh_type::element_type::topological_face_type::numVertices;++p )
                    {
                        //uint16_type ptIdInElt = ghostelt.fToP( f, p );
                        if ( ghostelt.point( ghostelt.fToP( f, p ) ).isGhostCell() )
                        {
                            allFaceVerticesAreInActiveMeshPart = false;
                            break;
                        }
                    }
                    if ( allFaceVerticesAreInActiveMeshPart )
                    {
                        myelts->push_back(boost::cref(ghostelt));
                        break;
                    }
                }

            }

            // generate a range object
            auto myrange = range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                              myelts->begin(),myelts->end(),myelts ), _mesh=mesh );
            DVLOG(2) << "ghost element in doftable : nelements(myrange) ["<<mesh.worldComm().rank()<<"] : " << nelements(myrange) << "\n";

            this->buildGhostDofMapExtended( mesh, myrange, myrangeActive );
        }
    }
}

template<typename MeshType, typename FEType, typename PeriodicityType,typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType,MortarType>::buildGlobalProcessToGlobalClusterDofMapOthersMeshNonBlockingComm( mesh_type& mesh,
                                                                                                                          std::vector< std::map<size_type,std::vector< std::vector<std::pair<uint16_type,size_type> > > > > const& listToSend )
{
    typedef std::vector< boost::tuple<std::vector<uint16_type>, ublas::vector<double> > > mpidofs_subcontainer_type; // tuple : ldof of each component + node
    typedef boost::tuple<size_type, mpidofs_subcontainer_type > mpidofs_container_type; // tuple : eltId + localDofs
    typedef std::vector< mpidofs_container_type > dofs_container_to_send_type;

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    //--------------------------------------------------------------------------------------------------------//
    std::map< rank_type,  dofs_container_to_send_type> dataToSend;
    // prepare container to send
    std::map< rank_type, std::vector< std::vector<std::vector<size_type> > > > memoryInitialRequest;
    std::map< rank_type, int > nDataInVecToSendBis;
    for ( rank_type proc=0; proc<this->worldComm().size(); ++proc )
    {
        if ( listToSend[proc].size() == 0 ) continue;

        auto itElements = listToSend[proc].begin();
        auto const enElements = listToSend[proc].end();
        const int nEltToSend = std::distance(itElements,enElements);
        memoryInitialRequest[proc].resize(nEltToSend);

        dataToSend[proc].resize(nEltToSend);
        for ( int cptElt=0 ; itElements!=enElements ; ++itElements, ++cptElt)
        {
            auto itDof = itElements->second.begin();
            auto const enDof = itElements->second.end();
            const int nDofsInEltToSearch = std::distance(itDof,enDof);

            mpidofs_subcontainer_type dofsInEltContainer(nDofsInEltToSearch);
            memoryInitialRequest[proc][cptElt].resize(nDofsInEltToSearch);
            for (int cptDof=0; itDof!=enDof ; ++itDof,++cptDof)
            {
                uint16_type nCompOfLocalDof = itDof->size();
                memoryInitialRequest[proc][cptElt][cptDof].resize(nCompOfLocalDof);
                for (uint16_type comp=0; comp<nCompOfLocalDof ; ++comp )
                {
                    const uint16_type thelocdof = itDof->operator[](comp).first;
                    const size_type theglobdof = itDof->operator[](comp).second;
                    // save the tag of mpi send
                    memoryInitialRequest[proc][cptElt][cptDof][comp] = theglobdof;
                    //------------------------------------------------------------------------------//
                    boost::get<0>( dofsInEltContainer[cptDof] ).push_back( thelocdof );
                    if (comp==0)
                    {
                        // get the dof node
                        ublas::vector<double> nodeDofToSend( nRealDim );
                        auto itFindDofPoint = M_dof_points.find( theglobdof );
                        CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                        nodeDofToSend[0]=itFindDofPoint->second.template get<0>()[0];
                        if ( nRealDim>1 )
                            nodeDofToSend[1]=itFindDofPoint->second.template get<0>()[1];
                        if ( nRealDim>2 )
                            nodeDofToSend[2]=itFindDofPoint->second.template get<0>()[2];
                        boost::get<1>( dofsInEltContainer[cptDof] ) = nodeDofToSend;
                    }
                }
            }
            dataToSend[proc][cptElt] = boost::make_tuple(itElements->first,dofsInEltContainer);
        } // for ( int cptElt=0 ; itElements!=enElements ; ++itElements, ++cptElt)
    }

    //--------------------------------------------------------------------------------------------------------//
    // counter of request
    int nbRequest = 2*mesh.neighborSubdomains().size();
    if ( nbRequest ==0 ) return;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    //--------------------------------------------------------------------------------------------------------//
    // first send/recv
    std::map<rank_type,dofs_container_to_send_type> dataToRecv;
    for ( const rank_type procNeigborId : mesh.neighborSubdomains() )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( procNeigborId , 0, dataToSend[procNeigborId] );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( procNeigborId , 0, dataToRecv[procNeigborId] );
    }
    //--------------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //--------------------------------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector< std::vector< std::vector<size_type> > > > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        auto itEltRecv = itDataRecv->second.begin();
        auto const enEltRecv = itDataRecv->second.end();
        const int nEltRecv=  itDataRecv->second.size();
        dataToReSend[idProc].resize( nEltRecv );
        for ( int cptElt=0 ; itEltRecv!=enEltRecv ; ++itEltRecv,++cptElt )
        {
            auto const idEltInMyPartition = itEltRecv->template get<0>();
            DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                    << "idEltInMyPartition: " << idEltInMyPartition << "\n";

            auto const& theelt = mesh.element( idEltInMyPartition );
            double dofPtCompareTol = std::max(1e-15,theelt.hMin()*1e-5);

            // prepare local dofs in element for localization process
            std::map<uint16_type,std::tuple<node_type, std::map<uint16_type,std::tuple<uint16_type,size_type> > > > mapLocalDofByCompToLocalDofAllCompWithNode;
            for ( auto const& ldof : this->localDof( idEltInMyPartition ) )
            {
                uint16_type ldofId = ldof.first.localDof();
                uint16_type ldofParentId = this->fe().dofParent( ldofId );
                size_type gdofId = ldof.second.index();
                uint16_type comp = this->fe().component(ldofId);

                if ( mapLocalDofByCompToLocalDofAllCompWithNode.find( ldofParentId ) == mapLocalDofByCompToLocalDofAllCompWithNode.end() )
                {
                    auto itFindDofPoint = M_dof_points.find( gdofId  );
                    CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built :" << gdofId;
                    auto const& thedofPtInElt = itFindDofPoint->second.template get<0>();
                    std::get<0>( mapLocalDofByCompToLocalDofAllCompWithNode[ldofParentId] ) = thedofPtInElt;
                }
                std::get<1>( mapLocalDofByCompToLocalDofAllCompWithNode[ldofParentId] )[comp] = std::make_tuple( ldofId,gdofId );
            }

            // loop over each dof in the mesh element
            auto itDofInElt = itEltRecv->template get<1>().begin();
            auto const enDofInElt = itEltRecv->template get<1>().end();
            const int nDofInElt = distance(itDofInElt,enDofInElt);
            dataToReSend[idProc][cptElt].resize( nDofInElt );
            for ( int cptDofInElt=0 ; itDofInElt != enDofInElt ; ++itDofInElt,++cptDofInElt )
            {
                std::vector<uint16_type> const& locDofToSearch = boost::get<0>( *itDofInElt );
                CHECK( !locDofToSearch.empty() ) << "wrong data ";
                uint16_type locDofToSearchType = this->fe().dofType( locDofToSearch.front() );
                std::map<uint16_type,uint16_type> mapComponentToCommContainerId;
                for ( uint16_type k=0;k<locDofToSearch.size();++k)
                    mapComponentToCommContainerId[ this->fe().component( locDofToSearch[k] ) ] = k;

                auto const& nodeDofRecv = boost::get<1>( *itDofInElt );
                //------------------------------------------------------------------------------//
                // search dof on elt id recv
                uint16_type locDofFound = invalid_uint16_type_value;
                for ( auto const& localDofDatas : mapLocalDofByCompToLocalDofAllCompWithNode )
                {
                    uint16_type l = localDofDatas.first;
                    // ignore local dof with different type
                    if ( locDofToSearchType != this->fe().dofType( l ) )
                         continue;
                    // dof point in element
                    auto const& thedofPtInElt = std::get<0>( localDofDatas.second );
                    DVLOG(3) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                             << "thedofPtInElt: " << thedofPtInElt << "nodeDofRecv: " << nodeDofRecv << "\n";
                    // test equatlity of dofs point
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                    {
                        find2 = find2 && (std::abs( thedofPtInElt[d]-nodeDofRecv[d] )<dofPtCompareTol);
                    }
                    // if find else save local dof
                    if (find2)
                    {
                        locDofFound = l;
                        break;
                    }
                }
                //------------------------------------------------------------------------------//
                // check
                CHECK( locDofFound != invalid_uint16_type_value ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";
                //------------------------------------------------------------------------------//
                auto const& locDofFoundAllComp = std::get<1>( mapLocalDofByCompToLocalDofAllCompWithNode.find( locDofFound )->second );
                dataToReSend[idProc][cptElt][cptDofInElt].resize( mapComponentToCommContainerId.size(),invalid_v<size_type> );
                for ( auto const& locDofByComp : locDofFoundAllComp )
                {
                    uint16_type thecomp = locDofByComp.first;
                    auto itFindComp = mapComponentToCommContainerId.find( thecomp );
                    if ( itFindComp == mapComponentToCommContainerId.end() )
                        continue;
                    uint16_type compContainerId = itFindComp->second;
                    //uint16_type ldofIdAsked = std::get<0>( locDofByComp.second );
                    const size_type dofGlobFound = std::get<1>( locDofByComp.second );
                    // save response
                    dataToReSend[idProc][cptElt][cptDofInElt][compContainerId] = this->M_mapGlobalProcessToGlobalCluster[dofGlobFound];
                    this->M_activeDofSharedOnCluster[dofGlobFound].insert(idProc);
                }
            }
        } // for ( int cptElt=0 ... )
    } // for ( ; itDataRecv ... )

    //--------------------------------------------------------------------------------------------------------//
    // send/recv respond to the request
    cptRequest=0;
    std::map<rank_type, std::vector<std::vector<std::vector<size_type> > > > finalDataToRecv;
    for ( const rank_type procNeigborId : mesh.neighborSubdomains() )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( procNeigborId, 0, dataToReSend[procNeigborId] );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( procNeigborId, 0, finalDataToRecv[procNeigborId] );
    }
    //--------------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //--------------------------------------------------------------------------------------------------------//
    // update datamap for ghost dof
    for ( auto const& finalDataByProc : finalDataToRecv )
    {
        const rank_type idProc = finalDataByProc.first;
        auto itEltRecv = finalDataByProc.second.begin();
        auto const enEltRecv = finalDataByProc.second.end();
        for ( int cptElt=0 ; itEltRecv!=enEltRecv ; ++itEltRecv,++cptElt )
        {
            for ( uint16_type ll=0;ll<itEltRecv->size();++ll )
            {
                auto const& gdofRecvAllComp = itEltRecv->operator[](ll);
                for ( uint16_type cc=0;cc<gdofRecvAllComp.size();++cc)
                {
                    const size_type myGlobProcessDof = memoryInitialRequest[idProc][cptElt][ll][cc];
                    const size_type dofGlobRecv = gdofRecvAllComp[cc];
                    DCHECK( dofGlobRecv != invalid_v<size_type> ) << "a gdof is missing from mpi recv";
                    //update data map
                    this->M_mapGlobalProcessToGlobalCluster[myGlobProcessDof] = dofGlobRecv;
                }
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

    // extract range of elements
    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myActiveEltsTouchInterProcess( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );
    typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype myGhostEltsExtended( new typename MeshTraits<mesh_type>::elements_reference_wrapper_type );

    std::set<size_type> dofdoneActive, dofdoneGhost;
#if 0
    auto rangeInterProcessFaces = mesh.interProcessFaces();
    auto face_it = std::get<0>( rangeInterProcessFaces );
    auto const face_en = std::get<1>( rangeInterProcessFaces );
#endif
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

    this->buildGhostDofMapExtended( mesh, myrangeGhost, myrangeActive );

}


template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostDofMapExtended( mesh_type& mesh,
                                                                                   Range<mesh_type,MESH_ELEMENTS> const& ghostEltRange,
                                                                                   Range<mesh_type,MESH_ELEMENTS> const& activeEltTouchInterProcessRange )
{
    DVLOG(2) << "[buildGhostDofMap] call buildGhostDofMapExtended on rank "<<  this->worldComm().rank() << "\n";

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    size_type start_next_free_dof = this->M_last_df[myRank]+1;
    //------------------------------------------------------------------------------//
    // build extended dof table
    size_type next_free_dof = start_next_free_dof;
    DofFromElement<self_type,fe_type> dfe( this, *M_fe );
    DofFromBoundary<self_type, fe_type> dfb( this, *M_fe );
    std::set<size_type> faceGhostDone;
    for ( auto const& ghostEltWrap : ghostEltRange )
    {
        auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
        // elements doftable
        if ( !this->isElementDone( ghostElt.id() ) )
            dfe.add( ghostEltWrap, next_free_dof, myRank );
        // faces doftable
        for ( size_type f = 0; f < ghostElt.nTopologicalFaces(); f++ )
        {
            if ( !ghostElt.facePtr(f) )
                continue;
            auto const& theface = ghostElt.face(f);
            if ( theface.isGhostCell() && faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
            {
                M_face_l2g[ theface.id()].resize( nLocalDofOnFace() );
                dfb.add( theface );
                faceGhostDone.insert( theface.id() );
            }
        }
    }
    this->M_nGhostDofAddedInExtendedDofTable = next_free_dof-start_next_free_dof;
    //------------------------------------------------------------------------------//
    // update local datamap
    bool hasNoElt = this->M_n_localWithGhost_df[myRank] == 0;
    const size_type thelastDof = ( !hasNoElt )?next_free_dof-1:0;

    std::vector<size_type> dataRecvFromGather;
    mpi::all_gather( this->worldComm().localComm(),
                     thelastDof,
                     dataRecvFromGather );

    for (rank_type p=0;p<nProc;++p)
    {
        this->M_last_df[p] = dataRecvFromGather[p];

        size_type mynDofWithGhost = ( this->M_n_localWithGhost_df[p]>0 )?
            this->M_last_df[p] - this->M_first_df[p] + 1 : 0;
        this->M_n_localWithGhost_df[p] = mynDofWithGhost;
    }

    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_v<size_type> );

    //------------------------------------------------------------------------------//
    // generate dof point on extended part
    this->generateDofPoints( ghostEltRange );

    //------------------------------------------------------------------------------//

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

    //
    std::map< rank_type, std::map<size_type,std::vector< std::vector<std::pair<uint16_type,size_type> > > > >  dofsDataToUpdate;
    std::unordered_map<size_type,size_type> mapGlobalClusterToGlobalProcessAroundInterProcess;
    std::vector<bool> dofdone( this->M_nGhostDofAddedInExtendedDofTable,false);

    // get dofs in extended part
    for ( auto const& ghostEltWrap : ghostEltRange )
    {
        auto const& ghostElt = boost::unwrap_ref( ghostEltWrap );
        const rank_type processIdOfGhost = ghostElt.processId();
        const size_type eltIdOfGhost = ghostElt.id();
        const size_type eltIdInOtherPartOfGhost = ghostElt.idInOthersPartitions(processIdOfGhost);

        auto & dofsDataToUpdateOnCurrentElement = dofsDataToUpdate[processIdOfGhost][eltIdInOtherPartOfGhost];
        // prepare local dofs to analyse
        std::map<uint16_type,std::map<uint16_type,std::tuple<uint16_type,size_type> > > mapLocalDofByCompToLocalDofAllComp;
        for ( auto const& ldof : this->localDof( eltIdOfGhost ) )
        {
            uint16_type ldofId = ldof.first.localDof();

            if ( is_tensor2symm )
                if ( localDofUsedForTensor2symm.find( ldofId ) == localDofUsedForTensor2symm.end() )
                    continue;

            size_type gdofId = ldof.second.index();
            if ( gdofId < start_next_free_dof )
                continue;

            if ( dofdone[gdofId-start_next_free_dof] )
                continue;
            dofdone[gdofId-start_next_free_dof] = true;

            uint16_type ldofParentId = this->fe().dofParent( ldofId );
            uint16_type comp = this->fe().component(ldofId);
            mapLocalDofByCompToLocalDofAllComp[ldofParentId][comp] = std::make_tuple( ldofId,gdofId );
        }

        // loop over local dof for the detection of ghost dofs
        for (auto const& localDofDatas : mapLocalDofByCompToLocalDofAllComp )
        {
            DCHECK( !localDofDatas.second.empty() ) << "no localdof data is empty";
            std::vector<std::pair<uint16_type,size_type> > compglobdofs;
            for ( auto const& localDofAllComp : localDofDatas.second )
            {
                const size_type theglobdof = std::get<1>( localDofAllComp.second );
                const uint16_type thelocdof = std::get<0>( localDofAllComp.second );
                compglobdofs.push_back( std::make_pair(thelocdof,theglobdof) );
            }
            dofsDataToUpdateOnCurrentElement.push_back( compglobdofs );
        }
    }

    // store mapping between cluster to process
    for( auto const& eltWrap : activeEltTouchInterProcessRange )
    {
        auto const& elt = boost::unwrap_ref( eltWrap );
        if ( elt.isGhostCell() )
            continue;

        for ( auto const& ldof : this->localDof( elt.id() ) )
        {
            size_type gdofId = ldof.second.index();
            mapGlobalClusterToGlobalProcessAroundInterProcess[ this->mapGlobalProcessToGlobalCluster()[gdofId] ] = gdofId;

        }
    }

    typedef std::vector< boost::tuple<std::vector<uint16_type>, ublas::vector<double> > > mpidofs_subcontainer_type; // tuple : ldof of each component + node
    typedef boost::tuple<size_type, mpidofs_subcontainer_type > mpidofs_container_type; // tuple : eltId + localDofs
    typedef std::vector< mpidofs_container_type > dofs_container_to_send_type;

    // prepare container to send
    std::map< rank_type, dofs_container_to_send_type > dataToSend, dataToRecv;
    std::map< rank_type, std::vector< std::vector<std::vector<size_type> > > > memoryInitialRequest;

    for ( auto const& dofsDataToUpdateByProc : dofsDataToUpdate )
    {
        const rank_type idProc = dofsDataToUpdateByProc.first;
        auto const& dofsDataToUpdateEltContainer = dofsDataToUpdateByProc.second;
        if ( dofsDataToUpdateEltContainer.empty() )
            continue;

        dataToSend[idProc].resize( dofsDataToUpdateEltContainer.size() );
        memoryInitialRequest[idProc].resize( dofsDataToUpdateEltContainer.size() );
        auto itEltGhost = dofsDataToUpdateEltContainer.begin();
        auto const enEltGhost = dofsDataToUpdateEltContainer.end();
        for ( size_type cptElt = 0 ; itEltGhost!=enEltGhost; ++itEltGhost,++cptElt )
        {
            const size_type idElt = itEltGhost->first;

            auto const& dofsDataToUpdateLocalDofContainer = itEltGhost->second;
            auto itDof = itEltGhost->second.begin();
            auto const enDof = itEltGhost->second.end();
            const int nDofsInEltToSearch = dofsDataToUpdateLocalDofContainer.size();
            mpidofs_subcontainer_type dofsInEltContainer(nDofsInEltToSearch);
            memoryInitialRequest[idProc][cptElt].resize(nDofsInEltToSearch);

            for ( uint16_type cptDof=0 ; itDof!=enDof ; ++itDof, ++cptDof )
            {
                uint16_type nCompOfLocalDof = itDof->size();
                memoryInitialRequest[idProc][cptElt][cptDof].resize(nCompOfLocalDof);

                for ( uint16_type comp=0; comp<nCompOfLocalDof ; ++comp )
                {
                    const uint16_type thelocdof = itDof->operator[](comp).first;
                    const size_type theglobdof = itDof->operator[](comp).second;
                    // save the tag of mpi send
                    memoryInitialRequest[idProc][cptElt][cptDof][comp] = theglobdof;
                    //------------------------------------------------------------------------------//
                    boost::get<0>( dofsInEltContainer[cptDof] ).push_back( thelocdof );
                    if (comp==0)
                    {
                        // get the dof node
                        ublas::vector<double> nodeDofToSend( nRealDim );
                        // nodeDofToSend.resize( nRealDim );
                        auto itFindDofPoint = M_dof_points.find( theglobdof );
                        CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built";
                        nodeDofToSend[0]=itFindDofPoint->second.template get<0>()[0];
                        if ( nRealDim>1 )
                            nodeDofToSend[1]=itFindDofPoint->second.template get<0>()[1];
                        if ( nRealDim>2 )
                            nodeDofToSend[2]=itFindDofPoint->second.template get<0>()[2];

                        boost::get<1>( dofsInEltContainer[cptDof] ) = nodeDofToSend;
                    }
                }
            }

            dataToSend[idProc][cptElt] = boost::make_tuple(idElt,dofsInEltContainer);
        }
    }

    //------------------------------------------------------------------------------//
    // counter of request
    int nbRequest=2*mesh.neighborSubdomains().size();
    if ( nbRequest == 0 ) return;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    // first send/recv
    for ( const rank_type procNeigborId : mesh.neighborSubdomains() )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( procNeigborId , 0, dataToSend[procNeigborId] );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( procNeigborId , 0, dataToRecv[procNeigborId] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector< std::vector< std::vector<size_type> > > > dataToReSend;
    std::map<rank_type, std::vector< boost::tuple<int,size_type> > > dataToSendNewNeigbor;
    for ( rank_type p : this->neighborSubdomains() )
        dataToSendNewNeigbor[p].clear();

    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        auto itEltRecv = itDataRecv->second.begin();
        auto const enEltRecv = itDataRecv->second.end();
        dataToReSend[idProc].resize( std::distance(itEltRecv,enEltRecv) );
        for ( int cptElt=0 ; itEltRecv!=enEltRecv ; ++itEltRecv,++cptElt )
        {
            const size_type idEltToSearch = itEltRecv->template get<0>();

            auto const& theelt = mesh.element( idEltToSearch );
            double dofPtCompareTol = std::max(1e-15,theelt.hMin()*1e-5);

            // prepare local dofs in element for localization process
            std::map<uint16_type,std::tuple<node_type, std::map<uint16_type,std::tuple<uint16_type,size_type> > > > mapLocalDofByCompToLocalDofAllCompWithNode;
            for ( auto const& ldof : this->localDof( idEltToSearch ) )
            {
                uint16_type ldofId = ldof.first.localDof();
                uint16_type ldofParentId = this->fe().dofParent( ldofId );
                size_type gdofId = ldof.second.index();
                uint16_type comp = this->fe().component(ldofId);

                if ( mapLocalDofByCompToLocalDofAllCompWithNode.find( ldofParentId ) == mapLocalDofByCompToLocalDofAllCompWithNode.end() )
                {
                    auto itFindDofPoint = M_dof_points.find( gdofId  );
                    CHECK( itFindDofPoint != M_dof_points.end() ) << "dof point is not built :" << gdofId;
                    auto const& thedofPtInElt = itFindDofPoint->second.template get<0>();
                    std::get<0>( mapLocalDofByCompToLocalDofAllCompWithNode[ldofParentId] ) = thedofPtInElt;
                }
                std::get<1>( mapLocalDofByCompToLocalDofAllCompWithNode[ldofParentId] )[comp] = std::make_tuple( ldofId,gdofId );
            }

            // loop over each dof in the mesh element
            auto itDofInElt = itEltRecv->template get<1>().begin();
            auto const enDofInElt = itEltRecv->template get<1>().end();
            const int nDofInElt = distance(itDofInElt,enDofInElt);
            dataToReSend[idProc][cptElt].resize( nDofInElt );
            for ( int cptDof=0; itDofInElt!=enDofInElt ; ++itDofInElt, ++cptDof )
            {
                std::vector<uint16_type> const& locDofToSearch = boost::get<0>( *itDofInElt );
                CHECK( !locDofToSearch.empty() ) << "wrong data ";
                uint16_type locDofToSearchType = this->fe().dofType( locDofToSearch.front() );

                std::map<uint16_type,uint16_type> mapComponentToCommContainerId;
                for ( uint16_type k=0;k<locDofToSearch.size();++k)
                    mapComponentToCommContainerId[ this->fe().component( locDofToSearch[k] ) ] = k;


                auto const& nodeDofRecv = boost::get<1>( *itDofInElt );

                uint16_type locDofFound = invalid_uint16_type_value;
                for ( auto const& localDofDatas : mapLocalDofByCompToLocalDofAllCompWithNode )
                {
                    uint16_type l = localDofDatas.first;
                    // ignore local dof with different type
                    if ( locDofToSearchType != this->fe().dofType( l ) )
                        continue;
                    // dof point in element
                    auto const& thedofPtInElt = std::get<0>( localDofDatas.second );
                    DVLOG(3) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                             << "thedofPtInElt: " << thedofPtInElt << "nodeDofRecv: " << nodeDofRecv << "\n";
                    // test equatlity of dofs point
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                    {
                        find2 = find2 && (std::abs( thedofPtInElt[d]-nodeDofRecv[d] )<dofPtCompareTol);
                    }
                    // if find else save local dof
                    if (find2)
                    {
                        locDofFound = l;
                        break;
                    }
                }
                CHECK( locDofFound != invalid_uint16_type_value ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";

                auto const& locDofFoundAllComp = std::get<1>( mapLocalDofByCompToLocalDofAllCompWithNode.find( locDofFound )->second );
                dataToReSend[idProc][cptElt][cptDof].resize( mapComponentToCommContainerId.size(),invalid_v<size_type> );
                for ( auto const& locDofByComp : locDofFoundAllComp )
                {
                    uint16_type thecomp = locDofByComp.first;
                    auto itFindComp = mapComponentToCommContainerId.find( thecomp );
                    if ( itFindComp == mapComponentToCommContainerId.end() )
                        continue;
                    uint16_type compContainerId = itFindComp->second;
                    //uint16_type ldofId = std::get<0>( locDofByComp.second );
                    const size_type dofGlobFound = std::get<1>( locDofByComp.second );

                    size_type gcdofFound = this->M_mapGlobalProcessToGlobalCluster[dofGlobFound];
                    // save response
                    dataToReSend[idProc][cptElt][cptDof][compContainerId] = gcdofFound;
                    if ( !this->dofGlobalProcessIsGhost( dofGlobFound ) )
                        this->M_activeDofSharedOnCluster[dofGlobFound].insert(idProc);
                    else
                    {
                        rank_type activeProcId = this->procOnGlobalCluster( gcdofFound );
                        if ( activeProcId != idProc )
                            dataToSendNewNeigbor[activeProcId].push_back( boost::make_tuple(idProc,gcdofFound) );
                    }
                }
            } // for ( int cptDof=0; itDofInElt ... )
        }
    } // for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )

    //------------------------------------------------------------------------------//
    // send respond to the request
    cptRequest=0;
    // second send/recv
    std::map<rank_type, std::vector<std::vector<std::vector<size_type> > > > finalDataToRecv;

    for ( const rank_type procNeigborId : mesh.neighborSubdomains() )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( procNeigborId , 0, dataToReSend[procNeigborId] );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( procNeigborId , 0, finalDataToRecv[procNeigborId] );
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //------------------------------------------------------------------------------//
    // update datamap for ghost dof
    for ( auto const& finalDataByProc : finalDataToRecv )
    {
        const rank_type idProc = finalDataByProc.first;
        auto itEltRecv = finalDataByProc.second.begin();
        auto const enEltRecv = finalDataByProc.second.end();
        for ( int cptElt=0 ; itEltRecv!=enEltRecv ; ++itEltRecv, ++cptElt )
        {
            for ( uint16_type ll=0;ll<itEltRecv->size();++ll )
            {
                auto const& gdofRecvAllComp = itEltRecv->operator[](ll);
                for ( uint16_type cc=0;cc<gdofRecvAllComp.size();++cc)
                {
                    const size_type myGlobProcessDof = memoryInitialRequest[idProc][cptElt][ll][cc];
                    const size_type dofGlobRecv = gdofRecvAllComp[cc];
                    DCHECK( dofGlobRecv != invalid_v<size_type> ) << "a gdof is missing from mpi recv";
                    //update data map
                    this->M_mapGlobalProcessToGlobalCluster[myGlobProcessDof] = dofGlobRecv;
                    rank_type activeProcId = this->procOnGlobalCluster( dofGlobRecv );
                    if ( activeProcId != myRank )
                        this->addNeighborSubdomain( activeProcId );
                }
            }
        }
    }
    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    // others isend/irecv between neigboor part : get new neigbor due to the extended part
    std::map<rank_type, std::vector< boost::tuple<int,size_type> > > dataToRecvNewNeigbor;
    int nbRequestNewNeigbor = 2*dataToSendNewNeigbor.size();
    mpi::request * reqsNewNeigbor = new mpi::request[nbRequestNewNeigbor];
    cptRequest=0;
    for ( auto const& mydataSend : dataToSendNewNeigbor )
    {
        const rank_type idProc = mydataSend.first;
        reqsNewNeigbor[cptRequest++] = this->worldComm().localComm().isend( idProc, 0, mydataSend.second );
        reqsNewNeigbor[cptRequest++] = this->worldComm().localComm().irecv( idProc, 0, dataToRecvNewNeigbor[idProc] );
    }
    // wait all requests
    mpi::wait_all(reqsNewNeigbor, reqsNewNeigbor + nbRequestNewNeigbor);
    // delete reqs because finish comm
    delete [] reqsNewNeigbor;
    //------------------------------------------------------------------------------//
    // update info recv about newNeighbor
    for ( auto const& mydataRecv : dataToRecvNewNeigbor )
    {
        for ( auto const& newDofNeigbor : mydataRecv.second )
        {
            int idProcGhost = boost::get<0>( newDofNeigbor );
            if ( idProcGhost != myRank )
            {
                size_type gcdof = boost::get<1>( newDofNeigbor );
                CHECK( this->dofGlobalClusterIsOnProc( gcdof ) ) << "must be an active dof";
                auto itFindGpDof = mapGlobalClusterToGlobalProcessAroundInterProcess.find( gcdof );
                CHECK( itFindGpDof != mapGlobalClusterToGlobalProcessAroundInterProcess.end() ) << "gcdof not register : " << gcdof;
                size_type gpdof = itFindGpDof->second;
                this->M_activeDofSharedOnCluster[gpdof].insert(idProcGhost);
                this->addNeighborSubdomain( idProcGhost );
            }
        }
    }

}

} // namespace Feel

#endif /* FEELPP_DOFTABLE_MPI_HPP */
