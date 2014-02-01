/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
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
#if 0 //old

    // dofs map on interprocess faces : ( dofOnMyRank, tuple< otherProcRank, dofOnOtherRank> )
    std::map<size_type,boost::tuple<size_type,size_type> > mapInterProcessDof;
    // dofs no present on interprocess faces : ( dofOnMyRank, tuple< otherProcRank, dofOnOtherRank> )
    std::map<size_type,boost::tuple<size_type,size_type> > setInterProcessDofNotPresent;

    if (is_continuous)
    {
        DVLOG(2) << "[buildGhostDofMap] call buildGhostInterProcessDofMap() with god rank "<<  this->worldComm().godRank() << "\n";
        buildGhostInterProcessDofMap( mesh,mapInterProcessDof );

        DVLOG(2) << "[buildGhostDofMap] call buildDofNotPresent() with rank "<<  this->worldComm().rank() << "\n";
        buildDofNotPresent( mapInterProcessDof,setInterProcessDofNotPresent );
    }

    DVLOG(2) << "[buildGhostDofMap] call buildGlobalProcessToGlobalClusterDofMap() with rank "<<  this->worldComm().rank() << "\n";
    buildGlobalProcessToGlobalClusterDofMap( mesh,setInterProcessDofNotPresent );

#else

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

#endif


    DVLOG(2) << "[buildGhostDofMap] call localtoglobalOnCluster() with rank "<<  this->worldComm().rank() << "\n";
    auto it_elt = mesh.beginElementWithProcessId( this->comm().rank() );
    auto en_elt = mesh.endElementWithProcessId( this->comm().rank() );

    for ( ; it_elt != en_elt; ++it_elt )
    {
        size_type elid= it_elt->id();

        for ( int i = 0; i < FEType::nLocalDof; ++i )
        {
            int nc1 = ( is_product?nComponents1:1 );

            for ( int c1 =0; c1 < nc1; ++c1 )
            {
                int ind = FEType::nLocalDof*c1+i;
                auto const& dof = localToGlobalOnCluster( elid, i, c1 );

                M_locglobOnCluster_indices[elid][ind] = dof.index();
                M_locglobOnCluster_signs[elid][ind] = dof.sign();
            }
        }
    }

    DVLOG(2) << "[buildGhostDofMap] finish () with god rank "<< this->worldComm().godRank() << "\n";

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
    size_type nbFaceDof = invalid_size_type_value;
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
    std::set<int> procRecvData;
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

template <typename MeshType>
boost::tuple<int,size_type >
updateDofOnVertices( MeshType const& mesh, typename MeshType::face_type const& theface, const int myIdProcess,
                     const int IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                     const uint16_type locDof, std::set<int> & procRecvData )
{
    int procMin = IdProcessOfGhost;
    size_type idFaceMin = idFaceInPartition;

    uint16_type iFaEl = ( theface.processId() == theface.proc_first() )? theface.pos_first():theface.pos_second();
    //local point number (in element)
    uint16_type iPtEl = MeshType::element_type::fToP( iFaEl, locDof );

    auto const& thept = eltOnProc.point(iPtEl);
    auto const theptId = thept.id();

    auto itprocghost=thept.elementsGhost().begin();
    auto const enprocghost=thept.elementsGhost().end();
    for ( ; itprocghost!=enprocghost ; ++itprocghost)
    {
        if (procMin>itprocghost->first)
        {
            const int theprocGhost=itprocghost->first;
            bool findFace=false;
            DCHECK(itprocghost->second.size()>0) << "need to have at least one ghost element\n";
            auto iteltghost = itprocghost->second.begin();
            auto const& eltGhost = mesh.element(*iteltghost,theprocGhost);
            for ( uint16_type f = 0; f < MeshType::element_type::numTopologicalFaces && !findFace; ++f )
            {
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
            const int procIdGhost = itprocghost->first;
            procRecvData.insert( procIdGhost );
        }
    }

    return boost::make_tuple(procMin,idFaceMin);
}

//--------------------------------------------------------------------------------------------------------//

template <typename MeshType>
boost::tuple<int,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface,const int myIdProcess,
                  const int IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<int> & procRecvData, mpl::int_<0> /**/ )
{
    return boost::make_tuple(0,0);
}
template <typename MeshType>
boost::tuple<int,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const int myIdProcess,
                  const int IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<int> & procRecvData, mpl::int_<1> /**/ )
{
    return boost::make_tuple(0,0);
}
template <typename MeshType>
boost::tuple<int,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const int myIdProcess,
                  const int IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<int> & procRecvData, mpl::int_<2> /**/ )
{
    return boost::make_tuple(0,0);
}
template <typename MeshType>
boost::tuple<int,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const int myIdProcess,
                  const int IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<int> & procRecvData, mpl::int_<3> /**/ )
{
    int procMin = IdProcessOfGhost;
    size_type idFaceMin = idFaceInPartition;

    uint16_type iFaEl = ( theface.processId() == theface.proc_first() )? theface.pos_first():theface.pos_second();
    //local edge number (in element)
    uint16_type iEdEl = MeshType::element_type::fToE(  iFaEl, idEdgesInFace );

    auto const& theedge = eltOnProc.edge(iEdEl);
    auto const theedgeId = theedge.id();

    auto itprocghost=theedge.elementsGhost().begin();
    auto const enprocghost=theedge.elementsGhost().end();
    for ( ; itprocghost!=enprocghost ; ++itprocghost)
    {
        if ( procMin>itprocghost->first )
        {
            const int theprocGhost=itprocghost->first;

            DCHECK(itprocghost->second.size()>0) << "need to have at least one ghost element\n";
            auto iteltghost = itprocghost->second.begin();
            auto const& eltGhost = mesh.element(*iteltghost,theprocGhost);
            bool findFace=false;
            for ( uint16_type f = 0; f < MeshType::element_type::numTopologicalFaces && !findFace; ++f )
            {
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
            const int procIdGhost = itprocghost->first;
            procRecvData.insert( procIdGhost );
        }

    }

    return boost::make_tuple(procMin,idFaceMin);
}

//--------------------------------------------------------------------------------------------------------//

template <typename MeshType>
boost::tuple<int,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const int myIdProcess,
                  const int IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<int> & procRecvData )
{
    return updateDofOnEdges( mesh,theface,myIdProcess,IdProcessOfGhost,idFaceInPartition,eltOnProc,idEdgesInFace,procRecvData,mpl::int_<MeshType::nDim>() );
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
                                                                                                        std::set<int> & procRecvData )
{
    // goal init container listToSend
    // std::vector< std::map<size_type,std::set<size_type> > > ( proc,( idFace,(globDof,..)), ...   )) )
    const int myRank = this->worldComm().rank();
    //------------------------------------------------------------------------------//
    // get nbFaceDof
    size_type nbFaceDof = invalid_size_type_value;
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
    auto face_it = mesh.interProcessFaces().first;
    auto const face_en = mesh.interProcessFaces().second;
    DVLOG(2) << "[buildGhostInterProcessDofMap] nb interprocess faces: " << std::distance( face_it, face_en ) << "\n";
    for ( ; face_it != face_en ; ++face_it )
    {
        DVLOG(2) << "[buildGhostInterProcessDofMap] face id: " << face_it->id() << "\n";
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") eltOnProc id: "  << eltOnProc.id()  << "G(): " << eltOnProc.G() << "\n";
        DVLOG(2) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") eltOffProc id: " << eltOffProc.id() << "G(): " << eltOffProc.G() << "\n";

        //------------------------------------------------------------------------------//

        const int IdProcessOfGhostIP = eltOffProc.processId();
        const size_type idFaceInPartitionIP = face_it->idInOthersPartitions( IdProcessOfGhostIP );
        int IdProcessOfGhost = IdProcessOfGhostIP;
        size_type idFaceInPartition = idFaceInPartitionIP;

        //------------------------------------------------------------------------------//
        // for each dof in face
        for ( uint16_type locDof = 0; locDof < nbFaceDof; ++locDof )
        {
            // check only component 0
            const size_type theglobdoftest = faceLocalToGlobal( face_it->id(),locDof, 0 ).template get<0>();
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

                boost::tie( IdProcessOfGhost, idFaceInPartition ) = Feel::detail::updateDofOnVertices<mesh_type>(mesh,*face_it, myRank, IdProcessOfGhost, idFaceInPartition, eltOnProc, pointGetLocDof,
                                                                                                           procRecvData );
            }
            else if ( nDim == 3 && locDof < (face_type::numVertices*fe_type::nDofPerVertex + face_type::numEdges*fe_type::nDofPerEdge) )
            {
                int locDofInEgde = locDof - face_type::numVertices*fe_type::nDofPerVertex;
                const int nDofPerEdgeTemp = mpl::if_<boost::is_same<mpl::int_<fe_type::nDofPerEdge>,mpl::int_<0> >,
                                                     mpl::int_<1>,
                                                     mpl::int_<fe_type::nDofPerEdge> >::type::value;
                int edgeGetLocDof = locDofInEgde / nDofPerEdgeTemp;

                boost::tie( IdProcessOfGhost, idFaceInPartition ) = Feel::detail::updateDofOnEdges<mesh_type>(mesh,*face_it, myRank, IdProcessOfGhost, idFaceInPartition, eltOnProc, edgeGetLocDof,
                                                                                                        procRecvData);
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
                    const size_type theglobdof = faceLocalToGlobal( face_it->id(),locDof,c ).template get<0>();
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
    const size_type mynDofWithoutGhost = this->M_n_localWithGhost_df[myRank] - nDofNotPresent;
    mpi::all_gather( this->worldComm(),
                     mynDofWithoutGhost,
                     this->M_n_localWithoutGhost_df );

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
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_size_type_value );
    this->M_mapGlobalClusterToGlobalProcess.resize( this->M_n_localWithoutGhost_df[myRank],invalid_size_type_value );
    //------------------------------------------------------------------------------//
    // add in map the dofs presents
    size_type firstGlobIndex = this->M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( size_type i=0; i< this->M_n_localWithGhost_df[myRank]; ++i )
    {
        if ( !dofIsGhost[i] )
        {
            this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
            this->M_mapGlobalClusterToGlobalProcess[nextGlobIndex-firstGlobIndex]=i;
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
                                                                       std::set<int> const& procRecvData )
{
    const int myRank = this->worldComm().rank();
    //--------------------------------------------------------------------------------------------------------//
    size_type nbFaceDof = invalid_size_type_value;
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
                    nodeDofToSend[0]=dofPoint( theglobdof ).template get<0>()[0];
                    if ( nRealDim>1 )
                        nodeDofToSend[1]=dofPoint( theglobdof ).template get<0>()[1];
                    if ( nRealDim>2 )
                        nodeDofToSend[2]=dofPoint( theglobdof ).template get<0>()[2];
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
                                                            << "paritioning data incorect "
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
                        const auto thedofPt = itdofpt->template get<0>();
                        if ( itdofpt->template get<2>() != comp ) continue;

                        DVLOG(3) << "[buildGhostInterProcessDofMap] (myRank:" <<  myRank << ") "
                                 << "thedofPt: " << thedofPt << "nodeDofRecv: " << nodeDofRecv << "\n";

                        // test equatlity of dofs point
                        bool find2=true;
                        for (uint16_type d=0;d<nRealDim;++d)
                        {
                            find2 = find2 && (std::abs( thedofPt[d]-nodeDofRecv[d] )<1e-9);
                        }
                        // if find else save local dof
                        if (find2) { locDof = itdofpt->template get<1>();find=true; }
                    }
                    // check
                    CHECK( find ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";
                    //------------------------------------------------------------------------------//
                    // get global dof
                    const auto dofGlobAsked = locDof;
                    // save response
                    resAskedWithMultiProcess[cptDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                }
                else
                {
                for ( uint16_type l = 0; ( l < nbFaceDof && !find ) ; ++l )
                {
                    // dof point in face
                    const auto thedofPtInFace = dofPoint(faceLocalToGlobal( idFaceInMyPartition, l, comp ).template get<0>()).template get<0>();
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
                const auto dofGlobAsked = thedof.template get<0>();
                // save response
                resAskedWithMultiProcess[cptDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
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
                                                                                                        std::set<int> const& procRecvData )
{
    typedef std::vector< boost::tuple<uint16_type, ublas::vector<double> > > dofs_in_face_subcontainer_type;
    typedef boost::tuple<size_type, dofs_in_face_subcontainer_type > dofs_in_face_container_type;
    typedef std::vector< dofs_in_face_container_type > dofs_container_to_send_type;

    const int myRank = this->worldComm().localRank();
    const int nProc = this->worldComm().localSize();
    //--------------------------------------------------------------------------------------------------------//
    size_type nbFaceDof = invalid_size_type_value;
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
    std::map< int, int > nDataInVecToSend;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        const int nFaceToSend = listToSend[proc].size();
        if ( nFaceToSend == 0 ) continue;
        nDataInVecToSend[proc] = nFaceToSend;
    }
    //--------------------------------------------------------------------------------------------------------//
    // init and resize the container to send
    std::map< int,  dofs_container_to_send_type> dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const int idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //--------------------------------------------------------------------------------------------------------//
    // prepare container to send
    std::map< int, std::vector< std::vector<size_type> > > memoryInitialRequest;
    std::map< int, int > nDataInVecToSendBis;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
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
                    // get info to send
                    ublas::vector<double> nodeDofToSend( nRealDim );
                    nodeDofToSend[0]=dofPoint( theglobdof ).template get<0>()[0];
                    if ( nRealDim>1 )
                        nodeDofToSend[1]=dofPoint( theglobdof ).template get<0>()[1];
                    if ( nRealDim>2 )
                        nodeDofToSend[2]=dofPoint( theglobdof ).template get<0>()[2];
                    // up container
                    if (!componentsAreSamePoint)
                        dofsInFaceContainer[cptDof] = boost::make_tuple(comp,nodeDofToSend);
                    else if (comp==0)
                        dofsInFaceContainer[cptDof2] = boost::make_tuple(0,nodeDofToSend);

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
    for ( int proc=0; proc<nProc; ++proc )
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
    std::map<int,dofs_container_to_send_type> dataToRecv;
    auto itProcRecvData = procRecvData.begin();
    auto const enProcRecvData = procRecvData.end();
    for ( ; itProcRecvData!=enProcRecvData ; ++itProcRecvData )
    {
        const int proc = *itProcRecvData;
        reqs[cptRequest] = this->worldComm().localComm().irecv( proc , 0, dataToRecv[proc] );
        ++cptRequest;
    }
    //--------------------------------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //--------------------------------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<int, std::vector< std::vector<size_type> > > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const int idProc = itDataRecv->first;
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
            dataToReSend[idProc][cptFace].resize( nDofInFaceOpt,invalid_size_type_value );
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
                    const auto thedofPtInFace = dofPoint(faceLocalToGlobal( idFaceInMyPartition, l, comp ).template get<0>()).template get<0>();
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
                    const auto dofGlobAsked = thedof.template get<0>();
                    // save response
                    dataToReSend[idProc][cptFace][cptDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                }
                else
                {
                    for (uint16_type comp2=0; comp2<ncdof ; ++comp2)
                    {
                        const auto thedof = faceLocalToGlobal( idFaceInMyPartition, locDof, comp2 );
                        const size_type dofGlobAsked = thedof.template get<0>();
                        const int indexDofInFace = comp2*nDofInFace + cptDofInFace;
                        // save response
                        dataToReSend[idProc][cptFace][indexDofInFace] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
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
    std::map<int, std::vector<std::vector<size_type> > > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const int idProc = itDataToSend->first;
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
        const int idProc = itFinalDataToRecv->first;
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
    const int myRank = this->worldComm().rank();
    //------------------------------------------------------------------------------//
    // update datamap info
    this->M_n_dofs=0;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        this->M_n_localWithoutGhost_df[proc] = this->M_n_localWithGhost_df[proc];
        this->M_n_dofs+=this->M_n_localWithoutGhost_df[proc];
    }

    this->M_first_df_globalcluster[0]=0;
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
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_size_type_value );
    this->M_mapGlobalClusterToGlobalProcess.resize( this->M_n_localWithoutGhost_df[myRank],invalid_size_type_value );
    //------------------------------------------------------------------------------//
    size_type firstGlobIndex = this->M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( size_type i=0; i< this->M_n_localWithGhost_df[myRank]; ++i )
    {
        this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
        this->M_mapGlobalClusterToGlobalProcess[nextGlobIndex-firstGlobIndex]=i;
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
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostDofMapExtended( mesh_type& mesh )
{
    DVLOG(2) << "[buildGhostDofMap] call buildGhostDofMapExtended on rank "<<  this->worldComm().rank() << "\n";

    const int myRank = this->worldComm().localRank();
    const int nProc = this->worldComm().localSize();
    const int ncdof  = is_product?nComponents:1;
    const bool componentsAreSamePoint=true;

    size_type start_next_free_dof = this->M_last_df[myRank]+1;
    //------------------------------------------------------------------------------//
    // build extended dof table
    size_type next_free_dof = start_next_free_dof;
    DofFromElement<self_type,fe_type> dfe( this, *M_fe );
    auto face_it = mesh.interProcessFaces().first;
    auto const face_en = mesh.interProcessFaces().second;
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        if ( !this->isElementDone( eltOffProc.id() ) )
        {
            dfe.add( eltOffProc, next_free_dof, myRank );
        }
    }
    //------------------------------------------------------------------------------//
    // update local datamap
    bool hasNoElt = this->M_n_localWithGhost_df[myRank] == 0;
    const size_type thelastDof = ( !hasNoElt )?next_free_dof-1:0;

    std::vector<size_type> dataRecvFromGather;
    mpi::all_gather( this->worldComm().localComm(),
                     thelastDof,
                     dataRecvFromGather );

    for (int p=0;p<nProc;++p)
    {
        this->M_last_df[p] = dataRecvFromGather[p];

        size_type mynDofWithGhost = ( this->M_n_localWithGhost_df[p]>0 )?
            this->M_last_df[p] - this->M_first_df[p] + 1 : 0;
        this->M_n_localWithGhost_df[p] = mynDofWithGhost;
    }

    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_size_type_value );

    //------------------------------------------------------------------------------//
    // update doftable for face
    DofFromBoundary<self_type, fe_type> dfb( this, *M_fe );
    face_it = mesh.interProcessFaces().first;
    std::set<size_type> faceGhostDone;
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        for ( size_type f = 0; f < mesh.numLocalFaces(); f++ )
        {
            auto const& theface = eltOffProc.face(f);
            if ( theface.isGhostCell() && faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
            {
                auto faceIt = mesh.faceIterator( theface.id() );
                dfb.add( faceIt );
                faceGhostDone.insert( theface.id() );
            }
        }
    }
    //------------------------------------------------------------------------------//
    // update M_locglob_indices and M_locglob_signs
    face_it = mesh.interProcessFaces().first;
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        size_type elid = eltOffProc.id();
        if ( is_mortar && eltOffProc.isOnBoundary() )
        {
            VLOG(1) << "resizing indices and signs for mortar...";
            auto const& ldof = this->localDof( elid );
            size_type ne = std::distance( ldof.first, ldof.second );
            VLOG(1) << "resizing indices and signs for mortar:  " << ne;
            M_locglob_indices[elid].resize( ne );
            M_locglob_signs[elid].resize( ne );
            for( auto const& dof: this->localDof( elid ) )
            {
                M_locglob_indices[elid][dof.first.localDof()] = dof.second.index();
                M_locglob_signs[elid][dof.first.localDof()] = dof.second.sign();
            }
        }
        else
            for ( int i = 0; i < FEType::nLocalDof; ++i )
            {
                int nc1 = ( is_product?nComponents1:1 );

                for ( int c1 =0; c1 < nc1; ++c1 )
                {
                    int ind = FEType::nLocalDof*c1+i;
                    auto const& dof = localToGlobal( elid, i, c1 );
                    M_locglob_indices[elid][ind] = dof.index();
                    M_locglob_signs[elid][ind] = dof.sign();
                }
            }
    }
    //------------------------------------------------------------------------------//
    // generate dof point on extended part
    this->generateDofPointsExtendedGhostMap( mesh);

    //------------------------------------------------------------------------------//
    //
    std::vector< std::map<size_type, std::set< std::vector<size_type>  > > > eltIdToSend(this->worldComm().localSize());
    face_it = mesh.interProcessFaces().first;
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;
        const int processIdOfGhost = eltOffProc.processId();
        const size_type eltIdOfGhost = eltOffProc.id();
        const size_type eltIdInOtherPartOfGhost = eltOffProc.idInOthersPartitions(processIdOfGhost);
        for ( uint16_type l =0; l < fe_type::nLocalDof; ++l )
        {
            std::vector<size_type> compglobdofs( ncdof );
            bool addThisDof=true;
            for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
            {
                const size_type thedof = boost::get<0>( localToGlobal( eltIdOfGhost, l, c1 ) );
                // the optimisation below not work because we can lose the symterie of comm
                //if ( thedof < start_next_free_dof ) { addThisDof=false; break; }

                compglobdofs[c1]=thedof;
            }
            if ( addThisDof )
                eltIdToSend[processIdOfGhost][eltIdInOtherPartOfGhost].insert( compglobdofs );
        }

    }

    //------------------------------------------------------------------------------//
    // compute size of container to send
    std::map< int, int > nDataInVecToSend;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        const int nDofToSend = eltIdToSend[proc].size();
        if ( nDofToSend == 0 ) continue;
        nDataInVecToSend[proc] = nDofToSend;
    }
    //------------------------------------------------------------------------------//
    // prepare container to send
    typedef std::vector< boost::tuple< size_type, std::vector< ublas::vector<double> > > > dofs_container_to_send_type;
    std::map< int, dofs_container_to_send_type > dataToSend;
    std::map< int, std::vector< std::vector<size_type> > > memoryInitialRequest;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const int idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );

        auto itEltGhost = eltIdToSend[idProc].begin();
        auto const enEltGhost = eltIdToSend[idProc].end();
        memoryInitialRequest[idProc].resize(std::distance(itEltGhost,enEltGhost));
        for ( size_type cptElt = 0 ; itEltGhost!=enEltGhost;++itEltGhost,++cptElt)
        {
            const size_type idElt = itEltGhost->first;
            auto itDof = itEltGhost->second.begin();
            auto const enDof = itEltGhost->second.end();
            std::vector< ublas::vector<double> > nodesToSend( std::distance(itDof,enDof) );
            memoryInitialRequest[idProc][cptElt].resize( ncdof*std::distance(itDof,enDof) );
            for ( uint16_type cptDof=0 ; itDof!=enDof ; ++itDof, ++cptDof )
            {
                for ( uint16_type comp=0; comp<ncdof ; ++comp )
                {
                    const size_type theglobdof = itDof->operator[](comp);
                    const int indexDof = (componentsAreSamePoint)? cptDof*ncdof+comp : cptDof;
                    memoryInitialRequest[idProc][cptElt][indexDof] = theglobdof;
                }
                const size_type globDofUsedForPoint = itDof->operator[](0);

                nodesToSend[cptDof].resize( nRealDim );
                nodesToSend[cptDof][0]=dofPoint( globDofUsedForPoint ).template get<0>()[0];
                if ( nRealDim>1 )
                    nodesToSend[cptDof][1]=dofPoint( globDofUsedForPoint ).template get<0>()[1];
                if ( nRealDim>2 )
                    nodesToSend[cptDof][2]=dofPoint( globDofUsedForPoint ).template get<0>()[2];
            }

            dataToSend[idProc][cptElt] = boost::make_tuple( idElt,nodesToSend );
        }
    }

    //------------------------------------------------------------------------------//
    // counter of request
    int nbRequest=0;
    for ( int proc=0; proc<nProc; ++proc )
    {
        if ( dataToSend.find(proc) != dataToSend.end() )
            nbRequest+=2;
    }
    if ( nbRequest == 0 ) return;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    //------------------------------------------------------------------------------//
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToSend->first , 0, itDataToSend->second );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------//
    // first recv
    std::map< int,dofs_container_to_send_type> dataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const int proc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( proc , 0, dataToRecv[proc] );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<int, std::vector< std::vector<size_type> > > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const int idProc = itDataRecv->first;
        auto itEltRecv = itDataRecv->second.begin();
        auto const enEltRecv = itDataRecv->second.end();
        dataToReSend[idProc].resize( std::distance(itEltRecv,enEltRecv) );
        for ( int cptElt=0 ; itEltRecv!=enEltRecv ; ++itEltRecv,++cptElt )
        {
            const size_type idEltToSearch = itEltRecv->template get<0>();
            auto itNodes = itEltRecv->template get<1>().begin();
            auto const enNodes = itEltRecv->template get<1>().end();
            int nDofInElt = ncdof*std::distance(itNodes,enNodes);
            dataToReSend[idProc][cptElt].resize( nDofInElt,invalid_size_type_value );
            for ( int cptDof=0; itNodes!=enNodes ; ++itNodes, cptDof+=ncdof )
            {
                auto nodeDofRecv = *itNodes;
                bool find=false;
                uint16_type locDof=0;
                for ( uint16_type l =0; l < fe_type::nLocalDof && !find; ++l )
                {
                    size_type thedof = boost::get<0>( localToGlobal( idEltToSearch, l, 0/*c1*/ ) );
                    const auto thedofPt = dofPoint( thedof ).template get<0>();

                    // test equatlity of dofs point
                    bool find2=true;
                    for (uint16_type d=0;d<nRealDim;++d)
                    {
                        find2 = find2 && (std::abs( thedofPt[d]-nodeDofRecv[d] )<1e-9);
                    }
                    // if find else save local dof
                    if (find2)
                    {
                        locDof = l;
                        find=true;
                    }
                } // for ( uint16_type l =0; l < fe_type::nLocalDof && !find; ++l )
                CHECK( find ) << "\nPROBLEM with parallel dof table construction : Dof point not find " << nodeDofRecv << "\n";

                for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                {
                    size_type dofGlobAsked = boost::get<0>( localToGlobal( idEltToSearch, locDof, c1 ) );
                    dataToReSend[idProc][cptElt][cptDof+c1] = this->M_mapGlobalProcessToGlobalCluster[dofGlobAsked];
                }
            } // for ( int cptDof=0; itNodes ... )
        }
    } // for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )

    //------------------------------------------------------------------------------//
    // send respond to the request
    cptRequest=0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend!=enDataToReSend ; ++itDataToReSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToReSend->first , 0, itDataToReSend->second );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------//
    // recv the initial request
    std::map<int, std::vector<std::vector<size_type> > > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const int idProc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc, 0, finalDataToRecv[idProc] );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //------------------------------------------------------------------------------//
    // update datamap for ghost dof
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const int idProc = itFinalDataToRecv->first;
        auto itEltRecv = itFinalDataToRecv->second.begin();
        auto const enEltRecv = itFinalDataToRecv->second.end();
        for ( int cptElt=0 ; itEltRecv!=enEltRecv ; ++itEltRecv, ++cptElt )
        {
            if (componentsAreSamePoint)
            {
                const int nDofsInElt = itEltRecv->size()/ncdof;
                for ( int cptDof=0 ; cptDof< nDofsInElt ; ++cptDof )
                {
                    for (uint16_type comp2=0; comp2<ncdof ; ++comp2)
                    {
                        const int myindexDof = cptDof*ncdof + comp2; //comp2*nDofsInElt + cptDof;
                        const size_type myGlobProcessDof = memoryInitialRequest[idProc][cptElt][myindexDof];
                        const size_type dofGlobRecv = itEltRecv->operator[](myindexDof);
                        //update data map
                        this->M_mapGlobalProcessToGlobalCluster[myGlobProcessDof] = dofGlobRecv;
                    }
                }
            }
            else
            {
                CHECK(false) << "not implement\n";
            }
        }
    }
    //------------------------------------------------------------------------------//
    // update M_locglobOnCluster_indices and M_locglobOnCluster_signs
    face_it = mesh.interProcessFaces().first;
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

        size_type elid = eltOffProc.id();

        for ( int i = 0; i < FEType::nLocalDof; ++i )
        {
            int nc1 = ( is_product?nComponents1:1 );

            for ( int c1 =0; c1 < nc1; ++c1 )
            {
                int ind = FEType::nLocalDof*c1+i;
                auto const& dof = localToGlobalOnCluster( elid, i, c1 );

                M_locglobOnCluster_indices[elid][ind] = dof.index();
                M_locglobOnCluster_signs[elid][ind] = dof.sign();
            }
        }
    }
    //------------------------------------------------------------------------------//

}

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::generateDofPointsExtendedGhostMap( mesh_type& mesh )
{
    if ( fe_type::is_modal )
        return;

    DVLOG(2) << "[Dof::generateDofPoints] generating dof coordinates\n";
    typedef typename gm_type::template Context<vm::POINT, element_type> gm_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;

    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    gm_ptrtype gm( new gm_type );
    fe_type fe;
    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, fe.points() ) );


    gm_context_ptrtype __c;

    std::vector<bool> dof_done( nLocalDofWithGhost(), false );
    M_dof_points.resize( nLocalDofWithGhost() );

    auto face_it = mesh.interProcessFaces().first;
    auto const face_en = mesh.interProcessFaces().second;
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

        if ( __c )
            __c->update( eltOffProc );
        else
            __c = gm_context_ptrtype( new gm_context_type( gm, eltOffProc, __geopc ) );

        for ( uint16_type l =0; l < fe_type::nLocalDof; ++l )
        {
            int ncdof  = is_product?nComponents:1;

            for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
            {
                size_type thedof = boost::get<0>( localToGlobal( eltOffProc.id(), l, c1 ) );
                //std::cout << "aaaa pt with thedof " << thedof << std::endl;

                if ( ( thedof >= firstDof() ) && ( thedof <= lastDof() ) )
                {
                    // get only the local dof
                    //size_type thedofonproc = thedof - firstDof();
                    thedof -= firstDof();
                    DCHECK( thedof < nLocalDofWithGhost() )
                        << "invalid local dof index "
                        <<  thedof << ", " << nLocalDofWithGhost() << "," << firstDof()  << ","
                        <<  lastDof() << "," << eltOffProc.id() << "," << l << "," <<  c1;

                    if ( dof_done[ thedof ] == false )
                    {
                        //std::cout << "add pt with thedof " << thedof << std::endl;
                        //M_dof_points[dof_id] = boost::make_tuple( thedof, __c->xReal( l ) );
                        M_dof_points[thedof] = boost::make_tuple( __c->xReal( l ), firstDof()+thedof, c1 );
                        dof_done[thedof] = true;
                    }
                }
            }
        }

    }

}



























//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostInterProcessDofMap( mesh_type& mesh,
        std::map<size_type,boost::tuple<size_type,size_type> > & mapInterProcessDof )
{
    //------------------------------------------------------------------------------//
    // if nbFaceDof == 0 else return
    size_type nbFaceDof = invalid_size_type_value;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );
    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    if ( nbFaceDof == 0 ) return;

    //------------------------------------------------------------------------------//
    // init dof list to search with interprocess dof
    std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > listToSend(this->worldComm().size());
    this->buildGhostInterProcessDofMapInit(mesh,listToSend);
    //------------------------------------------------------------------------------//
    // search real and ghost dofs and build mapping between them
    // the search is iterative beacause some dofs are connected over 2 proc
    // and we need to search all proc around these dofs
    std::vector< std::set<size_type > > memoryFace(this->worldComm().size());
    bool allDofAreFind =false;
    while (!allDofAreFind)
    {
        auto const& res = this->buildGhostInterProcessDofMapRecursive(mesh,listToSend,
                                                                      mapInterProcessDof,
                                                                      memoryFace);

        allDofAreFind = res.template get<0>();
        if (!allDofAreFind) listToSend = res.template get<1>();
    }
    //------------------------------------------------------------------------------//
}

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//


template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostInterProcessDofMapInit( mesh_type& mesh,
                                                                               std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > & listToSend )
{

    // goal init container listToSend
    // std::vector< std::map<size_type,std::set<size_type> > > ( proc,( idFace,(globDof,..)), ...   )) )
    const int myRank = this->worldComm().rank();
    //------------------------------------------------------------------------------//
    // get nbFaceDof
    size_type nbFaceDof = invalid_size_type_value;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );
    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    if ( nbFaceDof == 0 ) return;
    //------------------------------------------------------------------------------//
    std::vector<bool> dofdone(this->M_n_localWithGhost_df[myRank],false);
    // iteration on all interprocessfaces in order to send requests to the near proc
    auto face_it = mesh.interProcessFaces().first;
    auto const face_en = mesh.interProcessFaces().second;
    DVLOG(2) << "[buildGhostInterProcessDofMap] nb interprocess faces: " << std::distance( face_it, face_en ) << "\n";
    for ( ; face_it != face_en ; ++face_it )
    {
        DVLOG(2) << "[buildGhostInterProcessDofMap] face id: " << face_it->id() << "\n";
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
        const bool elt0isGhost = elt0.isGhostCell();
        auto const& eltOnProc = (elt0isGhost)?elt1:elt0;
        auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

#if 0
        element_type eltOnProc,eltOffProc;
        if ( elt0.processId()!=myRank )
        {
            eltOffProc = elt0;
            eltOnProc=elt1;
        }
        else if ( elt1.processId()!=myRank )
        {
            eltOffProc = elt1;
            eltOnProc=elt0;
        }
        else std::cout << "\n[buildGhostInterProcessDofMapInit] PROBLEME1!!!!" << std::endl;
#endif
        //------------------------------------------------------------------------------//

        const int IdProcessOfGhost = eltOffProc.processId();
        const size_type idFaceInPartition = face_it->idInOthersPartitions( IdProcessOfGhost );
        //------------------------------------------------------------------------------//
        // for each dof in face
        const uint16_type ncdof = is_product?nComponents:1;
        // subcontainer
        std::set<boost::tuple<size_type,uint16_type> > dofSetAdd;

        for ( uint16_type l = 0; l < nbFaceDof; ++l )
        {
            for ( uint16_type c = 0; c < ncdof; ++c )
            {
                // add dof in subcontainer
                const size_type theglobdof = faceLocalToGlobal( face_it->id(),l,c ).template get<0>();
                if (!dofdone[theglobdof])
                {
                    dofSetAdd.insert(boost::make_tuple(theglobdof,c));
                    dofdone[theglobdof]=true;
                }
            }
            // add composante dofs in result container
        }
        listToSend[IdProcessOfGhost].insert( std::make_pair( idFaceInPartition/*face_it->id()*/, dofSetAdd  ) );


    } // for ( ; face_it != face_en ; ++face_it )
}

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

namespace detail {

template<typename element_type>
void
searchPartitionAroundNode( const int myRank, ublas::vector<double> const& nodeSearched,
                           element_type const& eltOnProc, std::set<size_type> & memoryIdsFaces,
                           std::set<boost::tuple<int,size_type> > & multiProcessRes )
{
    for ( uint16_type f = 0; f < element_type::numTopologicalFaces; ++f )
        {
            auto const& theFace = eltOnProc.face(f);
            size_type theFaceId = theFace.id();
            if (memoryIdsFaces.find(theFaceId) != memoryIdsFaces.end()) continue;

            // save face
            memoryIdsFaces.insert( theFaceId );
            // face contains node?
            bool find=false;
            auto const& facevertices = theFace.vertices();
            for (uint16_type v=0;v<element_type::face_type::numVertices && !find ;++v)
                {
                    bool find2=true;
                    for (uint16_type d=0;d<element_type::nDim;++d)
                        {
                            find2 = find2 && ( std::abs( facevertices(d,v)-nodeSearched[d] )<1e-9 );
                        }
                    if (find2) find=true;
                }

            if (find)
                {
                    if ( theFace.isConnectedTo0() )
                        {
                            auto const& eltOnProc0 = theFace.element0();
                            if (eltOnProc0.processId()!=myRank)
                                multiProcessRes.insert( boost::make_tuple(eltOnProc0.processId(),theFace.idInOthersPartitions(eltOnProc0.processId()) ));
                            else if (eltOnProc.id()!=eltOnProc0.id())
                                searchPartitionAroundNode(myRank,nodeSearched,eltOnProc0,memoryIdsFaces,multiProcessRes );
                        }
                    if ( theFace.isConnectedTo1() )
                        {
                            auto const& eltOnProc1 = theFace.element1();
                            if (eltOnProc1.processId()!=myRank)
                                multiProcessRes.insert( boost::make_tuple(eltOnProc1.processId(),theFace.idInOthersPartitions(eltOnProc1.processId()) ));
                            else if (eltOnProc.id()!=eltOnProc1.id())
                                searchPartitionAroundNode(myRank,nodeSearched,eltOnProc1,memoryIdsFaces,multiProcessRes );
                        }
                } // if (find)

        } // for ( uint16_type f = 0; f < element_type::numTopologicalFaces; ++f )

} // searchPartitionAroundNode


template<typename mesh_type>
void
searchPartitionAroundEdge( const int myRank,
                           typename mesh_type::face_type const& faceContainEdge,const uint16_type idEdgesInFace,
                           typename mesh_type::element_type const& eltOnProc, std::set<size_type> & memoryIdsFaces,
                           std::set<boost::tuple<int,size_type> > & multiProcessRes, bool firstExp, mpl::int_<0> /**/ )
{}

template<typename mesh_type>
void
searchPartitionAroundEdge( const int myRank,
                           typename mesh_type::face_type const& faceContainEdge,const uint16_type idEdgesInFace,
                           typename mesh_type::element_type const& eltOnProc, std::set<size_type> & memoryIdsFaces,
                           std::set<boost::tuple<int,size_type> > & multiProcessRes, bool firstExp, mpl::int_<1> /**/ )
{}

template<typename mesh_type>
void
searchPartitionAroundEdge( const int myRank,
                           typename mesh_type::face_type const& faceContainEdge,const uint16_type idEdgesInFace,
                           typename mesh_type::element_type const& eltOnProc, std::set<size_type> & memoryIdsFaces,
                           std::set<boost::tuple<int,size_type> > & multiProcessRes, bool firstExp, mpl::int_<2> /**/ )
{}

template<typename mesh_type>
void
searchPartitionAroundEdge( const int myRank,
                           typename mesh_type::face_type const& faceContainEdge,const uint16_type idEdgesInFace,
                           typename mesh_type::element_type const& eltOnProc, std::set<size_type> & memoryIdsFaces,
                           std::set<boost::tuple<int,size_type> > & multiProcessRes, bool firstExp, mpl::int_<3> /**/ )
{

    // only usefull whith first pass
    uint16_type iFaEl=0,iEdEl=0;
    if (firstExp)
    {
        iFaEl = ( faceContainEdge.processId() == faceContainEdge.proc_first() )? faceContainEdge.pos_first():faceContainEdge.pos_second();
        //local edge number (in element)
        iEdEl = mesh_type::element_type::fToE(  iFaEl, idEdgesInFace );
    }

    auto const& theedge = (firstExp)?eltOnProc.edge(iEdEl):faceContainEdge.edge(idEdgesInFace);
    //auto const& theEdgesVertices = theedge.vertices();

    for ( uint16_type f = 0; f < mesh_type::element_type::numTopologicalFaces; ++f )
    {
        auto const& theFace = eltOnProc.face(f);

        size_type theFaceId = theFace.id();
        //auto const& theneighborface = mesh.face( eltOnProc.face(f).id() );

        if (memoryIdsFaces.find(theFaceId) != memoryIdsFaces.end()) continue;

        // save face
        memoryIdsFaces.insert( theFaceId );

        // face contains edge?
        bool find=false;uint16_type newEdgeIdInFace=0;
        for (uint16_type ed=0;ed<mesh_type::element_type::face_type::numEdges && !find ;++ed)
        {
            //auto const& otheredges = theFace.edge(ed);
            if (theFace.edge(ed).id() == theedge.id())
            {
                find=true; newEdgeIdInFace=ed;
#if 0
                // check
                auto const& newEdgesVertices = theFace.edge(ed).vertices();
                //auto const& theEdgesVertices = theedge.vertices();
                bool find2=true;
                for (uint16_type d=0;d<mesh_type::element_type::nDim;++d)
                    {
                        find2 = find2 && ( std::abs( newEdgesVertices(d,0)-theEdgesVertices(d,0) )<1e-9 );
                    }
                if (!find2) std::cout << "\n IAOIOIOAIOIO"<<std::endl;
                //if (newEdgesVertices!=theEdgesVertices) std::cout << "\IAOIOIOAIOIO"<<std::endl;
#endif
            }
        }

        if (find)
        {
            if ( theFace.isConnectedTo0() )
            {
                auto const& eltOnProc0 = theFace.element0();
                if (eltOnProc0.processId()!=myRank)
                    multiProcessRes.insert( boost::make_tuple(eltOnProc0.processId(),theFace.idInOthersPartitions(eltOnProc0.processId()) ));
                else if (eltOnProc.id()!=eltOnProc0.id())
                    searchPartitionAroundEdge<mesh_type>(myRank,theFace,newEdgeIdInFace,eltOnProc0,memoryIdsFaces,multiProcessRes,false,mpl::int_<3>());
            }
            if ( theFace.isConnectedTo1() )
            {
                auto const& eltOnProc1 = theFace.element1();
                if (eltOnProc1.processId()!=myRank)
                    multiProcessRes.insert( boost::make_tuple(eltOnProc1.processId(),theFace.idInOthersPartitions(eltOnProc1.processId()) ));
                else if (eltOnProc.id()!=eltOnProc1.id())
                    searchPartitionAroundEdge<mesh_type>(myRank,theFace,newEdgeIdInFace,eltOnProc1,memoryIdsFaces,multiProcessRes,false,mpl::int_<3>());
            }
        } // if (find)

    } // for ( uint16_type f = 0; f < element_type::numTopologicalFaces; ++f )

}

template<typename mesh_type>
void
searchPartitionAroundEdge( const int myRank,
                           typename mesh_type::face_type const& faceContainEdge,const uint16_type idEdgesInFace,
                           typename mesh_type::element_type const& eltOnProc, std::set<size_type> & memoryIdsFaces,
                           std::set<boost::tuple<int,size_type> > & multiProcessRes )
{
    searchPartitionAroundEdge<mesh_type>(myRank,faceContainEdge,idEdgesInFace,eltOnProc,memoryIdsFaces,multiProcessRes,true,mpl::int_<mesh_type::nDim>());
}

} // namespace detail

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
boost::tuple<bool, std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > >
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGhostInterProcessDofMapRecursive( mesh_type& mesh,
                                                                                    std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > const& listToSend,
                                                                                    std::map<size_type,boost::tuple<size_type,size_type> > & mapInterProcessDof,
                                                                                    std::vector< std::set<size_type > > & memoryFace )
{
    const int myRank = this->worldComm().rank();
    //--------------------------------------------------------------------------------------------------------//
    size_type nbFaceDof = invalid_size_type_value;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;
    //--------------------------------------------------------------------------------------------------------//

    std::vector<int> nbMsgToSend( this->worldComm().size(), 0 );

    std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > listNeedMoreSearch(this->worldComm().size());
    int counterListNeedMoreSearch = 0, counterListNeedMoreSearchLoc = 0;

    std::vector< std::vector< std::vector<boost::tuple<uint16_type,int> > > > memoryInitialRequest( this->worldComm().size() );

    typedef std::vector< boost::tuple<uint16_type, ublas::vector<double> > > dofs_in_face_subcontainer_type;
    typedef boost::tuple<size_type, dofs_in_face_subcontainer_type > dofs_in_face_container_type;

    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        auto itFaces = listToSend[proc].begin();
        auto const enFaces = listToSend[proc].end();
        memoryInitialRequest[proc].resize(std::distance(itFaces,enFaces));
        for ( int cptFaces=0 ; itFaces!=enFaces ; ++itFaces, ++cptFaces)
        {
            auto itDof = itFaces->second.begin();
            auto const enDof = itFaces->second.end();
            auto const nDofsInFace = std::distance(itDof,enDof);
            dofs_in_face_subcontainer_type dofsInFaceContainer(nDofsInFace);
            memoryInitialRequest[proc][cptFaces].resize(nDofsInFace);
            for (int cptDof=0 ; itDof!=enDof ; ++itDof,++cptDof)
            {
                auto const theglobdof = itDof->get<0>();
                auto const comp = itDof->get<1>();
                // save the tag of mpi send
                memoryInitialRequest[proc][cptFaces][cptDof] = boost::make_tuple(comp,theglobdof);
                //------------------------------------------------------------------------------//
                // get info to send
                ublas::vector<double> nodeDofToSend( nDim );
                nodeDofToSend[0]=dofPoint( theglobdof ).template get<0>()[0];
                if ( nDim>1 )
                    nodeDofToSend[1]=dofPoint( theglobdof ).template get<0>()[1];
                if ( nDim>2 )
                    nodeDofToSend[2]=dofPoint( theglobdof ).template get<0>()[2];
                // up container
                dofsInFaceContainer[cptDof] = boost::make_tuple(comp,nodeDofToSend);
                //------------------------------------------------------------------------------//
            }

            this->worldComm().send( proc , nbMsgToSend[proc], boost::make_tuple(itFaces->first,dofsInFaceContainer) );
            ++nbMsgToSend[proc];
        } // for ( int cptFaces=0 ; itFaces!=enFaces ; ++itFaces, ++cptFaces)
    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )

    //--------------------------------------------------------------------------------------------------------//
    // counter of msg received for each process
    std::vector<int> nbMsgToRecv;
    mpi::all_to_all( this->worldComm(),
                     nbMsgToSend,
                     nbMsgToRecv );
    //--------------------------------------------------------------------------------------------------------//
    // recv dof asked and re-send
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
        {
            dofs_in_face_container_type dataToRecvVec;
            this->worldComm().recv( proc, cpt, dataToRecvVec );

            auto const idFaceInMyPartition = dataToRecvVec.template get<0>();

            auto itDofInFace = dataToRecvVec.template get<1>().begin();
            auto const enDofInFace = dataToRecvVec.template get<1>().end();
            std::vector< boost::tuple<size_type,std::vector<boost::tuple<int,size_type> > > > resAskedWithMultiProcess(std::distance(itDofInFace,enDofInFace));
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

#if 0
                element_type eltOnProc, eltOffProc;
                uint16_type faceIdInEltOnProc = invalid_uint16_type_value;
                if ( elt0.processId()!=myRank )
                {
                    eltOffProc = elt0;
                    eltOnProc=elt1;
                    faceIdInEltOnProc = theface.pos_second();
                }
                else if ( elt1.processId()!=myRank )
                {
                    eltOffProc = elt1;
                    eltOnProc=elt0;
                    faceIdInEltOnProc = theface.pos_first();
                }
                else
                {
                    CHECK( (elt0.processId()==myRank) ||
                           (elt1.processId()==myRank) )
                        << "\nPROBLEM with parallel dof table construction\n"
                        << " elt0.processId() " << elt0.processId()
                        << " elt1.processId() " << elt1.processId()
                        << "\n";
                }
#endif

                //------------------------------------------------------------------------------//
                // search dof on face recv
                int locDof = nbFaceDof;
                bool find=false;

                for ( uint16_type l = 0; ( l < nbFaceDof && !find ) ; ++l )
                {
                    // dof point in face
                    const auto thedofPtInFace = dofPoint(faceLocalToGlobal( idFaceInMyPartition, l, comp ).template get<0>()).template get<0>();
                    // test equatlity of dofs point
                    bool find2=true;
                    for (uint16_type d=0;d<nDim;++d)
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
                //if ( !find ) std::cout<<"\n [buildGhostInterProcessDofMap] : Dof point not find on interprocess face " << nodeDofRecv << std::endl;
                CHECK( find ) << "\nPROBLEM with parallel dof table construction : Dof point not find on interprocess face " << nodeDofRecv << "\n";
                //------------------------------------------------------------------------------//
                // get global dof
                const auto thedof = faceLocalToGlobal( idFaceInMyPartition, locDof, comp );
                const auto dofGlobAsked = thedof.template get<0>();
                //------------------------------------------------------------------------------//
                //------------------------------------------------------------------------------//
                //------------------------------------------------------------------------------//
                //------------------------------------------------------------------------------//
                // search maybe with neighboors elt
                std::set<size_type> memoryIdsFaces;
                std::set<boost::tuple<int,size_type> > multiProcessRes;
                if ( locDof < face_type::numVertices*fe_type::nDofPerVertex)
                {
                        Feel::detail::searchPartitionAroundNode( myRank, nodeDofRecv, eltOnProc, memoryIdsFaces, multiProcessRes );
                }
                else if ( nDim == 3 && locDof < (face_type::numVertices*fe_type::nDofPerVertex + face_type::numEdges*fe_type::nDofPerEdge) )
                {
                    int locDofInEgde = locDof - face_type::numVertices*fe_type::nDofPerVertex;

                    const int nDofPerEdgeTemp = mpl::if_<boost::is_same<mpl::int_<fe_type::nDofPerEdge>,mpl::int_<0> >,
                                                         mpl::int_<1>,
                                                         mpl::int_<fe_type::nDofPerEdge> >::type::value;

                    int edgeGetLocDof = locDofInEgde / nDofPerEdgeTemp;

                    Feel::detail::searchPartitionAroundEdge<mesh_type>( myRank, theface, edgeGetLocDof, eltOnProc, memoryIdsFaces, multiProcessRes );
                }
                //if (proc==0 ) std::cout << "myRank " << myRank << " multiProcessRes.size() " << multiProcessRes.size() << " " << nodeDofRecv << std::endl;
                //------------------------------------------------------------------------------//
                // convert in vector for send
                auto const nbMultiProcessRes = multiProcessRes.size();
                std::vector<boost::tuple<int,size_type> > multiProcessResVec( nbMultiProcessRes );
                std::copy( multiProcessRes.begin(),multiProcessRes.end(),
                           multiProcessResVec.begin() );
                //------------------------------------------------------------------------------//
                //------------------------------------------------------------------------------//
                //------------------------------------------------------------------------------//
                // send response
                resAskedWithMultiProcess[cptDofInFace] = boost::make_tuple(dofGlobAsked,multiProcessResVec);
                //------------------------------------------------------------------------------//
            }
            this->worldComm().send( proc, cpt, resAskedWithMultiProcess );

        } // for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )

    //--------------------------------------------------------------------------------------------------------//
    // get response to initial request and update Feel::Mesh::Faces data
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        {
            //------------------------------------------------------------------------------//
            // recv response
            std::vector< boost::tuple<size_type,std::vector<boost::tuple<int,size_type> > > > resultRecvWithMultiProcess;
            this->worldComm().recv( proc, cpt, resultRecvWithMultiProcess );
            //------------------------------------------------------------------------------//
            // iterate on dofs
            auto itDofRes = resultRecvWithMultiProcess.begin();
            auto const enDofRes = resultRecvWithMultiProcess.end();
            for ( int cptDofRes=0 ; itDofRes != enDofRes ; ++itDofRes,++cptDofRes )
            {

                auto const mycomp = memoryInitialRequest[proc][cpt][cptDofRes].template get<0>();
                auto const myGlobProcessDof = memoryInitialRequest[proc][cpt][cptDofRes].template get<1>();
                const auto dofGlobRecv = itDofRes->template get<0>();
                const auto multiProcessRes = itDofRes->template get<1>();

                //------------------------------------------------------------------------------//
                // update mapInterProcessDof data
                if ( mapInterProcessDof.find( myGlobProcessDof ) == mapInterProcessDof.end() )
                {
                    mapInterProcessDof.insert( std::make_pair( myGlobProcessDof, boost::make_tuple( proc,dofGlobRecv ) ) );
                }
                else if ( ( int )mapInterProcessDof.find( myGlobProcessDof )->second.template get<0>() > proc )
                {
                    mapInterProcessDof[ myGlobProcessDof ] = boost::make_tuple( proc,dofGlobRecv );
                }
                //------------------------------------------------------------------------------//
                // save global dof with search on proc
                memoryFace[proc].insert(myGlobProcessDof);
                //------------------------------------------------------------------------------//
                // add maybe more search
                auto itMP = multiProcessRes.begin();
                auto const enMP = multiProcessRes.end();
                for ( ; itMP!=enMP ;++itMP )
                {
                    auto const eltMultiProcessProcId = itMP->template get<0>();
                    auto const eltMultiProcessFaceId = itMP->template get<1>();
                    if ( memoryFace[eltMultiProcessProcId].find( myGlobProcessDof ) == memoryFace[eltMultiProcessProcId].end() )
                    {
                        listNeedMoreSearch[eltMultiProcessProcId][eltMultiProcessFaceId].insert(boost::make_tuple(myGlobProcessDof,mycomp));
                        ++counterListNeedMoreSearchLoc;
                    }
                }
            }
        } // for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )

    //--------------------------------------------------------------------------------------------------------//
    // need to apply an all_reduce
    mpi::all_reduce( this->worldComm().globalComm(),
                     counterListNeedMoreSearchLoc,
                     counterListNeedMoreSearch,
                     std::plus<int>() );
    //std::cout<< " counterListNeedMoreSearch "<< counterListNeedMoreSearch << " counterListNeedMoreSearchLoc " << counterListNeedMoreSearchLoc << std::endl;
    //--------------------------------------------------------------------------------------------------------//
    return boost::make_tuple( (bool) (counterListNeedMoreSearch == 0),
                              listNeedMoreSearch);


} // buildGhostInterProcessDofMapRecursive

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildDofNotPresent( std::map<size_type,boost::tuple<size_type,size_type> > const & mapInterProcessDof,
                                                                 std::map<size_type,boost::tuple<size_type,size_type> > & setInterProcessDofNotPresent )

{
    const int myRank = this->worldComm().rank();

    auto it_mapIPDof = mapInterProcessDof.begin();
    auto const en_mapIPDof = mapInterProcessDof.end();
    for ( ; it_mapIPDof!=en_mapIPDof; ++it_mapIPDof )
    {
        if ( it_mapIPDof->second.template get<0>() < myRank )
            setInterProcessDofNotPresent[it_mapIPDof->first] = it_mapIPDof->second;
    }
}

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::buildGlobalProcessToGlobalClusterDofMap( mesh_type& mesh,
        std::map<size_type,boost::tuple<size_type,size_type> > const& setInterProcessDofNotPresent )
{
    const int myRank = this->worldComm().rank();

    //------------------------------------------------------------------------------//
    // update datamap info
    //const size_type mynDofWithoutGhost = this->M_last_df[myRank] - this->M_first_df[myRank] + 1 - setInterProcessDofNotPresent.size();
    const size_type mynDofWithoutGhost = this->M_n_localWithGhost_df[myRank] - setInterProcessDofNotPresent.size();
    mpi::all_gather( this->worldComm(),
                     mynDofWithoutGhost,
                     this->M_n_localWithoutGhost_df );

    this->M_n_dofs=0;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        this->M_n_dofs+=this->M_n_localWithoutGhost_df[proc];
    }

    this->M_first_df_globalcluster=this->M_first_df;
    this->M_last_df_globalcluster=this->M_last_df;
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
    this->M_mapGlobalProcessToGlobalCluster.resize( this->M_n_localWithGhost_df[myRank],invalid_size_type_value );
    this->M_mapGlobalClusterToGlobalProcess.resize( this->M_n_localWithoutGhost_df[myRank],invalid_size_type_value );

    //------------------------------------------------------------------------------//
    // add in map the dofs presents
    size_type firstGlobIndex = this->M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( size_type i=0; i< this->M_n_localWithGhost_df[myRank]; ++i )
    {
        if ( setInterProcessDofNotPresent.find( i )==setInterProcessDofNotPresent.end() )
        {
            this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
            this->M_mapGlobalClusterToGlobalProcess[nextGlobIndex-firstGlobIndex]=i;
            ++nextGlobIndex;
        }
    }

    //------------------------------------------------------------------------------//
    // update datamap with dofs non presents (ghosts)
    this->updateGhostGlobalDof( setInterProcessDofNotPresent );

} // buildGlobalProcessToGlobalClusterDofMap

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType, typename MortarType>
void
DofTable<MeshType, FEType, PeriodicityType, MortarType>::updateGhostGlobalDof( std::map<size_type,boost::tuple<size_type,size_type> > const& setInterProcessDofNotPresent )
{
    std::vector<int> nbMsgToSend2( this->worldComm().size(), 0 );
    std::vector< std::map<int,size_type> > mapMsg2( this->worldComm().size() );

    auto dofNP_it = setInterProcessDofNotPresent.begin();
    auto dofNP_en = setInterProcessDofNotPresent.end();
    for ( ; dofNP_it!=dofNP_en ; ++dofNP_it )
    {
        auto const globalProcessDof = dofNP_it->first;
        auto const IdProcessOfGhost = dofNP_it->second.template get<0>();
        auto const dofFaceInPartition = dofNP_it->second.template get<1>();

        this->worldComm().send( IdProcessOfGhost , nbMsgToSend2[IdProcessOfGhost], dofFaceInPartition );
        mapMsg2[IdProcessOfGhost].insert( std::make_pair( nbMsgToSend2[IdProcessOfGhost],globalProcessDof ) );
        ++nbMsgToSend2[IdProcessOfGhost];
    }

    // counter of msg received for each process
    std::vector<int> nbMsgToRecv2;
    mpi::all_to_all( this->worldComm(),
                     nbMsgToSend2,
                     nbMsgToRecv2 );

    // recv dof asked and re-send dof in this proc
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToRecv2[proc]; ++cpt )
        {
            //recv
            size_type dofRecv;
            this->worldComm().recv( proc, cpt, dofRecv );
            // send
            size_type const dofToSend= this->M_mapGlobalProcessToGlobalCluster[dofRecv];
            this->worldComm().send( proc, cpt,dofToSend );
        }
    }

    // get response to initial request and update Feel::Mesh::Faces data
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend2[proc]; ++cpt )
        {
            //recv
            size_type dofGlobClusterRecv;
            this->worldComm().recv( proc, cpt, dofGlobClusterRecv );
            //update data
            this->M_mapGlobalProcessToGlobalCluster[mapMsg2[proc][cpt]]=dofGlobClusterRecv;
        }
    }

}



} // namespace Feel

#endif /* FEELPP_DOFTABLE_MPI_HPP */
