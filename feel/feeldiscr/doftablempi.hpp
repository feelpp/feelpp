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

template <typename MeshType>
boost::tuple<rank_type,size_type >
updateDofOnVertices( MeshType const& mesh, typename MeshType::face_type const& theface, const rank_type myIdProcess,
                     const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                     const uint16_type locDof, std::set<rank_type> & procRecvData )
{
    rank_type procMin = IdProcessOfGhost;
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
            const rank_type theprocGhost=itprocghost->first;
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
            const rank_type procIdGhost = itprocghost->first;
            procRecvData.insert( procIdGhost );
        }
    }

    return boost::make_tuple(procMin,idFaceMin);
}

//--------------------------------------------------------------------------------------------------------//

template <typename MeshType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface,const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData, mpl::int_<0> /**/ )
{
    return boost::make_tuple(0,0);
}
template <typename MeshType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData, mpl::int_<1> /**/ )
{
    return boost::make_tuple(0,0);
}
template <typename MeshType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData, mpl::int_<2> /**/ )
{
    return boost::make_tuple(0,0);
}
template <typename MeshType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData, mpl::int_<3> /**/ )
{
    rank_type procMin = IdProcessOfGhost;
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
            const rank_type theprocGhost=itprocghost->first;

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
            const rank_type procIdGhost = itprocghost->first;
            procRecvData.insert( procIdGhost );
        }

    }

    return boost::make_tuple(procMin,idFaceMin);
}

//--------------------------------------------------------------------------------------------------------//

template <typename MeshType>
boost::tuple<rank_type,size_type>
updateDofOnEdges( MeshType const& mesh, typename MeshType::face_type const& theface, const rank_type myIdProcess,
                  const rank_type IdProcessOfGhost, const size_type idFaceInPartition, typename MeshType::element_type const& eltOnProc,
                  const uint16_type idEdgesInFace, std::set<rank_type> & procRecvData )
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
                                                                                                        std::set<rank_type> & procRecvData )
{
    // goal init container listToSend
    // std::vector< std::map<size_type,std::set<size_type> > > ( proc,( idFace,(globDof,..)), ...   )) )
    const rank_type myRank = this->worldComm().rank();
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

        const rank_type IdProcessOfGhostIP = eltOffProc.processId();
        const size_type idFaceInPartitionIP = face_it->idInOthersPartitions( IdProcessOfGhostIP );
        rank_type IdProcessOfGhost = IdProcessOfGhostIP;
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
                                                                       std::set<rank_type> const& procRecvData )
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
                    this->M_activeDofSharedOnCluster[dofGlobAsked].insert(proc);
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
                    this->M_activeDofSharedOnCluster[dofGlobAsked].insert(idProc);
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

    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();
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
                M_face_l2g[ faceIt->id()].resize( nLocalDofOnFace() );
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
    std::map< rank_type, int > nDataInVecToSend;
    for ( rank_type proc=0; proc<this->worldComm().size(); ++proc )
    {
        const int nDofToSend = eltIdToSend[proc].size();
        if ( nDofToSend == 0 ) continue;
        nDataInVecToSend[proc] = nDofToSend;
    }
    //------------------------------------------------------------------------------//
    // prepare container to send
    typedef std::vector< boost::tuple< size_type, std::vector< ublas::vector<double> > > > dofs_container_to_send_type;
    std::map< rank_type, dofs_container_to_send_type > dataToSend;
    std::map< rank_type, std::vector< std::vector<size_type> > > memoryInitialRequest;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
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
    for ( rank_type proc=0; proc<nProc; ++proc )
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
    std::map< rank_type,dofs_container_to_send_type> dataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type proc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( proc , 0, dataToRecv[proc] );
        ++cptRequest;
    }
    //------------------------------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //------------------------------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector< std::vector<size_type> > > dataToReSend;
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
                    this->M_activeDofSharedOnCluster[dofGlobAsked].insert(idProc);
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
    std::map<rank_type, std::vector<std::vector<size_type> > > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
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
        const rank_type idProc = itFinalDataToRecv->first;
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
DofTable<MeshType, FEType, PeriodicityType, MortarType>::generateDofPointsExtendedGhostMap( mesh_type& mesh ) const
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






} // namespace Feel

#endif /* FEELPP_DOFTABLE_MPI_HPP */
