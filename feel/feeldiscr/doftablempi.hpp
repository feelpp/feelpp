/*
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

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::buildGhostDofMap( mesh_type& mesh )
{
    // dofs map on interprocess faces : ( dofOnMyRank, tuple< otherProcRank, dofOnOtherRank> )
    std::map<size_type,boost::tuple<size_type,size_type> > mapInterProcessDof;
    // dofs no present on interprocess faces
    std::set<int> setInterProcessDofNotPresent;

    if (is_continuous)
    {
        DVLOG(2) << "[buildGhostDofMap] call buildGhostInterProcessDofMap() with god rank "<<  this->worldComm().godRank() << "\n";
        buildGhostInterProcessDofMap( mesh,mapInterProcessDof );

        DVLOG(2) << "[buildGhostDofMap] call buildDofNotPresent() with rank "<<  this->worldComm().rank() << "\n";
        buildDofNotPresent( mapInterProcessDof,setInterProcessDofNotPresent );
    }

    DVLOG(2) << "[buildGhostDofMap] call buildGlobalProcessToGlobalClusterDofMap() with rank "<<  this->worldComm().rank() << "\n";
    buildGlobalProcessToGlobalClusterDofMap( mesh,mapInterProcessDof,setInterProcessDofNotPresent );

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
                boost::tie( M_locglobOnCluster_indices[elid][ind],
                            M_locglobOnCluster_signs[elid][ind], boost::tuples::ignore ) =
                                localToGlobalOnCluster( elid, i, c1 );
            }
        }
    }

    DVLOG(2) << "[buildGhostDofMap] finish () with god rank "<< this->worldComm().godRank() << "\n";

}



//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::buildGhostInterProcessDofMapInit( mesh_type& mesh,
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
    // iteration on all interprocessfaces in order to send requests to the near proc
    auto face_it = mesh.interProcessFaces().first;
    auto face_en = mesh.interProcessFaces().second;
    DVLOG(2) << "[buildGhostInterProcessDofMap] nb interprocess faces: " << std::distance( face_it, face_en ) << "\n";
    for ( ; face_it != face_en ; ++face_it )
    {
        DVLOG(2) << "[buildGhostInterProcessDofMap] face id: " << face_it->id() << "\n";
        element_type eltOnProc,eltOffProc;
        auto const& elt0 = face_it->element0();
        auto const& elt1 = face_it->element1();
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
        //------------------------------------------------------------------------------//

        const int IdProcessOfGhost = eltOffProc.processId();
        const size_type idFaceInPartition = face_it->idInPartition( IdProcessOfGhost );
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
                dofSetAdd.insert(boost::make_tuple(theglobdof,c));
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

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::buildGhostInterProcessDofMap( mesh_type& mesh,
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

template<typename MeshType, typename FEType, typename PeriodicityType>
boost::tuple<bool, std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > >
DofTable<MeshType, FEType, PeriodicityType>::buildGhostInterProcessDofMapRecursive( mesh_type& mesh,
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

    std::vector<int> nbMsgToSend( this->worldComm().size() );
    std::fill( nbMsgToSend.begin(),nbMsgToSend.end(),0 );
    std::vector< std::map<int,int> > mapMsg( this->worldComm().size() );
    std::vector< std::map<int,uint16_type> > mapMsgComp( this->worldComm().size() );

    std::vector< std::map<size_type,std::set<boost::tuple<size_type,uint16_type> > > > listNeedMoreSearch(this->worldComm().size());
    int counterListNeedMoreSearch = 0, counterListNeedMoreSearchLoc = 0;

    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        //std::cout << "size list " << listToSend[proc].size() << std::endl;
        auto itFaces = listToSend[proc].begin();
        auto const enFaces = listToSend[proc].end();
        for ( ; itFaces!=enFaces ; ++itFaces)
        {
            auto itDof = itFaces->second.begin();
            auto const enDof = itFaces->second.end();
            for ( ; itDof!=enDof ; ++itDof)
            {
                auto const theglobdof = itDof->get<0>();
                auto const comp = itDof->get<1>();
                // save the tag of mpi send
                mapMsg[proc].insert( std::make_pair( nbMsgToSend[proc], theglobdof ) );
                mapMsgComp[proc].insert( std::make_pair( nbMsgToSend[proc], comp ) );
                //------------------------------------------------------------------------------//
                // get info to send
                //ublas::vector<double> dataToSend( 2+nDim );
                //dataToSend[0]=itFaces->first;//idFaceInPartition;
                //dataToSend[1]=comp;// composante
                ublas::vector<double> nodeDofToSend( nDim );
                nodeDofToSend[0]=dofPoint( theglobdof ).template get<0>()[0];
                if ( nDim>1 )
                    nodeDofToSend[1]=dofPoint( theglobdof ).template get<0>()[1];
                if ( nDim>2 )
                    nodeDofToSend[2]=dofPoint( theglobdof ).template get<0>()[2];
                //boost::tuple< size_type, uint16_type
                auto dataToSend= boost::make_tuple(itFaces->first,comp,nodeDofToSend);
                //------------------------------------------------------------------------------//
                // send
                this->worldComm().send( proc , nbMsgToSend[proc], dataToSend );
                //------------------------------------------------------------------------------//
                // update nb send
                ++nbMsgToSend[proc];
            }
        }
    }


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
            //recv
            //ublas::vector<double> dataToRecv( 2+nDim );
            //ublas::vector<double> dataToRecv( nDim );
            boost::tuple<size_type,uint16_type,ublas::vector<double> > dataToRecv;
            this->worldComm().recv( proc, cpt, dataToRecv );
            //------------------------------------------------------------------------------//

            auto const idFaceInMyPartition = dataToRecv.get<0>();
            auto const comp = dataToRecv.get<1>();
            auto const nodeDofRecv = dataToRecv.get<2>();

            auto const& theface = mesh.face( dataToRecv.get<0>() );
            auto const& elt0 = theface.element0();
            auto const& elt1 = theface.element1();

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
            if ( !find ) std::cout<<"\n [buildGhostInterProcessDofMap] : Dof point not find on interprocess face " << nodeDofRecv << std::endl;
            //------------------------------------------------------------------------------//
            // search maybe an another neighbor
            bool eltHasMultiProcess = false;
            int eltMultiProcessProcId = 0;
            size_type eltMultiProcessFaceId = 0;
            int eltOnProc0ProcId=0,eltOnProc1ProcId=0;
            for ( uint16_type f = 0; f < face_type::numTopologicalFaces; ++f )
            {
                if (f==faceIdInEltOnProc) continue;

                auto const& theneighborface = mesh.face( eltOnProc.face(f).id() );
                //auto const& eltOnProc0 = mesh.face( eltOnProc.face(f).id() ).element0();
                //auto const& eltOnProc1 = mesh.face( eltOnProc.face(f).id() ).element1();

                if ( theneighborface.isConnectedTo0() )
                    eltOnProc0ProcId = theneighborface.element0().processId();
                else eltOnProc0ProcId = myRank;

                if ( theneighborface.isConnectedTo1() )
                    eltOnProc1ProcId = theneighborface.element1().processId();
                else eltOnProc1ProcId = myRank;

                //if( (eltOnProc0.processId()!=myRank && eltOnProc0.processId()!=proc) ||
                //(eltOnProc1.processId()!=myRank && eltOnProc1.processId()!=proc) )
                if( (eltOnProc0ProcId!=myRank && eltOnProc0ProcId!=proc) ||
                    (eltOnProc1ProcId!=myRank && eltOnProc1ProcId!=proc) )
                {
                    auto const facevertices = eltOnProc.face(f).vertices();
                    for (uint16_type v=0;v<face_type::numVertices;++v)
                        if (true)//locDof<fe_type::nDofPerVertex) // dof on vertex
                    {
                        bool find2=true;
                        for (uint16_type d=0;d<nDim;++d)
                        {
                            find2 = find2 && ( std::abs( facevertices(d,v)-nodeDofRecv[d] )<1e-9 );
                        }

                        if (find2)
                        {
                            eltHasMultiProcess=true;
                            if (eltOnProc0ProcId!=myRank)
                                eltMultiProcessProcId = eltOnProc0ProcId;
                            else
                                eltMultiProcessProcId = eltOnProc1ProcId;

                            eltMultiProcessFaceId = eltOnProc.face(f).idInPartition( eltMultiProcessProcId );
                        }

                    } // if (locDof<fe_type::nDofPerVertex)
                } // if (locDof<fe_type::nDofPerVertex)
            } // for ( uint16_type f = 0; f < face_type::numTopologicalFaces; ++f )

            //------------------------------------------------------------------------------//
            // get global dof
            const auto thedof = faceLocalToGlobal( idFaceInMyPartition, locDof, comp );
            const int dofGlobAsked = thedof.template get<0>();
            //------------------------------------------------------------------------------//
            // prepare data to send
            boost::tuple<size_type,bool,int,size_type> eltAskedMultiProcess = boost::make_tuple(thedof.template get<0>(),
                                                                                                eltHasMultiProcess,
                                                                                                eltMultiProcessProcId,
                                                                                                eltMultiProcessFaceId );
            //------------------------------------------------------------------------------//
            // send response
            this->worldComm().send( proc, cpt, eltAskedMultiProcess );
            //------------------------------------------------------------------------------//
        } // for ( int cpt=0; cpt<nbMsgToRecv[proc]; ++cpt )
    } // for ( int proc=0; proc<this->worldComm().size(); ++proc )


    //--------------------------------------------------------------------------------------------------------//
    // get response to initial request and update Feel::Mesh::Faces data
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        {
            auto const myGlobProcessDof = mapMsg[proc][cpt];
            //------------------------------------------------------------------------------//
            // recv data
            boost::tuple<size_type,bool,int,size_type> eltAskedMultiProcess;
            this->worldComm().recv( proc, cpt, eltAskedMultiProcess );
            //------------------------------------------------------------------------------//
            // get data
            const int dofGlobRecv = eltAskedMultiProcess.template get<0>();
            const bool eltHasMultiProcess = eltAskedMultiProcess.template get<1>();
            const size_type eltMultiProcessProcId = eltAskedMultiProcess.template get<2>();
            const size_type eltMultiProcessFaceId = eltAskedMultiProcess.template get<3>();
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
            if ( eltHasMultiProcess && memoryFace[eltMultiProcessProcId].find( myGlobProcessDof ) == memoryFace[eltMultiProcessProcId].end() )
            {
                listNeedMoreSearch[eltMultiProcessProcId][eltMultiProcessFaceId].insert(boost::make_tuple(myGlobProcessDof,mapMsgComp[proc][cpt]));
                ++counterListNeedMoreSearchLoc;
            }
            //------------------------------------------------------------------------------//
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

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::buildDofNotPresent( std::map<size_type,boost::tuple<size_type,size_type> > const & mapInterProcessDof,
                                                                 std::set<int> & setInterProcessDofNotPresent )

{
    const int myRank = this->worldComm().rank();

    auto it_mapIPDof = mapInterProcessDof.begin();
    auto const en_mapIPDof = mapInterProcessDof.end();
    for ( ; it_mapIPDof!=en_mapIPDof; ++it_mapIPDof )
    {
        if ( it_mapIPDof->second.template get<0>() < myRank )
            setInterProcessDofNotPresent.insert(it_mapIPDof->first);
    }
}

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::buildGlobalProcessToGlobalClusterDofMap( mesh_type& mesh,
        std::map<size_type,boost::tuple<size_type,size_type> > const& mapInterProcessDof,
        std::set<int> const& setInterProcessDofNotPresent )
{
    const int myRank = this->worldComm().rank();

    //------------------------------------------------------------------------------//
    // update datamap info

    size_type mynDofWithoutGhost = this->_M_last_df[myRank] - this->_M_first_df[myRank] + 1 - setInterProcessDofNotPresent.size();
    mpi::all_gather( this->worldComm(),
                     mynDofWithoutGhost,
                     this->_M_n_localWithoutGhost_df );

    this->_M_n_dofs=0;
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        this->_M_n_dofs+=this->_M_n_localWithoutGhost_df[proc];
    }

    this->_M_first_df_globalcluster=this->_M_first_df;
    this->_M_last_df_globalcluster=this->_M_last_df;
    for ( int i=1; i<this->worldComm().size(); ++i )
    {
        this->_M_first_df_globalcluster[i]=this->_M_last_df_globalcluster[i-1]+1;
        this->_M_last_df_globalcluster[i]=this->_M_first_df_globalcluster[i]+this->_M_n_localWithoutGhost_df[i]-1;
    }


    //------------------------------------------------------------------------------//
    // init map
    this->M_mapGlobalProcessToGlobalCluster.resize( this->_M_n_localWithGhost_df[myRank] );
    this->M_mapGlobalClusterToGlobalProcess.resize( this->_M_n_localWithoutGhost_df[myRank] );
    std::fill( this->M_mapGlobalProcessToGlobalCluster.begin(),this->M_mapGlobalProcessToGlobalCluster.end(),invalid_size_type_value );
    std::fill( this->M_mapGlobalClusterToGlobalProcess.begin(),this->M_mapGlobalClusterToGlobalProcess.end(),invalid_size_type_value );

    //------------------------------------------------------------------------------//
    // add in map the dofs presents
    size_type firstGlobIndex = this->_M_first_df_globalcluster[myRank];
    size_type nextGlobIndex = firstGlobIndex;
    for ( int i=0; i< ( int )this->_M_n_localWithGhost_df[myRank]; ++i )
    {
        if ( setInterProcessDofNotPresent.find( i )==setInterProcessDofNotPresent.end() )
        {
            this->M_mapGlobalProcessToGlobalCluster[i]=nextGlobIndex;
            this->M_mapGlobalClusterToGlobalProcess[nextGlobIndex-firstGlobIndex]=i;
            ++nextGlobIndex;
        }
    }

    //------------------------------------------------------------------------------//
    // add in map the dofs non presents (ghosts)
    // this is important to build properly the interprocess dofs
    for ( int p = 1 ; p<this->comm().size() ; ++p )
    {
        this->updateGhostGlobalDof( mapInterProcessDof,
                                    setInterProcessDofNotPresent,
                                    p );
    }

} // buildGlobalProcessToGlobalClusterDofMap

//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------//

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::updateGhostGlobalDof( std::map<size_type,boost::tuple<size_type,size_type> > const& mapInterProcessDof,
        std::set<int> const& setInterProcessDofNotPresent,
        int procToUpdate )
{
    std::vector<int> nbMsgToSend2( this->worldComm().size() );
    std::fill( nbMsgToSend2.begin(),nbMsgToSend2.end(),0 );

    std::vector< std::map<int,int> > mapMsg2( this->worldComm().size() );

    // send msg only on proc
    if ( procToUpdate == this->comm().rank() )
    {
        auto dofNP_it = setInterProcessDofNotPresent.begin();
        auto dofNP_en = setInterProcessDofNotPresent.end();

        for ( ; dofNP_it!=dofNP_en ; ++dofNP_it )
        {
            int IdProcessOfGhost = mapInterProcessDof.find( *dofNP_it )->second.template get<0>();
            int dofFaceInPartition = mapInterProcessDof.find( *dofNP_it )->second.template get<1>();

            this->worldComm().send( IdProcessOfGhost , nbMsgToSend2[IdProcessOfGhost], dofFaceInPartition );

            mapMsg2[IdProcessOfGhost].insert( std::make_pair( nbMsgToSend2[IdProcessOfGhost],*dofNP_it ) );

            ++nbMsgToSend2[IdProcessOfGhost];
        }
    }

    // counter of msg received for each process
    std::vector<int> nbMsgToRecv2;
    mpi::all_to_all( this->worldComm(),
                     nbMsgToSend2,
                     nbMsgToRecv2 );

    //NEWBARRIER this->worldComm().barrier();

    // recv id asked and re-send set of face id
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToRecv2[proc]; ++cpt )
        {
            int dofRecv;
            //recv
            this->worldComm().recv( proc, cpt, dofRecv );

            // response // attention pas tout le temps exact : A mediter vincent (vrai pour 2 proc sur!)
            // je pense qu'il faut faire une boucle tant que pas tout est OK
            int dofToSend= this->M_mapGlobalProcessToGlobalCluster[dofRecv];

            if ( dofToSend==-1 ) std::cout << "\n AIEIAIIE rank "<< this->worldComm().rank() << "recv du proc " << proc <<  std::endl;

            if ( proc != procToUpdate ) std::cout << "\n[updateGhostGlobalDof] IL Y A UN SOUCIS !!!";

            this->worldComm().send( proc, cpt,dofToSend );
        }
    }

    // get response to initial request and update Feel::Mesh::Faces data
    //size_type firstGlobIndex = this->_M_first_df_globalcluster[this->comm().rank()];
    for ( int proc=0; proc<this->worldComm().size(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend2[proc]; ++cpt )
        {
            int dofGlobClusterRecv;
            //recv
            this->worldComm().recv( proc, cpt, dofGlobClusterRecv );
            //update data

            this->M_mapGlobalProcessToGlobalCluster[mapMsg2[proc][cpt]]=dofGlobClusterRecv;
            //this->M_mapGlobalClusterToGlobalProcess[dofGlobClusterRecv-firstGlobIndex]=mapMsg2[proc][cpt];
        }
    }

}

} // namespace Feel

#endif /* FEELPP_DOFTABLE_MPI_HPP */
