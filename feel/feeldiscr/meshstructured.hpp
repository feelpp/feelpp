/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Jun 2016

 Copyright (C) 2016 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_MESHSTRUCTURED_HPP
#define FEELPP_MESHSTRUCTURED_HPP 1

namespace Feel {

/**
 *
 */
class MeshStructured: public Mesh<Hypercube<2>>
{

  public:
    using super = Mesh<Hypercube<2>>;
    using point_type = super::point_type;
    using element_type = super::element_type;
    using node_type = super::node_type;
    MeshStructured() = default;
    MeshStructured( MeshStructured const& ) = default;
    MeshStructured( MeshStructured && ) = default;
    MeshStructured& operator=( MeshStructured const& ) = default;
    MeshStructured& operator=( MeshStructured && ) = default;
                 
    //!
    //! 
    //! 
    MeshStructured( int nx, int ny, double pixelsize, WorldComm const& );

    void updateGhostCellInfoByUsingNonBlockingComm( 
        MeshStructured* mesh, 
        std::map<int,int> const& __idGmshToFeel, 
        std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
        std::vector<int> const& nbMsgToRecv );
    
   void updateGhostCellInfoByUsingBlockingComm( 
        MeshStructured* mesh, 
        std::map<int,int> const& __idGmshToFeel, 
        std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
        std::vector<int> const& nbMsgToRecv );

    /*
     * Create an element
     * g_i = global x index 
     * g_j = global y index 
     */
   element_type newElement(int g_i, int g_j)
   {
       element_type e;
       return e;
   }


  private:
    int M_nx; // Global X number of elements
    int M_ny; // Global Y number of elements
    int M_l_nx; // local X number of elements (ghost excluded!)
    int M_l_ny; // local Y number of elements
    int M_s_x; // local first x index (0 for first element)
    int M_s_y; // local first y index (0 for first element)

    double M_pixelsize;
    std::map<int,boost::tuple<int,rank_type> > mapGhostElt; 
    std::vector<rank_type> ghosts; 
    std::map<int,int> __idGmshToFeel;


    int localToGlobal(int ii, int jj , int rank)
    {
       return ii*M_ny+(jj+ (rank*M_ny/this->worldComm().godSize()));
       
    }
    int globalToLocal(int i, int j, int rank)
    {
        int LM_ny= (M_ny/this->worldComm().godSize())*(rank+1)-(M_ny/this->worldComm().godSize())*rank;
        return i*LM_ny+j-rank*M_ny/this->worldComm().godSize();
    }

};

MeshStructured::MeshStructured( int nx, int ny, double pixelsize, WorldComm const& wc )
    :
    super( wc ),
    M_nx( nx ),
    M_ny( ny ),
    M_pixelsize( pixelsize )
{
    tic();
    // origin at (0,0)
    node_type coords( 2 );
    rank_type partId = this->worldComm().localRank(); 
    std::vector<int> nbMsgToRecv(M_nx*M_ny,0);
    int procSize = this->worldComm().godSize();
    int cx[procSize+1];
    rank_type id; 
    std::vector<int>::iterator lowProc;
    
    for (int tmpProc=0;tmpProc<=procSize;tmpProc++)
        cx[tmpProc]=(tmpProc)*M_ny/procSize;



    std::vector<int> vx(cx,cx+(procSize+1));
    
    
    /*
     * Nodes creations
     * NO ghost points to avoid costly tests at this moment
     */
#pragma omp parallel
    for( int j = cx[partId]; j <= cx[partId+1]; ++j )
    {
        ghosts.clear();
        if (j == cx[partId] && j != 0)
            ghosts.push_back(partId-1);
        else if (j == cx[partId+1] && j != M_ny)
            ghosts.push_back(partId+1);

        for( int i = 0; i <= M_nx; ++i )
        {
            // point
            int ptid = (M_ny+1)*i+j; 
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid, 
                           coords, 
                           i == 0 || i == M_nx || j == 0 || j == M_ny ); // Is on boundary ?
            // Actual rank ID 
            pt.setProcessId( partId );
            // NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId );
            pt.setNeighborPartitionIds( ghosts );
#if 0
            /// Indicates which elements the current point belongs to
            int id_0 = -1, id_1=-1, id_2=-1, id_3=-1;
            if( i == 0 )
             {
                 only id_3 & id_0
             }
             else if( i == M_nx )
             {
                 only id_1 & id_2
             }
            (id_0 >= 0) ? pt.addElement( id_0, 0) : false; 
            (id_1 >= 0) ? pt.addElement( id_1, 1) : false; 
            (id_2 >= 0) ? pt.addElement( id_2, 2) : false; 
            (id_3 >= 0) ? pt.addElement( id_3, 3) : false; 
#endif
            this->addPoint( pt );
        }
    }
    toc("MeshStructured Nodes", FLAGS_v>0);
    tic();
    int eltid = 0;
    /*
     * Generates elements whithout ghosts
     */
    tic(); 
    for( int j = cx[partId]; j < cx[partId+1]; ++j ) 
    {
        ghosts.clear();
        if (j == cx[partId] && j != 0)
            ghosts.push_back(partId-1);
        else if (j == cx[partId]-1 && j != M_ny-1)
            ghosts.push_back(partId+1);
    
        for( int i = 0; i < M_nx; ++i )
        {
            element_type e;   
            int eid = (M_ny)*i+j; // GMSH Id
            e.setId( eid ); // eid == gmsh id
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );
            //int l_j = j-cx[procId];
            //int l_i = i;
            //int l_eid = (cx[procId+1]-cx[procId])*l_i+l_j;
            //std::cout << eid << "\t" << partId << "\t" << partId << "\n";
            e.setNeighborPartitionIds( ghosts );
               
            // points definition
            int ptid[4] = { (M_ny+1)*(i+1)+j,   // 0
                            (M_ny+1)*(i+1)+j+1, // 1
                            (M_ny+1)*i+j+1,     // 2
                            (M_ny+1)*i+j        // 3
                          };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
            e.setOnBoundary( i==0 || j==0 || i==(M_nx-1) || j==(M_ny-1) );
                // true -> id = global elt id
                // false -> id = local elt id 
            this->addElement( e, true ); // e.id() is modified to Feel++ Id
            __idGmshToFeel.insert( std::make_pair(eid,e.id())); // e.id() == Feel++ Id
        }
    }
    toc("setElementInPart", FLAGS_v>0);
    tic();

/* Ghost management */
    int j_inf =cx[partId]; 
    int j_sup =cx[partId+1]; 
    bool is_first = (partId == 0); /* proc 0 ? */
    bool is_last = (partId == procSize-1);
    int loc_m_ny = cx[partId+1]-1-cx[partId];
    if(! is_last ) // look at the right
    {
        /*
         * Indicate I am a ghost to neighborhood partition
         */
        ghosts.clear();
        ghosts.push_back(partId+1);  
        for( int i = 0; i <= M_nx; ++i )
        {
            int ptId = (M_ny+1)*i+j_sup+1; 
            coords[0]=M_pixelsize*j_sup+1;
            coords[1]=M_pixelsize*i;
            point_type pt( ptId, coords,
                       i == 0 || i == M_nx  ); // Is on boundary ?

            pt.setProcessIdInPartition( partId );
            pt.setProcessId( partId+1 );  
            //pt.setNeighborPartitionIds( ghosts );

            this->addPoint( pt );
        } 
        for( int i = 0; i < M_nx; ++i )
        {    
            element_type e;   
            int eid = (M_ny)*i+j_sup;
            e.setId( eid );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId+1 );
            //e.setNeighborPartitionIds( ghosts );

            int ptid[4] = { (M_ny+1)*(i+1)+j_sup,   // 0
                            (M_ny+1)*(i+1)+j_sup+1, // 1
                            (M_ny+1)*i+j_sup+1,     // 2
                            (M_ny+1)*i+j_sup        // 3
                          };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
                
            e.setOnBoundary( i==0 || i==(M_nx-1) );
            nbMsgToRecv[eid-1]+=1; 
            //std::cout << partId << std::endl;
            //__idGmshToFeel.insert( std::make_pair(eid,eid)); 
            // true -> id = global elt id
            // false -> id = local elt id 
            /*std::cout << "*******************************\n";
            std::cout << "id \t"    << e.id() << std::endl;
            std::cout << "pid\t"    << e.processId() << std::endl;
            std::cout << "pidinp\t" << e.pidInPartition() << std::endl;
            std::cout << "refd\t"   << e.refDim() << std::endl;
            std::cout << "np\t"     << e.nPoints() << std::endl;
            auto tete =  e.idInOthersPartitions();
            std::cout << tete.size() << std::endl;
            for (auto it : tete )
                std::cout << it.first << "\t" << it.second << std::endl;*/
            auto toto = this->addElement( e, true );
            /*std::cout <<"id \t"    <<  toto.id() << std::endl;
            std::cout <<"pid\t"    <<  toto.processId() << std::endl;
            std::cout <<"pidinp\t" <<  toto.pidInPartition() << std::endl;
            std::cout <<"refd\t"   <<  toto.refDim() << std::endl;
            std::cout <<"np\t"     <<  toto.nPoints() << std::endl;
            auto titi =  toto.idInOthersPartitions();
            std::cout << titi.size() << std::endl;
            for (auto it : titi )
                std::cout << it.first << "\t" << it.second << std::endl;
            */
                
            mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( e.id(), partId+1 ) ) );


        }
    } // is_last
    ghosts.clear();
    if(! is_first) // look at the left
    {
        ghosts.push_back(partId-1);  
        for( int i = 0; i <= M_nx; ++i )
        {
            int ptId = (M_ny+1)*i+j_inf-1; 
            coords[0]=M_pixelsize*j_inf-1;
            coords[1]=M_pixelsize*i;
            point_type pt( ptId, coords,
                    i == 0 || i == M_nx  ); // Is on boundary ?
            pt.setProcessIdInPartition( partId );
            pt.setProcessId( partId-1 );  
            //pt.setNeighborPartitionIds( ghosts );

            this->addPoint( pt );
        }

        for( int i = 0; i < M_nx; ++i )
        {
            element_type e;   
            int eid = (M_ny)*i+j_inf-1;
            e.setId( eid );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId-1 );
            //e.setNeighborPartitionIds( ghosts );

            int ptid[4] = { (M_ny+1)*(i+1)+j_inf-1,   // 0
                (M_ny+1)*(i+1)+j_inf, // 1
                (M_ny+1)*i+j_inf,     // 2
                (M_ny+1)*i+j_inf-1        // 3
            };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
            e.setOnBoundary( i==0 || i==(M_nx-1) );
            nbMsgToRecv[eid+1]+=1; 
            // true -> id = global elt id
            // false -> id = local elt id
            /*std::cout << "*******************************\n";
            std::cout << "id \t"    << e.id() << std::endl;
            std::cout << "pid\t"    << e.processId() << std::endl;
            std::cout << "pidinp\t" << e.pidInPartition() << std::endl;
            std::cout << "refd\t"   << e.refDim() << std::endl;
            std::cout << "np\t"     << e.nPoints() << std::endl;
            auto tete =  e.idInOthersPartitions();
            std::cout << tete.size() << std::endl;
            for (auto it : tete )
                std::cout << it.first << "\t" << it.second << std::endl;*/
            //this->addElement( e, true );
            auto toto = this->addElement( e, true );
            /*std::cout <<"id \t"    <<  toto.id() << std::endl;
            std::cout <<"pid\t"    <<  toto.processId() << std::endl;
            std::cout <<"pidinp\t" <<  toto.pidInPartition() << std::endl;
            std::cout <<"refd\t"   <<  toto.refDim() << std::endl;
            std::cout <<"np\t"     <<  toto.nPoints() << std::endl;
            auto titi =  toto.idInOthersPartitions();
            std::cout << titi.size() << std::endl;
            for (auto it : titi )
                std::cout << it.first << "\t" << it.second << std::endl;
            */    
                
            mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( e.id(), partId-1 ) ) ); 

        }
    } // is_first
    this->worldComm().barrier();
    updateGhostCellInfoByUsingNonBlockingComm(this,__idGmshToFeel,mapGhostElt,nbMsgToRecv);
    toc("Ghost Management", FLAGS_v>0);
    toc("MeshStructured Elements", FLAGS_v>0);

}



//template<typename MeshType>
    void
MeshStructured::updateGhostCellInfoByUsingNonBlockingComm( MeshStructured* mesh, 
        std::map<int,int> const& __idGmshToFeel, 
        std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
        std::vector<int> const& nbMsgToRecv )
{

    DVLOG(1) << "updateGhostCellInfoNonBlockingComm : start on rank "<< this->worldComm().localRank() << "\n";

    const int nProc = this->worldComm().localSize();

    //std::cout << nProc << std::endl;

    //-----------------------------------------------------------//
    // compute size of container to send
    std::map< rank_type, int > nDataInVecToSend;
    auto it_map = __mapGhostElt.begin();
    auto const en_map = __mapGhostElt.end();
    for ( ; it_map!=en_map ; ++it_map )
    {
        const rank_type idProc = it_map->second.template get<1>();
        if ( nDataInVecToSend.find(idProc) == nDataInVecToSend.end() )
            nDataInVecToSend[idProc]=0;
        nDataInVecToSend[idProc]++;
    }
    //-----------------------------------------------------------//
    // init and resize the container to send
    std::map< rank_type, std::vector<int> > dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend!=enNDataInVecToSend ; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //-----------------------------------------------------------//
    // prepare container to send
    std::map< rank_type, std::map<int,int> > memoryMsgToSend;
    std::map< rank_type, int > nDataInVecToSendBis;
    it_map = __mapGhostElt.begin();
    for ( ; it_map!=en_map ; ++it_map )
    {
        const int idGmsh = it_map->first;
        const int idFeel = it_map->second.template get<0>();
        const rank_type idProc = it_map->second.template get<1>();

        if ( nDataInVecToSendBis.find(idProc) == nDataInVecToSendBis.end() )
            nDataInVecToSendBis[idProc]=0;
        // save request
        memoryMsgToSend[idProc][nDataInVecToSendBis[idProc]] = idFeel;
        // update container
        dataToSend[idProc][nDataInVecToSendBis[idProc]] = idGmsh;
        // update counter
        nDataInVecToSendBis[idProc]++;
        // std::cout << idProc << std::endl;

    }
    //-----------------------------------------------------------//
    // counter of request
    int nbRequest=0;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( dataToSend.find(proc) != dataToSend.end() )
            ++nbRequest;
        if ( nbMsgToRecv[proc] > 0 )
            ++nbRequest;
    }
    if ( nbRequest ==0 ) return;

    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;
    //-----------------------------------------------------------//
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToSend->first , 0, itDataToSend->second );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // first recv
    std::map<rank_type,std::vector<int> > dataToRecv;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( nbMsgToRecv[proc] > 0 )
        {
            reqs[cptRequest] = this->worldComm().localComm().irecv( proc , 0, dataToRecv[proc] );
            ++cptRequest;
        }
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    //-----------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector<int> > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        dataToReSend[idProc].resize( nDataRecv );
        //store the idFeel corresponding
        for ( int k=0; k<nDataRecv; ++k )
            dataToReSend[idProc][k] = __idGmshToFeel.find( itDataRecv->second[k] )->second;
    }
    //-----------------------------------------------------------//
    // send respond to the request
    cptRequest=0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend!=enDataToReSend ; ++itDataToReSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToReSend->first , 0, itDataToReSend->second );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // recv the initial request
    std::map<rank_type, std::vector<int> > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc, 0, finalDataToRecv[idProc] );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //-----------------------------------------------------------//
    // update mesh : id in other partitions for the ghost cells
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const rank_type idProc = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        std::cout << idProc << ":" << nDataRecv << std::endl;
        for ( int k=0; k<nDataRecv; ++k )
        {
            auto eltToUpdate = mesh->elementIterator( memoryMsgToSend[idProc][k],idProc );
            /*std::cout << "k = " << k << std::endl;
            std::cout << itFinalDataToRecv->second[k] << std::endl;
            std::cout << eltToUpdate->id() << std::endl;
            std::cout << eltToUpdate->processId() << std::endl;
            std::cout << eltToUpdate->pidInPartition() << std::endl;
            std::cout << eltToUpdate->refDim() << std::endl;
            std::cout << eltToUpdate->nPoints() << std::endl;
            auto toto = eltToUpdate->idInOthersPartitions();
            std::cout << toto.size() << std::endl;
            for (auto it : toto )
                std::cout << it.first << "\t" << it.second << std::endl;*/
            mesh->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( idProc, itFinalDataToRecv->second[k] ) );
        }
    }
    //-----------------------------------------------------------//
    DVLOG(1) << "updateGhostCellInfoNonBlockingComm : finish on rank "<< this->worldComm().localRank() << "\n";
}


    void
MeshStructured::updateGhostCellInfoByUsingBlockingComm( MeshStructured* mesh, std::map<int,int> const& __idGmshToFeel, std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
        std::vector<int> const& nbMsgToRecv )
{
    // counter of msg sent for each process
    std::vector<int> nbMsgToSend( this->worldComm().localSize(),0 );
    // map usefull to get final result
    std::vector< std::map<int,int> > mapMsg( this->worldComm().localSize() );

    // iterate over ghost elt
    auto it_map = __mapGhostElt.begin();
    auto const en_map = __mapGhostElt.end();
    for ( ; it_map!=en_map ; ++it_map )
    {
        auto const idGmsh = it_map->first;
        auto const idProc = it_map->second.template get<1>();

        // send
        this->worldComm().localComm().send( idProc , nbMsgToSend[idProc], idGmsh );
        // save tag of request
        mapMsg[idProc].insert( std::make_pair( nbMsgToSend[idProc],it_map->second.template get<0>() ) );
        // update nb send
        nbMsgToSend[idProc]++;
    }

#if !defined( NDEBUG )
    // check nbMsgToRecv computation
    std::vector<int> nbMsgToRecv2;
    mpi::all_to_all( this->worldComm().localComm(),
            nbMsgToSend,
            nbMsgToRecv2 );
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        CHECK( nbMsgToRecv[proc]==nbMsgToRecv2[proc] ) << "paritioning data incorect "
            << "myrank " << this->worldComm().localRank() << " proc " << proc
            << " nbMsgToRecv[proc] " << nbMsgToRecv[proc]
            << " nbMsgToRecv2[proc] " << nbMsgToRecv2[proc] << "\n";
    }
#endif

    // get gmsh id asked and re-send the correspond id Feel
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0 ; cpt<nbMsgToRecv[proc] ; ++cpt )
        {
            int idGmsh;
            //reception idGmsh
            this->worldComm().localComm().recv( proc, cpt, idGmsh );
            this->worldComm().localComm().send( proc, cpt, __idGmshToFeel.find( idGmsh )->second );
        }
    }

    // get response to initial request and update Feel::Mesh data
    for ( int proc=0; proc<this->worldComm().localSize(); ++proc )
    {
        for ( int cpt=0; cpt<nbMsgToSend[proc]; ++cpt )
        {
            int idFeel;
            // receive idFeel
            this->worldComm().localComm().recv( proc, cpt, idFeel );
            // update data
            auto elttt = mesh->elementIterator( mapMsg[proc][cpt],proc );
            mesh->elements().modify( elttt, Feel::detail::updateIdInOthersPartitions( proc, idFeel ) );
#if 0
            std::cout << "[updateGhostCellInfo]----3---\n"
                << "END! I am the proc" << this->worldComm().localRank()<<" I receive of the proc " << proc
                <<" with tag "<< cpt
                << " for the G " << mesh->element( mapMsg[proc][cpt], proc ).G()
                << " with idFeel Classic " << mapMsg[proc][cpt]
                << " with idFeel " << idFeel
                << " and modif " << mesh->element( mapMsg[proc][cpt] , proc ).idInPartition( proc )
                << std::endl;
#endif
        }
    }

}




}



#endif

