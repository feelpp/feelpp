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
    Feel::cout << "nx x ny = " << nx << " x " << ny << "\t" << nx*ny << std::endl;
    // origin at (0,0)
    node_type coords( 2 );
    rank_type partId = wc.localRank(); 
    std::vector<int> nbMsgToRecv(wc.godSize(),0);
    //std::vector<int> nbMsgToRecv((M_nx-1)*(M_ny-1),0);
    int procSize = wc.godSize();
    // cx[rank] = first point Id of partition
    int cx[procSize+1];
    rank_type id; 
    std::vector<int>::iterator lowProc;

    GmshOrdering<element_type> ordering;
   
    for (int tmpProc=0;tmpProc<=procSize;tmpProc++)
    {
        cx[tmpProc]=(tmpProc)*(M_ny-1)/procSize;
        Feel::cout << "cx[" << tmpProc << "] = " << cx[tmpProc] << "\t";
    }
    Feel::cout << "\n";
    
    /*
     * Nodes creations
     * NO ghost points to avoid costly tests at this moment
     */
#pragma omp parallel
#if 1
    // First column
   if(partId > 0)
   {
       int j = cx[partId];
       std::cout << "j = " << j << std::endl;
       ghosts.clear();
       ghosts.push_back(partId-1);
       for( int i = 0; i < M_nx; ++i )
       {
            // point
            int ptid = (M_ny)*i+j; 
            //std::cout << "ptId = " << ptid << std::endl;
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid, 
                           coords, 
                           i == 0 || i == M_nx-1 || j%((M_ny-1)) == 0  ); // Is on boundary ?
                           //i == 0 || i == M_nx-1 || j == 0 || j == M_ny-1 ); // Is on boundary ?
            // Actual rank ID 
            pt.setProcessId( partId-1);
            // NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId-1 );
            pt.setNeighborPartitionIds( ghosts );
            this->addPoint( pt );
        }
    }
#endif
    int start = (partId == 0) ?  0  : cx[partId]+1; 
    int stop = (partId == wc.size()-1) ?  M_ny : cx[partId+1]; 
    for( int j = start ; j < stop; ++j )
    {
        std::cout << "j = " << j << std::endl;
        ghosts.clear();
        if (j == cx[partId] && j != 0)
            ghosts.push_back(partId-1);
        else if (j == cx[partId+1]-1 && j != M_ny-1)
            ghosts.push_back(partId+1);

        for( int i = 0; i < M_nx; ++i )
        {
            // point
            int ptid = (M_ny)*i+j; 
            //std::cout << "ptId = " << ptid << std::endl;
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid, 
                           coords,
                           i == 0 || i == M_nx-1 || j%((M_ny-1)) == 0  ); // Is on boundary ?
                           //i == 0 || i == M_nx-1 || j == 0 || j == M_ny-1 ); // Is on boundary ?
            //// Actual rank ID 
            pt.setProcessId( partId );
            //// NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId );
            pt.setNeighborPartitionIds( ghosts );
            this->addPoint( pt );
        }
    }
#if 1
    // Last column
   if(partId < wc.size()-1)
   {
       int j = cx[partId+1];
       //std::cout << "j = " << j << std::endl;
       ghosts.clear();
       ghosts.push_back(partId+1);
       for( int i = 0; i < M_nx; ++i )
       {
            // point
            int ptid = (M_ny)*i+j; 
            //std::cout << "ptId = " << ptid << std::endl;
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid, 
                           coords,
                           i == 0 || i == M_nx-1 || j%((M_ny-1)) == 0  ); // Is on boundary ?
                           //i == 0 || i == M_nx-1 || j == 0 || j == M_ny-1 ); // Is on boundary ?
            // Actual rank ID 
            pt.setProcessId( partId );
            // NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId );
            pt.setNeighborPartitionIds( ghosts );
            this->addPoint( pt );
        }
    }
#endif
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
        else if (j == cx[partId+1]-1 && j != M_ny-1)
            ghosts.push_back(partId+1);
    
        for( int i = 0; i < M_nx-1; ++i )
        {
            element_type e;   
            int eid = (M_ny-1)*i+j; // GMSH Id
            //int eid =((M_ny-1)/procSize)*i+j%((M_ny-1)/procSize)+M_nx-1;
            e.setId( eid ); // eid == gmsh id
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );
            //int l_j = j-cx[procId];
            //int l_i = i;
            //int l_eid = (cx[procId+1]-cx[procId])*l_i+l_j;
            //std::cout << eid << "\t" << partId << "\t" << partId << "\n";
            e.setNeighborPartitionIds( ghosts );
               
            // points definition
            int ptid[4] = { (M_ny)*(i+1)+j,   // 0
                            (M_ny)*(i+1)+j+1, // 1
                            (M_ny)*i+j+1,     // 2
                            (M_ny)*i+j        // 3
                          };
           /* std::cout << "eid = " << eid << std::endl;
            for(auto it : ptid)
                std::cout << it << "\t" << std::endl;
            std::cout << std::endl;
      */
            
            for( int k = 0; k < 4; ++k )
            {
                //this->points().modify( this->pointIterator( ptid[k] ), Feel::detail::UpdateProcessId(e.processId()) );
                e.setPoint( ordering.fromGmshId(k), this->point( ptid[k]  ) );
            }
                //e.setPoint( k, this->point( ptid[k]  ) );
            //e.setOnBoundary( i==0 || j==0 || i==(M_nx-2) || j==(M_ny-2) );
            e.setOnBoundary( i==0 || i==(M_nx-2) || j%((M_ny-2))==0 );
                // true -> id = global elt id
                // false -> id = local elt id 
            this->addElement( e, true ); // e.id() is modified to Feel++ Id
            //std::cout << partId << "\t" << i << " : " << j << "\t" << e.id() << std::endl;
            __idGmshToFeel.insert( std::make_pair(eid,e.id())); // e.id() == Feel++ Id
            //__idGmshToFeel.insert( std::make_pair(eid,((M_ny-1)/procSize)*i+j%((M_ny-1)/procSize)+M_nx-1)); // e.id() == Feel++ Id
        }
    }
    toc("setElementInPart", FLAGS_v>0);
    tic();


#if 1
/* Ghost management */
    int j_inf =cx[partId]; 
    int j_sup =cx[partId+1]; 
    bool is_first = (partId == 0); /* proc 0 ? */
    bool is_last = (partId == procSize-1);
    if(! is_last ) // look at the right
    {
        /*
         * Indicate I am a ghost to neighborhood partition
         */
        ghosts.clear();
        ghosts.push_back(partId+1);  
        for( int i = 0; i < M_nx; ++i )
        {
            int ptId = (M_ny)*i+j_sup+1; 
            coords[0]=M_pixelsize*j_sup+1;
            coords[1]=M_pixelsize*i;
            point_type pt( ptId, coords,
                       i == 0 || i == M_nx-1  ); // Is on boundary ?

            pt.setProcessIdInPartition( partId );
            pt.setProcessId( partId+1 );  
            //pt.setNeighborPartitionIds( ghosts );

            this->addPoint( pt );
        } 
        for( int i = 0; i < M_nx-1; ++i )
        {    
            element_type e;   
            int eid = (M_ny-1)*i+j_sup;
            int eidF = i;
            //int eidF =(cx[partId]-cx[partId-1])*(M_nx-1)+i+M_nx-1 ;
            e.setId( eidF );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId+1 );
            //e.setIdInOtherPartitions(partId+1,i+M_nx-1);
            //e.setNeighborPartitionIds( ghosts );

            int ptid[4] = { (M_ny)*(i+1)+j_sup,   // 0
                            (M_ny)*(i+1)+j_sup+1, // 1
                            (M_ny)*i+j_sup+1,     // 2
                            (M_ny)*i+j_sup        // 3
                          };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
                
            e.setOnBoundary( i==0 || i==(M_nx-2) );
            nbMsgToRecv[partId+1]+=1; 
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
            this->addElement( e, false );
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
            mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( eidF, partId+1 )));
            /*mapGhostElt.insert( std::make_pair( eid,boost::make_tuple(
            (nx-1)*cx[partId+1]+i*(cx[partId+1]-cx[partId]), 
            partId+1 ) ) );*/
            __idGmshToFeel.insert( std::make_pair(eid,eidF)); 
            //__idGmshToFeel.insert( std::make_pair(eid,(nx-1)*cx[partId+1]+i*(cx[partId+1]-cx[partId]))); 


        }
    } // is_last
    ghosts.clear();
    if(! is_first) // look at the left
    {
        ghosts.push_back(partId-1);  
        for( int i = 0; i < M_nx; ++i )
        {
            int ptId = (M_ny)*i+j_inf-1; 
            coords[0]=M_pixelsize*j_inf-1;
            coords[1]=M_pixelsize*i;
            point_type pt( ptId, coords,
                    i == 0 || i == M_nx-1  ); // Is on boundary ?
            pt.setProcessIdInPartition( partId );
            pt.setProcessId( partId-1 );  
            //pt.setNeighborPartitionIds( ghosts );

            this->addPoint( pt );
        }

        for( int i = 0; i < M_nx-1; ++i )
        {
            element_type e;   
            int eid = (M_ny-1)*i+j_inf-1;
            int eidF =(cx[partId]-cx[partId-1]-1)*(M_nx-1)+i ;
            //int eidF = i ;
            e.setId( eidF );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId-1 );
            //e.setIdInOtherPartitions(partId-1,(cx[partId]-cx[partId-1])*(M_nx-1)+i);
            //e.setNeighborPartitionIds( ghosts );
            

            int ptid[4] = { (M_ny)*(i+1)+j_inf-1,   // 0
                            (M_ny)*(i+1)+j_inf, // 1
                            (M_ny)*i+j_inf,     // 2
                            (M_ny)*i+j_inf-1        // 3
            };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
            e.setOnBoundary( i==0 || i==(M_nx-2) );
            nbMsgToRecv[partId-1]+=1; 
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
            this->addElement( e, false ); // DO NOT modify e.id !!!
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
            mapGhostElt.insert( std::make_pair( eid,boost::make_tuple(eidF,partId-1)));

            /*mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( 
            (nx-1)*cx[partId-1]+(i+1)*(cx[partId]-cx[partId-1])-1, 
            partId-1 ) ) );*/
            __idGmshToFeel.insert( std::make_pair(eid,eidF)); 
            //__idGmshToFeel.insert( std::make_pair(eid,(nx-1)*cx[partId-1]+(i+1)*(cx[partId]-cx[partId-1])-1)); 

        }
    } // is_first
    
    for (std::map<int,boost::tuple<int,rank_type>>::iterator it=mapGhostElt.begin(); it!=mapGhostElt.end(); ++it)
        std::cout << "test : " << it->first << " ( " << it->second.get<0>()<< "," << it->second.get<1>() << " )" << std::endl;
    this->worldComm().barrier();
    updateGhostCellInfoByUsingNonBlockingComm(this,__idGmshToFeel,mapGhostElt,nbMsgToRecv);
    toc("Ghost Management", FLAGS_v>0);
    toc("MeshStructured Elements", FLAGS_v>0);
#endif
}



//template<typename MeshType>
    void
MeshStructured::updateGhostCellInfoByUsingNonBlockingComm( MeshStructured* mesh, 
        std::map<int,int> const& __idGmshToFeel, 
        std::map<int,boost::tuple<int,rank_type> > const& __mapGhostElt,
        std::vector<int> const& nbMsgToRecv )
{
    for(auto toto : nbMsgToRecv)
        std::cout << this->worldComm().localRank() << " : " << toto << "\n";
        

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
    LOG(INFO) << "TOTO\n";
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
    LOG(INFO) << "TOTO\n";
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
    LOG(INFO) << "TOTO\n";
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

    LOG(INFO) << "TOTO\n";
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
    LOG(INFO) << "TOTO\n";
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
    LOG(INFO) << "TOTO\n";
    //-----------------------------------------------------------//
    // wait all requests
    std::cout << "My, processor " << this->worldComm().localRank() << ", is waiting for " << nbRequest << " to be finished\n";
    mpi::wait_all(reqs, reqs + nbRequest);
    LOG(INFO) << "TOTO\n";
    //-----------------------------------------------------------//
    // build the container to ReSend
    std::map<rank_type, std::vector<int> > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        LOG(INFO) << "TOTO\n";
        const rank_type idProc = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        LOG(INFO) << "Resize to : " << nDataRecv << "\n";
        dataToReSend[idProc].resize( nDataRecv );
        //store the idFeel corresponding
        for ( int k=0; k<nDataRecv; ++k )
        {
            LOG(INFO) << "k = " << k << "\n";
            dataToReSend[idProc][k] = __idGmshToFeel.find( itDataRecv->second[k] )->second;
        }
    }
    LOG(INFO) << "TOTO\n";
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
    LOG(INFO) << "TOTO\n";
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
    LOG(INFO) << "TOTO\n";
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;
    //-----------------------------------------------------------//
    // update mesh : id in other partitions for the ghost cells
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    //for(auto toto = mesh->beginElement(); toto != mesh->endElement(); toto++)
    //    std::cout<< toto->processId()<< " : "   << toto->id() << std::endl;

    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const rank_type idProc = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        //std::cout << idProc << ":" << nDataRecv << std::endl;
        for ( int k=0; k<nDataRecv; ++k )
        {
           /* std::cout << "I want element " << memoryMsgToSend[idProc][k] << ": " << idProc << std::endl;*/
            auto eltToUpdate = mesh->elementIterator( memoryMsgToSend[idProc][k],idProc );
            std::cout << "k = " << k << std::endl;
            std::cout << "itFinalDataToRecv->second[k]  " << itFinalDataToRecv->second[k]  << std::endl;
            std::cout << "eltToUpdate->id()             " << eltToUpdate->id()             << std::endl;
            std::cout << "eltToUpdate->processId()      " << eltToUpdate->processId()      << std::endl;
            std::cout << "eltToUpdate->pidInPartition() " << eltToUpdate->pidInPartition() << std::endl;
            std::cout << "eltToUpdate->refDim()         " << eltToUpdate->refDim()         << std::endl;
            std::cout << "eltToUpdate->nPoints()        " << eltToUpdate->nPoints()        << std::endl;
            //auto toto = eltToUpdate->idInOthersPartitions();
            //std::cout << toto.size() << std::endl;
            //for (auto it : toto )
            //    std::cout << it.first << "\t" << it.second << std::endl;
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

