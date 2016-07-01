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

  private:
    int M_nx;
    int M_ny;
    double M_pixelsize;
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
    int procSize = this->worldComm().godSize();
    int cx[procSize+1];
    rank_type id; 
    std::vector<int>::iterator lowProc;

    for (int tmpProc=0;tmpProc<=procSize;tmpProc++)
        cx[tmpProc]=(tmpProc)*M_ny/procSize;
    //for (int tmpProc=0;tmpProc<procSize;tmpProc++)
       // cx[tmpProc]=(tmpProc+1)*M_ny/procSize;


    std::vector<int> vx(cx,cx+(procSize+1));
    
 

    for( int i = 0; i <= M_nx; ++i )
        for( int j = cx[partId]; j <= cx[partId+1]; ++j )
    //for( int i = 0; i <= M_nx; ++i )
        //for( int j = 0; j <= M_ny; ++j )
        {
            // point
            int ptid = (M_ny+1)*i+j; 
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid, coords, false );
            //lowProc=std::lower_bound(vx.begin(),vx.end(),j);
            //id = lowProc - vx.begin();
            //pt.setProcessIdInPartition( id );
            //pt.setProcessId( id );
            pt.setProcessIdInPartition( partId );
            pt.setProcessId( partId );  

            /*
            if ( (ptid%(M_ny/procSize)==0) && (ptid%M_ny !=0))
            {
                //ghosts.push_back(partId-1);
                //e.setNeighborPartitionIds(ghosts);
                //ghosts.clear();
                pt.addElementGhost(partId,ptid-1);

            }
            else if ( ((ptid+1)%(M_ny/procSize)==0) && ((ptid+1)%M_ny !=0))
            {
                //ghosts.push_back(partId+1);
                //e.setNeighborPartitionIds(ghosts);
                //ghosts.clear();
                pt.addElementGhost(partId,ptid+1);

            }
            */

            this->addPoint( pt );
        }
    toc("MeshStructured Nodes", FLAGS_v>0);
    tic();
    int eltid = 0;
    element_type e;
    std::vector<rank_type> ghosts;
   
   tic(); 
    for( int i = 0; i < M_nx; ++i )
        for( int j = cx[partId]; j < cx[partId+1]; ++j )
        //for( int j = 0; j < M_ny; ++j )
        {
            int eid = (M_ny)*i+j;
            e.setId( eid );

            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );

/*
       if ( partId < procSize )
            ghosts.push_back(partId+1);
        if ( partId > 0 )
            ghosts.push_back(partId-1);
            //lowProc=std::lower_bound(vx.begin(),vx.end(),j);
            //id = lowProc - vx.begin();
            //e.setProcessIdInPartition( id );
            //e.setProcessId( id );

         
            
            if ( (eid%(M_ny/procSize)==0) && (eid%M_ny !=0))
            {
                if (j!=0 && j!=M_ny)
                {
                    //ghosts.push_back((M_ny)*(i-1)+j-1);
                    //ghosts.push_back(eid-1);
                    //ghosts.push_back((M_ny)*(i+1)+j-1);
                    e.addElementGhost(partId,(M_ny)*(i-1)+j-1);
                    e.addElementGhost(partId,eid-1);
                    e.addElementGhost(partId,(M_ny)*(i+1)+j-1);
                }
                else 
                {
                    e.addElementGhost(partId,eid-1);
                    //ghosts.push_back(eid-1);
                }
                e.setNeighborPartitionIds(ghosts);
                ghosts.clear();
                //e.addElementGhost(partId,eid-1);

            }
            else if ( ((eid+1)%(M_ny/procSize)==0) && ((eid+1)%M_ny !=0))
            {
                ghosts.push_back(partId+1);
                if (j!=0 && j!=M_ny)
                {
                    //ghosts.push_back((M_ny)*(i-1)+j+1);
                    //ghosts.push_back(eid+11);
                    //ghosts.push_back((M_ny)*(i+1)+j+1);
                    e.addElementGhost(partId,(M_ny)*(i-1)+j+1);
                    e.addElementGhost(partId,eid+1);
                    e.addElementGhost(partId,(M_ny)*(i+1)+j+1);

                }
                else 
                {
                     e.addElementGhost(partId,eid+1);
                    //ghosts.push_back(eid+1);
                }

                e.setNeighborPartitionIds(ghosts);
                ghosts.clear();
                //e.addElementGhost(partId,eid+1);

            }*/
           /* else 
            {
                ghosts.push_back(partId);
                e.setNeighborPartitionIds(ghosts);
                ghosts.clear();
            }
*/
            //
            int ptid[4] = { (M_ny+1)*(i+1)+j,   // 0
                            (M_ny+1)*(i+1)+j+1, // 1
                            (M_ny+1)*i+j+1,     // 2
                            (M_ny+1)*i+j        // 3
            };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
                // true -> id = global elt id
                // false -> id = local elt id
            this->addElement( e, true );
        }
        toc("setElementInPart", FLAGS_v>0);
        tic();

/* Ghost management */
            int j_inf =cx[partId]; 
            int j_sup =cx[partId+1]-1; 
            bool is_first = (partId == 0); /* proc 0 ? */
            bool is_last = (partId == procSize-1);
            int loc_m_ny = cx[partId+1]-1-cx[partId];
#pragma omp parallel
{
            if(! is_last ) // look at the right
            {
            /*
             * Get my Id : 
             * Get my right neighbour ID
             * Set it as ghost
             */ 
    for( int i = 1; i < M_nx-1; ++i )
        {
            auto el = const_cast<element_type&>(this->element(globalToLocal(i,j_sup,partId)));
            el.addElementGhost(partId,(M_ny)*(i+1)+j_sup+1);
            el.addElementGhost(partId,(M_ny)*(i+0)+j_sup+1);
            el.addElementGhost(partId,(M_ny)*(i-1)+j_sup+1);
        }
        {
            int i = 0;
            auto el = const_cast<element_type&>(this->element(globalToLocal(i,j_sup,partId)));
            el.addElementGhost(partId,(M_ny)*(i+1)+j_sup+1);
            el.addElementGhost(partId,(M_ny)*(i+0)+j_sup+1);
            //el.addElementGhost(partId,(M_ny)*(i-1)+j_sup+1);
        }
        {
            int i = M_nx-1;
            auto el = const_cast<element_type&>(this->element(globalToLocal(i,j_sup,partId)));
            //el.addElementGhost(partId,(M_ny)*(i+1)+j_sup+1);
            el.addElementGhost(partId,(M_ny)*(i+0)+j_sup+1);
            el.addElementGhost(partId,(M_ny)*(i-1)+j_sup+1);
        }
            }
            if(! is_first) // look at the left
            {
    for( int i = 1; i < M_nx-1; ++i )
        {
            auto el = const_cast<element_type&>(this->element(globalToLocal(i,j_inf,partId)));
            el.addElementGhost(partId,(M_ny)*(i+1)+j_inf-1);
            el.addElementGhost(partId,(M_ny)*(i+0)+j_inf-1);
            el.addElementGhost(partId,(M_ny)*(i-1)+j_inf-1);
        }
        {
            int i = 0;
            auto el = const_cast<element_type&>(this->element(globalToLocal(i,j_inf,partId)));
            el.addElementGhost(partId,(M_ny)*(i+1)+j_inf-1);
            el.addElementGhost(partId,(M_ny)*(i+0)+j_inf-1);
            //el.addElementGhost(partId,(M_ny)*(i-1)+j_inf-1);
        }
        {
            int i = M_nx-1;
            auto el = const_cast<element_type&>(this->element(globalToLocal(i,j_inf,partId)));
            //el.addElementGhost(partId,(M_ny)*(i+1)+j_inf-1);
            el.addElementGhost(partId,(M_ny)*(i+0)+j_inf-1);
            el.addElementGhost(partId,(M_ny)*(i-1)+j_inf-1);
        }
            }
        }
    toc("Ghost Management", FLAGS_v>0);
    toc("MeshStructured Elements", FLAGS_v>0);
}
}

#endif
