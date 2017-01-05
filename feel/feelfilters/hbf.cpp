/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2014-06-04

Copyright (C) 2014 Feel++ Consortium

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
#include <fstream>
#include <iostream>
#include <feel/feelcore/environment.hpp>
//#include <fftw3.h>
#include "hbf.hpp"


namespace Feel
{
holo3_image<float>
readHBF( std::string const& s )
{
    #if 0
        std::ifstream in( Environment::expand(s).c_str(), std::ios::binary );
    #else
        std::ifstream in( s, std::ios::binary );
    #endif  
    if ( !in )
    {
        std::cout << "Error opening file " << s.c_str() << std::endl;
        exit(0);
    }

    std::string header;
    char ch;
    size_t count = 0;
    while ((ch = in.get()) != '\0')
    {
        header += ch;
        ++count;
    }
    if ( Environment::isMasterRank() )
    {
        std::cout << "header: " << header << "\n";
        std::cout << count << std::endl;
    }
    int32_t version = 0;
    in.read( (char*)&version, sizeof( int32_t ) );
    int32_t rows=0, cols=0;
    in.read( (char*)&rows, sizeof( int32_t ) );
    in.read( (char*)&cols, sizeof( int32_t ) );
    if ( Environment::isMasterRank() )
        std::cout << "rows: " << rows << " , cols: " << cols << " , size: " << rows*cols << std::endl;
    Eigen::MatrixXf x( cols, rows  );
    in.read( (char*)x.data(), x.size()*sizeof(float) );
    if(x.rows() <= 6 && x.cols() <= 6)
        std::cout << x << std::endl;
   // if ( Environment::isMasterRank() )
       // std::cout << "x.rows: " << x.rows() << " , x.cols(): " << x.cols() << std::endl;
    return x.transpose();
}
   
void
writeHBF( std::string const& s, holo3_image<float> const& x )
{
    #if 0  
        std::ofstream out( Environment::expand(s).c_str(), std::ios::binary );
    #else 
        std::ofstream out( s, std::ios::binary );
    #endif    
    if ( !out )
    {
        std::cout << "Error opening file " << s.c_str() << std::endl;
        exit(0);
    }
    const char* header = "class CH3Array2D<float> *";
    out.write( header, strlen(header) );
    out.write( "\0",sizeof(char) );

    if ( Environment::isMasterRank() )
    {
        std::cout << "writeHBF: write header " << header << "\n";
    }
    int32_t version = 100;
    out.write( (char*)&version, sizeof(int32_t) );
    int32_t rows=x.rows(), cols=x.cols();
    out.write( (char*)&rows, sizeof(int32_t) );
    out.write( (char*)&cols, sizeof(int32_t) );
    if ( Environment::isMasterRank() )
    {
        std::cout << "writeHBF: rows: " << rows << " , cols: " << cols << " , size: " << rows*cols << std::endl;
    }


    out.write( (char*)x.data(),x.size()*sizeof(float) );
    if ( Environment::isMasterRank() )
    {
        std::cout << "writeHBF: array written to disk" << std::endl;
    }
    if(x.rows() <= 6 && x.cols() <= 6)
        std::cout << x << std::endl;
    if ( Environment::isMasterRank() )
        std::cout << "x.rows: " << x.rows() << " , x.cols(): " << x.cols() << std::endl;
}



holo3_image<float> cutHbf ( holo3_image<float> const& im, int n, int start )
{
    holo3_image<float> res (n+1,n+1);
    for (int i=start;i<=start+n;i++)
        for (int j=start;j<=start+n;j++)
            res(i-start,j-start)=im(i,j);
    return res;
}
//
// Hbf2Feelpp
//

Hbf2Feelpp::Hbf2Feelpp( int nx, int ny, q1_space_ptrtype Yh ):M_rows(ny), M_cols(nx), M_Xh( Yh )
{
    tic();
    auto relation = Yh->dof()->pointIdToDofRelation();
    LOG(INFO) << "relation.size() = " << relation.first.size() << "\t" << relation.second.size()<< std::endl;
    for( int i = 0; i < ny; ++i )
    {
        for( int j = 0; j < nx; ++j )
        {
            int grid_vid = (nx+1)*i+j;
            int mesh_vid = -1;
            
            // vertices
            if ( ( i == 0    ) && ( j == 0    ) ) mesh_vid = 4;
            else if ( ( i == ny-1 ) && ( j == 0    ) ) mesh_vid = 1;
            else if ( ( i == ny-1 ) && ( j == nx-1 ) ) mesh_vid = 2;
            else if ( ( i == 0    ) && ( j == nx-1 ) ) mesh_vid = 3;

            // edges
            if ( ( j == 0    ) && ( ( i > 0 ) && ( i < ny-1 ) ) ) mesh_vid = 4 + i;
            else if ( ( i == ny-1 ) && ( ( j > 0 ) && ( j < nx-1 ) ) ) mesh_vid = 4 + ny-2 + j;
            else if ( ( j == nx-1 ) && ( ( i > 0 ) && ( i < ny-1 ) ) ) mesh_vid = 4 + ny-2 + nx-2 + (ny-1-i);
            else if ( ( i == 0    ) && ( ( j > 0 ) && ( j < nx-1 ) ) ) mesh_vid = 4 + ny-2 + nx-2 + ny-2 + (nx-1-j);

            // interior
            if ( ( i > 0 ) && ( i < ny-1 ) && ( j > 0 ) && ( j < nx -1 ) ) mesh_vid = 4+2*(nx-2)+2*(ny-2) + (ny-2)*(j-1) + (ny-1-i);

            if ( mesh_vid != -1 )
            {
                int dofid = relation.second[mesh_vid];
                M_relation.insert( dof_relation( std::make_pair(i,j), dofid ) );
            }
        }
    }
    toc("structured 2 feelpp relation");
}

Hbf2Feelpp::q1_element_type
Hbf2Feelpp::operator()( holo3_image<float> const& x )
{
    q1_element_type u( M_Xh );
    for( auto dof : M_relation.left )
    {
        u( dof.second ) = x(dof.first.first,dof.first.second);
    }
    return u;
}

//
//Hbf2FeelppStruc
//

/*
 * Nx = number of columns (= number of dots in P1 in X-holo3 direction)
 * Ny = number of lines
 */
Hbf2FeelppStruc::Hbf2FeelppStruc( int nx, int ny, q1_space_ptrtype Yh ):M_rows(ny), M_cols(nx), M_Xh( Yh )
{
    tic();
    std::cout << "Hbf2FeelppStruc relation creation \n";
    std::cout << nx << " : " << ny <<std::endl;
    
    int procSize = Environment::worldComm().godSize();
    rank_type partId = Environment::worldComm().localRank(); 
    int cx[procSize+1];
    for (int tmpProc=0;tmpProc<=procSize;tmpProc++)
        cx[tmpProc]=(tmpProc)*(nx-1)/procSize;

    auto relation = Yh->dof()->pointIdToDofRelation();
    //LOG(INFO) << "relation.size() = " << relation.first.size() << "\t" << relation.second.size()<< std::endl;
//#pragma omp parallel for
    for( int i = 0; i < ny; ++i )
    {
        for( int j = std::max(cx[partId]-1,0); j <= std::min(cx[partId+1]+1,nx-1); ++j )
        //for( int j = cx[partId]; j < cx[partId+1]; ++j )
        {
            //LOG(INFO) << i << " : " << j << "\t->\t" << (nx)*i+j << std::endl;
            //int grid_vid = (nx+1)*i+j;
            //int dofid = relation.second[grid_vid];
            //int dofid = (nx+1)*i+j;
            M_relation.insert( dof_relation( std::make_pair(i,j), relation.second[(nx)*i+j] ));
           //int dofid = (ny/procSize)*(j%((nx)/procSize))+i;
           //int dofid = (nx)*(i%(ny/procSize))+j;
           //std::cout << "r : " << relation.second[(nx)*i+j] <<  std::endl;
           //std::cout << "r : " << relation.second[(nx)*i+j] << " dof : "<< dofid << std::endl;
           //if (relation.second[nx*i+j] != dofid) 
           //     std::cout << "STOOOOOOOOOOOP ------------------------------------"<< "r : " << relation.second[(nx)*i+j] << " dof : "<< dofid  << std::endl;
           //M_relation.insert( dof_relation( std::make_pair(i,j), dofid ));
       }
    }
    //std::cout << "M_relation.size() = " << M_relation.size() << std::endl;
    toc("structured 2 feelpp relation");
}

Hbf2FeelppStruc::q1_element_type
Hbf2FeelppStruc::operator()( holo3_image<float> const& x )
{
    q1_element_type u( M_Xh );
    //std::cout << "u.nDof() = " << u.nDof() << std::endl;
    //std::cout << "M_relation.size() = " << M_relation.size() << std::endl;
//#pragma omp parallel for
    for( auto const & dof : M_relation.left )
    {
        //LOG(INFO) << "dof.second = " << dof.second << " = " << dof.first.first << " : " << dof.first.second << std::endl;
        //tic();
        u( dof.second ) = x(dof.first.first,dof.first.second);
        //toc("hbfTOFeel::operator()");
    }
    sync(u,"=");
    u.close();
    return u;

}

Hbf2FeelppStruc::q1_element_type
Hbf2FeelppStruc::operator()( holo3_image<float> const& x, q1_element_type u)
{
    //q1_element_type u( M_Xh );
    //std::cout << "u.nDof() = " << u.nDof() << std::endl;
    //std::cout << "M_relation.size() = " << M_relation.size() << std::endl;
//#pragma omp parallel for
    for( auto const & dof : M_relation.left )
    {
        //LOG(INFO) << "dof.second = " << dof.second << " = " << dof.first.first << " : " << dof.first.second << std::endl;
        //tic();
        u( dof.second ) = x(dof.first.first,dof.first.second);
        //toc("hbfTOFeel::operator()");
    }
    sync(u,"=");
    u.close();
    return u;

}


//
//HbfFineToCoarse
// 

HbfFineToCoarse::HbfFineToCoarse(int nx, int ny,int m): M_rows(ny), M_cols(nx),SizeElem(m)
{
    tic(); 
    int itilde=-1;
    int jtilde=-1;
    for( int i = 0; i < ny; ++i )
    {
        if(i%m==0)
            for( int j = 0; j < nx; ++j )
             {
                if(j%m==0)
                {
                    itilde=i/m;
                    jtilde=j/m;
                }
                if ( itilde != -1 )
                {
                    M_relation.insert( dof_relation( std::make_pair(i,j), std::make_pair(itilde,jtilde) ) ); 
                }
                itilde=-1;
                jtilde=-1;

            }
    }
    toc("structured thin 2 structured rough");
}

holo3_image<float>
HbfFineToCoarse::operator()( holo3_image<float> const& u, std::string way )
{
    if ( way == std::string("F2C") )
    {
        holo3_image<float> x(M_rows/SizeElem,M_cols/SizeElem);
        float tmp=0;    
        for( auto dof : M_relation.left )
        {
            for( int i = 0; i < SizeElem; ++i )
            {
                for( int j = 0; j < SizeElem; ++j )
                {
                    if (dof.first.first+i<M_rows && dof.first.second+j<M_cols)
                    tmp+=u(dof.first.first+i,dof.first.second+j);
                }
            }
            if ((dof.second.first< M_rows/SizeElem) && (dof.second.second<M_cols/SizeElem))
                x( dof.second.first, dof.second.second ) = tmp/(SizeElem*SizeElem);
                //x( dof.second.first, dof.second.second ) = u(dof.first.first,dof.first.second );
                tmp=0;
         }
         return x;
    }
    else if ( way == std::string("C2F"))
         {
            holo3_image<float> x( M_rows, M_cols ); 
            for( int i = 0; i < M_rows; ++i )
            {
                for( int j = 0; j < M_cols; ++j )
                { 
                    x(i,j)=u(i/SizeElem,j/SizeElem);        
                    //x( dof.first.first, dof.first.second ) = u(dof.second.first,dof.second.second);
                }
             }

             return x;
         }
         else
         {
            printf("No known way");
            holo3_image<float> x(1,1);    
            return x;
         }
};

holo3_image<float>
HbfFineToCoarse::integ( holo3_image<float> const& u, std::string way, double pixel )
{ 
    holo3_image<float> xt(M_rows,M_cols);
    holo3_image<float> x(M_rows/SizeElem,M_cols/SizeElem);
    int i=0;
    int j=0;
    double tmp=0;
    for(i=0;i<M_cols;i++)
    {
        for(j=0;j<M_rows;j++)
        {
            xt(i,j)=u(i,j)*pixel*pixel;
        }
    }
    for( auto dof : M_relation.left )
    {
        for( int i = 0; i < SizeElem; ++i )
        {
            for( int j = 0; j < SizeElem; ++j )
            {
                if (dof.first.first+i<M_rows && dof.first.second+j<M_cols)
                    tmp+=xt(dof.first.first+i,dof.first.second+j);
            }
         }
         if ((dof.second.first< M_rows/SizeElem) && (dof.second.second<M_cols/SizeElem))
            x( dof.second.first, dof.second.second ) = tmp/(SizeElem*SizeElem);
            //x( dof.second.first, dof.second.second ) = u(dof.first.first,dof.first.second );
            tmp=0;
     }
     return x;
}


int TransImage::T(holo3_image<float> im, std::pair<double,double> c)
{
    double x = c.first;
    double y = c.second;

    int i = x/dx;
    int j = y/dy;
    /*
    auto Xhc = Pch<1>();
    Hbf2Feelpp h2f(im.cols(),im.rows,Xhc);
    */
    return im(i,j);
}

int TransImage::T(holo3_image<float> im, std::pair<double,double> c, int L)
{
    double x = c.first;
    double y = c.second;

    int i = x/dx;
    int j = y/dy;
    /*
    auto Xhc = Pch<1>();
    Hbf2Feelpp h2f(im.cols(),im.rows,Xhc);
    */
    return im(L*i,L*j);
}
    //Test table

ElemFineToCoarse::ElemFineToCoarse(int nx, int ny,int N): M_rows(ny), M_cols(nx),SizeElem(N)
{
    /*
    tic(); 
    int itilde=-1;
    int jtilde=-1;
    for( int i = 0; i < ny; ++i )
    {
        itilde=i%N;
        for( int j = 0; j < nx; ++j )
        {
            jtilde=j%N;
            if ( itilde != -1 )
            {
                M_relation.insert( dof_relation( std::make_pair(i,j), std::make_pair(itilde,jtilde) ) ); 
            }

            jtilde=-1;

        }
        itilde=-1;
    }
    toc("structured super-pixel to fine grid ");
    */
    tic();
    int i=0;
    int j=0;
    int m=1;
    for (j=0;j<M_cols;j+=SizeElem)
    {
        for (i=0;i<M_rows;i+=SizeElem)
        {
            M_relation.insert( dof_relation( std::make_pair(i,j), m ) ); 
            m++;
        }        
    }
    toc("structured super-pixel to fine grid ");

}

std::pair<double,double> ElemFineToCoarse::operator()( std::pair<int,int> c)
{
    return std::make_pair(-1.+c.first%SizeElem*(2./SizeElem),-1.+c.second%SizeElem*(2./SizeElem));
}

std::pair<int,int> ElemFineToCoarse::operator()( std::pair<double,double> c, int num)
{
    auto kx=M_cols/SizeElem;
    auto ky=M_rows/SizeElem;
    auto ki=(num-1)%kx;
    auto kj=(num-1)/kx;
    return std::make_pair(c.first+ki*SizeElem,c.second+kj*SizeElem);
}

std::pair<int,int> ElemFineToCoarse::operator()( int num)
{
    auto kx=M_cols/SizeElem;
    auto ky=M_rows/SizeElem;
    auto ki=(num-1)%kx;
    auto kj=(num-1)/kx;

    return std::make_pair(ki*SizeElem,kj*SizeElem);
}


}// Feel
