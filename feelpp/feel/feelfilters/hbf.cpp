/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2014-06-04

Copyright (C) 2014-2020 Feel++ Consortium

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
#include <feel/feelfilters/hbf.hpp>
#include <map>

namespace Feel
{

std::tuple<int32_t, int32_t, std::string>
readHBFHeaderAndSizes( std::string const& s )
{
    LOG(INFO) << "Reading Header and sizes " << s << std::endl;
    std::ifstream in( s, std::ios::binary );
    if ( !in )
    {
        using namespace std::string_literals;
        throw std::invalid_argument ( "readHBFHeaderAndSizes - Error opening file "s + s.c_str() );
    }
    return readHBFHeaderAndSizes( in );
}
std::tuple<int32_t, int32_t, std::string>
readHBFHeaderAndSizes( std::istream& in )
{
    LOG(INFO) << "Reading Header and sizes" << std::endl;
    std::string header;
    std::string data_type;
    char ch;
    size_t count = 0;
    bool inside_type = false;
    while ((ch = in.get()) != '\0')
    {
        header += ch;

        if ( ch == '>' || ch == 'c' )
            inside_type = false;
        if ( inside_type )
            data_type += ch;
        if ( ch == '<' )
            inside_type = true;

        ++count;
    }
    LOG(INFO) << "header: " << header << std::endl;
    LOG(INFO) << "data type:" << data_type << std::endl;
    LOG(INFO) << count << std::endl;

    int32_t version = 0;
    in.read( (char*)&version, sizeof( int32_t ) );
    LOG(INFO) << "version: " << version;
    int32_t rows=0, cols=0;
    in.read( (char*)&rows, sizeof( int32_t ) );
    in.read( (char*)&cols, sizeof( int32_t ) );
    LOG(INFO) << "rows: " << rows << " , cols: " << cols << " , size: " << rows*cols << std::endl;
    return {rows, cols, data_type};
}


template <typename T>
holo3_image<T>
readHBF( std::istream& in )
{
    auto [rows,cols, data_type] = readHBFHeaderAndSizes( in );

    holo3_image<T> x ( rows, cols );
    in.read( (char*)x.data(), x.size()*sizeof(T) );
    if (x.rows() <= 6 && x.cols() <= 6)
        LOG(INFO) << x << std::endl;
    LOG(INFO) << "[readhbf] rows: " << x.rows() 
                << ", cols: " << x.cols()
                << ", data type: " << data_type << std::endl;
    return x;
}
template <typename T>
holo3_image<T>
readHBF( std::string const& s )
{
    std::ifstream in( s, std::ios::binary );
    if ( !in )
    {
        using namespace std::string_literals;
        throw std::invalid_argument ( "ReadHBF - Error opening file "s + s.c_str() );
    }
    return readHBF<T>( in );
}
template <typename T>
void writeHBF( std::ostream& out, holo3_image<T> const& x )
{
    std::map<std::string, std::string> data_sizes;
    data_sizes[typeid( int8_t ).name()] = "__int8";
    data_sizes[typeid( int16_t ).name()] = "__int16";
    data_sizes[typeid( int32_t ).name()] = "__int32";
    data_sizes[typeid( int64_t ).name()] = "__int64";
    data_sizes[typeid( long ).name()] = "long";
    data_sizes[typeid( uint8_t ).name()] = "unsigned __int8";
    data_sizes[typeid( uint16_t ).name()] = "unsigned __int16";
    data_sizes[typeid( uint32_t ).name()] = "unsigned __int32";
    data_sizes[typeid( uint64_t ).name()] = "unsigned __int64";
    data_sizes[typeid( float ).name()] = "float";
    data_sizes[typeid( double ).name()] = "double";

    std::string data_type = data_sizes[typeid( T ).name()];
    std::string header = "class CH3Array2D<" + data_type + "> *";

    out.write( header.c_str(), header.size() );
    out.write( "\0", sizeof( char ) );

    LOG( INFO ) << "[writeHBF] write header " << header << "\n";

    int32_t version = 100;
    out.write( (char*)&version, sizeof( int32_t ) );
    int32_t rows = x.rows(), cols = x.cols();
    out.write( (char*)&rows, sizeof( int32_t ) );
    out.write( (char*)&cols, sizeof( int32_t ) );
    LOG( INFO ) << "[writeHBF] rows: " << rows << " , cols: " << cols << " , size: " << rows * cols << std::endl;
    //holo3_image<T> y ( cols, rows );
    //y = x.transpose();
    out.write( (char*)x.data(), x.size() * sizeof( T ) );

    LOG( INFO ) << "[writeHBF]: array written to disk" << std::endl;

    if ( x.rows() <= 6 && x.cols() <= 6 )
        LOG( INFO ) << x << std::endl;
    LOG( INFO ) << "[writehbf] x.rows: " << x.rows() << " , x.cols(): " << x.cols() << std::endl;
}
template <typename T>
void
writeHBF( std::string const& s, holo3_image<T> const& x, worldcomm_ptr_t wc )
{
    if ( wc->localRank() == 0 )
    {
#if 0  
        std::ofstream out( Environment::expand(s).c_str(), std::ios::binary );
#else 
        std::ofstream out( s, std::ios::binary );
#endif    
        if ( !out )
        {
            using namespace std::string_literals;
            throw std::invalid_argument ( "writeHBF - Error opening file "s + s.c_str() );
        }
        writeHBF( out, x );
    }
    wc->barrier();
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
    toc("structured 2 feelpp relation",FLAGS_v>0);
}

Hbf2Feelpp::q1_element_type
Hbf2Feelpp::operator()( holo3_image<float> const& x )
{
    q1_element_type u( M_Xh );
    for( auto const& dof : M_relation.left )
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
Hbf2FeelppStruc::Hbf2FeelppStruc( int nx, int ny, q1_space_ptrtype Yh )
    :
    M_rows(ny), 
    M_cols(nx), 
    M_Xh( Yh )
{
    tic();
    LOG(INFO) << "Hbf2FeelppStruc relation creation \n";
    LOG(INFO) << nx << " : " << ny <<std::endl;
    
    int procSize = Yh->worldComm().localSize();
    rank_type partId = Yh->worldComm().localRank(); 
    int *cx = new int[procSize+1];
    for (int tmpProc=0;tmpProc<=procSize;tmpProc++)
        cx[tmpProc]=(tmpProc)*(nx-1)/procSize;

    tic();
    auto [dof2pid, pid2dof] = Yh->dof()->pointIdToDofRelation("", false, true );
    toc("pidToDof relation",FLAGS_v>0);
    for( int i = std::max(cx[partId]-1,0); i <= std::min(cx[partId+1]+1,nx-1); ++i )
    {
        for( int j = 0; j < ny; ++j )
        {
            M_relation.push_back( dof_relation( std::make_pair(j,i), pid2dof[(ny)*i+j] ));
        }
    }
    toc("structured to feelpp relation",FLAGS_v>0);
    delete [] cx;
}

Hbf2FeelppStruc::q1_element_type
Hbf2FeelppStruc::operator()( holo3_image<float> const& x )
{
    tic();
    q1_element_type u = M_Xh->element();
    toc("h2f", FLAGS_v>0);
    tic();
    for( auto const & dof : M_relation.left )
    {
        u( dof.second ) = x(dof.first.first,dof.first.second);
    }
    toc("h2f dof", FLAGS_v>0);
    tic();
    sync(u,"=");
    u.close();
    toc("h2f sync+close()", FLAGS_v>0);
    return u;


}

holo3_image<float>
Hbf2FeelppStruc::operator()( q1_element_type const& u )
{
    tic();
    holo3_image<float> y( M_rows, M_cols );
    y = holo3_image<float>::Zero( M_rows, M_cols );
    std::vector<int> sizes;
    mpi::all_gather( u.worldComm(), (int)u.dof()->nLocalDofWithoutGhost(), sizes );
    
    for( auto const& dof : M_relation.left )
    {
        DCHECK( dof.first.first < M_rows ) << "invalid row index " << dof.first.first;
        DCHECK( dof.first.second < M_cols ) << "invalid col index " << dof.first.second;
        if ( dof.second < u.size() )
            y( dof.first.first, dof.first.second ) = u(dof.second);
    }
    int p = u.worldComm().rank();
    int pm1 = (p==0)?0:p-1;
    int s = 0;
    for( int i = 0; i < p; ++i )
        s += sizes[i];
    holo3_image<float> x = y.transpose();
    holo3_image<float> yy= y.transpose();
    mpi::gatherv( u.worldComm(), yy.data()+s, sizes[p], 
                  x.data(), sizes, 0 );
    toc("H2F feelpp to holo3_image",FLAGS_v>0);
    return x.transpose();
}

Hbf2FeelppStruc::q1_element_type
Hbf2FeelppStruc::operator()( holo3_image<float> const& x, q1_element_type& u)
{
    //q1_element_type u( M_Xh );
    //std::cout << "u.nDof() = " << u.nDof() << std::endl;
    //std::cout << "M_relation.size() = " << M_relation.size() << std::endl;
//#pragma omp parallel for
    for( auto const & dof : M_relation.left )
    {
        //LOG(INFO) << "dof.second = " << dof.second << " = " << dof.first.first << " : " << dof.first.second << std::endl;
        //tic();
        if ( dof.second < u.size() )
            u( dof.second ) = x(dof.first.first,dof.first.second);
        //toc("hbfTOFeel::operator()");
    }
    sync(u,"=");
    u.close();
    return u;

}

void
Hbf2FeelppStruc::cellToQ1( holo3_image<float> const& x, q1_element_type& u)
{
    //q1_element_type u( M_Xh );
    //std::cout << "u.nDof() = " << u.nDof() << std::endl;
    //std::cout << "M_relation.size() = " << M_relation.size() << std::endl;
//#pragma omp parallel for
    for( auto const & dof : M_relation.left )
    {
        //LOG(INFO) << "dof.second = " << dof.second << " = " << dof.first.first << " : " << dof.first.second << std::endl;
        //tic();
        if ( dof.second < u.size() )
            u( dof.second ) = x(dof.first.first,dof.first.second);
        //toc("hbfTOFeel::operator()");
    }
    sync(u,"=");
    u.close();

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


template holo3_image<double>   readHBF(std::string const&);
template holo3_image<float>    readHBF(std::string const&);
template holo3_image<int8_t>   readHBF(std::string const&);
template holo3_image<int16_t>  readHBF(std::string const&);
template holo3_image<int32_t>  readHBF(std::string const&);
template holo3_image<int64_t>  readHBF(std::string const&);
template holo3_image<uint8_t>  readHBF(std::string const&);
template holo3_image<uint16_t> readHBF(std::string const&);
template holo3_image<uint32_t> readHBF(std::string const&);
template holo3_image<uint64_t> readHBF(std::string const&);

template holo3_image<double>   readHBF(std::istream&);
template holo3_image<float>    readHBF(std::istream&);
template holo3_image<int8_t>   readHBF(std::istream&);
template holo3_image<int16_t>  readHBF(std::istream&);
template holo3_image<int32_t>  readHBF(std::istream&);
template holo3_image<int64_t>  readHBF(std::istream&);
template holo3_image<uint8_t>  readHBF(std::istream&);
template holo3_image<uint16_t> readHBF(std::istream&);
template holo3_image<uint32_t> readHBF(std::istream&);
template holo3_image<uint64_t> readHBF(std::istream&);

template void writeHBF( std::string const& , holo3_image<double> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<float> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<int8_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<int16_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<int32_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<int64_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<uint8_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<uint16_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<uint32_t> const&, worldcomm_ptr_t wc );
template void writeHBF( std::string const& , holo3_image<uint64_t> const&, worldcomm_ptr_t wc );

template void writeHBF( std::ostream& , holo3_image<double>   const& );
template void writeHBF( std::ostream& , holo3_image<float>    const& );
template void writeHBF( std::ostream& , holo3_image<int8_t>   const& );
template void writeHBF( std::ostream& , holo3_image<int16_t>  const& );
template void writeHBF( std::ostream& , holo3_image<int32_t>  const& );
template void writeHBF( std::ostream& , holo3_image<int64_t>  const& );
template void writeHBF( std::ostream& , holo3_image<uint8_t>  const& );
template void writeHBF( std::ostream& , holo3_image<uint16_t> const& );
template void writeHBF( std::ostream& , holo3_image<uint32_t> const& );
template void writeHBF( std::ostream& , holo3_image<uint64_t> const& );

}// Feel
