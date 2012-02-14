/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-10-29

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file graphcsr.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-10-29
 */
#include <boost/timer.hpp>

#include <feel/feelcore/application.hpp>
#include <feel/feelalg/graphcsr.hpp>

namespace Feel
{
GraphCSR::GraphCSR( size_type n,
                    size_type first_row_entry_on_proc,
                    size_type last_row_entry_on_proc,
                    size_type first_col_entry_on_proc,
                    size_type last_col_entry_on_proc )
    :
    M_first_row_entry_on_proc( first_row_entry_on_proc ),
    M_last_row_entry_on_proc( last_row_entry_on_proc ),
    M_first_col_entry_on_proc( first_col_entry_on_proc ),
    M_last_col_entry_on_proc( last_col_entry_on_proc ),
    M_max_nnz( 0 ),
    M_n_total_nz( n, 0 ),
    M_n_nz( n, 0 ),
    M_n_oz( n, 0 ),
    M_storage(),
    M_is_closed( false )
{
    //std::cout << "creating graph " << this << "\n";
}

GraphCSR::GraphCSR( GraphCSR const & g )
    :
    M_first_row_entry_on_proc( g.M_first_row_entry_on_proc ),
    M_last_row_entry_on_proc( g.M_last_row_entry_on_proc ),
    M_first_col_entry_on_proc( g.M_first_col_entry_on_proc ),
    M_last_col_entry_on_proc( g.M_last_col_entry_on_proc ),
    M_max_nnz( g.M_max_nnz ),
    M_n_total_nz( g.M_n_total_nz ),
    M_n_nz( g.M_n_nz ),
    M_n_oz( g.M_n_oz ),
    M_storage( g.M_storage ),
    M_graphT( g.M_graphT),
    M_is_closed( g.M_is_closed )
{
}

GraphCSR::~GraphCSR()
{
    M_storage.clear();
}

GraphCSR&
GraphCSR::operator=( GraphCSR const& g )
{
    if ( this != &g )
        {
            M_first_row_entry_on_proc = g.M_first_row_entry_on_proc;
            M_last_row_entry_on_proc = g.M_last_row_entry_on_proc;
            M_first_col_entry_on_proc = g.M_first_col_entry_on_proc;
            M_last_col_entry_on_proc = g.M_last_col_entry_on_proc;
            M_max_nnz = g.M_max_nnz;
            M_n_total_nz= g.M_n_total_nz;
            M_n_nz = g.M_n_nz;
            M_n_oz = g.M_n_oz;
            M_storage = g.M_storage;
            M_graphT = g.M_graphT;
            M_is_closed = g.M_is_closed;
        }
    return *this;
}


void
GraphCSR::zero()
{
#if 0
    auto nbDof = M_last_row_entry_on_proc-M_first_row_entry_on_proc+1;
    for (size_type i=0 ; i<nbDof ; ++i)
        {
            row_type& row = this->row(i);
            row.get<0>() = 0;//proc
            row.get<1>() = i; //local index : warning false in parallel!!!!
            row.get<2>().clear(); //all is zero
        }
#else
    auto nbDof = M_last_row_entry_on_proc-M_first_row_entry_on_proc+1;
    for (size_type i=M_first_row_entry_on_proc ; i<=M_last_row_entry_on_proc ; ++i)
        {
            row_type& row = this->row(i);
            row.get<0>() = 0;//proc
            row.get<1>() = i-M_first_row_entry_on_proc; //local index
            row.get<2>().clear(); //all is zero
        }

#endif
}

GraphCSR::self_ptrtype
GraphCSR::transpose()
{

    if (M_graphT) return M_graphT;
    this->close();
    M_graphT = self_ptrtype(new self_type( M_n_total_nz.size(),
                                           M_first_col_entry_on_proc,
                                           M_last_col_entry_on_proc,
                                           M_first_row_entry_on_proc,
                                           M_last_row_entry_on_proc ) );

    for( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
        {
            // Get the row of the sparsity pattern
            row_type const& irow = it->second;
            if ( boost::get<0>( irow ) == M_comm.rank() )
                {
                    // num line
                    size_type globalindex = it->first;

                    size_type localindex = boost::get<1>( irow );
                    for ( auto colit = boost::get<2>( irow ).begin(), colen=boost::get<2>( irow ).end() ; colit!=colen ; ++colit )
                        {
                            self_type::row_type& row = M_graphT->row(*colit);
                            row.get<0>()=irow.get<0>();
                            // Warning : wrong in parallele
                            row.get<1>()=*colit;//globalindex;
                            row.get<2>().insert(globalindex);
                        }
                }
        }

    M_graphT->close();

    return M_graphT;
}

void
GraphCSR::addMissingZeroEntriesDiagonal()
{
    size_type m = this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1;
    size_type n = this->lastColEntryOnProc()-this->firstColEntryOnProc()+1;

    for ( size_type i = this->firstRowEntryOnProc() ; i< this->firstRowEntryOnProc()+std::min(m,n) ; ++i)
        {
            if (this->storage().find(i)!=this->end())
                {
                    if (this->row(i).get<2>().find(i) == this->row(i).get<2>().end())
                        {
                            this->row(i).get<2>().insert(i);
                        }
                }
            else
                {
                    this->row(i).get<0>() = 0;//rank
                    this->row(i).get<1>() = i;//loc
                    this->row(i).get<2>().insert(i);
                }
        }
}


void
GraphCSR::close()
{
    if ( M_is_closed )
    {
       //std::cout << "already closed graph " << this << "...\n";
        return ;
    }
    M_is_closed = true;

    //std::cout << "closing graph " << this << "...\n";
    boost::timer ti;
    Debug(5050) << "[close] nrows=" << this->size() << "\n";
    Debug(5050) << "[close] firstRowEntryOnProc()=" << this->firstRowEntryOnProc() << "\n";
    Debug(5050) << "[close] lastRowEntryOnProc()=" << this->lastRowEntryOnProc() << "\n";
    Debug(5050) << "[close] firstColEntryOnProc()=" << this->firstColEntryOnProc() << "\n";
    Debug(5050) << "[close] lastColEntryOnProc()=" << this->lastColEntryOnProc() << "\n";
    Debug(5050) << "[close] M_n_total_nz=" << M_n_total_nz.size() << "\n";
    Debug(5050) << "[close] M_storage size=" << M_storage.size() << "\n";
    Debug(5050) << "[close] nrows=" << this->size() << "\n";
#if !defined(FEEL_ENABLE_MPI_MODE)
    M_n_total_nz.resize( M_last_row_entry_on_proc+1/*M_storage.size()*/ );
    M_n_nz.resize( M_last_row_entry_on_proc+1/*M_storage.size()*/ );
    M_n_oz.resize( M_last_row_entry_on_proc+1/*M_storage.size()*/ );
#else // MPI
    M_n_total_nz.resize( this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1/*M_storage.size()*/ );
    M_n_nz.resize( this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1/*M_storage.size()*/ );
    M_n_oz.resize( this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1/*M_storage.size()*/ );

    std::vector<int> nbMsgToSend(this->comm().size());
    std::fill(nbMsgToSend.begin(),nbMsgToSend.end(),0);
    std::vector< std::map<int,int> > mapMsg(this->comm().size());
#endif

    std::fill( M_n_nz.begin(), M_n_nz.end(), 0 );
    std::fill( M_n_oz.begin(), M_n_oz.end(), 0 );

    size_type sum_nz = 0;
    //size_type sum_n_nz = 0;
    //size_type sum_n_oz = 0;

    M_max_nnz = 0;
    for( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
        {
            // Get the row of the sparsity pattern
            row_type const& irow = it->second;
            if ( boost::get<0>( irow ) == M_comm.rank() )
                {
                    //std::vector<size_type> const& ivec = boost::get<2>( irow );
                    size_type vec_size = boost::get<2>( irow ).size();
                    size_type globalindex = it->first;
                    size_type localindex = boost::get<1>( irow );


                    FEEL_ASSERT( globalindex >= firstRowEntryOnProc() )
                        ( globalindex <= lastRowEntryOnProc() )
                        ( globalindex )( firstRowEntryOnProc() )
                        ( lastRowEntryOnProc() ).error ( "invalid local/global index" );
                    FEEL_ASSERT( globalindex >= 0 )( globalindex < M_n_total_nz.size() )
                        ( globalindex )
                        ( M_n_total_nz.size() ).error ( "invalid local/global index for M_n_total_nz" );
                    M_n_total_nz[localindex] = vec_size;
                    sum_nz += vec_size;
                    for( auto vecit = boost::get<2>( irow ).begin(), vecen = boost::get<2>( irow ).end(); vecit != vecen; ++vecit )
                        {
#if defined(FEEL_ENABLE_MPI_MODE) // MPI
                            if ( (*vecit < firstColEntryOnProc()) ||
                                 (*vecit > lastColEntryOnProc() ))
                                {
                                    //Debug() << "globalindex=" << globalindex << " localindex="
                                    //<< localindex << " off-block diag: " << M_n_oz[localindex] << "\n";
                                    //std::cout << "rank " << M_comm.rank() << " globalindex=" << globalindex << " localindex="
                                    //          << localindex << " off-block diag: " << M_n_oz[localindex] << std::endl;
                                    // entry is off block-diagonal
                                    ++M_n_oz[localindex];
                                    //vincent++sum_n_oz;
                                }
                            else
#endif
                                {
                                    //Debug() << "globalindex=" << globalindex << " localindex="
                                    //<< localindex << " on-block diag: " << M_n_nz[localindex] << "\n";
                                    //std::cout << "rank " << M_comm.rank() << " globalindex=" << globalindex << " localindex="
                                    //          << localindex << " on-block diag: " << M_n_nz[localindex] << std::endl;
                                    // entry is in block-diagonal
                                    ++M_n_nz[localindex];
                                    //vincent++sum_n_nz;
                                }
                        }


#if !defined( NDEBUG )
                    Debug( 5050 ) << "M_total_nz [  " << localindex << "]=" << M_n_total_nz[localindex] << "\n";

                    Debug( 5050 ) << "M_nz [  " << localindex << "]=" << M_n_nz[localindex] << "\n";
                    Debug( 5050 ) << "M_oz [  " << localindex << "]=" << M_n_oz[localindex] << "\n";
#endif // NDEBUG

                    M_max_nnz = std::max( M_n_total_nz[localindex], M_max_nnz );
                }
            else
                {
#if defined(FEEL_ENABLE_MPI_MODE) // MPI

                    //size_type globalindex = it->first;
                    //size_type localindex = boost::get<1>( irow );
                    //size_type vec_size = boost::get<2>( irow ).size();

                    auto dofOnGlobalCluster = it->first;

                    // Get the row of the sparsity pattern
                    row_type const& irow = it->second;

                    auto procOnGlobalCluster = irow.get<0>();

                    auto dofOnProc = irow.get<1>();
                    //std::cout << " OOOOdofOnProc " << dofOnProc << " vecsize" << boost::get<2>( irow ).size() << std::endl;
                    std::vector<size_type> ivec( irow.get<2>().size()+2 );
                    ivec[0]=dofOnGlobalCluster;
                    ivec[1]=dofOnProc;

                    auto icol_it = irow.get<2>().begin();
                    auto icol_en = irow.get<2>().end();
                    for (int i=0; icol_it!=icol_en ; ++i,++icol_it)
                        {
                            //std::cout << " *icol_it " << *icol_it << std::endl;
                            ivec[i+2] = *icol_it;//-this->firstColEntryOnProc(); //ivec[i+1] = *icol_it;
                        }

#if 0
                    std::cout << "/n I am proc " << world.rank() << " I send to " << procOnGlobalCluster
                              << " with dofGlobCluster " << dofOnGlobalCluster
                              << " with dofLoc " << dofOnProc
                              << " size of ivec " << ivec.size()
                              << std::endl;
#endif
                    this->comm().send(procOnGlobalCluster, nbMsgToSend[procOnGlobalCluster],ivec );

                    mapMsg[procOnGlobalCluster].insert(std::make_pair(nbMsgToSend[procOnGlobalCluster], irow.get<1>()  ) );

                    ++nbMsgToSend[procOnGlobalCluster];
#endif // MPI
                }


        }

#if defined(FEEL_ENABLE_MPI_MODE) // MPI
    // counter of msg received for each process
    std::vector<int> nbMsgToRecv;
    mpi::all_to_all(this->comm(),
                    nbMsgToSend,
                    nbMsgToRecv);

    // recv id asked and re-send set of face id
    for (int proc=0; proc<this->comm().size();++proc)
        {
            for ( int cpt=0;cpt<nbMsgToRecv[proc];++cpt)
                {
                    //size_type localindex = mapMsg[proc][cpt];
                    std::vector<size_type> ivec;

                    this->comm().recv(proc, cpt,ivec );

                    size_type globalindex = ivec[0];
                    size_type localindex = ivec[1];
                    row_type& row = this->row(globalindex);
                    row.get<0>()=this->comm().rank();
                    row.get<1>()=localindex;

                    //std::cout << " OOOOlocalindex " << localindex << std::endl;
                    for (size_type j=1; j< ivec.size(); j++)
                        {
                            bool isInserted = false;
                            if (row.get<2>().find(ivec[j]) == row.get<2>().end())
                                {
                                    row.get<2>().insert(ivec[j]);
                                    isInserted = true;
                                }

                            if ( (ivec[j] < firstColEntryOnProc()) ||
                                 (ivec[j] > lastColEntryOnProc() ))
                                {
                                    // SUREMENT UTILE!!!!!
                                    // M_colindices_nz[localindex].insert(ivec[j]);
                                    //std::cout << "AYTYATYATYAT"<<std::endl;
                                    //std::cout << "globalindex=" << globalindex << " localindex="
                                    //          << localindex << " off-block diag: " << M_n_oz[localindex] << std::endl;
                                    // entry is off block-diagonal
                                    //if (this->storage().find(localindex)==this->end())
                                    //    if (this->storage().find(localindex)->get<2>().find(ivec[j]) ==this->end())
                                    if (isInserted) { ++M_n_oz[localindex]; ++sum_nz; }
                                }
                            else
                                {
                                    //std::cout << "rank " << this->comm().rank() << " AIAIAIAI localindex"<< localindex << std::endl;
                                    //if (this->storage().find(localindex)==this->end())
                                    if (isInserted) { ++M_n_nz[localindex]; ++sum_nz; }

#if 0
                                    auto isInserted = M_colindices_nz[localindex].insert(ivec[j]).second;
                                    if (isInserted)
                                        ++M_n_nz[localindex];
#endif
                                }


                        }
                    //std::sort(M_colindices_nz[localindex].begin(), M_colindices_nz[localindex].end());

                }

        } // for (int proc=0; proc<M_comm.size();++proc)


#endif


#if 1
#if 0
    M_ia.resize( M_storage.size()+1 );
    M_ja.resize( sum_nz );
    M_a.resize(  sum_nz, 0. );
    size_type col_cursor = 0;
    auto jait = M_ja.begin();
    for( auto it = M_storage.begin(), en = M_storage.end()  ; it != en; ++it )
    {
        row_type const& irow = it->second;
        size_type localindex = boost::get<1>( irow );
        M_ia[localindex] = col_cursor;
        jait = std::copy( boost::get<2>( irow ).begin(), boost::get<2>( irow ).end(), jait );
        col_cursor+=boost::get<2>( irow ).size();
    }
    M_ia[M_storage.size()] = sum_nz;
#else

#if !defined(FEEL_ENABLE_MPI_MODE) // NOT MPI
    M_ia.resize( M_last_row_entry_on_proc+2,0 );
    M_ja.resize( sum_nz );
    M_a.resize(  sum_nz, 0. );
    size_type col_cursor = 0;
    auto jait = M_ja.begin();
    //for( auto it = M_storage.begin(), en = M_storage.end()  ; it != en; ++it )

    for( int i = 0 ; i<(M_last_row_entry_on_proc+1); ++i )
    {
        if (M_storage.find(i)!=M_storage.end())
            {
                row_type const& irow = this->row(i);
                size_type localindex = boost::get<1>( irow );
                M_ia[localindex] = col_cursor;
                jait = std::copy( boost::get<2>( irow ).begin(), boost::get<2>( irow ).end(), jait );
                col_cursor+=boost::get<2>( irow ).size();
            }
        else
            {
                M_ia[i] = col_cursor;
            }
    }
    M_ia[M_last_row_entry_on_proc+1] = sum_nz;

#else // MPI
    size_type nRowLoc = this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1;
    M_ia.resize( nRowLoc+1,0 );
    M_ja.resize( /*sum_n_nz*/sum_nz );
    M_a.resize(  /*sum_n_nz*/sum_nz, 0. );
    size_type col_cursor = 0;
    auto jait = M_ja.begin();
    for( int i = 0 ; i<nRowLoc; ++i )
    {
        if (M_storage.find(this->firstRowEntryOnProc()+i)!=M_storage.end())
            {
                row_type const& irow = this->row(this->firstRowEntryOnProc()+i);
                size_type localindex = boost::get<1>( irow );
                M_ia[i/*localindex*/] = col_cursor;
                jait = std::copy( boost::get<2>( irow ).begin(), boost::get<2>( irow ).end(), jait );
                /*for( auto vecit = boost::get<2>( irow ).begin(), vecen = boost::get<2>( irow ).end(); vecit != vecen; ++vecit )
                    {
                        //if ( (*vecit >= firstColEntryOnProc()) && (*vecit <= lastColEntryOnProc() ))
                        {
                            M_ja[col_cursor]=*vecit;std::cout << *vecit << std::endl;
                            ++col_cursor;
                        }
                        }*/
                col_cursor+=boost::get<2>( irow ).size();
            }
        else
            {
                //std::cout << "\n STRANGE " << std::endl;
                M_ia[i] = col_cursor;
            }
    }
    M_ia[nRowLoc] = /*sum_n_nz*/sum_nz;
#endif // MPI
#endif // 0
#endif // 1

} // close


void
GraphCSR::showMe( std::ostream& __out ) const
{
    for (int proc = 0;proc<this->comm().size();++proc)
        {
            if (proc==this->comm().rank())
                {
                    __out << "first_row_entry_on_proc " << M_first_row_entry_on_proc << std::endl;
                    __out << "last_row_entry_on_proc " << M_last_row_entry_on_proc << std::endl;
                    __out << "first_col_entry_on_proc " << M_first_col_entry_on_proc << std::endl;
                    __out << "last_col_entry_on_proc " << M_last_col_entry_on_proc << std::endl;
                    __out << "max_nnz " << M_max_nnz << std::endl;

                    for( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
                        {
                            // Get the row of the sparsity pattern
                            row_type const& row = it->second;

                            __out << " proc " << row.get<0>()
                                  << " globalindex " << it->first
                                  << " localindex " << row.get<1>()
                                  << "(nz " << M_n_nz[row.get<1>()]
                                  << " oz " << M_n_oz[row.get<1>()]
                                  << ") : ";
                            //size_type vec_size = boost::get<2>( irow ).size();
                            //size_type globalindex = it->first;
                            //size_type localindex = boost::get<1>( irow );
                            for (auto it = row.get<2>().begin(), en= row.get<2>().end() ; it!=en ; ++it)
                                __out << *it << " ";

                            __out << std::endl;
                        }
                } // if (proc==this->comm().rank())
            this->comm().barrier();
        } // for (int proc = 0;proc<nProc;++proc)

}

void
GraphCSR::printPython( std::string const& nameFile) const
{

    std::ofstream graphFile(nameFile, std::ios::out /*| std::ios::app*/);
#if 0
    std::cout << "first_row_entry_on_proc " << M_first_row_entry_on_proc << std::endl;
    std::cout << "last_row_entry_on_proc " << M_last_row_entry_on_proc << std::endl;
    std::cout << "first_col_entry_on_proc " << M_first_col_entry_on_proc << std::endl;
    std::cout << "last_col_entry_on_proc " << M_last_col_entry_on_proc << std::endl;
    std::cout << "max_nnz " << M_max_nnz << std::endl;
#endif

    graphFile << "import numpy" << std::endl
              << "from scipy.sparse import * " << std::endl
              << "from scipy import * " << std::endl
              << "from pylab import * " << std::endl;

    graphFile << "nRow=" << M_last_row_entry_on_proc-M_first_row_entry_on_proc+1 << std::endl
              << "nCol=" << M_last_col_entry_on_proc-M_first_col_entry_on_proc+1 << std::endl;

    graphFile << "mattt = array([" << std::endl;
    for( auto it = M_storage.begin(), en = --M_storage.end() ; it != en; ++it )
        {
            auto const& row = it->second;

            for (auto it2 = row.get<2>().begin(), en2= row.get<2>().end() ; it2!=en2 ; ++it2)
                graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ],";// << std::endl;
        }

    auto it = --M_storage.end();
    auto const& row = it->second;

    if (row.get<2>().size()>0)
        {
            if (row.get<2>().size()>1)
                {
                    for (auto it2 = row.get<2>().begin(), en2= --row.get<2>().end() ; it2!=en2 ; ++it2)
                        graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ],";
                }
            auto it2 = --row.get<2>().end();
            graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ] ])" << std::endl;
        }
    else { /*???*/ }

    graphFile << "row = array(mattt[:,0],dtype=int);" << std::endl
              << "col = array(mattt[:,1],dtype=int);" << std::endl
              << "data = array(mattt[:,2]);" << std::endl
              << "A = csr_matrix( (data,(row,col)), shape=(nRow,nCol) );" << std::endl
              << "fig = plt.figure();" << std::endl
              << "matplotlib.pyplot.spy(A,precision=1e-8,aspect='equal');" << std::endl
              << "plt.show();" << std::endl;


} // printPython


} // Feel
