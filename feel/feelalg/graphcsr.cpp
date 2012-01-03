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



GraphCSR::self_ptrtype
GraphCSR::transpose()
{
    if (M_graphT) return M_graphT;

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
GraphCSR::close()
{
    if ( M_is_closed )
    {
        std::cout << "already closed graph " << this << "...\n";
        return ;
    }
    M_is_closed = true;

    std::cout << "closing graph " << this << "...\n";
    boost::timer ti;
    Debug(5050) << "[close] nrows=" << this->size() << "\n";
    Debug(5050) << "[close] firstRowEntryOnProc()=" << this->firstRowEntryOnProc() << "\n";
    Debug(5050) << "[close] lastRowEntryOnProc()=" << this->lastRowEntryOnProc() << "\n";
    Debug(5050) << "[close] firstColEntryOnProc()=" << this->firstColEntryOnProc() << "\n";
    Debug(5050) << "[close] lastColEntryOnProc()=" << this->lastColEntryOnProc() << "\n";
    Debug(5050) << "[close] M_n_total_nz=" << M_n_total_nz.size() << "\n";
    Debug(5050) << "[close] M_storage size=" << M_storage.size() << "\n";
    M_n_total_nz.resize( M_storage.size() );
    M_n_nz.resize( M_storage.size() );
    M_n_oz.resize( M_storage.size() );

    std::fill( M_n_nz.begin(), M_n_nz.end(), 0 );
    std::fill( M_n_oz.begin(), M_n_oz.end(), 0 );

    size_type sum_nz = 0;
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
// parallel version coming soon....
#if 0
                            if ( (*vecit < firstColEntryOnProc()) ||
                                 (*vecit > lastColEntryOnProc() ))
                                {
                                    //Debug() << "globalindex=" << globalindex << " localindex="
                                    //<< localindex << " off-block diag: " << M_n_oz[localindex] << "\n";
                                    // entry is off block-diagonal
                                    ++M_n_oz[localindex];
                                }
                            else
#endif
                                {
                                    //Debug() << "globalindex=" << globalindex << " localindex="
                                    //<< localindex << " on-block diag: " << M_n_nz[localindex] << "\n";
                                    // entry is in block-diagonal
                                    ++M_n_nz[localindex];
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
                    // something to do here ?
                }


        }
#if 1
    //std::cout << "sum_nz = " << sum_nz << "\n";
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
#endif // 0

} // close


void
GraphCSR::showMe( std::ostream& __out ) const
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
}


} // Feel
