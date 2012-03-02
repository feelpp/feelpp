/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2007-08-14

   Copyright (C) 2008, 2009 Christophe Prud'homme
   Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)
   Copyright (C) 2007 école Polytechnique Fédérale de Lausanne (EPFL)

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
   \file matrixepetra.cpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2007-08-14
*/
#include <feel/feelalg/vectorepetra.hpp>
#include <feel/feelalg/matrixepetra.hpp>

#if defined( HAVE_TRILINOS_EPETRA )
#include <Epetra_FECrsGraph.h>
#include <Epetra_Time.h>
#include <Epetra_RowMatrixTransposer.h>
#endif


namespace Feel
{
#if defined ( HAVE_TRILINOS_EPETRA )

namespace detail
{
Epetra_Map epetraMap( DataMap const& dmap )
{
    std::vector<int> e( dmap.nMyElements() );
    std::copy( dmap.myGlobalElements().begin(),
               dmap.myGlobalElements().end(),
               e.begin() );
    return Epetra_Map( -1, dmap.nMyElements(), e.data(), 0, Epetra_MpiComm( dmap.comm() ) );
}
}
MatrixEpetra::real_type
MatrixEpetra::energy ( vector_type const& v1, vector_type const& v2, bool transpose ) const
{
    this->close();

    real_type res;
    if ( dynamic_cast<epetra_vector_type const*>( &v1 ) != (epetra_vector_type const*)0 )
    {
        epetra_vector_type const& ev1( dynamic_cast<epetra_vector_type const&>( v1 ) );
        epetra_vector_type const& ev2( dynamic_cast<epetra_vector_type const&>( v2 ) );
        epetra_vector_type ev3( this->getRowMap() );

        _M_mat->Multiply( transpose, ev2.vec(), ev3.vec() );
        ev3.vec().Dot( ev1.vec(), &res );
    }
    else
    {
        Epetra_BlockMap bmap1( detail::epetraMap( v1.map() ) );
        VectorEpetra<value_type> u( bmap1 );
        {
            size_type s = u.localSize();
            size_type start = u.firstLocalIndex();
            for( size_type i = 0; i < s; ++i )
                u.set( start + i, v1( start + i ) );
        }
        Epetra_BlockMap bmap2( detail::epetraMap( v2.map() ) );
        VectorEpetra<value_type> v( bmap2 );
        {
            size_type s = v.localSize();
            size_type start = v.firstLocalIndex();
            for( size_type i = 0; i < s; ++i )
                v.set( start + i, v2( start + i ) );
        }
        VectorEpetra<value_type> z( bmap1 );
        _M_mat->Multiply( transpose, v.vec(), z.vec() );
        z.vec().Dot( u.vec(), &res );
    }

    return res;
}

void
MatrixEpetra::init ( const size_type m,
                     const size_type n,
                     const size_type m_l,
                     const size_type /*n_l*/,
                     graph_ptrtype const& graph )
{
    if ((m==0) || (n==0))
        return;

    this->setGraph( graph );

    {
#if 1
        // Clear initialized matrices
        if (this->isInitialized())
            this->clear();

        /*
          _M_emap = Epetra_Map( m, m_l, 0, Epetra_MpiComm(M_comm));
          _M_col_emap = Epetra_Map( n, n, 0, Epetra_MpiComm(M_comm));
          _M_dom_map = _M_emap;
          _M_range_map = _M_emap;
        */
#endif
        //_M_emap.Print(std::cout );
        //_M_col_emap.Print(std::cout );
        //_M_dom_map.Print(std::cout );
        //_M_range_map.Print(std::cout );
        //Epetra_FECrsGraph epetra_graph( Copy, _M_emap, _M_col_emap, (int*)this->graph()->nNz().data() );
        //std::cout << "_M_row_map " << _M_emap << "\n";
        //std::cout << "_M_col_map " << _M_col_emap << "\n";

        boost::shared_ptr<Epetra_Map> graph_map ( new Epetra_Map( _M_emap ) );


        //find which map (between col and row map has more lines)
        if ( _M_emap.MaxAllGID() < _M_col_emap.MaxAllGID() )
            graph_map = boost::shared_ptr<Epetra_Map>(new Epetra_Map(_M_col_emap));

        std::vector<int> nnz( this->graph()->nNz().size() );
        std::copy( this->graph()->nNz().begin(), this->graph()->nNz().end(), nnz.begin() );
        //Epetra_FECrsGraph epetra_graph( Copy, _M_emap, _M_col_emap, (int*)this->graph()->nNz().data() );
        Epetra_FECrsGraph epetra_graph( Copy, _M_emap, _M_col_emap, nnz.data() );

        //Epetra_FECrsGraph epetra_graph( Copy, _M_emap, (int*)this->graph()->nNz().data() );
        //int start = _M_emap.MinMyGID();
        //int stop = _M_emap.MaxMyGID()+1;

        // insert global indices, even the non-local ones
        graph_type::const_iterator it = this->graph()->begin();
        graph_type::const_iterator en = this->graph()->end();
        for( ; it != en; ++it )
            {
                Debug( 10010 ) << "[MatrixEpetra::init] row with gid/lid=("
                               << it->first << "/" << it->second.get<1>() << ")"
                               << " on proc : " << it->second.get<0>()
                               << " with nnz: " << it->second.get<2>().size() << "\n";

                // Insert global indices line by line
                const int numRows = 1;
                const int ii = it->first;
                int rows[1];
                rows[0] =  ii;
                const int numCols = it->second.get<2>().size();

                std::vector<int> cols( it->second.get<2>().size() );
                std::copy( it->second.get<2>().begin(), it->second.get<2>().end(), cols.begin() );
                //const int *cols = (int*)it->second.get<2>().data();

                //std::cout << "Inserted row number: " << ii << " and column: " << *cols << "\n";

                int ierr = epetra_graph.InsertGlobalIndices( numRows, rows,
                                                             numCols, cols.data() );


                Feel::detail::ignore_unused_variable_warning(ierr);

                FEELPP_ASSERT( ierr == 0 )( ierr )( it->first )( it->second.get<0>() )( it->second.get<1>() )( it->second.get<2>().size() ).warn( "problem with Epetra_FECrsGraph::InsertGlobalIndices" );
                //                 Debug( 10010 ) << "row = " << ii << " irow.size " << irow.size() <<  " ierr = " << ierr << "\n";
            }
        //std::cout << "------------------------------------------------------------\n";
        //std::cout << "Epetra graph: " << epetra_graph << "\n";
        //std::cout << "------------------------------------------------------------\n";

        //int ierr = epetra_graph.GlobalAssemble( _M_dom_map, _M_range_map );
        //std::cout << "_M_dom_map " << _M_dom_map << "\n";
        //std::cout << "_M_range_map " << _M_range_map << "\n";

        //std::cout << "Call GlobalAssemble...\n";
        //int ierr = epetra_graph.GlobalAssemble( _M_dom_map, _M_range_map, true );
        int ierr = epetra_graph.GlobalAssemble( _M_range_map, _M_dom_map, true );
        //int ierr = epetra_graph.GlobalAssemble( true );
        //int ierr = epetra_graph.GlobalAssemble();
        //std::cout << "Epetra graph: " << epetra_graph << "\n";
        Feel::detail::ignore_unused_variable_warning(ierr);
        FEELPP_ASSERT( ierr == 0 )( ierr ).warn ( "[MatrixEpetra::init] GlobalAssemble failed" );
        Debug( 10010 ) << "Global assemble  ierr = " << ierr << "\n";
        //epetra_graph.Print( std::cout );
        _M_mat = boost::shared_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix(Copy, epetra_graph ) );

        //_M_mat->Print( std::cout );
        this->setInitialized( true );
    }
}




void
MatrixEpetra::add (const size_type i,
                   const size_type j,
                   const value_type& value)
{
    FEELPP_ASSERT (this->isInitialized()).error( "MatrixEpetra<> not properly initialized" );

    int i_val = static_cast<int>(i);
    int j_val = static_cast<int>(j);
    value_type epetra_value = static_cast<value_type>(value);

    int ierr=0;

#if 0
    //works only for values that already exist in the matrix. Can not be used to put new values into.
    ierr = _M_mat->SumIntoGlobalValues(1, &i_val, 1,  &j_val, &epetra_value);

    if (ierr != 0)
        {
            ierr=_M_mat->InsertGlobalValues(1, &i_val, 1,  &j_val, &epetra_value);

            Debug(10010) << "ERRORCODE SumIntoGlobalValues: " << ierr
                         <<  " in M(" << i_val << "," << j_val << ") for value "
                         << epetra_value << ".\n";
        }
#else
    //ierr=_M_mat->InsertGlobalValues(1, &i_val, 1,  &j_val, &epetra_value);
    ierr = _M_mat->SumIntoGlobalValues(1, &i_val, 1,  &j_val, &epetra_value);

#if 0
    Debug( 10010 ) << "ERRORCODE SumIntoGlobalValues: " << ierr
                   <<  " in M(" << i_val << "," << j_val << ") for value "
                   << epetra_value << ".\n";
#endif
    //ierr=_M_mat->InsertMyValues(i_val, 1,  &j_val, &epetra_value);
#endif
}


void
MatrixEpetra::multiplyMatrix( const MatrixEpetra& A, const MatrixEpetra& B )
{
    EpetraExt::MatrixMatrix::Multiply( A.mat(), false, B.mat(), false, (*this).mat() );
}

void
MatrixEpetra::diagonal (Vector<double>& dest) const
{
    // TBD
}
void
MatrixEpetra::transpose( MatrixSparse<value_type>& Mt ) const
{
    FEELPP_ASSERT( 0 ).warn( "not implemented yet" );
#if 0
 	Epetra_FECrsMatrix* Atrans;

    Epetra_RowMatrixTransposer transposer ( &*_M_mat );
    transposer.CreateTranspose( true, Atrans );

    Mt =  new MatrixEpetra( *Atrans );
    bool isSymmetric;

    if (isSymmetric) {
    } else {

    }
#endif

}


// print into Matlab sparse Matrix
void
MatrixEpetra::printMatlab (const std::string name) const
{
    FEELPP_ASSERT (this->isInitialized()).error( "epetra matrix not properly initialized" );

    FEELPP_ASSERT (this->closed()).warn("epetra matrix not closed");
    if ( !this->closed() )
        const_cast<MatrixEpetra*>(this)->close();

    Debug( 10010 ) << "[printMatlab] print matrix in matlab file " << name << "\n";
    //std::cout << "[printMatlab] print matrix in matlab file " << name << "\n";

    //this->printKonsole();
    //std::cout << "[printMatlab] print matrix in matlab file done\n";
    int ret = EpetraExt::RowMatrixToMatlabFile( name.c_str(), *_M_mat);
    //int ret = EpetraExt::RowMatrixToMatrixMarketFile( name.c_str(), *_M_mat, "toto", "tutu" );
    if ( ret != 0 )
        std::cout << "error in printMatlab\n";
}


// print to console
void
MatrixEpetra::printKonsole () const
{
    FEELPP_ASSERT (this->isInitialized()).error( "epetra matrix not properly initialized" );

    std::cout << "\n+---------------Information about the Matrix---------------+\n" <<  std::endl;
    cout << *_M_mat;
    std::cout << "+----------------------------------------------------------+\n" <<  std::endl;
    std::cout << "\n+---------------Information about the RowMap------------------+\n" <<  std::endl;
    cout << _M_mat->RowMap();
    std::cout << "+----------------------------------------------------------+\n" <<  std::endl;
    std::cout << "\n+---------------Information about the ColMap------------------+\n" <<  std::endl;
    cout << _M_mat->ColMap();
    std::cout << "+----------------------------------------------------------+\n" <<  std::endl;
}



void
MatrixEpetra::zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context )
{
    FEELPP_ASSERT (this->isInitialized()).error( "MatrixEpetra<> not properly initialized" );

    const Epetra_Map& rowMap(_M_mat->RowMatrixRowMap());
    const Epetra_Map& colMap(_M_mat->RowMatrixColMap());


    std::vector<int>::const_iterator rit = rows.begin();
    std::vector<int>::const_iterator ren = rows.end();
    int start = rowMap.MinMyGID();
    int stop = rowMap.MaxMyGID()+1;

    Debug(10010) << "[MatrixEpetra::zeroRows] Number of rows to zero out (except diagonal) : " << rows.size() << "\n";
    Debug(10010) << "[MatrixEpetra::zeroRows] start : " << start << " stop  : " << stop << "\n";

    Debug(10010) << "[MatrixEpetra::zeroRows] keep diag ? " << on_context.test( ON_ELIMINATION_KEEP_DIAGONAL ) << "\n";
    Debug(10010) << "[MatrixEpetra::zeroRows] on symmetric ? " << on_context.test( ON_ELIMINATION_SYMMETRIC ) << "\n";

    Epetra_Vector Diagonal( rowMap );
    _M_mat->ExtractDiagonalCopy( Diagonal );


    for( size_type i = 0; rit != ren; ++rit, ++i )
        {
            int myRow = rowMap.LID( *rit );

            if (myRow >= 0)
                {
                    int    NumEntries;
                    double* Values;
                    int* Indices;

                    Debug(10010) << "extract row gid: " << *rit << "( " << i << ") lid: " << myRow  << "\n";
                    //When doing ExtractMyRowView, Indices contain the local indices
                    int ierr = _M_mat->ExtractMyRowView( myRow, NumEntries, Values, Indices);
                    Debug( 10010 ) << "ExtractGlobalRowView1  ierr = " << ierr << "\n";
                    FEELPP_ASSERT( ierr == 0 )( ierr )( myRow )( NumEntries ).warn ("error in ExtractGlobalRowView");
                    if ( on_context.test( ON_ELIMINATION_SYMMETRIC ) )
                        {
                            // we suppose here that the pattern is
                            // symmetric (not necessarily the matrix)

                            // do a loop on all indices, j, in Indices and
                            // zero out corresponding columns
                            // (Indices[j],myRow)
                            for( int j = 0; j < NumEntries; ++j )
                                {
                                    if ( Indices[j] == myRow )
                                        continue;

                                    int NumRowEntries;
                                    double *RowValues;
                                    int *RowIndices;
                                    Debug( 10010 ) << "[zeroRows] working with row lid: " << Indices[j]
                                                   << " (" << *rit << ", " << i << ")"
                                                   << " gid: " << rowMap.LID( Indices[j] )
                                                   << " is_local: " << rowMap.MyLID( Indices[j] ) << "\n";
                                    if ( rowMap.MyLID( Indices[j] ) )
                                        {
                                            int gid = Indices[j];
                                            Debug( 10010 ) << "[zeroRows] local, gid =" << gid << " lid = " << Indices[j] << "\n";
                                            ierr = _M_mat->ExtractMyRowView( Indices[j], NumRowEntries, RowValues, RowIndices);
                                            Debug( 10010 ) << "ExtractMyRowView ierr = " << ierr << "\n";
                                            FEELPP_ASSERT( ierr == 0 )( ierr )( Indices[j] )( NumRowEntries ).warn ("error in ExtractMyRowView/symm");
                                            bool found = false;
                                            for( int k = 0; k < NumRowEntries; ++k )
                                                if ( RowIndices[k] == myRow )
                                                    {
                                                        found=true;
                                                        rhs.add( gid, -values[i]*RowValues[k] );
                                                        RowValues[k] = 0.0;
                                                        break;
                                                    }
                                            FEELPP_ASSERT( found == true )( myRow )( Indices[j] ).error ( "matrix graph is not symmetric" );
                                        }
                                    else
                                        {
                                            FEELPP_ASSERT( 0 ).error ( "zeroRows() not available in parallel computation" );
                                        }
                                }
                        }

                    std::fill( Values, Values+NumEntries, 0.0 );

                    NumEntries = 1;
                    if ( on_context.test( ON_ELIMINATION_KEEP_DIAGONAL ) )
                        {
                            ierr = _M_mat->ReplaceMyValues(myRow, NumEntries, &Diagonal[myRow], &myRow  );
                            FEELPP_ASSERT( ierr == 0 )( ierr )( myRow )( NumEntries ).warn ("error in ReplaceMyValues()/diag");
                        }
                    else
                        {
                            /* put diagonal to 1 */
                            Diagonal[myRow] = 1.0;
                            ierr = _M_mat->ReplaceMyValues(myRow, NumEntries, &Diagonal[myRow], &myRow  );
                            FEELPP_ASSERT( ierr == 0 )( ierr )( myRow )( NumEntries ).warn ("error in ReplaceMyValues()/1");
                        }

                    M_bc_index.push_back(*rit);

                    // warning: a row index may belong to another
                    // processor, so make sure that we access only the
                    // rows that belong to this processor
                    rhs.set( *rit, values[i]*Diagonal[myRow] );

                }
        }

    if ( on_context.test( ON_ELIMINATION_SYMMETRIC ) )
        {

        }
}

void
MatrixEpetra::addMatrix ( int* rows, int nrows,
                          int* cols, int ncols,
                          value_type* data )
{
    FEELPP_ASSERT (this->isInitialized()).error( "MatrixEpetra<> not properly initialized" );

    _M_mat->SumIntoGlobalValues( nrows, rows, ncols, cols, data, Epetra_FECrsMatrix::ROW_MAJOR );
}


//
// Explicit instantiation
//

#endif // HAVE_TRILINOS_EPETRA
}
