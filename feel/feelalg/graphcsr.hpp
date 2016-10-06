/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \file graphcsr.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-10-29
 */

#ifndef __GraphCSR_H
#define __GraphCSR_H 1

#include <vector>
#include <map>
#include <set>
#include <boost/tuple/tuple.hpp>
//#include <boost/mpi/communicator.hpp>

#include <boost/shared_ptr.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/datamap.hpp>
#include <feel/feelvf/block.hpp>

namespace Feel
{
namespace mpi=boost::mpi;

/**
 * \class GraphCSR
 * \brief Graph representation of the Compressed Sparse Row format
 *
 * @author Christophe Prud'homme
 * @see
 */
class GraphCSR
{
public:


    /** @name Typedefs
     */
    //@{

    typedef GraphCSR self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef std::vector<size_type> nz_type;
    typedef boost::shared_ptr<nz_type> nz_ptrtype;

    //typedef boost::tuple<size_type, size_type, std::vector<size_type> > row_type;
    typedef boost::tuple<size_type, size_type, std::set<size_type> > row_type;
    typedef std::map<size_type, row_type > storage_type;
    typedef boost::shared_ptr<storage_type> storage_ptrtype;

    typedef storage_type::iterator iterator;
    typedef storage_type::const_iterator const_iterator;

    typedef DataMap datamap_type;
    typedef boost::shared_ptr<DataMap> datamap_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     * \param n number of rows in the graph
     */
    GraphCSR( size_type n = 0,
              size_type first_row_entry_on_proc = 0,
              size_type last_row_entry_on_proc = 0,
              size_type first_col_entry_on_proc = 0,
              size_type last_col_entry_on_proc = 0,
              WorldComm const& worldcomm = Environment::worldComm() );

    GraphCSR( datamap_ptrtype const& mapRow,
              datamap_ptrtype const& mapCol );

    GraphCSR( datamap_type const& mapRow,
              datamap_type const& mapCol );

    GraphCSR( vf::BlocksBase<self_ptrtype> const & blockSet,
              bool diagIsNonZero=true,
              bool close=true );

    /**
     * copy constructor
     */
    GraphCSR( GraphCSR const & g );

    /**
     * destructor
     */
    ~GraphCSR();

    //@}

    /** @name Operator overloads
     */
    //@{

    /** copy operator */
    GraphCSR& operator=( GraphCSR const& g );

    //@}

    /** @name Accessors
     */
    //@{

    size_type nRows() const
    {
        return M_last_row_entry_on_proc[this->worldComm().globalRank()]+1;
    }
    size_type nCols() const
    {
        return M_last_col_entry_on_proc[this->worldComm().globalRank()]+1;
    }

    /**
     * \return the first entry index on proc
     */
    size_type firstRowEntryOnProc() const
    {
        return M_first_row_entry_on_proc[this->worldComm().globalRank()];
    }

    /**
     * \return the last entry index on proc
     */
    size_type lastRowEntryOnProc() const
    {
        return M_last_row_entry_on_proc[this->worldComm().globalRank()];
    }
    /**
     * \return the first entry index on proc
     */
    size_type firstColEntryOnProc() const
    {
        return M_first_col_entry_on_proc[this->worldComm().globalRank()];
    }

    /**
     * \return the last entry index on proc
     */
    size_type lastColEntryOnProc() const
    {
        return M_last_col_entry_on_proc[this->worldComm().globalRank()];
    }

    /**
     * \return the number of rows in the pattern
     */
    size_type size() const
    {
        return M_storage.size();
    }
    /**
     * \return true if the graph is empty, false otherwise
     */
    bool empty() const
    {
        return M_storage.empty();
    }
    /**
     * \return begin (rw) iterator on graph
     */
    iterator begin()
    {
        return M_storage.begin();
    }

    /**
     * \return end (rw) iterator on graph
     */
    iterator end()
    {
        return M_storage.end();
    }

    /**
     * \return begin (ro) iterator on graph
     */
    const_iterator begin() const
    {
        return M_storage.begin();
    }
    /**
     * \return end (ro) iterator on graph
     */
    const_iterator end() const
    {
        return M_storage.end();
    }
    /**
     * get the i-th row
     */
    row_type& row( size_type i )
    {
        return M_storage[i];
    }

    /**
     * get the i-th row (const)
     */
    row_type const& row( size_type i ) const
    {
        return M_storage.find( i )->second;
    }

    /**
     * return storage
     */
    storage_type const& storage() const
    {
        return M_storage;
    }


    /**
     * \return the maximum number of non-zero entries per row
     */
    size_type maxNnz() const
    {
        return M_max_nnz;
    }

    /**
     * \return the array containing the number of non-zero entries per
     * row that the current processor will deal with (rows that belong
     * to the proc or not)
     */
    std::vector<size_type> const&
    nNz() const
    {
        return M_n_total_nz;
    }

    /**
     * \return the array containing the number of non-zero entries per
     * row on the current processor
     */
    std::vector<size_type> const&
    nNzOnProc() const
    {
        return M_n_nz;
    }

    /**
     * \return the array containing the number of non-zero entries per
     * row on other processors
     */
    std::vector<size_type> const&
    nNzOffProc() const
    {
        return M_n_oz;
    }

    /**
     * \return the communicator
     */
    //mpi::communicator const& comm() const { return M_comm; }
    WorldComm const& worldComm() const
    {
        return M_worldComm;
    }

    nz_type const& ia() const
    {
        return M_ia;
    }
    nz_type const& ja() const
    {
        return M_ja;
    }
    std::vector<double> const& a() const
    {
        return M_a;
    }
    std::vector<double>& a()
    {
        return M_a;
    }


    //@}

    /** @name  Mutators
     */
    //@{


    void setFirstRowEntryOnProc( size_type entry )
    {
        M_first_row_entry_on_proc[this->worldComm().globalRank()] = entry;
    }
    void setFirstColEntryOnProc( size_type entry )
    {
        M_first_col_entry_on_proc[this->worldComm().globalRank()] = entry;
    }

    void setLastRowEntryOnProc( size_type entry )
    {
        M_last_row_entry_on_proc[this->worldComm().globalRank()] = entry;
    }
    void setLastColEntryOnProc( size_type entry )
    {
        M_last_col_entry_on_proc[this->worldComm().globalRank()] = entry;
    }

    datamap_type const& mapRow() const { return *M_mapRow; }
    datamap_type const& mapCol() const { return *M_mapCol; }
    datamap_ptrtype const& mapRowPtr() const { return M_mapRow; }
    datamap_ptrtype const& mapColPtr() const { return M_mapCol; }
    datamap_ptrtype mapRowPtr() { return M_mapRow; }
    datamap_ptrtype mapColPtr() { return M_mapCol; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * graph has not entries
     */
    void zero();

    /**
     * close the graph: compute some information per row (ie number of
     * non-zero entries per row )
     */
    void close();

    /**
     * transpose graph
     */
    self_ptrtype transpose( bool doClose = true );

    /**
     * create subgraph from a list of global process index
     */
    self_ptrtype createSubGraph( std::vector<size_type> const& rows, std::vector<size_type> const& cols,
                                 datamap_ptrtype const& submapRow = datamap_ptrtype(),
                                 datamap_ptrtype const& submapCol = datamap_ptrtype(),
                                 bool useSameDataMap=false, bool checkAndFixRange=true ) const;

    /**
     * add missing zero entries diagonal
     */
    void addMissingZeroEntriesDiagonal();

    /**
     * showMe
     */
    void showMe( std::ostream& __out = std::cout ) const;

    void printPython( std::string const& nameFile ) const;

    //@}

private :

    void mergeBlockGraph( self_ptrtype const& g,
                          size_type start_i, size_type start_j );

    void mergeBlockGraphMPI( self_ptrtype const& g, vf::BlocksBase<self_ptrtype> const & blockSet, int i, int j,
                             size_type start_i, size_type start_j );

    void updateDataMap( vf::BlocksBase<self_ptrtype> const & blockSet );

    size_type nLocalDofWithoutGhostOnProcStartRow( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const;
    size_type nLocalDofWithoutGhostOnProcStartCol( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const;
    size_type nLocalDofWithGhostOnProcStartRow( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const;
    size_type nLocalDofWithGhostOnProcStartCol( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const;

protected:

private:
    bool M_is_closed;
    //mpi::communicator M_comm;
    WorldComm M_worldComm;

    std::vector<size_type> M_first_row_entry_on_proc;
    std::vector<size_type> M_last_row_entry_on_proc;
    std::vector<size_type> M_first_col_entry_on_proc;
    std::vector<size_type> M_last_col_entry_on_proc;
    size_type M_max_nnz;
    nz_type M_n_total_nz;
    nz_type M_n_nz;
    nz_type M_n_oz;
    storage_type M_storage;
    nz_type M_ia, M_ja;
    std::vector<double> M_a;

    self_ptrtype M_graphT;

    datamap_ptrtype M_mapRow, M_mapCol;
};


class BlocksBaseGraphCSR : public vf::BlocksBase<boost::shared_ptr<GraphCSR> >
{
public :
    typedef vf::BlocksBase<boost::shared_ptr<GraphCSR> > super_type;
    typedef super_type::index_type index_type;
    typedef BlocksBaseGraphCSR self_type;
    typedef boost::shared_ptr<GraphCSR> graph_ptrtype;

    BlocksBaseGraphCSR( index_type nr=0,index_type nc=0 )
        :
        super_type( nr,nc ),
        M_isClosed( false )
    {}

    BlocksBaseGraphCSR( self_type const & b )
        :
        super_type( b ),
        M_isClosed( b.M_isClosed )
    {}

    BlocksBaseGraphCSR( super_type const & b )
        :
        super_type( b ),
        M_isClosed( false )
    {}

    self_type
    operator<<( graph_ptrtype const& g ) const
    {
        return super_type::operator<<( g );
    }

    void close();

    bool isClosed() const { return M_isClosed; }

    /**
     * termination function for fill
     */
    template<typename Arg1>
    void  fill( int& n, int& c, const Arg1& arg1 ) {}

    /**
     * fill row and col
     */
    template<typename Arg1, typename Arg2, typename ...Args>
    void
    fill( int& n, int& c, const Arg1& arg1, const Arg2& arg2, const Args&... args )
        {
            this->operator()( n,c ) = stencil( _test=arg1,_trial=arg2, _diag_is_nonzero=false, _close=false)->graph();
            if ( c > n )
                this->operator()( c,n ) = stencil( _test=arg2,_trial=arg1, _diag_is_nonzero=false, _close=false)->graph();
            DCHECK( this->operator()( n,c ) != nullptr ) << "Invalid stencil " << n << "," << c;
            DCHECK( this->operator()( c,n ) != nullptr ) << "Invalid stencil " << c << "," << n;
            cout << "stencil for (" << n << "," << c << ") and (" << c << "," << n << ")\n";
            // go to next submatrix
            fill( n, ++c, arg1, args... );
        }


    /**
     * termination function for rowcol recusion
     */
    void rowcol( int r  ) {}

    /**
     * recurse through first row and column of submatrix of size (n,n)
     */
    template<typename Arg1, typename ...Args>
    void
    rowcol( int n, const Arg1& arg1, const Args&... args )
        {
            int c = n;
            // fill row and column (n,n) with csr graph
            fill( n, c, arg1, arg1, args... );
            // do submatrix (n+1,n+1)
            rowcol( ++n, args...);
        }

private :
    bool M_isClosed;


};

/**
 * Build blocks of CSR graphs with function spaces \p
 * (args1,args2,argn). Variadic template is used to handle an arbitrary number of
 * function spaces.
 * The blocks are organized then matrix wise with the stencil associated of pairs of function spaces in \p (arg1,...,argn)
 *
 */
template<typename Arg1, typename ...Args>
BlocksBaseGraphCSR
csrGraphBlocks( const Arg1& arg1, const Args&... args )
{
    const int size = sizeof...(Args)+1;
    BlocksBaseGraphCSR g( size, size );
    int n = 0;
    g.rowcol( n, arg1, args...);
    return g;
}


} // Feel
#endif /* __GraphCSR_H */
