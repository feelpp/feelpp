/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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

    GraphCSR( boost::shared_ptr<DataMap> const& mapRow,
              boost::shared_ptr<DataMap> const& mapCol );

    GraphCSR( DataMap const& mapRow,
              DataMap const& mapCol );

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

    DataMap const& mapRow() const { return *M_mapRow; }
    DataMap const& mapCol() const { return *M_mapCol; }
    boost::shared_ptr<DataMap> const& mapRowPtr() const { return M_mapRow; }
    boost::shared_ptr<DataMap> const& mapColPtr() const { return M_mapCol; }
    boost::shared_ptr<DataMap> mapRowPtr() { return M_mapRow; }
    boost::shared_ptr<DataMap> mapColPtr() { return M_mapCol; }

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

    boost::shared_ptr<DataMap> M_mapRow, M_mapCol;
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

private :
    bool M_isClosed;


};



} // Feel
#endif /* __GraphCSR_H */
