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
   \file graphcsr.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-10-29
 */

#ifndef __GraphCSR_H
#define __GraphCSR_H 1

#include <vector>
#include <map>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/shared_ptr.hpp>

#include <feel/feelcore/feel.hpp>

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
              size_type last_col_entry_on_proc = 0 );

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

    /**
     * \return the first entry index on proc
     */
    size_type firstRowEntryOnProc() const
    {
        return M_first_row_entry_on_proc;
    }

    /**
     * \return the last entry index on proc
     */
    size_type lastRowEntryOnProc() const
    {
        return M_last_row_entry_on_proc;
    }
    /**
     * \return the first entry index on proc
     */
    size_type firstColEntryOnProc() const
    {
        return M_first_col_entry_on_proc;
    }

    /**
     * \return the last entry index on proc
     */
    size_type lastColEntryOnProc() const
    {
        return M_last_col_entry_on_proc;
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
    iterator begin() { return M_storage.begin(); }

    /**
     * \return end (rw) iterator on graph
     */
    iterator end() { return M_storage.end(); }

    /**
     * \return begin (ro) iterator on graph
     */
    const_iterator begin() const { return M_storage.begin(); }
    /**
     * \return end (ro) iterator on graph
     */
    const_iterator end() const { return M_storage.end(); }
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
     * \return the maximum number of non-zero entries per row
     */
    size_type maxNnz() const { return M_max_nnz; }

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
    mpi::communicator const& comm() const { return M_comm; }

    nz_type const& ia() const { return M_ia; }
    nz_type const& ja() const { return M_ja; }
    std::vector<double> const& a() const { return M_a; }


    //@}

    /** @name  Mutators
     */
    //@{


    void setFirstRowEntryOnProc( size_type entry ) { M_first_row_entry_on_proc = entry; }
    void setFirstColEntryOnProc( size_type entry ) { M_first_col_entry_on_proc = entry; }

    void setLastRowEntryOnProc( size_type entry ) { M_last_row_entry_on_proc = entry; }
    void setLastColEntryOnProc( size_type entry ) { M_last_col_entry_on_proc = entry; }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * close the graph: compute some information per row (ie number of
     * non-zero entries per row )
     */
    void close();

    //@}



protected:

private:

    mpi::communicator M_comm;

    size_type M_first_row_entry_on_proc;
    size_type M_last_row_entry_on_proc;
    size_type M_first_col_entry_on_proc;
    size_type M_last_col_entry_on_proc;
    size_type M_max_nnz;
    nz_type M_n_total_nz;
    nz_type M_n_nz;
    nz_type M_n_oz;
    storage_type M_storage;
    nz_type M_ia, M_ja;
    std::vector<double> M_a;
};

} // Feel
#endif /* __GraphCSR_H */
