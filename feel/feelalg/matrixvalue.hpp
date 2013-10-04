/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-13

  Copyright (C) 2005,2006 EPFL

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
   \file matrixgmm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-13
 */
#ifndef __MatrixValue_H
#define __MatrixValue_H 1

#include <set>

#include <boost/numeric/ublas/vector.hpp>


namespace Feel
{
/*!
 * \class MatrixValue
 * \brief interface to matrix
 *
 * \code
 * MatrixValue<T> m;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T>
class MatrixValue
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;

    typedef value_type matrix_type;
    typedef std::vector<std::set<size_type> > pattern_type;

    static const bool is_row_major = true;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixValue( value_type acc = value_type( 0 ) )
        :
        M_mat( acc )
    {}
    MatrixValue( MatrixValue const & m )
        :
        M_mat( m.M_mat )
    {}

    ~MatrixValue()
    {
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    unsigned int size1 () const
    {
        return 1;
    }

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    unsigned int size2 () const
    {
        return 1;
    }

    /**
     * \return the number of non-zeros entries in the matrix
     */
    size_type nnz() const
    {
        return 1;
    }

    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    unsigned int rowStart () const
    {
        return 0;
    }

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    unsigned int rowStop () const
    {
        return 0;
    }

    /**
     * \return true if matrix is initialized/usable, false otherwise
     */
    bool isInitialized() const
    {
        return true;
    }

    /**
     * \c close the gmm matrix, that will copy the content of write
     * optimized matrix into a read optimized matrix
     */
    void close () const;


    /**
     * see if Gmm matrix has been closed
     * and fully assembled yet
     */
    bool closed() const
    {
        return true;
    }


    /**
     * Returns the read optimized gmm matrix.
     */
    matrix_type const& mat () const
    {
        return M_mat;
    }

    /**
     * Returns the read optimized gmm matrix.
     */
    matrix_type & mat ()
    {
        return M_mat;
    }

    /**
     * Returns the write optimized gmm matrix.
     */
    matrix_type const& wmat () const
    {
        return M_mat;
    }

    /**
     * Returns the write optimized gmm matrix.
     */
    matrix_type & wmat ()
    {
        return M_mat;
    }



    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Initialize a Value matrix that is of global
     * dimension \f$ m \times  n \f$ with local dimensions
     * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
     * nonzeros per row (defaults to 30).
     * \p noz is the number of on-processor
     * nonzeros per row (defaults to 30).
     */
    void init ( const unsigned int /*m*/,
                const unsigned int /*n*/,
                const unsigned int /*m_l*/,
                const unsigned int /*n_l*/,
                const unsigned int /*nnz*/=30,
                const unsigned int /*noz*/=10 )
    {
        this->zero();
    }

    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear ()
    {
        M_mat = 0;
    }

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ()
    {
        M_mat = 0;
    }

    void zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
    {
        M_mat = 0;
    }

    /**
     * Add \p value to the value already accumulated
     */
    void add ( const unsigned int /*i*/,
               const unsigned int /*j*/,
               const value_type value )
    {
        M_mat += value;
    }

    /**
     * set to \p value
     */
    void set ( const unsigned int /*i*/,
               const unsigned int /*j*/,
               const value_type value )
    {
        M_mat = value;
    }



    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL" ) const;


    /**
     * fill sparse matrix with non zero entries
     */
    void fill( pattern_type const& );

    void resize( size_type /* nr*/, size_type /*nc*/, bool /*preserve*/ = false )
    {

    }

    /**
     * In this case the energy is the value that has been accumulated
     */
    value_type
    energy( ublas::vector<value_type> const& /*__v*/,
            ublas::vector<value_type> const& /*__u*/ ) const
    {
        return M_mat;
    }

    /**
     *
     */
    void diagonalize( size_type );
    //@}



protected:

private:

    /**
     * the gmm sparse matrix data structure
     */
    mutable matrix_type M_mat;

};

template<typename T>
void
MatrixValue<T>::diagonalize( size_type __dof_index )
{
    FEELPP_ASSERT( 0 ).error( "diagonalize is undefined for this matrix type" );
}
template<typename T>
void
MatrixValue<T>::fill( pattern_type const& /*__pattern*/ )
{
}

template<typename T>
void
MatrixValue<T>::close() const
{
}

template<typename T>
void
MatrixValue<T>::printMatlab( const std::string /*filename*/ ) const
{
}

} // Feel
#endif /* __MatrixValue_H */
