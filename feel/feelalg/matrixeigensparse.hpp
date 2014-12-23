/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-13

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \file matrixeigen.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-13
 */
#ifndef __MatrixEigenSparse_H
#define __MatrixEigenSparse_H 1

#include <set>

#include <boost/version.hpp>
#if (BOOST_VERSION >= 103400)
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vectorublas.hpp>


namespace Feel
{
template<typename T, typename Storage> class VectorUblas;

/*!
 *
 * \brief interface to eigen sparse matrix
 *
 * this class is a wrapper around \c csr_matrix<> and \c csc_matrix<>
 * data type from \c eigen:: .
 *
 *
 * \code
 * // csr matrix
 * MatrixEigen<T,eigen::row_major> m;
 * // csc matrix
 * MatrixEigen<T,eigen::col_major> m;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T>
class MatrixEigenSparse : public MatrixSparse<T>
{
    typedef MatrixSparse<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename type_traits<value_type>::real_type real_type;
    typedef Eigen::SparseMatrix<T> matrix_type;
    typedef typename super::graph_type graph_type;
    typedef typename super::graph_ptrtype graph_ptrtype;
    typedef Eigen::Triplet<T> triplet;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixEigenSparse();

    MatrixEigenSparse( size_type r, size_type c, WorldComm const& worldComm=Environment::worldComm() );

    MatrixEigenSparse( MatrixEigenSparse const & m );

    ~MatrixEigenSparse();


    //@}

    /** @name Operator overloads
     */
    //@{

    MatrixEigenSparse<T> & operator = ( MatrixSparse<value_type> const& M )
    {
        return *this;
    }


    value_type  operator()( size_type i, size_type j ) const
    {
        //return M_mat.row(i).col(j);
        return 0.;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size1 () const
    {
        return M_mat.rows();
    }

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size2 () const
    {
        return M_mat.cols();
    }

    /**
     * \return the number of non-zeros entries in the matrix
     */
    size_type nnz() const
    {
        return M_mat.rows()*M_mat.cols();
    }

    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    size_type rowStart () const
    {
        return 0;
    }

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    size_type rowStop () const
    {
        return 0;
    }

    /**
     * \return true if matrix is initialized/usable, false otherwise
     */
    bool isInitialized() const
    {
        return M_is_initialized;
    }

    /**
     * \c close the eigen matrix, that will copy the content of write
     * optimized matrix into a read optimized matrix
     */
    void close () const;


    /**
     * see if Eigen matrix has been closed
     * and fully assembled yet
     */
    bool closed() const
    {
        return M_is_closed;
    }


    /**
     * Returns the read optimized eigen matrix.
     */
    matrix_type const& mat () const
    {
        return M_mat;
    }

    /**
     * Returns the read optimized eigen matrix.
     */
    matrix_type & mat ()
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
     * Initialize a Eigen matrix that is of global
     * dimension \f$ m \times  n \f$ with local dimensions
     * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
     * nonzeros per row (defaults to 30).
     * \p noz is the number of on-processor
     * nonzeros per row (defaults to 30).
     */
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                const size_type nnz=30,
                const size_type noz=10 );

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                graph_ptrtype const& graph );

    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear ()
    {
        //eigen::resize( M_mat, 0, 0 );
        M_mat.setZero();
    }

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ()
    {
        M_mat.setZero();
    }

    void zero ( size_type start1, size_type stop1, size_type start2, size_type stop2 )
    {
    }

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    void add ( const size_type i,
               const size_type j,
               const value_type& value )
    {
        M_tripletList.push_back(triplet(i, j, value) );
    }

    /**
      * set \p value to the element
      * \p (i,j).  Throws an error if
      * the entry does not
      * exist. Still, it is allowed to
      * store zero values in
      * non-existent fields.
      */
    void set ( const size_type i,
               const size_type j,
               const value_type& value )
    {
        M_tripletList.push_back(triplet(i, j, value));
    }


    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL" ) const;



    void resize( size_type nr, size_type nc, bool /*preserve*/ = false );

    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    void diagonal ( Vector<T>& dest ) const;

    /**
     * \return \f$ v^T M u \f$
     */
    real_type
    energy( Vector<value_type> const& __v,
            Vector<value_type> const& __u, bool transpose = false ) const;

    /**
     * eliminates row without change pattern, and put 1 on the diagonal
     * entry
     *
     *\warning if the matrix was symmetric before this operation, it
     * won't be afterwards. So use the proper solver (nonsymmetric)
     */
    void zeroRows( std::vector<int> const& rows, Vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context );

    void init() {}

    /**
     * Add the full matrix to the
     * Petsc matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix( const ublas::matrix<T, ublas::row_major>&,
                    const std::vector<size_type>&,
                    const std::vector<size_type>& ) {}

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    void addMatrix( const boost::numeric::ublas::matrix<T, ublas::row_major>&, const std::vector<size_type>& ) {}

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    void addMatrix( value_type v, MatrixSparse<value_type>& _m );

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix ( int* rows, int nrows,
                     int* cols, int ncols,
                     value_type* data );

    void scale( const T a );

    /**
     * Returns the transpose of a matrix
     *
     * \param M the matrix to transpose
     * \param Mt the matrix transposed
     */
    void transpose( MatrixSparse<value_type>& Mt, size_type options ) const;

    /**
     * Return the l1-norm of the matrix, that is
     * \f$|M|_1=max_{all columns j}\sum_{all rows i} |M_ij|\f$, (max. sum of columns).
     *
     * This is the natural matrix norm that is compatible to the
     * l1-norm for vectors, i.e.  \f$|Mv|_1\leq |M|_1 |v|_1\f$.
     */
    real_type l1Norm() const
    {
        return real_type( 0 );
    }

    /**
     * Return the linfty-norm of the matrix, that is
     *
     * \f$|M|_\infty=max_{all rows i}\sum_{all columns j} |M_ij|\f$,
     *
     * (max. sum of rows).
     * This is the natural matrix norm that is
     * compatible to the linfty-norm of vectors, i.e.
     * \f$|Mv|_\infty \leq |M|_\infty |v|_\infty\f$.
     */
    real_type linftyNorm() const
    {
        return real_type( 0 );
    }

    /**
     * update a block matrix
     */
    void updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m, std::vector<size_type> start_i, std::vector<size_type> start_j );

    //@}



protected:

private:

    bool M_is_initialized;
    mutable bool M_is_closed;

    /**
     * the eigen sparse matrix data structure
     */
    mutable matrix_type M_mat;
    mutable std::vector<triplet> M_tripletList;
};

#if !defined( FEELPP_INSTANTIATE_MATRIXEIGENSPARSE )
extern template class MatrixEigenSparse<double>;
extern template class MatrixEigenSparse<std::complex<double>>;
#endif


} // Feel
#endif /* __MatrixEigenSparse_H */
