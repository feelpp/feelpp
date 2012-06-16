/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-13
 */
#ifndef __MatrixEigen_H
#define __MatrixEigen_H 1

#include <set>

#include <boost/version.hpp>
#if (BOOST_VERSION >= 103400)
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <Eigen/Core>



#include <eigen_matrix.h>
#include <eigen_sub_matrix.h>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vectorublas.hpp>


namespace eigen
{
namespace Feel
{
template<typename T, typename Storage> class VectorUblas;

/*!
 * \class MatrixEigen
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
class MatrixEigen : public MatrixSparse<T>
{
    typedef MatrixSparse<T> super;
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename type_traits<value_type>::real_type real_type;

    typedef typename mpl::if_<boost::is_same<LayoutType, eigen::row_major>,
            mpl::identity<eigen::csr_matrix<value_type> >,
            typename mpl::if_<boost::is_same<LayoutType, eigen::col_major>,
            mpl::identity<eigen::csc_matrix<value_type> >,
            mpl::identity<boost::none_t> >::type>::type::type matrix_type;


    static const bool is_row_major = boost::is_same<LayoutType,eigen::row_major>::value;


    typedef std::vector<std::set<size_type> > pattern_type;

    typedef eigen::row_matrix<eigen::wsvector<value_type> > write_matrix_type;

    typedef typename super::graph_type graph_type;
    typedef typename super::graph_ptrtype graph_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixEigen();

    MatrixEigen( size_type r, size_type c );

    MatrixEigen( MatrixEigen const & m );

    ~MatrixEigen();


    //@}

    /** @name Operator overloads
     */
    //@{

    MatrixEigen<T,LayoutType> & operator = ( MatrixSparse<value_type> const& M )
    {
        return *this;
    }


    value_type  operator()( size_type i, size_type j ) const
    {
        return _M_mat( i, j );
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
        return _M_mat.nrows();
    }

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size2 () const
    {
        return _M_mat.ncols();
    }

    /**
     * \return the number of non-zeros entries in the matrix
     */
    size_type nnz() const
    {
        return eigen::nnz( _M_mat );
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
        return _M_is_initialized;
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
        return _M_is_closed;
    }


    /**
     * Returns the read optimized eigen matrix.
     */
    matrix_type const& mat () const
    {
        return _M_mat;
    }

    /**
     * Returns the read optimized eigen matrix.
     */
    matrix_type & mat ()
    {
        return _M_mat;
    }

    /**
     * Returns the write optimized eigen matrix.
     */
    write_matrix_type const& wmat () const
    {
        return _M_wmat;
    }

    /**
     * Returns the write optimized eigen matrix.
     */
    write_matrix_type & wmat ()
    {
        return _M_wmat;
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
        //eigen::resize( _M_mat, 0, 0 );
        eigen::resize( _M_wmat, 0, 0 );
    }

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ()
    {
        eigen::clear( _M_wmat );
    }

    void zero ( size_type start1, size_type stop1, size_type start2, size_type stop2 )
    {
        eigen::clear( eigen::sub_matrix( _M_wmat,
                                     eigen::sub_interval( start1, stop1-start1 ),
                                     eigen::sub_interval( start2, stop2-start2 ) ) );
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
        _M_wmat( i, j ) += value;
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
        _M_wmat( i, j ) = value;
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

    void resize( size_type nr, size_type nc, bool /*preserve*/ = false );

    /**
     * \return \f$ v^T M u \f$
     */
    value_type
    energy( Vector<value_type> const& __v,
            Vector<value_type> const& __u, bool transpose = false ) const;

    /**
     * eliminates row without change pattern, and put 1 on the diagonal
     * entry
     *
     *\warning if the matrix was symmetric before this operation, it
     * won't be afterwards. So use the proper solver (nonsymmetric)
     */
    void zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context );

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
                     value_type* data )
    {
        // NOT IMPLEMENTED YET (eigen support should get dropped in fact)
    }

    void scale( const T a ) {}

    /**
     * Returns the transpose of a matrix
     *
     * \param M the matrix to transpose
     * \param Mt the matrix transposed
     */
    void transpose( MatrixSparse<value_type>& Mt ) const;

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
    void updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m, size_type start_i, size_type start_j )
    {
#warning todo!
    }

    //@}



protected:

private:

    bool _M_is_initialized;
    mutable bool _M_is_closed;

    /**
     * the eigen sparse matrix data structure
     */
    mutable matrix_type _M_mat;

    /**
     * write optimized matrix
     */
    mutable write_matrix_type _M_wmat;
};


template<typename T, typename LayoutType>
void
MatrixEigen<T, LayoutType>::zeroRows( std::vector<int> const& rows,
                                    std::vector<value_type> const& vals,
                                    Vector<value_type>& rhs,
                                    Context const& on_context )
{
    Feel::detail::ignore_unused_variable_warning( rhs );
    Feel::detail::ignore_unused_variable_warning( vals );

    eigen::resize( _M_wmat, eigen::mat_nrows( _M_mat ), eigen::mat_ncols( _M_mat ) );
    eigen::copy( _M_mat, _M_wmat );

    for ( size_type i = 0; i < rows.size(); ++i )
    {
        value_type value = 1.0;

        if ( on_context.test( ON_ELIMINATION_KEEP_DIAGONAL ) )
            value = _M_wmat.row( rows[i] ).r( rows[i] );

        eigen::clear( eigen::mat_row( _M_wmat, rows[i] ) );

        // set diagonal
        _M_wmat.row( rows[i] ).w( rows[i], value );

        // multiply rhs by value of the diagonal entry value
        rhs.set( rows[i], value * vals[i] );
    }

    eigen::copy( _M_wmat, _M_mat );

}

} // Feel
#endif /* __MatrixEigen_H */
