/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2008-01-03

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file matrixblock.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-06-10
 */

#ifndef __MatrixBlock_H
#define __MatrixBlock_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelvf/block.hpp>

namespace Feel
{

template <typename T,typename SizeT>
class Backend;

template <typename T = double>
class BlocksBaseSparseMatrix : public vf::BlocksBase<std::shared_ptr<MatrixSparse<T>>>
{
  public:
    typedef vf::BlocksBase<std::shared_ptr<MatrixSparse<T>>> super_type;
    typedef typename super_type::index_type index_type;
    typedef BlocksBaseSparseMatrix<T> self_type;
    typedef std::shared_ptr<MatrixSparse<T>> matrix_sparse_ptrtype;

    BlocksBaseSparseMatrix( index_type nr = 0, index_type nc = 0 )
        : super_type( nr, nc ),
          M_isClosed( false )
    {
    }

    BlocksBaseSparseMatrix( self_type const& b )
        : super_type( b ),
          M_isClosed( b.M_isClosed )
    {
    }

    BlocksBaseSparseMatrix( super_type const& b )
        : super_type( b ),
          M_isClosed( false )
    {
    }

    self_type
    operator<<( matrix_sparse_ptrtype const& m ) const
    {
        return super_type::operator<<( m );
    }

    void close();

    bool isClosed() const { return M_isClosed; }

  private:
    bool M_isClosed;
};

template <int NR, int NC, typename T = double>
class BlocksSparseMatrix : public BlocksBaseSparseMatrix<T>
{
  public:
    typedef typename BlocksBaseSparseMatrix<T>::index_type index_type;
    static const index_type NBLOCKROWS = NR;
    static const index_type NBLOCKCOLS = NC;

    typedef BlocksBaseSparseMatrix<T> super_type;

    BlocksSparseMatrix()
        : super_type( NBLOCKROWS, NBLOCKCOLS )
    {
    }
};

/**
 * \class MatrixBlock
 * \brief block of matrices
 *
 * <code>
 * auto myBlocks = Blocks<2,2,double>()<< A11 << A12
 *                                     << A21 << A22;
 *
 * auto A = backend->newBlockMatrix(myBlocks);
 * </code>
 *
 * @author Vincent Chabannes
 */

template <typename T>
class MatrixBlockBase : public MatrixSparse<T>
{
    typedef MatrixSparse<T> super;

  public:
    /** @name Typedefs
     */
    //@{
    using size_type = typename super::size_type;
    typedef MatrixBlockBase<T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;

    typedef Backend<value_type,size_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    typedef super matrix_type;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;

    typedef std::vector<matrix_ptrtype> vector_matrix_ptrtype;

    typedef typename super::graph_type graph_type;
    typedef typename super::graph_ptrtype graph_ptrtype;

    typedef typename super::indexsplit_type indexsplit_type;
    typedef typename super::indexsplit_ptrtype indexsplit_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixBlockBase( vf::BlocksBase<matrix_ptrtype> const& blockSet,
                     backend_ptrtype backend,
                     bool copy_values = true,
                     bool diag_is_nonzero = true );

    MatrixBlockBase( vf::BlocksBase<graph_ptrtype> const& graph,
                     backend_ptrtype backend,
                     bool diag_is_nonzero = true );

    MatrixBlockBase( vf::BlocksBase<matrix_ptrtype> const& blockSet,
                     backend_type& backend,
                     bool copy_values = true,
                     bool diag_is_nonzero = true )
        : MatrixBlockBase( blockSet, backend.shared_from_this(), copy_values, diag_is_nonzero )
    {
    }

    MatrixBlockBase( vf::BlocksBase<graph_ptrtype> const& graph,
                     backend_type& backend,
                     bool diag_is_nonzero = true )
        : MatrixBlockBase( graph, backend.shared_from_this(), diag_is_nonzero )
    {
    }

    explicit MatrixBlockBase( worldcomm_ptr_t const& wc )
        : super( wc ),
          M_backend(),
          M_mat()
        {}
    
    MatrixBlockBase( MatrixBlockBase const& mb )
        : super( mb ),
          M_backend( mb.M_backend ),
          M_mat( mb.M_mat )
    {
    }

    ~MatrixBlockBase() override
    {
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    MatrixBlockBase& operator=( MatrixBlockBase const& mb )
    {
        if ( this != &mb )
        {
            M_mat = mb.M_mat;
        }

        return *this;
    }
    using clone_ptrtype = typename super::clone_ptrtype;
    clone_ptrtype clone() const override { return std::make_shared<MatrixBlockBase<T>>( *this ); }
    //@}

    /** @name Accessors
     */
    //@{
    matrix_ptrtype getSparseMatrix()
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
     * Initialize a Petsc matrix that is of global
     * dimension \f$ m \times  n \f$ with local dimensions
     * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
     * nonzeros per row (defaults to 30).
     * \p noz is the number of on-processor
     * nonzeros per row (defaults to 30).
     */
    void init( const size_type m,
               const size_type n,
               const size_type m_l,
               const size_type n_l,
               const size_type nnz = 30,
               const size_type noz = 10 ) override;

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    void init( const size_type m,
               const size_type n,
               const size_type m_l,
               const size_type n_l,
               graph_ptrtype const& graph ) override;

    //!
    //! get the backend read-only
    //!
    backend_ptrtype const& backend() const { return M_backend; }

    //!
    //! get the backend read-write
    //!
    backend_ptrtype backend()  { return M_backend; }

    //!
    //! set the backend
    //!
    void setBackend( backend_ptrtype b ) { M_backend = b; }
    
    /**
     * Release all memory and return to a state just like after having
     * called the default constructor.
     */
    void clear() override;

    /**
     * Set all entries to 0.
     */
    void zero() override;

    /**
     * Set entries between to 0.
     */
    void zero( size_type start1, size_type size1,
               size_type start2, size_type size2 ) override;

    /**
     * Call the Sparse assemble routines.  sends necessary messages to
     * other processors
     */
    void close() const override;

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size1() const override;

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size2() const override;

    //!
    //! @return the number of non-zero entries
    //!
    size_type nnz() const override { return M_mat->nnz(); }
    
    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    size_type rowStart() const override;

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    size_type rowStop() const override;

    /**
     * Set the element \p (i,j) to \p value.
     * Throws an error if the entry does
     * not exist. Still, it is allowed to store
     * zero values in non-existent fields.
     */
    void set( const size_type i,
              const size_type j,
              const value_type& value ) override;

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    void add( const size_type i,
              const size_type j,
              const value_type& value ) override;

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix( const ublas::matrix<value_type>& dm,
                    const std::vector<size_type>& rows,
                    const std::vector<size_type>& cols ) override;

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix( int* rows, int nrows,
                    int* cols, int ncols,
                    value_type* data,
                    size_type K = 0,
                    size_type K2 = invalid_v<size_type> ) override;

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    void addMatrix( const ublas::matrix<value_type>& dm,
                    const std::vector<size_type>& dof_indices ) override
    {
        this->addMatrix( dm, dof_indices, dof_indices );
    }

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    void addMatrix( const value_type, MatrixSparse<value_type> const&, Feel::MatrixStructure matStruc = Feel::SAME_NONZERO_PATTERN ) override;

    void scale( const value_type ) override;

    /**
     * Return the value of the entry \p (i,j).  This may be an
     * expensive operation and you should always take care where to
     * call this function.  In order to avoid abuse, this function
     * throws an exception if the required element does not exist in
     * the matrix.
     *
     * In case you want a function that returns zero instead (for
     * entries that are not in the sparsity pattern of the matrix),
     * use the \p el function.
     */
    value_type operator()( const size_type i,
                           const size_type j ) const override;

    /**
     *
     */
    self_type& operator=( MatrixSparse<value_type> const& M ) override;

    /**
     * Returns the diagonal of the block matrix
     *
     * \param out the vector to store the diagonal
     */
    void diagonal( Vector<value_type>& out ) const override;

    /**
     * Returns the transpose of a matrix
     *
     * \param Mt the matrix transposed
     */
    void transpose( MatrixSparse<value_type>& Mt, size_type options ) const override;

    /**
     * \return \f$ v^T M u \f$
     */
    real_type
    energy( Vector<value_type> const& __v,
            Vector<value_type> const& __u,
            bool transpose = false ) const override;

    /**
     * Return the l1-norm of the matrix, that is
     * \f$|M|_1=max_{all columns j}\sum_{all rows i} |M_ij|\f$, (max. sum of columns).
     *
     * This is the natural matrix norm that is compatible to the
     * l1-norm for vectors, i.e.  \f$|Mv|_1\leq |M|_1 |v|_1\f$.
     */
    real_type l1Norm() const override;

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
    real_type linftyNorm() const override;

    /**
     * see if Sparse matrix has been closed
     * and fully assembled yet
     */
    bool closed() const override;

    /**
     * Print the contents of the matrix to the screen
     * in a uniform style, regardless of matrix/solver
     * package being used.
     */
    void print( std::ostream& os = std::cout ) const;

    /**
     * Same as the print method above, but allows you
     * to print to a stream in the standard syntax.
     */
    template <typename U>
    friend std::ostream& operator<<( std::ostream& os, const MatrixSparse<U>& m );

    /**
     * Print the contents of the matrix to the screen
     * in a package-personalized style, if available.
     */
    void printPersonal( std::ostream& /*os*/ = std::cout ) const override
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        FEELPP_ASSERT( 0 )
            .error( "invalid call" );
    }

    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name = "NULL" ) const override;

    /**
     * This function creates a matrix called "submatrix" which is defined
     * by the row and column indices given in the "rows" and "cols" entries.
     * Currently this operation is only defined for the PetscMatrix type.
     */
    void createSubmatrix( MatrixSparse<value_type>& submatrix,
                          const std::vector<size_type>& rows,
                          const std::vector<size_type>& cols ) const override
    {
        this->_get_submatrix( submatrix,
                              rows,
                              cols,
                              false ); // false means DO NOT REUSE submatrix
    }

    /**
     * This function is similar to the one above, but it allows you to reuse
     * the existing sparsity pattern of "submatrix" instead of reallocating
     * it again.  This should hopefully be more efficient if you are frequently
     * extracting submatrices of the same size.
     */
    void reinitSubmatrix( MatrixSparse<value_type>& submatrix,
                          const std::vector<size_type>& rows,
                          const std::vector<size_type>& cols ) const override
    {
        this->_get_submatrix( submatrix,
                              rows,
                              cols,
                              true ); // true means REUSE submatrix
    }

    /**
     * eliminate rows without change pattern, and put 1 on the diagonal
     * entry
     *
     *\warning if the matrix was symmetric before this operation, it
     * won't be afterwards. So use the proper solver (nonsymmetric)
     */
    void zeroRows( std::vector<int> const& rows, Vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context, value_type value_on_diagonal ) override;

    void updateBlockMat( std::shared_ptr<MatrixSparse<value_type>> const& m, std::vector<size_type> const& start_i, std::vector<size_type> const& start_j ) override;

    //@}o

  protected:
  private:
    backend_ptrtype M_backend;

    std::shared_ptr<MatrixSparse<value_type>> M_mat;
};

template <int NR, int NC, typename T>
class MatrixBlock : public MatrixBlockBase<T>
{
    typedef MatrixBlockBase<T> super_type;

  public:
    static const uint16_type NBLOCKROWS = NR;
    static const uint16_type NBLOCKCOLS = NC;
    static const uint16_type NBLOCKSIZE = NR * NC;

    typedef typename super_type::value_type value_type;
    typedef typename super_type::matrix_ptrtype matrix_ptrtype;
    typedef typename super_type::backend_type backend_type;
    typedef vf::Blocks<NBLOCKROWS, NBLOCKCOLS, matrix_ptrtype> blocks_type;
    typedef vf::BlocksBase<matrix_ptrtype> blocksbase_type;

    template<typename BackendT>
    MatrixBlock( blocksbase_type const& blockSet,
                 BackendT && b,
                 bool copy_values = true,
                 bool diag_is_nonzero = true )
        : super_type( blockSet, std::forward<BackendT>(b), copy_values, diag_is_nonzero )
    {
    }

    MatrixBlock( MatrixBlock const& mb )
        : super_type( mb )
    {
    }

    MatrixBlock operator=( MatrixBlock const& mb )
    {
        super_type::operator=( mb );
        return *this;
    }

    MatrixBlock& operator=( matrix_ptrtype const& M )
    {
        super_type::operator=( M );
        return *this;
    }

}; // MatrixBlock

} // Feel


#endif /* __MatrixBlock_H */
