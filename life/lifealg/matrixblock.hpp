/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-10

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-10
 */

#ifndef __MatrixBlock_H
#define __MatrixBlock_H 1

#include <boost/fusion/vector.hpp>


namespace Life
{
/**
 * \class MatrixBlock
 * \brief block of matrices
 *
 * <code>
 * Matrix<2,2, MatrixType> MB( fusion::make_vector( A, B, C, D ) );
 * </code>
 *
 * it is likely that MatrixBlock<> won't be used directly but will be
 * constructed by a free-function, e.g. block<NR,NC>( make_vector() ),
 * to ease the matrix assembly
 *
 * <code>
 * form2( Xh, Xh, block<2,2>( A, B, C, D  ) ) = integrate(...);
 * </code>
 *
 * or alternatively
 *
 * <code>
 *  // store the object MatrixBlock, we can reuse it
 *  AUTO( myblock, block<2,2>( A, B, C, D  ) );
 *  form( Xh, Xh , myblock ) = ...;
 *  form( Xh, Xh , myblock ) += ...;
 * </code>
 *
 * @author Christophe Prud'homme
 */
template<int NR, int NC, typename MatrixType>
class MatrixBlock : public MatrixSparse<typename MatrixSparse::value_type>
{
    typedef MatrixSparse<typename MatrixSparse::value_type> super;
public:

    /** @name Constants
     */
    //@{

    //! number of rows
    static const uint16_type NBLOCKROWS = NR;

    //! number of columns
    static const uint16_type NBLOCKCOLS = NC;

    static const uint16_type NBLOCKSIZE = NR * NC;
    //@}


    /** @name Typedefs
     */
    //@{

    typedef MatrixType matrix_type;
    typedef typename MatrixSparse::value_type value_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    typedef fusion::vector<NBLOCKSIZE, matrix_type> vector_matrix_type;
    typedef fusion::vector<NBLOCKSIZE, matrix_ptrtype> vector_matrix_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixBlock( vector_matrix_ptrtype const & v )
        :
        super(),
        M_v( v )
    {}

    MatrixBlock( MatrixBlock const & mb )
        :
        super( mb ),
        M_v( mb.M_v )
    {}

    ~MatrixBlock()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    MatrixBlock operator=( MatrixBlock const& mb )
    {
        if ( this != &mb )
            {
                M_v = mb.M_v;
            }
        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                graph_ptrtype const& graph );

    /**
     * \return true if matrix has a graph, false otherwise
     */
    bool hasGraph() const { return _M_graph != 0; }

    /**
     * \return the graph associated to the sparse matrix
     */
    graph_ptrtype const& graph() const { return _M_graph; }

    /**
     * set the graph associated to the sparse matrix
     */
    void setGraph( graph_ptrtype const& graph ) { _M_graph = graph; }

    /**
     * Release all memory and return to a state just like after having
     * called the default constructor.
     */
    void clear ();

    /**
     * Set all entries to 0.
     */
    void zero ();

    /**
     * Call the Sparse assemble routines.  sends necessary messages to
     * other processors
     */
    void close () const;

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size1 () const;

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size2 () const;

    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    size_type rowStart () const;

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    size_type rowStop () const;

    /**
     * Set the element \p (i,j) to \p value.
     * Throws an error if the entry does
     * not exist. Still, it is allowed to store
     * zero values in non-existent fields.
     */
    void set (const size_type i,
                      const size_type j,
                      const value_type& value);

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    void add (const size_type i,
                      const size_type j,
                      const value_type& value);

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix (const ublas::matrix<value_type> &dm,
                            const std::vector<size_type> &rows,
                            const std::vector<size_type> &cols);

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    void addMatrix (const ublas::matrix<value_type> &dm,
                            const std::vector<size_type> &dof_indices);

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    void addMatrix (const T, MatrixSparse<T> &);

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
    T operator () (const size_type i,
                           const size_type j) const;

    /**
     * Return the l1-norm of the matrix, that is
     * \f$|M|_1=max_{all columns j}\sum_{all rows i} |M_ij|\f$, (max. sum of columns).
     *
     * This is the natural matrix norm that is compatible to the
     * l1-norm for vectors, i.e.  \f$|Mv|_1\leq |M|_1 |v|_1\f$.
     */
    real_type l1Norm () const;

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
    real_type linftyNorm () const;

    /**
     * see if Sparse matrix has been closed
     * and fully assembled yet
     */
    bool closed() const;

    /**
     * Print the contents of the matrix to the screen
     * in a uniform style, regardless of matrix/solver
     * package being used.
     */
    void print(std::ostream& os=std::cout) const;

    /**
     * Same as the print method above, but allows you
     * to print to a stream in the standard syntax.
     */
    template <typename U>
    friend std::ostream& operator << (std::ostream& os, const MatrixSparse<U>& m);

    /**
     * Print the contents of the matrix to the screen
     * in a package-personalized style, if available.
     */
    void printPersonal(std::ostream& /*os*/=std::cout) const
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        LIFE_ASSERT( 0 ).error("invalid call");
    }

    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab(const std::string name="NULL") const
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        std::cerr << "ERROR writing MATLAB file " << name << std::endl;
        LIFE_ASSERT( 0 ).error("invalid call");
    }

    /**
     * This function creates a matrix called "submatrix" which is defined
     * by the row and column indices given in the "rows" and "cols" entries.
     * Currently this operation is only defined for the PetscMatrix type.
     */
    void createSubmatrix(MatrixSparse<T>& submatrix,
                                 const std::vector<size_type>& rows,
                                 const std::vector<size_type>& cols) const
    {
        this->_get_submatrix(submatrix,
                             rows,
                             cols,
                             false); // false means DO NOT REUSE submatrix
    }

    /**
     * This function is similar to the one above, but it allows you to reuse
     * the existing sparsity pattern of "submatrix" instead of reallocating
     * it again.  This should hopefully be more efficient if you are frequently
     * extracting submatrices of the same size.
     */
    void reinitSubmatrix(MatrixSparse<T>& submatrix,
                                 const std::vector<size_type>& rows,
                                 const std::vector<size_type>& cols) const
    {
        this->_get_submatrix(submatrix,
                             rows,
                             cols,
                             true); // true means REUSE submatrix
    }

    /**
     * eliminate rows without change pattern, and put 1 on the diagonal
     * entry
     *
     *\warning if the matrix was symmetric before this operation, it
     * won't be afterwards. So use the proper solver (nonsymmetric)
     */
    void zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context );


    //@}



protected:

private:

};
}
} // Life
#endif /* __MatrixBlock_H */
