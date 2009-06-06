/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Universit� Joseph Fourier (Grenoble I)

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
   \file matrixsparse.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-11-27
 */
// $Id: sparse_matrix.h,v 1.12 2005/05/24 12:54:57 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __sparse_matrix_h__
#define __sparse_matrix_h__

// C++ includes
#include <iomanip>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>



#include <life/lifecore/life.hpp>
#include <life/lifecore/traits.hpp>
#include <life/lifecore/context.hpp>

#include <life/lifealg/enums.hpp>
#include <life/lifealg/graphcsr.hpp>
#include <life/lifealg/vector.hpp>

namespace Life
{
namespace ublas = boost::numeric::ublas;

// forward declarations
template <typename T> class MatrixSparse;
template <typename T> inline std::ostream& operator << (std::ostream& os, const MatrixSparse<T>& m);


/**
 * Generic sparse matrix. This class contains
 * pure virtual members that must be overloaded
 * in derived classes.  Using a derived class
 * allows for uniform access to sparse matrices
 * from various different solver packages in
 * different formats.
 *
 * @author Benjamin S. Kirk, 2003
 * @author Christophe Prud'homme, 2005
 */
template <typename T>
class MatrixSparse
{
public:

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;
    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    typedef Vector<T> vector_type;
    typedef boost::shared_ptr<Vector<T> > vector_ptrtype;

    /**
     * Constructor; initializes the matrix to be empty, without any
     * structure, i.e.  the matrix is not usable at all. This
     * constructor is therefore only useful for matrices which are
     * members of a class. All other matrices should be created at a
     * point in the data flow where all necessary information is
     * available.
     *
     * You have to initialize
     * the matrix before usage with
     * \p init(...).
     */
    MatrixSparse ();

    /**
     * Destructor. Free all memory, but do not release the memory of
     * the sparsity structure.
     */
    virtual ~MatrixSparse ();

    /**
     * @returns true if the matrix has been initialized,
     * false otherwise.
     */
    virtual bool isInitialized() const { return _M_is_initialized; }

    /**
     * Updates the matrix sparsity pattern. When your \p
     * MatrixSparse<T> implementation does not need this data simply
     * do not overload this method.
     */
    virtual void updateSparsityPattern (const std::vector<std::vector<size_type> >&) {}

    /**
     * Initialize a Sparse matrix that is of global dimension \f$ m
     * \times n \f$ with local dimensions \f$ m_l \times n_l \f$.  \p
     * nnz is the number of on-processor nonzeros per row (defaults to
     * 30).  \p noz is the number of on-processor nonzeros per row
     * (defaults to 10).
     */
    virtual void init (const size_type m,
                       const size_type n,
                       const size_type m_l,
                       const size_type n_l,
                       const size_type nnz=30,
                       const size_type noz=10) = 0;

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    virtual void init ( const size_type m,
                        const size_type n,
                        const size_type m_l,
                        const size_type n_l,
                        graph_ptrtype const& graph ) = 0;

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
    virtual void clear () = 0;

    /**
     * Set all entries to 0.
     */
    virtual void zero () = 0;

    /**
     * Set entries between to 0.
     */
    virtual void zero ( size_type start1, size_type size1,
                        size_type start2, size_type size2 ) = 0;

    /**
     * Call the Sparse assemble routines.  sends necessary messages to
     * other processors
     */
    virtual void close () const = 0;

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    virtual size_type size1 () const = 0;

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    virtual size_type size2 () const = 0;

    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    virtual size_type rowStart () const = 0;

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    virtual size_type rowStop () const = 0;

    /**
     * Set the element \p (i,j) to \p value.
     * Throws an error if the entry does
     * not exist. Still, it is allowed to store
     * zero values in non-existent fields.
     */
    virtual void set (const size_type i,
                      const size_type j,
                      const value_type& value) = 0;

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    virtual void add (const size_type i,
                      const size_type j,
                      const value_type& value) = 0;

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    virtual void addMatrix (const ublas::matrix<value_type> &dm,
                            const std::vector<size_type> &rows,
                            const std::vector<size_type> &cols) = 0;

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    virtual void addMatrix (const ublas::matrix<value_type> &dm,
                            const std::vector<size_type> &dof_indices) = 0;

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    virtual void addMatrix (const T, MatrixSparse<T> &) = 0;

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    void addMatrix (const T& s, boost::shared_ptr<MatrixSparse<T> > & m )
    {
        this->addMatrix( s, *m );
    }

    virtual void scale ( const T ) = 0;

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
    virtual T operator () (const size_type i,
                           const size_type j) const = 0;

    virtual MatrixSparse<T> & operator = ( MatrixSparse<value_type> const& M ) = 0;
    MatrixSparse<T> & operator = ( boost::shared_ptr<MatrixSparse<value_type> > const& M )
    {
        *this = *M;
        return *this;
    }
    /**
     * compute the A scalar product \f$v^T A u\f$
     *
     * \param u a vector
     * \param v a vector
     * \param transpose true to compute \f$v^T A^T u\f$
     * \return the energy \f$v^T A u\f$
     */
    virtual real_type energy ( vector_type const& v,
                               vector_type const& u,
                               bool transpose = false ) const = 0;

    /**
     * Compute the scalar product \f$(Au, v)= v^T A u\f$
     *
     * \param u a vector
     * \param v a vector
     * \param transpose true to compute \f$v^T A^T u\f$ instead, false otherwise
     * \return the energy \f$v^T A u\f$
     */
    virtual real_type energy ( vector_ptrtype const& v,
                               vector_ptrtype const& u,
                               bool transpose = false ) const
    {
        return this->energy( *v, *u, transpose );
    }

    /**
     * Return the l1-norm of the matrix, that is
     * \f$|M|_1=max_{all columns j}\sum_{all rows i} |M_ij|\f$, (max. sum of columns).
     *
     * This is the natural matrix norm that is compatible to the
     * l1-norm for vectors, i.e.  \f$|Mv|_1\leq |M|_1 |v|_1\f$.
     */
    virtual real_type l1Norm () const = 0;

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
    virtual real_type linftyNorm () const = 0;

    /**
     * see if Sparse matrix has been closed
     * and fully assembled yet
     */
    virtual bool closed() const = 0;

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
    virtual void printPersonal(std::ostream& /*os*/=std::cout) const
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
    virtual void printMatlab(const std::string name="NULL") const
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
    virtual void createSubmatrix(MatrixSparse<T>& submatrix,
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
    virtual void reinitSubmatrix(MatrixSparse<T>& submatrix,
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
    virtual void zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context ) = 0;



    /**
     * set initialized only for subclasses
     */
    void setInitialized( bool init )
    {
        _M_is_initialized = init;
    }
protected:
    /**
     * Protected implementation of the create_submatrix and reinit_submatrix
     * routines.  Note that this function must be redefined in derived classes
     * for it to work properly!
     */
    virtual void _get_submatrix(MatrixSparse<T>& ,
                                const std::vector<size_type>& ,
                                const std::vector<size_type>& ,
                                const bool) const
    {
        std::cerr << "Error! This function is not yet implemented in the base class!"
                  << std::endl;
        LIFE_ASSERT( 0 ).error("invalid call");
    }

    /**
     * Flag indicating whether or not the matrix
     * has been initialized.
     */
    bool _M_is_initialized;

    graph_ptrtype _M_graph;
};



//-----------------------------------------------------------------------
// MatrixSparse inline members
template <typename T>
inline
MatrixSparse<T>::MatrixSparse () :
    _M_is_initialized(false)
{}



template <typename T>
inline
MatrixSparse<T>::~MatrixSparse ()
{}



template <typename T>
inline
void MatrixSparse<T>::print(std::ostream& os) const
{
    assert (this->isInitialized());

    for (size_type i=0; i<this->m(); i++)
    {
        for (size_type j=0; j<this->n(); j++)
            os << std::setw(8) << (*this)(i,j) << " ";
        os << std::endl;
    }
}



#if 0
// Full specialization for Complex datatypes
template <>
inline
void MatrixSparse<std::complex<double> >::print(std::ostream& os) const
{
    // std::complex<>::operator<<() is defined, but use this form

    std::cout << "Real part:" << std::endl;
    for (size_type i=0; i<this->size1(); i++)
    {
        for (size_type j=0; j<this->size2(); j++)
            os << std::setw(8) << (*this)(i,j).real() << " ";
        os << std::endl;
    }

    os << std::endl << "Imaginary part:" << std::endl;
    for (size_type i=0; i<this->size1(); i++)
    {
        for (size_type j=0; j<this->size2(); j++)
            os << std::setw(8) << (*this)(i,j).imag() << " ";
        os << std::endl;
    }
}
#endif


// For SGI MIPSpro this implementation must occur after
// the partial specialization of the print() member.
template <typename T>
inline
std::ostream& operator << (std::ostream& os, const MatrixSparse<T>& m)
{
    m.print(os);
    return os;
}

} // Life

#endif // #ifndef __sparse_matrix_h__

