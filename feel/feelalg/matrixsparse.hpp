/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
#include <boost/fusion/include/fold.hpp>

#include <Eigen/Core>



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelcore/context.hpp>

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/graphcsr.hpp>
#include <feel/feelalg/vector.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

// forward declarations
template <typename T> class MatrixSparse;
template <typename T> inline std::ostream& operator << ( std::ostream& os, const MatrixSparse<T>& m );

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

    typedef DataMap datamap_type;
    typedef boost::shared_ptr<datamap_type> datamap_ptrtype;

    typedef typename datamap_type::indexsplit_type indexsplit_type;
    typedef typename datamap_type::indexsplit_ptrtype indexsplit_ptrtype;

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
    MatrixSparse( WorldComm const& worldComm=Environment::worldComm() );

    MatrixSparse( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, WorldComm const& worldComm=Environment::worldComm() );

    /**
     * Destructor. Free all memory, but do not release the memory of
     * the sparsity structure.
     */
    virtual ~MatrixSparse ();

    /**
     * Return datamap for rows
     */
    datamap_type const& mapRow() const
    {
        return *M_mapRow;
    }

    /**
     * Return datamap for cols
     */
    datamap_type const& mapCol() const
    {
        return *M_mapCol;
    }

    /**
     * Return datamap for rows
     */
    datamap_ptrtype const& mapRowPtr() const
    {
        return M_mapRow;
    }

    /**
     * Return datamap for cols
     */
    datamap_ptrtype const& mapColPtr() const
    {
        return M_mapCol;
    }

    void setMapRow( datamap_ptrtype const& d )
    {
        M_mapRow=d;
    }
    void setMapCol( datamap_ptrtype const& d )
    {
        M_mapCol=d;
    }

    /**
     * @returns true if the matrix has been initialized,
     * false otherwise.
     */
    virtual bool isInitialized() const
    {
        return M_is_initialized;
    }

    /**
     * Updates the matrix sparsity pattern. When your \p
     * MatrixSparse<T> implementation does not need this data simply
     * do not overload this method.
     */
    virtual void updateSparsityPattern ( const std::vector<std::vector<size_type> >& ) {}

    /**
     * Initialize a Sparse matrix that is of global dimension \f$ m
     * \times n \f$ with local dimensions \f$ m_l \times n_l \f$.  \p
     * nnz is the number of on-processor nonzeros per row (defaults to
     * 30).  \p noz is the number of on-processor nonzeros per row
     * (defaults to 10).
     */
    virtual void init ( const size_type m,
                        const size_type n,
                        const size_type m_l,
                        const size_type n_l,
                        const size_type nnz=30,
                        const size_type noz=10 ) = 0;

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    virtual void init ( const size_type m,
                        const size_type n,
                        const size_type m_l,
                        const size_type n_l,
                        graph_ptrtype const& graph ) = 0;

    /**
     * set the indexSplit associated to the sparse matrix
     */
    virtual void setIndexSplit( indexsplit_ptrtype const& is )
    {
        M_indexSplit = is;
    }

    /**
     * \return the indexSplit associated to the sparse matrix
     */
    indexsplit_ptrtype const& indexSplit() const
    {
        return M_indexSplit;
    }

    /**
     * \return true if matrix has a graph, false otherwise
     */
    bool hasGraph() const
    {
        return M_graph != 0;
    }

    /**
     * \return the graph associated to the sparse matrix
     */
    graph_ptrtype const& graph() const
    {
        return M_graph;
    }

    /**
     * set the graph associated to the sparse matrix
     */
    void setGraph( graph_ptrtype const& graph )
    {
        M_graph = graph;
    }

    /**
     * set matrix properties, @see MatrixProperties
     */
    void setMatrixProperties( size_type p )
    {
        M_mprop = ( size_type )p;
        checkProperties();
    }
    /**
     * \return true if matrix is hermitian, false otherwise
     */
    bool isHermitian() const
    {
        checkProperties();
        return M_mprop.test( HERMITIAN );
    }

    /**
     * \return true if matrix is non hermitian, false otherwise
     */
    bool isNonHermitian() const
    {
        checkProperties();
        return M_mprop.test( NON_HERMITIAN );
    }

    /**
     * \return true if matrix is positive definite, false otherwise
     */
    bool isHermitianPositiveDefinite() const
    {
        checkProperties();
        return M_mprop.test( HERMITIAN | POSITIVE_DEFINITE );
    }

    /**
     * \return true if matrix is singular, false otherwise
     */
    bool isSingular() const
    {
        checkProperties();
        return M_mprop.test( SINGULAR );
    }

    /**
     * \return true if matrix is singular, false otherwise
     */
    bool isPositiveDefinite() const
    {
        checkProperties();
        return M_mprop.test( POSITIVE_DEFINITE );
    }

    bool haveConsistentProperties() const
    {
        bool p1 = M_mprop.test( SINGULAR ) && M_mprop.test( POSITIVE_DEFINITE );
        bool p2 = M_mprop.test( HERMITIAN ) && M_mprop.test( NON_HERMITIAN );
        return ( p1 == false ) && ( p2 == false );
    }
    bool isDense() const
    {
        return M_mprop.test( DENSE );
    }

    /**
     * \return true if matrix is symmetric, false otherwise
     */
    virtual bool isSymmetric() const;

    /**
     * \return true if \p this is the transpose of Trans, false otherwise
     */
    virtual bool isTransposeOf ( MatrixSparse<value_type> &Trans ) const;


    void checkProperties() const
    {
        if ( !haveConsistentProperties() )
        {
            std::ostringstream ostr;
            ostr << "Invalid matrix properties:\n"
                 << "           HERMITIAN: " << isHermitian() << "\n"
                 << "       NON_HERMITIAN: " << isNonHermitian() << "\n"
                 << "            SINGULAR: " << isSingular() << "\n"
                 << "   POSITIVE_DEFINITE: " << isPositiveDefinite() << "\n"
                 << "               DENSE: " << isDense() << "\n";
            throw std::logic_error( ostr.str() );
        }
    }
    /**
     * \return the communicator
     */
    //mpi::communicator const& comm() const { return M_comm; }
    WorldComm const& comm() const
    {
        return M_worldComm;
    }

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
    virtual void set ( const size_type i,
                       const size_type j,
                       const value_type& value ) = 0;

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    virtual void add ( const size_type i,
                       const size_type j,
                       const value_type& value ) = 0;

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    virtual void addMatrix ( const ublas::matrix<value_type> &dm,
                             const std::vector<size_type> &rows,
                             const std::vector<size_type> &cols ) = 0;

    /**
     * Add the full matrix to the
     * Sparse matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    virtual void addMatrix ( int* rows, int nrows,
                             int* cols, int ncols,
                             value_type* data ) = 0;

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    virtual void addMatrix ( const ublas::matrix<value_type> &dm,
                             const std::vector<size_type> &dof_indices ) = 0;

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    virtual void addMatrix ( const T, MatrixSparse<T> & ) = 0;

    /**
     * Add a Sparse matrix \p _X, scaled with \p _a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = \_a*\_X + \texttt{this} \f$.
     */
    void addMatrix ( const T& s, boost::shared_ptr<MatrixSparse<T> > & m )
    {
        this->addMatrix( s, *m );
    }

    virtual void scale ( const T ) = 0;

    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    void multVector ( const Vector<T>& arg,
                      Vector<T>& dest ) const;

    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    void multVector ( const boost::shared_ptr<Vector<T> >& arg,
                      boost::shared_ptr<Vector<T> >& dest ) const
    {
        this->multVector( *arg, *dest );
    }

    /**
     * Multiplies the matrix with \p arg and adds the result to \p dest.
     */
    void multAddVector ( const Vector<T>& arg,
                         Vector<T>& dest ) const;


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
    virtual T operator () ( const size_type i,
                            const size_type j ) const = 0;

    virtual MatrixSparse<T> & operator = ( MatrixSparse<value_type> const& M ) = 0;
    MatrixSparse<T> & operator = ( boost::shared_ptr<MatrixSparse<value_type> > const& M )
    {
        *this = *M;
        return *this;
    }
    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    virtual void diagonal ( Vector<T>& dest ) const = 0;

    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    void diagonal ( boost::shared_ptr<Vector<T> >& dest ) const
    {
        diagonal( *dest );
    }

    /**
     * Returns the transpose of a matrix
     *
     * \param Mt the matrix transposed
     */
    virtual void transpose( MatrixSparse<value_type>& Mt, size_type options = MATRIX_TRANSPOSE_ASSEMBLED ) const = 0;

    /**
     * \return the transpose of the matrix
     */
    boost::shared_ptr<MatrixSparse<T> > transpose( size_type options = MATRIX_TRANSPOSE_ASSEMBLED ) const
    {
        boost::shared_ptr<MatrixSparse<T> > Mt;
        transpose( *Mt, options );
        return Mt;
    }

    /**
     * Returns the transpose of a matrix
     *
     * \param M the matrix to transpose
     * \param Mt the matrix transposed
     */
    void transpose( boost::shared_ptr<MatrixSparse<value_type> >& Mt, size_type options = MATRIX_TRANSPOSE_ASSEMBLED ) const
    {
        this->transpose( *Mt, options );
    }

    /**
     * Returns the symmetric part of the matrix
     */
    virtual void symmetricPart( MatrixSparse<value_type>& Ms ) const {}

    /**
     * Returns the symmetric part of the matrix
     */
    void symmetricPart( boost::shared_ptr<MatrixSparse<value_type> >& Ms ) const
    {
        symmetricPart( *Ms );
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
    real_type energy ( vector_ptrtype const& v,
                       vector_ptrtype const& u,
                       bool _transpose = false ) const
    {
        return this->energy( *v, *u, _transpose );
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
    void print( std::ostream& os=std::cout ) const;

    /**
     * Same as the print method above, but allows you
     * to print to a stream in the standard syntax.
     */
    template <typename U>
    friend std::ostream& operator << ( std::ostream& os, const MatrixSparse<U>& m );

    /**
     * Print the contents of the matrix to the screen
     * in a package-personalized style, if available.
     */
    virtual void printPersonal( std::ostream& /*os*/=std::cout ) const
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        FEELPP_ASSERT( 0 ).error( "invalid call" );
    }

    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    virtual void printMatlab( const std::string name="NULL" ) const
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        std::cerr << "ERROR writing MATLAB file " << name << std::endl;
        FEELPP_ASSERT( 0 ).error( "invalid call" );
    }

    /**
     * This function creates a matrix called "submatrix" which is defined
     * by the row and column indices given in the "rows" and "cols" entries.
     * useSameDataMap : (opimisation) put at true if dataMapCol == dataMapRow and rows == cols for each proc
     * checkAndFixRange : add missing dof entries in // ( typically a ghost dof present but not active dof associated )
     */
    virtual
    boost::shared_ptr<MatrixSparse<T> >
    createSubMatrix( std::vector<size_type> const& rows,
                     std::vector<size_type> const& cols,
                     bool useSameDataMap=false,
                     bool checkAndFixRange=true ) const
    {
        CHECK( false ) << "invalid call : Not Implemented in base class";
        boost::shared_ptr<MatrixSparse<T> > res;
        return res;
    }

    /**
     * This function creates a matrix called "submatrix" which is defined
     * by the row and column indices given in the "rows" and "cols" entries.
     * Currently this operation is only defined for the PetscMatrix type.
     */
    virtual void createSubmatrix( MatrixSparse<T>& submatrix,
                                  const std::vector<size_type>& rows,
                                  const std::vector<size_type>& cols ) const
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
    virtual void reinitSubmatrix( MatrixSparse<T>& submatrix,
                                  const std::vector<size_type>& rows,
                                  const std::vector<size_type>& cols ) const
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
    virtual void zeroRows( std::vector<int> const& rows, Vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context ) = 0;

    /**
     * update a block matrix
     */
    virtual void  updateBlockMat( boost::shared_ptr<MatrixSparse<T> > m, std::vector<size_type> start_i, std::vector<size_type> start_j ) = 0;


    /**
     * set initialized only for subclasses
     */
    void setInitialized( bool _init )
    {
        M_is_initialized = _init;
    }
#if 0
    template<typename DomainSpace, typename ImageSpace>
    void updateIndexSplit( DomainSpace const& dm, ImageSpace const& im )
    {
        auto nSpace = DomainSpace::element_type::nSpaces;

        if ( nSpace>1 )
        {
            //std::cout << "\n Debug : nSpace " << nSpace << "\n";
            std::vector < std::vector<int> > is( nSpace );
            uint cptSpaces=0;
            uint start=0;
            auto result = boost::make_tuple( cptSpaces,start );

            std::vector < std::vector<int> > indexSplit( nSpace );
            //detail::computeNDofForEachSpace cndof(nSpace);
            detail::computeNDofForEachSpace cndof( indexSplit );
            boost::fusion::fold( dm->functionSpaces(), result,  cndof );

            this->setIndexSplit( indexSplit );
        }
    }
#endif
    virtual void sqrt( MatrixSparse<value_type>& _m ) const;

    void sqrt( boost::shared_ptr<MatrixSparse<value_type> >& _m ) const
    {
        sqrt(*_m);
    }

    virtual void matMatMult ( MatrixSparse<value_type> const& In, MatrixSparse<value_type> &Res );

    virtual void matInverse ( MatrixSparse<value_type> &Inv );

    virtual void applyInverseSqrt( Vector<value_type>& vec_in, Vector<value_type>& vec_out );

protected:
    /**
     * Protected implementation of the create_submatrix and reinit_submatrix
     * routines.  Note that this function must be redefined in derived classes
     * for it to work properly!
     */
    virtual void _get_submatrix( MatrixSparse<T>& ,
                                 const std::vector<size_type>& ,
                                 const std::vector<size_type>& ,
                                 const bool ) const
    {
        std::cerr << "Error! This function is not yet implemented in the base class!"
                  << std::endl;
        FEELPP_ASSERT( 0 ).error( "invalid call" );
    }

    //! mpi communicator
    //mpi::communicator M_comm;
    WorldComm M_worldComm;

    /**
     * Flag indicating whether or not the matrix
     * has been initialized.
     */
    bool M_is_initialized;

    graph_ptrtype M_graph;

    Context M_mprop;

    indexsplit_ptrtype M_indexSplit;

    /**
     * data distribution map of the vector over the processors
     */
    datamap_ptrtype M_mapRow;
    datamap_ptrtype M_mapCol;


};

typedef MatrixSparse<double> d_sparse_matrix_type;
typedef boost::shared_ptr<d_sparse_matrix_type> d_sparse_matrix_ptrtype;
typedef boost::shared_ptr<d_sparse_matrix_type> sparse_matrix_ptrtype;


//-----------------------------------------------------------------------
// MatrixSparse inline members
template <typename T>
inline
MatrixSparse<T>::MatrixSparse( WorldComm const& worldComm ) :
    M_worldComm( worldComm ),
    M_is_initialized( false ),
    M_mprop( NON_HERMITIAN )
{}

template <typename T>
inline
MatrixSparse<T>::MatrixSparse ( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, WorldComm const& worldComm ) :
    M_worldComm( worldComm ),
    M_is_initialized( false ),
    M_mprop( NON_HERMITIAN ),
    M_mapRow( dmRow ),
    M_mapCol( dmCol )
{}

template <typename T>
inline
MatrixSparse<T>::~MatrixSparse ()
{}



template <typename T>
inline
void MatrixSparse<T>::print( std::ostream& os ) const
{
    assert ( this->isInitialized() );

    for ( size_type i=0; i<this->size1(); i++ )
    {
        for ( size_type j=0; j<this->size2(); j++ )
            os << std::setw( 8 ) << ( *this )( i,j ) << " ";

        os << std::endl;
    }
}

template <typename T>
void MatrixSparse<T>::multVector ( const Vector<T>& arg,
                                   Vector<T>& dest ) const
{
    dest.zero();
    this->multAddVector( arg,dest );
}



template <typename T>
void MatrixSparse<T>::multAddVector ( const Vector<T>& arg,
                                      Vector<T>& dest ) const
{
    /* This functionality is actually implemented in the \p
       Vector class.  */
    dest.addVector( arg,*this );
}


#if 0
// Full specialization for Complex datatypes
template <>
inline
void MatrixSparse<std::complex<double> >::print( std::ostream& os ) const
{
    // std::complex<>::operator<<() is defined, but use this form

    std::cout << "Real part:" << std::endl;

    for ( size_type i=0; i<this->size1(); i++ )
    {
        for ( size_type j=0; j<this->size2(); j++ )
            os << std::setw( 8 ) << ( *this )( i,j ).real() << " ";

        os << std::endl;
    }

    os << std::endl << "Imaginary part:" << std::endl;

    for ( size_type i=0; i<this->size1(); i++ )
    {
        for ( size_type j=0; j<this->size2(); j++ )
            os << std::setw( 8 ) << ( *this )( i,j ).imag() << " ";

        os << std::endl;
    }
}
#endif


// For SGI MIPSpro this implementation must occur after
// the partial specialization of the print() member.
template <typename T>
inline
std::ostream& operator << ( std::ostream& os, const MatrixSparse<T>& m )
{
    m.print( os );
    return os;
}

namespace detail
{
template <class MatrixType>
struct is_matrix_ptr : mpl::false_ {};

template <class MatrixType>
struct is_matrix_ptr<boost::shared_ptr<MatrixType> >
        :
        boost::is_base_of<MatrixSparse<typename MatrixType::value_type>,
        MatrixType>
{};
}

template <typename T>
void MatrixSparse<T>::sqrt( MatrixSparse<value_type>& _m ) const
{
    std::cerr << "Error! This function is not yet implemented in the base class!"
              << std::endl;
    FEELPP_ASSERT( 0 ).error( "invalid call" );
}

template <typename T>
void MatrixSparse<T>::matMatMult ( MatrixSparse<value_type> const& In, MatrixSparse<value_type> &Res )
{
    std::cerr << "Error! This function is not yet implemented in the base class!"
              << std::endl;
    FEELPP_ASSERT( 0 ).error( "invalid call" );
}

template <typename T>
void MatrixSparse<T>::matInverse ( MatrixSparse<value_type> &Inv )
{
    std::cerr << "Error! This function is not yet implemented in the base class!"
              << std::endl;
    FEELPP_ASSERT( 0 ).error( "invalid call" );
}

template <typename T>
void MatrixSparse<T>::applyInverseSqrt( Vector<value_type>& vec_in, Vector<value_type>& vec_out )
{
    std::cerr << "Error! This function is not yet implemented in the base class!"
              << std::endl;
    FEELPP_ASSERT( 0 ).error( "invalid call" );
}

template <typename T>
bool MatrixSparse<T>::isSymmetric () const
{
    std::cerr << "Error! This function is not yet implemented in the base class!"
              << std::endl;
    FEELPP_ASSERT( 0 ).error( "invalid call" );

    return 0;
}

template <typename T>
bool MatrixSparse<T>::isTransposeOf ( MatrixSparse<value_type> &Trans ) const
{
    std::cerr << "Error! This function is not yet implemented in the base class!"
              << std::endl;
    FEELPP_ASSERT( 0 ).error( "invalid call" );

    return 0;
}

} // Feel

#endif // #ifndef __sparse_matrix_h__
