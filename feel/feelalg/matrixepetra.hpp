/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s):
  Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Klaus Sapelza <klaus.sapelza@epfl.ch>
       Date: 2006-09-14

  Copyright (C) 2006,2007 EPFL
  Copyright (C) 2006,2007 Universit� Joseph Fourier (Grenoble I)

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
   \file matrixepetra.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Klaus.Sapelza <klaus.sapelza@epfl.ch>
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2005-09-14
 */
#ifndef __MatrixEpetra_H
#define __MatrixEpetra_H 1

#include <feel/feelconfig.h>


#include <feel/feelcore/application.hpp>

#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vectorepetra.hpp>

#if defined( FEELPP_HAS_TRILINOS_EPETRA )
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#if defined(FEELPP_HAS_MPI)
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif /* FEELPP_HAS_MPI */
#include <Epetra_Map.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION

namespace Feel
{
template<typename T> class VectorEpetra;

/**
 * \class MatrixEpetra
 * \brief Wrapper for epetra crs matrices
 *
 * Epetra matrix. Provides a nice interface to the
 * Epetra C-based data structures for parallel,
 * sparse matrices.
 *
 * @author Christophe Prud'homme
 * @author Klaus Sapelza
 * @author Goncalo Pena
 * @see
 */
class MatrixEpetra : public MatrixSparse<double>
{
    typedef MatrixSparse<double> super;


public:
    /** @name Typedefs
     */

    //@{
    static const bool is_row_major = true;
    typedef super::value_type value_type;
    typedef super::real_type real_type;
    typedef std::vector<std::set<size_type> > pattern_type;

    typedef super::graph_type graph_type;
    typedef super::graph_ptrtype graph_ptrtype;

    typedef Vector<value_type> vector_type;
    typedef boost::shared_ptr<Vector<value_type> > vector_ptrtype;

    typedef VectorEpetra<value_type> epetra_vector_type;
    typedef boost::shared_ptr<epetra_vector_type> epetra_vector_ptrtype;

    typedef MatrixEpetra epetra_sparse_matrix_type;

    //@}


    /** @name Constructors, destructor
     */
    //@{

    /**
     * Constructor.  Creates a EpetraMatrix assuming you already
     * have a valid Epetra_FECrsMatrix object.  In this case, m is NOT destroyed
     * by the EpetraMatrix destructor when this object goes out of scope.
     * This allows ownership of m to remain with the original creator,
     * and to simply provide additional functionality with the EpetraMatrix.
     */

    MatrixEpetra ( Epetra_FECrsMatrix const& m );

    /**
     * Copy constructor
     */
    MatrixEpetra ( MatrixEpetra const& T );

    /**
     * Other constructors
     */
    MatrixEpetra ( Epetra_Map const& emap, int nnz = 50 );

    MatrixEpetra ( Epetra_Vector const& x );

    MatrixEpetra ( Epetra_Map const& rowmap, Epetra_Map const& colmap );

    MatrixEpetra( Epetra_Map const& row_emap, Epetra_Map const& col_emap,
                  Epetra_Map const& dom_map, Epetra_Map const& range_map );


    MatrixEpetra ( const size_type m,
                   const size_type n,
                   const size_type m_l,
                   const size_type n_l,
                   const size_type nnz,
                   const size_type /*noz*/ );

    /**
     * Destructor. Free all memory, but do not
     * release the memory of the sparsity
     * structure.
     */

    ~MatrixEpetra();

    //@}

    /** @name Operator overloads
     */
    //@{

    MatrixEpetra & operator = ( MatrixSparse<value_type> const& M )
    {
        return *this;
    }
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
    value_type operator () ( const size_type i,
                             const size_type j ) const;


    /**
     * Overload of operator =
     * In this way, the pointer associated with MatrixEpetra doesn't point to the
     * same copied matrix, but to a new one
    */
    MatrixEpetra& operator=( const MatrixEpetra& T )
    {
        if ( &T != this )
        {
            M_emap = T.getRowMap();
            M_col_emap = T.getColMap();
            M_dom_map = T.mat().DomainMap();
            M_range_map = T.mat().RangeMap();
            *M_mat = T.mat();

            this->setInitialized( T.isInitialized() );
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

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
     * returns the row map
     */
    Epetra_Map getRowMap() const
    {
        return M_emap;
    }

    /**
     * returns the column map
     */
    Epetra_Map getColMap() const
    {
        return M_col_emap;
    }

    /**
     * returns the domain map
     */
    Epetra_Map getDomainMap() const
    {
        return M_mat->DomainMap();
    }

    /**
     * returns the range map
     */
    Epetra_Map getRangeMap() const
    {
        return M_mat->RangeMap();
    }
    //@}


    std::vector<size_type> bcIndices()
    {
        return M_bc_index;
    }

    boost::shared_ptr<Epetra_FECrsMatrix> matrix()
    {
        return M_mat;
    }

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Initialize a Epetra matrix that is of global
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
                const size_type nnz=10,
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
     * Initialize using sparsity structure computed by \p dof_map.
     */
    void init ();

    /**
     * reinitialize the matrix
     */
    //     void resize( size_type m, size_type n, bool /*preserve*/ )
    //     {
    //         this->setInitialized( false );
    //         init( m, n, m, n );
    //     }

    //     void fill( pattern_type const& patt )
    //     {}

    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear ();


    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ();

    /**
     * Set all entries to 0 in the
     * block[start1,stop1]x[start2,stop2]. This method retains
     * sparsity structure.
     */
    void zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
    {
    }

    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    void diagonal ( Vector<double>& dest ) const;

    void setDiagonal( Epetra_Vector const& x );

    /**
     * Call the Epetra assemble routines.
     * sends necessary messages to other
     * processors
     * GlobalAssemble calls FillComplete().
     * This does really an assemble: Gather
     * any overlapping /shared data into the
     * non-overlapping partitioning defied by
     * the map that was passed to this matrix
     * at construiction time. Data imported
     * from other processors is stored on the
     * owning processor with a "sumInto"
     * or accumulate operation.
     * This is a collective operation -- every processor must enter it before any will complete it.
     */

    void close () const;


    /**
     * see if Epetra matrix has been closed
     * and fully assembled yet
     */
    bool closed() const;

    bool EpetraIndicesAreLocal() const;

    bool HaveColMap() const;

    /**
     * Returns the transpose of a matrix
     *
     * \param Mt the matrix transposed
     */
    void transpose( MatrixSparse<value_type>& Mt ) const;

    /**
     * Return the energy v1^T M v2
     */
    real_type energy ( vector_type const& v1, vector_type const& v2, bool tranpose = false ) const;

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
     * Set the element \p (i,j) to \p value.
     * Throws an error if the entry does
     * not exist. Still, it is allowed to store
     * zero values in non-existent fields.
     */
    void set ( const size_type i,
               const size_type j,
               const value_type& value );

    /**
     * Multiplies matrix A by B and stores the result
     */
    void multiplyMatrix ( const MatrixEpetra& A, const MatrixEpetra& B );

    template<typename T>
    void multiply ( bool trans, const Vector<T>& v,  Vector<T>& r ) const
    {
        epetra_vector_type const& ev( dynamic_cast<epetra_vector_type const&>( v ) );
        epetra_vector_type& er( dynamic_cast<epetra_vector_type&>( r ) );
        M_mat->Multiply( trans, ev.vec(), er.vec() );
    }

    void multiply ( bool trans, const Epetra_MultiVector &X, Epetra_MultiVector &Y ) const
    {
        M_mat->Multiply( trans, X, Y );
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
               const value_type& value );

    /**
     * Add the full matrix to the
     * Epetra matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix ( double coeff, const MatrixEpetra& _m )
    {
        EpetraExt::MatrixMatrix::Add( _m.mat(), false, coeff, ( *this ).mat(), 1. );
    }

    void addMatrix ( const ublas::matrix<value_type> &/*dm*/,
                     const std::vector<size_type> &/*rows*/,
                     const std::vector<size_type> &/*cols*/ )
    {
    }

    /**
     * Add the full matrix to the
     * Petsc matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     *
     * \warning the \p data array is a matrix and must be ROW_MAJOR
     */
    void addMatrix ( int* rows, int nrows,
                     int* cols, int ncols,
                     value_type* data );

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    void addMatrix ( const ublas::matrix<value_type> &dm,
                     const std::vector<size_type> &dof_indices )
    {
        this->addMatrix ( dm, dof_indices, dof_indices );
    }


    /**
     * Add a Sparse matrix \p X, scaled with \p a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = a*X + \texttt{this} \f$.
     * Use this with caution, the sparse matrices need to have the
     * same nonzero pattern, otherwise \p Epetra will crash!
     * It is advisable to not only allocate appropriate memory with
     * \p init() , but also explicitly zero the terms of \p this
     * whenever you add a non-zero value to \p X.  Note: \p X will
     * be closed, if not already done, before performing any work.
     */
    void addMatrix ( const value_type a, super &X );


    /**
     * Multiplies all the entries in the matrix by scalar
     */
    void scale ( double const scalar );


    /**
     * Returns the raw Epetra matrix context pointer.
     */
    Epetra_FECrsMatrix& mat ()
    {
        //FEELPP_ASSERT (M_mat != NULL).error("null epetra matrix");
        return *M_mat;
    }

    /**
     * Returns the raw Epetra matrix context pointer.
     */
    Epetra_FECrsMatrix const& mat () const
    {
        //FEELPP_ASSERT (M_mat != NULL).error("null epetra matrix");
        return *M_mat;
    }

    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL" ) const;
    void printKonsole() const;


    /**
     * eliminate row without change pattern, and put 1 on the diagonal
     * entry
     */

    //@}

    /**
     * eliminate rows without change pattern, and put 1 on the diagonal
     * entry
     *
     *\warning if the matrix was symmetric before this operation, it
     * won't be afterwards. So use the proper solver (nonsymmetric)
     */
    void zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context );


    /**
     * update a block matrix
     */
    void updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m, std::vector<size_type> start_i, std::vector<size_type> start_j );

protected:

private:

    mpi::communicator M_comm;

    /**
     * Constructor; initializes the matrix to
     * be empty, without any structure, i.e.
     * the matrix is not usable at all. This
     * constructor is therefore only useful
     * for matrices which are members of a
     * class. All other matrices should be
     * created at a point in the data flow
     * where all necessary information is
     * available.
     *
     * You have to initialize
     * the matrix before usage with
     * \p init(...).
     */
    MatrixEpetra();

    /**
     * Epetra matrix datatype to store values
     */

    Epetra_Map M_emap;
    Epetra_Map M_col_emap;
    Epetra_Map M_dom_map;
    Epetra_Map M_range_map;

    mutable boost::shared_ptr<Epetra_FECrsMatrix> M_mat;

    std::vector<size_type> M_bc_index;
};



inline
MatrixEpetra::MatrixEpetra()
    :
    super(),

#ifdef FEELPP_HAS_MPI
    M_emap( Epetra_Map( -1, 0, 0, Epetra_MpiComm( M_comm ) ) ),
    M_col_emap( Epetra_Map( -1, 0, 0, Epetra_MpiComm( M_comm ) ) ),
    M_dom_map( Epetra_Map( -1, 0, 0, Epetra_MpiComm( M_comm ) ) ),
    M_range_map( Epetra_Map( -1, 0, 0, Epetra_MpiComm( M_comm ) ) ),
    M_mat( new Epetra_FECrsMatrix( Copy, M_emap, 0 ) )
#else
    M_emap( Epetra_Map( -1, 0, 0, Epetra_SerialComm ) ),
    M_col_emap( Epetra_Map( -1, 0, 0, Epetra_SerialComm ) ),
    M_dom_map( Epetra_Map( -1, 0, 0, Epetra_SerialComm ) ),
    M_range_map( Epetra_Map( -1, 0, 0, Epetra_SerialComm ) ),
    M_mat( new Epetra_FECrsMatrix( Copy, M_emap, 0 ) )

#endif
{}

inline
MatrixEpetra::MatrixEpetra( Epetra_Map const& emap, int nnz )
    :
    super(),
    M_emap( emap ),
    M_col_emap( emap ),
    M_dom_map( emap ),
    M_range_map( emap ),
    M_mat( new Epetra_FECrsMatrix( Copy, emap, 0 ) )
{
    //this->setInitialized( true );
}


inline
MatrixEpetra::MatrixEpetra( Epetra_Map const& row_emap, Epetra_Map const& col_emap )
    :
    super(),
    M_emap( row_emap ),
    M_col_emap( col_emap ),
    M_dom_map( col_emap ),
    M_range_map( row_emap ),
    M_mat( new Epetra_FECrsMatrix( Copy, row_emap, col_emap, 0 ) )
{
    //this->setInitialized( true );

}

inline
MatrixEpetra::MatrixEpetra( Epetra_Map const& row_emap, Epetra_Map const& col_emap,
                            Epetra_Map const& dom_map, Epetra_Map const& range_map )
    :
    super(),
    M_emap( row_emap ),
    M_col_emap( col_emap ),
    M_dom_map( dom_map ),
    M_range_map( range_map ),
    M_mat( new Epetra_FECrsMatrix( Copy, M_emap, M_col_emap, 0 ) )
{
    //this->setInitialized( true );

}


inline
MatrixEpetra::MatrixEpetra( Epetra_FECrsMatrix const& M )
    :
    super(),
    M_emap( M.RowMap() ),
    M_col_emap( M.ColMap() ),
    M_dom_map( M.DomainMap() ),
    M_range_map( M.RangeMap() ),
    M_mat( new Epetra_FECrsMatrix( Copy, M_emap, M_col_emap, 0 ) )
{
}


inline
MatrixEpetra::MatrixEpetra( MatrixEpetra const& T )
    :
    super(),
    M_emap( T.M_emap ),
    M_col_emap( T.M_col_emap ),
    M_dom_map( T.M_dom_map ),
    M_range_map( T.M_range_map ),
    M_mat( new Epetra_FECrsMatrix( T.mat() ) )

{
}


inline
MatrixEpetra::MatrixEpetra( const size_type m,
                            const size_type n,
                            const size_type m_l,
                            const size_type /*n_l*/,
                            const size_type nnz,
                            const size_type /*noz*/ )
    :
    super(),
    M_emap( Epetra_Map( m, m_l, 0, Epetra_MpiComm( M_comm ) ) ),
    M_col_emap( Epetra_Map( n, n, 0, Epetra_MpiComm( M_comm ) ) ),
    M_dom_map( M_emap ),
    M_range_map( M_emap ),
    M_mat( new Epetra_FECrsMatrix( Copy, M_emap, M_col_emap, nnz ) )
{
}


inline
MatrixEpetra::~MatrixEpetra()
{

}


inline
void MatrixEpetra::init ( const size_type m,
                          const size_type n,
                          const size_type /*m_l*/,
                          const size_type /*n_l*/,
                          const size_type /*nnz*/,
                          const size_type /*noz*/ )
{
    if ( ( m==0 ) || ( n==0 ) )
        return;

    {
        // Clear initialized matrices
        if ( this->isInitialized() )
            this->clear();

        this->setInitialized( true );
    }
}





inline
void MatrixEpetra::init ()
{
    if ( this->isInitialized() )
        this->clear();

    this->setInitialized( true );
}

inline
void MatrixEpetra::zero ()
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "epetra matrix not properly initialized" ) ;

    M_bc_index.resize( 0 );
    M_mat->PutScalar( 0.0 );

}


inline
void MatrixEpetra::clear ()
{
    if ( this->isInitialized() )
    {
        M_mat = boost::shared_ptr<Epetra_FECrsMatrix>( new Epetra_FECrsMatrix( Copy, M_emap, M_col_emap, 50 ) );

        this->setInitialized( false );
    }
}




inline
void MatrixEpetra::close () const
{
    int ierr=0;

    if ( !this->closed() )
        ierr =  M_mat->GlobalAssemble( M_dom_map, M_range_map, true );

    if ( ierr != 0 )
    {
        DVLOG(2) << "ERRCODE GlobalAssemble: " << ierr << "\n";
    }
}


inline
void MatrixEpetra::setDiagonal ( Epetra_Vector const& x )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );;

    const Epetra_Map& rowMap( M_mat->RowMatrixRowMap() );
    const Epetra_Map& colMap( M_mat->RowMatrixColMap() );

    size_type L = x.MyLength();

    value_type zero_value = static_cast<value_type>( 0.0 );

    for ( size_type i=0; i< L; i++ )
    {
        int i_val = static_cast<int>( rowMap.GID( i ) );
        int j_val = static_cast<int>( colMap.GID( i ) );

        M_mat->InsertGlobalValues( 1, &i_val, 1,  &j_val, &zero_value );
    }

    this->close();

    M_mat->ReplaceDiagonalValues( x );
}



inline
size_type MatrixEpetra::size1 () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );;

    DVLOG(2) << "Size in size1(): " << M_mat->NumGlobalRows() << "\n";

    int epetra_m = M_mat->NumGlobalRows();

    return static_cast<size_type>( epetra_m );
}



inline
size_type MatrixEpetra::size2 () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );;

    int epetra_n = M_mat->NumGlobalCols();

    return static_cast<size_type>( epetra_n );
}



inline
size_type MatrixEpetra::rowStart () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );;

    int start = M_emap.MinMyGID();

    return static_cast<size_type>( start );
}




inline
size_type MatrixEpetra::rowStop () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );;

    int stop = M_emap.MaxMyGID();

    return static_cast<size_type>( stop );
}


inline
void MatrixEpetra::set ( const size_type i,
                         const size_type j,
                         const value_type& value )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );;

    int i_val = static_cast<int>( i );
    int j_val = static_cast<int>( j );

    int ierr = 0;

    value_type epetra_value = static_cast<value_type>( value );

    ierr = M_mat->ReplaceGlobalValues( 1, &i_val, 1,  &j_val, &epetra_value );

    if ( ierr )
    {
        ierr = M_mat->InsertGlobalValues( 1, &i_val, 1,  &j_val, &epetra_value );

        if ( ierr != 0 )
        {
            DVLOG(2) << "ERRORCODE InsertGlobalValues: " << ierr <<  " in M(" << i_val << "," << j_val << ") for value "<< epetra_value << "." << "\n";
        }
    }

    //FEELPP_ASSERT( ierr == 0 )( ierr ).warn ( "invalid MatrixEpetra::set operation" );
}




inline
bool MatrixEpetra::closed() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );

    bool filled;
    filled = M_mat->Filled();

    return static_cast<bool>( filled );
}


inline
bool MatrixEpetra::EpetraIndicesAreLocal() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );

    bool IsLocal;
    IsLocal = M_mat->IndicesAreLocal();

    return static_cast<bool>( IsLocal );
}


inline
bool MatrixEpetra::HaveColMap() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixEpetra<> not properly initialized" );

    bool HaveMap;
    HaveMap = M_mat->HaveColMap();

    return static_cast<bool>( HaveMap );
}

inline
void
MatrixEpetra::addMatrix ( const value_type coeff, super &X )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "epetra matrix not initialized" );

    epetra_sparse_matrix_type* A_ptr = dynamic_cast< epetra_sparse_matrix_type*>( &X );

    this->addMatrix( coeff, *A_ptr );
}

inline
void
MatrixEpetra::scale ( double scalar )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "epetra matrix not initialized" );
    this->close();
    M_mat->Scale( scalar );
}

inline
MatrixEpetra::real_type
MatrixEpetra::l1Norm() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "epetra matrix not initialized" );
    this->close();

    real_type NormOne = M_mat->NormOne();

    return static_cast<real_type>( NormOne );
}

inline
MatrixEpetra::real_type
MatrixEpetra::linftyNorm() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "epetra matrix not initialized" );
    this->close();

    real_type NormInf = M_mat->NormInf();

    return static_cast<real_type>( NormInf );
}


// This is for a GLOBAL matrix!

inline
MatrixEpetra::value_type   //doesnt work in parallel yet
MatrixEpetra::operator () ( const size_type i,
                            const size_type j ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "epetra matrix not initialized" );


    int i_val=static_cast<int>( i ),
        j_val=static_cast<int>( j );

    int    NumEntries;
    double* Values;
    int* Indices;


    // int ierr = M_mat->ExtractMyRowView( i_val, NumEntries, Values, Indices);
    M_mat->ExtractMyRowView( i_val, NumEntries, Values, Indices );

    for ( int k = 0; k < NumEntries; ++k )
    {
        if ( Indices[k] == j_val )
            return static_cast<double> ( Values[k] );
    }

    /*
        int *epetra_cols;
        double *epetra_row;

        int
            ncols=0,
            i_val=static_cast<int>(i),
            j_val=static_cast<int>(j);

        // the matrix must NOT be closed
        //    this->close();


        M_mat->ExtractGlobalRowView(i_val, ncols, epetra_row, epetra_cols);

        // Perform a binary search to find the contiguous index in
        // petsc_cols (resp. petsc_row) corresponding to global index j_val
        std::pair<const int*, const int*> p =  std::equal_range (&epetra_cols[0], &epetra_cols[0] + ncols, j_val);

        std::cout << "Epetra row: " << *epetra_row << "\n";
        std::cout << "Epetra col: " << *epetra_cols << "\n";


        // Found an entry for j_val
        if (p.first != p.second)
            {
                // The entry in the contiguous row corresponding
                // to the j_val column of interest
                const int j = std::distance (const_cast<int*>(&epetra_cols[0]),
                                             const_cast<int*>(p.first));

                assert (j < ncols);
                assert (epetra_cols[j] == j_val);

                return static_cast<double> (epetra_row[j]);

                //return value;
            }

        std::cout << "Element " << epetra_row[j] << "\n";*/

    // Otherwise the entry is not in the sparse matrix,
    // i.e. it is 0.
    return 0.;

    //return M_mat[i_val][j_val]; //no good: locally indexed!
}

inline
void
MatrixEpetra::updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m, std::vector<size_type> start_i, std::vector<size_type> start_j )
{
    LOG(ERROR) << "Invalid call to updateBlockMat, not yet implemented\n";
}


inline
NdebugStream&
operator<<( NdebugStream& __os, MatrixEpetra const& /*__n*/ )
{
    return __os;
}

} // Feel

#endif /* FEELPP_HAS_TRILINOS_EPETRA */
#endif /* __MatrixEpetra_H */
