/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-03

  Copyright (C) 2008, 2009 Christophe Prud'homme
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
   \file matrixpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-01-03
 */
#include <life/lifealg/vectorpetsc.hpp>
#include <life/lifealg/matrixpetsc.hpp>

#if defined( HAVE_PETSC_H )




namespace Life
{

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc()
    :
    _M_destroy_mat_on_exit(true)
{}




template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc(Mat m)
    :
    _M_destroy_mat_on_exit(false)
{
    this->_M_mat = m;
#if (PETSC_VERSION_MAJOR >= 3)
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS);
#endif
    this->setInitialized( true );
}




template <typename T>
inline
MatrixPetsc<T>::~MatrixPetsc()
{
    this->clear();
}


template <typename T>
void MatrixPetsc<T>::init (const size_type m,
                           const size_type n,
                           const size_type m_l,
                           const size_type n_l,
                           const size_type nnz,
                           const size_type /*noz*/)
{
    if ((m==0) || (n==0))
        return;

    {
        // Clear initialized matrices
        if (this->isInitialized())
            this->clear();

        this->setInitialized( true );
    }


    int ierr     = 0;
    int m_global = static_cast<int>(m);
    int n_global = static_cast<int>(n);
    int m_local  = static_cast<int>(m_l);
    int n_local  = static_cast<int>(n_l);
    int n_nz     = static_cast<int>(nnz);
    //int n_oz     = static_cast<int>(noz);

    // create a sequential matrix on one processor
    if ((m_l == m) && (n_l == n))
    {
        // Create matrix.  Revisit later to do preallocation and make more efficient
        ierr = MatCreateSeqAIJ (Application::COMM_WORLD, m_global, n_global,
                                n_nz, PETSC_NULL, &_M_mat);
        CHKERRABORT(Application::COMM_WORLD,ierr);



    }

    else
    {

        ierr = MatCreateMPIAIJ (Application::COMM_WORLD, m_local, n_local, m_global, n_global,
                                PETSC_DECIDE, PETSC_NULL, PETSC_DECIDE, PETSC_NULL, &_M_mat);
        //ierr = MatCreateMPIAIJ (Application::COMM_WORLD, m_local, n_local, m_global, n_global,
        ///n_nz, PETSC_NULL, n_oz, PETSC_NULL, &_M_mat);
        //MatCreate(Application::COMM_WORLD,m_local,n_local,m_global,n_global, &_M_mat);
        //MatCreate(Application::COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m_global,n_global, &_M_mat);
        //MatSetSizes(_M_mat,m_local,n_local,m_global,n_global);
        CHKERRABORT(Application::COMM_WORLD,ierr);

    }
    ierr = MatSetFromOptions (_M_mat);
    CHKERRABORT(Application::COMM_WORLD,ierr);

#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    ierr = MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS );
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);
#if 0
    // additional insertions will not be allowed if they generate
    // a new nonzero
    ierr = MatSetOption (_M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption (_M_mat, MAT_NEW_NONZERO_LOCATION_ERR);
    CHKERRABORT(Application::COMM_WORLD,ierr);
#endif // 0

    this->zero ();
}

template <typename T>
void MatrixPetsc<T>::init (const size_type m,
                           const size_type n,
                           const size_type m_l,
                           const size_type n_l,
                           graph_ptrtype const& graph )
{
    this->setGraph( graph );

    {
        // Clear initialized matrices
        if (this->isInitialized())
            this->clear();

        this->setInitialized(  true );
    }


    int proc_id = 0;

    MPI_Comm_rank (Application::COMM_WORLD, &proc_id);

    Debug( 7013 ) << "[MatrixPetsc::init()] m   = " << m << "\n";
    Debug( 7013 ) << "[MatrixPetsc::init()] n   = " << n << "\n";
    Debug( 7013 ) << "[MatrixPetsc::init()] m_l = " << m_l << "\n";
    Debug( 7013 ) << "[MatrixPetsc::init()] n_l = " << n_l << "\n";

    // Make sure the sparsity pattern isn't empty
    LIFE_ASSERT (this->graph()->size() == n_l)( this->graph()->size() )( n_l ).warn( "incompatible diagonal non zero pattern" );
    Debug( 7013 ) << "[MatrixPetsc::init()] graph size   = " << this->graph()->size() << "\n";
    Debug( 7013 ) << "[MatrixPetsc::init()] graph first row entry on proc   = " << this->graph()->firstRowEntryOnProc() << "\n";
    Debug( 7013 ) << "[MatrixPetsc::init()] graph last row entry on proc   = " << this->graph()->lastRowEntryOnProc() << "\n";

    if (m==0)
        return;

    int ierr     = 0;
    int m_global = static_cast<int>(m);
    int n_global = static_cast<int>(n);
    int m_local  = static_cast<int>(m_l);
    int n_local  = static_cast<int>(n_l);


    // create a sequential matrix on one processor
    if ((m_l == m) && (n_l == n))
    {
#if 0
        PetscInt nrows = m_local;
        PetscInt ncols = n_local;
        PetscInt *dnz;
        PetscInt *onz;
        MatPreallocateInitialize( Application::COMM_WORLD, nrows, ncols, dnz, onz );
        typename graph_type::iterator it = this->graph()->begin();
        typename graph_type::iterator en = this->graph()->end();
        for( ; it != en; ++it )
            {
                std::vector<PetscInt> nzcol( boost::get<2>( it->second ).size() );
                std::copy( boost::get<2>( it->second ).begin(),
                           boost::get<2>( it->second ).end(),
                           nzcol.begin() );


                MatPreallocateSet( it->first, 1, nzcol.data(), dnz, onz );
            }
        MatPreallocateFinalize( dnz, onz );
#endif
        PetscInt *dnz;
        dnz = new PetscInt[ this->graph()->nNzOnProc().size() ];
        std::copy( this->graph()->nNzOnProc().begin(),
                   this->graph()->nNzOnProc().end(),
                   dnz );

        //std::copy( dnz, dnz+this->graph()->nNzOnProc().size(), std::ostream_iterator<PetscInt>( std::cout, "\n" ) );
        ierr = MatCreateSeqAIJ (Application::COMM_WORLD, m_global, n_global,
                                0,
                                dnz,
                                //(int*) this->graph()->nNzOnProc().data(),
                                &_M_mat);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        //ierr = MatSeqAIJSetPreallocation( _M_mat, 0, (int*)this->graph()->nNzOnProc().data() );
        ierr = MatSeqAIJSetPreallocation( _M_mat, 0, dnz );
        CHKERRABORT(Application::COMM_WORLD,ierr);
    }

    else
    {
        ierr = MatCreateMPIAIJ (Application::COMM_WORLD,
                                m_local, n_local,
                                m_global, n_global,
                                0, (int*) this->graph()->nNzOnProc().data(),
                                0, (int*) this->graph()->nNzOffProc().data(), &_M_mat);
        CHKERRABORT(Application::COMM_WORLD,ierr);


    }
    ierr = MatSetFromOptions (_M_mat);
    CHKERRABORT(Application::COMM_WORLD,ierr);

#if (PETSC_VERSION_MAJOR >= 3)
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS);
#endif
#if 0
    // additional insertions will not be allowed if they generate
    // a new nonzero
    ierr = MatSetOption (_M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption (_M_mat, MAT_NEW_NONZERO_LOCATION_ERR);
    CHKERRABORT(Application::COMM_WORLD,ierr);
#endif

    //MatShift( _M_mat, 1 );
    //printMatlab( "shift.m" );

    //this->zero();
}


template <typename T>
void MatrixPetsc<T>::zero ()
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    ierr = MatZeroEntries(_M_mat);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}

template <typename T>
void MatrixPetsc<T>::zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    ierr = MatZeroEntries(_M_mat);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}


template <typename T>
void MatrixPetsc<T>::clear ()
{
    int ierr=0;

    if ((this->isInitialized()) && (this->_M_destroy_mat_on_exit))
    {
        ierr = MatDestroy (_M_mat);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        this->setInitialized( false );
    }
}




template <typename T>
inline
void MatrixPetsc<T>::close () const
{
    // BSK - 1/19/2004
    // strictly this check should be OK, but it seems to
    // fail on matrix-free matrices.  Do they falsely
    // state they are assembled?  Check with the developers...
//   if (this->closed())
//     return;

    int ierr=0;

    ierr = MatAssemblyBegin (_M_mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(Application::COMM_WORLD,ierr);
    ierr = MatAssemblyEnd   (_M_mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}



template <typename T>
inline
size_type MatrixPetsc<T>::size1 () const
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize (_M_mat, &petsc_m, &petsc_n);

    return static_cast<size_type>(petsc_m);
}



template <typename T>
inline
size_type MatrixPetsc<T>::size2 () const
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize (_M_mat, &petsc_m, &petsc_n);

    return static_cast<size_type>(petsc_n);
}



template <typename T>
inline
size_type MatrixPetsc<T>::rowStart () const
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    return static_cast<size_type>(start);
}



template <typename T>
inline
size_type MatrixPetsc<T>::rowStop () const
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    return static_cast<size_type>(stop);
}



template <typename T>
inline
void MatrixPetsc<T>::set (const size_type i,
                          const size_type j,
                          const value_type& value)
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>(value);
    ierr = MatSetValues(_M_mat, 1, &i_val, 1, &j_val,
                        &petsc_value, INSERT_VALUES);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}



template <typename T>
inline
void MatrixPetsc<T>::add (const size_type i,
                          const size_type j,
                          const value_type& value)
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );

    int ierr=0, i_val=i, j_val=j;
    //Debug( 7013 ) << "[MatrixPetsc<>::add] adding value " << value << " at (" << i << "," << j << ")\n";
    PetscScalar petsc_value = static_cast<PetscScalar>(value);


    ierr = MatSetValues(_M_mat, 1, &i_val, 1, &j_val,
                        &petsc_value, ADD_VALUES);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}
/*                                   */

template <typename T>
inline
bool MatrixPetsc<T>::closed() const
{
    LIFE_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );


    int ierr=0;
    PetscTruth assembled;

    ierr = MatAssembled(_M_mat, &assembled);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    return (assembled == PETSC_TRUE);
}
template <typename T>
void
MatrixPetsc<T>::addMatrix(const ublas::matrix<value_type>& dm,
                          const std::vector<size_type>& rows,
                          const std::vector<size_type>& cols)
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    const size_type m = dm.size1();
    const size_type n = dm.size2();

    LIFE_ASSERT (rows.size() == size1()).error( "invalid row size" );
    LIFE_ASSERT (cols.size() == size2()).error( "invalid column size" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues(_M_mat,
                        m, (int*) boost::addressof( rows[0] ),
                        n, (int*) boost::addressof( cols[0] ),
                        (PetscScalar*) dm.data().begin(),
                        ADD_VALUES);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}

// print
template <typename T>
void
MatrixPetsc<T>::printMatlab (const std::string name) const
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" );

    // assert (this->closed());
    this->close();

    int ierr=0;
    PetscViewer petsc_viewer;


    ierr = PetscViewerCreate (Application::COMM_WORLD,
                              &petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    /**
     * Create an ASCII file containing the matrix
     * if a filename was provided.
     */
    if (name != "NULL")
    {
        ierr = PetscViewerASCIIOpen( Application::COMM_WORLD,
                                     name.c_str(),
                                     &petsc_viewer);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        ierr = PetscViewerSetFormat (petsc_viewer,
                                     PETSC_VIEWER_ASCII_MATLAB);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        ierr = MatView (_M_mat, petsc_viewer);
        CHKERRABORT(Application::COMM_WORLD,ierr);
    }

    /**
     * Otherwise the matrix will be dumped to the screen.
     */
    else
    {
        ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                     PETSC_VIEWER_ASCII_MATLAB);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        ierr = MatView (_M_mat, PETSC_VIEWER_STDOUT_WORLD);
        CHKERRABORT(Application::COMM_WORLD,ierr);
    }


    /**
     * Destroy the viewer.
     */
    ierr = PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}

template <typename T>
inline
void
MatrixPetsc<T>::addMatrix (const T a_in, MatrixSparse<T> &X_in)
{
    LIFE_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    // sanity check. but this cannot avoid
    // crash due to incompatible sparsity structure...
    LIFE_ASSERT( this->size1() == X_in.size1())( this->size1() )( X_in.size1() ).error( "incompatible dimension" );
    LIFE_ASSERT( this->size2() == X_in.size2())( this->size2() )( X_in.size2() ).error( "incompatible dimension" );

    PetscScalar     a = static_cast<PetscScalar>      (a_in);
    MatrixPetsc<T>* X = dynamic_cast<MatrixPetsc<T>*> (&X_in);

    LIFE_ASSERT (X != 0).error( "invalid petsc matrix" );

    int ierr=0;

    // the matrix from which we copy the values has to be assembled/closed
    X->close ();

// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = MatAXPY(&a,  X->_M_mat, _M_mat, (MatStructure)SAME_NONZERO_PATTERN);
    CHKERRABORT(Application::COMM_WORLD,ierr);

// 2.3.x & newer
#else

    //ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)SAME_NONZERO_PATTERN);
    ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)SUBSET_NONZERO_PATTERN );
    CHKERRABORT(Application::COMM_WORLD,ierr);

#endif
}

template <typename T>
void
MatrixPetsc<T>::scale( T const a )
{
    int ierr = MatScale( _M_mat, a );
    CHKERRABORT(Application::COMM_WORLD,ierr);
}
template <typename T>
inline
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::l1Norm() const
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    int ierr=0;
    double petsc_value;
    real_type value;

    assert (this->closed());

    ierr = MatNorm(_M_mat, NORM_1, &petsc_value);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    value = static_cast<real_type>(petsc_value);

    return value;
}

template <typename T>
inline
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::linftyNorm() const
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    int ierr=0;
    double petsc_value;
    real_type value;

    assert (this->closed());

    ierr = MatNorm(_M_mat, NORM_INFINITY, &petsc_value);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    value = static_cast<real_type>(petsc_value);

    return value;
}

template <typename T>
MatrixPetsc<T> &
MatrixPetsc<T>::operator = ( MatrixSparse<value_type> const& M )
{
    MatrixPetsc<T> const* X = dynamic_cast<MatrixPetsc<T> const*> (&M);

    //MatConvert(X->mat(), MATSAME, MAT_INITIAL_MATRIX, &_M_mat);
    MatDuplicate(X->mat(),MAT_COPY_VALUES,&_M_mat);

    return *this;
}
template <typename T>
inline
typename MatrixPetsc<T>::value_type
MatrixPetsc<T>::operator () (const size_type i,
                             const size_type j) const
{
    LIFE_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

#if (((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR >= 2) && (PETSC_VERSION_SUBMINOR >= 1)) || \
    ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR >= 3))  || ( PETSC_VERSION_MAJOR >= 3 ) )
    // PETSc 2.2.1 & newer
    const PetscScalar *petsc_row;
    const PetscInt    *petsc_cols;

#else
    // PETSc 2.2.0 & older
    PetscScalar *petsc_row;
    int* petsc_cols;

#endif

    T value=0.;

    int
        ierr=0,
        ncols=0,
        i_val=static_cast<int>(i),
        j_val=static_cast<int>(j);


    // the matrix needs to be closed for this to work
    this->close();

    ierr = MatGetRow(_M_mat, i_val, &ncols, &petsc_cols, &petsc_row);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Perform a binary search to find the contiguous index in
    // petsc_cols (resp. petsc_row) corresponding to global index j_val
    std::pair<const int*, const int*> p =
        std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);

    // Found an entry for j_val
    if (p.first != p.second)
    {
        // The entry in the contiguous row corresponding
        // to the j_val column of interest
        const int j = std::distance (const_cast<int*>(&petsc_cols[0]),
                                     const_cast<int*>(p.first));

        assert (j < ncols);
        assert (petsc_cols[j] == j_val);

        value = static_cast<T> (petsc_row[j]);

        ierr  = MatRestoreRow(_M_mat, i_val,
                              &ncols, &petsc_cols, &petsc_row);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        return value;
    }

    // Otherwise the entry is not in the sparse matrix,
    // i.e. it is 0.
    return 0.;
}


template<typename T>
void
MatrixPetsc<T>::zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context )
{
#if (PETSC_VERSION_MAJOR >= 3)
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS);
#endif
    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    if ( on_context.test( ON_ELIMINATION_KEEP_DIAGONAL ) )
        {
            VectorPetsc<value_type> diag( this->size1(), stop-start );
            MatGetDiagonal( _M_mat, diag.vec() );
            MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0);
            MatDiagonalSet( _M_mat, diag.vec(), INSERT_VALUES );
            for( size_type i = 0; i < rows.size(); ++i )
                {
                    // eliminate column

                    // warning: a row index may belong to another
                    // processor, so make sure that we access only the
                    // rows that belong to this processor
                    if ( rows[i] >= start && rows[i] < stop )
                        rhs.set( rows[i], values[i]*diag(rows[i]) );
                }
        }
    else
        {




            MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0);
            for( size_type i = 0; i < rows.size(); ++i )
                {
                    // eliminate column

                    // warning: a row index may belong to another
                    // processor, so make sure that we access only the
                    // rows that belong to this processor
                    if ( rows[i] >= start && rows[i] < stop )
                        rhs.set( rows[i], values[i] );
                }
        }
    //if ( on_context.test( ON_ELIMINATION_SYMMETRIC ) )
    if ( 0 )
        {
            Debug() << "symmetrize zero-out  operation\n";
            for( size_type i = 0; i < rows.size(); ++i )
                {
                    int myRow = rows[i];
                    int ncols;
                    const PetscInt* cols;
                    const PetscScalar *vals;
                    MatGetRow( _M_mat, myRow, &ncols, &cols, &vals );
                    // we suppose here that the pattern is
                    // symmetric (not necessarily the matrix)

                    // do a loop on all indices, j, in Indices and
                    // zero out corresponding columns
                    // (Indices[j],myRow)
                    for( int j = 0; j < ncols; ++j )
                        {
                            if ( cols[j] == myRow )
                                continue;

                            int ncols2;
                            const PetscInt* cols2;
                            const PetscScalar *vals2;

                            MatGetRow( _M_mat, cols[j], &ncols2, &cols2, &vals2 );
                            PetscScalar* v = new PetscScalar[ncols2];
                            std::copy( vals2, vals2+ncols, v );
                            bool found = false;
                            for( int k = 0; k < ncols2; ++k )
                                if ( cols2[k] == myRow )
                                    {
                                        found=true;
                                        rhs.add( cols[j], -values[i]*vals2[k] );
                                        //int j_val = cols2[k];
                                        v[k] = 0;
                                        MatSetValuesRow( _M_mat, myRow,v);

                                        break;
                                    }
                            MatRestoreRow( _M_mat, cols[j], &ncols2, &cols2, &vals2 );
                            LIFE_ASSERT( found == true )( myRow )( cols[j] ).error ( "matrix graph is not symmetric" );
                        }
                    MatRestoreRow( _M_mat, myRow, &ncols, &cols, &vals );
                } // i
        }


}

template<typename T>
void
MatrixPetsc<T>::transpose( MatrixSparse<value_type>& Mt ) const
{
    MatrixPetsc<T>* Atrans = dynamic_cast<MatrixPetsc<T>*> (&Mt);

  
#if (PETSC_VERSION_MAJOR >= 3)
    int ierr = MatTranspose( _M_mat, MAT_INITIAL_MATRIX,&Atrans->_M_mat );
#else
    int ierr = MatTranspose( _M_mat, &Atrans->_M_mat );
#endif

    CHKERRABORT(Application::COMM_WORLD,ierr);

    PetscTruth isSymmetric;
    MatEqual( _M_mat, Atrans->_M_mat, &isSymmetric);
    if (isSymmetric) {
#if (PETSC_VERSION_MAJOR >= 3)
    MatSetOption(_M_mat,MAT_SYMMETRIC,PETSC_TRUE);
#else
    MatSetOption(_M_mat,MAT_SYMMETRIC );
#endif
    } else {
        PetscPrintf(PETSC_COMM_WORLD,"Warning: Petsc matrix is non-symmetric \n");
    }
}

template<typename T>
void
MatrixPetsc<T>::symmetricPart( MatrixSparse<value_type>& Mt ) const
{
    int ierr = 0;

    // first check if the matrix is symmetric
    Mat Atrans;
    MatDuplicate(_M_mat,MAT_COPY_VALUES,&Atrans);

#if (PETSC_VERSION_MAJOR >= 3)
    MatTranspose( _M_mat, MAT_INITIAL_MATRIX, &Atrans );
#else
    MatTranspose( _M_mat, &Atrans );
#endif
    PetscTruth isSymmetric;
    MatEqual( _M_mat, Atrans, &isSymmetric);
    MatDestroy( Atrans );

    if (isSymmetric) {
        Log() << "[MatrixPetsc<T>::symmetricPart] Matrix is already symmetric, don't do anything\n";
#if (PETSC_VERSION_MAJOR >= 3)
	MatSetOption(_M_mat,MAT_SYMMETRIC,PETSC_TRUE);
#else
	MatSetOption(_M_mat,MAT_SYMMETRIC );
#endif


        // duplicate the matrix
        MatrixPetsc<T>* B = dynamic_cast<MatrixPetsc<T>*> (&Mt);
        ierr = MatDuplicate( _M_mat, MAT_COPY_VALUES, &B->_M_mat );
#if (PETSC_VERSION_MAJOR >= 3)
	MatSetOption( B->_M_mat,MAT_SYMMETRIC,PETSC_TRUE);
#else
	MatSetOption( B->_M_mat,MAT_SYMMETRIC );
#endif
        return ;
    }

    // The matrix was not symmetric, compute 1/2 (A + A^T )
    Mat A[2];
    ierr = MatDuplicate(_M_mat,MAT_COPY_VALUES,&A[0]);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = MatDuplicate(_M_mat,MAT_COPY_VALUES,&A[1]);
    CHKERRABORT(Application::COMM_WORLD,ierr);

#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( A[1], MAT_INITIAL_MATRIX, &A[1] );
#else
    MatTranspose( A[1], &A[1] );
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);

    MatrixPetsc<T>* B = dynamic_cast<MatrixPetsc<T>*> (&Mt);
    LIFE_ASSERT( B ).error ( "[symmetricPart] invalid dynamic_cast<MatrixPetsc>" );

    ierr = MatCreateComposite(PETSC_COMM_WORLD,2,A,&B->_M_mat);
    CHKERRABORT(Application::COMM_WORLD,ierr);


//A Completer pour PETSC < 3 !!!!!!!!
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatCompositeSetType(B->_M_mat,MAT_COMPOSITE_ADDITIVE);
#endif
    
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = MatCompositeMerge(B->_M_mat);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = MatScale( B->_M_mat, 0.5 );
    CHKERRABORT(Application::COMM_WORLD,ierr);


    // check that B is now symmetric
    Mat Btrans;
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( B->_M_mat, MAT_INITIAL_MATRIX, &Btrans );
#else
    ierr = MatTranspose( B->_M_mat, &Btrans );
#endif
    
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = MatEqual( B->_M_mat, Btrans, &isSymmetric);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    if (isSymmetric) {
        Log() << "[MatrixPetsc<T>::symmetricPart] symmetric part is really symmetric\n";
#if (PETSC_VERSION_MAJOR >= 3)
        ierr = MatSetOption(B->_M_mat,MAT_SYMMETRIC,PETSC_TRUE);
#else
	ierr = MatSetOption(B->_M_mat,MAT_SYMMETRIC);
#endif
        CHKERRABORT(Application::COMM_WORLD,ierr);

    } else {
        PetscPrintf(PETSC_COMM_WORLD,"Warning: Petsc matrix is non-symmetric \n");
    }
    ierr = MatDestroy (Btrans);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}

template class MatrixPetsc<double>;
} // Life

#endif // HAVE_PETSC_H
