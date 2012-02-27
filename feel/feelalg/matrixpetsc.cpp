/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-03

  Copyright (C) 2008-2011 Christophe Prud'homme
  Copyright (C) 2008-2011 Universite Joseph Fourier (Grenoble I)

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
#include <boost/timer.hpp>
#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>

#if defined( HAVE_PETSC_H )




namespace Feel
{

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc()
    :
    _M_destroy_mat_on_exit(true)
{}

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc(DataMap const& dmRow, DataMap const& dmCol, WorldComm const& worldComm)
    :
    super(dmRow,dmCol,worldComm),
    _M_destroy_mat_on_exit(true)
{}

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc(Mat m)
    :
    _M_destroy_mat_on_exit(false)
{
    this->_M_mat = m;
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption(_M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
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
    //std::cout << "\nSEQUENTIAL : init without graph"<<std::endl;

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
        ierr = MatCreateSeqAIJ (this->comm(), m_global, n_global,
                                n_nz, PETSC_NULL, &_M_mat);
        CHKERRABORT(this->comm(),ierr);



    }

    else
    {

        ierr = MatCreateMPIAIJ (this->comm(), m_local, n_local, m_global, n_global,
                                PETSC_DECIDE, PETSC_NULL, PETSC_DECIDE, PETSC_NULL, &_M_mat);
        //ierr = MatCreateMPIAIJ (this->comm(), m_local, n_local, m_global, n_global,
        ///n_nz, PETSC_NULL, n_oz, PETSC_NULL, &_M_mat);
        //MatCreate(this->comm(),m_local,n_local,m_global,n_global, &_M_mat);
        //MatCreate(this->comm(),PETSC_DECIDE,PETSC_DECIDE,m_global,n_global, &_M_mat);
        //MatSetSizes(_M_mat,m_local,n_local,m_global,n_global);
        CHKERRABORT(this->comm(),ierr);

    }
    ierr = MatSetFromOptions (_M_mat);
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    ierr = MatSetOption(_M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(_M_mat,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    ierr = MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    ierr = MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS );
#endif
    CHKERRABORT(this->comm(),ierr);
#if 1
    // additional insertions will not be allowed if they generate
    // a new nonzero
    //ierr = MatSetOption (_M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    //CHKERRABORT(this->comm(),ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption (_M_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    CHKERRABORT(this->comm(),ierr);
#endif // 0

    this->zero ();
    //this->zeroEntriesDiagonal();
}

template <typename T>
void MatrixPetsc<T>::init (const size_type m,
                           const size_type n,
                           const size_type m_l,
                           const size_type n_l,
                           graph_ptrtype const& graph )
{
    //std::cout << "\n SEQUENTIAL : init with graph"<<std::endl;
    this->setGraph( graph );

    {
        // Clear initialized matrices
        if (this->isInitialized())
            this->clear();

        this->setInitialized(  true );
    }

    int proc_id = 0;

    MPI_Comm_rank (this->comm(), &proc_id);

    Debug( 7013 ) << "[MatrixPETSc::init()] m   = " << m << "\n";
    Debug( 7013 ) << "[MatrixPETSc::init()] n   = " << n << "\n";
    Debug( 7013 ) << "[MatrixPETSc::init()] m_l = " << m_l << "\n";
    Debug( 7013 ) << "[MatrixPETSc::init()] n_l = " << n_l << "\n";

    // Make sure the sparsity pattern isn't empty
    FEEL_ASSERT (this->graph()->size() == n_l)( this->graph()->size() )( n_l ).warn( "incompatible diagonal non zero pattern" );
    Debug( 7013 ) << "[MatrixPETSc::init()] graph size   = " << this->graph()->size() << "\n";
    Debug( 7013 ) << "[MatrixPETSc::init()] graph first row entry on proc   = " << this->graph()->firstRowEntryOnProc() << "\n";
    Debug( 7013 ) << "[MatrixPETSc::init()] graph last row entry on proc   = " << this->graph()->lastRowEntryOnProc() << "\n";

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
#if 1 // MatCreateSeqAIJ
#if 0
        PetscInt nrows = m_local;
        PetscInt ncols = n_local;
        PetscInt *dnz;
        PetscInt *onz;
        MatPreallocateInitialize( this->comm(), nrows, ncols, dnz, onz );
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
        ierr = MatCreateSeqAIJ (this->comm(), m_global, n_global,
                                0,
                                dnz,
                                //(int*) this->graph()->nNzOnProc().data(),
                                &_M_mat);
        CHKERRABORT(this->comm(),ierr);

        //ierr = MatSeqAIJSetPreallocation( _M_mat, 0, (int*)this->graph()->nNzOnProc().data() );
#if 0
        ierr = MatSeqAIJSetPreallocation( _M_mat, 0, dnz );
#else
        this->graph()->close();
        //std::cout << "sizes:" << this->graph()->ia().size() << "," << this->graph()->ja().size()  << std::endl;
        _M_ia.resize( this->graph()->ia().size() );
        _M_ja.resize( this->graph()->ja().size() );
        std::copy( this->graph()->ia().begin(), this->graph()->ia().end(), _M_ia.begin() );
        std::copy( this->graph()->ja().begin(), this->graph()->ja().end(), _M_ja.begin() );
        //std::for_each( ia.begin(), ia.end(), [](const int& i ){ std::cout << i << std::endl; } );

        ierr = MatSeqAIJSetPreallocationCSR( _M_mat,
                                             _M_ia.data(), _M_ja.data(),
                                             this->graph()->a().data() );
#endif
        CHKERRABORT(this->comm(),ierr);

#else // MatCreateSeqAIJWithArrays

        //std::cout << "matrix csr creating...\n" << std::endl;
        this->graph()->close();
        //std::cout << "sizes:" << this->graph()->ia().size() << "," << this->graph()->ja().size()  << std::endl;
        _M_ia.resize( this->graph()->ia().size() );
        _M_ja.resize( this->graph()->ja().size() );
        std::copy( this->graph()->ia().begin(), this->graph()->ia().end(), _M_ia.begin() );
        std::copy( this->graph()->ja().begin(), this->graph()->ja().end(), _M_ja.begin() );
        //std::for_each( ia.begin(), ia.end(), [](const int& i ){ std::cout << i << std::endl; } );

        ierr = MatCreateSeqAIJWithArrays( this->comm(),
                                          m_global, n_global,
                                          _M_ia.data(), _M_ja.data(), this->graph()->a().data(), &_M_mat );
        CHKERRABORT(this->comm(),ierr);
        //std::cout << "matrix csr created\n" << std::endl;
#endif
    }

    else
    {
        ierr = MatCreateMPIAIJ (this->comm(),
                                m_local, n_local,
                                m_global, n_global,
                                0, (int*) this->graph()->nNzOnProc().data(),
                                0, (int*) this->graph()->nNzOffProc().data(), &_M_mat);
        CHKERRABORT(this->comm(),ierr);


    }

    ierr = MatSetFromOptions (_M_mat);
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption(_M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(_M_mat,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS);
#endif


#if 0
    // additional insertions will not be allowed if they generate
    // a new nonzero
    //ierr = MatSetOption (_M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    //CHKERRABORT(this->comm(),ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption (_M_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    CHKERRABORT(this->comm(),ierr);
#endif

    //MatShift( _M_mat, 1 );
    //printMatlab( "shift.m" );
    //this->close();
    //this->zero();
    //this->close();

    //this->zeroEntriesDiagonal();
}


template <typename T>
void MatrixPetsc<T>::setIndexSplit(std::vector< std::vector<int> > const &indexSplit )
{
    this->M_IndexSplit=indexSplit;

    _M_petscIS.resize(indexSplit.size());

    int ierr=0;

    for (uint i = 0 ; i < indexSplit.size(); ++i)
        {
            PetscInt nDofForThisField = indexSplit[i].size();
            //std::cout << "\n setIndexSplit " << i << " ndof:" << nDofForThisField << "\n";

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            ierr = ISCreateGeneral(this->comm(),nDofForThisField,this->M_IndexSplit[i].data(),PETSC_COPY_VALUES,&_M_petscIS[i] );
#else
            ierr = ISCreateGeneral(this->comm(),nDofForThisField,this->M_IndexSplit[i].data(),&_M_petscIS[i] );
#endif
            CHKERRABORT(this->comm(),ierr);

#if 0
            ISView(_M_petscIS[i],PETSC_VIEWER_STDOUT_SELF);

            PetscInt n;
            /*
              Get the number of indices in the set
            */
            ISGetLocalSize(_M_petscIS[i],&n);
            std::cout << "Local size: " << n << "\n";
            const PetscInt *nindices;

            /*
              Get the indices in the index set
            */
            ISGetIndices(_M_petscIS[i],&nindices);
            for(int j = 0;j < n; ++j )
            {
                std::cout << nindices[j] << " ";
            }
            std::cout << "\n";
#endif
        }
    //std::cout << "\n setIndexSplit done\n";
}

template <typename T>
void MatrixPetsc<T>::updatePCFieldSplit(PC & pc)
{
#if 1
    int ierr=0;
    if (_M_mapPC.find(&pc)==_M_mapPC.end() )
        {
            _M_mapPC[&pc]=false;
        }

    if (_M_mapPC[&pc]==false)
        {
            const PCType pcName;
            ierr = PCGetType(pc,&pcName);
            CHKERRABORT(this->comm(),ierr);


            if ( std::string( PCFIELDSPLIT ) == std::string(pcName) )
                {
                    //std::cout << "\n updatePCFieldSplit " << _M_petscIS.size() << "\n";
                    _M_mapPC[&pc]=true;
                    for (uint i = 0 ; i < _M_petscIS.size(); ++i)
                        {
                            //std::cout << "\n split " << i << "\n";
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
                            //std::cout << "\n version >= 3.2 \n";
                            std::ostringstream os;
                            os << i;
                            ierr=PCFieldSplitSetIS(pc,os.str().c_str(),_M_petscIS[i]);
#else
                            //std::cout << "\n version < 3.2 \n";
                            ierr=PCFieldSplitSetIS(pc,_M_petscIS[i]);
#endif
                            //std::cout << "\n split " << i << "done\n" << std::endl;

                            CHKERRABORT(this->comm(),ierr);
                        }
                }
        }
#else
    int ierr=0;

    const PCType pcName;
    ierr = PCGetType(pc,&pcName);
    CHKERRABORT(this->comm(),ierr);

    if ( std::string( PCFIELDSPLIT ) == std::string(pcName) )
        {
            //std::cout << "\n updatePCFieldSplit \n";
            for (uint i = 0 ; i < _M_petscIS.size(); ++i)
                {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
                    ierr=PCFieldSplitSetIS(pc,PETSC_NULL,_M_petscIS[i]);
#else
                    ierr=PCFieldSplitSetIS(pc,_M_petscIS[i]);
#endif
                    CHKERRABORT(this->comm(),ierr);
                }
        }

#endif
}


template <typename T>
void MatrixPetsc<T>::zero ()
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    PetscBool is_assembled;
    MatAssembled( _M_mat, &is_assembled );
    if ( is_assembled )
    {
        ierr = MatZeroEntries(_M_mat);
        CHKERRABORT(this->comm(),ierr);

        //this->zeroEntriesDiagonal();
    }
    else
    {
        if ( this->graph() )
        {
            std::vector<PetscInt> cols( this->graph()->nCols(), 0 );
            std::vector<PetscScalar> v( this->graph()->nCols(), 0. );
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                PetscInt row = it->second.get<1>();
                std::copy( it->second.get<2>().begin(), it->second.get<2>().end(), cols.begin() );
                MatSetValues( _M_mat, 1, &row, it->second.get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
            }
        }
    }

}

template <typename T>
void MatrixPetsc<T>::zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{

    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    PetscBool is_assembled;
    MatAssembled( _M_mat, &is_assembled );
    if ( is_assembled )
    {
        ierr = MatZeroEntries(_M_mat);
        CHKERRABORT(this->comm(),ierr);

        //this->zeroEntriesDiagonal();
    }
    else
    {
        if ( this->graph() )
        {
            std::vector<PetscInt> cols( this->graph()->nCols(), 0 );
            std::vector<PetscScalar> v( this->graph()->nCols(), 0. );
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                PetscInt row = it->second.get<1>();
                std::copy( it->second.get<2>().begin(), it->second.get<2>().end(), cols.begin() );
                MatSetValues( _M_mat, 1, &row, it->second.get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
            }
        }
     }

}


template <typename T>
void MatrixPetsc<T>::clear ()
{
    int ierr=0;

    if ((this->isInitialized()) && (this->_M_destroy_mat_on_exit))
    {
        ierr = PETSc::MatDestroy (_M_mat);
        CHKERRABORT(this->comm(),ierr);

        for (int i=0;i< _M_petscIS.size(); ++i)
            {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
                ierr = ISDestroy(&_M_petscIS[i]);
                CHKERRABORT(this->comm(),ierr);
#else
                ierr = ISDestroy(_M_petscIS[i]);
                CHKERRABORT(this->comm(),ierr);
#endif
            }


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
    CHKERRABORT(this->comm(),ierr);
    ierr = MatAssemblyEnd   (_M_mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(this->comm(),ierr);
}



template <typename T>
inline
size_type MatrixPetsc<T>::size1 () const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize (_M_mat, &petsc_m, &petsc_n);

    return static_cast<size_type>(petsc_m);
}



template <typename T>
inline
size_type MatrixPetsc<T>::size2 () const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize (_M_mat, &petsc_m, &petsc_n);

    return static_cast<size_type>(petsc_n);
}



template <typename T>
inline
size_type MatrixPetsc<T>::rowStart () const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(this->comm(),ierr);

    return static_cast<size_type>(start);
}



template <typename T>
inline
size_type MatrixPetsc<T>::rowStop () const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(this->comm(),ierr);

    return static_cast<size_type>(stop);
}



template <typename T>
inline
void MatrixPetsc<T>::set (const size_type i,
                          const size_type j,
                          const value_type& value)
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>(value);
    ierr = MatSetValues(_M_mat, 1, &i_val, 1, &j_val,
                        &petsc_value, INSERT_VALUES);
    CHKERRABORT(this->comm(),ierr);
}



template <typename T>
inline
void MatrixPetsc<T>::add (const size_type i,
                          const size_type j,
                          const value_type& value)
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );

    int ierr=0, i_val=i, j_val=j;
    //Debug( 7013 ) << "[MatrixPetsc<>::add] adding value " << value << " at (" << i << "," << j << ")\n";
    PetscScalar petsc_value = static_cast<PetscScalar>(value);


    ierr = MatSetValues(_M_mat, 1, &i_val, 1, &j_val,
                        &petsc_value, ADD_VALUES);
    CHKERRABORT(this->comm(),ierr);
}
/*                                   */

template <typename T>
inline
bool MatrixPetsc<T>::closed() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );


    int ierr=0;
    PetscTruth assembled;

    ierr = MatAssembled(_M_mat, &assembled);
    CHKERRABORT(this->comm(),ierr);

    return (assembled == PETSC_TRUE);
}
template <typename T>
void
MatrixPetsc<T>::addMatrix(const ublas::matrix<value_type>& dm,
                          const std::vector<size_type>& rows,
                          const std::vector<size_type>& cols)
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    const size_type m = dm.size1();
    const size_type n = dm.size2();

    FEEL_ASSERT (rows.size() == size1()).error( "invalid row size" );
    FEEL_ASSERT (cols.size() == size2()).error( "invalid column size" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues(_M_mat,
                        m, (int*) boost::addressof( rows[0] ),
                        n, (int*) boost::addressof( cols[0] ),
                        (PetscScalar*) dm.data().begin(),
                        ADD_VALUES);
    CHKERRABORT(this->comm(),ierr);
}
template <typename T>
void
MatrixPetsc<T>::addMatrix ( int* rows, int nrows,
                            int* cols, int ncols,
                            value_type* data )
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues(_M_mat,
                        nrows, (int*) rows,
                        ncols, (int*) cols,
                        (PetscScalar*) data,
                        ADD_VALUES);
    CHKERRABORT(this->comm(),ierr);
}


#if 0
template <typename T>
void
MatrixPetsc<T>::printPython(const std::string name) const
{
    this->close();
    auto nbRow = this->size1();

    int ierr = 0;

    std::ofstream graphFile(name, std::ios::out /*| std::ios::app*/);


    for (size_type i = 0 ; i < nbRow ; ++i)
        {
            PetscInt row = i;
            PetscInt ncols;
            const PetscInt *cols;
            const PetscScalar *vals;

            PetscInt nrow2=1;
            PetscInt *row2 = new PetscInt[1]; row2[0] = i;

            ierr = MatGetRow(this->mat(), row, &ncols, &cols, &vals);
            CHKERRABORT(PETSC_COMM_WORLD,ierr);

            for (size_type jj = 0 ; jj < ncols ; ++jj)
                {
                    /*PetscInt *cols2 = new PetscInt[1];cols2[0]=cols[jj];
                    const PetscScalar *vals3;
                    ierr = MatGetValues(this->mat(), nrow2, row2, 1, cols2, vals3);
                    CHKERRABORT(PETSC_COMM_WORLD,ierr);
                    */

                    if(std::abs(vals[jj])>1.e-8 )graphFile << "[" << row <<","<< cols[jj] << "," << vals[jj] << "],";
                    //std::cout << "\n [updateBlockMat] i "<< start_i+row << " j " << start_j+cols[jj] << " val " << vals[jj] << std::endl;
                    //this->set(start_i+row,
                    //      start_j+cols[jj],
                    //      vals[jj]);
                }

            ierr = MatRestoreRow(this->mat(),row,&ncols,&cols,&vals);
            CHKERRABORT(PETSC_COMM_WORLD,ierr);
        }
    graphFile.close();
}
#endif

// print
template <typename T>
void
MatrixPetsc<T>::printMatlab (const std::string name) const
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" );

    // assert (this->closed());
    this->close();

    int ierr=0;
    PetscViewer petsc_viewer;


    ierr = PetscViewerCreate (this->comm(),
                              &petsc_viewer);
    CHKERRABORT(this->comm(),ierr);

    /**
     * Create an ASCII file containing the matrix
     * if a filename was provided.
     */
    if (name != "NULL")
    {
        ierr = PetscViewerASCIIOpen( this->comm(),
                                     name.c_str(),
                                     &petsc_viewer);
        CHKERRABORT(this->comm(),ierr);

        ierr = PetscViewerSetFormat (petsc_viewer,
                                     PETSC_VIEWER_ASCII_MATLAB);
                                     //PETSC_VIEWER_ASCII_PYTHON );
        CHKERRABORT(this->comm(),ierr);

        ierr = MatView (_M_mat, petsc_viewer);
        CHKERRABORT(this->comm(),ierr);
    }

    /**
     * Otherwise the matrix will be dumped to the screen.
     */
    else
    {
        //ierr = PetscViewerSetFormat (petsc_viewer,//PETSC_VIEWER_STDOUT_SELF,//PETSC_VIEWER_STDOUT_WORLD,
        //                             PETSC_VIEWER_ASCII_INFO_DETAIL);//PETSC_VIEWER_ASCII_MATLAB);
        //CHKERRABORT(this->comm(),ierr);

        //ierr = MatView (_M_mat,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_SELF);//
        //CHKERRABORT(this->comm(),ierr);
    }


    /**
     * Destroy the viewer.
     */
    ierr = PETSc::PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(this->comm(),ierr);
}

template <typename T>
inline
void
MatrixPetsc<T>::addMatrix (const T a_in, MatrixSparse<T> &X_in)
{
    FEEL_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    // sanity check. but this cannot avoid
    // crash due to incompatible sparsity structure...
    FEEL_ASSERT( this->size1() == X_in.size1())( this->size1() )( X_in.size1() ).error( "incompatible dimension" );
    FEEL_ASSERT( this->size2() == X_in.size2())( this->size2() )( X_in.size2() ).error( "incompatible dimension" );

    PetscScalar     a = static_cast<PetscScalar>      (a_in);
    MatrixPetsc<T>* X;
    if (this->comm().size()>1)
        {
            X = dynamic_cast<MatrixPetscMPI<T>*> (&X_in);
        }
    else
        {
            X = dynamic_cast<MatrixPetsc<T>*> (&X_in);
        }
    FEEL_ASSERT (X != 0).error( "invalid petsc matrix" );

    int ierr=0;

    // the matrix from which we copy the values has to be assembled/closed
    X->close ();

// 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = MatAXPY(&a,  X->_M_mat, _M_mat, (MatStructure)SAME_NONZERO_PATTERN);
    CHKERRABORT(this->comm(),ierr);

// 2.3.x & newer
#else
    if (this->comm().size()>1)
        {
            ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)DIFFERENT_NONZERO_PATTERN);
            //ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)SAME_NONZERO_PATTERN);
        }
    else
        {
            //this->close();
            //ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)SAME_NONZERO_PATTERN);
            //ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)SUBSET_NONZERO_PATTERN );
            ierr = MatAXPY(_M_mat, a, X->_M_mat, (MatStructure)DIFFERENT_NONZERO_PATTERN);
            //ierr = MatDuplicate(X->mat(),MAT_COPY_VALUES,&_M_mat);
        }
    CHKERRABORT(this->comm(),ierr);
#endif
}

template <typename T>
void
MatrixPetsc<T>::scale( T const a )
{
    int ierr = MatScale( _M_mat, a );
    CHKERRABORT(this->comm(),ierr);
}
template <typename T>
inline
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::l1Norm() const
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    int ierr=0;
    double petsc_value;
    real_type value;

    this->close();

    ierr = MatNorm(_M_mat, NORM_1, &petsc_value);
    CHKERRABORT(this->comm(),ierr);

    value = static_cast<real_type>(petsc_value);

    return value;
}

template <typename T>
inline
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::linftyNorm() const
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    int ierr=0;
    double petsc_value;
    real_type value;

    this->close();

    ierr = MatNorm(_M_mat, NORM_INFINITY, &petsc_value);
    CHKERRABORT(this->comm(),ierr);

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
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

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
    CHKERRABORT(this->comm(),ierr);

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
        CHKERRABORT(this->comm(),ierr);

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
    // the matrix needs to be closed for this to work
    this->close();

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption(_M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(_M_mat,MAT_KEEP_ZEROED_ROWS);
#endif
    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(this->comm(),ierr);

    if ( on_context.test( ON_ELIMINATION_KEEP_DIAGONAL ) )
        {
            VectorPetsc<value_type> diag( this->size1(), stop-start );
            MatGetDiagonal( _M_mat, diag.vec() );
            // in Petsc 3.2, we might want to look at the new interface so that
            // right hand side is automatically changed wrt to zeroing out the
            // matrix entries
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
            MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0,PETSC_NULL,PETSC_NULL);
#else
            MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0);
#endif
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



#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
            MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0,PETSC_NULL,PETSC_NULL);
#else
            MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0);
#endif
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
                            FEEL_ASSERT( found == true )( myRow )( cols[j] ).error ( "matrix graph is not symmetric" );
                        }
                    MatRestoreRow( _M_mat, myRow, &ncols, &cols, &vals );
                } // i
        }


}
template<typename T>
void
MatrixPetsc<T>::diagonal( Vector<value_type>& out ) const
{
    this->close();
    VectorPetsc<T>* v = dynamic_cast<VectorPetsc<T>*> (&out);
    int ierr = MatGetDiagonal( _M_mat, v->vec());
    CHKERRABORT(this->comm(),ierr);
}
template<typename T>
void
MatrixPetsc<T>::transpose( MatrixSparse<value_type>& Mt ) const
{
    this->close();

    MatrixPetsc<T>* Atrans = dynamic_cast<MatrixPetsc<T>*> (&Mt);

    int ierr = PETSc::MatDestroy( Atrans->_M_mat );
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( _M_mat, MAT_INITIAL_MATRIX,&Atrans->_M_mat );
#else
    ierr = MatTranspose( _M_mat, &Atrans->_M_mat );
#endif

    CHKERRABORT(this->comm(),ierr);

    if (this->size1()==this->size2())
        {
            PetscTruth isSymmetric;
            MatEqual( _M_mat, Atrans->_M_mat, &isSymmetric);
            if (isSymmetric) {
#if (PETSC_VERSION_MAJOR >= 3)
                MatSetOption(_M_mat,MAT_SYMMETRIC,PETSC_TRUE);
#else
                MatSetOption(_M_mat,MAT_SYMMETRIC );
#endif
            } else {
                Debug(7013) << "[MatrixPETSc::transpose] Petsc matrix is non-symmetric \n";
            }
        }
    Mt.setGraph(this->graph()->transpose());

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
    PETSc::MatDestroy( Atrans );

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
    CHKERRABORT(this->comm(),ierr);

    ierr = MatDuplicate(_M_mat,MAT_COPY_VALUES,&A[1]);
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( A[1], MAT_INITIAL_MATRIX, &A[1] );
#else
    MatTranspose( A[1], &A[1] );
#endif
    CHKERRABORT(this->comm(),ierr);

    MatrixPetsc<T>* B = dynamic_cast<MatrixPetsc<T>*> (&Mt);
    FEEL_ASSERT( B ).error ( "[symmetricPart] invalid dynamic_cast<MatrixPetsc>" );

    ierr = MatCreateComposite(PETSC_COMM_WORLD,2,A,&B->_M_mat);
    CHKERRABORT(this->comm(),ierr);


//A Completer pour PETSC < 3 !!!!!!!!
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatCompositeSetType(B->_M_mat,MAT_COMPOSITE_ADDITIVE);
#endif

    CHKERRABORT(this->comm(),ierr);

    ierr = MatCompositeMerge(B->_M_mat);
    CHKERRABORT(this->comm(),ierr);

    ierr = MatScale( B->_M_mat, 0.5 );
    CHKERRABORT(this->comm(),ierr);


    // check that B is now symmetric
    Mat Btrans;
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( B->_M_mat, MAT_INITIAL_MATRIX, &Btrans );
#else
    ierr = MatTranspose( B->_M_mat, &Btrans );
#endif

    CHKERRABORT(this->comm(),ierr);

    ierr = MatEqual( B->_M_mat, Btrans, &isSymmetric);
    CHKERRABORT(this->comm(),ierr);

    if (isSymmetric) {
        Log() << "[MatrixPetsc<T>::symmetricPart] symmetric part is really symmetric\n";
#if (PETSC_VERSION_MAJOR >= 3)
        ierr = MatSetOption(B->_M_mat,MAT_SYMMETRIC,PETSC_TRUE);
#else
	ierr = MatSetOption(B->_M_mat,MAT_SYMMETRIC);
#endif
        CHKERRABORT(this->comm(),ierr);

    } else {
        PetscPrintf(PETSC_COMM_WORLD,"Warning: Petsc matrix is non-symmetric \n");
    }
    ierr = PETSc::MatDestroy (Btrans);
    CHKERRABORT(this->comm(),ierr);
}


template<typename T>
typename MatrixPetsc<T>::value_type
MatrixPetsc<T>::energy( Vector<value_type> const& __v,
                        Vector<value_type> const& __u,
                        bool transpose ) const
{
    this->close();

    PetscScalar e;
    if ( dynamic_cast<VectorPetsc<T> const*>( &__v ) != (VectorPetsc<T> const*)0 )
    {
        VectorPetsc<T> const& v   = dynamic_cast<VectorPetsc<T> const&>( __v );
        VectorPetsc<T> const& u   = dynamic_cast<VectorPetsc<T> const&>( __u );
        VectorPetsc<value_type> z( __v.size(), __v.localSize() );
        if ( !transpose )
            MatMult( _M_mat, u.vec(), z.vec() );
        else
            MatMultTranspose( _M_mat, u.vec(), z.vec() );
        VecDot( v.vec(), z.vec(), &e );
    }
    else
    {
        VectorPetsc<value_type> u( __u.size(), __u.localSize() );
        {
            size_type s = u.localSize();
            size_type start = u.firstLocalIndex();
            for( size_type i = 0; i < s; ++i )
                u.set( start + i, __u( start + i ) );
        }
        VectorPetsc<value_type> v( __v.size(), __v.localSize() );
        {
            size_type s = v.localSize();
            size_type start = v.firstLocalIndex();
            for( size_type i = 0; i < s; ++i )
                v.set( start + i, __v( start + i ) );
        }
        VectorPetsc<value_type> z( __v.size(), __v.localSize() );
        if ( !transpose )
            MatMult( _M_mat, u.vec(), z.vec() );
        else
            MatMultTranspose( _M_mat, u.vec(), z.vec() );

        VecDot( v.vec(), z.vec(), &e );
    }
    return e;
}

template<typename T>
void
MatrixPetsc<T>::updateBlockMat(boost::shared_ptr<MatrixSparse<T> > m, size_type start_i, size_type start_j)
{
#if 1
        auto blockMatrix = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>(&*(m) ) );

        auto nbRowInBlock = blockMatrix->size1();

        int ierr = 0;
#if 1
        for (size_type ii = 0 ; ii < nbRowInBlock ; ++ii)
            {

                PetscInt row = ii;
                PetscInt ncols;
                const PetscInt *cols;
                const PetscScalar *vals;

                ierr = MatGetRow(blockMatrix->mat(), row, &ncols, &cols, &vals);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);

                for (size_type jj = 0 ; jj < ncols ; ++jj)
                    {
                        //std::cout << "\n [updateBlockMat] i "<< start_i+row << " j " << start_j+cols[jj] << " val " << vals[jj] << std::endl;
                        this->set(start_i+row,
                                  start_j+cols[jj],
                                  vals[jj]);
                    }

                ierr = MatRestoreRow(blockMatrix->mat(),row,&ncols,&cols,&vals);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);

            }

#else

                int row = ii;
                int ncols;
                const int *cols;
                const double *vals;

                ierr = MatGetRow(blockMatrix->mat(), row, &ncols, &cols, &vals);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);

                //PetscScalar petsc_value = static_cast<PetscScalar>(value);

                int *aaa = new int[nbRowInBlock];
                for (size_type ii = 0 ; ii < nbRowInBlock ; ++ii) aaa[start_i+ii];
                int bbb=nbRowInBlock;

                //int *aaa = new int[1];aaa[0]=start_i+row;
                //int bbb=1;

                //this->addMatrix ( aaa, 1,
                //                  cols, ncols,
                //                  vals );
                ierr = MatSetValues(_M_mat,
                                    bbb, aaa/*(int*) rows*/,
                                    ncols, /*(int*)*/ cols,
                                    /*(PetscScalar*) data*/vals,
                                    INSERT_VALUES/*ADD_VALUES*/);

                ierr = MatRestoreRow(blockMatrix->mat(),row,&ncols,&cols,&vals);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);


#endif

#if 0




        PetscInt shift=0;
        PetscBool symmetric = PETSC_FALSE;
        PetscBool inodecompressed = PETSC_FALSE;
        PetscInt n;
        PetscInt* ia;
        PetscInt* ja;
        PetscBool done;
        ierr = MatGetRowIJ(blockMatrix->mat(),
                           shift,
                           symmetric,
                           inodecompressed,
                           &n,&ia,&ja,
                           &done);


        this->addMatrix ( int* rows, int nrows,
                          int* cols, int ncols,
                          value_type* data )

                    this->set(start_i+row,
                              start_j+cols[jj],
                              vals[jj]);


        ierr = MatRestoreRowIJ(blockMatrix->mat(),
                               shift,
                               symmetric,
                               inodecompressed,
                               &n,&ia,&ja,
                               &done );


#endif

#endif





}

template <typename T>
inline
void MatrixPetsc<T>::zeroEntriesDiagonal()
{
#if 1
#if 0 //trop lent!
    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    CHKERRABORT(this->comm(),ierr);

    //VectorPetsc<value_type> diag( std::min(this->size1(),this->size2()), stop-start );
    //VectorPetsc<value_type> diag( std::max(this->size1(),this->size2()), stop-start );
    VectorPetsc<value_type> diag( this->size1(), stop-start );
    //VectorPetsc<value_type> diag( this->size1(), this->size1() );

    //MatGetDiagonal( _M_mat, diag.vec() );
    //diag.zero();

    //std:: cout << diag.size();
    //MatZeroRows( _M_mat, rows.size(), rows.data(), 1.0);
    ierr = MatDiagonalSet( _M_mat, diag.vec(), INSERT_VALUES );
    CHKERRABORT(this->comm(),ierr);
#else
    //std::cout << "\n this->size1(),this->size2() " << this->size1() << " " << this->size2() << std::endl;
#if 0
    for (uint i = 0;i <std::min(this->size1(),this->size2());++i)
        this->set(i,i,0.);
#else
    //std::cout << "zeroEntriesDiagonal()"<< std::endl;
    for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
        {
            size_type index = it->first;std::cout << "\n index " << index;
            if (index < std::min(this->size1(),this->size2()))
                {
                    //if (!it->second.get<2>().empty())
                    if (it->second.get<2>().find(index)!=it->second.get<2>().end()  )
                        this->set(index,index,0.);
                }
        }
#endif
#endif
#endif
}


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI()
    :
    super()
{}

//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI(DataMap const& dmRow, DataMap const& dmCol, WorldComm const& worldComm )
    :
    super(dmRow,dmCol,worldComm)
{}

//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI(Mat m, DataMap const& dmRow, DataMap const& dmCol)
    :
    super(m)
{
    this->setMapRow(dmRow);
    this->setMapCol(dmCol);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void MatrixPetscMPI<T>::init( const size_type m,
                              const size_type n,
                              const size_type m_l,
                              const size_type n_l,
                              const size_type nnz,
                              const size_type /*noz*/)
{
    //std::cout << "\n MatrixPetscMPI<T>::init without graph " << std::endl;
}

template <typename T>
void MatrixPetscMPI<T>::init( const size_type m,
                              const size_type n,
                              const size_type m_l,
                              const size_type n_l,
                              graph_ptrtype const& graph )
{
    //std::cout << "\n MatrixPetscMPI<T>::init with graph " << std::endl;

    this->setGraph( graph );

    // Clear initialized matrices
    if (this->isInitialized())
        this->clear();

    this->setInitialized(  true );

    if (m==0)
        return;

    int ierr     = 0;
    int m_global = static_cast<int>(m);
    int n_global = static_cast<int>(n);
    int m_local  = static_cast<int>(m_l);
    int n_local  = static_cast<int>(n_l);


    PetscInt *dnz;
    PetscInt n_dnz = this->graph()->nNzOnProc().size();
    //std::cout << "\n n_ndz = " << n_dnz << std::endl;
    dnz = new PetscInt[n_dnz];
    std::copy( this->graph()->nNzOnProc().begin(),
               this->graph()->nNzOnProc().end(),
               dnz );


    PetscInt *dnzOffProc;
    PetscInt n_dnzOffProc = this->graph()->nNzOffProc().size();
    dnzOffProc = new PetscInt[ n_dnzOffProc ];
    std::copy( this->graph()->nNzOffProc().begin(),
               this->graph()->nNzOffProc().end(),
               dnzOffProc );

    if (n_dnzOffProc==0)
        dnzOffProc = PETSC_NULL;

    //  Mat _M_matttt;
#if 1
    ierr = MatCreateMPIAIJ (this->comm(),
                            m_local, n_local,
                            m_global, n_global,
                            /*PETSC_DECIDE*//*n_dnz*/0, /*PETSC_NULL*/dnz,
                            /*PETSC_DECIDE*/0/*n_dnzOffProc*/, dnzOffProc,
                            //&_M_matttt);
                            &(this->mat()) );//(&this->_M_mat));
#else
    ierr = MatCreateMPIAIJ (this->comm(),
                            m_local, n_local,
                            m_global, n_global,
                            /*PETSC_DECIDE*//*n_dnz*/0, PETSC_NULL,
                            /*PETSC_DECIDE*/0/*n_dnzOffProc*/,PETSC_NULL,
                            //&_M_matttt);
                            &(this->mat()) );//(&this->_M_mat));
#endif

//ierr = MatSetFromOptions(_M_matttt);//this->mat());
//    CHKERRABORT(this->comm(),ierr);

    CHKERRABORT(this->comm(),ierr);

    this->graph()->close();
    std::vector<PetscInt> ia( this->graph()->ia().size() );
    std::vector<PetscInt> ja( this->graph()->ja().size() );
    std::copy( this->graph()->ia().begin(), this->graph()->ia().end(), ia.begin() );
    std::copy( this->graph()->ja().begin(), this->graph()->ja().end(), ja.begin() );
#if 0
    ierr = MatMPIAIJSetPreallocation(this->mat(), 0, dnz, 0, dnzOffProc);
    //ierr = MatMPIAIJSetPreallocation(this->mat(), 0/*n_dnz*/ , PETSC_NULL, /*n_dnzOffProc*/0, PETSC_NULL);
#else
    ierr = MatMPIAIJSetPreallocationCSR(this->mat(), ia.data() , ja.data(), this->graph()->a().data() );
#endif
    CHKERRABORT(this->comm(),ierr);
    //----------------------------------------------------------------------------------//
    // localToGlobalMapping
    IS isRow;
    IS isCol;
    ISLocalToGlobalMapping isLocToGlobMapRow;
    ISLocalToGlobalMapping isLocToGlobMapCol;
#if 0
    auto idxRow = this->mapRow().mapGlobalProcessToGlobalCluster();
    auto idxCol = this->mapCol().mapGlobalProcessToGlobalCluster();
    //std::cout << "idxRow.size()" << idxRow.size() << std::endl;
    ierr = ISCreateGeneral(this->comm(), idxRow.size(), &idxRow[0], PETSC_COPY_VALUES, &isRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISCreateGeneral(this->comm(), idxCol.size(), &idxCol[0], PETSC_COPY_VALUES,&isCol);
    CHKERRABORT(this->comm(),ierr);
#else
    PetscInt *idxRow;
    PetscInt *idxCol;
    PetscInt n_idxRow =  this->mapRow().mapGlobalProcessToGlobalCluster().size();
    PetscInt n_idxCol =  this->mapCol().mapGlobalProcessToGlobalCluster().size();
    idxRow = new PetscInt[n_idxRow];
    idxCol = new PetscInt[n_idxCol];
    std::copy( this->mapRow().mapGlobalProcessToGlobalCluster().begin(),
               this->mapRow().mapGlobalProcessToGlobalCluster().end(),
               idxRow );
    std::copy( this->mapCol().mapGlobalProcessToGlobalCluster().begin(),
               this->mapCol().mapGlobalProcessToGlobalCluster().end(),
               idxCol );

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(), n_idxRow, idxRow, PETSC_COPY_VALUES, &isRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISCreateGeneral(this->comm(), n_idxCol, idxCol, PETSC_COPY_VALUES, &isCol);
    CHKERRABORT(this->comm(),ierr);
#else
    ierr = ISCreateGeneral(this->comm(), n_idxRow, idxRow, &isRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISCreateGeneral(this->comm(), n_idxCol, idxCol, &isCol);
    CHKERRABORT(this->comm(),ierr);

#endif
#endif
    ierr=ISLocalToGlobalMappingCreateIS(isRow, &isLocToGlobMapRow);
    CHKERRABORT(this->comm(),ierr);

    ierr=ISLocalToGlobalMappingCreateIS(isCol, &isLocToGlobMapCol);
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = MatSetLocalToGlobalMapping(this->mat(),isLocToGlobMapRow,isLocToGlobMapCol);
#else
    ierr = MatSetLocalToGlobalMapping(this->mat(),isLocToGlobMapRow);
#endif
    CHKERRABORT(this->comm(),ierr);

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy(&isRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISDestroy(&isCol);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISLocalToGlobalMappingDestroy(&isLocToGlobMapRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISLocalToGlobalMappingDestroy(&isLocToGlobMapCol);
    CHKERRABORT(this->comm(),ierr);
#else
    ierr = ISDestroy(isRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISDestroy(isCol);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISLocalToGlobalMappingDestroy(isLocToGlobMapRow);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISLocalToGlobalMappingDestroy(isLocToGlobMapCol);
    CHKERRABORT(this->comm(),ierr);
#endif
    delete idxRow;
    delete idxCol;

    //----------------------------------------------------------------------------------//
    // options
    //    Mat _M_matttt;

    ierr = MatSetFromOptions(this->mat());
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    //ATTENTION remetre
    //MatSetOption(this->mat(),MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
    MatSetOption(this->mat(),MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption(this->mat(),MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(this->mat(),MAT_KEEP_ZEROED_ROWS);
#endif



#if 0
    // additional insertions will not be allowed if they generate
    // a new nonzero
    //ierr = MatSetOption (_M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    //CHKERRABORT(this->comm(),ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption (this->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    CHKERRABORT(this->comm(),ierr);
#endif

    //this->zero();
    //this->zeroEntriesDiagonal();
    //this->printMatlab("NULL");
    //std::cout << "finish init with graph " << std::endl;

}


//----------------------------------------------------------------------------------------------------//



template <typename T>
inline
size_type MatrixPetscMPI<T>::size1() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    //ierr = MatGetLocalSize (this->mat(), &petsc_m, &petsc_n);
    //CHKERRABORT(this->comm(),ierr);
    petsc_m = this->mapRow().nLocalDofWithGhost();

    return static_cast<size_type>(petsc_m);
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
size_type MatrixPetscMPI<T>::size2() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    //ierr = MatGetLocalSize (this->mat(), &petsc_m, &petsc_n);
    //CHKERRABORT(this->comm(),ierr);
    petsc_n = this->mapCol().nLocalDofWithGhost();

    return static_cast<size_type>(petsc_n);
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
size_type MatrixPetscMPI<T>::rowStart() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    //ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    //CHKERRABORT(this->comm(),ierr);
    start=0; stop=this->mapRow().nLocalDofWithGhost();
    return static_cast<size_type>(start);
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
size_type MatrixPetscMPI<T>::rowStop() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    //ierr = MatGetOwnershipRange(_M_mat, &start, &stop);
    //CHKERRABORT(this->comm(),ierr);
    start=0; stop=this->mapRow().nLocalDofWithGhost();
    return static_cast<size_type>(stop);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
size_type MatrixPetscMPI<T>::colStart() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0;
    start=0; stop=this->mapCol().nLocalDofWithGhost();
    return static_cast<size_type>(start);
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
size_type MatrixPetscMPI<T>::colStop() const
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0;
    start=0; stop=this->mapCol().nLocalDofWithGhost();
    return static_cast<size_type>(stop);
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
void MatrixPetscMPI<T>::set(const size_type i,
                            const size_type j,
                            const value_type& value)
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );;

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>(value);
    ierr = MatSetValuesLocal(this->mat(), 1, &i_val, 1, &j_val,
                             &petsc_value, INSERT_VALUES);

    CHKERRABORT(this->comm(),ierr);
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
void MatrixPetscMPI<T>::add (const size_type i,
                             const size_type j,
                             const value_type& value)
{
    FEEL_ASSERT (this->isInitialized()).error( "MatrixPetsc<> not properly initialized" );

    int ierr=0, i_val=i, j_val=j;
    //Debug( 7013 ) << "[MatrixPetsc<>::add] adding value " << value << " at (" << i << "," << j << ")\n";
    PetscScalar petsc_value = static_cast<PetscScalar>(value);

    ierr = MatSetValuesLocal(this->mat(), 1, &i_val, 1, &j_val,
                             &petsc_value, ADD_VALUES);
    CHKERRABORT(this->comm(),ierr);
}


//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::addMatrix(const ublas::matrix<value_type>& dm,
                             const std::vector<size_type>& rows,
                             const std::vector<size_type>& cols)
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );

    const size_type m = dm.size1();
    const size_type n = dm.size2();

    FEEL_ASSERT (rows.size() == size1()).error( "invalid row size" );
    FEEL_ASSERT (cols.size() == size2()).error( "invalid column size" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValuesLocal(this->mat(),
                             m, (int*) boost::addressof( rows[0] ),
                             n, (int*) boost::addressof( cols[0] ),
                             (PetscScalar*) dm.data().begin(),
                             ADD_VALUES);


    CHKERRABORT(this->comm(),ierr);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::addMatrix( int* rows, int nrows,
                              int* cols, int ncols,
                              value_type* data )
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not initialized" );
    /*for (int k=0;k<nrows;++k)
        for (int q=0;q<ncols;++q)
            { std::cout << "rank " << this->comm().rank()
                        << "add mat("<< rows[k] <<"," << cols[q]
                        << ") = "//<< data[k][q]
                        << " " << this->mapCol().mapGlobalProcessToGlobalCluster()[cols[q]]
                        << std::endl; }*/
    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValuesLocal(this->mat(),
                             nrows, (int*) rows,
                             ncols, (int*) cols,
                             (PetscScalar*) data,
                             ADD_VALUES);

    CHKERRABORT(this->comm(),ierr);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::zero()
{
    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    PetscBool is_assembled;
    MatAssembled( this->mat(), &is_assembled );
    if( is_assembled )
    {
        //std::cout << "MPI is_assembled " << std::endl;
        ierr = MatZeroEntries(this->mat());
        CHKERRABORT(this->comm(),ierr);

        //this->zeroEntriesDiagonal();
    }
    else
    {
        if ( this->graph() )
        {
#if 1
            //std::vector<PetscInt> cols( this->graph()->nCols(), 0 );
            //std::vector<PetscScalar> v( this->graph()->nCols(), 0. );
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                if ( it->second.get<0>() == this->comm().rank() )
                    {

                        std::vector<PetscInt> cols(  it->second.get<2>().size(), 0 );
                        //std::set<PetscInt> cols;

                        //PetscInt row = it->second.get<1>();
                        PetscInt row =  it->first;


                        //MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_NULL, PETSC_NULL);

                        //this->mapRow()->firstDofGlobalCluster( this->comm().rank() );
                        auto it2=it->second.get<2>().begin();
                        auto  en2=it->second.get<2>().end();
                        size_type cpt=0;
                        for ( ; it2!=en2 ; ++it2)
                            {
                                if ( (*it2 >= this->graph()->firstColEntryOnProc() ) &&
                                     (*it2 <= this->graph()->lastColEntryOnProc() ) )
                                    {
                                        cols[cpt] = this->mapRow().mapGlobalClusterToGlobalProcess()[*it2];
                                        ++cpt;
                                    }
                            }
                        cols.resize(cpt);
                        std::vector<PetscScalar> v(  cpt,0 );

                        //std::copy( it->second.get<2>().begin(), it->second.get<2>().end(), cols.begin() );
                        MatSetValuesLocal( this->mat(), 1, &row, it->second.get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
                    }
            }
#endif
        }
    }

}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::zero( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{

    FEEL_ASSERT (this->isInitialized()).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    PetscBool is_assembled;
    MatAssembled( this->mat(), &is_assembled );
    if ( is_assembled )
    {
        ierr = MatZeroEntries(this->mat());
        CHKERRABORT(this->comm(),ierr);

        //this->zeroEntriesDiagonal();
    }
    else
    {
        if ( this->graph() )
        {
#if 0
            std::vector<PetscInt> cols( this->graph()->nCols(), 0 );
            std::vector<PetscScalar> v( this->graph()->nCols(), 0. );
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                PetscInt row = it->second.get<1>();
                std::copy( it->second.get<2>().begin(), it->second.get<2>().end(), cols.begin() );
                MatSetValuesLocal( this->mat(), 1, &row, it->second.get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
            }
#else
            //std::vector<PetscInt> cols( this->graph()->nCols(), 0 );
            //std::vector<PetscScalar> v( this->graph()->nCols(), 0. );
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                if ( it->second.get<0>() == this->comm().rank() )
                    {

                        std::vector<PetscInt> cols(  it->second.get<2>().size(), 0 );
                        //std::set<PetscInt> cols;

                        //PetscInt row = it->second.get<1>();
                        PetscInt row =  it->first;


                        //MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_NULL, PETSC_NULL);

                        //this->mapRow()->firstDofGlobalCluster( this->comm().rank() );
                        auto it2=it->second.get<2>().begin();
                        auto  en2=it->second.get<2>().end();
                        size_type cpt=0;
                        for ( ; it2!=en2 ; ++it2)
                            {
                                if ( (*it2 >= this->graph()->firstColEntryOnProc() ) &&
                                     (*it2 <= this->graph()->lastColEntryOnProc() ) )
                                    {
                                        cols[cpt] = this->mapRow().mapGlobalClusterToGlobalProcess()[*it2];
                                        ++cpt;
                                    }
                            }
                        cols.resize(cpt);
                        std::vector<PetscScalar> v(  cpt,0 );

                        //std::copy( it->second.get<2>().begin(), it->second.get<2>().end(), cols.begin() );
                        MatSetValuesLocal( this->mat(), 1, &row, it->second.get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
                    }
            }

#endif
        }
     }

}

//----------------------------------------------------------------------------------------------------//

template<typename T>
void
MatrixPetscMPI<T>::zeroRows( std::vector<int> const& rows,
                             std::vector<value_type> const& values,
                             Vector<value_type>& rhs,
                             Context const& on_context )
{
    // the matrix doesn't be closed because not all processors are present here with composite spaces(this call must be done after)
    // this->close();

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption(this->mat(),MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
    MatSetOption(this->mat(),MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption(this->mat(),MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
#else
    MatSetOption(this->mat(),MAT_KEEP_ZEROED_ROWS);
#endif

    int start=0, stop=this->mapRow().nLocalDofWithGhost(), ierr=0;
    //ierr = MatGetOwnershipRange(_M_mat, &start, &stop);

    if (false)// on_context.test( ON_ELIMINATION_KEEP_DIAGONAL ) )
        {
            VectorPetscMPI<value_type> diag( this->mapRow() );

            //VectorPetsc<value_type> diag( this->mapRow().nLocalDofWithoutGhost(),this->mapRow().worldComm() );
            //diag( this->mapRow().nLocalDofWithGhost(),this->mapRow().worldComm().subWorldComm(this->mapRow().worldComm().mapColorWorld()[this->mapRow().worldComm().globalRank()  ] ));
            ierr =MatGetDiagonal( this->mat(), diag.vec() );
            CHKERRABORT(this->comm(),ierr);

            // in Petsc 3.2, we might want to look at the new interface so that
            // right hand side is automatically changed wrt to zeroing out the
            // matrix entries
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
            ierr = MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_NULL, PETSC_NULL);
            //CHKERRABORT(this->comm(),ierr);
#else
            MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0);
#endif
            // doesn't work with composite space
            ierr=MatDiagonalSet( this->mat(), diag.vec(), INSERT_VALUES );
            //CHKERRABORT(this->comm(),ierr);

            // important close
            diag.close();

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
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
            MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_NULL, PETSC_NULL);
            //CHKERRABORT(this->comm(),ierr);
#else
            MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0);
#endif
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

    // rsh doesn't be closed because not all processors are present here with composite spaces(this call must be done after)
    // rhs.close();

    //reset MatOption (assemble with communication)
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption(this->mat(),MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_FALSE);
#else
    // ???
#endif
} // zeroRows

//----------------------------------------------------------------------------------------------------//

template class MatrixPetsc<double>;
template class MatrixPetscMPI<double>;
} // Feel

#endif // HAVE_PETSC_H
