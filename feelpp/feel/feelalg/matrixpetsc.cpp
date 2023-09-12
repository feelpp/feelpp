/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-01-03
 */

#include <boost/smart_ptr/make_shared.hpp>

#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/functionspetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feeltiming/tic.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/string.hpp> // Needed to send/receive strings!

#if defined( FEELPP_HAS_PETSC_H )

extern "C"
{
#include "petscsys.h"
}



namespace Feel
{

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    M_destroy_mat_on_exit( true )
{}
template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol )
    :
    super( dmRow,dmCol ),
    M_destroy_mat_on_exit( true )
{}
template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, worldcomm_ptr_t const& worldComm )
    :
    super( dmRow,dmCol,worldComm ),
    M_destroy_mat_on_exit( true )
{}

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( MatrixSparse<value_type> const& M, IS& isrow, IS& iscol )
    :
    super(),
    M_destroy_mat_on_exit( true )
{
    MatrixPetsc<T> const* A = dynamic_cast<MatrixPetsc<T> const*> ( &M );
    int ierr=0;
    PetscInt nrow;
    PetscInt ncol;
    ISGetSize(isrow,&nrow);
    ISGetSize(iscol,&ncol);
    datamap_ptrtype dmrow( new datamap_type(nrow, nrow) );
    datamap_ptrtype dmcol( new datamap_type(ncol, ncol) );
    this->setMapRow(dmrow);
    this->setMapCol(dmcol);
#if PETSC_VERSION_LESS_THAN( 3,8,0 )
    ierr = MatGetSubMatrix(A->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &this->M_mat);
#else
    ierr = MatCreateSubMatrix(A->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &this->M_mat);
#endif
    CHKERRABORT( this->comm(),ierr );
    this->setInitialized( true );
    this->close();
}

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( MatrixSparse<value_type> const& M, std::vector<int> const& rowIndex, std::vector<int> const& colIndex )
    :
    super(),
    M_destroy_mat_on_exit( true )
{
    MatrixPetsc<T> const* A = dynamic_cast<MatrixPetsc<T> const*> ( &M );
    int ierr=0;
    IS isrow;
    IS iscol;
    PetscInt *rowMap;
    PetscInt *colMap;
    int nrow = rowIndex.size();
    int ncol = colIndex.size();

    PetscMalloc(nrow*sizeof(PetscInt),&rowMap);
    PetscMalloc(ncol*sizeof(PetscInt),&colMap);

    for (int i=0; i<nrow; i++) rowMap[i] = rowIndex[i];
    for (int i=0; i<ncol; i++) colMap[i] = colIndex[i];

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(Environment::worldComm(),nrow,rowMap,PETSC_COPY_VALUES,&isrow);
    CHKERRABORT( this->comm(),ierr );
    ierr = ISCreateGeneral(Environment::worldComm(),ncol,colMap,PETSC_COPY_VALUES,&iscol);
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISCreateGeneral(Environment::worldComm(),nrow,rowMap,&isrow);
    CHKERRABORT( this->comm(),ierr );
    ierr = ISCreateGeneral(Environment::worldComm(),ncol,colMap,&iscol);
    CHKERRABORT( this->comm(),ierr );
#endif
    PetscFree(rowMap);
    PetscFree(colMap);

    datamap_ptrtype dmrow( new datamap_type(nrow, nrow) );
    datamap_ptrtype dmcol( new datamap_type(ncol, ncol) );
    this->setMapRow(dmrow);
    this->setMapCol(dmcol);
#if PETSC_VERSION_LESS_THAN( 3,8,0 )
    ierr = MatGetSubMatrix(A->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &this->M_mat);
#else
    ierr = MatCreateSubMatrix(A->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &this->M_mat);
#endif
    CHKERRABORT( this->comm(),ierr );
    this->setInitialized( true );
    this->close();
}


template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( Mat m )
    :
    M_destroy_mat_on_exit( false )
{
    this->M_mat = m;
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS );
#endif
    this->setInitialized( true );
    this->setIsClosed( true );
}

template <typename T>
inline
MatrixPetsc<T>::MatrixPetsc( Mat m, datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, bool destroyMatOnExit )
    :
    super( dmRow,dmCol ),
    M_destroy_mat_on_exit( destroyMatOnExit )
{
    this->M_mat = m;
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS );
#endif
    this->setInitialized( true );
    this->setIsClosed( true );
}




template <typename T>
inline
MatrixPetsc<T>::~MatrixPetsc()
{
    this->clear();
}


template <typename T>
typename MatrixPetsc<T>::clone_ptrtype
MatrixPetsc<T>::clone () const
{
    clone_ptrtype cloned_matrix;
    if ( auto m = dynamic_cast<MatrixPetscMPI<T> const*>( this ); m )
    {
        Mat M;
        MatDuplicate( m->mat(), MAT_COPY_VALUES, &M );
        cloned_matrix.reset( new MatrixPetscMPI<T>( M, this->mapRowPtr(), this->mapColPtr(), false, true ) );
    }  
    else if ( auto m = dynamic_cast<MatrixPetsc<T> const*>( this ); m )
    {
        Mat M;
        MatDuplicate( m->mat(), MAT_COPY_VALUES, &M );
        cloned_matrix.reset( new MatrixPetsc<T>( M, this->mapRowPtr(), this->mapColPtr(), true ) );
    }
    return cloned_matrix;
}

template <typename T>
void MatrixPetsc<T>::init ( const size_type m,
                            const size_type n,
                            const size_type m_l,
                            const size_type n_l,
                            const size_type nnz,
                            const size_type /*noz*/ )
{
    //std::cout << "\nSEQUENTIAL : init without graph"<<std::endl;

    if ( ( m==0 ) || ( n==0 ) )
        return;

    {
        // Clear initialized matrices
        if ( this->isInitialized() )
            this->clear();

        this->setInitialized( true );
    }


    int ierr     = 0;
    int m_global = static_cast<int>( m );
    int n_global = static_cast<int>( n );
    int m_local  = static_cast<int>( m_l );
    int n_local  = static_cast<int>( n_l );
    int n_nz     = static_cast<int>( nnz );
    //int n_oz     = static_cast<int>(noz);

    // create a sequential matrix on one processor
    if ( ( m_l == m ) && ( n_l == n ) )
    {
        // Create matrix.  Revisit later to do preallocation and make more efficient
        ierr = MatCreateSeqAIJ ( this->comm(), m_global, n_global,
                                 n_nz, PETSC_IGNORE, &M_mat );
        CHKERRABORT( this->comm(),ierr );



    }

    else
    {
#if PETSC_VERSION_LESS_THAN(3,3,0)
        ierr = MatCreateMPIAIJ ( this->comm(), m_local, n_local, m_global, n_global,
                                 PETSC_DECIDE, PETSC_IGNORE, PETSC_DECIDE, PETSC_IGNORE, &M_mat );
#else
        ierr = MatCreateAIJ ( this->comm(), m_local, n_local, m_global, n_global,
                                 PETSC_DECIDE, PETSC_IGNORE, PETSC_DECIDE, PETSC_IGNORE, &M_mat );
#endif
        //ierr = MatCreateMPIAIJ (this->comm(), m_local, n_local, m_global, n_global,
        ///n_nz, PETSC_IGNORE, n_oz, PETSC_IGNORE, &M_mat);
        //MatCreate(this->comm(),m_local,n_local,m_global,n_global, &M_mat);
        //MatCreate(this->comm(),PETSC_DECIDE,PETSC_DECIDE,m_global,n_global, &M_mat);
        //MatSetSizes(M_mat,m_local,n_local,m_global,n_global);
        CHKERRABORT( this->comm(),ierr );

    }

    ierr = MatSetFromOptions ( M_mat );
    CHKERRABORT( this->comm(),ierr );

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    ierr = MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
    MatSetOption( M_mat,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    ierr = MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    ierr = MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS );
#endif
    CHKERRABORT( this->comm(),ierr );
#if 1
    // additional insertions will not be allowed if they generate
    // a new nonzero
    //ierr = MatSetOption (M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    //CHKERRABORT(this->comm(),ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption ( M_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
    CHKERRABORT( this->comm(),ierr );
#endif // 0

    //this->zero ();
    //this->zeroEntriesDiagonal();
}

template <typename T>
void MatrixPetsc<T>::init ( const size_type /*m*/,
                            const size_type /*n*/,
                            const size_type /*m_l*/,
                            const size_type /*n_l*/,
                            graph_ptrtype const& graph )
{
    VLOG(1) << "MatrixPetsc init with graph";
    this->setGraph( graph );

    // Clear initialized matrices
    if ( this->isInitialized() )
        this->clear();

    this->setInitialized(  true );
    this->setIsClosed( true );

    int proc_id = 0;
    MPI_Comm_rank ( this->comm(), &proc_id );

    int m_global = static_cast<int>( this->graph()->mapRow().nDof()/*m*/ );
    int n_global = static_cast<int>( this->graph()->mapCol().nDof()/*n*/ );
    int m_local  = static_cast<int>( this->graph()->mapRow().nLocalDofWithoutGhost()/*m_l*/ );
    int n_local  = static_cast<int>( this->graph()->mapCol().nLocalDofWithoutGhost()/*n_l*/ );
    if ( m_global==0 )
        return;
    DVLOG(1) << "[MatrixPETSc::init()] proc_id   = " << proc_id << "\n";
    DVLOG(1) << "[MatrixPETSc::init()] m   = " << m_global << "\n";
    DVLOG(1) << "[MatrixPETSc::init()] n   = " << n_global << "\n";
    DVLOG(1) << "[MatrixPETSc::init()] m_l = " << m_local << "\n";
    DVLOG(1) << "[MatrixPETSc::init()] n_l = " << n_local << "\n";
    // Make sure the sparsity pattern isn't empty
    //FEELPP_ASSERT ( this->graph()->size() == n_l )( this->graph()->size() )( n_l ).warn( "incompatible diagonal non zero pattern" );
    DVLOG(1) << "[MatrixPETSc::init()] graph size   = " << this->graph()->size() << "\n";
    DVLOG(1) << "[MatrixPETSc::init()] graph first row entry on proc   = " << this->graph()->firstRowEntryOnProc() << "\n";
    DVLOG(1) << "[MatrixPETSc::init()] graph last row entry on proc   = " << this->graph()->lastRowEntryOnProc() << "\n";

    int ierr     = 0;

    // create a sequential matrix on one processor
    if ( ( m_local == m_global ) && ( n_local == n_global ) )
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

        for ( ; it != en; ++it )
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
        for( int i = 0; i < this->graph()->nNzOnProc().size(); ++i )
        {
            DLOG_IF(ERROR, dnz[i] == n_global ) << "row " << i << " is full : number of non zero entries = " << dnz[i] << " == cols " << n_global;
            LOG_IF(FATAL, dnz[i] > n_global ) << "row " << i << " has invalid data : number of non zero entries = " << dnz[i] << " == cols " << n_global;
        }

        ierr = MatCreateSeqAIJ ( this->comm(), m_global, n_global,
                                 0,
                                 dnz,
                                 //(int*) this->graph()->nNzOnProc().data(),
                                 &M_mat );
        CHKERRABORT( this->comm(),ierr );
        delete[] dnz;
        //ierr = MatSeqAIJSetPreallocation( M_mat, 0, (int*)this->graph()->nNzOnProc().data() );
#if 0
        ierr = MatSeqAIJSetPreallocation( M_mat, 0, dnz );
#else
        this->graph()->close();
        //std::cout << "sizes:" << this->graph()->ia().size() << "," << this->graph()->ja().size()  << std::endl;
        M_ia.resize( this->graph()->ia().size() );
        M_ja.resize( this->graph()->ja().size() );
        std::copy( this->graph()->ia().begin(), this->graph()->ia().end(), M_ia.begin() );
        std::copy( this->graph()->ja().begin(), this->graph()->ja().end(), M_ja.begin() );
        //std::for_each( ia.begin(), ia.end(), [](const int& i ){ std::cout << i << std::endl; } );

        ierr = MatSeqAIJSetPreallocationCSR( M_mat,
                                             M_ia.data(), M_ja.data(),
                                             this->graph()->a().data() );
#endif
        CHKERRABORT( this->comm(),ierr );

#else // MatCreateSeqAIJWithArrays

        //std::cout << "matrix csr creating...\n" << std::endl;
        this->graph()->close();
        //std::cout << "sizes:" << this->graph()->ia().size() << "," << this->graph()->ja().size()  << std::endl;
        M_ia.resize( this->graph()->ia().size() );
        M_ja.resize( this->graph()->ja().size() );
        std::copy( this->graph()->ia().begin(), this->graph()->ia().end(), M_ia.begin() );
        std::copy( this->graph()->ja().begin(), this->graph()->ja().end(), M_ja.begin() );
        //std::for_each( ia.begin(), ia.end(), [](const int& i ){ std::cout << i << std::endl; } );

        ierr = MatCreateSeqAIJWithArrays( this->comm(),
                                          m_global, n_global,
                                          M_ia.data(), M_ja.data(), this->graph()->a().data(), &M_mat );
        CHKERRABORT( this->comm(),ierr );
        //std::cout << "matrix csr created\n" << std::endl;
#endif
    }

    else
    {
#if PETSC_VERSION_LESS_THAN(3,3,0)
        ierr = MatCreateMPIAIJ ( this->comm(),
                                 m_local, n_local,
                                 m_global, n_global,
                                 0, ( int* ) this->graph()->nNzOnProc().data(),
                                 0, ( int* ) this->graph()->nNzOffProc().data(), &M_mat );
#else
        ierr = MatCreateAIJ ( this->comm(),
                                 m_local, n_local,
                                 m_global, n_global,
                                 0, ( int* ) this->graph()->nNzOnProc().data(),
                                 0, ( int* ) this->graph()->nNzOffProc().data(), &M_mat );
#endif
        CHKERRABORT( this->comm(),ierr );


    }

    ierr = MatSetFromOptions ( M_mat );
    CHKERRABORT( this->comm(),ierr );

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
    MatSetOption( M_mat,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS );
#endif


#if 1
    // additional insertions will not be allowed if they generate
    // a new nonzero
    //ierr = MatSetOption (M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    //CHKERRABORT(this->comm(),ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption ( M_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
    CHKERRABORT( this->comm(),ierr );
#endif

    if ( this->isSPD() )
    {
        ierr = MatSetOption ( M_mat, MAT_SPD, PETSC_TRUE );
        CHKERRABORT( this->comm(),ierr );
    }
    else if ( this->isSymmetric() )
    {
        ierr = MatSetOption ( M_mat, MAT_SYMMETRIC, PETSC_TRUE );
        CHKERRABORT( this->comm(),ierr );
    }
    else if ( this->isHermitian() )
    {
        ierr = MatSetOption ( M_mat, MAT_HERMITIAN, PETSC_TRUE );
        CHKERRABORT( this->comm(),ierr );
    }
    else if ( this->isStructurallySymmetric() )
    {
        ierr = MatSetOption ( M_mat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE );
        CHKERRABORT( this->comm(),ierr );
    }

    //MatShift( M_mat, 1 );
    //printMatlab( "shift.m" );
    //this->close();
    //this->zero();
    //this->close();

    //this->zeroEntriesDiagonal();
}


template <typename T>
void MatrixPetsc<T>::setIndexSplit( indexsplit_ptrtype const& indexSplit )
{
    //this->M_IndexSplit=indexSplit;
    super::setIndexSplit( indexSplit );

    PetscConvertIndexSplit( M_petscIS, *indexSplit, this->comm() );
    DVLOG(1) << "setIndexSplit done\n";
}

template <typename T>
void MatrixPetsc<T>::updatePCFieldSplit( PC & pc, indexsplit_ptrtype const& is )
{
    PetscConvertIndexSplit( M_petscIS, *is, this->comm() );
    updatePCFieldSplit( pc );
}

template <typename T>
void MatrixPetsc<T>::updatePCFieldSplit( PC & pc )
{
#if 1
    int ierr=0;

    if ( M_mapPC.find( &pc )==M_mapPC.end() )
    {
        M_mapPC[&pc]=false;
    }

    if ( M_mapPC[&pc]==false )
    {
#if PETSC_VERSION_LESS_THAN(3,4,0)
        const PCType pcName;
#else
        PCType pcName;
#endif
        ierr = PCGetType( pc,&pcName );
        CHKERRABORT( this->comm(),ierr );


        if ( std::string( PCFIELDSPLIT ) == std::string( pcName ) )
        {


            //std::cout << "\n updatePCFieldSplit " << M_petscIS.size() << "\n";
            M_mapPC[&pc]=true;

            for ( uint i = 0 ; i < M_petscIS.size(); ++i )
            {
                //std::cout << "\n split " << i << "\n";
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
                //std::cout << "\n version >= 3.2 \n";
                std::ostringstream os;
                os << i;
                ierr=PCFieldSplitSetIS( pc,os.str().c_str(),M_petscIS[i] );
#else
                //std::cout << "\n version < 3.2 \n";
                ierr=PCFieldSplitSetIS( pc,M_petscIS[i] );
#endif
                //std::cout << "\n split " << i << "done\n" << std::endl;

                CHKERRABORT( this->comm(),ierr );
            }
#if 0
            std::cout << " matrix petsc\n";
            ierr = PCFieldSplitSetBlockSize(pc,3);
            //CHKERRABORT( this->comm(),ierr );
            const PetscInt ufields[] = {0,2},pfields[] = {1};
            ierr = PCFieldSplitSetFields( pc ,"u", 2, ufields,ufields);
            CHKERRABORT( this->comm(),ierr );
            ierr = PCFieldSplitSetFields( pc , "p", 1, pfields,pfields);
            CHKERRABORT( this->comm(),ierr );
#endif

        }

    }

#else
    int ierr=0;

    const PCType pcName;
    ierr = PCGetType( pc,&pcName );
    CHKERRABORT( this->comm(),ierr );

    if ( std::string( PCFIELDSPLIT ) == std::string( pcName ) )
    {
        //std::cout << "\n updatePCFieldSplit \n";
        for ( uint i = 0 ; i < M_petscIS.size(); ++i )
        {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            ierr=PCFieldSplitSetIS( pc,PETSC_IGNORE,M_petscIS[i] );
#else
            ierr=PCFieldSplitSetIS( pc,M_petscIS[i] );
#endif
            CHKERRABORT( this->comm(),ierr );
        }



    }

#endif
}


template <typename T>
void MatrixPetsc<T>::zero ()
{
    CHECK( this->isInitialized() ) << "petsc matrix not properly initialized";

    if ( !this->closed() )
        this->close();

    int ierr=0;
    ierr = MatZeroEntries( M_mat );
    CHKERRABORT( this->comm(),ierr );

#if 0
    PetscBool is_assembled;
    MatAssembled( M_mat, &is_assembled );
    VLOG(2) << "Matrix is assembled : " << (is_assembled?"true":"false");
    if ( is_assembled )
    {
        ierr = MatZeroEntries( M_mat );
        CHKERRABORT( this->comm(),ierr );

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
                PetscInt row = it->second.template get<1>();
                std::copy( it->second.template get<2>().begin(), it->second.template get<2>().end(), cols.begin() );
                MatSetValues( M_mat, 1, &row, it->second.template get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
            }
        }
    }
#endif

}

template <typename T>
void MatrixPetsc<T>::zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{
    this->zero();
#if 0
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not properly initialized" ) ;

    int ierr=0;

    PetscBool is_assembled;
    MatAssembled( M_mat, &is_assembled );

    if ( is_assembled )
    {
        ierr = MatZeroEntries( M_mat );
        CHKERRABORT( this->comm(),ierr );

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
                PetscInt row = it->second.template get<1>();
                std::copy( it->second.template get<2>().begin(), it->second.template get<2>().end(), cols.begin() );
                MatSetValues( M_mat, 1, &row, it->second.template get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
            }
        }
    }
#endif
}


template <typename T>
void MatrixPetsc<T>::clear ()
{
    int ierr=0;
    PetscBool pinit;
    PetscInitialized( &pinit );
    if ( pinit && ( this->isInitialized() ) && ( this->M_destroy_mat_on_exit ) )
    {
        ierr = PETSc::MatDestroy ( M_mat );
        CHKERRABORT( this->comm(),ierr );

        for ( int i=0; i< ( int ) M_petscIS.size(); ++i )
        {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            ierr = ISDestroy( &M_petscIS[i] );
            CHKERRABORT( this->comm(),ierr );
#else
            ierr = ISDestroy( M_petscIS[i] );
            CHKERRABORT( this->comm(),ierr );
#endif
        }


        this->setInitialized( false );
        this->setIsClosed( false );
    }
}

template <typename T>
inline
void MatrixPetsc<T>::close ()  const
{
    if ( !this->isInitialized() )
        return;

    int ierr=0;
    PetscBool assembled = PETSC_FALSE;
    ierr = MatAssembled(M_mat,&assembled);
    VLOG(1) << "Matrix assembled ? " << assembled;
    //if ( !assembled )
    if ( 1 )
    {
        tic();
        // BSK - 1/19/2004
        // strictly this check should be OK, but it seems to
        // fail on matrix-free matrices.  Do they falsely
        // state they are assembled?  Check with the developers...
        //   if (this->closed())
        //     return;


        CHECK( M_mat ) << "invalid matrix";
        ierr = MatAssemblyBegin ( M_mat, MAT_FINAL_ASSEMBLY );
        CHKERRABORT( this->comm(),ierr );
        ierr = MatAssemblyEnd   ( M_mat, MAT_FINAL_ASSEMBLY );
        CHKERRABORT( this->comm(),ierr );
        toc("MatrixPETSc::close",FLAGS_v>0);
    }
    this->setIsClosed( true );
    //const_cast<MatrixPetsc<T>*>( this )->setIsClosed( true );
}



template <typename T>
inline
typename MatrixPetsc<T>::size_type
MatrixPetsc<T>::size1 () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize ( M_mat, &petsc_m, &petsc_n );
    CHKERRABORT( this->comm(),ierr );
    return static_cast<size_type>( petsc_m );
}



template <typename T>
inline
typename MatrixPetsc<T>::size_type
MatrixPetsc<T>::size2 () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    int petsc_m=0, petsc_n=0, ierr=0;

    ierr = MatGetSize ( M_mat, &petsc_m, &petsc_n );
    CHKERRABORT( this->comm(),ierr );
    return static_cast<size_type>( petsc_n );
}



template <typename T>
inline
typename MatrixPetsc<T>::size_type
MatrixPetsc<T>::rowStart () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange( M_mat, &start, &stop );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<size_type>( start );
}



template <typename T>
inline
typename MatrixPetsc<T>::size_type
MatrixPetsc<T>::rowStop () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange( M_mat, &start, &stop );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<size_type>( stop );
}

template <typename T>
inline
typename MatrixPetsc<T>::size_type
MatrixPetsc<T>::nnz () const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    MatInfo info;
    int ierr = MatGetInfo( M_mat, MAT_GLOBAL_SUM, &info);

    CHKERRABORT( this->comm(),ierr );
    return info.nz_allocated;
}



template <typename T>
inline
void MatrixPetsc<T>::set ( const size_type i,
                           const size_type j,
                           const value_type& value )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>( value );
    ierr = MatSetValues( M_mat, 1, &i_val, 1, &j_val,
                         &petsc_value, INSERT_VALUES );
    CHKERRABORT( this->comm(),ierr );
}



template <typename T>
inline
void MatrixPetsc<T>::add ( const size_type i,
                           const size_type j,
                           const value_type& value )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );

    int ierr=0, i_val=i, j_val=j;
    //DVLOG(2) << "[MatrixPetsc<>::add] adding value " << value << " at (" << i << "," << j << ")\n";
    PetscScalar petsc_value = static_cast<PetscScalar>( value );


    ierr = MatSetValues( M_mat, 1, &i_val, 1, &j_val,
                         &petsc_value, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}
/*                                   */

template <typename T>
inline
bool MatrixPetsc<T>::closed() const
{
    return super::closed();
#if 0
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );
    int ierr=0;
    PetscTruth assembled;
    ierr = MatAssembled( M_mat, &assembled );
    CHKERRABORT( this->comm(),ierr );
    return ( assembled == PETSC_TRUE );
#endif
}
template <typename T>
void
MatrixPetsc<T>::addMatrix( const ublas::matrix<value_type>& dm,
                           const std::vector<size_type>& rows,
                           const std::vector<size_type>& cols )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    const size_type m = dm.size1();
    const size_type n = dm.size2();

    FEELPP_ASSERT ( rows.size() == this->size1() ).error( "invalid row size" );
    FEELPP_ASSERT ( cols.size() == this->size2() ).error( "invalid column size" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues( M_mat,
                         m, ( int* ) boost::addressof( rows[0] ),
                         n, ( int* ) boost::addressof( cols[0] ),
                         ( PetscScalar* ) dm.data().begin(),
                         ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}
template <typename T>
void
MatrixPetsc<T>::addMatrix ( int* rows, int nrows,
                            int* cols, int ncols,
                            value_type* data,
                            size_type K,
                            size_type K2 )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValues( M_mat,
                         nrows, ( int* ) rows,
                         ncols, ( int* ) cols,
                         ( PetscScalar* ) data,
                         ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
MatrixPetsc<T>::setDiagonal( const Vector<T>& vecDiag )
{
    if ( !vecDiag.closed() )
        const_cast<Vector<T>*>( &vecDiag )->close();

    this->setIsClosed( false );

    // all data are in petsc type
    VectorPetsc<T> const* vecDiag_petsc = dynamic_cast<VectorPetsc<T> const*>( &vecDiag );
    if ( vecDiag_petsc )
    {
        int ierr = 0;
        ierr = MatDiagonalSet(this->mat(),vecDiag_petsc->vec(),INSERT_VALUES);
        CHKERRABORT( this->comm(),ierr );
        return;
    }
    // others vector type
    auto vecPetscCastPair = Feel::detail::toPETScPairPtr( vecDiag, true );
    VectorPetsc<T> const* vecDiag_petscCast = vecPetscCastPair.first;
    if ( vecDiag_petscCast )
    {
        this->setDiagonal( *vecDiag_petscCast );
        return;
    }
    CHECK( false ) << "invalid vector type";
}

template <typename T>
void
MatrixPetsc<T>::addDiagonal( const Vector<T>& vecDiag )
{
    if ( !vecDiag.closed() )
        const_cast<Vector<T>*>( &vecDiag )->close();

    this->setIsClosed( false );

    // all data are in petsc type
    VectorPetsc<T> const* vecDiag_petsc = dynamic_cast<VectorPetsc<T> const*>( &vecDiag );
    if ( vecDiag_petsc )
    {
        int ierr = 0;
        ierr = MatDiagonalSet(this->mat(),vecDiag_petsc->vec(),ADD_VALUES);
        CHKERRABORT( this->comm(),ierr );
        return;
    }
    // others vector type
    auto vecPetscCastPair = Feel::detail::toPETScPairPtr( vecDiag, true );
    VectorPetsc<T> const* vecDiag_petscCast = vecPetscCastPair.first;
    if ( vecDiag_petscCast )
    {
        this->addDiagonal( *vecDiag_petscCast );
        return;
    }
    CHECK( false ) << "invalid vector type";

}

template <typename T>
void
MatrixPetsc<T>::multVector( const Vector<T>& arg, Vector<T>& dest, bool transpose ) const
{
    CHECK( this->isInitialized() ) << "is not initialized";
    CHECK( arg.isInitialized() ) << "is not initialized";

    if ( !this->closed() )
        this->close();
    if ( !arg.closed() )
        const_cast<Vector<T>*>( &arg )->close();
    if ( !dest.isInitialized() )
        dest.init( (transpose)? this->mapRowPtr() :this->mapColPtr() );
    else if ( !dest.closed() )
        dest.close();

    // all data are in petsc type
    VectorPetsc<T> const* vec_petsc_arg = dynamic_cast<VectorPetsc<T> const*>( &arg );
    VectorPetsc<T> * vec_petsc_dest = dynamic_cast<VectorPetsc<T> *>( &dest );
    if ( vec_petsc_arg && vec_petsc_dest )
    {
        int ierr = 0;
        if ( !transpose )
        {
            if ( this->mapCol().worldComm().globalSize() == arg.map().worldComm().globalSize() )
            {
                ierr = MatMult( this->mat(), vec_petsc_arg->vec(), vec_petsc_dest->vec() );
                CHKERRABORT( this->comm(),ierr );
            }
            else
            {
                auto x_convert = VectorPetscMPI<T>(this->mapColPtr());
                x_convert.duplicateFromOtherPartition(arg);
                x_convert.close();
                ierr = MatMult( this->mat(), x_convert.vec(), vec_petsc_dest->vec() );
                CHKERRABORT( this->comm(),ierr );
            }
        }
        else
        {
            if ( this->mapRow().worldComm().globalSize() == arg.map().worldComm().globalSize() )
            {
                ierr = MatMultTranspose( this->mat(), vec_petsc_arg->vec(), vec_petsc_dest->vec() );
                CHKERRABORT( this->comm(),ierr );
            }
            else
            {
                auto x_convert = VectorPetscMPI<T>(this->mapRowPtr());
                x_convert.duplicateFromOtherPartition(arg);
                x_convert.close();
                ierr = MatMultTranspose( this->mat(), x_convert.vec(), vec_petsc_dest->vec() );
                CHKERRABORT( this->comm(),ierr );
            }
        }
        // update ghost if parallel vector
        vec_petsc_dest->localize();
        return;
    }

    // others vector type
    auto vecPetscCastPair_arg = Feel::detail::toPETScPairPtr( arg, true );
    VectorPetsc<T> const* vecPetscCast_arg = vecPetscCastPair_arg.first;
    auto vecPetscCastPair_dest = Feel::detail::toPETScPairPtr( dest );
    VectorPetsc<T> * vecPetscCast_dest = vecPetscCastPair_dest.first;
    CHECK( vecPetscCast_dest ) << "vector dest type can not be convert into a petsc vector";
    this->multVector( *vecPetscCast_arg, *vecPetscCast_dest, transpose );
    return;
}

template <typename T>
void
MatrixPetsc<T>::matMatMult ( MatrixSparse<T> const& matIn, MatrixSparse<T> &matRes ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );
    FEELPP_ASSERT( this->size2() == matIn.size1() )( this->size2() )( matIn.size1() ).error( "incompatible dimension" );

    if ( !this->closed() )
        this->close();
    if ( matIn.closed() )
        matIn.close();
    if ( !matRes.closed() )
        matRes.close();

    MatrixPetsc<T> const* matInPetsc = dynamic_cast<MatrixPetsc<T> const*> ( &matIn );
    MatrixPetsc<T>* matResPetsc = dynamic_cast<MatrixPetsc<T>*> ( &matRes );
    if ( matInPetsc && matResPetsc )
    {
        int ierr=0;
        if ( matResPetsc->isInitialized() )
        {
            ierr = MatMatMult(this->M_mat, matInPetsc->mat(), MAT_REUSE_MATRIX, PETSC_DEFAULT, &matResPetsc->mat());
            CHKERRABORT( this->comm(),ierr );
        }
        else
        {
            ierr = MatMatMult(this->M_mat, matInPetsc->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &matResPetsc->mat());
            CHKERRABORT( this->comm(),ierr );

            matRes.setMapRow( this->mapRowPtr() );
            matRes.setMapCol( matIn.mapColPtr() );
            MatrixPetscMPI<T>* matResPetscMPI = dynamic_cast<MatrixPetscMPI<T>*> ( &matRes );
            if ( matResPetscMPI )
                matResPetscMPI->initLocalToGlobalMapping();
            matRes.setInitialized( true );
            matRes.setIsClosed( true );
            // TODO graph
        }
    }
    else
    {
        CHECK( false ) << "TODO other kind of matrix";
    }

}


template <typename T>
void
MatrixPetsc<T>::PtAP( MatrixSparse<value_type> const& matP, MatrixSparse<value_type> & matC ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );
    FEELPP_ASSERT( this->size1() == matP.size1() )( this->size1() )( matP.size1() ).error( "incompatible dimension" );
    FEELPP_ASSERT( this->size2() == matP.size1() )( this->size2() )( matP.size1() ).error( "incompatible dimension" );

    if ( !this->closed() )
        this->close();
    if ( matP.closed() )
        matP.close();
    if ( !matC.closed() )
        matC.close();

    MatrixPetsc<T> const* matP_petsc = dynamic_cast<MatrixPetsc<T> const*> ( &matP );
    MatrixPetsc<T>* matC_petsc = dynamic_cast<MatrixPetsc<T>*> ( &matC );
    if ( matP_petsc && matC_petsc )
    {
        int ierr=0;
        if ( matC.isInitialized() )
        {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,3,0)
            ierr = MatPtAP( this->mat(), matP_petsc->mat(), MAT_REUSE_MATRIX, 1.0, &matC_petsc->mat() );
            CHKERRABORT( this->comm(),ierr );
#else
            CHECK( false ) << "PtAP requiert a PETSc version  >= 3.3";
#endif
        }
        else
        {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,3,0)
            ierr = MatPtAP( this->mat(), matP_petsc->mat(), MAT_INITIAL_MATRIX, 1.0, &matC_petsc->mat() );
            CHKERRABORT( this->comm(),ierr );
#else
            CHECK( false ) << "PtAP requiert a PETSc version  >= 3.3";
#endif
            matC.setMapRow( matP.mapColPtr() );
            matC.setMapCol( matP.mapColPtr() );
            MatrixPetscMPI<T>* matC_petscMPI = dynamic_cast<MatrixPetscMPI<T>*> ( &matC );
            if ( matC_petscMPI )
                matC_petscMPI->initLocalToGlobalMapping();
            matC.setInitialized( true );
            matC.setIsClosed( true );
            // TODO graph
        }
        return;
    }

    CHECK( false ) << "TODO other kind of matrix";
}

template <typename T>
void
MatrixPetsc<T>::PAPt( MatrixSparse<value_type> const& matP, MatrixSparse<value_type> & matC ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );
    FEELPP_ASSERT( this->size1() == matP.size2() )( this->size1() )( matP.size2() ).error( "incompatible dimension" );
    FEELPP_ASSERT( this->size2() == matP.size2() )( this->size2() )( matP.size2() ).error( "incompatible dimension" );

    if ( !this->closed() )
        this->close();
    if ( matP.closed() )
        matP.close();
    if ( !matC.closed() )
        matC.close();

    MatrixPetsc<T> const* matP_petsc = dynamic_cast<MatrixPetsc<T> const*> ( &matP );
    MatrixPetsc<T>* matC_petsc = dynamic_cast<MatrixPetsc<T>*> ( &matC );
    if ( matP_petsc && matC_petsc )
    {
        int ierr=0;
        if ( matC.isInitialized() )
        {
            if ( this->comm().size() == 1 )
            {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 3)
                ierr = MatRARt( this->mat(), matP_petsc->mat(), MAT_REUSE_MATRIX, 1.0, &matC_petsc->mat() );
                CHKERRABORT( this->comm(),ierr );
#else
                CHECK( false ) << "PtAP requiert a PETSc version  >= 3.3";
#endif
            }
            else
            {
                MatrixPetscMPI<T> matTemp_petscMPI( this->mapRowPtr(),matP.mapRowPtr() );
                matP.matMatMult( *this, matTemp_petscMPI );
                matTemp_petscMPI.matMatMult( *matP.transpose(), matC );
            }
        }
        else
        {
            if ( this->comm().size() == 1 )
            {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 3)
                ierr = MatRARt( this->mat(), matP_petsc->mat(), MAT_INITIAL_MATRIX, 1.0, &matC_petsc->mat() );
                CHKERRABORT( this->comm(),ierr );
#else
                CHECK( false ) << "PtAP requiert a PETSc version  >= 3.3";
#endif
            }
            else
            {
                MatrixPetscMPI<T> matTemp_petscMPI( this->mapRowPtr(),matP.mapRowPtr() );
                matP.matMatMult( *this, matTemp_petscMPI );
                matTemp_petscMPI.matMatMult( *matP.transpose(), matC );
            }

            matC.setMapRow( matP.mapRowPtr() );
            matC.setMapCol( matP.mapRowPtr() );
            MatrixPetscMPI<T>* matC_petscMPI = dynamic_cast<MatrixPetscMPI<T>*> ( &matC );
            if ( matC_petscMPI )
                matC_petscMPI->initLocalToGlobalMapping();
            matC.setInitialized( true );
            matC.setIsClosed( true );
            // TODO graph
        }
        return;
    }

    CHECK( false ) << "TODO other kind of matrix";
}




#if 0
template <typename T>
void
MatrixPetsc<T>::printPython( const std::string name ) const
{
    this->close();
    auto nbRow = this->size1();

    int ierr = 0;

    std::ofstream graphFile( name, std::ios::out /*| std::ios::app*/ );


    for ( size_type i = 0 ; i < nbRow ; ++i )
    {
        PetscInt row = i;
        PetscInt ncols;
        const PetscInt *cols;
        const PetscScalar *vals;

        PetscInt nrow2=1;
        PetscInt *row2 = new PetscInt[1];
        row2[0] = i;

        ierr = MatGetRow( this->mat(), row, &ncols, &cols, &vals );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        for ( size_type jj = 0 ; jj < ncols ; ++jj )
        {
            /*PetscInt *cols2 = new PetscInt[1];cols2[0]=cols[jj];
            const PetscScalar *vals3;
            ierr = MatGetValues(this->mat(), nrow2, row2, 1, cols2, vals3);
            CHKERRABORT(PETSC_COMM_WORLD,ierr);
            */

            if ( std::abs( vals[jj] )>1.e-8 )graphFile << "[" << row <<","<< cols[jj] << "," << vals[jj] << "],";

            //std::cout << "\n [updateBlockMat] i "<< start_i+row << " j " << start_j+cols[jj] << " val " << vals[jj] << std::endl;
            //this->set(start_i+row,
            //      start_j+cols[jj],
            //      vals[jj]);
        }

        ierr = MatRestoreRow( this->mat(),row,&ncols,&cols,&vals );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }

    graphFile.close();
}
#endif

// print
template <typename T>
void
MatrixPetsc<T>::printMatlab ( const std::string name ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not properly initialized" );

    // assert (this->closed());
    if ( !this->closed() )
        this->close();
    PetscObjectSetName((PetscObject)M_mat,fs::path("var_"+name).stem().string().c_str());
    int ierr=0;
    PetscViewer petsc_viewer;


    ierr = PetscViewerCreate ( this->comm(),
                               &petsc_viewer );
    CHKERRABORT( this->comm(),ierr );

    /**
     * Create an ASCII file containing the matrix
     * if a filename was provided.
     */
    if ( name != "NULL" )
    {
#if 0
        ierr = PetscViewerBinaryOpen( this->comm(),
                                     name.c_str(),
                                     FILE_MODE_WRITE,
                                     &petsc_viewer );
#else
        ierr = PetscViewerASCIIOpen( this->comm(),
                                     name.c_str(),
                                     &petsc_viewer );
#endif
        CHKERRABORT( this->comm(),ierr );
#if 0
#if PETSC_VERSION_LESS_THAN(3,7,0)
        ierr = PetscViewerSetFormat ( petsc_viewer,
                                      PETSC_VIEWER_BINARY_MATLAB );
#else
        ierr = PetscViewerPushFormat ( petsc_viewer,
                                      PETSC_VIEWER_BINARY_MATLAB );
#endif
#else
        ierr = PetscViewerPushFormat ( petsc_viewer,
                                       PETSC_VIEWER_ASCII_MATLAB );
#endif
        //PETSC_VIEWER_ASCII_PYTHON );

        CHKERRABORT( this->comm(),ierr );

        ierr = MatView ( M_mat, petsc_viewer );
        CHKERRABORT( this->comm(),ierr );
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
        ierr = PetscViewerPopFormat ( petsc_viewer );
        CHKERRABORT( this->comm(),ierr );
#endif
    }

    /**
     * Otherwise the matrix will be dumped to the screen.
     */
    else
    {
        //ierr = PetscViewerSetFormat (petsc_viewer,//PETSC_VIEWER_STDOUT_SELF,//PETSC_VIEWER_STDOUT_WORLD,
        //                             PETSC_VIEWER_ASCII_INFO_DETAIL);//PETSC_VIEWER_ASCII_MATLAB);
        //CHKERRABORT(this->comm(),ierr);

        //ierr = MatView (M_mat,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_SELF);//
        //CHKERRABORT(this->comm(),ierr);
    }


    /**
     * Destroy the viewer.
     */
    ierr = PETSc::PetscViewerDestroy ( petsc_viewer );
    CHKERRABORT( this->comm(),ierr );
}


template <typename T>
std::shared_ptr<MatrixSparse<T> >
MatrixPetsc<T>::createSubMatrix( std::vector<size_type> const& _rows,
                                 std::vector<size_type> const& _cols,
                                 bool useSameDataMap, bool checkAndFixRange ) const
{
    if ( useSameDataMap )
        CHECK( _rows.size() == _cols.size() ) << "if useSameDataMap we must have equality of index set ( _rows == _cols )";

    // update maybe input index set
    std::vector<size_type> rows = ( checkAndFixRange )?
        this->mapRowPtr()->buildIndexSetWithParallelMissingDof( _rows ) : _rows;
    std::vector<size_type> cols = ( checkAndFixRange && !useSameDataMap )?
        this->mapColPtr()->buildIndexSetWithParallelMissingDof( _cols ) : ( useSameDataMap )? rows : _cols;

    // build subdatamap row and col
    datamap_ptrtype subMapRow = this->mapRowPtr()->createSubDataMap( rows, false );
    datamap_ptrtype subMapCol = ( useSameDataMap )? subMapRow : this->mapColPtr()->createSubDataMap( cols, false );

    // build submatrix petsc
    Mat subMatPetsc = NULL;
    this->getSubMatrixPetsc( rows,cols,subMatPetsc );

    // build matrixsparse object
    std::shared_ptr<MatrixSparse<T> > subMat;
    if ( this->comm().size()>1 )
        subMat.reset( new MatrixPetscMPI<T>( subMatPetsc,subMapRow,subMapCol,true,true ) );
    else
        subMat.reset( new MatrixPetsc<T>( subMatPetsc,subMapRow,subMapCol,true ) );

    bool computeSubGraph = true;
    if ( computeSubGraph && this->hasGraph() )
    {
        auto subgraph = this->graph()->createSubGraph( rows, cols, subMapRow, subMapCol, useSameDataMap, false );
        subMat->setGraph( subgraph );
    }

    return subMat;
}

template <typename T>
void
MatrixPetsc<T>::updateSubMatrix( std::shared_ptr<MatrixSparse<T> > & submatrix,
                                 std::vector<size_type> const& rows,
                                 std::vector<size_type> const& cols, bool doClose )
{
    CHECK( submatrix ) << "submatrix is not init";
    std::shared_ptr<MatrixPetsc<T> > submatrixPetsc = std::dynamic_pointer_cast<MatrixPetsc<T> >( submatrix );
    this->getSubMatrixPetsc( rows,cols, submatrixPetsc->mat() , doClose );
}


template <typename T>
void
MatrixPetsc<T>::getSubMatrixPetsc( std::vector<size_type> const& rows,
                                   std::vector<size_type> const& cols,
                                   Mat &submat, bool doClose ) const
{
    if ( !this->closed() )
        this->close();
    //if(doClose) // with close(), petsc believe the matrix has changed. That can be confusing
          //    const_cast<MatrixPetsc<T>*>( this )->close();
    int ierr=0;
    IS isrow;
    IS iscol;
    PetscInt *rowMap;
    PetscInt *colMap;

    std::set<size_type> rowMapOrdering, colMapOrdering;
    if ( this->comm().size()>1 )
    {
        // convert global process ids into global cluster ids, remove ghost dofs
        // and build ordering row and col map
        for (int i=0; i<rows.size(); i++)
        {
            if ( !this->mapRow().dofGlobalProcessIsGhost( rows[i] ) )
                rowMapOrdering.insert( this->mapRow().mapGlobalProcessToGlobalCluster( rows[i] ) );
        }
        for (int i=0; i<cols.size(); i++)
        {
            if ( !this->mapCol().dofGlobalProcessIsGhost( cols[i] ) )
                colMapOrdering.insert( this->mapCol().mapGlobalProcessToGlobalCluster( cols[i] ) );
        }
    }
    else
    {
        // build ordering row and col map
        rowMapOrdering.insert( rows.begin(), rows.end() );
        colMapOrdering.insert( cols.begin(), cols.end() );
    }

    // copying into PetscInt vector
    int nrow = rowMapOrdering.size();
    int ncol = colMapOrdering.size();
    rowMap = new PetscInt[nrow];
    colMap = new PetscInt[ncol];
    size_type curId=0;
    for ( auto& rowId : rowMapOrdering )
    {
        rowMap[curId] = rowId;
        ++curId;
    }
    curId=0;
    for ( auto& colId : colMapOrdering )
    {
        colMap[curId] = colId;
        ++curId;
    }


#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(),nrow,rowMap,PETSC_COPY_VALUES,&isrow);
    CHKERRABORT( this->comm(),ierr );
    ierr = ISCreateGeneral(this->comm(),ncol,colMap,PETSC_COPY_VALUES,&iscol);
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISCreateGeneral(this->comm(),nrow,rowMap,&isrow);
    CHKERRABORT( this->comm(),ierr );
    ierr = ISCreateGeneral(this->comm(),ncol,colMap,&iscol);
    CHKERRABORT( this->comm(),ierr );
#endif
#if PETSC_VERSION_LESS_THAN( 3,8,0 )
    if ( submat == NULL )
        ierr = MatGetSubMatrix(this->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &submat);
    else
        ierr = MatGetSubMatrix(this->mat(), isrow, iscol, MAT_REUSE_MATRIX, &submat);
#else
    if ( submat == NULL )
        ierr = MatCreateSubMatrix(this->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &submat);
    else
        ierr = MatCreateSubMatrix(this->mat(), isrow, iscol, MAT_REUSE_MATRIX, &submat);
#endif
    CHKERRABORT( this->comm(),ierr );

    ierr = PETSc::ISDestroy( isrow );
    CHKERRABORT( this->comm(),ierr );
    ierr = PETSc::ISDestroy( iscol );
    CHKERRABORT( this->comm(),ierr );

    delete[] rowMap;
    delete[] colMap;
}

// previous createSubmatrix
template <typename T>
void
MatrixPetsc<T>::createSubmatrix( MatrixSparse<T>& submatrix,
                                 const std::vector<size_type>& rows,
                                 const std::vector<size_type>& cols ) const
{
    //auto sub_graph = GraphCSR(rows.size()*cols.size(),0,rows.size()-1,0,cols.size()-1,this->comm());
    auto sub_graph = graph_ptrtype(new graph_type(0,0,rows.size(),0,cols.size(),this->worldCommPtr()));

    MatrixPetsc<T>* A = dynamic_cast<MatrixPetsc<T>*> ( &submatrix );
    A->setGraph( sub_graph );

    int ierr=0;
    IS isrow;
    IS iscol;
    PetscInt *rowMap;
    PetscInt *colMap;
    int nrow = rows.size();
    int ncol = cols.size();

    rowMap = new PetscInt[nrow];
    colMap = new PetscInt[ncol];

    for (int i=0; i<nrow; i++) rowMap[i] = rows[i];
    for (int i=0; i<ncol; i++) colMap[i] = cols[i];

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(),nrow,rowMap,PETSC_COPY_VALUES,&isrow);
    CHKERRABORT( this->comm(),ierr );
    ierr = ISCreateGeneral(this->comm(),ncol,colMap,PETSC_COPY_VALUES,&iscol);
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISCreateGeneral(this->comm(),nrow,rowMap,&isrow);
    CHKERRABORT( this->comm(),ierr );
    ierr = ISCreateGeneral(this->comm(),ncol,colMap,&iscol);
    CHKERRABORT( this->comm(),ierr );
#endif
#if PETSC_VERSION_LESS_THAN( 3,8,0 )
    ierr = MatGetSubMatrix(this->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &A->mat());
#else
    ierr = MatCreateSubMatrix(this->mat(), isrow, iscol, MAT_INITIAL_MATRIX, &A->mat());
#endif
    CHKERRABORT( this->comm(),ierr );

    PETSc::ISDestroy( isrow );
    PETSc::ISDestroy( iscol );
    delete[] rowMap;
    delete[] colMap;
}

template <typename T>
inline
void
MatrixPetsc<T>::addMatrix ( const T a_in, MatrixSparse<T> const&X_in, Feel::MatrixStructure matStruc )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    // sanity check. but this cannot avoid
    // crash due to incompatible sparsity structure...
    FEELPP_ASSERT( this->size1() == X_in.size1() )( this->size1() )( X_in.size1() ).error( "incompatible dimension" );
    FEELPP_ASSERT( this->size2() == X_in.size2() )( this->size2() )( X_in.size2() ).error( "incompatible dimension" );

    // the matrix have to be assembled/closed
    if ( !this->closed() )
        this->close();
    if ( !X_in.closed() )
        X_in.close();


    PetscScalar     a = static_cast<PetscScalar>      ( a_in );
    MatrixPetsc<T> const* matPetscIn = dynamic_cast<MatrixPetsc<T> const*> ( &X_in );

    if ( matPetscIn )
    {
        int ierr=0;

        // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
        ierr = MatAXPY( &a,  matPetscIn->M_mat, M_mat, MatStructure::SAME_NONZERO_PATTERN );
        // 2.3.x & newer
#else
        ierr = MatAXPY(M_mat, a, matPetscIn->M_mat, PetscGetMatStructureEnum( matStruc ) /*MatStructure::SAME_NONZERO_PATTERN*/);
        //ierr = MatAXPY(M_mat, a, X->M_mat, MatStructure::SUBSET_NONZERO_PATTERN );
        //ierr = MatAXPY( M_mat, a, X->M_mat, MatStructure::DIFFERENT_NONZERO_PATTERN );
        //ierr = MatDuplicate(X->mat(),MAT_COPY_VALUES,&M_mat);
#endif
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        CHECK( false ) << "TODO other kind of matrix";
    }

}

template <typename T>
void
MatrixPetsc<T>::scale( T const a )
{
    if ( !this->closed() )
        this->close();
    int ierr = MatScale( M_mat, a );
    CHKERRABORT( this->comm(),ierr );
}
template <typename T>
inline
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::l1Norm() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );
    if ( !this->closed() )
        this->close();

    int ierr=0;
    double petsc_value;
    real_type value;

    ierr = MatNorm( M_mat, NORM_1, &petsc_value );
    CHKERRABORT( this->comm(),ierr );

    value = static_cast<real_type>( petsc_value );

    return value;
}

template <typename T>
inline
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::linftyNorm() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );
    if ( !this->closed() )
        this->close();

    int ierr=0;
    double petsc_value;
    real_type value;

    ierr = MatNorm( M_mat, NORM_INFINITY, &petsc_value );
    CHKERRABORT( this->comm(),ierr );

    value = static_cast<real_type>( petsc_value );

    return value;
}

template <typename T>
MatrixPetsc<T> &
MatrixPetsc<T>::operator = ( MatrixSparse<value_type> const& M )
{
    if ( !this->closed() )
        this->close();
    if ( !M.closed() )
        M.close();

    MatrixPetsc<T> const* X = dynamic_cast<MatrixPetsc<T> const*> ( &M );
    this->setGraph( M.graph() );

    //MatConvert(X->mat(), MATSAME, MAT_INITIAL_MATRIX, &M_mat);
    int ierr= MatDuplicate( X->mat(),MAT_COPY_VALUES,&M_mat );
    CHKERRABORT( this->comm(),ierr );

    return *this;
}
template <typename T>
inline
typename MatrixPetsc<T>::value_type
MatrixPetsc<T>::operator () ( const size_type i,
                              const size_type j ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

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
    i_val=static_cast<int>( i ),
    j_val=static_cast<int>( j );


    // the matrix needs to be closed for this to work
    if ( !this->closed() )
        this->close();

    if ( this->comm().size()>1 )
    {
        i_val = this->mapRow().mapGlobalProcessToGlobalCluster( i );
        j_val = this->mapCol().mapGlobalProcessToGlobalCluster( j );
        ierr = MatGetRow( M_mat, i_val, &ncols, &petsc_cols, &petsc_row );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = MatGetRow( M_mat, i_val, &ncols, &petsc_cols, &petsc_row );
        CHKERRABORT( this->comm(),ierr );
    }

    // Perform a binary search to find the contiguous index in
    // petsc_cols (resp. petsc_row) corresponding to global index j_val
    std::pair<const int*, const int*> p =
        std::equal_range ( &petsc_cols[0], &petsc_cols[0] + ncols, j_val );

    // Found an entry for j_val
    if ( p.first != p.second )
    {
        // The entry in the contiguous row corresponding
        // to the j_val column of interest
        const int j = std::distance ( const_cast<int*>( &petsc_cols[0] ),
                                      const_cast<int*>( p.first ) );

        assert ( j < ncols );
        assert ( petsc_cols[j] == j_val );

        value = static_cast<T> ( petsc_row[j] );

        ierr  = MatRestoreRow( M_mat, i_val,
                               &ncols, &petsc_cols, &petsc_row );
        CHKERRABORT( this->comm(),ierr );

        return value;
    }

    ierr  = MatRestoreRow( M_mat, i_val,
                           &ncols, &petsc_cols, &petsc_row );
    CHKERRABORT( this->comm(),ierr );

    // Otherwise the entry is not in the sparse matrix,
    // i.e. it is 0.
    return 0.;
}


template<typename T>
void
MatrixPetsc<T>::zeroRows( std::vector<int> const& rows,
                          Vector<value_type> const& values,
                          Vector<value_type>& rhs,
                          Context const& on_context,
                          value_type value_on_diagonal )
{
    // the matrix needs to be closed for this to work
    if ( !this->closed() )
        this->close();
    if ( !rhs.closed() )
        rhs.close();

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    MatSetOption( M_mat,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( M_mat,MAT_KEEP_ZEROED_ROWS );
#endif

    if ( on_context.test( ContextOn::ELIMINATION) )
    {
        LOG(INFO) << "MatrixPETSc:: zeroRows seq elimination";
        rhs.setIsClosed( false );
        VectorPetsc<T>* prhs = dynamic_cast<VectorPetsc<T>*> ( &rhs );
        const VectorPetsc<T>* pvalues = dynamic_cast<const VectorPetsc<T>*> ( &values );

        int start, stop;
        int ierr = MatGetOwnershipRange( M_mat, &start, &stop );
        CHKERRABORT( this->comm(),ierr );

        VectorPetsc<value_type> diag( this->mapColPtr() );
        if ( on_context.test( ContextOn::KEEP_DIAGONAL ) )
        {
            LOG(INFO) << "MatrixPETSc:: zeroRows seq getdiag";
            MatGetDiagonal( M_mat, diag.vec() );
        }

        if ( on_context.test( ContextOn::SYMMETRIC ) )
        {
            LOG(INFO) << "MatrixPETSc:: zeroRows seq symmetric";
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
            MatZeroRowsColumns(M_mat, rows.size(), rows.data(), value_on_diagonal, pvalues->vec(), prhs->vec() );
#else
            MatZeroRows( M_mat, rows.size(), rows.data(), value_on_diagonal );
#endif
            if ( on_context.test( ContextOn::KEEP_DIAGONAL ) )
            {
                LOG(INFO) << "MatrixPETSc:: zeroRows seq setdiag";
                MatDiagonalSet( M_mat, diag.vec(), INSERT_VALUES );

                for ( size_type i = 0; i < rows.size(); ++i )
                {
                    // warning: a row index may belong to another
                    // processor, so make sure that we access only the
                    // rows that belong to this processor
                    if ( rows[i] >= start && rows[i] < stop )
                        rhs.set( rows[i], values(rows[i])*diag( rows[i] ) );
                }
            }
        }
        else // non symmetric case
        {
            LOG(INFO) << "MatrixPETSc:: zeroRows seq unsymmetric";
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
            MatZeroRows( M_mat, rows.size(), rows.data(), value_on_diagonal,PETSC_IGNORE,PETSC_IGNORE );
#else
            MatZeroRows( M_mat, rows.size(), rows.data(), value_on_diagonal );
#endif
            if ( on_context.test( ContextOn::KEEP_DIAGONAL ) )
            {
                LOG(INFO) << "MatrixPETSc:: zeroRows seq set diag";
                MatDiagonalSet( M_mat, diag.vec(), INSERT_VALUES );
                for ( size_type i = 0; i < rows.size(); ++i )
                {
                    // warning: a row index may belong to another
                    // processor, so make sure that we access only the
                    // rows that belong to this processor
                    if ( rows[i] >= start && rows[i] < stop )
                        rhs.set( rows[i], values(rows[i])*diag( rows[i] ) );
                }
            }
            else
            {
                for ( size_type i = 0; i < rows.size(); ++i )
                {
                    // eliminate column

                    // warning: a row index may belong to another
                    // processor, so make sure that we access only the
                    // rows that belong to this processor
                    if ( rows[i] >= start && rows[i] < stop )
                        rhs.set( rows[i], values(rows[i]) );
                }
            }
        }
    }

    if ( !rhs.closed() )
        rhs.close();

}
template<typename T>
void
MatrixPetsc<T>::diagonal( Vector<value_type>& out ) const
{
    if ( !this->closed() )
        this->close();
    if ( !out.isInitialized() )
        out.init( this->mapRowPtr() );
    else if ( !out.closed() )
        const_cast<Vector<T>*>( &out )->close();

    VectorPetsc<T>* vecOutPetsc = dynamic_cast<VectorPetsc<T>*> ( &out );
    if ( vecOutPetsc )
    {
        int ierr = MatGetDiagonal( M_mat, vecOutPetsc->vec() );
        CHKERRABORT( this->comm(),ierr );
        // update ghost in //
        vecOutPetsc->localize();
        return;
    }

    auto vecPetscCastPair = Feel::detail::toPETScPairPtr( out );
    VectorPetsc<T> * vecPetscCast = vecPetscCastPair.first;
    if ( vecPetscCast )
    {
        this->diagonal( *vecPetscCast );
        return;
    }

    CHECK( false ) << "TODO other kind of vector";
}
template<typename T>
std::shared_ptr<Vector<T> >
MatrixPetsc<T>::diagonal() const
{
    std::shared_ptr<Vector<T> > vecRes;
    const MatrixPetscMPI<T>* matPetscMpi = dynamic_cast<const MatrixPetscMPI<T>*> ( this );
    if ( matPetscMpi )
        vecRes.reset( new VectorPetscMPI<T>( this->mapRowPtr() ) );
    else
        vecRes.reset( new VectorPetsc<T>( this->mapRowPtr() ) );
    this->diagonal( *vecRes );

    return vecRes;
}
template<typename T>
void
MatrixPetsc<T>::transpose( MatrixSparse<value_type>& Mt, size_type options ) const
{
    Context ctx( options );
    tic();
    if ( !this->closed() )
        this->close();
    if ( !Mt.closed() )
        Mt.close();
    toc("transpose: close()", FLAGS_v > 0);
    tic();
    MatrixPetsc<T>* Atrans = dynamic_cast<MatrixPetsc<T>*> ( &Mt );
    CHECK( Atrans ) << "support only petsc matrix";
    int ierr = 0;
    if ( ctx.test( MATRIX_TRANSPOSE_ASSEMBLED ) )
    {
        if ( Atrans->isInitialized() )
        {
#if ( PETSC_VERSION_MAJOR >= 3 && PETSC_VERSION_MINOR >= 18 )
            ierr = MatTransposeSetPrecursor( M_mat, Atrans->M_mat );
            CHKERRABORT( this->comm(), ierr );
#endif
            ierr = MatTranspose( M_mat, MAT_REUSE_MATRIX,&Atrans->M_mat );
            CHKERRABORT( this->comm(),ierr );
        }
        else
        {
            ierr = MatTranspose( M_mat, MAT_INITIAL_MATRIX,&Atrans->M_mat );
            CHKERRABORT( this->comm(),ierr );

            Mt.setMapRow( this->mapColPtr() );
            Mt.setMapCol( this->mapRowPtr() );
            MatrixPetscMPI<T>* matTransposedPetscMPI = dynamic_cast<MatrixPetscMPI<T>*> ( &Mt );
            if ( matTransposedPetscMPI )
                matTransposedPetscMPI->initLocalToGlobalMapping();
            Mt.setInitialized( true );
            Mt.setIsClosed( true );
            if ( this->hasGraph() )
                Mt.setGraph( this->graph()->transpose() );
        }
    }
    else if ( ctx.test( MATRIX_TRANSPOSE_UNASSEMBLED ) )
    {
        ierr = MatCreateTranspose( M_mat, &Atrans->M_mat );
        CHKERRABORT( this->comm(),ierr );
    }

    toc("transpose: create mat transpose", FLAGS_v > 0);

    if ( ctx.test( MATRIX_TRANSPOSE_CHECK ) )
    {
        tic();
        if ( this->size1()==this->size2() )
        {
            PetscTruth isSymmetric;
            MatEqual( M_mat, Atrans->M_mat, &isSymmetric );

            if ( isSymmetric )
            {
#if (PETSC_VERSION_MAJOR >= 3)
                MatSetOption( M_mat,MAT_SYMMETRIC,PETSC_TRUE );
#else
                MatSetOption( M_mat,MAT_SYMMETRIC );
#endif
            }

            else
            {
                DVLOG(2) << "[MatrixPETSc::transpose] Petsc matrix is non-symmetric \n";
            }
        }
        toc("transpose: init done", FLAGS_v > 0);
    }

}

template<typename T>
std::shared_ptr<MatrixSparse<T> >
MatrixPetsc<T>::transpose( size_type options ) const
{
    std::shared_ptr<MatrixSparse<T> > matRes;
    const MatrixPetscMPI<T>* matPetscMpi = dynamic_cast<const MatrixPetscMPI<T>*> ( this );
    if ( matPetscMpi )
        matRes = std::make_shared<MatrixPetscMPI<T>>( this->mapColPtr(),this->mapRowPtr(), this->worldCommPtr() );
    else
        matRes = std::make_shared<MatrixPetsc<T>>( this->mapColPtr(),this->mapRowPtr(), this->worldCommPtr() );

    this->transpose( *matRes, options );

    return matRes;
}


template<typename T>
void
MatrixPetsc<T>::symmetricPart( MatrixSparse<value_type>& Mt ) const
{
    this->transpose( Mt, Feel::MATRIX_TRANSPOSE_ASSEMBLED );
    Mt.addMatrix( 1.,*this );
    Mt.scale( 0.5 );
#if 0
    int ierr = 0;
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatSetOption( Mt.M_mat,MAT_SYMMETRIC,PETSC_TRUE );
#else
    ierr = MatSetOption( Mt.M_mat,MAT_SYMMETRIC );
#endif
    CHKERRABORT( this->comm(),ierr );
#endif
    return;

#if 0
    if ( !this->closed() )
        this->close();
    int ierr = 0;

    // first check if the matrix is symmetric
    Mat Atrans;
    MatDuplicate( M_mat,MAT_COPY_VALUES,&Atrans );

#if (PETSC_VERSION_MAJOR >= 3)
    MatTranspose( M_mat, MAT_INITIAL_MATRIX, &Atrans );
#else
    MatTranspose( M_mat, &Atrans );
#endif
    PetscTruth isSymmetric;
    MatEqual( M_mat, Atrans, &isSymmetric );
    PETSc::MatDestroy( Atrans );

    if ( isSymmetric )
    {
        LOG(INFO) << "[MatrixPetsc<T>::symmetricPart] Matrix is already symmetric, don't do anything\n";
#if (PETSC_VERSION_MAJOR >= 3)
        MatSetOption( M_mat,MAT_SYMMETRIC,PETSC_TRUE );
#else
        MatSetOption( M_mat,MAT_SYMMETRIC );
#endif


        // duplicate the matrix
        MatrixPetsc<T>* B = dynamic_cast<MatrixPetsc<T>*> ( &Mt );
        ierr = MatDuplicate( M_mat, MAT_COPY_VALUES, &B->M_mat );
#if (PETSC_VERSION_MAJOR >= 3)
        MatSetOption( B->M_mat,MAT_SYMMETRIC,PETSC_TRUE );
#else
        MatSetOption( B->M_mat,MAT_SYMMETRIC );
#endif
        return ;
    }

    // The matrix was not symmetric, compute 1/2 (A + A^T )
    Mat A[2];
    ierr = MatDuplicate( M_mat,MAT_COPY_VALUES,&A[0] );
    CHKERRABORT( this->comm(),ierr );

    ierr = MatDuplicate( M_mat,MAT_COPY_VALUES,&A[1] );
    CHKERRABORT( this->comm(),ierr );

#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( A[1], MAT_INITIAL_MATRIX, &A[1] );
#else
    MatTranspose( A[1], &A[1] );
#endif
    CHKERRABORT( this->comm(),ierr );

    MatrixPetsc<T>* B = dynamic_cast<MatrixPetsc<T>*> ( &Mt );
    FEELPP_ASSERT( B ).error ( "[symmetricPart] invalid dynamic_cast<MatrixPetsc>" );

    ierr = MatCreateComposite( PETSC_COMM_WORLD,2,A,&B->M_mat );
    CHKERRABORT( this->comm(),ierr );


    //A Completer pour PETSC < 3 !!!!!!!!
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatCompositeSetType( B->M_mat,MAT_COMPOSITE_ADDITIVE );
#endif

    CHKERRABORT( this->comm(),ierr );

    ierr = MatCompositeMerge( B->M_mat );
    CHKERRABORT( this->comm(),ierr );

    ierr = MatScale( B->M_mat, 0.5 );
    CHKERRABORT( this->comm(),ierr );


    // check that B is now symmetric
    Mat Btrans;
#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( B->M_mat, MAT_INITIAL_MATRIX, &Btrans );
#else
    ierr = MatTranspose( B->M_mat, &Btrans );
#endif

    CHKERRABORT( this->comm(),ierr );

    ierr = MatEqual( B->M_mat, Btrans, &isSymmetric );
    CHKERRABORT( this->comm(),ierr );

    if ( isSymmetric )
    {
        LOG(INFO) << "[MatrixPetsc<T>::symmetricPart] symmetric part is really symmetric\n";
#if (PETSC_VERSION_MAJOR >= 3)
        ierr = MatSetOption( B->M_mat,MAT_SYMMETRIC,PETSC_TRUE );
#else
        ierr = MatSetOption( B->M_mat,MAT_SYMMETRIC );
#endif
        CHKERRABORT( this->comm(),ierr );

    }

    else
    {
        PetscPrintf( PETSC_COMM_WORLD,"Warning: Petsc matrix is non-symmetric \n" );
    }

    ierr = PETSc::MatDestroy ( Btrans );
    CHKERRABORT( this->comm(),ierr );
#endif
}


template<typename T>
typename MatrixPetsc<T>::real_type
MatrixPetsc<T>::energy( Vector<value_type> const& __v,
                        Vector<value_type> const& __u,
                        bool transpose ) const
{
    if ( !this->closed() )
        this->close();
    if ( !__v.closed() )
        const_cast<Vector<value_type>*>( &__v )->close();
    if ( !__u.closed() )
        const_cast<Vector<value_type>*>( &__u )->close();


    // case all vector are petsc
    VectorPetsc<T> const* vec_petsc_v = dynamic_cast<VectorPetsc<T> const*>( &__v );
    VectorPetsc<T> const* vec_petsc_u = dynamic_cast<VectorPetsc<T> const*>( &__u );
    if ( vec_petsc_v && vec_petsc_u )
    {
        int ierr = 0;
        std::shared_ptr<VectorPetsc<value_type> > z;
        if ( this->comm().size() > 1 )
            z = std::make_shared< VectorPetscMPI<value_type> >( __v.mapPtr() );
        else
            z = std::make_shared< VectorPetsc<value_type> >( __v.mapPtr() );

        if ( !transpose )
            ierr = MatMult( M_mat, vec_petsc_u->vec(), z->vec() );
        else
            ierr = MatMultTranspose( M_mat, vec_petsc_u->vec(), z->vec() );
        CHKERRABORT( this->comm(),ierr );

        PetscScalar e;
        ierr = VecDot( vec_petsc_v->vec(), z->vec(), &e );
        CHKERRABORT( this->comm(),ierr );

        return e;
    }

    // others vector type
    auto vecPetscCastPair_v = Feel::detail::toPETScPairPtr(  __v, true );
    VectorPetsc<T> const* vecPetscCast_v = vecPetscCastPair_v.first;
    auto vecPetscCastPair_u = Feel::detail::toPETScPairPtr(  __u, true );
    VectorPetsc<T> const* vecPetscCast_u = vecPetscCastPair_u.first;
    CHECK( vecPetscCast_v && vecPetscCast_u ) << "invalid vector type";
    return this->energy( *vecPetscCast_v, *vecPetscCast_u, transpose );
}

template<typename T>
void
MatrixPetsc<T>::updateBlockMat( std::shared_ptr<MatrixSparse<T> > const& m, std::vector<size_type> const& start_i, std::vector<size_type> const& start_j )
{
    if ( !m->closed() )
        m->close();

    auto const& mapRowBlock = m->mapRow();
    auto const& mapColBlock = m->mapCol();
    this->setIsClosed( false );

    auto blockMatrix = std::dynamic_pointer_cast< MatrixPetsc<T> >( m );
    CHECK( blockMatrix ) << "support only PetscMatrix";

    const size_type firstDofGcRowBlock = mapRowBlock.firstDofGlobalCluster();
    const size_type firstDofGcColBlock = mapColBlock.firstDofGlobalCluster();
    const size_type nLocalDofWithoutGhostRowBlock = mapRowBlock.nLocalDofWithoutGhost();
    const rank_type myrank = mapRowBlock.worldComm().globalRank();
    int ierr = 0;
    std::vector<PetscInt> idcXShift;
    auto const& gp2gcRowBlock = mapRowBlock.mapGlobalProcessToGlobalCluster();
    for ( size_type k=0;k<nLocalDofWithoutGhostRowBlock;++k )
    {
        const PetscInt gDof = gp2gcRowBlock[k];
        PetscInt ncolsX;
        const PetscInt *idcX;
        const PetscScalar *valX;
        ierr = MatGetRow( blockMatrix->mat(), gDof, &ncolsX, &idcX, &valX );
        CHKERRABORT( this->comm(),ierr );

        const PetscInt gDofShift = start_i[myrank]+ (gDof-firstDofGcRowBlock);
        idcXShift.resize(ncolsX,0);
        for (int c=0;c<ncolsX;++c)
        {
            if ( mapColBlock.dofGlobalClusterIsOnProc(idcX[c]) )
            {
                idcXShift[c] = start_j[myrank]+ ( idcX[c]-firstDofGcColBlock );
            }
            else
            {
                const rank_type realproc = mapColBlock.procOnGlobalCluster(idcX[c]);
                idcXShift[c] = start_j[realproc] + ( idcX[c]-mapColBlock.firstDofGlobalCluster(realproc) );
            }
        }
        ierr = MatSetValues( this->mat(),1, &gDofShift,ncolsX,idcXShift.data(),valX, INSERT_VALUES );
        CHKERRABORT( this->comm(),ierr );

        // apply this when finish with MatGetRow
        ierr = MatRestoreRow( blockMatrix->mat(), gDof, &ncolsX, &idcX, &valX );
        CHKERRABORT( this->comm(),ierr );

    }

}

template<typename T>
bool
MatrixPetsc<T>::isSymmetric( bool check ) const
{
    PetscBool b = super::isSymmetric()?PETSC_TRUE:PETSC_FALSE;
    if ( check )
        MatIsSymmetric( M_mat, 1e-13, &b );
    return b;
}

template<typename T>
bool
MatrixPetsc<T>::isTransposeOf ( MatrixSparse<T> &Trans ) const
{
    FEELPP_ASSERT( this->size2() == Trans.size1() )( this->size2() )( Trans.size1() ).error( "incompatible dimension" );

    MatrixPetsc<T>* In = dynamic_cast<MatrixPetsc<T>*> ( &Trans );

    Mat Atrans;
    PetscTruth isSymmetric;
    int ierr = 0;

#if (PETSC_VERSION_MAJOR >= 3)
    ierr = MatTranspose( M_mat, MAT_INITIAL_MATRIX, &Atrans );
#else
    ierr = MatTranspose( M_mat, &Atrans );
#endif
    CHKERRABORT( this->comm(),ierr );

    ierr = MatEqual( In->M_mat, Atrans, &isSymmetric );
    CHKERRABORT( this->comm(),ierr );

    ierr = PETSc::MatDestroy ( Atrans );
    CHKERRABORT( this->comm(),ierr );

    return isSymmetric;
}

template <typename T>
inline
void MatrixPetsc<T>::zeroEntriesDiagonal()
{
#if 1
#if 0 //trop lent!
    int start=0, stop=0, ierr=0;

    ierr = MatGetOwnershipRange( M_mat, &start, &stop );
    CHKERRABORT( this->comm(),ierr );

    //VectorPetsc<value_type> diag( std::min(this->size1(),this->size2()), stop-start );
    //VectorPetsc<value_type> diag( std::max(this->size1(),this->size2()), stop-start );
    VectorPetsc<value_type> diag( this->size1(), stop-start );
    //VectorPetsc<value_type> diag( this->size1(), this->size1() );

    //MatGetDiagonal( M_mat, diag.vec() );
    //diag.zero();

    //std:: cout << diag.size();
    //MatZeroRows( M_mat, rows.size(), rows.data(), 1.0);
    ierr = MatDiagonalSet( M_mat, diag.vec(), INSERT_VALUES );
    CHKERRABORT( this->comm(),ierr );
#else
    //std::cout << "\n this->size1(),this->size2() " << this->size1() << " " << this->size2() << std::endl;
#if 0

    for ( uint i = 0; i <std::min( this->size1(),this->size2() ); ++i )
        this->set( i,i,0. );

#else

    //std::cout << "zeroEntriesDiagonal()"<< std::endl;
    for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
    {
        size_type index = it->first;
        std::cout << "\n index " << index;

        if ( index < std::min( this->size1(),this->size2() ) )
        {
            //if (!it->second.template get<2>().empty())
            if ( it->second.template get<2>().find( index )!=it->second.template get<2>().end()  )
                this->set( index,index,0. );
        }
    }

#endif
#endif
#endif
}

template <typename T>
inline
void MatrixPetsc<T>::getMatInfo(std::vector<double> &vec)
{
    MatGetInfo(this->M_mat,MAT_GLOBAL_SUM,&M_info);
    vec.push_back(M_info.block_size);
    vec.push_back(M_info.nz_allocated);
    vec.push_back(M_info.nz_used);
    vec.push_back(M_info.nz_unneeded);
    vec.push_back(M_info.memory);
    vec.push_back(M_info.assemblies);
    vec.push_back(M_info.mallocs);
    vec.push_back(M_info.fill_ratio_given);
    vec.push_back(M_info.fill_ratio_needed);
    vec.push_back(M_info.factor_mallocs);
}
template <typename T>
void MatrixPetsc<T>::threshold(void)
{
    int ierr=0;
    int cpt = 0;
    this->close();
    PetscInt m, n, tot, maxTot;
    ierr = MatGetOwnershipRange(this->M_mat,&m,&n);CHKERRABORT( this->comm(),ierr );
    std::cout << Environment::worldComm().globalRank()  << "\t" << m << " : " << n << std::endl;
    tot = n-m;
    boost::mpi::all_reduce(Environment::worldComm(), tot, maxTot ,  mpi::maximum<int>());
    bool doIt = true;

    for(int i = m ; i < n; i++)
    {
        const PetscScalar *v;
        PetscInt ncols, ncols2;
        const PetscInt *cols;
        PetscInt *cols2;
        PetscScalar *toPut;
        ierr = MatGetRow(this->M_mat, i, &ncols, &cols, &v); CHKERRABORT( this->comm(),ierr );

        ncols2=ncols;
        toPut = new PetscScalar[ncols];
        cols2 = new PetscInt[ncols];
        //if(doIt)
        {
            for(int j = 0; j < ncols ; j++)
            {
                cols2[j] = cols[j];
                if(fabs(v[j]) < 0.8 )
                    toPut[j] = 0.;
                else if(v[j] < 0.)
                    toPut[j] = -1.;
                else
                    toPut[j] =  1.;
            }
        }
        // apply this when finish with MatGetRow
        ierr = MatRestoreRow( this->M_mat, i, &ncols, &cols , &v );CHKERRABORT( this->comm(),ierr );

        //if(doIt)
        {
            ierr = MatSetValues(this->M_mat, 1, &i, ncols2, cols2, toPut, INSERT_VALUES); CHKERRABORT( this->comm(),ierr );
            this->close();
            delete [] toPut;
            delete [] cols2;
        }
        cpt++;
        if( (i == n-1) && (cpt < maxTot))
        {
            doIt = false;
            i--;
        }
    }
    Environment::worldComm().barrier();

    this->close();
}

template <typename T>
void
MatrixPetsc<T>::save( std::string const& filename, std::string const& format )
{
    std::string filenameUsed = boost::str( boost::format("%1%_%2%_%3%") %filename %format %this->comm().globalRank() );
    std::ofstream ofs( filenameUsed );
    if (ofs)
    {
        if ( format=="binary" )
        {
            boost::archive::binary_oarchive oa(ofs);
            oa << *this;
        }
        else if ( format=="xml")
        {
            boost::archive::xml_oarchive oa(ofs);
            oa << boost::serialization::make_nvp("vectorpetsc", *this );
        }
        else if ( format=="text")
        {
            boost::archive::text_oarchive oa(ofs);
            oa << *this;
        }
        else
            CHECK( false ) << "MatrixPetsc save() function : error with unknown format " << format;
    }
    else
    {
        CHECK( false ) << "MatrixPetsc save() function : error opening ofstream with name " << filenameUsed;
    }
    ofs.close();
}

template <typename T>
void
MatrixPetsc<T>::load( std::string const& filename, std::string const& format )
{
    std::string filenameUsed = boost::str( boost::format("%1%_%2%_%3%") %filename %format %this->comm().globalRank() );
    std::ifstream ifs( filenameUsed );
    if ( ifs )
    {
        if ( format=="binary" )
        {
            boost::archive::binary_iarchive ia(ifs);
            ia >> *this;
        }
        else if ( format=="xml")
        {
            boost::archive::xml_iarchive ia(ifs);
            ia >> boost::serialization::make_nvp("vectorpetsc", *this );
        }
        else if ( format=="text")
        {
            boost::archive::text_iarchive ia(ifs);
            ia >> *this;
        }
        else
            CHECK( false ) << "MatrixPetsc save() function : error with unknown format " << format;
    }
    else
    {
        CHECK( false ) << "MatrixPetsc load() function : error opening ofstream with name " << filenameUsed;
    }
    ifs.close();
}

template <typename T>
template<class Archive>
void MatrixPetsc<T>::serialize(Archive & ar, const unsigned int version )
{
    if ( !this->closed() )
        this->close();

    CHECK( this->isInitialized() ) << "TODO : matrix not initialized";
    int ierr=0;
#if 0
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(super);
#endif

    if ( Archive::is_saving::value )
    {
        CHECK( this->hasGraph() ) << "graph is required";

        std::vector<double> matValues;
        matValues.reserve( this->graph()->ja().size() );
        auto itVal = matValues.begin();
        PetscInt rowIds[1];
        std::vector<PetscInt> colIds;
        for ( auto const& rowDataEntries : this->graph()->storage() )
        {
            if ( !this->mapRow().dofGlobalClusterIsOnProc( rowDataEntries.first ) )
                continue;
            rowIds[0] = rowDataEntries.first;
            auto const& colData = boost::get<2>( rowDataEntries.second );
            if ( colData.empty() )
                continue;
            colIds.assign( colData.begin(),colData.end() );
            itVal = matValues.end();
            matValues.resize( matValues.size() + colIds.size() );
            ierr = MatGetValues( this->mat(), 1, rowIds, colIds.size(), colIds.data(), &(*itVal)/*matValues.data()*/ );
            CHKERRABORT( this->comm(),ierr );
        }
        ar & boost::serialization::make_nvp("values", matValues );
    }
    else
    {
        CHECK( this->hasGraph() ) << "graph is required";

        std::vector<double> matValues;
        ar & boost::serialization::make_nvp("values", matValues );
        auto itVal = matValues.begin();
        PetscInt rowIds[1];
        std::vector<PetscInt> colIds;
        for ( auto const& rowDataEntries : this->graph()->storage() )
        {
            if ( !this->mapRow().dofGlobalClusterIsOnProc( rowDataEntries.first ) )
                continue;
            rowIds[0] = rowDataEntries.first;
            auto const& colData = boost::get<2>( rowDataEntries.second );
            if ( colData.empty() )
                continue;
            colIds.assign( colData.begin(),colData.end() );
            ierr = MatSetValues( this->mat(), 1, rowIds, colIds.size(), colIds.data(), &(*itVal), INSERT_VALUES );
            CHKERRABORT( this->comm(),ierr );
            itVal += colIds.size();
        }
        this->close();
    }
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm )
{}
template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol )
    :
    super( dmRow,dmCol )
{}
template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, worldcomm_ptr_t const& worldComm )
    :
    super( dmRow,dmCol,worldComm )
{}

//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
MatrixPetscMPI<T>::MatrixPetscMPI( Mat m, datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol,bool initLocalToGlobalMapping, bool destroyMatOnExit )
    :
    super( m, dmRow, dmCol, destroyMatOnExit )
{
    if ( initLocalToGlobalMapping )
        this->initLocalToGlobalMapping();
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void MatrixPetscMPI<T>::init( const size_type m,
                              const size_type n,
                              const size_type m_l,
                              const size_type n_l,
                              const size_type nnz,
                              const size_type /*noz*/ )
{
    //std::cout << "\n MatrixPetscMPI<T>::init without graph " << std::endl;
}

template <typename T>
void MatrixPetscMPI<T>::init( const size_type /*m*/,
                              const size_type /*n*/,
                              const size_type /*m_l*/,
                              const size_type /*n_l*/,
                              graph_ptrtype const& graph )
{
    //this->comm().globalComm().barrier();
    VLOG(1) << "MatrixPetscMPI<T>::init with graph start on proc"<< this->comm().globalRank() << "("<<this->comm().godRank() <<")" << std::endl;

    this->setGraph( graph );
    this->graph()->close();

    // Clear initialized matrices
    if ( this->isInitialized() )
        this->clear();

    this->setInitialized(  true );
    this->setIsClosed( true );

    int ierr     = 0;
    int m_global = static_cast<int>( this->graph()->mapRow().nDof()/*m*/ );
    int n_global = static_cast<int>( this->graph()->mapCol().nDof()/*n*/ );
    int m_local  = static_cast<int>( this->graph()->mapRow().nLocalDofWithoutGhost()/*m_l*/ );
    int n_local  = static_cast<int>( this->graph()->mapCol().nLocalDofWithoutGhost()/*n_l*/ );
    if ( m_global==0 )
        return;


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

    if ( n_dnzOffProc==0 )
        dnzOffProc = PETSC_IGNORE;

#if PETSC_VERSION_LESS_THAN(3,3,0)
    ierr = MatCreateMPIAIJ ( this->comm(),
                             m_local, n_local,
                             m_global, n_global,
                             /*PETSC_DECIDE*//*n_dnz*/0, /*PETSC_IGNORE*/dnz,
                             /*PETSC_DECIDE*/0/*n_dnzOffProc*/, dnzOffProc,
                             //&M_matttt);
                             &this->M_mat);
                             //&( this->mat() ) ); //(&this->M_mat));
#else
    ierr = MatCreateAIJ ( this->comm(),
                          m_local, n_local,
                          m_global, n_global,
                          /*PETSC_DECIDE*//*n_dnz*/0, /*PETSC_IGNORE*/dnz,
                          /*PETSC_DECIDE*/0/*n_dnzOffProc*/, dnzOffProc,
                          //&M_matttt);
                          &this->M_mat);
    //&( this->mat() ) ); //(&this->M_mat));

#endif
    CHKERRABORT( this->comm(),ierr );

    // free
    delete[] dnz;
    delete[] dnzOffProc;

    std::vector<PetscInt> ia( this->graph()->ia().size() );
    std::vector<PetscInt> ja( this->graph()->ja().size() );
    std::copy( this->graph()->ia().begin(), this->graph()->ia().end(), ia.begin() );
    std::copy( this->graph()->ja().begin(), this->graph()->ja().end(), ja.begin() );
#if 0
    ierr = MatMPIAIJSetPreallocation( this->mat(), 0, dnz, 0, dnzOffProc );
#else
    ierr = MatMPIAIJSetPreallocationCSR( this->mat(), ia.data() , ja.data(), this->graph()->a().data() );
    //ierr = MatMPIAIJSetPreallocationCSR( this->mat(), ia.data() , ja.data(),NULL );
#endif
    CHKERRABORT( this->comm(),ierr );

    //----------------------------------------------------------------------------------//
    // localToGlobal mapping
    this->initLocalToGlobalMapping();

    //----------------------------------------------------------------------------------//
    // options

    ierr = MatSetFromOptions( this->mat() );
    CHKERRABORT( this->comm(),ierr );

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 0)
    //MatSetOption(this->mat(),MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
    MatSetOption( this->mat(),MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
#elif (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR == 0)
    MatSetOption( this->mat(),MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( this->mat(),MAT_KEEP_ZEROED_ROWS );
#endif



#if 1
    // additional insertions will not be allowed if they generate
    // a new nonzero
    //ierr = MatSetOption (M_mat, MAT_NO_NEW_NONZERO_LOCATIONS);
    //CHKERRABORT(this->comm(),ierr);

    // generates an error for new matrix entry
    ierr = MatSetOption ( this->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
    CHKERRABORT( this->comm(),ierr );
#endif

    //this->zero();
    //this->zeroEntriesDiagonal();
    //this->printMatlab("NULL");
    //std::cout << "MatrixPetscMPI<T>::init with graph finish on proc"
    //          << this->comm().globalRank() << "("<<this->comm().godRank() <<")" << std::endl;
    //this->comm().globalComm().barrier();
}


//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::initLocalToGlobalMapping()
{
    int ierr = 0;
    IS isRow;
    IS isCol;
    ISLocalToGlobalMapping isLocToGlobMapRow;
    ISLocalToGlobalMapping isLocToGlobMapCol;

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
    ierr = ISCreateGeneral( this->comm(), n_idxRow, idxRow, PETSC_COPY_VALUES, &isRow );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISCreateGeneral( this->comm(), n_idxCol, idxCol, PETSC_COPY_VALUES, &isCol );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISCreateGeneral( this->comm(), n_idxRow, idxRow, &isRow );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISCreateGeneral( this->comm(), n_idxCol, idxCol, &isCol );
    CHKERRABORT( this->comm(),ierr );
#endif

    ierr=ISLocalToGlobalMappingCreateIS( isRow, &isLocToGlobMapRow );
    CHKERRABORT( this->comm(),ierr );

    ierr=ISLocalToGlobalMappingCreateIS( isCol, &isLocToGlobMapCol );
    CHKERRABORT( this->comm(),ierr );

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = MatSetLocalToGlobalMapping( this->mat(),isLocToGlobMapRow,isLocToGlobMapCol );
#else
    ierr = MatSetLocalToGlobalMapping( this->mat(),isLocToGlobMapRow );
#endif
    CHKERRABORT( this->comm(),ierr );

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy( &isRow );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISDestroy( &isCol );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISLocalToGlobalMappingDestroy( &isLocToGlobMapRow );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISLocalToGlobalMappingDestroy( &isLocToGlobMapCol );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy( isRow );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISDestroy( isCol );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISLocalToGlobalMappingDestroy( isLocToGlobMapRow );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISLocalToGlobalMappingDestroy( isLocToGlobMapCol );
    CHKERRABORT( this->comm(),ierr );
#endif
    delete[] idxRow;
    delete[] idxCol;

}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
void MatrixPetscMPI<T>::set( const size_type i,
                             const size_type j,
                             const value_type& value )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );;

    int ierr=0, i_val=i, j_val=j;

    PetscScalar petsc_value = static_cast<PetscScalar>( value );
    ierr = MatSetValuesLocal( this->mat(), 1, &i_val, 1, &j_val,
                              &petsc_value, INSERT_VALUES );

    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
inline
void MatrixPetscMPI<T>::add ( const size_type i,
                              const size_type j,
                              const value_type& value )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "MatrixPetsc<> not properly initialized" );

    int ierr=0, i_val=i, j_val=j;
    //DVLOG(2) << "[MatrixPetsc<>::add] adding value " << value << " at (" << i << "," << j << ")\n";
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr = MatSetValuesLocal( this->mat(), 1, &i_val, 1, &j_val,
                              &petsc_value, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}


//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::addMatrix( const ublas::matrix<value_type>& dm,
                              const std::vector<size_type>& rows,
                              const std::vector<size_type>& cols )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    const size_type m = dm.size1();
    const size_type n = dm.size2();

    FEELPP_ASSERT ( rows.size() == this->size1() ).error( "invalid row size" );
    FEELPP_ASSERT ( cols.size() == this->size2() ).error( "invalid column size" );

    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValuesLocal( this->mat(),
                              m, ( int* ) boost::addressof( rows[0] ),
                              n, ( int* ) boost::addressof( cols[0] ),
                              ( PetscScalar* ) dm.data().begin(),
                              ADD_VALUES );


    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::addMatrix( int* rows, int nrows,
                              int* cols, int ncols,
                              value_type* data,
                              size_type K,
                              size_type K2)
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );
    /*for (int k=0;k<nrows;++k)
        for (int q=0;q<ncols;++q)
            { std::cout << "rank " << this->comm().rank()
                        << "add mat("<< rows[k] <<"," << cols[q]
                        << ") = "//<< data[k][q]
                        << " " << this->mapCol().mapGlobalProcessToGlobalCluster()[cols[q]]
                        << std::endl; }*/
    int ierr=0;

    // These casts are required for PETSc <= 2.1.5
    ierr = MatSetValuesLocal( this->mat(),
                              nrows, ( int* ) rows,
                              ncols, ( int* ) cols,
                              ( PetscScalar* ) data,
                              ADD_VALUES );

    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//
#if 0
template <typename T>
inline
void
MatrixPetscMPI<T>::addMatrix( const T a_in, MatrixSparse<T> const&X_in )
{
#if 0
    if (this->hasGraph() && X_in.hasGraph() &&
        static_cast<void*>( const_cast<graph_type*>(this->graph().get()) ) == static_cast<void*>( const_cast<graph_type*>(X_in.graph().get())) )
        {
            this->addMatrixSameNonZeroPattern(a_in,X_in);
        }
    else
        {
            super::addMatrix(a_in,X_in);
        }
#else
    super::addMatrix(a_in,X_in);
#endif
}
#endif
//----------------------------------------------------------------------------------------------------//

template <typename T>
inline
void
MatrixPetscMPI<T>::addMatrixSameNonZeroPattern( const T a_in, MatrixSparse<T> &X_in )
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "petsc matrix not initialized" );

    // sanity check. but this cannot avoid
    // crash due to incompatible sparsity structure...
    FEELPP_ASSERT( this->size1() == X_in.size1() )( this->size1() )( X_in.size1() ).error( "incompatible dimension" );
    FEELPP_ASSERT( this->size2() == X_in.size2() )( this->size2() )( X_in.size2() ).error( "incompatible dimension" );

    PetscScalar a = static_cast<PetscScalar>( a_in );
    MatrixPetscMPI<T>* X = dynamic_cast<MatrixPetscMPI<T>*> ( &X_in );
    FEELPP_ASSERT ( X != 0 ).error( "invalid petsc matrix" );

    int ierr=0;

    // the matrix from which we copy the values has to be assembled/closed
    X->close ();

    const size_type nLocalDofWithGhost = this->mapRow().nLocalDofWithGhost();

    if (a==1)
        {
            //std::cout << "case a==1 " << std::endl;
            for ( size_type k=0;k<nLocalDofWithGhost;++k )
                {
                    if (!this->mapRow().dofGlobalProcessIsGhost(k))
                        {
                            const PetscInt gDof = X->mapRow().mapGlobalProcessToGlobalCluster(k);
                            PetscInt ncolsX;
                            const PetscInt *idcX;
                            const PetscScalar *valX;
                            ierr = MatGetRow( X->mat(), gDof, &ncolsX, &idcX, &valX );
                            CHKERRABORT( this->comm(),ierr );

                            // set new values
                            ierr = MatSetValues(this->mat(),1, &gDof,ncolsX,idcX,valX,ADD_VALUES);
                            CHKERRABORT( this->comm(),ierr );

                            // apply this when finish with MatGetRow
                            ierr = MatRestoreRow( X->mat(), gDof, &ncolsX, &idcX, &valX );
                            CHKERRABORT( this->comm(),ierr );
                        }
                }
        }
    else // case a!=1
        {
            for ( size_type k=0;k<nLocalDofWithGhost;++k )
                {
                    if (!this->mapRow().dofGlobalProcessIsGhost(k))
                        {
                            const PetscInt gDof = X->mapRow().mapGlobalProcessToGlobalCluster(k);
                            PetscInt ncolsX;
                            const PetscInt *idcX;
                            const PetscScalar *valX;
                            ierr = MatGetRow( X->mat(), gDof, &ncolsX, &idcX, &valX );
                            CHKERRABORT( this->comm(),ierr );

                            //get new values in row
                            PetscScalar *valNewRow = new PetscScalar[ncolsX];
                            for (int col=0;col<ncolsX;++col)
                                valNewRow[col]=a*valX[col];

                            // set new values
                            ierr = MatSetValues(this->mat(),1, &gDof,ncolsX,idcX,valNewRow,ADD_VALUES);
                            CHKERRABORT( this->comm(),ierr );

                            // apply this when finish with MatGetRow
                            ierr = MatRestoreRow( X->mat(), gDof, &ncolsX, &idcX, &valX );
                            CHKERRABORT( this->comm(),ierr );
                            // clean
                            delete[] valNewRow;
                        }
                }

        } // case : a!=1

    this->close();

} // addMatrixSameNonZeroPattern


//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::zero()
{
    CHECK ( this->isInitialized() ) <<  "petsc matrix not properly initialized";

    if ( !this->closed() )
        this->close();

    int ierr=0;
    ierr = MatZeroEntries( this->mat() );
    CHKERRABORT( this->comm(),ierr );

#if 0
    int ierr=0;
    LOG(INFO) << "assembly status of the matrix";
    PetscBool is_assembled;
    MatAssembled( this->mat(), &is_assembled );
    VLOG(2) << "Matrix is assembled : " << (is_assembled?"true":"false");
    if ( is_assembled )
    {
        ierr = MatZeroEntries( this->mat() );
        CHKERRABORT( this->comm(),ierr );
    }

    else
    {
        if ( this->graph() )
        {
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                if ( ( int )it->second.template get<0>() == this->comm().globalRank() || this->mapRow().mapGlobalClusterToGlobalProcess().size()==0 )
                {

                    // Work in progress (but normally this part is useless because we use the CSR prealocation)
#if 0
                    std::vector<PetscInt> cols(  it->second.template get<2>().size(), 0 );
                    //PetscInt row = it->second.template get<1>();

                    PetscInt row = 0;
                    if (this->mapRow().mapGlobalClusterToGlobalProcess().size()==0)
                        row=0;//this->mapRow().mapGlobalClusterToGlobalProcess()[it->first-this->mapRow().firstDofGlobalCluster()];
                    else
                        row=it->first;//this->mapRow().mapGlobalClusterToGlobalProcess()[it->first-this->mapRow().firstDofGlobalCluster()];

                    //MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE);

                    //this->mapRow()->firstDofGlobalCluster( this->comm().rank() );
                    auto it2=it->second.template get<2>().begin();
                    auto  en2=it->second.template get<2>().end();
                    size_type cpt=0;

                    for ( ; it2!=en2 ; ++it2 )
                    {
                        //if ( ( *it2 >= this->graph()->firstColEntryOnProc() ) &&
                        //( *it2 <= this->graph()->lastColEntryOnProc() ) )
                        {
                            //cols[cpt] = this->mapRow().mapGlobalClusterToGlobalProcess()[*it2];
                            if (this->mapCol().mapGlobalClusterToGlobalProcess().size()>0)
                                cols[cpt] = *it2;//this->mapCol().mapGlobalClusterToGlobalProcess()[*it2-this->mapCol().firstDofGlobalCluster() ];
                            else cols[cpt] = 0;
                            ++cpt;
                        }
                    }

                    if (cpt>0)
                        {
                            //if (this->mapRow().mapGlobalClusterToGlobalProcess().size()>0) {
                            cols.resize( cpt );
                            std::vector<PetscScalar> v(  cpt,0 );
                            //for (int k=0;k<cpt;++k)
                            //    std::cout<< "zero : row " << row << " col" << cols[k] << std::endl;
                            //std::copy( it->second.template get<2>().begin(), it->second.template get<2>().end(), cols.begin() );
                            //MatSetValuesLocal( this->mat(), 1, &row, it->second.template get<2>().size(), cols.data(), v.data(), INSERT_VALUES );

                            //MatSetValuesLocal( this->mat(), 1, &row, cpt, cols.data(), v.data(), INSERT_VALUES );
                            MatSetValues( this->mat(), 1, &row, cpt, cols.data(), v.data(), INSERT_VALUES );
                        }
#else
                    std::vector<PetscInt> cols(  it->second.template get<2>().size(), 0 );
                    //PetscInt row = it->second.template get<1>();

                    PetscInt row = 0;
                    if (this->mapRow().mapGlobalClusterToGlobalProcess().size()==0)
                        row=0;//this->mapRow().mapGlobalClusterToGlobalProcess()[it->first-this->mapRow().firstDofGlobalCluster()];
                    else
                        row=it->first;//this->mapRow().mapGlobalClusterToGlobalProcess()[it->first-this->mapRow().firstDofGlobalCluster()];

                    //MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE);

                    //this->mapRow()->firstDofGlobalCluster( this->comm().rank() );
                    auto it2=it->second.template get<2>().begin();
                    auto  en2=it->second.template get<2>().end();
                    size_type cpt=0;

                    for ( ; it2!=en2 ; ++it2 )
                    {
                        //if ( ( *it2 >= this->graph()->firstColEntryOnProc() ) &&
                        //( *it2 <= this->graph()->lastColEntryOnProc() ) )
                        {
                            //cols[cpt] = this->mapRow().mapGlobalClusterToGlobalProcess()[*it2];
                            if (this->mapCol().mapGlobalClusterToGlobalProcess().size()>0)
                                cols[cpt] = this->mapCol().mapGlobalClusterToGlobalProcess()[*it2-this->mapCol().firstDofGlobalCluster() ];
                            else cols[cpt] = 0;
                            ++cpt;
                        }
                    }

                    if (cpt>0)
                        {
                            //if (this->mapRow().mapGlobalClusterToGlobalProcess().size()>0) {
                            cols.resize( cpt );
                            std::vector<PetscScalar> v(  cpt,0 );
                            for (int k=0;k<cpt;++k)
                                std::cout<< "zero : row " << row << " col" << cols[k] << std::endl;
                            //std::copy( it->second.template get<2>().begin(), it->second.template get<2>().end(), cols.begin() );
                            //MatSetValuesLocal( this->mat(), 1, &row, it->second.template get<2>().size(), cols.data(), v.data(), INSERT_VALUES );

                            MatSetValuesLocal( this->mat(), 1, &row, cpt, cols.data(), v.data(), INSERT_VALUES );
                            //MatSetValues( this->mat(), 1, &row, cpt, cols.data(), v.data(), INSERT_VALUES );
                        }

#endif
                }
            }

            //this->close();
        }
    }
#endif
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
MatrixPetscMPI<T>::zero( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{
    this->zero();
#if 0
    CHECK ( this->isInitialized() ) <<  "petsc matrix not properly initialized";

    int ierr=0;

    PetscBool is_assembled;
    MatAssembled( this->mat(), &is_assembled );

    if ( is_assembled )
    {
        ierr = MatZeroEntries( this->mat() );
        CHKERRABORT( this->comm(),ierr );

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
                PetscInt row = it->second.template get<1>();
                std::copy( it->second.template get<2>().begin(), it->second.template get<2>().end(), cols.begin() );
                MatSetValuesLocal( this->mat(), 1, &row, it->second.template get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
            }

#else

            //std::vector<PetscInt> cols( this->graph()->nCols(), 0 );
            //std::vector<PetscScalar> v( this->graph()->nCols(), 0. );
            for ( auto it=this->graph()->begin(), en=this->graph()->end() ; it!=en ; ++it )
            {
                if ( ( int )it->second.template get<0>() == this->comm().rank() )
                {

                    std::vector<PetscInt> cols(  it->second.template get<2>().size(), 0 );
                    //std::set<PetscInt> cols;

                    //PetscInt row = it->second.template get<1>();
                    PetscInt row =  it->first;


                    //MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE);

                    //this->mapRow()->firstDofGlobalCluster( this->comm().rank() );
                    auto it2=it->second.template get<2>().begin();
                    auto  en2=it->second.template get<2>().end();
                    size_type cpt=0;

                    for ( ; it2!=en2 ; ++it2 )
                    {
                        if ( ( *it2 >= this->graph()->firstColEntryOnProc() ) &&
                                ( *it2 <= this->graph()->lastColEntryOnProc() ) )
                        {
                            cols[cpt] = this->mapRow().mapGlobalClusterToGlobalProcess()[*it2];
                            ++cpt;
                        }
                    }

                    cols.resize( cpt );
                    std::vector<PetscScalar> v(  cpt,0 );

                    //std::copy( it->second.template get<2>().begin(), it->second.template get<2>().end(), cols.begin() );
                    MatSetValuesLocal( this->mat(), 1, &row, it->second.template get<2>().size(), cols.data(), v.data(), INSERT_VALUES );
                }
            }

#endif
        }
    }
#endif
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
void
MatrixPetscMPI<T>::zeroRows( std::vector<int> const& rows,
                             Vector<value_type> const& values,
                             Vector<value_type>& rhs,
                             Context const& on_context,
                             value_type value_on_diagonal )
{

    // specific treatment if not all process going here
    bool hasAllProcess = true;
    if ( hasAllProcess )
    {
        if ( !this->closed() )
            this->close();
        if ( !rhs.closed() )
            rhs.close();
    }
    else
    {
#if !PETSC_VERSION_LESS_THAN(3,2,0)
        MatSetOption( this->mat(),MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE );
#endif
    }


#if PETSC_VERSION_LESS_THAN(3,0,0)
    MatSetOption( this->mat(),MAT_KEEP_ZEROED_ROWS );
#elif PETSC_VERSION_LESS_THAN(3,1,0)
    MatSetOption( this->mat(),MAT_KEEP_ZEROED_ROWS,PETSC_TRUE );
#else
    MatSetOption( this->mat(),MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE );
#endif

    int start=0, stop=this->mapRow().nLocalDofWithGhost(), ierr=0;
    VectorPetscMPI<T>* prhs = dynamic_cast<VectorPetscMPI<T>*> ( &rhs );
    CHECK( prhs != 0 ) << "dynamic cast from Vector to VectorPetscMPI failed for rhs";
    const VectorPetscMPI<T>* pvalues = dynamic_cast<const VectorPetscMPI<T>*> ( &values );
    CHECK( pvalues != 0 ) << "dynamic cast from Vector to VectorPetscMPI failed for values";

    if ( on_context.test( ContextOn::ELIMINATION ) )
    {
        rhs.setIsClosed( false );
        std::shared_ptr<VectorPetscMPI<value_type>> diag;
        if ( on_context.test( ContextOn::KEEP_DIAGONAL ) )
        {
            diag = std::make_shared<VectorPetscMPI<value_type>>( this->mapColPtr() );
            MatGetDiagonal( this->M_mat, diag->vec() );
            diag->close();
        }

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
        if ( on_context.test(ContextOn::SYMMETRIC ) )
        {
            MatZeroRowsColumnsLocal(this->M_mat, rows.size(), rows.data(), value_on_diagonal, pvalues->vec(), prhs->vec() );
        }
        else
        {
            MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), value_on_diagonal, pvalues->vec(), prhs->vec() );
            //CHKERRABORT(this->comm(),ierr);
        }
        if ( on_context.test( ContextOn::KEEP_DIAGONAL ) )
        {
            MatDiagonalSet( this->M_mat, diag->vec(), INSERT_VALUES );
            for ( size_type i = 0; i < rows.size(); ++i )
            {
                // warning: a row index may belong to another
                // processor, so make sure that we access only the
                // rows that belong to this processor
                if ( rows[i] >= start && rows[i] < stop )
                    rhs.set( rows[i], values(rows[i])*diag->operator()( rows[i] ) );
            }
        }

#else
        MatZeroRowsLocal( this->mat(), rows.size(), rows.data(), value_on_diagonal );

        for ( size_type i = 0; i < rows.size(); ++i )
        {
            // eliminate column

            // warning: a row index may belong to another
            // processor, so make sure that we access only the
            // rows that belong to this processor
            if ( rows[i] >= start && rows[i] < stop )
                rhs.set( rows[i], values(rows[i]) );
        }
#endif
    }
    // rsh doesn't be closed because not all processors are present here with composite spaces(this call must be done after)
    if ( hasAllProcess )
    {
        if ( !rhs.closed() )
            rhs.close();
    }
    else
    {
        //reset MatOption (assemble with communication)
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 1)
        MatSetOption( this->mat(),MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_FALSE );
#else
        // ???
#endif
    }

} // zeroRows

//----------------------------------------------------------------------------------------------------//

template class MatrixPetsc<double>;
template class MatrixPetscMPI<double>;

template void MatrixPetsc<double>::serialize(boost::archive::text_oarchive & ar, const unsigned int version );
template void MatrixPetsc<double>::serialize(boost::archive::text_iarchive & ar, const unsigned int version );
template void MatrixPetsc<double>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version );
template void MatrixPetsc<double>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version );
template void MatrixPetsc<double>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version );
template void MatrixPetsc<double>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version );

} // Feel



//BOOST_CLASS_EXPORT_IMPLEMENT( Feel::MatrixPetsc<double> )
//BOOST_CLASS_EXPORT_IMPLEMENT( Feel::MatrixPetscMPI<double> )


#endif // FEELPP_HAS_PETSC_H
