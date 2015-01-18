/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-02

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file vectorpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-02
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>

#if defined( FEELPP_HAS_PETSC_H )

extern "C"
{
#include "petscsys.h"
}



namespace Feel
{
template <typename T>
typename VectorPetsc<T>::clone_ptrtype
VectorPetsc<T>::clone () const
{
    clone_ptrtype cloned_vector ( new VectorPetsc<T>( this->mapPtr() ) );
    CHECK( cloned_vector->size() == this->size() ) << "Invalid cloned vector size : " << cloned_vector->size()
                                                   << " expected size : " << this->size() ;
    //*cloned_vector = *this;
    CHECK( this->closed() ) << "VectorPETSc is closed and should not";
    return cloned_vector;
}

template <typename T>
inline
void
VectorPetsc<T>::init ( const size_type n,
                       const size_type n_local,
                       const bool fast )
{
    int ierr=0;
    int petsc_n=static_cast<int>( n );
    int petsc_n_local=static_cast<int>( n_local );


    // Clear initialized vectors
    if ( this->isInitialized() )
        this->clear();


    // create a sequential vector if on only 1 processor
    if ( n_local == n )
    {
        ierr = VecCreateSeq ( PETSC_COMM_SELF, petsc_n, &M_vec );
        CHKERRABORT( PETSC_COMM_SELF,ierr );

        ierr = VecSetFromOptions ( M_vec );
        CHKERRABORT( PETSC_COMM_SELF,ierr );
    }

    // otherwise create an MPI-enabled vector
    else
    {
        DCHECK( n_local < n ) << "invalid local size : " << n_local << " is not less than global size " <<  n;

        ierr = VecCreateMPI ( this->comm(), petsc_n_local, petsc_n,
                              &M_vec );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecSetFromOptions ( M_vec );
        CHKERRABORT( this->comm(),ierr );
    }

    this->M_is_initialized = true;


    if ( fast == false )
        this->zero ();
}
template <typename T>
void
VectorPetsc<T>::set ( const value_type& value )
{
    int ierr=0;
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr = VecSet ( M_vec, petsc_value );
    CHKERRABORT( this->comm(),ierr );
}
template <typename T>
void
VectorPetsc<T>::set ( size_type i, const value_type& value )
{
    DCHECK( i<size() ) << "invalid index " << i <<  " size : " << size();


    int ierr=0;
    int i_val = static_cast<int>( i );
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr = VecSetValues ( M_vec, 1, &i_val, &petsc_value, INSERT_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetsc<T>::setVector ( int* i, int n, value_type* v )
{
    //FEELPP_ASSERT(n<=size())( n )( size() ).error( "invalid local index array size" );

#if PETSC_VERSION_LESS_THAN(3,5,3)
    if ( n == 0 ) return;
#endif
    int ierr=0;
    ierr = VecSetValues ( M_vec, n, i, v, INSERT_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetsc<T>::add ( const size_type i, const value_type& value )
{
    DCHECK( i<size() ) << "invalid index " << i <<  " size : " << size();

    int ierr=0;
    int i_val = static_cast<int>( i );
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr = VecSetValues ( M_vec, 1, &i_val, &petsc_value, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetsc<T>::addVector ( int* i, int n, value_type* v )
{
    //FEELPP_ASSERT(n<=size())( n )( size() ).error( "invalid local index array size" );

#if PETSC_VERSION_LESS_THAN(3,5,3)
    if ( n == 0 ) return;
#endif
    int ierr=0;
    ierr = VecSetValues ( M_vec, n, i, v, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );

}
template <typename T>
typename VectorPetsc<T>::value_type
VectorPetsc<T>::operator() ( const size_type i ) const
    {
        DCHECK( this->isInitialized() ) << "VectorPETSc not initialized";
        DCHECK ( ( ( i >= this->firstLocalIndex() ) &&
                   ( i <  this->lastLocalIndex() ) ) ) << "invalid vector index " <<  i
                                                       << " first local index: "  << this->firstLocalIndex()
                                                       << " last local index:  " << this->lastLocalIndex();

        int ierr=0;
        PetscScalar *values, value=0.;


        ierr = VecGetArray( M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        value = values[i - this->firstLocalIndex()];

        ierr = VecRestoreArray ( M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        return static_cast<value_type>( value );
    }
template <typename T>
typename VectorPetsc<T>::value_type&
VectorPetsc<T>::operator() ( const size_type i )
    {
        DCHECK( this->isInitialized() ) << "VectorPETSc not initialized";
        DCHECK ( ( ( i >= this->firstLocalIndex() ) &&
                   ( i <  this->lastLocalIndex() ) ) ) << "invalid vector index " <<  i
                                                       << " first local index: "  << this->firstLocalIndex()
                                                       << " last local index:  " << this->lastLocalIndex();

        int ierr=0;
        PetscScalar *values;

        ierr = VecGetArray( M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        PetscScalar& value = values[i - this->firstLocalIndex()];

        ierr = VecRestoreArray ( M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        return static_cast<value_type&>( value );
    }

template <typename T>
void
VectorPetsc<T>::pointwiseMult ( Vector<T> const& xx, Vector<T> const& yy )
{    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( this->comm().size()>1 )
    {
        const_cast<VectorPetscMPI<T>*>( dynamic_cast<const VectorPetscMPI<T>*>( &xx ) )->close();
        const_cast<VectorPetscMPI<T>*>( dynamic_cast<const VectorPetscMPI<T>*>( &yy ) )->close();
    }
    else
    {
        const_cast<VectorPetsc<T>*>( dynamic_cast<const VectorPetsc<T>*>( &xx ) )->close();
        const_cast<VectorPetsc<T>*>( dynamic_cast<const VectorPetsc<T>*>( &yy ) )->close();
    }

    const VectorPetsc<T>* x = dynamic_cast<const VectorPetsc<T>*>( &xx );
    const VectorPetsc<T>* y = dynamic_cast<const VectorPetsc<T>*>( &yy );
    CHECK( x != 0 && y != 0 ) << "invalid Vector<> types";
    auto ierr =  VecPointwiseMult(this->vec(), x->vec(), y->vec() );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetsc<T>::pointwiseDivide ( Vector<T> const& xx, Vector<T> const& yy )
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( this->comm().size()>1 )
    {
        const_cast<VectorPetscMPI<T>*>( dynamic_cast<const VectorPetscMPI<T>*>( &xx ) )->close();
        const_cast<VectorPetscMPI<T>*>( dynamic_cast<const VectorPetscMPI<T>*>( &yy ) )->close();
    }
    else
    {
        const_cast<VectorPetsc<T>*>( dynamic_cast<const VectorPetsc<T>*>( &xx ) )->close();
        const_cast<VectorPetsc<T>*>( dynamic_cast<const VectorPetsc<T>*>( &yy ) )->close();
    }

    const VectorPetsc<T>* x = dynamic_cast<const VectorPetsc<T>*>( &xx );
    const VectorPetsc<T>* y = dynamic_cast<const VectorPetsc<T>*>( &yy );
    CHECK( x != 0 && y != 0 ) << "invalid Vector<> types";
    auto ierr =  VecPointwiseDivide(this->vec(), x->vec(), y->vec() );
    CHKERRABORT( this->comm(),ierr );
}
template <typename T>
void
VectorPetsc<T>::zero()
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";

    int ierr=0;

    PetscScalar z=0.;
    this->close();
    // 2.2.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)
    ierr = VecSet ( &z, M_vec );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = VecSet ( M_vec, z );
    CHKERRABORT( this->comm(),ierr );
#endif
}

template <typename T>
int
VectorPetsc<T>::reciprocal()
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    return VecReciprocal( M_vec );
}


template <typename T>
void
VectorPetsc<T>::clear ()
{
    if ( ( this->isInitialized() ) && ( this->M_destroy_vec_on_exit ) )
    {
        int ierr=0;

        ierr = PETSc::VecDestroy( M_vec );
        CHKERRABORT( this->comm(),ierr );
    }

    this->M_is_closed = this->M_is_initialized = false;
}

template <typename T>
void
VectorPetsc<T>::localize(const Vector<T>& /*V*/)
{
    CHECK( 0 ) << "invalid call, not implemented yet";
}

template <typename T>
void
VectorPetsc<T>::insert ( const Vector<T>& /*V*/,
                         const std::vector<size_type>& /*dof_indices*/ )
{
    CHECK( 0 ) << "invalid call, not implemented yet";
}


template <typename T>
void
VectorPetsc<T>::insert ( const ublas::vector<T>& /*V*/,
                         const std::vector<size_type>& /*dof_indices*/ )
{
    CHECK( 0 ) << "invalid call, not implemented yet";
}

template <typename T>
void
VectorPetsc<T>::scale ( T factor_in )
{
    int ierr = 0;
    PetscScalar factor = static_cast<PetscScalar>( factor_in );

    // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = VecScale( &factor, M_vec );
    CHKERRABORT( this->comm(),ierr );

    // 2.3.x & later style
#else

    ierr = VecScale( M_vec, factor );
    CHKERRABORT( this->comm(),ierr );

#endif
}
template <typename T>
void
VectorPetsc<T>::add ( const value_type& v_in )
{
    int ierr=0;
    PetscScalar* values;
    const PetscScalar v = static_cast<PetscScalar>( v_in );
    const int n   = static_cast<int>( this->localSize() );
    const int fli = static_cast<int>( this->firstLocalIndex() );

    for ( int i=0; i<n; i++ )
    {
        ierr = VecGetArray ( M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        int ig = fli + i;

        PetscScalar value = ( values[ig] + v );

        ierr = VecRestoreArray ( M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecSetValues ( M_vec, 1, &ig, &value, INSERT_VALUES );
        CHKERRABORT( this->comm(),ierr );
    }
}
template <typename T>
void
VectorPetsc<T>::add ( const Vector<value_type>& v )
{
    this->add ( 1., v );
}
template <typename T>
void
VectorPetsc<T>::add ( const value_type& a_in, const Vector<value_type>& v_in )
{
    int ierr = 0;
    PetscScalar a = static_cast<PetscScalar>( a_in );

    if ( this->comm().size()>1 )
    {
        const_cast<VectorPetscMPI<T>*>( dynamic_cast<const VectorPetscMPI<T>*>( &v_in ) )->close();
    }
    else
    {
        const_cast<VectorPetsc<T>*>( dynamic_cast<const VectorPetsc<T>*>( &v_in ) )->close();
    }

    const VectorPetsc<T>* v = dynamic_cast<const VectorPetsc<T>*>( &v_in );

    CHECK ( v != NULL ) << "dynamic cast failed";
    CHECK( this->size() == v->size() ) << "invalid vector this.size : " << this->size() << " != v.size: " << v->size();

    // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = VecAXPY( &a, v->M_vec, M_vec );
    CHKERRABORT( this->comm(),ierr );

    // 2.3.x & later style
#else

    ierr = VecAXPY( M_vec, a, v->M_vec );
    CHKERRABORT( this->comm(),ierr );

#endif
}


template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::min () const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";

    int index=0, ierr=0;
    PetscReal min=0.;

    ierr = VecMin ( M_vec, &index, &min );
    CHKERRABORT( this->comm(),ierr );

    // this return value is correct: VecMin returns a PetscReal
    return static_cast<Real>( min );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::max() const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";

    int index=0, ierr=0;
    PetscReal max=0.;

    ierr = VecMax ( M_vec, &index, &max );
    CHKERRABORT( this->comm(),ierr );

    // this return value is correct: VecMax returns a PetscReal
    return static_cast<Real>( max );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>:: l1Norm () const
{
    CHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr=0;
    double value=0.;

    ierr = VecNorm ( M_vec, NORM_1, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::l2Norm () const
{
    CHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr=0;
    double value=0.;

    ierr = VecNorm ( M_vec, NORM_2, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::linftyNorm () const
{
    CHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr=0;
    double value=0.;

    ierr = VecNorm ( M_vec, NORM_INFINITY, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
typename VectorPetsc<T>::value_type
VectorPetsc<T>:: sum () const
{
    CHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr=0;
    double value=0.;

    ierr = VecSum ( M_vec, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
void VectorPetsc<T>::printMatlab ( const std::string name, bool renumber ) const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    LOG_IF( WARNING, ! this->closed() ) <<  "vector is not closed";

    if ( !this->closed() )
    {
        VLOG(1) << "closing vector\n";
        const_cast<VectorPetsc<T>*>( this )->close();
    }

    const_cast<VectorPetsc<T>*>( this )->close();
    PetscObjectSetName((PetscObject)M_vec,fs::path("var_"+name).stem().string().c_str());
    //this->close();
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
        ierr = PetscViewerASCIIOpen( this->comm(),
                                     name.c_str(),
                                     &petsc_viewer );
        CHKERRABORT( this->comm(),ierr );

        ierr = PetscViewerSetFormat ( petsc_viewer,
                                      PETSC_VIEWER_ASCII_MATLAB );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecView ( const_cast<Vec>( M_vec ), petsc_viewer );
        CHKERRABORT( this->comm(),ierr );
    }

    /**
     * Otherwise the matrix will be dumped to the screen.
     */
    else
    {
        ierr = PetscViewerSetFormat ( PETSC_VIEWER_STDOUT_WORLD,
                                      PETSC_VIEWER_ASCII_MATLAB );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecView ( const_cast<Vec>( M_vec ), PETSC_VIEWER_STDOUT_WORLD );
        CHKERRABORT( this->comm(),ierr );
    }


    /**
     * Destroy the viewer.
     */
    ierr = PETSc::PetscViewerDestroy ( petsc_viewer );
    CHKERRABORT( this->comm(),ierr );
}


template <typename T>
typename VectorPetsc<T>::value_type
VectorPetsc<T>::dot( Vector<T> const& __v )
{
    this->close();
    PetscScalar e;

    VectorPetsc<value_type> v( __v.size(), __v.localSize() );
    {
        size_type s = v.localSize();
        size_type start = v.firstLocalIndex();

        for ( size_type i = 0; i < s; ++i )
            v.set( start + i, __v( start + i ) );
    }

    VecDot( this->vec(), v.vec(), &e );

    return e;
}

template <typename T>
void
VectorPetsc<T>::addVector ( const Vector<value_type>& V_in,
                            const MatrixSparse<value_type>& A_in )

    {
        const VectorPetsc<T>* V = dynamic_cast<const VectorPetsc<T>*>( &V_in );
        const MatrixPetsc<T>* A = dynamic_cast<const MatrixPetsc<T>*>( &A_in );

        CHECK ( A != 0 ) << "Invalid PETSc matrix\n";
        this->close();
        A->close();
        int ierr=0;


        if ( !V )
        {
            if ( this->comm().size()>1 )
            {
                VectorPetscMPI<T> tmp( V_in.mapPtr() );
                dynamic_cast<Vector<T>&>( tmp ) = V_in;
                ierr = MatMultAdd( const_cast<MatrixPetsc<T>*>( A )->mat(), tmp.M_vec, M_vec, M_vec );
            }
            else
            {
                VectorPetsc<T> tmp( V_in.mapPtr() );
                dynamic_cast<Vector<T>&>( tmp ) = V_in;
                ierr = MatMultAdd( const_cast<MatrixPetsc<T>*>( A )->mat(), tmp.M_vec, M_vec, M_vec );
            }
        }
        else
        {
            // The const_cast<> is not elegant, but it is required since PETSc
            // is not const-correct.
            ierr = MatMultAdd( const_cast<MatrixPetsc<T>*>( A )->mat(), V->M_vec, M_vec, M_vec );

        }
        CHKERRABORT( this->comm(),ierr );
    }

//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//

template <typename T>
typename VectorPetscMPI<T>::clone_ptrtype
VectorPetscMPI<T>::clone () const
{
    clone_ptrtype cloned_vector ( new VectorPetscMPI<T>( this->mapPtr() ) );
    CHECK( cloned_vector->size() == this->size() ) << "Invalid cloned vector size : " << cloned_vector->size()
                                                   << " expected size : " << this->size() ;
    //*cloned_vector = *this;
    CHECK( this->closed() ) << "VectorPETSc is closed and should not";
    return cloned_vector;
}


template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( datamap_ptrtype const& dm )
    :
    super( dm,false ) //false for not init
{
    this->init( dm->nDof(), dm->nLocalDofWithoutGhost() );
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( Vec v, datamap_ptrtype const& dm )
    :
    super( v,dm )
{
    int ierr=0;
    int petsc_n_localWithGhost=static_cast<int>( this->map().nLocalDofWithGhost() );

    ierr = VecCreateSeq ( PETSC_COMM_SELF, petsc_n_localWithGhost, &  M_vecLocal );
    CHKERRABORT( this->comm(),ierr );

    IS isGlob;
    IS isLoc;

    // create IS for vecScatter
    PetscInt *idx;
    PetscInt n_idx =  this->map().mapGlobalProcessToGlobalCluster().size();
    idx = new PetscInt[n_idx];
    std::copy( this->map().mapGlobalProcessToGlobalCluster().begin(),
               this->map().mapGlobalProcessToGlobalCluster().end(),
               idx );

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral( this->comm(), n_idx, idx, PETSC_COPY_VALUES, &isGlob );
#else
    ierr = ISCreateGeneral( this->comm(), n_idx, idx, &isGlob );
#endif
    CHKERRABORT( this->comm(),ierr );

    ierr = ISCreateStride( PETSC_COMM_SELF,n_idx,0,1,&isLoc );
    CHKERRABORT( this->comm(),ierr );

    // create vecScatter
    ierr = VecScatterCreate( this->vec(), isGlob,
                             M_vecLocal, isLoc,
                             &M_vecScatter );
    CHKERRABORT( this->comm(),ierr );

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy ( &isGlob );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( &isLoc );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy ( isGlob );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( isLoc );
    CHKERRABORT( this->comm(),ierr );
#endif
    delete[] idx;

    this->M_is_initialized = true;

    this->close();
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
void
VectorPetscMPI<T>::init( const size_type n,
                         const size_type n_localWithoutGhost,
                         const bool fast )
{
    //std::cout << "MPI init start" << std::endl;
    int ierr=0;
    int petsc_n=static_cast<int>( n );
    int petsc_n_localWithoutGhost=static_cast<int>( n_localWithoutGhost );
    int petsc_n_localWithGhost=static_cast<int>( this->map().nLocalDofWithGhost() );
    //std::cout << "petsc_n_localWithoutGhost "<< petsc_n_localWithoutGhost << std::endl;
    //std::cout << "petsc_n_localWithGhost "<< petsc_n_localWithGhost << std::endl;

    // Clear initialized vectors
    if ( this->isInitialized() )
        this->clear();

    FEELPP_ASSERT( n_localWithoutGhost < n )( n_localWithoutGhost )( n ).warn( "invalid local size" );

    ierr = VecCreateMPI ( this->comm(), petsc_n_localWithoutGhost, petsc_n,
                          &this->M_vec );
    CHKERRABORT( this->comm(),ierr );

    // localToGlobalMapping
    IS is;
    ISLocalToGlobalMapping isLocToGlobMap;

    PetscInt *idx;
    PetscInt n_idx =  this->map().mapGlobalProcessToGlobalCluster().size();
    idx = new PetscInt[n_idx];
    std::copy( this->map().mapGlobalProcessToGlobalCluster().begin(),
               this->map().mapGlobalProcessToGlobalCluster().end(),
               idx );
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral( this->comm(), n_idx, idx, PETSC_COPY_VALUES, &is );
#else
    ierr = ISCreateGeneral( this->comm(), n_idx, idx, &is );
#endif
    CHKERRABORT( this->comm(),ierr );

    // create LocalToGlobalMapping
    ierr=ISLocalToGlobalMappingCreateIS( is, &isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
    ierr=VecSetLocalToGlobalMapping( this->vec(),isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );

    // create local vector
    ierr = VecCreateSeq ( PETSC_COMM_SELF, petsc_n_localWithGhost, &  M_vecLocal );
    CHKERRABORT( this->comm(),ierr );

    // create vecScatter
    IS isLoc;
    ierr = ISCreateStride( PETSC_COMM_SELF,n_idx,0,1,&isLoc );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecScatterCreate( this->vec(), is,
                             M_vecLocal, isLoc,
                             &M_vecScatter );
    CHKERRABORT( this->comm(),ierr );



    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy( &is );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISLocalToGlobalMappingDestroy( &isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( &isLoc );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy( is );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISLocalToGlobalMappingDestroy( isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( isLoc );
    CHKERRABORT( this->comm(),ierr );
#endif

    delete[] idx;

    ierr = VecSetFromOptions( this->vec() );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecSetFromOptions( M_vecLocal );
    CHKERRABORT( this->comm(),ierr );


    this->M_is_initialized = true;

    if ( fast == false )
        this->zero ();
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
typename VectorPetscMPI<T>::value_type
VectorPetscMPI<T>::operator() ( const size_type i ) const
{
    int ierr=0;
    PetscScalar *values, value=0.;
    ierr = VecGetArray( M_vecLocal, &values );
    CHKERRABORT( this->comm(),ierr );
    //std::cout << "\n operator MPI ";
    value =  values[i /*- this->firstLocalIndex()*/ ];

    ierr = VecRestoreArray( M_vecLocal, &values );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<value_type>( value );
}
template <typename T>
typename VectorPetscMPI<T>::value_type&
VectorPetscMPI<T>::operator() ( const size_type i )
{
    int ierr=0;
    PetscScalar *values;
    ierr = VecGetArray( M_vecLocal, &values );
    CHKERRABORT( this->comm(),ierr );
    //std::cout << "\n operator MPI ";
    PetscScalar& value =  values[i /*- this->firstLocalIndex()*/ ];

    ierr = VecRestoreArray( M_vecLocal, &values );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<value_type&>( value );
}
//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::set( size_type i, const value_type& value )
{
    //FEELPP_ASSERT(i<size())( i )( size() ).error( "invalid index" );
    //if ( this->map().dofGlobalProcessIsGhost(i) ) return;

    int ierr=0;
    int i_val = static_cast<int>( i );
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr=VecSetValuesLocal( this->vec(),1,&i_val,&petsc_value,INSERT_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetscMPI<T>::setVector ( int* i, int n, value_type* v )
{
    DCHECK( this->isInitialized() ) << "vector not initialized";
    DCHECK(n<=this->size()) << "invalid local index array size: " << n << " > " << this->size();

#if PETSC_VERSION_LESS_THAN(3,5,3)
    if ( n == 0 ) return;
#endif
    int ierr=0;
    ierr=VecSetValuesLocal( this->vec(), n, i, v, INSERT_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
void
VectorPetscMPI<T>::add ( const size_type i, const value_type& value )
{
    //FEELPP_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int ierr=0;
    int i_val = static_cast<int>( i );
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr=VecSetValuesLocal( this->vec(), 1, &i_val, &petsc_value, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::addVector ( int* i, int n, value_type* v )
{
    DCHECK( this->isInitialized() ) << "vector not initialized";
    DCHECK(n<=this->size()) << "invalid local index array size: " << n << " > " << this->size();

#if PETSC_VERSION_LESS_THAN(3,5,3)
    if ( n == 0 ) return;
#endif
    int ierr=0;
    ierr=VecSetValuesLocal( this->vec(), n, i, v, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::addVector( const Vector<value_type>& V_in,
                              const MatrixSparse<value_type>& A_in )
{
    super::addVector( V_in,A_in );
    this->localize();
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::clear()
{
    if ( this->isInitialized() )
    {
        super::clear();

       int ierr=0;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = VecDestroy( &M_vecLocal );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecScatterDestroy( &M_vecScatter );
        CHKERRABORT( this->comm(),ierr );
#else
        ierr = VecDestroy( M_vecLocal );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecScatterDestroy( M_vecScatter );
        CHKERRABORT( this->comm(),ierr );
#endif
    }
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void VectorPetscMPI<T>::localize()
{
    int ierr = 0;

    // Perform the scatter
    ierr = VecScatterBegin( M_vecScatter, this->vec(), M_vecLocal, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecScatterEnd  ( M_vecScatter, this->vec(), M_vecLocal, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::close()
{
    //FEELPP_ASSERT (this->isInitialized()).error( "VectorPetsc<> not initialized" );
    //std::cout << "\n MPI CLOSE "<<std::endl;;
    super::close();

    this->localize();

}

//----------------------------------------------------------------------------------------------------//

template <typename T>
size_type
VectorPetscMPI<T>::firstLocalIndex() const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";

    int petsc_first=0;

    return static_cast<size_type>( petsc_first );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
size_type
VectorPetscMPI<T>::lastLocalIndex() const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";

    int petsc_last=this->map().nLocalDofWithGhost();

    return static_cast<size_type>( petsc_last );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::duplicateFromOtherPartition( Vector<T> const& vecInput)
{
    auto testCommActivities_this=this->map().worldComm().hasMultiLocalActivity();

    if (testCommActivities_this.template get<0>())
        {
            //std::cout << "VectorPetscMPI<T>::duplicateFromOtherPartition hasMultiLocalActivity " << std::endl;
            // save initial activities
            std::vector<int> saveActivities_this = this->map().worldComm().activityOnWorld();
            // iterate on each local activity
            const auto colorWhichIsActive = testCommActivities_this.template get<1>();
            auto it_color=colorWhichIsActive.begin();
            auto const en_color=colorWhichIsActive.end();
            for ( ;it_color!=en_color;++it_color )
                {
                    this->map().worldComm().applyActivityOnlyOn( *it_color );
                    this->duplicateFromOtherPartition_run( vecInput );
                }
            // revert initial activities
            this->map().worldComm().setIsActive(saveActivities_this);
        }
    else
        {
            this->duplicateFromOtherPartition_run( vecInput );
        }
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::duplicateFromOtherPartition_run( Vector<T> const& vecInput)
{
    std::list<boost::tuple<size_type,size_type> > memory_dofCluster;
    std::vector<size_type> dofClusterMissing;
    std::vector<size_type> originaldofClusterMissing;
    std::vector<size_type> originaldofClusterMissing_recv;
    std::vector<double> dofClusterMissing_RequestVal;
    std::vector<int> dofClusterMissing_RequestIsFind;

    if (this->map().worldComm().isActive())
    {
        const size_type mynWithGhost = this->map().mapGlobalProcessToGlobalCluster().size();
        for (size_type k=0;k<mynWithGhost;++k)
        {
            const size_type convert_dofProcess = k;
            const size_type convert_dofCluster = this->map().mapGlobalProcessToGlobalCluster()[k];
            const size_type original_dofCluster = convert_dofCluster;
            const size_type original_firstDofCluster = vecInput.map().firstDofGlobalCluster();
            if ( original_dofCluster>=vecInput.map().firstDofGlobalCluster() &&
                 original_dofCluster<=vecInput.map().lastDofGlobalCluster() )
            {
                const size_type original_dofProcess = vecInput.map().mapGlobalClusterToGlobalProcess()[original_dofCluster-original_firstDofCluster];
                this->set( convert_dofProcess, vecInput(original_dofProcess) );
            }
            else
            {
                memory_dofCluster.push_back( boost::make_tuple(k,original_dofCluster) );
            }
        }
        // init data to send
        dofClusterMissing.resize(memory_dofCluster.size());
        originaldofClusterMissing.resize(memory_dofCluster.size());
        auto it_dof = memory_dofCluster.begin();
        for (int k=0;k<dofClusterMissing.size();++k,++it_dof)
        {
            dofClusterMissing[k]=it_dof->get<0>();
            originaldofClusterMissing[k]=it_dof->get<1>();
        }
    }

    auto worldCommFusion = this->map().worldComm()+vecInput.map().worldComm();
    std::vector<rank_type> globalRankToFusionRank_this(this->map().worldComm().globalSize());
    mpi::all_gather( this->map().worldComm().globalComm(),
                     worldCommFusion.globalRank(),
                     globalRankToFusionRank_this );
    std::vector<rank_type> globalRankToFusionRank_input(vecInput.map().worldComm().globalSize());
    mpi::all_gather( vecInput.map().worldComm().globalComm(),
                     worldCommFusion.globalRank(),
                     globalRankToFusionRank_input );

    std::vector<int> thisProcIsActive_fusion(worldCommFusion.globalSize());
    mpi::all_gather( worldCommFusion.globalComm(),
                     (int)this->map().worldComm().isActive(),
                     thisProcIsActive_fusion );
    std::vector<int> inputProcIsActive_fusion(worldCommFusion.globalSize());
    mpi::all_gather( worldCommFusion.globalComm(),
                     (int)vecInput.map().worldComm().isActive(),
                     inputProcIsActive_fusion );

    int firstActiveProc_this=0;
    bool findFirstActive_this=false;
    while (!findFirstActive_this)
    {
        if (thisProcIsActive_fusion[firstActiveProc_this])
        {
            findFirstActive_this=true;
        }
        else ++firstActiveProc_this;
    }
    int firstActiveProc_input=0;
    bool findFirstActive_input=false;
    while (!findFirstActive_input)
    {
        if (inputProcIsActive_fusion[firstActiveProc_input])
        {
            findFirstActive_input=true;
        }
        else ++firstActiveProc_input;
    }


    for (int p=0;p<globalRankToFusionRank_this.size(); ++p)
    {
        if (!this->map().worldComm().isActive()) globalRankToFusionRank_this[p]=p%this->map().worldComm().globalSize()+firstActiveProc_this; // FAIRE COMMMUNICATION!!!!!
    }
    for (int p=0;p<globalRankToFusionRank_input.size(); ++p)
    {
        if (!vecInput.map().worldComm().isActive()) globalRankToFusionRank_input[p]=p%vecInput.map().worldComm().globalSize()+firstActiveProc_input; // FAIRE COMMMUNICATION!!!!!
    }

    std::vector<std::list<int> > searchDistribution(this->map().worldComm().globalSize());
    for (int p=0;p<this->map().worldComm().globalSize();++p)
    {
        searchDistribution[p].clear();
        for (int q=0;q<vecInput.map().worldComm().globalSize();++q)
        {
            //if (q!=p)
            if( (globalRankToFusionRank_this[p])!=globalRankToFusionRank_input[q] )
            {
                searchDistribution[p].push_back(q);
            }
        }
    }

#if 0
    vecInput.map().worldComm().globalComm().barrier();
    for (int p=0;p<vecInput.map().worldComm().globalSize();++p)
    {
        if (p==vecInput.map().worldComm().globalRank())
        {
            std::cout << "I am proc " << p << "\n";
            for (int q=0;q<this->map().worldComm().globalSize();++q)
            {
                auto it_list = searchDistribution[q].begin();
                auto en_list = searchDistribution[q].end();
                for ( ; it_list!=en_list;++it_list) { std::cout << *it_list <<" "; }
                std::cout << std::endl;
            }
        }
        vecInput.map().worldComm().globalComm().barrier();
    }
#endif


    for (int proc=0;proc<this->map().worldComm().globalSize();++proc)
    {
        for (auto it_rankLocalization=searchDistribution[proc].begin(),en_rankLocalization=searchDistribution[proc].end();
             it_rankLocalization!=en_rankLocalization;++it_rankLocalization)
        {
            const int rankLocalization = *it_rankLocalization;
            if ( this->map().worldComm().globalRank() == proc  && thisProcIsActive_fusion[worldCommFusion.globalRank()] )  // send info to rankLocalization
            {
                const int rankToSend = globalRankToFusionRank_input[rankLocalization];
                worldCommFusion.globalComm().send(rankToSend,0,originaldofClusterMissing );
            }
            else if ( vecInput.map().worldComm().globalRank()==rankLocalization && inputProcIsActive_fusion[worldCommFusion.globalRank()] )
            {
                const int rankToRecv = globalRankToFusionRank_this[proc];
                worldCommFusion.globalComm().recv(rankToRecv,0,originaldofClusterMissing_recv );

                const size_type nDataRecv = originaldofClusterMissing_recv.size();
                dofClusterMissing_RequestVal.resize(nDataRecv);
                dofClusterMissing_RequestIsFind.resize(nDataRecv);
                for (size_type k=0;k<nDataRecv;++k)
                {
                    const size_type original_firstDofCluster = vecInput.map().firstDofGlobalCluster();
                    const size_type original_dofCluster = originaldofClusterMissing_recv[k];
                    if (original_dofCluster >=vecInput.map().firstDofGlobalCluster() && original_dofCluster<=vecInput.map().lastDofGlobalCluster())
                    {
                        const size_type original_dofProcess = vecInput.map().mapGlobalClusterToGlobalProcess()[original_dofCluster-original_firstDofCluster];
                        dofClusterMissing_RequestVal[k]=vecInput(original_dofProcess);
                        dofClusterMissing_RequestIsFind[k]=1;
                    }
                    else dofClusterMissing_RequestIsFind[k]=0;
                }
                worldCommFusion.globalComm().send( rankToRecv, 1, dofClusterMissing_RequestVal );
                worldCommFusion.globalComm().send( rankToRecv, 2, dofClusterMissing_RequestIsFind );
            }

            if ( this->map().worldComm().globalRank() == proc && thisProcIsActive_fusion[worldCommFusion.globalRank()]  )
            {
                const int rankToRecv = globalRankToFusionRank_input[rankLocalization];
                worldCommFusion.globalComm().recv( rankToRecv, 1, dofClusterMissing_RequestVal );
                worldCommFusion.globalComm().recv( rankToRecv, 2, dofClusterMissing_RequestIsFind );

                const size_type nDataRecv = dofClusterMissing_RequestVal.size();
                for (size_type k=0;k<nDataRecv;++k)
                {
                    if (dofClusterMissing_RequestIsFind[k])
                    {
                        const size_type convert_dofProcess = dofClusterMissing[k];
                        this->set( convert_dofProcess,dofClusterMissing_RequestVal[k]);
                    }
                }
            }
            //---------------------------------------
            worldCommFusion.globalComm().barrier();
            //---------------------------------------
        }
    }
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
typename VectorPetsc<T>::value_type
VectorPetscMPI<T>::dot( Vector<T> const& __v )
{
    this->close();
    PetscScalar e;

    VectorPetscMPI<value_type> v( this->mapPtr() );
    {
        size_type s = v.map().nLocalDofWithGhost();
        size_type start = v.firstLocalIndex();

        for ( size_type i = 0; i < s; ++i )
            v.set( start + i, __v( start + i ) );
    }

    v.close();

    VecDot( this->vec(), v.vec(), &e );

    return e;
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
size_type
VectorPetscMPI<T>::localSize() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "VectorPetsc not initialized" );

    int petsc_size=0;
    int ierr = VecGetLocalSize( M_vecLocal, &petsc_size );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<size_type>( petsc_size );
}


template class VectorPetsc<double>;
template class VectorPetscMPI<double>;

} // Feel

#endif // FEELPP_HAS_PETSC_H
