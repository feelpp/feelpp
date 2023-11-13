/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#include <feel/feelalg/topetsc.hpp>
#include <feel/feeltiming/tic.hpp>
#if BOOST_VERSION < 105900
#include <boost/smart_ptr/make_shared.hpp>
#endif

BOOST_CLASS_EXPORT_IMPLEMENT( Feel::VectorPetsc<double> )
BOOST_CLASS_EXPORT_IMPLEMENT( Feel::VectorPetscMPI<double> )

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
    clone_ptrtype cloned_vector;
    if ( dynamic_cast<VectorPetscMPIRange<T> const*>( this ) )
        cloned_vector.reset( new VectorPetscMPIRange<T>( this->mapPtr() ) );
    else if ( dynamic_cast<VectorPetscMPI<T> const*>( this ) )
        cloned_vector.reset( new VectorPetscMPI<T>( this->mapPtr() ) );
    else
        cloned_vector.reset( new VectorPetsc<T>( this->mapPtr() ) );
    CHECK( cloned_vector->size() == this->size() ) << "Invalid cloned vector size : " << cloned_vector->size()
                                                   << " expected size : " << this->size() ;
    //*cloned_vector = *this;
    //CHECK( this->closed() ) << "VectorPETSc is closed and should not";
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
    int petsc_n = static_cast<PetscInt>( n );
    int petsc_n_local=static_cast<PetscInt>( n_local );

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

    this->setInitialized( true );
    this->setIsClosed( true );

    if ( !fast )
        this->zero ();
}

template<typename T>
void
VectorPetsc<T>::init( datamap_ptrtype const& dm )
{
    if ( !this->map().isCompatible( *dm ) )
        this->setMap( dm );
    this->init( this->map().nLocalDofWithoutGhost(), this->map().nDof() );
}

template <typename T>
void
VectorPetsc<T>::set( const value_type& value )
{
    int ierr=0;
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr = VecSet ( M_vec, petsc_value );
    CHKERRABORT( this->comm(),ierr );
}
template <typename T>
void
VectorPetsc<T>::set( const size_type i, const value_type& value )
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
VectorPetsc<T>::addVector ( int* i, int n, value_type* v, size_type K, size_type K2 )
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
const typename VectorPetsc<T>::value_type &
VectorPetsc<T>::operator() ( const size_type i ) const
{
    DCHECK( this->isInitialized() ) << "VectorPETSc not initialized";
    DCHECK ( ( ( i >= this->firstLocalIndex() ) &&
               ( i <  this->lastLocalIndex() ) ) ) << "invalid vector index " <<  i
                                                   << " first local index: "  << this->firstLocalIndex()
                                                   << " last local index:  " << this->lastLocalIndex();

    int ierr=0;
    const PetscScalar *values;


    ierr = VecGetArrayRead( M_vec, &values );
    CHKERRABORT( this->comm(),ierr );

    const PetscScalar& value = values[i - this->firstLocalIndex()];

    ierr = VecRestoreArrayRead ( M_vec, &values );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<const value_type&>( value );
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
VectorPetsc<T>::close()
{
    if ( !this->isInitialized() )
        return;

    int ierr=0;

    ierr = VecAssemblyBegin( M_vec );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecAssemblyEnd( M_vec );
    CHKERRABORT( this->comm(),ierr );

    this->M_is_closed = true;
}

template <typename T>
Vector<typename VectorPetsc<T>::value_type>&
VectorPetsc<T>::operator= ( const Vector<value_type> &V )
{
    if ( !this->closed() )
        this->close();
    if ( !V.closed() )
        const_cast<Vector<T>*>( &V )->close();

    const VectorPetsc<T>* vecPetsc =  dynamic_cast<const VectorPetsc<T>*>( &V );
    if ( vecPetsc )
    {
        if ( !this->map().isCompatible( V.map() ) )
            this->setMap( V.mapPtr() );

        int ierr=0;
        ierr = VecCopy( vecPetsc->vec(), M_vec );
        CHKERRABORT( this->comm(),ierr );
        return *this;
    }

    auto vecPetscCastPair = Feel::detail::toPETScPairPtr( V, false );
    VectorPetsc<T> const* vecPetscCast = vecPetscCastPair.first;
    if ( vecPetscCast )
        return this->operator=( *vecPetscCast );

    return super::operator=( V );
}

template <typename T>
Vector<typename VectorPetsc<T>::value_type>&
VectorPetsc<T>::operator= ( const VectorPetsc<value_type> &V )
{
    return this->operator=( *dynamic_cast< Vector<value_type> const* >( &V ) );
}

template <typename T>
void
VectorPetsc<T>::pointwiseOperationsImpl( Vector<T> const& xx, Vector<T> const& yy, int op )
{
   if ( !this->closed() )
       this->close();
   if ( !xx.closed() )
        const_cast<Vector<T>*>( &xx )->close();
   if ( !yy.closed() )
        const_cast<Vector<T>*>( &yy )->close();

    const VectorPetsc<T>* vecxPetsc = dynamic_cast<const VectorPetsc<T>*>( &xx );
    const VectorPetsc<T>* vecyPetsc = dynamic_cast<const VectorPetsc<T>*>( &yy );
    if ( vecxPetsc && vecyPetsc )
    {
        int ierr = 0;
        if ( op == 0 )
            ierr =  VecPointwiseMult( this->vec(), vecxPetsc->vec(), vecyPetsc->vec() );
        else
            ierr =  VecPointwiseDivide( this->vec(), vecxPetsc->vec(), vecyPetsc->vec() );
        CHKERRABORT( this->comm(),ierr );
        return;
    }

    // try to use a view petsc format else create a new vector (values copied)
    auto vecPetscCastPair_xx = Feel::detail::toPETScPairPtr( xx, true );
    VectorPetsc<T> const* vecPetscCast_xx = vecPetscCastPair_xx.first;
    auto vecPetscCastPair_yy = Feel::detail::toPETScPairPtr( yy, true );
    VectorPetsc<T> const* vecPetscCast_yy = vecPetscCastPair_yy.first;
    if ( vecPetscCast_xx && vecPetscCast_yy )
    {
        if ( op == 0 )
            this->pointwiseMult( *vecPetscCast_xx,*vecPetscCast_yy );
        else
            this->pointwiseDivide( *vecPetscCast_xx,*vecPetscCast_yy );
        return;
    }

    CHECK( false ) << "TODO other kind of vector";

}

template <typename T>
void
VectorPetsc<T>::pointwiseMult ( Vector<T> const& xx, Vector<T> const& yy )
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    this->pointwiseOperationsImpl( xx,yy,0 );
}

template <typename T>
void
VectorPetsc<T>::pointwiseDivide ( Vector<T> const& xx, Vector<T> const& yy )
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    this->pointwiseOperationsImpl( xx,yy,1 );
}
template <typename T>
void
VectorPetsc<T>::zero()
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        this->close();

    int ierr=0;
    PetscScalar z=0.;
    //this->close();
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
    if ( !this->closed() )
        this->close();

    int ierr=0;
    ierr = VecReciprocal( M_vec );
    CHKERRABORT( this->comm(),ierr );
    return ierr;
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
    if ( !this->closed() )
        this->close();

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
    if ( !this->closed() )
        this->close();

    const PetscScalar v = static_cast<PetscScalar>( v_in );
    int ierr=0;
    ierr = VecShift( M_vec, v );
    CHKERRABORT( this->comm(),ierr );
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
    if ( !this->closed() )
        this->close();
    if ( !v_in.closed() )
        const_cast<Vector<T>*>( &v_in )->close();

    int ierr = 0;
    PetscScalar a = static_cast<PetscScalar>( a_in );

    const VectorPetsc<T>* vecPetsc =  dynamic_cast<const VectorPetsc<T>*>( &v_in );
    if ( vecPetsc )
    {
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
        ierr = VecAXPY( &a, vecPetsc->M_vec, M_vec );
        CHKERRABORT( this->comm(),ierr );
#else
        ierr = VecAXPY( M_vec, a, vecPetsc->M_vec );
        CHKERRABORT( this->comm(),ierr );
#endif
        return;
    }

    // try to use a view petsc format else create a new vector (values copied)
    auto vecPetscCastPairIn = Feel::detail::toPETScPairPtr( v_in, true );
    VectorPetsc<T> const* vecPetscCastIn = vecPetscCastPairIn.first;
    if ( vecPetscCastIn )
    {
        this->add( a_in, *vecPetscCastIn );
        return;
    }
}


template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::min () const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

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
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int index, ierr=0;
    PetscReal max=0.;

    ierr = VecMax ( M_vec, &index, &max );
    CHKERRABORT( this->comm(),ierr );

    // this return value is correct: VecMax returns a PetscReal
    return static_cast<Real>( max );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::maxWithIndex( int* index ) const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr=0;
    PetscReal max=0.;

    ierr = VecMax ( M_vec, index, &max );
    CHKERRABORT( this->comm(),ierr );

    // this return value is correct: VecMax returns a PetscReal
    return static_cast<Real>( max );
}


template <typename T>
void
VectorPetsc<T>::abs()
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr = 0;
    ierr = VecAbs( M_vec );
    CHKERRABORT( this->comm(),ierr );
}


template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::l1Norm () const
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
VectorPetsc<T>::sum () const
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

#if PETSC_VERSION_LESS_THAN( 3, 7, 0 )
        ierr = PetscViewerSetFormat( petsc_viewer,
                                     PETSC_VIEWER_ASCII_MATLAB );
#else
        ierr = PetscViewerPushFormat( petsc_viewer,
                                      PETSC_VIEWER_ASCII_MATLAB );
#endif
        CHKERRABORT( this->comm(),ierr );

        ierr = VecView ( const_cast<Vec>( M_vec ), petsc_viewer );
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
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
        ierr = PetscViewerPushFormat ( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB );
#else
        ierr = PetscViewerSetFormat ( PETSC_VIEWER_STDOUT_WORLD,
                                      PETSC_VIEWER_ASCII_MATLAB );
#endif
        CHKERRABORT( this->comm(),ierr );

        ierr = VecView ( const_cast<Vec>( M_vec ), PETSC_VIEWER_STDOUT_WORLD );
        CHKERRABORT( this->comm(),ierr );
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
        ierr = PetscViewerPopFormat ( PETSC_VIEWER_STDOUT_WORLD);
        CHKERRABORT( this->comm(),ierr );
#endif
    }


    /**
     * Destroy the viewer.
     */
    ierr = PETSc::PetscViewerDestroy ( petsc_viewer );
    CHKERRABORT( this->comm(),ierr );
}


template <typename T>
typename VectorPetsc<T>::value_type
VectorPetsc<T>::dot( Vector<T> const& __v ) const
{
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();


    const VectorPetsc<T>* vecPetsc =  dynamic_cast<const VectorPetsc<T>*>( &__v );
    if ( vecPetsc )
    {
        int ierr = 0;
        PetscScalar e;
        ierr = VecDot( this->vec(), vecPetsc->vec(), &e );
        CHKERRABORT( this->comm(),ierr );
        return e;
    }

    // try to use a view petsc format else create a new vector (values copied)
    auto vecPetscCastPairIn = Feel::detail::toPETScPairPtr( __v, true );
    VectorPetsc<T> const* vecPetscCastIn = vecPetscCastPairIn.first;
    if ( vecPetscCastIn )
        return this->dot( *vecPetscCastIn );

    return value_type(0);
}

template <typename T>
void
VectorPetsc<T>::addVector ( const Vector<value_type>& V_in,
                            const MatrixSparse<value_type>& A_in )
{
    if ( !this->closed() )
        this->close();
    if ( !V_in.closed() )
        const_cast<Vector<T>*>( &V_in )->close();
    if ( !A_in.closed() )
        A_in.close();

    const MatrixPetsc<T>* A = dynamic_cast<const MatrixPetsc<T>*>( &A_in );
    CHECK( A != 0 ) << "Invalid PETSc matrix\n";

    int ierr = 0;
    const VectorPetsc<T>* vecPetsc =  dynamic_cast<const VectorPetsc<T>*>( &V_in );
    if ( vecPetsc )
    {
        ierr = MatMultAdd( const_cast<MatrixPetsc<T>*>( A )->mat(), vecPetsc->M_vec, M_vec, M_vec );
        CHKERRABORT( this->comm(),ierr );
        // update ghost values (do nothing in sequential)
        this->localize();
        return;
    }

    // try to use a view petsc format else create a new vector (values copied)
    auto vecPetscCastPairIn = Feel::detail::toPETScPairPtr( V_in, true );
    VectorPetsc<T> const* vecPetscCastIn = vecPetscCastPairIn.first;
    if ( vecPetscCastIn )
    {
        this->addVector( *vecPetscCastIn, A_in );
        return;
    }

}

template <typename T>
std::shared_ptr<Vector<T> >
VectorPetsc<T>::createSubVector( std::vector<size_type> const& _rows,
                                 bool checkAndFixRange ) const
{
    // update maybe input index set
    std::vector<size_type> rows = ( checkAndFixRange )?
        this->mapPtr()->buildIndexSetWithParallelMissingDof( _rows ) : _rows;

    // build subdatamap row
    datamap_ptrtype subMapRow = this->mapPtr()->createSubDataMap( rows, false );

    // build subvector petsc
    Vec subVecPetsc = NULL;
    this->getSubVectorPetsc( rows, subVecPetsc );

    // build vectorsparse object
    std::shared_ptr<Vector<T> > subVec;
    if ( this->comm().size()>1 )
        subVec.reset( new VectorPetscMPI<T>( subVecPetsc,subMapRow,true ) );
    else
        subVec.reset( new VectorPetsc<T>( subVecPetsc,subMapRow,true ) );

    return subVec;
}

template <typename T>
void
VectorPetsc<T>::updateSubVector( std::shared_ptr<Vector<T> > & subvector,
                                 std::vector<size_type> const& rows,
                                 bool init )
{
    CHECK( subvector ) << "subvector is not init";
    std::shared_ptr<VectorPetsc<T> > subvectorPetsc = std::dynamic_pointer_cast<VectorPetsc<T> >( subvector );
    this->getSubVectorPetsc( rows, subvectorPetsc->vec(), init );
}

template <typename T>
void
VectorPetsc<T>::getSubVectorPetsc( std::vector<size_type> const& rows,
                                   Vec &subvec,
                                   bool init ) const
{
    if ( !this->closed() )
        const_cast<VectorPetsc<T>*>( this )->close();

    int ierr=0;
    IS isrow;
    PetscInt *rowMap;

    std::set<size_type> rowMapOrdering;
    if ( this->comm().size()>1 )
    {
        // convert global process ids into global cluster ids, remove ghost dofs
        // and build ordering row map
        for (int i=0; i<rows.size(); i++)
        {
            if ( !this->map().dofGlobalProcessIsGhost( rows[i] ) )
                rowMapOrdering.insert( this->map().mapGlobalProcessToGlobalCluster( rows[i] ) );
        }
    }
    else
    {
        // build ordering row map
        rowMapOrdering.insert( rows.begin(), rows.end() );
    }

    // copying into PetscInt vector
    int nrow = rowMapOrdering.size();
    rowMap = new PetscInt[nrow];
    size_type curId=0;
    for ( auto& rowId : rowMapOrdering )
    {
        rowMap[curId] = rowId;
        ++curId;
    }


#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(),nrow,rowMap,PETSC_COPY_VALUES,&isrow);
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISCreateGeneral(this->comm(),nrow,rowMap,&isrow);
    CHKERRABORT( this->comm(),ierr );
#endif

    if( subvec == NULL ) //createSubVector
    {
        ierr = VecGetSubVector(this->vec(), isrow, &subvec);
        CHKERRABORT( this->comm(),ierr );
    }
    else //updateSubVector
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 5)
        //ierr = VecRestoreSubVector(this->vec(), isrow, &subvec);
        if( init )
        {
            ierr = VecISSet(this->vec(), isrow, 0); //re-init isrow indices to zero
            CHKERRABORT( this->comm(),ierr );
        }
        ierr = VecISAXPY(this->vec(), isrow, 1, subvec); //vec[isrow[i]] += alpha*subvec[i] with alpha=1
        CHKERRABORT( this->comm(),ierr );
#else
        std::cerr << "ERROR : update of subvectors requires petsc version >= 3.5" << std::endl;
#endif
    }

    ierr = PETSc::ISDestroy( isrow );
    CHKERRABORT( this->comm(),ierr );

    delete[] rowMap;
}
template <typename T>
void
VectorPetsc<T>::save( std::string const& filename, std::string const& format )
{
    std::string fullfilename = boost::str( boost::format("%1%_%2%_%3%") %filename %format %this->comm().rank() );
    std::ofstream ofs( fullfilename );
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
            CHECK( false ) << "VectorPetsc save() function : error with unknown format " << format;
    }
    else
    {
        CHECK( false ) << "VectorPetsc save() function : error opening ofstream with name " << filename;
    }
}
template <typename T>
void
VectorPetsc<T>::load( std::string const& filename, std::string const& format )
{
    std::string fullfilename = boost::str( boost::format("%1%_%2%_%3%") %filename %format %this->comm().rank() );
    std::ifstream ifs( fullfilename );
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
            CHECK( false ) << "VectorPetsc save() function : error with unknown format " << format;
    }
    else
    {
        CHECK( false ) << "VectorPetsc load() function : error opening ofstream with name "<< filename;
    }
}


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//


template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( datamap_ptrtype const& dm, bool doInit )
    :
    super( dm,false ) //false for not init
{
    if ( doInit )
        this->initImpl();

    //this->init( dm->nDof(), dm->nLocalDofWithoutGhost() );
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( Vec v, datamap_ptrtype const& dm, bool duplicate )
    :
    super( dm,false )
{
    this->M_destroy_vec_on_exit = duplicate;

    if ( duplicate )
    {
        int ierr = 0;
        ierr = VecDuplicate( v, &this->M_vec );
        CHKERRABORT( this->comm(),ierr );

        Vec lx;
        ierr = VecGhostGetLocalForm(this->vec(),&lx);
        CHKERRABORT( this->comm(),ierr );
        if ( !lx )
        {
            int petsc_n_localWithoutGhost = static_cast<PetscInt>( this->map().nLocalDofWithoutGhost() );
            int petsc_n_localGhost = static_cast<PetscInt>( this->map().nLocalGhosts() );
            PetscInt *idx = NULL;
            if ( petsc_n_localGhost > 0 )
            {
                idx = new PetscInt[petsc_n_localGhost];
                std::copy( this->map().mapGlobalProcessToGlobalCluster().begin()+petsc_n_localWithoutGhost,
                           this->map().mapGlobalProcessToGlobalCluster().end(),
                           idx );
            }

            ierr = VecMPISetGhost( this->vec(), petsc_n_localGhost, idx );
            CHKERRABORT( this->comm(),ierr );
        }
        ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
        CHKERRABORT( this->comm(),ierr );

        ierr = VecCopy( v, this->vec() );
        CHKERRABORT( this->comm(),ierr );

        // update ghosts
        this->M_is_initialized = true;
        this->localize();
    }
    else
    {
        this->M_vec = v;
        this->M_is_initialized = true;

        // make sure that ghosts are updated
        this->localize();
    }

    this->setIsClosed( true );
}

//----------------------------------------------------------------------------------------------------//
template<typename T>
void
VectorPetscMPI<T>::init( const size_type n,
                         const size_type n_localWithoutGhost,
                         const bool fast )
{
    CHECK( false ) << "not allowed";
}

template<typename T>
void
VectorPetscMPI<T>::init( datamap_ptrtype const& dm )
{
    if ( !this->map().isCompatible( *dm ) )
        this->setMap( dm );
    this->initImpl();
}

template<typename T>
void
VectorPetscMPI<T>::initImpl( const bool fast )
{
    int ierr=0;
    int petsc_n = static_cast<PetscInt>( this->map().nDof() );
    int petsc_n_localWithoutGhost = static_cast<PetscInt>( this->map().nLocalDofWithoutGhost() );
    int petsc_n_localWithGhost = static_cast<PetscInt>( this->map().nLocalDofWithGhost() );
    int petsc_n_localGhost = static_cast<PetscInt>( this->map().nLocalGhosts() );
    FEELPP_ASSERT( petsc_n_localWithoutGhost < petsc_n )( petsc_n_localWithoutGhost )( petsc_n ).warn( "invalid local size" );

    // clear initialized vectors
    if ( this->isInitialized() )
        this->clear();

    PetscInt *idx = NULL;
    if ( petsc_n_localGhost > 0 )
    {
        idx = new PetscInt[petsc_n_localGhost];
        std::copy( this->map().mapGlobalProcessToGlobalCluster().begin()+petsc_n_localWithoutGhost,
                   this->map().mapGlobalProcessToGlobalCluster().end(),
                   idx );
    }
    ierr = VecCreateGhostWithArray( this->comm(),
                                    petsc_n_localWithoutGhost, petsc_n,
                                    petsc_n_localGhost,
                                    idx,
                                    NULL, &this->M_vec );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecSetFromOptions( this->vec() );
    CHKERRABORT( this->comm(),ierr );

    if ( petsc_n_localGhost > 0 )
        delete[] idx;

    this->setInitialized( true );
    this->setIsClosed( true );

    if ( !fast )
        this->zero ();
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
const typename VectorPetscMPI<T>::value_type&
VectorPetscMPI<T>::operator() ( const size_type i ) const
{
    int ierr=0;
    Vec lx;
    ierr = VecGhostGetLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
#if 0
    CHECK( lx ) << "is not a GhostGetLocalForm";
    PetscBool hola;
    ierr = VecGhostIsLocalForm(this->vec(),lx, &hola);
    CHKERRABORT( this->comm(),ierr );
    CHECK( hola ) << "is not a GhostGetLocalForm2";
#endif
    const PetscScalar *values;
    ierr = VecGetArrayRead(lx,&values);
    CHKERRABORT( this->comm(),ierr );
    const PetscScalar& value =  values[i];
    ierr = VecRestoreArrayRead( lx, &values );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    return static_cast<const value_type&>( value );
}
template <typename T>
typename VectorPetscMPI<T>::value_type&
VectorPetscMPI<T>::operator() ( const size_type i )
{
    int ierr=0;
    Vec lx;
    ierr = VecGhostGetLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
#if 0
    CHECK( lx ) << "is not a GhostGetLocalForm";
    PetscBool hola;
    ierr = VecGhostIsLocalForm(this->vec(),lx, &hola);
    CHKERRABORT( this->comm(),ierr );
    CHECK( hola ) << "is not a GhostGetLocalForm2";
#endif
    PetscScalar *values;
    ierr = VecGetArray(lx,&values);
    CHKERRABORT( this->comm(),ierr );
    PetscScalar& value =  values[i];
    ierr = VecRestoreArray( lx, &values );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    return static_cast<value_type&>( value );
}


template <typename T>
Vector<typename VectorPetscMPI<T>::value_type>&
VectorPetscMPI<T>::operator= ( const Vector<value_type> &V )
{
    if ( !this->closed() )
        this->close();
    if ( !V.closed() )
        const_cast<Vector<T>*>( &V )->close();

    int ierr=0;
    VectorPetscMPI<T>* vecOutPetscMPIRange =  dynamic_cast<VectorPetscMPIRange<T>*>( &(*this) );
    if ( !vecOutPetscMPIRange )
    {
        if ( !this->map().isCompatible( V.map() ) )
            this->setMap( V.mapPtr() );

        const VectorPetscMPIRange<T>* vecPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &V );
        if ( vecPetscMPIRange )
        {
            ierr = VecCopy( vecPetscMPIRange->vec(), this->vec() );
            CHKERRABORT( this->comm(),ierr );

            size_type nLocalDofWithoutGhost = vecPetscMPIRange->map().nLocalDofWithoutGhost();
            size_type nLocalGhost = vecPetscMPIRange->map().nLocalGhosts();

            Vec lxOut;
            ierr = VecGhostGetLocalForm(this->vec(),&lxOut);
            CHKERRABORT( this->comm(),ierr );

            PetscScalar *valuesOut;
            PetscScalar *valuesInGhost;
            ierr = VecGetArray( lxOut, &valuesOut );
            CHKERRABORT( this->comm(),ierr );
            ierr = VecGetArray( vecPetscMPIRange->vecGhost(), &valuesInGhost );
            CHKERRABORT( this->comm(),ierr );
            for ( size_type k=0;k<nLocalGhost;++k )
                valuesOut[nLocalDofWithoutGhost+k] = valuesInGhost[k];
            ierr = VecRestoreArray( lxOut, &valuesOut );
            CHKERRABORT( this->comm(),ierr );
            ierr = VecRestoreArray( vecPetscMPIRange->vecGhost(), &valuesInGhost );
            CHKERRABORT( this->comm(),ierr );

            ierr = VecGhostRestoreLocalForm(this->vec(),&lxOut);
            CHKERRABORT( this->comm(),ierr );

            return *this;
        }
        const VectorPetscMPI<T>* vecPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &V );
        if ( vecPetscMPI )
        {
            Vec lx;
            ierr = VecGhostGetLocalForm(this->vec(),&lx);
            CHKERRABORT( this->comm(),ierr );
            Vec lxIn;
            ierr = VecGhostGetLocalForm(vecPetscMPI->vec(),&lxIn);
            CHKERRABORT( this->comm(),ierr );
            ierr = VecCopy( lxIn, lx );
            CHKERRABORT( this->comm(),ierr );
            ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
            CHKERRABORT( this->comm(),ierr );
            ierr = VecGhostRestoreLocalForm(vecPetscMPI->vec(),&lxIn);
            CHKERRABORT( this->comm(),ierr );
            return *this;
        }
    }
    return super::operator=( V );
}

template <typename T>
Vector<typename VectorPetsc<T>::value_type>&
VectorPetscMPI<T>::operator= ( const VectorPetscMPI<value_type> &V )
{
    return this->operator=( *dynamic_cast< Vector<value_type> const* >( &V ) );
}

template <typename T>
void
VectorPetscMPI<T>::set( const value_type& value )
{
    if ( !this->closed() )
        this->close();
    int ierr=0;
    const PetscScalar val = static_cast<PetscScalar>( value );
    Vec lx;
    ierr = VecGhostGetLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    ierr = VecSet ( lx, val );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetscMPI<T>::add( const value_type& v_in )
{
    if ( !this->closed() )
        this->close();
    int ierr=0;
    const PetscScalar v = static_cast<PetscScalar>( v_in );
    Vec lx;
    ierr = VecGhostGetLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    ierr = VecShift( lx, v );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetscMPI<T>::add( const value_type& a_in, const Vector<value_type>& v_in )
{
    if ( !this->closed() )
        this->close();
    if ( !v_in.closed() )
        const_cast<Vector<T>*>( &v_in )->close();

    int ierr = 0;
    PetscScalar a = static_cast<PetscScalar>( a_in );

    VectorPetscMPI<T>* vecOutPetscMPIRange =  dynamic_cast<VectorPetscMPIRange<T>*>( &(*this) );
    if ( !vecOutPetscMPIRange )
    {
        const VectorPetscMPIRange<T>* vecPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &v_in );
        if ( vecPetscMPIRange )
        {
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
            ierr = VecAXPY( &a, vecPetscMPIRange->vec(), this->vec() );
            CHKERRABORT( this->comm(),ierr );
#else
            ierr = VecAXPY( this->vec(), a, vecPetscMPIRange->vec() );
            CHKERRABORT( this->comm(),ierr );
#endif
            size_type nLocalDofWithoutGhost = vecPetscMPIRange->map().nLocalDofWithoutGhost();
            size_type nLocalGhost = vecPetscMPIRange->map().nLocalGhosts();

            Vec lxOut;
            ierr = VecGhostGetLocalForm(this->vec(),&lxOut);
            CHKERRABORT( this->comm(),ierr );

            PetscScalar *valuesOut;
            PetscScalar *valuesInGhost;
            ierr = VecGetArray( lxOut, &valuesOut );
            CHKERRABORT( this->comm(),ierr );
            ierr = VecGetArray( vecPetscMPIRange->vecGhost(), &valuesInGhost );
            CHKERRABORT( this->comm(),ierr );
            for ( size_type k=0;k<nLocalGhost;++k )
                valuesOut[nLocalDofWithoutGhost+k] += a*valuesInGhost[k];
            ierr = VecRestoreArray( lxOut, &valuesOut );
            CHKERRABORT( this->comm(),ierr );
            ierr = VecRestoreArray( vecPetscMPIRange->vecGhost(), &valuesInGhost );
            CHKERRABORT( this->comm(),ierr );

            ierr = VecGhostRestoreLocalForm(this->vec(),&lxOut);
            CHKERRABORT( this->comm(),ierr );

            return;
        }
        const VectorPetscMPI<T>* vecPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &v_in );
        if ( vecPetscMPI )
        {
            //CHECK( this->size() == v->size() ) << "invalid vector this.size : " << this->size() << " != v.size: " << v->size();
            Vec lx;
            ierr = VecGhostGetLocalForm(this->vec(),&lx);
            CHKERRABORT( this->comm(),ierr );
            Vec lxIn;
            ierr = VecGhostGetLocalForm(vecPetscMPI->vec(),&lxIn);
            CHKERRABORT( this->comm(),ierr );
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
            ierr = VecAXPY( &a, lxIn, lx );
            CHKERRABORT( this->comm(),ierr );
#else
            ierr = VecAXPY( lx, a, lxIn );
            CHKERRABORT( this->comm(),ierr );
#endif
            ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
            CHKERRABORT( this->comm(),ierr );
            ierr = VecGhostRestoreLocalForm(vecPetscMPI->vec(),&lxIn);
            CHKERRABORT( this->comm(),ierr );
            return;
        }
    }

    // try others casts
    super::add( a_in, v_in );
}


template <typename T>
void
VectorPetscMPI<T>::set( size_type i, const value_type& value )
{
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
    int ierr=0;
    int i_val = static_cast<int>( i );
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr=VecSetValuesLocal( this->vec(), 1, &i_val, &petsc_value, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::addVector ( int* i, int n, value_type* v, size_type K, size_type K2 )
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
VectorPetscMPI<T>::pointwiseMult( Vector<T> const& x, Vector<T> const& y )
{
    if ( !this->closed() )
        this->close();
    if ( !x.closed() )
        const_cast<Vector<T>*>( &x )->close();
    if ( !y.closed() )
        const_cast<Vector<T>*>( &y )->close();

    int ierr = 0;
    const VectorPetscMPI<T>* vecxPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &x );
    const VectorPetscMPI<T>* vecyPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &y );

    if( vecxPetscMPI && vecyPetscMPI )
    {
        VectorPetscMPI<T>* vecOutPetscMPIRange =  dynamic_cast<VectorPetscMPIRange<T>*>( &(*this) );
        if ( !vecOutPetscMPIRange )
        {
            const VectorPetscMPIRange<T>* vecxPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &x );
            const VectorPetscMPIRange<T>* vecyPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &y );
            if ( vecxPetscMPI && vecyPetscMPI && !vecxPetscMPIRange && !vecyPetscMPIRange )
            {
                Vec lx;
                ierr = VecGhostGetLocalForm(this->vec(),&lx);
                CHKERRABORT( this->comm(),ierr );
                Vec lxInx;
                ierr = VecGhostGetLocalForm(vecxPetscMPI->vec(),&lxInx);
                CHKERRABORT( this->comm(),ierr );
                Vec lxIny;
                ierr = VecGhostGetLocalForm(vecyPetscMPI->vec(),&lxIny);
                CHKERRABORT( this->comm(),ierr );

                ierr =  VecPointwiseMult(lx, lxInx, lxIny );
                CHKERRABORT( this->comm(),ierr );

                ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
                CHKERRABORT( this->comm(),ierr );
                ierr = VecGhostRestoreLocalForm(vecxPetscMPI->vec(),&lxInx);
                CHKERRABORT( this->comm(),ierr );
                ierr = VecGhostRestoreLocalForm(vecyPetscMPI->vec(),&lxIny);
                CHKERRABORT( this->comm(),ierr );
                return;
            }
        }

        this->pointwiseOperationOthersPetscImpl(x,y,0);
        return;
    }
    // default method
    super::pointwiseMult( x,y );
}

template <typename T>
void
VectorPetscMPI<T>::pointwiseDivide( Vector<T> const& x, Vector<T> const& y )
{
    if ( !this->closed() )
        this->close();
    if ( !x.closed() )
        const_cast<Vector<T>*>( &x )->close();
    if ( !y.closed() )
        const_cast<Vector<T>*>( &y )->close();

    int ierr = 0;
    const VectorPetscMPI<T>* vecxPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &x );
    const VectorPetscMPI<T>* vecyPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &y );
    if( vecxPetscMPI && vecyPetscMPI )
    {
        VectorPetscMPI<T>* vecOutPetscMPIRange =  dynamic_cast<VectorPetscMPIRange<T>*>( &(*this) );
        if ( !vecOutPetscMPIRange )
        {
            const VectorPetscMPIRange<T>* vecxPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &x );
            const VectorPetscMPIRange<T>* vecyPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &y );
            if ( !vecxPetscMPIRange && !vecyPetscMPIRange )
            {
                Vec lx;
                ierr = VecGhostGetLocalForm(this->vec(),&lx);
                CHKERRABORT( this->comm(),ierr );
                Vec lxInx;
                ierr = VecGhostGetLocalForm(vecxPetscMPI->vec(),&lxInx);
                CHKERRABORT( this->comm(),ierr );
                Vec lxIny;
                ierr = VecGhostGetLocalForm(vecyPetscMPI->vec(),&lxIny);
                CHKERRABORT( this->comm(),ierr );

                ierr =  VecPointwiseDivide(lx, lxInx, lxIny );
                CHKERRABORT( this->comm(),ierr );

                ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
                CHKERRABORT( this->comm(),ierr );
                ierr = VecGhostRestoreLocalForm(vecxPetscMPI->vec(),&lxInx);
                CHKERRABORT( this->comm(),ierr );
                ierr = VecGhostRestoreLocalForm(vecyPetscMPI->vec(),&lxIny);
                CHKERRABORT( this->comm(),ierr );
                return;
            }
        }

        this->pointwiseOperationOthersPetscImpl(x,y,1);
        return;
    }

    // default method
    super::pointwiseDivide( x,y );
}

template <typename T>
void
VectorPetscMPI<T>::pointwiseOperationOthersPetscImpl( Vector<T> const& x, Vector<T> const& y, int op )
{
    int ierr = 0;
    const VectorPetscMPIRange<T>* vecxPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &x );
    const VectorPetscMPIRange<T>* vecyPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &y );
    const VectorPetscMPI<T>* vecxPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &x );
    const VectorPetscMPI<T>* vecyPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &y );
    VectorPetscMPIRange<T>* vecOutPetscMPIRange =  dynamic_cast<VectorPetscMPIRange<T>*>( &(*this) );
    CHECK( vecxPetscMPI && vecyPetscMPI ) << "is not petsc vectors";

    // active part
    Vec thevecActiveInx = (vecxPetscMPIRange)? vecxPetscMPIRange->vec() : vecxPetscMPI->vec();
    Vec thevecActiveIny = (vecyPetscMPIRange)? vecyPetscMPIRange->vec() : vecyPetscMPI->vec();
    if ( op == 0 )
        ierr = VecPointwiseMult(this->vec(), thevecActiveInx, thevecActiveIny );
    else if ( op == 1 )
        ierr = VecPointwiseDivide(this->vec(), thevecActiveInx, thevecActiveIny );
    CHKERRABORT( this->comm(),ierr );

    // ghosts part
    size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
    size_type nLocalGhost = this->map().nLocalGhosts();

    PetscScalar *valuesOut;
    PetscScalar *valuesInxGhost;
    PetscScalar *valuesInyGhost;
    Vec lxOut = 0, lxInx = 0, lxIny = 0;
    size_type startOutGhostDof = 0, startxGhostDof = 0, startyGhostDof = 0;

    // get array in each context
    if ( vecOutPetscMPIRange )
    {
        ierr = VecGetArray( vecOutPetscMPIRange->vecGhost(), &valuesOut );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = VecGhostGetLocalForm(this->vec(),&lxOut);
        CHKERRABORT( this->comm(),ierr );
        CHECK( lxOut ) << "aie";
        ierr = VecGetArray( lxOut, &valuesOut );
        CHKERRABORT( this->comm(),ierr );
        startOutGhostDof = nLocalDofWithoutGhost;
    }
    if ( vecxPetscMPIRange )
    {
        ierr = VecGetArray( vecxPetscMPIRange->vecGhost(), &valuesInxGhost );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = VecGhostGetLocalForm(vecxPetscMPI->vec(),&lxInx);
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGetArray( lxInx, &valuesInxGhost );
        CHKERRABORT( this->comm(),ierr );
        startxGhostDof = nLocalDofWithoutGhost;
    }
    if ( vecyPetscMPIRange )
    {
        ierr = VecGetArray( vecyPetscMPIRange->vecGhost(), &valuesInyGhost );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = VecGhostGetLocalForm(vecyPetscMPI->vec(),&lxIny);
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGetArray( lxIny, &valuesInyGhost );
        CHKERRABORT( this->comm(),ierr );
        startyGhostDof = nLocalDofWithoutGhost;
    }

    // apply op on ghosts
    if ( op == 0)
    {
        // apply pointwiseMult on ghosts
        for ( size_type k=0;k<nLocalGhost;++k )
            valuesOut[startOutGhostDof+k] = valuesInxGhost[startxGhostDof+k]*valuesInyGhost[startyGhostDof+k];
    }
    else if ( op == 1 )
    {
        // apply pointwiseDiv on ghosts
        for ( size_type k=0;k<nLocalGhost;++k )
            valuesOut[startOutGhostDof+k] = valuesInxGhost[startxGhostDof+k]/valuesInyGhost[startyGhostDof+k];
    }

    // clean petsc object
    if ( vecOutPetscMPIRange )
    {
        ierr = VecRestoreArray( vecOutPetscMPIRange->vecGhost(), &valuesOut );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = VecRestoreArray( lxOut, &valuesOut );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGhostRestoreLocalForm(this->vec(),&lxOut);
        CHKERRABORT( this->comm(),ierr );
    }
    if ( vecxPetscMPIRange )
    {
        ierr = VecRestoreArray( vecxPetscMPIRange->vecGhost(), &valuesInxGhost );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = VecRestoreArray( lxInx, &valuesInxGhost );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGhostRestoreLocalForm(vecxPetscMPI->vec(),&lxInx);
        CHKERRABORT( this->comm(),ierr );
    }
    if ( vecyPetscMPIRange )
    {
        ierr = VecRestoreArray( vecyPetscMPIRange->vecGhost(), &valuesInyGhost );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        ierr = VecRestoreArray( lxIny, &valuesInyGhost );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGhostRestoreLocalForm(vecyPetscMPI->vec(),&lxIny);
        CHKERRABORT( this->comm(),ierr );
    }

}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::zero()
{
    if ( !this->closed() )
        this->close();
    int ierr = 0;
    //this->close();
    Vec lx;
    ierr = VecGhostGetLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    PetscScalar z=0.;
    ierr = VecSet ( lx, z );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::clear()
{
    if ( this->isInitialized() )
    {
        super::clear();
    }
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void VectorPetscMPI<T>::localize()
{
    if ( !this->isInitialized() )
        return;

    int ierr = 0;
    ierr = VecGhostUpdateBegin(this->vec(),INSERT_VALUES,SCATTER_FORWARD);
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostUpdateEnd(this->vec(),INSERT_VALUES,SCATTER_FORWARD);
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::close()
{
    tic();
    //FEELPP_ASSERT (this->isInitialized()).error( "VectorPetsc<> not initialized" );
    //std::cout << "\n MPI CLOSE "<<std::endl;;
    super::close();

    this->localize();
    toc("VectorPetscMPI::close",FLAGS_v>0);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
int
VectorPetscMPI<T>::reciprocal()
{
    if ( !this->closed() )
        this->close();
    int ierr = 0;
    Vec lx;
    ierr = VecGhostGetLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    ierr = VecReciprocal( lx );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecGhostRestoreLocalForm(this->vec(),&lx);
    CHKERRABORT( this->comm(),ierr );
    return ierr;
}


template <typename T>
typename VectorPetscMPI<T>::size_type
VectorPetscMPI<T>::firstLocalIndex() const
{
    DCHECK( this->isInitialized() ) << "VectorPetsc<> not initialized";

    int petsc_first=0;

    return static_cast<size_type>( petsc_first );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
typename VectorPetscMPI<T>::size_type
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
                const size_type original_dofProcess = original_dofCluster-original_firstDofCluster;
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
            dofClusterMissing[k]=it_dof->template get<0>();
            originaldofClusterMissing[k]=it_dof->template get<1>();
        }
    }

    auto worldCommFusion = this->map().worldComm()+vecInput.map().worldComm();
    std::vector<rank_type> globalRankToFusionRank_this(this->map().worldComm().globalSize());
    mpi::all_gather( this->map().worldComm().globalComm(),
                     worldCommFusion->globalRank(),
                     globalRankToFusionRank_this );
    std::vector<rank_type> globalRankToFusionRank_input(vecInput.map().worldComm().globalSize());
    mpi::all_gather( vecInput.map().worldComm().globalComm(),
                     worldCommFusion->globalRank(),
                     globalRankToFusionRank_input );

    std::vector<int> thisProcIsActive_fusion(worldCommFusion->globalSize());
    mpi::all_gather( worldCommFusion->globalComm(),
                     (int)this->map().worldComm().isActive(),
                     thisProcIsActive_fusion );
    std::vector<int> inputProcIsActive_fusion(worldCommFusion->globalSize());
    mpi::all_gather( worldCommFusion->globalComm(),
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
            if ( this->map().worldComm().globalRank() == proc  && thisProcIsActive_fusion[worldCommFusion->globalRank()] )  // send info to rankLocalization
            {
                const int rankToSend = globalRankToFusionRank_input[rankLocalization];
                worldCommFusion->globalComm().send(rankToSend,0,originaldofClusterMissing );
            }
            else if ( vecInput.map().worldComm().globalRank()==rankLocalization && inputProcIsActive_fusion[worldCommFusion->globalRank()] )
            {
                const int rankToRecv = globalRankToFusionRank_this[proc];
                worldCommFusion->globalComm().recv(rankToRecv,0,originaldofClusterMissing_recv );

                const size_type nDataRecv = originaldofClusterMissing_recv.size();
                dofClusterMissing_RequestVal.resize(nDataRecv);
                dofClusterMissing_RequestIsFind.resize(nDataRecv);
                for (size_type k=0;k<nDataRecv;++k)
                {
                    const size_type original_firstDofCluster = vecInput.map().firstDofGlobalCluster();
                    const size_type original_dofCluster = originaldofClusterMissing_recv[k];
                    if (original_dofCluster >=vecInput.map().firstDofGlobalCluster() && original_dofCluster<=vecInput.map().lastDofGlobalCluster())
                    {
                        const size_type original_dofProcess = original_dofCluster-original_firstDofCluster;
                        dofClusterMissing_RequestVal[k]=vecInput(original_dofProcess);
                        dofClusterMissing_RequestIsFind[k]=1;
                    }
                    else dofClusterMissing_RequestIsFind[k]=0;
                }
                worldCommFusion->globalComm().send( rankToRecv, 1, dofClusterMissing_RequestVal );
                worldCommFusion->globalComm().send( rankToRecv, 2, dofClusterMissing_RequestIsFind );
            }

            if ( this->map().worldComm().globalRank() == proc && thisProcIsActive_fusion[worldCommFusion->globalRank()]  )
            {
                const int rankToRecv = globalRankToFusionRank_input[rankLocalization];
                worldCommFusion->globalComm().recv( rankToRecv, 1, dofClusterMissing_RequestVal );
                worldCommFusion->globalComm().recv( rankToRecv, 2, dofClusterMissing_RequestIsFind );

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
            worldCommFusion->globalComm().barrier();
            //---------------------------------------
        }
    }
}

//----------------------------------------------------------------------------------------------------//


template <typename T>
typename VectorPetscMPI<T>::size_type
VectorPetscMPI<T>::localSize() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "VectorPetsc not initialized" );

#if 0
    int petsc_size=0;
    int ierr = VecGetLocalSize( M_vecLocal, &petsc_size );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<size_type>( petsc_size );
#else
    return this->map().nLocalDofWithGhost();
#endif
}



template <typename T>
VectorPetscMPIRange<T>::VectorPetscMPIRange( datamap_ptrtype const& dm )
    :
    super_type( dm, false ),
    M_destroyVecGhostOnExit( true ),
    M_destroyVecScatterGhostOnExit( true )
{
    const PetscScalar* arrayActive = NULL;
    const PetscScalar* arrayGhost = NULL;
    this->initRangeView( arrayActive,arrayGhost );
}

template<typename T>
VectorPetscMPIRange<T>::VectorPetscMPIRange( Vec v, datamap_ptrtype const& dm, bool duplicate )
    :
    super_type( dm, false ),
    M_destroyVecGhostOnExit( true ),
    M_destroyVecScatterGhostOnExit( true )
{

    this->M_destroy_vec_on_exit = duplicate;

    int ierr = 0;
    if ( duplicate )
    {
        ierr = VecDuplicate( v, &this->M_vec );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecCopy( v, this->M_vec );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        this->M_vec = v;
    }

    PetscInt petsc_n_localWithoutGhost=static_cast<PetscInt>( this->map().nLocalDofWithoutGhost() );
    PetscInt petsc_n_localGhost = static_cast<PetscInt>( this->map().nLocalGhosts() );
    // create ghost vector
    ierr = VecCreateSeq( PETSC_COMM_SELF, petsc_n_localGhost, &M_vecGhost );
    CHKERRABORT( this->comm(),ierr );

    // create ghosts mapping
    PetscInt *idxGhost = NULL;
    if ( petsc_n_localGhost > 0 )
    {
        idxGhost = new PetscInt[petsc_n_localGhost];
        std::copy( this->map().mapGlobalProcessToGlobalCluster().begin()+petsc_n_localWithoutGhost,
                   this->map().mapGlobalProcessToGlobalCluster().end(),
                   idxGhost );
    }

    // create vecScatter
    IS isScatter;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral( this->comm(), petsc_n_localGhost, idxGhost, PETSC_COPY_VALUES, &isScatter );
#else
    ierr = ISCreateGeneral( this->comm(), petsc_n_localGhost, idxGhost, &isScatter );
#endif
    CHKERRABORT( this->comm(),ierr );

    IS isLoc;
    ierr = ISCreateStride( PETSC_COMM_SELF,petsc_n_localGhost,0,1,&isLoc );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecScatterCreate( this->vec(), isScatter,
                             M_vecGhost, isLoc,
                             &M_vecScatterGhost );
    CHKERRABORT( this->comm(),ierr );

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy ( &isScatter );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( &isLoc );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy ( isScatter );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( isLoc );
    CHKERRABORT( this->comm(),ierr );
#endif
    delete[] idxGhost;

#if 0
    ierr = VecSetFromOptions( M_vecGhost );
    CHKERRABORT( this->comm(),ierr );
#endif
    this->setInitialized( true );
    this->localize();
    this->setIsClosed( true );
}

template<typename T>
VectorPetscMPIRange<T>::VectorPetscMPIRange( Vec v, Vec vGhost, datamap_ptrtype const& dm, bool duplicate )
    :
    super_type( dm, false ),
    M_destroyVecGhostOnExit( true ),
    M_destroyVecScatterGhostOnExit( true )
{
    this->M_destroy_vec_on_exit = duplicate;

    int ierr = 0;
    if ( duplicate )
    {
        ierr = VecDuplicate( v, &this->M_vec );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecCopy( v, this->M_vec );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecDuplicate( vGhost, &this->M_vecGhost );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecCopy( vGhost, this->M_vecGhost );
        CHKERRABORT( this->comm(),ierr );
    }
    else
    {
        this->M_vec = v;
        this->M_vecGhost = vGhost;
        M_destroyVecGhostOnExit = false;
    }

    this->initVecScatterGhost();

    this->setInitialized( true );
    this->localize();
    this->setIsClosed( true );
}

template<typename T>
VectorPetscMPIRange<T>::VectorPetscMPIRange( Vec v, Vec vGhost, VecScatter vecScatterGhost, datamap_ptrtype const& dm )
    :
    super_type( dm, false )
{
    this->M_destroy_vec_on_exit = false;
    this->M_destroyVecGhostOnExit = false;
    this->M_destroyVecScatterGhostOnExit = false;

    this->M_vec = v;
    this->M_vecGhost = vGhost;
    this->M_vecScatterGhost = vecScatterGhost;

    this->setInitialized( true );
    this->localize();
    this->setIsClosed( true );
}

template<typename T>
void
VectorPetscMPIRange<T>::init( datamap_ptrtype const& dm )
{
    if ( !this->map().isCompatible( *dm ) )
        this->setMap( dm );

    const PetscScalar* arrayActive = NULL;
    const PetscScalar* arrayGhost = NULL;
    this->initRangeView( arrayActive,arrayGhost );
}

template <typename T>
void
VectorPetscMPIRange<T>::initRangeView( const PetscScalar arrayActive[], const PetscScalar arrayGhost[] )
{
    // clear initialized vectors
    if ( this->isInitialized() )
        this->clear();

    int ierr=0;
    PetscInt petsc_n_dof=static_cast<PetscInt>( this->map().nDof() );
    PetscInt petsc_n_localWithoutGhost=static_cast<PetscInt>( this->map().nLocalDofWithoutGhost() );
    PetscInt petsc_n_localGhost=static_cast<PetscInt>( this->map().nLocalGhosts() );

    // active dof view
    if ( arrayActive )
        ierr = VecCreateMPIWithArray( this->comm(),1, petsc_n_localWithoutGhost, petsc_n_dof, arrayActive, &this->M_vec );
    else
        ierr = VecCreateMPI( this->comm(), petsc_n_localWithoutGhost, petsc_n_dof, &this->M_vec );
    CHKERRABORT( this->comm(),ierr );

    // ghost dof view
    if ( arrayGhost )
        ierr = VecCreateSeqWithArray( PETSC_COMM_SELF,1,petsc_n_localGhost,arrayGhost, &M_vecGhost );
    else
        ierr = VecCreateSeq( PETSC_COMM_SELF, petsc_n_localGhost, &M_vecGhost );
    CHKERRABORT( this->comm(),ierr );

    this->initVecScatterGhost();

    ierr = VecSetFromOptions( this->vec() );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecSetFromOptions( M_vecGhost );
    CHKERRABORT( this->comm(),ierr );

    this->setInitialized( true );
    this->setIsClosed( true );
}

template <typename T>
void
VectorPetscMPIRange<T>::initVecScatterGhost()
{

    M_destroyVecScatterGhostOnExit = false;


    int ierr=0;
    PetscInt petsc_n_localWithoutGhost=static_cast<PetscInt>( this->map().nLocalDofWithoutGhost() );
    PetscInt petsc_n_localGhost=static_cast<PetscInt>( this->map().nLocalGhosts() );
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

    // create vecScatter
    IS isScatter;
    PetscInt* idxGhost = ( petsc_n_localGhost > 0 )? std::addressof( idx[petsc_n_localWithoutGhost] ) : NULL;

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral( this->comm(), petsc_n_localGhost, idxGhost, PETSC_COPY_VALUES, &isScatter );
#else
    ierr = ISCreateGeneral( this->comm(), petsc_n_localGhost, idxGhost, &isScatter );
#endif
    CHKERRABORT( this->comm(),ierr );

    IS isLoc;
    ierr = ISCreateStride( PETSC_COMM_SELF,petsc_n_localGhost,0,1,&isLoc );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecScatterCreate( this->vec(), isScatter,
                             M_vecGhost, isLoc,
                             &M_vecScatterGhost );
    CHKERRABORT( this->comm(),ierr );


    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy( &is );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISLocalToGlobalMappingDestroy( &isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( &isScatter );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( &isLoc );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy( is );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISLocalToGlobalMappingDestroy( isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( isScatter );
    CHKERRABORT( this->comm(),ierr );
    ierr = ISDestroy ( isLoc );
    CHKERRABORT( this->comm(),ierr );
#endif

    delete[] idx;
}
template <typename T>
void
VectorPetscMPIRange<T>::clear()
{
    super_type::clear();

    if ( !this->isInitialized() )
        return;

    int ierr=0;

    if ( M_destroyVecGhostOnExit )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = VecDestroy( &M_vecGhost );
        CHKERRABORT( this->comm(),ierr );
#else
        ierr = VecDestroy( M_vecGhost );
        CHKERRABORT( this->comm(),ierr );
#endif
    }

    if ( M_destroyVecScatterGhostOnExit )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = VecScatterDestroy( &M_vecScatterGhost );
        CHKERRABORT( this->comm(),ierr );
#else
        ierr = VecScatterDestroy( M_vecScatterGhost );
        CHKERRABORT( this->comm(),ierr );
#endif
    }
}

template <typename T>
const typename VectorPetscMPIRange<T>::value_type&
VectorPetscMPIRange<T>::operator() ( const size_type i ) const
{
    int ierr=0;
    const PetscScalar *values;
    if ( i < this->map().nLocalDofWithoutGhost() )
    {
        ierr = VecGetArrayRead( this->vec(), &values );
        CHKERRABORT( this->comm(),ierr );
        const PetscScalar& value = values[i];
        ierr = VecRestoreArrayRead( this->vec(), &values );
        CHKERRABORT( this->comm(),ierr );
        return static_cast<const value_type&>( value );
    }
    else
    {
        ierr = VecGetArrayRead( M_vecGhost, &values );
        CHKERRABORT( this->comm(),ierr );
        const PetscScalar& value = values[i-this->map().nLocalDofWithoutGhost()];
        ierr = VecRestoreArrayRead( M_vecGhost, &values );
        CHKERRABORT( this->comm(),ierr );
        return static_cast<const value_type&>( value );
    }
}

template <typename T>
typename VectorPetscMPIRange<T>::value_type&
VectorPetscMPIRange<T>::operator() ( const size_type i )
{
    int ierr=0;
    PetscScalar *values;
    if ( i < this->map().nLocalDofWithoutGhost() )
    {
        ierr = VecGetArray( this->vec(), &values );
        CHKERRABORT( this->comm(),ierr );
        PetscScalar& value =  values[i];
        ierr = VecRestoreArray( this->vec(), &values );
        CHKERRABORT( this->comm(),ierr );
        return static_cast<value_type&>( value );
    }
    else
    {
        ierr = VecGetArray( M_vecGhost, &values );
        CHKERRABORT( this->comm(),ierr );
        PetscScalar& value =  values[i-this->map().nLocalDofWithoutGhost()];
        ierr = VecRestoreArray( M_vecGhost, &values );
        CHKERRABORT( this->comm(),ierr );
        return static_cast<value_type&>( value );
    }
}

template <typename T>
Vector<typename VectorPetscMPIRange<T>::value_type>&
VectorPetscMPIRange<T>::operator= ( const Vector<value_type> &V )
{
    if ( !this->closed() )
        this->close();
    if ( !V.closed() )
        const_cast<Vector<T>*>( &V )->close();

    int ierr=0;
    const VectorPetscMPIRange<T>* vecPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &V );
    if ( vecPetscMPIRange )
    {
        if ( !this->map().isCompatible( V.map() ) )
            this->setMap( V.mapPtr() );

        ierr = VecCopy( vecPetscMPIRange->vec(), this->vec() );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecCopy( vecPetscMPIRange->vecGhost(), this->vecGhost() );
        CHKERRABORT( this->comm(),ierr );
        return *this;
    }
    const VectorPetscMPI<T>* vecPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &V );
    if ( vecPetscMPI )
    {
        if ( !this->map().isCompatible( V.map() ) )
            this->setMap( V.mapPtr() );

        ierr = VecCopy( vecPetscMPI->vec(), this->vec() );
        CHKERRABORT( this->comm(),ierr );

        size_type nLocalDofWithoutGhost = vecPetscMPI->map().nLocalDofWithoutGhost();
        size_type nLocalGhost = vecPetscMPI->map().nLocalGhosts();

        Vec lxIn;
        ierr = VecGhostGetLocalForm(vecPetscMPI->vec(),&lxIn);
        CHKERRABORT( this->comm(),ierr );

        PetscScalar *valuesIn;
        PetscScalar *valuesOutGhost;
        ierr = VecGetArray( lxIn, &valuesIn );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGetArray( this->vecGhost(), &valuesOutGhost );
        CHKERRABORT( this->comm(),ierr );
        for ( size_type k=0;k<nLocalGhost;++k )
            valuesOutGhost[k] = valuesIn[nLocalDofWithoutGhost+k];

        ierr = VecRestoreArray( lxIn, &valuesIn );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecRestoreArray( this->vecGhost(), &valuesOutGhost );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecGhostRestoreLocalForm(vecPetscMPI->vec(),&lxIn);
        CHKERRABORT( this->comm(),ierr );

        return *this;
    }

    return super_type::operator= ( V );
}

template <typename T>
Vector<typename VectorPetscMPIRange<T>::value_type>&
VectorPetscMPIRange<T>::operator= ( const VectorPetscMPIRange<value_type> &V )
{
    return this->operator=( *dynamic_cast< Vector<value_type> const* >( &V ) );
}

template <typename T>
void
VectorPetscMPIRange<T>::set( const value_type& value )
{
    if ( !this->closed() )
        this->close();
    int ierr=0;
    const PetscScalar val = static_cast<PetscScalar>( value );
    ierr = VecSet ( this->vec(), val );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecSet ( M_vecGhost, val );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetscMPIRange<T>::add( const value_type& v_in )
{
    if ( !this->closed() )
        this->close();
    int ierr=0;
    const PetscScalar v = static_cast<PetscScalar>( v_in );
    ierr = VecShift( this->vec(), v );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecShift( M_vecGhost, v );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
void
VectorPetscMPIRange<T>::add( const value_type& a_in, const Vector<value_type>& v_in )
{
    if ( !this->closed() )
        this->close();
    if ( !v_in.closed() )
        const_cast<Vector<T>*>( &v_in )->close();

    int ierr = 0;
    PetscScalar a = static_cast<PetscScalar>( a_in );

    const VectorPetscMPIRange<T>* vecPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &v_in );
    if ( vecPetscMPIRange )
    {
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
        ierr = VecAXPY( &a, vecPetscMPIRange->vec(), this->vec() );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecAXPY( &a, vecPetscMPIRange->vecGhost(), this->vecGhost() );
        CHKERRABORT( this->comm(),ierr );

#else
        ierr = VecAXPY( this->vec(), a, vecPetscMPIRange->vec() );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecAXPY( this->vecGhost(), a, vecPetscMPIRange->vecGhost() );
        CHKERRABORT( this->comm(),ierr );
#endif
        return;
    }
    const VectorPetscMPI<T>* vecPetscMPI =  dynamic_cast<const VectorPetscMPI<T>*>( &v_in );
    if ( vecPetscMPI )
    {
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
        ierr = VecAXPY( &a, vecPetscMPI->vec(), this->vec() );
        CHKERRABORT( this->comm(),ierr );
#else
        ierr = VecAXPY( this->vec(), a, vecPetscMPI->vec() );
        CHKERRABORT( this->comm(),ierr );
#endif

        size_type nLocalDofWithoutGhost = vecPetscMPI->map().nLocalDofWithoutGhost();
        size_type nLocalGhost = vecPetscMPI->map().nLocalGhosts();

        Vec lxIn;
        ierr = VecGhostGetLocalForm(vecPetscMPI->vec(),&lxIn);
        CHKERRABORT( this->comm(),ierr );

        PetscScalar *valuesIn;
        PetscScalar *valuesOutGhost;
        ierr = VecGetArray( lxIn, &valuesIn );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecGetArray( this->vecGhost(), &valuesOutGhost );
        CHKERRABORT( this->comm(),ierr );
        for ( size_type k=0;k<nLocalGhost;++k )
            valuesOutGhost[k] += a*valuesIn[nLocalDofWithoutGhost+k];

        ierr = VecRestoreArray( lxIn, &valuesIn );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecRestoreArray( this->vecGhost(), &valuesOutGhost );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecGhostRestoreLocalForm(vecPetscMPI->vec(),&lxIn);
        CHKERRABORT( this->comm(),ierr );

        return;
    }

    super_type::add( a_in, v_in );
}

template <typename T>
void
VectorPetscMPIRange<T>::pointwiseMult( Vector<T> const& x, Vector<T> const& y )
{
    if ( !this->closed() )
        this->close();
    if ( !x.closed() )
        const_cast<Vector<T>*>( &x )->close();
    if ( !y.closed() )
        const_cast<Vector<T>*>( &y )->close();
    int ierr = 0;
    const VectorPetscMPIRange<T>* vecxPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &x );
    const VectorPetscMPIRange<T>* vecyPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &y );
    if ( vecxPetscMPIRange && vecyPetscMPIRange )
    {
        ierr = VecPointwiseMult(this->vec(), vecxPetscMPIRange->vec(), vecyPetscMPIRange->vec() );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecPointwiseMult(this->vecGhost(), vecxPetscMPIRange->vecGhost(), vecyPetscMPIRange->vecGhost() );
        CHKERRABORT( this->comm(),ierr );
        return;
    }

    super_type::pointwiseMult( x,y );
}
template <typename T>
void
VectorPetscMPIRange<T>::pointwiseDivide( Vector<T> const& x, Vector<T> const& y )
{
    if ( !this->closed() )
        this->close();
    if ( !x.closed() )
        const_cast<Vector<T>*>( &x )->close();
    if ( !y.closed() )
        const_cast<Vector<T>*>( &y )->close();
    int ierr = 0;
    const VectorPetscMPIRange<T>* vecxPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &x );
    const VectorPetscMPIRange<T>* vecyPetscMPIRange =  dynamic_cast<const VectorPetscMPIRange<T>*>( &y );
    if ( vecxPetscMPIRange && vecyPetscMPIRange )
    {
        ierr = VecPointwiseDivide(this->vec(), vecxPetscMPIRange->vec(), vecyPetscMPIRange->vec() );
        CHKERRABORT( this->comm(),ierr );
        ierr = VecPointwiseDivide(this->vecGhost(), vecxPetscMPIRange->vecGhost(), vecyPetscMPIRange->vecGhost() );
        CHKERRABORT( this->comm(),ierr );
        return;
    }
    super_type::pointwiseDivide( x,y );
}


template <typename T>
void
VectorPetscMPIRange<T>::zero()
{
    if ( !this->closed() )
        this->close();
    int ierr = 0;
    PetscScalar z=0.;
    ierr = VecSet( this->vec(), z );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecSet( M_vecGhost, z );
    CHKERRABORT( this->comm(),ierr );
}

template <typename T>
int
VectorPetscMPIRange<T>::reciprocal()
{
    if ( !this->closed() )
        this->close();
    int ierr = 0;
    ierr = VecReciprocal( this->vec() );
    CHKERRABORT( this->comm(),ierr );
    ierr = VecReciprocal( M_vecGhost );
    CHKERRABORT( this->comm(),ierr );
    return ierr;
}

template <typename T>
void
VectorPetscMPIRange<T>::localize()
{
    if ( !this->closed() )
        this->close();

    int ierr = 0;

    // Perform the scatter
    ierr = VecScatterBegin( M_vecScatterGhost, this->vec(), M_vecGhost, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecScatterEnd( M_vecScatterGhost, this->vec(), M_vecGhost, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT( this->comm(),ierr );
}



#if BOOST_VERSION < 105900
template<typename SizeT>
vector_ptrtype
vec( Vec v, datamap_ptrtype<SizeT> datamap )
{
    if ( datamap->worldComm().localSize() > 1 )
        return std::make_shared<Feel::VectorPetscMPI<double>>( v, datamap );
    else
        return std::make_shared<Feel::VectorPetsc<double>>( v, datamap );
}
#else
template<typename SizeT>
vector_uptrtype
vec( Vec v, datamap_ptrtype<SizeT> datamap )
{
    if ( datamap->worldComm().localSize() > 1 )
        return std::make_unique<Feel::VectorPetscMPI<double>>( v, datamap );
    else
        return std::make_unique<Feel::VectorPetsc<double>>( v, datamap );
    // using vector_ptrtype = std::shared_ptr<Feel::Vector<double> >;
}
#endif

template vector_uptrtype vec<uint32_type>( Vec v, datamap_ptrtype<uint32_type> datamap );
template class VectorPetsc<double>;
template class VectorPetscMPI<double>;
template class VectorPetscMPIRange<double>;

} // Feel


#endif // FEELPP_HAS_PETSC_H
