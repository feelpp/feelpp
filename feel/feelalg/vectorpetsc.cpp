/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-02
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>

#if defined( FEELPP_HAS_PETSC_H )




namespace Feel
{
/**
 * \p Utility::iota is a duplication of the SGI STL extension
 * \p std::iota.  It simply assigns sequentially increasing values
 * to a range. That is, it assigns \p value to \p *first, \p value + 1
 * to \p *(first + 1) and so on. In general, each iterator \p i in the
 * range [first, last) is assigned \p value + (i - \p first).
 */
template <typename ForwardIter, typename T>
void iota ( ForwardIter first, ForwardIter last, T value )
{
    while ( first != last )
    {
        *first = value++;
        ++first;
    }
}

template <typename T>
void
VectorPetsc<T>::zero()
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "VectorPetsc<> not initialized" );

    int ierr=0;

    PetscScalar z=0.;
    this->close();
    // 2.2.x & earlier style
#if PETSC_VERSION_LESS_THAN(2,2,0)
    ierr = VecSet ( &z, _M_vec );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = VecSet ( _M_vec, z );
    CHKERRABORT( this->comm(),ierr );
#endif
}


template <typename T>
void
VectorPetsc<T>::clear ()
{
    if ( ( this->isInitialized() ) && ( this->_M_destroy_vec_on_exit ) )
    {
        int ierr=0;

        ierr = PETSc::VecDestroy( _M_vec );
        CHKERRABORT( this->comm(),ierr );
    }

    this->M_is_closed = this->M_is_initialized = false;
}

template <typename T>
void
VectorPetsc<T>::insert ( const Vector<T>& /*V*/,
                         const std::vector<size_type>& /*dof_indices*/ )
{
    FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
}


template <typename T>
void
VectorPetsc<T>::insert ( const ublas::vector<T>& /*V*/,
                         const std::vector<size_type>& /*dof_indices*/ )
{
    FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
}

template <typename T>
void
VectorPetsc<T>::scale ( T factor_in )
{
    int ierr = 0;
    PetscScalar factor = static_cast<PetscScalar>( factor_in );

    // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = VecScale( &factor, _M_vec );
    CHKERRABORT( this->comm(),ierr );

    // 2.3.x & later style
#else

    ierr = VecScale( _M_vec, factor );
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
        ierr = VecGetArray ( _M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        int ig = fli + i;

        PetscScalar value = ( values[ig] + v );

        ierr = VecRestoreArray ( _M_vec, &values );
        CHKERRABORT( this->comm(),ierr );

        ierr = VecSetValues ( _M_vec, 1, &ig, &value, INSERT_VALUES );
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

    const VectorPetsc<T>* v = dynamic_cast<const VectorPetsc<T>*>( &v_in );

    assert ( v != NULL );
    assert( this->size() == v->size() );

    // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = VecAXPY( &a, v->_M_vec, _M_vec );
    CHKERRABORT( this->comm(),ierr );

    // 2.3.x & later style
#else

    ierr = VecAXPY( _M_vec, a, v->_M_vec );
    CHKERRABORT( this->comm(),ierr );

#endif
}


template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::min () const
{
    assert ( this->isInitialized() );

    int index=0, ierr=0;
    PetscReal min=0.;

    ierr = VecMin ( _M_vec, &index, &min );
    CHKERRABORT( this->comm(),ierr );

    // this return value is correct: VecMin returns a PetscReal
    return static_cast<Real>( min );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::max() const
{
    assert ( this->isInitialized() );

    int index=0, ierr=0;
    PetscReal max=0.;

    ierr = VecMax ( _M_vec, &index, &max );
    CHKERRABORT( this->comm(),ierr );

    // this return value is correct: VecMax returns a PetscReal
    return static_cast<Real>( max );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>:: l1Norm () const
{
    assert( this->closed() );

    int ierr=0;
    double value=0.;

    ierr = VecNorm ( _M_vec, NORM_1, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::l2Norm () const
{
    assert( this->closed() );

    int ierr=0;
    double value=0.;

    ierr = VecNorm ( _M_vec, NORM_2, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::linftyNorm () const
{
    assert( this->closed() );

    int ierr=0;
    double value=0.;

    ierr = VecNorm ( _M_vec, NORM_INFINITY, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
typename VectorPetsc<T>::value_type
VectorPetsc<T>:: sum () const
{
    assert( this->closed() );

    int ierr=0;
    double value=0.;

    ierr = VecSum ( _M_vec, &value );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<Real>( value );
}

template <typename T>
void VectorPetsc<T>::printMatlab ( const std::string name ) const
{
    assert ( this->isInitialized() );
    FEELPP_ASSERT ( this->closed() ).warn( "vector is not closed" );

    if ( !this->closed() )
    {
        Debug() << "closing vector\n";
        const_cast<VectorPetsc<T>*>( this )->close();
    }

    const_cast<VectorPetsc<T>*>( this )->close();

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

        ierr = VecView ( const_cast<Vec>( _M_vec ), petsc_viewer );
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

        ierr = VecView ( const_cast<Vec>( _M_vec ), PETSC_VIEWER_STDOUT_WORLD );
        CHKERRABORT( this->comm(),ierr );
    }


    /**
     * Destroy the viewer.
     */
    ierr = PETSc::PetscViewerDestroy ( petsc_viewer );
    CHKERRABORT( this->comm(),ierr );
}


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//

template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( DataMap const& dm )
    :
    super( dm,false ) //false for not init
{
    this->init( dm.nDof(), dm.nLocalDofWithoutGhost() );
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( Vec v, DataMap const& dm )
    :
    super( v,dm )
{
    ///HERE!!!
    int ierr=0;
    int petsc_n_localWithGhost=static_cast<int>( this->map().nLocalDofWithGhost()/*n_local*/ );

    ierr = VecCreateSeq ( PETSC_COMM_SELF, petsc_n_localWithGhost, &  _M_vecLocal );
    CHKERRABORT( this->comm(),ierr );
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
    int petsc_n_localWithGhost=static_cast<int>( this->map().nLocalDofWithGhost()/*n_local*/ );
    //std::cout << "petsc_n_localWithoutGhost "<< petsc_n_localWithoutGhost << std::endl;
    //std::cout << "petsc_n_localWithGhost "<< petsc_n_localWithGhost << std::endl;

    // Clear initialized vectors
    if ( this->isInitialized() )
        this->clear();

    FEELPP_ASSERT( n_localWithoutGhost < n )( n_localWithoutGhost )( n ).warn( "invalid local size" );

    ierr = VecCreateMPI ( this->comm(), petsc_n_localWithoutGhost, petsc_n,
                          &this->_M_vec );
    CHKERRABORT( this->comm(),ierr );

    //ierr = VecSetFromOptions (this->vec());
    //CHKERRABORT(this->comm(),ierr);

    // localToGlobalMapping
    IS is;
    ISLocalToGlobalMapping isLocToGlobMap;

    //auto idx = this->map().mapGlobalProcessToGlobalCluster();
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

    ierr=ISLocalToGlobalMappingCreateIS( is, &isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );

    ierr=VecSetLocalToGlobalMapping( this->vec(),isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );

    // local vector
    ierr = VecCreateSeq ( PETSC_COMM_SELF, petsc_n_localWithGhost, &  _M_vecLocal );
    CHKERRABORT( this->comm(),ierr );

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy( &is );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISLocalToGlobalMappingDestroy( &isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy( is );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISLocalToGlobalMappingDestroy( isLocToGlobMap );
    CHKERRABORT( this->comm(),ierr );
#endif

    delete idx;

    ierr = VecSetFromOptions( this->vec() );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecSetFromOptions( _M_vecLocal );
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
    ierr = VecGetArray( _M_vecLocal, &values );
    CHKERRABORT( this->comm(),ierr );
    //std::cout << "\n operator MPI ";
    value =  values[i /*- this->firstLocalIndex()*/ ];

    ierr = VecRestoreArray( _M_vecLocal, &values );
    CHKERRABORT( this->comm(),ierr );

    return static_cast<value_type>( value );
}
//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::set( size_type i, const value_type& value )
{
    //FEELPP_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int ierr=0;
    int i_val = static_cast<int>( i );
    PetscScalar petsc_value = static_cast<PetscScalar>( value );

    ierr=VecSetValuesLocal( this->vec(),1,&i_val,&petsc_value,INSERT_VALUES );
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
    //FEELPP_ASSERT(n<=size())( n )( size() ).error( "invalid local index array size" );

    int ierr=0;
    ierr=VecSetValuesLocal( this->vec(), n, i, v, ADD_VALUES );
    CHKERRABORT( this->comm(),ierr );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::clear()
{
    super::clear();

    if ( /*(this->isInitialized()) &&*/ ( this->destroy_vec_on_exit() ) )
    {
        int ierr=0;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = VecDestroy( &_M_vecLocal );
#else
        ierr = VecDestroy( _M_vecLocal );
#endif
        CHKERRABORT( this->comm(),ierr );

    }
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void VectorPetscMPI<T>::localize()
{
    //std::cout << "\n MPI LOCALIZE "<<std::endl;;

    int ierr = 0;

    IS isGlob;
    IS isLoc;

    VecScatter scatter;
#if 0
    auto idx = this->map().mapGlobalProcessToGlobalCluster();
    ierr = ISCreateGeneral( this->comm(), idx.size(), &idx[0], PETSC_COPY_VALUES, &isGlob );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISCreateStride( PETSC_COMM_SELF,idx.size(),0,1,&isLoc );
    CHKERRABORT( this->comm(),ierr );
#else
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
#endif

    // create scatter
    ierr = VecScatterCreate( this->vec(), isGlob,
                             _M_vecLocal, isLoc,
                             &scatter );
    CHKERRABORT( this->comm(),ierr );

    // Perform the scatter
    ierr = VecScatterBegin( scatter, this->vec(), _M_vecLocal, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecScatterEnd  ( scatter, this->vec(), _M_vecLocal, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT( this->comm(),ierr );

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy ( &isGlob );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISDestroy ( &isLoc );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecScatterDestroy( &scatter );
    CHKERRABORT( this->comm(),ierr );
#else
    ierr = ISDestroy ( isGlob );
    CHKERRABORT( this->comm(),ierr );

    ierr = ISDestroy ( isLoc );
    CHKERRABORT( this->comm(),ierr );

    ierr = VecScatterDestroy( scatter );
    CHKERRABORT( this->comm(),ierr );
#endif

    delete idx;
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
    assert ( this->isInitialized() );

    int petsc_first=0;

    return static_cast<size_type>( petsc_first );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
size_type
VectorPetscMPI<T>::lastLocalIndex() const
{
    assert ( this->isInitialized() );

    int petsc_last=this->map().nLocalDofWithGhost();

    return static_cast<size_type>( petsc_last );
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::duplicateFromOtherPartition( Vector<T> const& vecInput)
{
#if defined(FEELPP_ENABLE_MPI_MODE) // WITH MPI

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

#endif // MPI
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::duplicateFromOtherPartition_run( Vector<T> const& vecInput)
{
#if defined(FEELPP_ENABLE_MPI_MODE) // WITH MPI

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
    std::vector<int> globalRankToFusionRank_this(this->map().worldComm().globalSize());
    mpi::all_gather( this->map().worldComm().globalComm(),
                     worldCommFusion.globalRank(),
                     globalRankToFusionRank_this );
    std::vector<int> globalRankToFusionRank_input(vecInput.map().worldComm().globalSize());
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


#endif // MPI

}

//----------------------------------------------------------------------------------------------------//

template class VectorPetsc<double>;
template class VectorPetscMPI<double>;

} // Feel

#endif // FEELPP_HAS_PETSC_H
