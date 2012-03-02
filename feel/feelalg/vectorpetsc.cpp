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

#if defined( HAVE_PETSC_H )




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
void iota (ForwardIter first, ForwardIter last, T value)
{
    while (first != last)
        {
            *first = value++;
            ++first;
        }
}

template <typename T>
void
VectorPetsc<T>::clear ()
{
    if ((this->isInitialized()) && (this->_M_destroy_vec_on_exit))
    {
        int ierr=0;

        ierr = PETSc::VecDestroy( _M_vec );
        CHKERRABORT(this->comm(),ierr);
    }

    this->M_is_closed = this->M_is_initialized = false;
}

template <typename T>
void
VectorPetsc<T>::insert (const Vector<T>& /*V*/,
                        const std::vector<size_type>& /*dof_indices*/)
{
    FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
}


template <typename T>
void
VectorPetsc<T>::insert (const ublas::vector<T>& /*V*/,
                        const std::vector<size_type>& /*dof_indices*/)
{
    FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
}

template <typename T>
void
VectorPetsc<T>::scale ( T factor_in )
{
    int ierr = 0;
    PetscScalar factor = static_cast<PetscScalar>(factor_in);

    // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = VecScale(&factor, _M_vec);
    CHKERRABORT(this->comm(),ierr);

    // 2.3.x & later style
#else

    ierr = VecScale(_M_vec, factor);
    CHKERRABORT(this->comm(),ierr);

#endif
}
template <typename T>
void
VectorPetsc<T>::add (const value_type& v_in )
{
    int ierr=0;
    PetscScalar* values;
    const PetscScalar v = static_cast<PetscScalar>(v_in);
    const int n   = static_cast<int>(this->localSize());
    const int fli = static_cast<int>(this->firstLocalIndex());

    for (int i=0; i<n; i++)
        {
            ierr = VecGetArray (_M_vec, &values);
            CHKERRABORT(this->comm(),ierr);

            int ig = fli + i;

            PetscScalar value = (values[ig] + v);

            ierr = VecRestoreArray (_M_vec, &values);
            CHKERRABORT(this->comm(),ierr);

            ierr = VecSetValues (_M_vec, 1, &ig, &value, INSERT_VALUES);
            CHKERRABORT(this->comm(),ierr);
        }
}
template <typename T>
void
VectorPetsc<T>::add (const Vector<value_type>& v)
{
    this->add (1., v);
}
template <typename T>
void
VectorPetsc<T>::add (const value_type& a_in, const Vector<value_type>& v_in)
{
    int ierr = 0;
    PetscScalar a = static_cast<PetscScalar>(a_in);

    const VectorPetsc<T>* v = dynamic_cast<const VectorPetsc<T>*>(&v_in);

    assert (v != NULL);
    assert(this->size() == v->size());

    // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = VecAXPY(&a, v->_M_vec, _M_vec);
    CHKERRABORT(this->comm(),ierr);

    // 2.3.x & later style
#else

    ierr = VecAXPY(_M_vec, a, v->_M_vec);
    CHKERRABORT(this->comm(),ierr);

#endif
}


template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::min () const
{
    assert (this->isInitialized());

    int index=0, ierr=0;
    PetscReal min=0.;

    ierr = VecMin (_M_vec, &index, &min);
    CHKERRABORT(this->comm(),ierr);

    // this return value is correct: VecMin returns a PetscReal
    return static_cast<Real>(min);
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::max() const
{
    assert (this->isInitialized());

    int index=0, ierr=0;
    PetscReal max=0.;

    ierr = VecMax (_M_vec, &index, &max);
    CHKERRABORT(this->comm(),ierr);

    // this return value is correct: VecMax returns a PetscReal
    return static_cast<Real>(max);
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>:: l1Norm () const
{
    assert(this->closed());

    int ierr=0;
    double value=0.;

    ierr = VecNorm (_M_vec, NORM_1, &value);
    CHKERRABORT(this->comm(),ierr);

    return static_cast<Real>(value);
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::l2Norm () const
{
    assert(this->closed());

    int ierr=0;
    double value=0.;

    ierr = VecNorm (_M_vec, NORM_2, &value);
    CHKERRABORT(this->comm(),ierr);

    return static_cast<Real>(value);
}

template <typename T>
typename VectorPetsc<T>::real_type
VectorPetsc<T>::linftyNorm () const
{
    assert(this->closed());

    int ierr=0;
    double value=0.;

    ierr = VecNorm (_M_vec, NORM_INFINITY, &value);
    CHKERRABORT(this->comm(),ierr);

    return static_cast<Real>(value);
}

template <typename T>
typename VectorPetsc<T>::value_type
VectorPetsc<T>:: sum () const
{
    assert(this->closed());

    int ierr=0;
    double value=0.;

    ierr = VecSum (_M_vec, &value);
    CHKERRABORT(this->comm(),ierr);

    return static_cast<Real>(value);
}

template <typename T>
void
VectorPetsc<T>::localize (Vector<T>& v_local_in) const
{
    VectorPetsc<T>* v_local = dynamic_cast<VectorPetsc<T>*>(&v_local_in);

    assert (v_local != NULL);
    assert (v_local->localSize() == this->size());

    int ierr = 0;
    const int n = this->size();

    IS is;
    VecScatter scatter;

    // Create idx, idx[i] = i;
    std::vector<int> idx(n); Feel::iota (idx.begin(), idx.end(), 0);

    // Create the index set & scatter object
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(), n, &idx[0], PETSC_COPY_VALUES, &is);
#else
    ierr = ISCreateGeneral(this->comm(), n, &idx[0], &is);
#endif
    CHKERRABORT(this->comm(),ierr);

    ierr = VecScatterCreate(const_cast<Vec>(this->_M_vec), is,
                            v_local->_M_vec, is,
                            &scatter);
    CHKERRABORT(this->comm(),ierr);

    // Perform the scatter
#if ( (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) ) || ( PETSC_VERSION_MAJOR >= 3 )

    ierr = VecScatterBegin(scatter, _M_vec, v_local->_M_vec, INSERT_VALUES,
                           SCATTER_FORWARD );
#else
    ierr = VecScatterBegin( _M_vec, v_local->_M_vec, INSERT_VALUES,
                            SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3)|| ( PETSC_VERSION_MAJOR >= 3 )
    ierr = VecScatterEnd  (scatter, const_cast<Vec>(_M_vec), v_local->_M_vec, INSERT_VALUES,
                           SCATTER_FORWARD );
#else
    ierr = VecScatterEnd  ( _M_vec, v_local->_M_vec, INSERT_VALUES,
                            SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(this->comm(),ierr);

    // Clean up
    ierr = PETSc::ISDestroy (is);
    CHKERRABORT(this->comm(),ierr);

    ierr = PETSc::VecScatterDestroy(scatter);
    CHKERRABORT(this->comm(),ierr);
}



template <typename T>
void VectorPetsc<T>::localize (Vector<T>& v_local_in,
                               const std::vector<size_type>& send_list) const
{
    VectorPetsc<T>* v_local = dynamic_cast<VectorPetsc<T>*>(&v_local_in);

    assert (v_local != NULL);
    assert (v_local->localSize() == this->size());
    assert (send_list.size()     <= v_local->size());

    int ierr=0;
    const int n_sl = send_list.size();

    IS is;
    VecScatter scatter;

    std::vector<int> idx(n_sl);

    for (int i=0; i<n_sl; i++)
        idx[i] = static_cast<int>(send_list[i]);

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(), n_sl, &idx[0], PETSC_COPY_VALUES, &is);
#else
    // Create the index set & scatter object
    ierr = ISCreateGeneral(this->comm(), n_sl, &idx[0], &is);
#endif
    CHKERRABORT(this->comm(),ierr);

    ierr = VecScatterCreate(const_cast<Vec>(_M_vec),          is,
                            v_local->_M_vec, is,
                            &scatter);

    CHKERRABORT(this->comm(),ierr);


    // Perform the scatter
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
    ierr = VecScatterBegin(scatter, const_cast<Vec>(_M_vec), v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD );
#else
    ierr = VecScatterBegin( _M_vec, v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
    ierr = VecScatterEnd  (scatter, const_cast<Vec>(_M_vec), v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD);
#else
    ierr = VecScatterEnd  ( _M_vec, v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(this->comm(),ierr);

    // Clean up
    ierr = PETSc::ISDestroy (is);
    CHKERRABORT(this->comm(),ierr);

    ierr = PETSc::VecScatterDestroy(scatter);
    CHKERRABORT(this->comm(),ierr);
}


template <typename T>
void VectorPetsc<T>::localize (const size_type first_local_idx,
                               const size_type last_local_idx,
                               const std::vector<size_type>& send_list)
{
    // Only good for serial vectors.
    assert (this->size() == this->localSize());
    assert (last_local_idx > first_local_idx);
    assert (send_list.size() <= this->size());
    assert (last_local_idx < this->size());

    const size_type size       = this->size();
    const size_type local_size = (last_local_idx - first_local_idx + 1);
    int ierr=0;

    // Don't bother for serial cases
    if ((first_local_idx == 0) &&
        (local_size == size))
        return;


    // Build a parallel vector, initialize it with the local
    // parts of (*this)
    VectorPetsc<T> parallel_vec;

    parallel_vec.init (size, local_size, false );


    // Copy part of *this into the parallel_vec
    {
        IS is;
        VecScatter scatter;

        // Create idx, idx[i] = i+first_local_idx;
        std::vector<int> idx(local_size);
        Feel::iota (idx.begin(), idx.end(), first_local_idx);

        // Create the index set & scatter object
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(), local_size, &idx[0], PETSC_COPY_VALUES, &is);
#else
        ierr = ISCreateGeneral(this->comm(), local_size, &idx[0], &is);
#endif
        CHKERRABORT(this->comm(),ierr);

        ierr = VecScatterCreate(_M_vec,              is,
                                parallel_vec._M_vec, is,
                                &scatter);
        CHKERRABORT(this->comm(),ierr);

        // Perform the scatter
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
        ierr = VecScatterBegin(scatter, _M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD );
#else
        ierr = VecScatterBegin(_M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter);
#endif
        CHKERRABORT(this->comm(),ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
        ierr = VecScatterEnd  (scatter, _M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD );
#else
        ierr = VecScatterEnd  (_M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter);
#endif
        CHKERRABORT(this->comm(),ierr);

        // Clean up
        ierr = PETSc::ISDestroy (is);
        CHKERRABORT(this->comm(),ierr);

        ierr = PETSc::VecScatterDestroy(scatter);
        CHKERRABORT(this->comm(),ierr);
    }

    // localize like normal
    parallel_vec.close();
    parallel_vec.localize (*this, send_list);
    this->close();
}

// Full specialization for double datatypes
template <>
void VectorPetsc<double>::localize (std::vector<double>& v_local) const
{
    int ierr=0;
    const int n  = this->size();
    const int nl = this->localSize();
    PetscScalar *values;


    v_local.resize(n);


    for (int i=0; i<n; i++)
        v_local[i] = 0.;

    // only one processor
    if (n == nl)
        {
            ierr = VecGetArray (const_cast<Vec>(_M_vec), &values);
            CHKERRABORT(this->comm(),ierr);

            for (int i=0; i<n; i++)
                v_local[i] = static_cast<double>(values[i]);

            ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
            CHKERRABORT(this->comm(),ierr);
        }

    // otherwise multiple processors
    else
        {
            size_type ioff = firstLocalIndex();
            std::vector<double> local_values(n, 0.);

            {
                ierr = VecGetArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(this->comm(),ierr);

                for (int i=0; i<nl; i++)
                    local_values[i+ioff] = static_cast<double>(values[i]);

                ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(this->comm(),ierr);
            }

            MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
                           this->comm());
        }
}



#if 0
// Full specialization for Complex datatypes
template <>
void VectorPetsc<Complex>::localize (std::vector<Complex>& v_local) const
{
    int ierr=0;
    const int n  = size();
    const int nl = localSize();
    PetscScalar *values;

    v_local.resize(n);


    for (int i=0; i<n; i++)
        v_local[i] = 0.;

    // only one processor
    if (n == nl)
        {
            ierr = VecGetArray (_M_vec, &values);
            CHKERRABORT(this->comm(),ierr);

            for (int i=0; i<n; i++)
                v_local[i] = static_cast<Complex>(values[i]);

            ierr = VecRestoreArray (_M_vec, &values);
            CHKERRABORT(this->comm(),ierr);
        }

    // otherwise multiple processors
    else
        {
            size_type ioff = firstLocalIndex();

            /* in here the local values are stored, acting as send buffer for MPI
             * initialize to zero, since we collect using MPI_SUM
             */
            std::vector<Real> real_local_values(n, 0.);
            std::vector<Real> imag_local_values(n, 0.);

            {
                ierr = VecGetArray (_M_vec, &values);
                CHKERRABORT(this->comm(),ierr);

                // provide my local share to the real and imag buffers
                for (int i=0; i<nl; i++)
                    {
                        real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
                        imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
                    }

                ierr = VecRestoreArray (_M_vec, &values);
                CHKERRABORT(this->comm(),ierr);
            }

            /* have buffers of the real and imaginary part of v_local.
             * Once MPI_Reduce() collected all the real and imaginary
             * parts in these std::vector<double>, the values can be
             * copied to v_local
             */
            std::vector<Real> real_v_local(n);
            std::vector<Real> imag_v_local(n);

            // collect entries from other proc's in real_v_local, imag_v_local
            MPI_Allreduce (&real_local_values[0], &real_v_local[0], n,
                           MPI_DOUBLE, MPI_SUM, this->comm());

            MPI_Allreduce (&imag_local_values[0], &imag_v_local[0], n,
                           MPI_DOUBLE, MPI_SUM, this->comm());

            // copy real_v_local and imag_v_local to v_local
            for (int i=0; i<n; i++)
                v_local[i] = Complex(real_v_local[i], imag_v_local[i]);

        }
}
#endif
// Full specialization for Real datatypes
template <>
void VectorPetsc<Real>::localizeToOneProcessor (std::vector<Real>& v_local,
                                                const size_type pid) const
{
    int ierr=0;
    const int n  = size();
    const int nl = localSize();
    PetscScalar *values;


    v_local.resize(n);


    // only one processor
    if (n == nl)
        {
            ierr = VecGetArray (const_cast<Vec>(_M_vec), &values);
            CHKERRABORT(this->comm(),ierr);

            for (int i=0; i<n; i++)
                v_local[i] = static_cast<Real>(values[i]);

            ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
            CHKERRABORT(this->comm(),ierr);
        }

    // otherwise multiple processors
    else
        {
            size_type ioff = this->firstLocalIndex();
            std::vector<Real> local_values (n, 0.);

            {
                ierr = VecGetArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(this->comm(),ierr);

                for (int i=0; i<nl; i++)
                    local_values[i+ioff] = static_cast<Real>(values[i]);

                ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(this->comm(),ierr);
            }


            MPI_Reduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
                        pid, this->comm());
        }
}

template <typename T>
void VectorPetsc<T>::printMatlab (const std::string name) const
{
    assert (this->isInitialized());
    FEELPP_ASSERT (this->closed()).warn( "vector is not closed" );

    if ( !this->closed() )
        {
            Debug() << "closing vector\n";
            const_cast<VectorPetsc<T>*>(this)->close();
        }

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
            CHKERRABORT(this->comm(),ierr);

            ierr = VecView (const_cast<Vec>(_M_vec), petsc_viewer);
            CHKERRABORT(this->comm(),ierr);
        }

    /**
     * Otherwise the matrix will be dumped to the screen.
     */
    else
        {
            ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                         PETSC_VIEWER_ASCII_MATLAB);
            CHKERRABORT(this->comm(),ierr);

            ierr = VecView (const_cast<Vec>(_M_vec), PETSC_VIEWER_STDOUT_WORLD);
            CHKERRABORT(this->comm(),ierr);
        }


    /**
     * Destroy the viewer.
     */
    ierr = PETSc::PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(this->comm(),ierr);
}


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//

template<typename T>
VectorPetscMPI<T>::VectorPetscMPI( DataMap const& dm )
    :
    super(dm,false)//false for not init
{
    this->init(dm.nDof(), dm.nLocalDofWithoutGhost());
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
VectorPetscMPI<T>::VectorPetscMPI(Vec v, DataMap const& dm)
    :
    super(v,dm)
{
    ///HERE!!!
    int ierr=0;
    int petsc_n_localWithGhost=static_cast<int>(this->map().nLocalDofWithGhost()/*n_local*/);

    ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n_localWithGhost, &  _M_vecLocal);
    CHKERRABORT(this->comm(),ierr);
    this->close();
}

//----------------------------------------------------------------------------------------------------//

template<typename T>
void
VectorPetscMPI<T>::init(const size_type n,
                        const size_type n_localWithoutGhost,
                        const bool fast)
{
    //std::cout << "MPI init start" << std::endl;
    int ierr=0;
    int petsc_n=static_cast<int>(n);
    int petsc_n_localWithoutGhost=static_cast<int>(n_localWithoutGhost);
    int petsc_n_localWithGhost=static_cast<int>(this->map().nLocalDofWithGhost()/*n_local*/);
    //std::cout << "petsc_n_localWithoutGhost "<< petsc_n_localWithoutGhost << std::endl;
    //std::cout << "petsc_n_localWithGhost "<< petsc_n_localWithGhost << std::endl;

    // Clear initialized vectors
    if (this->isInitialized())
        this->clear();

    FEELPP_ASSERT(n_localWithoutGhost < n)( n_localWithoutGhost )( n ).error( "invalid local size" );

    ierr = VecCreateMPI (this->comm(), petsc_n_localWithoutGhost, petsc_n,
                         &this->vec() );//&_M_vec);
    CHKERRABORT(this->comm(),ierr);

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
    ierr = ISCreateGeneral(this->comm(), n_idx, idx, PETSC_COPY_VALUES, &is);
#else
    ierr = ISCreateGeneral(this->comm(), n_idx, idx, &is);
#endif
    CHKERRABORT(this->comm(),ierr);

    ierr=ISLocalToGlobalMappingCreateIS(is, &isLocToGlobMap);
    CHKERRABORT(this->comm(),ierr);

    ierr=VecSetLocalToGlobalMapping(this->vec(),isLocToGlobMap);
    CHKERRABORT(this->comm(),ierr);

    // local vector
    ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n_localWithGhost, &  _M_vecLocal);
    CHKERRABORT(this->comm(),ierr);

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy(&is);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISLocalToGlobalMappingDestroy(&isLocToGlobMap);
    CHKERRABORT(this->comm(),ierr);
#else
    ierr = ISDestroy(is);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISLocalToGlobalMappingDestroy(isLocToGlobMap);
    CHKERRABORT(this->comm(),ierr);
#endif

    delete idx;

    ierr = VecSetFromOptions(this->vec());
    CHKERRABORT(this->comm(),ierr);

    ierr = VecSetFromOptions(_M_vecLocal);
    CHKERRABORT(this->comm(),ierr);

    this->M_is_initialized = true;

    if (fast == false)
        this->zero ();
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
typename VectorPetscMPI<T>::value_type
VectorPetscMPI<T>::operator() (const size_type i) const
{
        int ierr=0;
        PetscScalar *values, value=0.;
        ierr = VecGetArray(_M_vecLocal, &values);
        CHKERRABORT(this->comm(),ierr);
        //std::cout << "\n operator MPI ";
        value =  values[i /*- this->firstLocalIndex()*/ ];

        ierr = VecRestoreArray(_M_vecLocal, &values);
        CHKERRABORT(this->comm(),ierr);

        return static_cast<value_type>(value);
}
//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::set( size_type i, const value_type& value)
{
    //FEELPP_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int ierr=0;
    int i_val = static_cast<int>(i);
    PetscScalar petsc_value = static_cast<PetscScalar>(value);

    ierr=VecSetValuesLocal(this->vec(),1,&i_val,&petsc_value,INSERT_VALUES );
    CHKERRABORT(this->comm(),ierr);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::add (const size_type i, const value_type& value)
{
    //FEELPP_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int ierr=0;
    int i_val = static_cast<int>(i);
    PetscScalar petsc_value = static_cast<PetscScalar>(value);

    ierr=VecSetValuesLocal(this->vec(), 1, &i_val, &petsc_value, ADD_VALUES);
    CHKERRABORT(this->comm(),ierr);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::addVector ( int* i, int n, value_type* v )
{
    //FEELPP_ASSERT(n<=size())( n )( size() ).error( "invalid local index array size" );

    int ierr=0;
    ierr=VecSetValuesLocal(this->vec(), n, i, v, ADD_VALUES);
    CHKERRABORT(this->comm(),ierr);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
void
VectorPetscMPI<T>::clear()
{
    super::clear();

    if (/*(this->isInitialized()) &&*/ (this->destroy_vec_on_exit()))
    {
        int ierr=0;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = VecDestroy(&_M_vecLocal);
#else
        ierr = VecDestroy(_M_vecLocal);
#endif
        CHKERRABORT(this->comm(),ierr);

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
    ierr = ISCreateGeneral(this->comm(), idx.size(), &idx[0], PETSC_COPY_VALUES, &isGlob);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISCreateStride(PETSC_COMM_SELF,idx.size(),0,1,&isLoc);
    CHKERRABORT(this->comm(),ierr);
#else
    PetscInt *idx;
    PetscInt n_idx =  this->map().mapGlobalProcessToGlobalCluster().size();
    idx = new PetscInt[n_idx];
    std::copy( this->map().mapGlobalProcessToGlobalCluster().begin(),
               this->map().mapGlobalProcessToGlobalCluster().end(),
               idx );

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISCreateGeneral(this->comm(), n_idx, idx, PETSC_COPY_VALUES, &isGlob);
#else
    ierr = ISCreateGeneral(this->comm(), n_idx, idx, &isGlob);
#endif
    CHKERRABORT(this->comm(),ierr);

    ierr = ISCreateStride(PETSC_COMM_SELF,n_idx,0,1,&isLoc);
    CHKERRABORT(this->comm(),ierr);
#endif

    // create scatter
    ierr = VecScatterCreate(this->vec(), isGlob,
                            _M_vecLocal, isLoc,
                            &scatter);
    CHKERRABORT(this->comm(),ierr);

    // Perform the scatter
    ierr = VecScatterBegin(scatter, this->vec(), _M_vecLocal, INSERT_VALUES, SCATTER_FORWARD );
    CHKERRABORT(this->comm(),ierr);

    ierr = VecScatterEnd  (scatter, this->vec(), _M_vecLocal, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRABORT(this->comm(),ierr);

    // Clean up
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = ISDestroy (&isGlob);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISDestroy (&isLoc);
    CHKERRABORT(this->comm(),ierr);

    ierr = VecScatterDestroy(&scatter);
    CHKERRABORT(this->comm(),ierr);
#else
    ierr = ISDestroy (isGlob);
    CHKERRABORT(this->comm(),ierr);

    ierr = ISDestroy (isLoc);
    CHKERRABORT(this->comm(),ierr);

    ierr = VecScatterDestroy(scatter);
    CHKERRABORT(this->comm(),ierr);
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
    assert (this->isInitialized());

    int petsc_first=0, petsc_last=0;

    petsc_first=0; petsc_last=this->map().nLocalDofWithGhost();

    return static_cast<size_type>(petsc_first);
}

//----------------------------------------------------------------------------------------------------//

template <typename T>
size_type
VectorPetscMPI<T>::lastLocalIndex() const
{
    assert (this->isInitialized());

    int petsc_first=0, petsc_last=0;

    petsc_first=0; petsc_last=this->map().nLocalDofWithGhost();

    return static_cast<size_type>(petsc_last);
}

//----------------------------------------------------------------------------------------------------//

template class VectorPetsc<double>;
template class VectorPetscMPI<double>;

} // Feel

#endif // HAVE_PETSC_H
