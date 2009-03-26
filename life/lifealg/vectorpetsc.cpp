/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-02

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

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
#include <life/lifecore/life.hpp>
#include <life/lifealg/vectorpetsc.hpp>
#include <life/lifealg/matrixpetsc.hpp>

#if defined( HAVE_PETSC_H )




namespace Life
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
VectorPetsc<T>::insert (const Vector<T>& /*V*/,
                        const std::vector<size_type>& /*dof_indices*/)
{
    LIFE_ASSERT( 0 ).error( "invalid call, not implemented yet" );
}


template <typename T>
void
VectorPetsc<T>::insert (const ublas::vector<T>& /*V*/,
                        const std::vector<size_type>& /*dof_indices*/)
{
    LIFE_ASSERT( 0 ).error( "invalid call, not implemented yet" );
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
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // 2.3.x & later style
#else

    ierr = VecScale(_M_vec, factor);
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
            CHKERRABORT(Application::COMM_WORLD,ierr);

            int ig = fli + i;

            PetscScalar value = (values[ig] + v);

            ierr = VecRestoreArray (_M_vec, &values);
            CHKERRABORT(Application::COMM_WORLD,ierr);

            ierr = VecSetValues (_M_vec, 1, &ig, &value, INSERT_VALUES);
            CHKERRABORT(Application::COMM_WORLD,ierr);
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
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // 2.3.x & later style
#else

    ierr = VecAXPY(_M_vec, a, v->_M_vec);
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    CHKERRABORT(Application::COMM_WORLD,ierr);

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
    std::vector<int> idx(n); Life::iota (idx.begin(), idx.end(), 0);

    // Create the index set & scatter object
    ierr = ISCreateGeneral(Application::COMM_WORLD, n, &idx[0], &is);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = VecScatterCreate(const_cast<Vec>(this->_M_vec), is,
                            v_local->_M_vec, is,
                            &scatter);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Perform the scatter
#if ( (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) ) || ( PETSC_VERSION_MAJOR >= 3 )

    ierr = VecScatterBegin(scatter, _M_vec, v_local->_M_vec, INSERT_VALUES,
                           SCATTER_FORWARD );
#else
    ierr = VecScatterBegin( _M_vec, v_local->_M_vec, INSERT_VALUES,
                            SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3)|| ( PETSC_VERSION_MAJOR >= 3 )
    ierr = VecScatterEnd  (scatter, const_cast<Vec>(_M_vec), v_local->_M_vec, INSERT_VALUES,
                           SCATTER_FORWARD );
#else
    ierr = VecScatterEnd  ( _M_vec, v_local->_M_vec, INSERT_VALUES,
                            SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Clean up
    ierr = ISDestroy (is);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = VecScatterDestroy(scatter);
    CHKERRABORT(Application::COMM_WORLD,ierr);
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

    // Create the index set & scatter object
    ierr = ISCreateGeneral(Application::COMM_WORLD, n_sl, &idx[0], &is);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = VecScatterCreate(const_cast<Vec>(_M_vec),          is,
                            v_local->_M_vec, is,
                            &scatter);

    CHKERRABORT(Application::COMM_WORLD,ierr);


    // Perform the scatter
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
    ierr = VecScatterBegin(scatter, const_cast<Vec>(_M_vec), v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD );
#else
    ierr = VecScatterBegin( _M_vec, v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
    ierr = VecScatterEnd  (scatter, const_cast<Vec>(_M_vec), v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD);
#else
    ierr = VecScatterEnd  ( _M_vec, v_local->_M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter );
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Clean up
    ierr = ISDestroy (is);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = VecScatterDestroy(scatter);
    CHKERRABORT(Application::COMM_WORLD,ierr);
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
        Life::iota (idx.begin(), idx.end(), first_local_idx);

        // Create the index set & scatter object
        ierr = ISCreateGeneral(Application::COMM_WORLD, local_size, &idx[0], &is);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        ierr = VecScatterCreate(_M_vec,              is,
                                parallel_vec._M_vec, is,
                                &scatter);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        // Perform the scatter
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
        ierr = VecScatterBegin(scatter, _M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD );
#else
        ierr = VecScatterBegin(_M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter);
#endif
        CHKERRABORT(Application::COMM_WORLD,ierr);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3) || ( PETSC_VERSION_MAJOR >= 3 )
        ierr = VecScatterEnd  (scatter, _M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD );
#else
        ierr = VecScatterEnd  (_M_vec, parallel_vec._M_vec, INSERT_VALUES, SCATTER_FORWARD, scatter);
#endif
        CHKERRABORT(Application::COMM_WORLD,ierr);

        // Clean up
        ierr = ISDestroy (is);
        CHKERRABORT(Application::COMM_WORLD,ierr);

        ierr = VecScatterDestroy(scatter);
        CHKERRABORT(Application::COMM_WORLD,ierr);
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
            CHKERRABORT(Application::COMM_WORLD,ierr);

            for (int i=0; i<n; i++)
                v_local[i] = static_cast<double>(values[i]);

            ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
            CHKERRABORT(Application::COMM_WORLD,ierr);
        }

    // otherwise multiple processors
    else
        {
            size_type ioff = firstLocalIndex();
            std::vector<double> local_values(n, 0.);

            {
                ierr = VecGetArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(Application::COMM_WORLD,ierr);

                for (int i=0; i<nl; i++)
                    local_values[i+ioff] = static_cast<double>(values[i]);

                ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(Application::COMM_WORLD,ierr);
            }

            MPI_Allreduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
                           Application::COMM_WORLD);
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
            CHKERRABORT(Application::COMM_WORLD,ierr);

            for (int i=0; i<n; i++)
                v_local[i] = static_cast<Complex>(values[i]);

            ierr = VecRestoreArray (_M_vec, &values);
            CHKERRABORT(Application::COMM_WORLD,ierr);
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
                CHKERRABORT(Application::COMM_WORLD,ierr);

                // provide my local share to the real and imag buffers
                for (int i=0; i<nl; i++)
                    {
                        real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
                        imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
                    }

                ierr = VecRestoreArray (_M_vec, &values);
                CHKERRABORT(Application::COMM_WORLD,ierr);
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
                           MPI_DOUBLE, MPI_SUM, Application::COMM_WORLD);

            MPI_Allreduce (&imag_local_values[0], &imag_v_local[0], n,
                           MPI_DOUBLE, MPI_SUM, Application::COMM_WORLD);

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
            CHKERRABORT(Application::COMM_WORLD,ierr);

            for (int i=0; i<n; i++)
                v_local[i] = static_cast<Real>(values[i]);

            ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
            CHKERRABORT(Application::COMM_WORLD,ierr);
        }

    // otherwise multiple processors
    else
        {
            size_type ioff = this->firstLocalIndex();
            std::vector<Real> local_values (n, 0.);

            {
                ierr = VecGetArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(Application::COMM_WORLD,ierr);

                for (int i=0; i<nl; i++)
                    local_values[i+ioff] = static_cast<Real>(values[i]);

                ierr = VecRestoreArray (const_cast<Vec>(_M_vec), &values);
                CHKERRABORT(Application::COMM_WORLD,ierr);
            }


            MPI_Reduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
                        pid, Application::COMM_WORLD);
        }
}

template <typename T>
void VectorPetsc<T>::printMatlab (const std::string name) const
{
    assert (this->isInitialized());
    LIFE_ASSERT (this->closed()).warn( "vector is not closed" );

    if ( !this->closed() )
        {
            Debug() << "closing vector\n";
            const_cast<VectorPetsc<T>*>(this)->close();
        }

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

            ierr = VecView (const_cast<Vec>(_M_vec), petsc_viewer);
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

            ierr = VecView (const_cast<Vec>(_M_vec), PETSC_VIEWER_STDOUT_WORLD);
            CHKERRABORT(Application::COMM_WORLD,ierr);
        }


    /**
     * Destroy the viewer.
     */
    ierr = PetscViewerDestroy (petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);
}




template class VectorPetsc<double>;

} // Life

#endif // HAVE_PETSC_H
