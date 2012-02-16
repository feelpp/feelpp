/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-10-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file vectorpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-10-18
 */
#ifndef __VectorPetsc_H
#define __VectorPetsc_H 1

#include <feel/feelconfig.h>

#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>

#if defined(HAVE_PETSC_H)
#include <feel/feelcore/application.hpp>


extern "C"
{
#if defined(MPICH_NAME)
#if !defined( MPICH_HAVE_MPI_WIN )
#define MPICH_HAVE_MPI_WIN
  struct MPI_Win {};
#endif
#endif
#include <petscmat.h>
}



namespace Feel
{
template<typename T> class MatrixPetsc;

/**
 * \class VectorPetsc
 * \brief Wrapper for petsc matrices
 *
 * Petsc vector. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 *
 * @author Benjamin S. Kirk, 2002
 * @author Christophe Prud'homme
 * @see
 */
template<typename T>
class VectorPetsc : public Vector<T>
{
    typedef Vector<T> super;

public:


    /** @name Typedefs
     */
    //@{

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef typename super::clone_ptrtype clone_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  Dummy-Constructor. Dimension=0
     */
    VectorPetsc ()
        :
        super(),
        M_comm(),
        _M_destroy_vec_on_exit( true )
    {
    }

    /**
     * Constructor. Set dimension to \p n and initialize all elements with zero.
     */
    VectorPetsc (const size_type n)
        :
        super( n ),
        M_comm(),
        _M_destroy_vec_on_exit( true )
    {
        this->init(n, n, false);
    }

    /**
     * Constructor. Set local dimension to \p n_local, the global dimension
     * to \p n, and initialize all elements with zero.
     */
    VectorPetsc (const size_type n,
                 const size_type n_local)
        :
        super( n, n_local ),
        M_comm(),
        _M_destroy_vec_on_exit( true )
    {
        this->init(n, n_local, false);
    }

    VectorPetsc ( DataMap const& dm, bool doInit=true )
        :
        super(dm),
        M_comm(),
        _M_destroy_vec_on_exit( true )
    {
        if (doInit)
            this->init(dm.nDof(), dm.nLocalDofWithoutGhost(), false);
    }


    /**
     * Constructor.  Creates a VectorPetsc assuming you already have a
     * valid PETSc Vec object.  In this case, v is NOT destroyed by the
     * VectorPetsc constructor when this object goes out of scope.
     * This allows ownership of v to remain with the original creator,
     * and to simply provide additional functionality with the VectorPetsc.
     */
    VectorPetsc(Vec v)
        :
        super(),
        M_comm(),
        _M_destroy_vec_on_exit( false )
    {
        this->_M_vec = v;
        this->M_is_initialized = true;
    }

    VectorPetsc(Vec v, DataMap const& dm)
        :
        super(dm),
        M_comm(),
        _M_destroy_vec_on_exit( false )
    {
        this->_M_vec = v;
        this->M_is_initialized = true;
    }

    /**
     * Destructor, deallocates memory. Made virtual to allow
     * for derived classes to behave properly.
     */
    ~VectorPetsc ()
    {
        this->clear();
    }

    /**
     * Creates a copy of this vector and returns it in an \p
     * shared_ptr<>.  This must be overloaded in the derived classes.
     */
    clone_ptrtype clone () const
    {
        clone_ptrtype cloned_vector (new VectorPetsc<T>);

        *cloned_vector = *this;

        return cloned_vector;
    }


    /**
     * Change the dimension of the vector to \p N. The reserved memory for
     * this vector remains unchanged if possible, to make things faster, but
     * this may waste some memory, so take this in the back of your head.
     * However, if \p N==0 all memory is freed, i.e. if you want to resize
     * the vector and release the memory not needed, you have to first call
     * \p init(0) and then \p init(N). This cited behaviour is analogous
     * to that of the STL containers.
     *
     * On \p fast==false, the vector is filled by
     * zeros.
     */
    void init (const size_type N,
               const size_type n_local,
               const bool         fast=false);

    /**
     * call init with n_local = N,
     */
    void init (const size_type N,
               const bool         fast=false)
    {
        this->init(N,N,fast);
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    value_type operator() (const size_type i) const
    {
        FEEL_ASSERT (this->isInitialized()).error( "vector not initialized" );
        FEEL_ASSERT ( ((i >= this->firstLocalIndex()) &&
                       (i <  this->lastLocalIndex())) )( i )( this->firstLocalIndex() )( this->lastLocalIndex() ).error( "invalid vector index" );

        int ierr=0;
        PetscScalar *values, value=0.;


        ierr = VecGetArray(_M_vec, &values);
        CHKERRABORT(M_comm,ierr);

        value = values[i - this->firstLocalIndex()];

        ierr = VecRestoreArray (_M_vec, &values);
        CHKERRABORT(M_comm,ierr);

        return static_cast<value_type>(value);
    }


    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    Vector<T> & operator += (const Vector<value_type> &V)
    {
        FEEL_ASSERT(this->closed()).error( "vector is not closed" );

        this->add(1., V);

        return *this;
    }

    /**
     * Subtraction operator.
     * Fast equivalent to \p U.add(-1, V).
     */
    Vector<T> & operator -= (const Vector<value_type> &V)
    {
        FEEL_ASSERT(this->closed()).error( "vector is not closed" );

        this->add(-1., V);

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    const bool destroy_vec_on_exit() const {return _M_destroy_vec_on_exit;}

    /**
     * @return dimension of the vector. This
     * function was formerly called \p n(), but
     * was renamed to get the \p PetscVector<T> class
     * closer to the C++ standard library's
     * \p std::vector container.
     */
    size_type size () const
    {
        FEEL_ASSERT (this->isInitialized()).error( "VectorPetsc not initialized" );


        if (!this->isInitialized())
            return 0;

        int petsc_size=0;
        int ierr = VecGetSize(_M_vec, &petsc_size);
        CHKERRABORT(M_comm,ierr);
        return static_cast<size_type>(petsc_size);
    }
    /**
     * @return the local size of the vector
     * (index_stop-index_start)
     */
    size_type localSize() const
    {
        FEEL_ASSERT (this->isInitialized()).error( "VectorPetsc not initialized" );

        int petsc_size=0;
        int ierr = VecGetLocalSize(_M_vec, &petsc_size);
        CHKERRABORT(M_comm,ierr);

        return static_cast<size_type>(petsc_size);
    }

    /**
     * Returns the raw PETSc vector context pointer.  Note this is generally
     * not required in user-level code. Just don't do anything crazy like
     * calling VecDestroy()!
     */
    Vec vec () const { FEEL_ASSERT (_M_vec != 0).error( "invalid petsc vector" ); return _M_vec; }
    Vec& vec ()  { FEEL_ASSERT (_M_vec != 0).error( "invalid petsc vector" ); return _M_vec; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Call the assemble functions
     */
    void close ()
    {
        FEEL_ASSERT (this->isInitialized()).error( "VectorPetsc<> not initialized" );

        int ierr=0;

        ierr = VecAssemblyBegin(_M_vec);
        CHKERRABORT(M_comm,ierr);
        ierr = VecAssemblyEnd(_M_vec);
        CHKERRABORT(M_comm,ierr);

        this->M_is_closed = true;
    }

    /**
     * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
     * faster.
     */
    void zero ()
    {
        FEEL_ASSERT (this->isInitialized()).error( "VectorPetsc<> not initialized" );

        int ierr=0;

        PetscScalar z=0.;

        // 2.2.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

        ierr = VecSet (&z, _M_vec);
        CHKERRABORT(M_comm,ierr);

        // 2.3.x & newer
#else

        ierr = VecSet (_M_vec, z);
        CHKERRABORT(M_comm,ierr);

#endif
    }

    /**
     * Set entries to zero between \p start and \p stop
     */
    void zero ( size_type /*start*/,  size_type /*stop*/ )
    {
        this->zero();
    }

    /**
     * set the entries to the constant \p v
     */
    void setConstant( value_type v ) { this->set( v ); }

    /**
     * @returns the \p VectorPetsc<T> to a pristine state.
     */
    FEEL_DONT_INLINE void clear ();

    /**
     * \f$ v(i) = \mathrm{value} \forall i\f$
     */
    void set ( const value_type& value);

    /**
     * v(i) = value
     */
    void set ( size_type i, const value_type& value);

    /**
     * v(i) += value
     */
    void add ( size_type i, const value_type& value);

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    void addVector ( int* i, int n, value_type* v );

    /**
     * \f$ U+=v \f$ where \p v is a std::vector<T>
     * and you
     * want to specify WHERE to add it
     */
    void addVector (const std::vector<value_type>& v,
                    const std::vector<size_type>& dof_indices)
    {
        FEEL_ASSERT (v.size() == dof_indices.size()).error( "invalid dof indices" );

        for (size_type i=0; i<v.size(); i++)
            this->add (dof_indices[i], v[i]);
    }

    /**
     * \f$ U+=V \f$ where U and V are type
     * \p NumericVector<T> and you
     * want to specify WHERE to add
     * the \p NumericVector<T> V
     */
    void addVector (const Vector<value_type>& V,
                    const std::vector<size_type>& dof_indices)
    {
        FEEL_ASSERT (V.size() == dof_indices.size()).error( "invalid dof indices" );

        for (size_type i=0; i<V.size(); i++)
            this->add (dof_indices[i], V(i));
    }


    /**
     * \f$ U+=A*V\f$, add the product of a \p MatrixSparse \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector (const Vector<value_type>& V_in,
                    const MatrixSparse<value_type>& A_in)
    {
        const VectorPetsc<T>* V = dynamic_cast<const VectorPetsc<T>*>(&V_in);
        const MatrixPetsc<T>* A = dynamic_cast<const MatrixPetsc<T>*>(&A_in);

        assert (V != 0);
        assert (A != 0);

        int ierr=0;

        A->close();

        // The const_cast<> is not elegant, but it is required since PETSc
        // is not const-correct.
        ierr = MatMultAdd(const_cast<MatrixPetsc<T>*>(A)->mat(), V->_M_vec, _M_vec, _M_vec);
        CHKERRABORT(M_comm,ierr);
    }


    /**
     * \f$U+=V \f$ where U and V are type
     * uvlas::vector<T> and you
     * want to specify WHERE to add
     * the DenseVector<T> V
     */
    void addVector (const ublas::vector<value_type>& V,
                    const std::vector<size_type>& dof_indices)
    {
        FEEL_ASSERT (V.size() == dof_indices.size()).error( "invalid dof indices" );

        for (size_type i=0; i<V.size(); i++)
            this->add (dof_indices[i], V(i));
    }

    /**
     * \f$ U=v \f$ where v is a DenseVector<T>
     * and you want to specify WHERE to insert it
     */
    void insert (const std::vector<T>& /*v*/,
                         const std::vector<size_type>& /*dof_indices*/)
    {
        FEEL_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * \f$U=V\f$, where U and V are type
     * Vector<T> and you
     * want to specify WHERE to insert
     * the Vector<T> V
     */
    void insert (const Vector<T>& V,
                 const std::vector<size_type>& dof_indices);


    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to insert
     * the DenseVector<T> V
     */
    void insert (const ublas::vector<T>& V,
                 const std::vector<size_type>& dof_indices);

    /**
     * Scale each element of the
     * vector by the given factor.
     */
    void scale (const T factor);


    /**
     * \f$ U(0-DIM)+=s\f$.
     * Addition of \p s to all components. Note
     * that \p s is a scalar and not a vector.
     */
    void add (const value_type& v_in );

    /**
     * \f$ U+=V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add (const Vector<value_type>& v);

    /**
     * \f$ U+=a*V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add (const value_type& a_in, const Vector<value_type>& v_in);

    /**
     * @return the minimum element in the vector.
     * In case of complex numbers, this returns the minimum
     * Real part.
     */
    real_type min () const;

    /**
     * @returns the maximum element in the vector.
     * In case of complex numbers, this returns the maximum
     * Real part.
     */
    real_type max() const;

    /**
     * @return the \f$l_1\f$-norm of the vector, i.e.
     * the sum of the absolute values.
     */
    real_type l1Norm () const;

    /**
     * @returns the \f$l_2\f$-norm of the vector, i.e.
     * the square root of the sum of the
     * squares of the elements.
     */
    real_type l2Norm () const;

    /**
     * @return the sum of the vector components, i.e.
     * the sum of the values.
     */
    value_type sum () const;

    /**
     * @returns the maximum absolute value of the
     * elements of this vector, which is the
     * \f$l_\infty\f$-norm of a vector.
     */
    real_type linftyNorm () const;

    /**
     * @return the index of the first vector element
     * actually stored on this processor
     */
    size_type firstLocalIndex () const
    {
        assert (this->isInitialized());

        int ierr=0, petsc_first=0, petsc_last=0;

        ierr = VecGetOwnershipRange (_M_vec, &petsc_first, &petsc_last);
        CHKERRABORT(M_comm,ierr);

        return static_cast<size_type>(petsc_first);
    }



    /**
     * @return the index of the last vector element
     * actually stored on this processor
     */
    size_type lastLocalIndex () const
    {
        assert (this->isInitialized());

        int ierr=0, petsc_first=0, petsc_last=0;

        ierr = VecGetOwnershipRange (_M_vec, &petsc_first, &petsc_last);
        CHKERRABORT(M_comm,ierr);

        return static_cast<size_type>(petsc_last);
    }

    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize (std::vector<T>& v_local) const;

    /**
     * Same, but fills a \p Vector<T> instead of
     * a \p std::vector.
     */
    void localize (Vector<T>& v_local) const;

    /**
     * Creates a local vector \p v_local containing
     * only information relevant to this processor, as
     * defined by the \p send_list.
     */
    void localize (Vector<T>& v_local,
                   const std::vector<size_type>& send_list) const;

    /**
     * Updates a local vector with selected values from neighboring
     * processors, as defined by \p send_list.
     */
    void localize (const size_type first_local_idx,
                   const size_type last_local_idx,
                   const std::vector<size_type>& send_list);

    /**
     * Creates a local copy of the global vector in
     * \p v_local only on processor \p proc_id.  By
     * default the data is sent to processor 0.  This method
     * is useful for outputting data from one processor.
     */
    void localizeToOneProcessor (std::vector<T>& v_local,
                                 const size_type proc_id=0) const;


    /**
     * Print the contents of the vector in Matlab's format. Optionally
     *  prints the vector to the file named \p name.  If \p name is
     *  not specified it is dumped to the screen.
     */
    void printMatlab(const std::string name="NULL") const;

    //@}



protected:

public:

    // disable
    VectorPetsc( VectorPetsc const & v)
        :
        super( v ),
        _M_destroy_vec_on_exit( true )
    {
        FEEL_ASSERT(v.closed()).error( "copied vector is not closed" );

        VecDuplicate(v._M_vec, &_M_vec);
        VecCopy(v._M_vec, _M_vec);
        this->M_is_initialized = true;
        this->close();
    }


private:

    mpi::communicator M_comm;

    /**
     * Petsc vector datatype to store values
     */
    Vec _M_vec;

    /**
     * This boolean value should only be set to false
     * for the constructor which takes a PETSc Vec object.
     */
    const bool _M_destroy_vec_on_exit;
};

template <typename T>
inline
void
VectorPetsc<T>::init (const size_type n,
                      const size_type n_local,
                      const bool fast)
{
    int ierr=0;
    int petsc_n=static_cast<int>(n);
    int petsc_n_local=static_cast<int>(n_local);


    // Clear initialized vectors
    if (this->isInitialized())
        this->clear();


    // create a sequential vector if on only 1 processor
    if (n_local == n)
        {
            ierr = VecCreateSeq (PETSC_COMM_SELF, petsc_n, &_M_vec);
            CHKERRABORT(PETSC_COMM_SELF,ierr);

            ierr = VecSetFromOptions (_M_vec);
            CHKERRABORT(PETSC_COMM_SELF,ierr);
        }
    // otherwise create an MPI-enabled vector
    else
        {
            FEEL_ASSERT(n_local < n)( n_local )( n ).error( "invalid local size" );

            ierr = VecCreateMPI (M_comm, petsc_n_local, petsc_n,
                                 &_M_vec);
            CHKERRABORT(M_comm,ierr);

            ierr = VecSetFromOptions (_M_vec);
            CHKERRABORT(M_comm,ierr);
        }

    this->M_is_initialized = true;


    if (fast == false)
        this->zero ();
}
template <typename T>
void
VectorPetsc<T>::set ( const value_type& value)
{
    int ierr=0;
    PetscScalar petsc_value = static_cast<PetscScalar>(value);

    ierr = VecSet (_M_vec, petsc_value );
    CHKERRABORT(M_comm,ierr);
}
template <typename T>
void
VectorPetsc<T>::set ( size_type i, const value_type& value)
{
    FEEL_ASSERT(i<size())( i )( size() ).error( "invalid index" );


    int ierr=0;
    int i_val = static_cast<int>(i);
    PetscScalar petsc_value = static_cast<PetscScalar>(value);

    ierr = VecSetValues (_M_vec, 1, &i_val, &petsc_value, INSERT_VALUES);
    CHKERRABORT(M_comm,ierr);
}

template <typename T>
void
VectorPetsc<T>::add (const size_type i, const value_type& value)
{
    FEEL_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int ierr=0;
    int i_val = static_cast<int>(i);
    PetscScalar petsc_value = static_cast<PetscScalar>(value);

    ierr = VecSetValues (_M_vec, 1, &i_val, &petsc_value, ADD_VALUES);
    CHKERRABORT(M_comm,ierr);
}

template <typename T>
void
VectorPetsc<T>::addVector ( int* i, int n, value_type* v )
{
    //FEEL_ASSERT(n<=size())( n )( size() ).error( "invalid local index array size" );

    int ierr=0;
    ierr = VecSetValues (_M_vec, n, i, v, ADD_VALUES);
    CHKERRABORT(M_comm,ierr);

}

//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//

template<typename T>
class VectorPetscMPI : public VectorPetsc<T>
{
    typedef VectorPetsc<T> super;
    typedef typename super::value_type value_type;

public:

    VectorPetscMPI()
        :
        super()
    {}

    VectorPetscMPI(Vec v, DataMap const& dm);

    VectorPetscMPI(DataMap const& dm );

    ~VectorPetscMPI() { this->clear(); }

    void init(const size_type N,
              const size_type n_local,
              const bool fast=false);

    value_type operator() (const size_type i) const;

    void set(size_type i, const value_type& value);

    void add(const size_type i, const value_type& value);

    void addVector(int* i, int n, value_type* v );

    void clear();

    void localize();

    void close();

    size_type firstLocalIndex() const;
    size_type lastLocalIndex() const;

private :

    Vec _M_vecLocal;

};

} // Feel
#endif /* HAVE_PETSC */
#endif /* __VectorPetsc_H */
