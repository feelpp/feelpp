/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */
#ifndef __VectorPetsc_H
#define __VectorPetsc_H 1

#include <feel/feelconfig.h>

#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vectorublas.hpp>

BOOST_CLASS_EXPORT_KEY(Feel::VectorPetsc<double>)
BOOST_CLASS_EXPORT_KEY(Feel::VectorPetscMPI<double>)

#if defined(FEELPP_HAS_PETSC_H)
#include <feel/feelcore/application.hpp>

extern "C"
{
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
class FEELPP_EXPORT VectorPetsc : public Vector<T>
{
    typedef Vector<T> super;

public:

    friend class boost::serialization::access;

    /** @name Typedefs
     */
    //@{

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef typename super::clone_ptrtype clone_ptrtype;
    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;
    using size_type = typename super::size_type;
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
        M_destroy_vec_on_exit( true )
    {
    }

    /**
     * Constructor. Set dimension to \p n and initialize all elements with zero.
     */
    VectorPetsc ( const size_type n, worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() )
        :
        super( n, _worldComm ),
        M_destroy_vec_on_exit( true )
    {
        this->init( n, n, false );
    }

    /**
     * Constructor. Set local dimension to \p n_local, the global dimension
     * to \p n, and initialize all elements with zero.
     */
    VectorPetsc ( const size_type n,
                  const size_type n_local,
                  worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() )
        :
        super( n, n_local, _worldComm ),
        M_destroy_vec_on_exit( true )
    {
        this->init( n, n_local, false );
    }

    VectorPetsc ( datamap_ptrtype const& dm, bool doInit=true )
        :
        super( dm ),
        M_destroy_vec_on_exit( true )
    {
        if ( doInit )
            this->init( dm->nDof(), dm->nLocalDofWithoutGhost(), false );
    }


    /**
     * Constructor.  Creates a VectorPetsc assuming you already have a
     * valid PETSc Vec object.  In this case, v is NOT destroyed by the
     * VectorPetsc constructor when this object goes out of scope.
     * This allows ownership of v to remain with the original creator,
     * and to simply provide additional functionality with the VectorPetsc.
     */
    VectorPetsc( Vec v, bool duplicate = false )
        :
        super(),
        M_destroy_vec_on_exit( duplicate )
    {
        if ( duplicate )
        {
            VecDuplicate( v, &M_vec );
            VecCopy( v, M_vec );
        }
        else
            this->M_vec = v;
        this->M_is_initialized = true;
        this->setIsClosed( true );
    }

    VectorPetsc( Vec v, datamap_ptrtype const& dm, bool duplicate = false )
        :
        super( dm ),
        M_destroy_vec_on_exit( duplicate )
    {
        if ( duplicate )
        {
            VecDuplicate( v, &M_vec );
            VecCopy( v, M_vec );
        }
        else
            this->M_vec = v;
        this->M_is_initialized = true;
        this->setIsClosed( true );
    }

    /**
     * @brief creates a VectorPetsc out of a VectorUblas
     *
     * there is no copy, PETSc will use the storage of the VectorUblas
     */
    template<typename Storage>
    VectorPetsc( VectorUblas<T,Storage> const& v )
        :
        super( v.mapPtr() ),
        M_destroy_vec_on_exit( true )
        {
            int ierr=0;
            PetscInt petsc_n_dof=static_cast<PetscInt>( this->map().nDof() );
            PetscInt petsc_n_localWithoutGhost=static_cast<PetscInt>( this->map().nLocalDofWithoutGhost() );
            PetscInt petsc_n_localGhost=static_cast<PetscInt>( this->map().nLocalGhosts() );
            const PetscScalar* thearray = ( this->map().nLocalDofWithGhost() > 0 )? std::addressof( *v.begin()/*v[0]*/ ) : NULL;
            PetscInt *idx = NULL;
            if ( petsc_n_localGhost > 0 )
            {
                idx = new PetscInt[petsc_n_localGhost];
                std::copy( this->map().mapGlobalProcessToGlobalCluster().begin()+petsc_n_localWithoutGhost,
                           this->map().mapGlobalProcessToGlobalCluster().end(),
                           idx );
            }
            ierr = VecCreateGhostWithArray( this->comm(),
                                            petsc_n_localWithoutGhost, petsc_n_dof,
                                            petsc_n_localGhost, idx,
                                            thearray, &M_vec );
            CHKERRABORT( this->comm(),ierr );

            if ( petsc_n_localGhost > 0 )
                delete[] idx;

            this->M_is_initialized = true;
            this->setIsClosed( true );
        }

    /**
     * Constructor,  extracts a subvector from 'v' using mapping 'is'
     * without copy.
     */
    VectorPetsc( VectorPetsc<value_type> &v, IS &is )
        :
        super(),
        M_destroy_vec_on_exit( false )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        /* map */
        PetscInt n;
        ISGetSize(is,&n);
        datamap_ptrtype dm( new datamap_type(n, n, v.worldCommPtr()) );
        this->setMap(dm);
        /* init */
        VecGetSubVector(v.vec(), is, &this->M_vec);
        this->M_is_initialized = true;
        /* close */
        this->close(); /* no // assembly required */
#endif
    }

    VectorPetsc( Vector<value_type> const& v, std::vector<int> const& index )
        :
        super(),
        //super(v,index),
        M_destroy_vec_on_exit( false )
    {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsequenced"
#endif

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)

        VectorPetsc<T> const* V = dynamic_cast<VectorPetsc<T> const*> ( &v );
        int ierr=0;
        IS is;
        PetscInt *map;
        int n = index.size();
        PetscMalloc(n*sizeof(PetscInt),&map);
        for (int i=0; i<n; i++) map[i] = index[i];
        ierr = ISCreateGeneral(Environment::worldComm(),n,map,PETSC_COPY_VALUES,&is);
        CHKERRABORT( this->comm(),ierr );
        PetscFree(map);

        datamap_ptrtype dm( new datamap_type(n, n, v.worldCommPtr()) );
        this->setMap(dm);
        /* init */
        ierr = VecGetSubVector(V->vec(), is, &this->M_vec);
        CHKERRABORT( this->comm(),ierr );
        this->M_is_initialized = true;
        /* close */
        this->close(); /* no // assembly required */
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#endif
    }

    /**
     * Destructor, deallocates memory. Made virtual to allow
     * for derived classes to behave properly.
     */
    ~VectorPetsc () override
    {
        this->clear();
    }

    /**
     * Creates a copy of this vector and returns it in an \p
     * shared_ptr<>.  This must be overloaded in the derived classes.
     */
    clone_ptrtype clone () const override;


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
    void init ( const size_type N,
                const size_type n_local,
                const bool         fast=false ) override;

    /**
     * call init with n_local = N,
     */
    void init ( const size_type N,
                const bool         fast=false ) override
    {
        this->init( N,N,fast );
    }

    /**
     * call init with datamap,
     */
    void init ( datamap_ptrtype const& dm ) override;

    //@}

    /** @name Operator overloads
     */
    //@{

    value_type operator() ( const size_type i ) const override;
    value_type& operator() ( const size_type i ) override;

    /**
     *  \f$U = V\f$: copy all components.
     */
    Vector<value_type>& operator= ( const Vector<value_type> &V ) override;
    Vector<value_type>& operator= ( const VectorPetsc<value_type> &V );
    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    VectorPetsc<T> & operator += ( const Vector<value_type> &V ) override
    {
        this->add( 1., V );
        return *this;
    }

    /**
     * Subtraction operator.
     * Fast equivalent to \p U.add(-1, V).
     */
    VectorPetsc<T> & operator -= ( const Vector<value_type> &V ) override
    {
        this->add( -1., V );
        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{

    bool destroy_vec_on_exit() const
    {
        return M_destroy_vec_on_exit;
    }

    /**
     * @return dimension of the vector. This
     * function was formerly called \p n(), but
     * was renamed to get the \p PetscVector<T> class
     * closer to the C++ standard library's
     * \p std::vector container.
     */
    size_type size () const override
    {
        DCHECK( this->isInitialized() ) << "VectorPetsc not initialized";


        if ( !this->isInitialized() )
            return 0;

        int petsc_size=0;
        int ierr = VecGetSize( M_vec, &petsc_size );
        CHKERRABORT( this->comm(),ierr );
        return static_cast<size_type>( petsc_size );
    }
    /**
     * @return the local size of the vector
     * (index_stop-index_start)
     */
    size_type localSize() const override
    {
        DCHECK( this->isInitialized() ) << "VectorPetsc not initialized";

        int petsc_size=0;
        int ierr = VecGetLocalSize( M_vec, &petsc_size );
        CHKERRABORT( this->comm(),ierr );

        return static_cast<size_type>( petsc_size );
    }

    /**
     * Returns the raw PETSc vector context pointer.  Note this is generally
     * not required in user-level code. Just don't do anything crazy like
     * calling VecDestroy()!
     */
    Vec vec () const
    {
        FEELPP_ASSERT ( M_vec != 0 ).error( "invalid petsc vector" );
        return M_vec;
    }
    Vec& vec ()
    {
        FEELPP_ASSERT ( M_vec != 0 ).error( "invalid petsc vector" );
        return M_vec;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     *  \f$v = x*y\f$: coefficient-wise multiplication
     */
    void pointwiseMult ( Vector<T> const& x, Vector<T> const& y ) override;

    /**
     *  \f$v = x/y\f$: coefficient-wise divide
     */
    void pointwiseDivide ( Vector<T> const& x, Vector<T> const& y ) override;

    /**
     * Call the assemble functions
     */
    void close() override;

    /**
     * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
     * faster.
     */
    void zero () override;

    /**
     * Set entries to zero between \p start and \p stop
     */
    void zero ( size_type /*start*/,  size_type /*stop*/ ) override
    {
        this->zero();
    }

    /**
     * set the entries to the constant \p v
     */
    void setConstant( value_type v ) override
    {
        this->set( v );
    }

    /**
     * @returns the \p VectorPetsc<T> to a pristine state.
     */
    FEELPP_DONT_INLINE void clear () override;

    /**
     * Update ghost values
     */
    virtual void localize() {}

    /**
     *
     */
    void localize( const Vector<T>& V);

    /**
     * \f$ v(i) = \mathrm{value} \forall i\f$
     */
    virtual void set( const value_type& value );

    /**
     * v(i) = value
     */
    void set( const size_type i, const value_type& value ) override;

    /**
     * v([i1,i2,...,in]) = [value1,...,valuen]
     */
    void setVector( int* i, int n, value_type* v ) override;

    /**
     * v(i) += value
     */
    void add( const size_type i, const value_type& value ) override;

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    void addVector( int* i, int n, value_type* v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override;

    /**
     * \f$ U+=v \f$ where \p v is a std::vector<T>
     * and you
     * want to specify WHERE to add it
     */
    void addVector ( const std::vector<value_type>& v,
                     const std::vector<size_type>& dof_indices ) override
    {
        FEELPP_ASSERT ( v.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<v.size(); i++ )
            this->add ( dof_indices[i], v[i] );
    }

    /**
     * \f$ U+=V \f$ where U and V are type
     * \p NumericVector<T> and you
     * want to specify WHERE to add
     * the \p NumericVector<T> V
     */
    void addVector ( const Vector<value_type>& V,
                     const std::vector<size_type>& dof_indices ) override
    {
        FEELPP_ASSERT ( V.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<V.size(); i++ )
            this->add ( dof_indices[i], V( i ) );
    }


    /**
     * \f$ U+=A*V\f$, add the product of a \p MatrixSparse \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const Vector<value_type>& V_in,
                     const MatrixSparse<value_type>& A_in ) override;


    /**
     * \f$U+=V \f$ where U and V are type
     * uvlas::vector<T> and you
     * want to specify WHERE to add
     * the DenseVector<T> V
     */
    void addVector ( const ublas::vector<value_type>& V,
                     const std::vector<size_type>& dof_indices )
    {
        FEELPP_ASSERT ( V.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<V.size(); i++ )
            this->add ( dof_indices[i], V( i ) );
    }

    /**
     * \f$ U=v \f$ where v is a DenseVector<T>
     * and you want to specify WHERE to insert it
     */
    void insert ( const std::vector<T>& /*v*/,
                  const std::vector<size_type>& /*dof_indices*/ ) override
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * \f$U=V\f$, where U and V are type
     * Vector<T> and you
     * want to specify WHERE to insert
     * the Vector<T> V
     */
    void insert ( const Vector<T>& V,
                  const std::vector<size_type>& dof_indices ) override;


    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to insert
     * the DenseVector<T> V
     */
    void insert ( const ublas::vector<T>& V,
                  const std::vector<size_type>& dof_indices ) override;

    /**
     * Scale each element of the
     * vector by the given factor.
     */
    void scale ( const T factor ) override;


    /**
     * \f$ U(0-DIM)+=s\f$.
     * Addition of \p s to all components. Note
     * that \p s is a scalar and not a vector.
     */
    void add ( const value_type& v_in ) override;

    /**
     * \f$ U+=V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add ( const Vector<value_type>& v ) override;

    /**
     * \f$ U+=a*V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add ( const value_type& a_in, const Vector<value_type>& v_in ) override;
    /**
     * Replaces each component of a vector by its reciprocal.
     */
    int reciprocal() override;

    /**
     * @return the minimum element in the vector.
     * In case of complex numbers, this returns the minimum
     * Real part.
     */
    real_type min () const override;

    /**
     * @returns the maximum element in the vector.
     * In case of complex numbers, this returns the maximum
     * Real part.
     */
    real_type max() const override;
    real_type maxWithIndex( int* index=nullptr ) const override;


    //! Replaces every element in a vector with its absolute value
    void abs() override;

    /**
     * @return the \f$l_1\f$-norm of the vector, i.e.
     * the sum of the absolute values.
     */
    real_type l1Norm () const override;

    /**
     * @returns the \f$l_2\f$-norm of the vector, i.e.
     * the square root of the sum of the
     * squares of the elements.
     */
    real_type l2Norm () const override;

    /**
     * @return the sum of the vector components, i.e.
     * the sum of the values.
     */
    value_type sum () const override;

    /**
     * @returns the maximum absolute value of the
     * elements of this vector, which is the
     * \f$l_\infty\f$-norm of a vector.
     */
    real_type linftyNorm () const override;

    /**
     * @return the index of the first vector element
     * actually stored on this processor
     */
    size_type firstLocalIndex () const override
    {
        assert ( this->isInitialized() );

        int ierr=0, petsc_first=0, petsc_last=0;

        ierr = VecGetOwnershipRange ( M_vec, &petsc_first, &petsc_last );
        CHKERRABORT( this->comm(),ierr );

        return static_cast<size_type>( petsc_first );
    }



    /**
     * @return the index of the last vector element
     * actually stored on this processor
     */
    size_type lastLocalIndex () const override
    {
        assert ( this->isInitialized() );

        int ierr=0, petsc_first=0, petsc_last=0;

        ierr = VecGetOwnershipRange ( M_vec, &petsc_first, &petsc_last );
        CHKERRABORT( this->comm(),ierr );

        return static_cast<size_type>( petsc_last );
    }


    /**
     * Print the contents of the vector in Matlab's format. Optionally
     *  prints the vector to the file named \p name.  If \p name is
     *  not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL", bool renumber = false ) const override;

    value_type dot( Vector<T> const& __v ) const override;

    /**
     * This function creates a vector which is defined
     * by the row indices given in the "rows" entries.
     */
    std::shared_ptr<Vector<T> >
    createSubVector( std::vector<size_type> const& rows,
                     bool checkAndFixRange=true ) const override;

    /**
     * Copy (default) or add (boolean init=false) entries of subvector (already built from a createSubVector)
     * into row indices given in the "rows" entries.
     */
    void
    updateSubVector( std::shared_ptr<Vector<T> > & subvector,
                     std::vector<size_type> const& rows,
                     bool init=true ) override;

    void save( std::string const& filename="default_archive_name", std::string const& format="binary" ) override;
    void load( std::string const& filename="default_archive_name", std::string const& format="binary" ) override;

    //@}


private:
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
    {
        CHECK( this->isInitialized() ) << "Vector not initialized";
        if ( !this->closed() )
            const_cast<VectorPetsc<T>*>( this )->close();

        int ierr = 0;

        double * array;
        ierr = VecGetArray(M_vec, &array);
        CHKERRABORT( this->comm(),ierr );

        size_type n = this->map().nLocalDofWithoutGhost();
        ar << BOOST_SERIALIZATION_NVP( n );
        ar << boost::serialization::make_array( array, n );

        ierr = VecRestoreArray(M_vec, &array );
        CHKERRABORT( this->comm(),ierr );
    }

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
    {
        CHECK( this->isInitialized() ) << "Vector not initialized";
        if ( !this->closed() )
            this->close();

        int ierr = 0;

        double * array;
        ierr = VecGetArray(M_vec, &array);
        CHKERRABORT( this->comm(),ierr );

        size_type n = invalid_v<size_type>;
        ar >> BOOST_SERIALIZATION_NVP( n );
        CHECK( n == this->map().nLocalDofWithoutGhost() ) << "vector is not compatible with array serialized";
        ar >>  boost::serialization::make_array( array, n );

        ierr = VecRestoreArray(M_vec, &array );
        CHKERRABORT( this->comm(),ierr );

        if ( this->comm().size()>1 )
        {
            VectorPetscMPI<T>* vecOutPetscMPI = dynamic_cast<VectorPetscMPI<T>*>( &(*this) );
            if ( vecOutPetscMPI )
                vecOutPetscMPI->localize();
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected:

public:

    // disable
    VectorPetsc( VectorPetsc const & v )
        :
        super( v ),
        M_destroy_vec_on_exit( true )
    {
        //FEELPP_ASSERT( v.closed() ).error( "copied vector is not closed" );
        if ( !v.closed() )
            const_cast<VectorPetsc<T>*>( &v )->close();

        VecDuplicate( v.M_vec, &M_vec );
        VecCopy( v.M_vec, M_vec );
        this->M_is_initialized = true;
        //this->close();
        this->setIsClosed( true );
    }

    void getSubVectorPetsc( std::vector<size_type> const& rows,
                            Vec &subvec,
                            bool init=true ) const;

private :
    void pointwiseOperationsImpl( Vector<T> const& x, Vector<T> const& y, int op );

protected:

    /**
     * Petsc vector datatype to store values
     */
    Vec M_vec;

    /**
     * This boolean value should only be set to false
     * for the constructor which takes a PETSc Vec object.
     */
    bool M_destroy_vec_on_exit;
};


//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------//

template<typename T>
class VectorPetscMPI : public VectorPetsc<T>
{
    typedef VectorPetsc<T> super;
public:
    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;
    
    friend class boost::serialization::access;

    typedef typename super::value_type value_type;
    using size_type = typename datamap_type::size_type;
    VectorPetscMPI()
        :
        super()
    {}

    VectorPetscMPI( Vec v, datamap_ptrtype const& dm, bool duplicate = false );

    VectorPetscMPI( datamap_ptrtype const& dm, bool doInit=true );

    template<typename Storage>
    VectorPetscMPI( VectorUblas<T,Storage> const& v )
        :
        super( v )
        {}

    ~VectorPetscMPI() override
    {
        this->clear();
    }


    FEELPP_DEPRECATED
    void init( const size_type N,
               const size_type n_local,
               const bool fast=false ) override;
    /**
     * call init with datamap,
     */
    void init( datamap_ptrtype const& dm ) override;

    value_type operator() ( const size_type i ) const override;
    value_type& operator() ( const size_type i ) override;

    /**
     *  \f$U = V\f$: copy all components.
     */
    Vector<value_type>& operator= ( const Vector<value_type> &V ) override;
    Vector<value_type>& operator= ( const VectorPetscMPI<value_type> &V );

    /**
     * \f$ v(i) = \mathrm{value} \forall i\f$
     */
    void set( const value_type& value ) override;

    /**
     * \f$ U(0-DIM)+=s\f$.
     * Addition of \p s to all components. Note
     * that \p s is a scalar and not a vector.
     */
    void add( const value_type& v_in ) override;

    /**
     * \f$ U+=a*V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add( const value_type& a_in, const Vector<value_type>& v_in ) override;

    /**
     * v(i) = value (i is global process index)
     */
    void set( size_type i, const value_type& value ) override;

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen] (i1,i2,... is global process index)
     */
    void setVector( int* i, int n, value_type* v ) override;

    /**
     * v(i) += value (i is global process index)
     */
    void add( const size_type i, const value_type& value ) override;

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen] (i1,i2,... is global process index)
     */
    void addVector( int* i, int n, value_type* v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override;

    /**
     *  \f$v = x*y\f$: coefficient-wise multiplication
     */
    void pointwiseMult( Vector<T> const& x, Vector<T> const& y ) override;

    /**
     *  \f$v = x/y\f$: coefficient-wise divide
     */
    void pointwiseDivide( Vector<T> const& x, Vector<T> const& y ) override;

    /**
     * Set all entries to zero. Equivalent to \p v = 0.
     */
    void zero() override;
    void zero( size_type /*start*/,  size_type /*stop*/ ) override
    {
        this->zero();
    }

    /**
     * Replaces each component of a vector by its reciprocal.
     */
    int reciprocal() override;

    /**
     * @returns the \p VectorPetsc<T> to a pristine state.
     */
    void clear() override;

    /**
     * Update ghost values
     */
    void localize() override;

    /**
     * Call the assemble functions and update ghost values
     */
    void close() override;

    size_type firstLocalIndex() const override;
    size_type lastLocalIndex() const override;
    size_type localSize() const override;

    void duplicateFromOtherPartition( Vector<T> const& vecInput );

private :

    void initImpl( const bool fast = false );

    void pointwiseOperationOthersPetscImpl( Vector<T> const& x, Vector<T> const& y, int op );

    void duplicateFromOtherPartition_run( Vector<T> const& vecInput );

};

template<typename T>
class VectorPetscMPIRange : public VectorPetscMPI<T>
{
    typedef VectorPetscMPI<T> super_type;
public:
    typedef typename super_type::value_type value_type;
    using datamap_ptrtype = typename super_type::datamap_ptrtype;
    using size_type = typename super_type::size_type;
    VectorPetscMPIRange( datamap_ptrtype const& dm );

    VectorPetscMPIRange( Vec v, datamap_ptrtype const& dm, bool duplicate = false );

    VectorPetscMPIRange( Vec v, Vec vGhost, datamap_ptrtype const& dm, bool duplicate = false );

    VectorPetscMPIRange( Vec v, Vec vGhost, VecScatter vecScatterGhost, datamap_ptrtype const& dm );

    template<typename Storage>
    VectorPetscMPIRange( VectorUblas<T,Storage> const& v )
        :
        super_type( v.mapPtr(), false ),
        M_destroyVecGhostOnExit( true ),
        M_destroyVecScatterGhostOnExit( true )
        {
            const PetscScalar* arrayActive = ( this->map().nLocalDofWithoutGhost() > 0 )? std::addressof( *v.begin() ) : NULL;
            const PetscScalar* arrayGhost = ( this->map().nLocalGhosts() > 0 )? std::addressof( *v.beginGhost() ) : NULL;
            this->initRangeView( arrayActive,arrayGhost );
        }
    ~VectorPetscMPIRange() override
    {
        this->clear();
    }

    /**
     * call init with datamap,
     */
    void init( datamap_ptrtype const& dm ) override;

    /**
     * @returns the \p VectorPetsc<T> to a pristine state.
     */
    void clear() override;

    value_type operator() ( const size_type i ) const override;
    value_type& operator() ( const size_type i ) override;

    /**
     *  \f$U = V\f$: copy all components.
     */
    Vector<value_type>& operator= ( const Vector<value_type> &V ) override;
    Vector<value_type>& operator= ( const VectorPetscMPIRange<value_type> &V );

    /**
     * \f$ v(i) = \mathrm{value} \forall i\f$
     */
    void set( const value_type& value ) override;

    /**
     * \f$ U(0-DIM)+=s\f$.
     * Addition of \p s to all components. Note
     * that \p s is a scalar and not a vector.
     */
    void add( const value_type& v_in ) override;

    /**
     * \f$ U+=a*V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add( const value_type& a_in, const Vector<value_type>& v_in ) override;

    /**
     *  \f$v = x*y\f$: coefficient-wise multiplication
     */
    void pointwiseMult( Vector<T> const& x, Vector<T> const& y ) override;

    /**
     *  \f$v = x/y\f$: coefficient-wise divide
     */
    void pointwiseDivide( Vector<T> const& x, Vector<T> const& y ) override;

    /**
     * Set all entries to zero. Equivalent to \p v = 0.
     */
    void zero() override;

    /**
     * Replaces each component of a vector by its reciprocal.
     */
    int reciprocal() override;

    /**
     * Update ghost values
     */
    void localize() override;

    /**
     * Returns the raw PETSc vector of ghosts in context pointer
     */
    Vec vecGhost() const
    {
        FEELPP_ASSERT ( M_vecGhost != 0 ).error( "invalid petsc vector" );
        return M_vecGhost;
    }
    Vec& vecGhost()
    {
        FEELPP_ASSERT ( M_vecGhost != 0 ).error( "invalid petsc vector" );
        return M_vecGhost;
    }

    VecScatter vecScatterGhost() const
    {
        FEELPP_ASSERT ( M_vecScatterGhost != 0 ).error( "invalid petsc vector" );
        return M_vecScatterGhost;
    }
    VecScatter& vecScatterGhost()
    {
        FEELPP_ASSERT ( M_vecScatterGhost != 0 ).error( "invalid vector scatter" );
        return M_vecScatterGhost;
    }

private :
    void initRangeView( const PetscScalar arrayActive[], const PetscScalar arrayGhost[] );
    void initVecScatterGhost();
private :
    Vec M_vecGhost;
    VecScatter M_vecScatterGhost;
    bool M_destroyVecGhostOnExit, M_destroyVecScatterGhostOnExit;


};
/**
 * @addtogroup FreeFunctions
 * @{
 */

/**
 * @brief create a shared pointer VectorPetsc from PETSc Vec \p v
 */
template<typename SizeT>
#if BOOST_VERSION < 105900
vector_ptrtype vec( Vec v, datamap_ptrtype<SizeT> d );
#else
vector_uptrtype vec( Vec v, datamap_ptrtype<SizeT> d );
#endif
/**
 * @}
 */


} // Feel


#endif /* FEELPP_HAS_PETSC */
#endif /* __VectorPetsc_H */
