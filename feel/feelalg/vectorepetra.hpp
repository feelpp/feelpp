/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s):
   Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Klaus Sapelza <klaus.sapelza@epfl.ch>
   Date: 2006-09-14

   Copyright (C) 2006,2007 EPFL
   Copyright (C) 2006,2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file vectorepetra.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Klaus.Sapelza <klaus.sapelza@epfl.ch>
   \date 2006-09-14
*/
#ifndef __VectorEpetra_H
#define __VectorEpetra_H 1

#include <feel/feelconfig.h>

#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelcore/application.hpp>
#if defined(FEELPP_HAS_TRILINOS_EPETRA)

#if defined(FEELPP_HAS_MPI)
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif /* FEELPP_HAS_MPI */

#include <EpetraExt_MultiVectorOut.h>
#include <Epetra_FEVector.h>
#include <Epetra_Vector.h>



namespace Feel
{
/**
 * \class VectorEpetra
 * \brief Wrapper for epetra vectors
 *
 * Epetra vector. Provides a nice interface to the Epetra data
 * structures for parallel vectors.
 *
 * @author Christophe Prud'homme
 * @author Klaus Sapelza
 * @see Trilinos, Vector
 */
template<typename T>
class VectorEpetra : public Vector<T>
{
    typedef Vector<T> super;


public:


    /** @name Typedefs
     */
    //@{

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    typedef typename super::clone_ptrtype clone_ptrtype;

    typedef VectorEpetra<value_type> epetra_vector_type;
    typedef boost::shared_ptr<epetra_vector_type> epetra_vector_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    VectorEpetra ();

    /**
     * Constructor. Set dimension to \p n and initialize all elements with zero.
     */
    VectorEpetra ( Epetra_BlockMap const& emap );

    /**
     * Constructor.  Creates a VectorEpetra assuming you already have a
     * valid Epetra_FEVector object.
     */
    VectorEpetra ( Epetra_FEVector const * v );

    /**
     * Constructor.  Creates a VectorEpetra assuming you already have a
     * valid Epetra_Vector object.
     */
    VectorEpetra ( Epetra_Vector const * v );

    /**
     * Constructor.  Creates a VectorEpetra assuming you already have a
     * valid Epetra_Vector object.
     */
    VectorEpetra ( VectorEpetra const& v );

    /**
     * Destructor, deallocates memory. Made virtual to allow
     * for derived classes to behave properly.
     */
    ~VectorEpetra ();

    /**
     * Returns the Epetra map
     */
    Epetra_BlockMap Map() const
    {
        return M_vec.Map();
    }


    /**
     * Creates a copy of this vector and returns it in an \p
     * shared_ptr<>.  This must be overloaded in the derived classes.
     */
    clone_ptrtype clone () const
    {
        clone_ptrtype cloned_vector ( new VectorEpetra<T> );

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
     * On \p fast = false the vector is filled by zeros.
     */

    void init ( Epetra_BlockMap const& emap, const bool fast = false  );

    void init ( const size_type N, const size_type n_local, const bool fast = false );

    /** @name Operator overloads
     */
    //@{
    //operator Epetra_Vector() { return M_vec; }



    value_type operator() ( const size_type i ) const
    {
        checkInvariants();
        FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
        FEELPP_ASSERT ( ( ( i >= this->firstLocalIndex() ) &&
                          ( i <  this->lastLocalIndex() ) ) )( i )( this->firstLocalIndex() )( this->lastLocalIndex() ).warn( "invalid vector index" );

        value_type value=0.;
        int ierr=0, dummy;
        double* values;

        ierr = M_vec.ExtractView( &values, &dummy );

        value = values[i - this->firstLocalIndex()];
        return static_cast<value_type>( value );

    }

    value_type& operator() ( const size_type i )
    {
        checkInvariants();
        FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
        FEELPP_ASSERT ( ( ( i >= this->firstLocalIndex() ) &&
                          ( i <  this->lastLocalIndex() ) ) )( i )( this->firstLocalIndex() )( this->lastLocalIndex() ).warn( "invalid vector index" );

        return M_vec[0][ ( int )i-this->firstLocalIndex() ];

    }


    VectorEpetra& operator=( VectorEpetra const& v )
    {
        if ( &v != this )
        {
            super::operator=( v );
            M_emap = v.Map();
            M_vec = v.vec();

            this->M_is_initialized = true;
        }

        checkInvariants();
        return *this;
    }

    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    Vector<T> & operator += ( const Vector<T> &V )
    {

        this->add( 1., V );

        return *this;
    }

    /**
     * Subtraction operator.
     * Fast equivalent to \p U.add(-1, V).
     */
    super & operator -= ( const super &V )
    {

        this->add( -1., V );

        return *this;
    }


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @return dimension of the vector. This
     * function was formerly called \p n(), but
     * was renamed to get the \p VectorEpetra class
     * closer to the C++ standard library's
     * \p std::vector container.
     */
    size_type size () const
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "VectorEpetra not initialized" );


        if ( !this->isInitialized() )
            return 0;

        int epetra_size=0;
        epetra_size = M_vec.GlobalLength();
        return static_cast<size_type>( epetra_size );

        return 0;
    }
    /**
     * @return the local size of the vector
     * (index_stop-index_start)
     */
    size_type localSize() const
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "VectorEpetra not initialized" );

        int epetra_size=0;
        epetra_size = M_vec.MyLength();
        DVLOG(2) << "[VectorEpetra::localSize] localSize= " << epetra_size  << "\n";
        return static_cast<size_type>( epetra_size );
    }

    /**
     * Returns the raw Epetra vector context pointer.  Note this is generally
     * not required in user-level code. Just don't do anything crazy like
     * calling VecDestroy()!
     */
    Epetra_FEVector const& vec () const
    {
        //FEELPP_ASSERT (M_vec != 0).error( "invalid epetra vector" );
        return M_vec;
    }

    /**
     * Returns the raw Epetra vector context pointer.  Note this is generally
     * not required in user-level code. Just don't do anything crazy like
     * calling VecDestroy()!
     */
    Epetra_FEVector& vec ()
    {
        //FEELPP_ASSERT (M_vec != 0).error( "invalid epetra vector" );
        return M_vec;
    }

    /**
     * @author Vielfaure Florent
     *
     * Returns the raw Epetra_Vector shared pointer.
     * The returned object is a view of the Epetra_FEVector contained
     * in the vectorEpetra
     */
    boost::shared_ptr<Epetra_Vector> epetraVector ()
    {
        double** V;
        M_vec.ExtractView( &V );
        boost::shared_ptr<Epetra_Vector> EV( new Epetra_Vector( View,M_vec.Map(),V[0] ) );
        return EV;
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
     * Call the assemble functions (necessary in PetSc, but not in Trilinos)
     */
    void close ()
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "VectorEpetra<> not initialized" );

        int ierr=0;
        ierr = M_vec.GlobalAssemble( Add );

        this->M_is_closed = true;
    }

    /**
     * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
     * faster.
     */
    void zero ()
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "VectorEpetra<> not initialized" );

        M_vec.PutScalar( 0.0 );

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
    void setConstant( value_type v )
    {
        this->set( v );
    }

    /**
     * @returns the \p VectorEpetra to a pristine state.
     */
    void clear ()
    {
        if ( this->isInitialized() ) //&& (this->M_destroy_vec_on_exit))
        {
            M_emap = Epetra_BlockMap( -1, 0, 0, M_vec.Comm() );
            M_vec.ReplaceMap( M_emap );
            M_vec.PutScalar( 0.0 );
        }

        this->M_is_closed = this->M_is_initialized = false;
    }

    /**
     * \f$ v(i) = \mathrm{value} \forall i\f$
     */
    void set ( const value_type& value );

    /**
     * v(i) = value
     */
    void set ( size_type i, const value_type& value );

    /**
     * v(i) += value
     */
    void add ( size_type i, const value_type& value );

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    void addVector ( int* i, int n, value_type* v );


    void addVector ( VectorEpetra& v )
    {
        M_vec.Update( 1.0, v.vec(), 1.0 );
    }


    /**
     * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const Vector<T>& _v,
                     const MatrixSparse<T>& _M );

    /**
     * \f$ U+=v \f$ where \p v is a std::vector<T>
     * and you want to specify WHERE to add it
     */
    void addVector ( const std::vector<value_type>& v,
                     const std::vector<size_type>& dof_indices )
    {
        FEELPP_ASSERT ( v.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<v.size(); i++ )
            this->add ( dof_indices[i], v[i] );
    }


    /**
     * \f$ U+=V \f$ where U and V are type
     * \p Vector<T> and you
     * want to specify WHERE to add
     * the \p Vector<T> V
     */
    void addVector ( const Vector<value_type>& V,
                     const std::vector<size_type>& dof_indices )
    {
        FEELPP_ASSERT ( V.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<V.size(); i++ )
            this->add ( dof_indices[i], V( i ) );
    }
    //
    //
    //     /**
    //      * \f$ U+=A*V\f$, add the product of a \p MatrixSparse \p A
    //      * and a \p Vector \p V to \p this, where \p this=U.
    //      */
    //     void addVector (const Vector<value_type>& V_in,
    //                     const MatrixSparse<value_type>& A_in)
    //     {
    //
    //
    //     }
    //

    /**
     * \f$U+=V \f$ where U and V are type
     * ublas::vector<T> and you
     * want to specify WHERE to add
     * the ublas::vector<T> V
     */
    void addVector ( const ublas::vector<value_type>& V,
                     const std::vector<size_type>& dof_indices )
    {
        FEELPP_ASSERT ( V.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<V.size(); i++ )
            this->add ( dof_indices[i], V( i ) );
    }

    /**
     * \f$ U(0-DIM)+=s\f$.
     * Addition of \p s to all components. Note
     * that \p s is a scalar and not a vector.
     */
    void add ( const value_type& v_in )
    {
        value_type v = static_cast<value_type>( v_in );
        const int n   = static_cast<int>( this->size() );

        for ( int i=0 ; i<n ; i++ )
        {
            M_vec.SumIntoGlobalValues( 1,&i,&v );  //indices are in global index space
            //M_vec.SumIntoMyValues(1,&v,&i);      //indices are in local index space
        }

    }

    /**
     * \f$ U+=V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add ( const Vector<value_type>& v )
    {
        VectorEpetra* v_ptr = const_cast<VectorEpetra*>( dynamic_cast< VectorEpetra const*>( &v ) );
        this->add ( 1., *v_ptr );
    }

    /**
     * \f$ U+=a*V \f$ .
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add ( const value_type& a_in, const Vector<value_type>& v_in )
    {

        const value_type a = static_cast<value_type>( a_in );
        const VectorEpetra<T>* v = dynamic_cast<const VectorEpetra<T>*>( &v_in );


        assert ( v != NULL );
        assert( this->size() == v->size() );

        M_vec.Update( a, v->M_vec,1. );

    }

    /**
     * \f$ U=v \f$ where v is a DenseVector<T>
     * and you want to specify WHERE to insert it
     */
    void insert ( const std::vector<T>& v,
                  const std::vector<size_type>& dof_indices )
    {
        //to be implemented
    }

    /**
     * \f$U=V\f$, where U and V are type
     * Vector<T> and you
     * want to specify WHERE to insert
     * the Vector<T> V
     */
    void insert ( const Vector<T>& V,
                  const std::vector<size_type>& dof_indices )
    {
        //to be implemented
    }

    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to insert
     * the DenseVector<T> V
     */
    void insert ( const ublas::vector<T>& /*V*/,
                  const std::vector<size_type>& /*dof_indices*/ )
    {
        //to be implemented
    }

    /**
     * Scale each element of the
     * vector by the given factor.
     */
    void scale ( const T /*factor*/ )
    {
        //to be implemented
    }


    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( std::vector<T>& /*v_local*/ ) const
    {
        //to be implemented
    }

    /**
     * Same, but fills a \p Vector<T> instead of
     * a \p std::vector.
     */
    void localize ( Vector<T>& /*v_local*/ ) const
    {
        //to be implemented
    }

    /**
     * Creates a local vector \p v_local containing
     * only information relevant to this processor, as
     * defined by the \p send_list.
     */
    void localize ( Vector<T>& /*v_local*/,
                    const std::vector<size_type>& /*send_list*/ ) const
    {
    }

    /**
     * Updates a local vector with selected values from neighboring
     * processors, as defined by \p send_list.
     */
    void localize ( const size_type /*first_local_idx*/,
                    const size_type /*last_local_idx*/,
                    const std::vector<size_type>& /*send_list*/ )
    {
        //to be implemented
    }
    /**
     * Creates a local copy of the global vector in
     * \p v_local only on processor \p proc_id.  By
     * default the data is sent to processor 0.  This method
     * is useful for outputting data from one processor.
     */
    void localizeToOneProcessor ( std::vector<T>& /*v_local*/,
                                  const size_type /*proc_id*/=0 ) const
    {
        //to be implemented
    }

    /**
     * @return the minimum element in the vector.
     * In case of complex numbers, this returns the minimum
     * Real part.
     */
    real_type min () const
    {
        assert ( this->isInitialized() );

        double min=0.;

        M_vec.MinValue( &min );

        // this return value is correct: VecMin returns a PetscReal
        return static_cast<double>( min );
    }

    /**
     * @returns the maximum element in the vector.
     * In case of complex numbers, this returns the maximum
     * Real part.
     */
    real_type max() const
    {
        assert ( this->isInitialized() );

        double max=0.;

        M_vec.MaxValue( &max );

        // this return value is correct: VecMin returns a PetscReal
        return static_cast<double>( max );
    }

    /**
     * @return the \f$l_1\f$-norm of the vector, i.e.
     * the sum of the absolute values.
     */
    real_type l1Norm () const
    {
        //assert(this->closed());

        double value=0.;

        M_vec.Norm1( &value );

        return static_cast<Real>( value );
    }

    /**
     * @returns the \f$l_2\f$-norm of the vector, i.e.
     * the square root of the sum of the
     * squares of the elements.
     */
    real_type l2Norm () const
    {
        //assert(this->closed());

        double value=0.;

        M_vec.Norm2( &value );

        return static_cast<Real>( value );

    }

    /**
     * @returns the maximum absolute value of the
     * elements of this vector, which is the
     * \f$l_\infty\f$-norm of a vector.
     */
    real_type linftyNorm () const
    {
        assert( this->closed() );

        double value=0.;

        M_vec.NormInf( &value );

        return static_cast<Real>( value );

    }

    /**
     * @return the sum of the vector.
     */
    value_type sum() const
    {
        //assert(this->closed());

        double value=0.;
        double global_sum=0;

        double const * pointers( M_vec[0] );

        for ( int i( 0 ); i < M_vec.MyLength(); ++i , ++pointers )
            value += *pointers;

        M_vec.Comm().SumAll( &value, &global_sum, 1 );

        return static_cast<Real>( global_sum );


    }


    /**
     * @return the global index of the first vector element
     * actually stored on this processor
     *
     * \return the minimum global Index owned by this processor
     */
    size_type firstLocalIndex () const
    {

        int epetra_first = 0;
        assert ( this->isInitialized() );
        //epetra_first = M_vec.Map().MinMyGID();
        epetra_first = M_vec.Map().MinLID();
        DVLOG(2) << "[VectorEpetra::firstLocalIndex] firstLocalIndex= " << epetra_first  << "\n";
        return static_cast<size_type>( epetra_first );
    }



    /**
     * @return the index of the last vector element
     * actually stored on this processor
     *
     * \return the maximum global Index owned by this processor + 1
     */
    size_type lastLocalIndex () const
    {
        int epetra_last = 0;
        assert ( this->isInitialized() );
        //epetra_last = M_vec.Map().MaxMyGID();
        epetra_last = M_vec.Map().MaxLID();
        DVLOG(2) << "[VectorEpetra::lastLocalIndex] lastLocalIndex= " << epetra_last+1  << "\n";
        return static_cast<size_type>( epetra_last )+1;
    }

    /**
     * print Epetra vector in a matlab file \p name
     * \param name filename of the matlab file
     * \sa MatrixEpetra::printMatlab
     */
    void printMatlab ( const std::string name = "", bool renumber = false ) const;

    //   @}



protected:

private:
    void checkInvariants() const;
private:
    Epetra_BlockMap M_emap;
    Epetra_FEVector M_vec;

    /**
     * This boolean value should only be set to false
     * for the constructor which takes a Epetra Vec object.
     */
    //const bool M_destroy_vec_on_exit;
};

template<typename T>
DebugStream&
operator<<( DebugStream& __os, VectorEpetra<T> const& __n );


template<typename T>
NdebugStream&
operator<<( NdebugStream& __os, VectorEpetra<T> const& __n );

DebugStream&
operator<<( DebugStream& __os, Epetra_BlockMap const& __n );


NdebugStream&
operator<<( NdebugStream& __os, Epetra_BlockMap const& __n );


} // Feel
#endif /* FEELPP_HAS_EPETRA */
#endif /* __VectorEpetra_H */
