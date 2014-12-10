/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-06-20

  Copyright (C) 2012 Universit� Joseph Fourier (Grenoble I)

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
   \file vectoreigen.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-06-20
 */
#ifndef __VectorEigen_H
#define __VectorEigen_H 1

#include <boost/serialization/complex.hpp>
#include <set>
#include <boost/operators.hpp>
#include <Eigen/Core>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/matrixeigendense.hpp>


namespace Feel
{
//template<typename T> class MatrixEigenDense;

template<typename T>
struct get_real_type
{};

template<>
struct get_real_type<double>
{
    typedef double value_type;
};

template<>
struct get_real_type<std::complex<double> >
{
    typedef double value_type;
};

/*!
 * \class VectorEigen
 * \brief interface to vector
 *
 * \code
 * VectorEigen<T> m;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T >
class VectorEigen
    : public Vector<T>
    , boost::addable<VectorEigen<T> >
    , boost::subtractable<VectorEigen<T> >

{
    typedef Vector<T> super1;
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename get_real_type<T>::value_type real_type;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vector_type;
    typedef VectorEigen<value_type> this_type;
    typedef typename super1::clone_ptrtype clone_ptrtype;

    typedef typename super1::datamap_type datamap_type;
    typedef typename super1::datamap_ptrtype datamap_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    VectorEigen();

    VectorEigen( size_type __s, WorldComm const& _worldComm = Environment::worldComm() );

    VectorEigen( datamap_ptrtype const& dm );

    VectorEigen( size_type __s, size_type __n_local, WorldComm const& _worldComm = Environment::worldComm() );

    VectorEigen( VectorEigen const & m );

    ~VectorEigen();

    /**
     * Creates a copy of this vector and returns it in an \p shared_ptr<>.
     * This must be overloaded in the derived classes.
     */
    clone_ptrtype clone () const ;

    /**
     * Change the dimension of the vector to \p N. The reserved memory
     * for this vector remains unchanged if possible, to make things
     * faster, but this may waste some memory, so take this in the
     * back of your head.  However, if \p N==0 all memory is freed,
     * i.e. if you want to resize the vector and release the memory
     * not needed, you have to first call \p init(0) and then \p
     * init(N). This cited behaviour is analogous to that of the STL
     * containers.
     *
     * On \p fast==false, the vector is filled by zeros.
     */
    void init ( const size_type N,
                const size_type n_local,
                const bool      fast=false );

    /**
     * call init with n_local = N,
     */
    void init ( const size_type n,
                const bool      fast=false );

    /**
     * init from a \p DataMap
     */
    void init( datamap_ptrtype const& dm );


    //@}
    /** @name Operator overloads
     */
    //@{

    /**
     *  \f$U = V\f$: copy all components.
     */
    Vector<value_type>& operator= ( const Vector<value_type> &V );

    /**
     * Access components, returns \p u(i).
     */
    T operator()( size_type i ) const
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
        FEELPP_ASSERT ( ( i >= this->firstLocalIndex() ) &&
                        ( i < this->lastLocalIndex() ) )
        ( i )
        ( this->firstLocalIndex() )
        ( this->lastLocalIndex() ).error( "vector invalid index" );

        return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Access components, returns \p u(i).
     */
    T& operator()( size_type i )
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
        FEELPP_ASSERT ( ( i >= this->firstLocalIndex() ) &&
                        ( i <  this->lastLocalIndex() ) )
        ( i )
        ( this->firstLocalIndex() )
        ( this->lastLocalIndex() ).error( "vector invalid index" );
        return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Access components, returns \p u(i).
     */
    T operator[]( size_type i ) const
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
        FEELPP_ASSERT ( ( i >= this->firstLocalIndex() ) &&
                        ( i < this->lastLocalIndex() ) )
        ( i )
        ( this->firstLocalIndex() )
        ( this->lastLocalIndex() ).error( "vector invalid index" );

        return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Access components, returns \p u(i).
     */
    T& operator[]( size_type i )
    {
        FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
        FEELPP_ASSERT ( ( i >= this->firstLocalIndex() ) &&
                        ( i <=  this->lastLocalIndex() ) )
        ( i )
        ( this->firstLocalIndex() )
        ( this->lastLocalIndex() ).error( "vector invalid index" );
        return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    Vector<T>& operator+=( const Vector<T>& v )
    {
        checkInvariant();
        add( 1., v );
        return *this;
    }

    /**
     * Subtraction operator.
     * Fast equivalent to \p U.add(-1, V).
     */
    Vector<T>& operator-=( const Vector<T>& v )
    {
        checkInvariant();
        add( -1., v );

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{


    /**
     * if the vector is a range, return the first index of the range,
     * otherwise returns 0
     */
    size_type start() const;

    /**
     * return row_start, the index of the first
     * vector row stored on this processor
     */
    unsigned int rowStart () const
    {
        checkInvariant();
        return 0;
    }

    /**
     * return row_stop, the index of the last
     * vector row (+1) stored on this processor
     */
    size_type rowStop () const
    {
        checkInvariant();
        return 0;
    }

    /**
     * \return true if vector is initialized/usable, false otherwise
     */
    bool isInitialized() const
    {
        return true;
    }

    /**
     * \c close the eigen vector, that will copy the content of write
     * optimized vector into a read optimized vector
     */
    void close () const;


    /**
     * see if vector has been closed
     * and fully assembled yet
     */
    bool closed() const
    {
        return true;
    }


    /**
     * Returns the read optimized eigen vector.
     */
    vector_type const& vec () const
    {
        return M_vec;
    }

    /**
     * Returns the read optimized eigen vector.
     */
    vector_type & vec ()
    {
        return M_vec;
    }


    //@


    /** @name  Mutators
     */
    //@{


    /**
     * set the entries to the constant \p v
     */
    void setConstant( value_type v )
        {
            M_vec.setConstant( v );
        }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * resize with size @p n
     */
    void resize( size_type n, bool preserve = true );

    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear ();

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ()
    {
        M_vec.setZero( M_vec.size() );
    }

    void zero ( size_type /*start1*/, size_type /*stop1*/ )
    {
        //eigen::project( (*this), eigen::range( start1, stop1 ) ) = eigen::zero_vector<value_type>( stop1 );
    }

    /**
     * Add \p value to the value already accumulated
     */
    void add ( const size_type i, const value_type& value )
    {
        checkInvariant();
        M_vec( i ) += value;
    }

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    void addVector ( int* i, int n, value_type* v )
    {
        for ( int j = 0; j < n; ++j )
            M_vec( i[j] ) += v[j];
    }

    /**
     * set to \p value
     */
    void set ( size_type i, const value_type& value )
    {
        checkInvariant();
        M_vec( i ) = value;
    }

    /**
     * v([i1,i2,...,in]) = [value1,...,valuen]
     */
    void setVector ( int* i, int n, value_type* v )
    {
        for ( int j = 0; j < n; ++j )
            M_vec( i[j] ) = v[j];
    }

    /**
     * \f$ U+=v \f$ where \p v is a std::vector<T>
     * and you
     * want to specify WHERE to add it
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
     * \p NumericVector<T> and you
     * want to specify WHERE to add
     * the \p NumericVector<T> V
     */
    void addVector ( const Vector<value_type>& V,
                     const std::vector<size_type>& dof_indices )
    {
        FEELPP_ASSERT ( V.size() == dof_indices.size() ).error( "invalid dof indices" );

        for ( size_type i=0; i<V.size(); i++ )
            this->add ( dof_indices[i], V( i ) );
    }


    /**
     * \f$ U+=A*V\f$, add the product of a \p MatrixSparse \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const Vector<value_type>& /*V_in*/,
                     const MatrixSparse<value_type>& /*A_in*/ );
    // {
    //     FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    // }

    /**
     * \f$U+=V \f$ where U and V are type
     * uvlas::vector<T> and you
     * want to specify WHERE to add
     * the DenseVector<T> V
     */
    void addVector ( const vector_type& V,
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
                  const std::vector<size_type>& /*dof_indices*/ )
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * \f$U=V\f$, where U and V are type
     * Vector<T> and you
     * want to specify WHERE to insert
     * the Vector<T> V
     */
    void insert ( const Vector<T>& /*V*/,
                  const std::vector<size_type>& /*dof_indices*/ )
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }


    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to insert
     * the DenseVector<T> V
     */
    void insert ( const vector_type& /*V*/,
                  const std::vector<size_type>& /*dof_indices*/ )
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to insert
     * the DenseVector<T> V
     */
    void insert ( const ublas::vector<T>& V,
                  const std::vector<size_type>& dof_indices );


    /**
     * Scale each element of the
     * vector by the given factor.
     */
    void scale ( const T factor )
    {
        M_vec *= factor;
    }

    /**
     * Print the contents of the vector in Matlab's
     * sparse vector forvec. Optionally prints the
     * vector to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL", bool renumber = false ) const;

    void close() {}

    /**
     * @return the minimum element in the vector.  In case of complex
     * numbers, this returns the minimum Real part.
     */
    real_type min() const
    {
        checkInvariant();

        real_type local_min = 0;//M_vec.minCoeff();

        real_type global_min = local_min;



#ifdef FEELPP_HAS_MPI

        if ( this->comm().size() > 1 )
        {
            MPI_Allreduce ( &local_min, &global_min, 1,
                            MPI_DOUBLE, MPI_MIN, this->comm() );
        }

#endif

        return global_min;
    }
    /**
     * @return the maximum element in the vector.  In case of complex
     * numbers, this returns the maximum Real part.
     */
    real_type max() const
    {
        checkInvariant();

        real_type local_max = 0;//M_vec.maxCoeff();

        real_type global_max = local_max;


#ifdef FEELPP_HAS_MPI

        if ( this->comm().size() > 1 )
        {
            MPI_Allreduce ( &local_max, &global_max, 1,
                            MPI_DOUBLE, MPI_MAX, this->comm() );
        }

#endif

        return global_max;
    }

    /**
     * @return the \f$l_1\f$-norm of the vector, i.e.  the sum of the
     * absolute values.
     */
    real_type l1Norm() const
    {
        checkInvariant();
        double local_l1 = M_vec.array().abs().sum();

        double global_l1 = local_l1;


#ifdef FEELPP_HAS_MPI

        if ( this->comm().size() > 1 )
        {
            mpi::all_reduce( this->comm(), local_l1, global_l1, std::plus<double>() );
        }

#endif

        return global_l1;

    }

    /**
     * @return the \f$l_2\f$-norm of the vector, i.e.  the square root
     * of the sum of the squares of the elements.
     */
    real_type l2Norm() const
    {
        checkInvariant();
        real_type local_norm2 = M_vec.squaredNorm();
        real_type global_norm2 = local_norm2;


#ifdef FEELPP_HAS_MPI

        if ( this->comm().size() > 1 )
        {
            mpi::all_reduce( this->comm(), local_norm2, global_norm2, std::plus<real_type>() );
        }

#endif
        return math::sqrt( global_norm2 );
    }

    /**
     * @return the maximum absolute value of the elements of this
     * vector, which is the \f$l_\infty\f$-norm of a vector.
     */
    real_type linftyNorm() const
    {
        checkInvariant();
        real_type local_norminf = M_vec.array().abs().maxCoeff();
        real_type global_norminf = local_norminf;


#ifdef FEELPP_HAS_MPI

        if ( this->comm().size() > 1 )
        {
            mpi::all_reduce( this->comm(), local_norminf, global_norminf, mpi::maximum<real_type>() );
        }

#endif
        return global_norminf;
    }


    /**
     * @return the sum of the vector.
     */
    value_type sum() const
    {
        checkInvariant();
        value_type local_sum = M_vec.array().sum();

        value_type global_sum = local_sum;


#ifdef FEELPP_HAS_MPI

        if ( this->comm().size() > 1 )
        {
            mpi::all_reduce( this->comm(), local_sum, global_sum, std::plus<value_type>() );
        }

#endif

        return global_sum;

    }

    /**
     * @compute sqrt on each element of the vector.
     */
    FEELPP_DONT_INLINE this_type sqrt() const;


    /**
     *@compute pow on each element of the vector.
     */
    this_type pow( int n ) const;


    /**
     * \f$U(0-DIM)+=s\f$.
     * Addition of \p s to all components.
     * \note \p s is a scalar and not a vector.
     */
    void add( const T& a )
    {
        checkInvariant();

        for ( size_type i = 0; i < this->localSize(); ++i )
            M_vec.operator[]( i ) += a;

        return;
    }

    /**
     * \f$U+=V\f$.
     * Simple vector addition, equal to the \p operator+=.
     */
    void add( const Vector<T>& v )
    {
        checkInvariant();
        add( 1., v );
        return;
    }

    /**
     * \f$U+=a*V\f$.
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add( const T& a, const Vector<T>& v )
    {
        checkInvariant();

        for ( size_type i = 0; i < this->localSize(); ++i )
            M_vec.operator()( i ) += a*v( v.firstLocalIndex() + i );

        return;
    }

    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( std::vector<value_type>& /*v_local*/ ) const
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * Updates a local vector with selected values from neighboring
     * processors, as defined by \p send_list.
     */
    void localize ( const size_type first_local_idx,
                    const size_type last_local_idx,
                    const std::vector<size_type>& send_list );
    /**
     * Same, but fills a \p Vector<T> instead of
     * a \p std::vector.
     */
    void localize ( Vector<T>& v_local ) const;
    void localize ( vector_type& v_local ) const;

    /**
     * Creates a local vector \p v_local containing
     * only information relevant to this processor, as
     * defined by the \p send_list.
     */
    void localize ( Vector<T>& v_local,
                    const std::vector<size_type>& send_list ) const;



    /**
     * Creates a local copy of the global vector in
     * \p v_local only on processor \p proc_id.  By
     * default the data is sent to processor 0.  This method
     * is useful for outputting data from one processor.
     */
    void localizeToOneProcessor ( std::vector<T>& v_local,
                                  const size_type proc_id = 0 ) const;
    void localizeToOneProcessor ( vector_type& v_local,
                                  const size_type proc_id = 0 ) const;


    value_type dot( Vector<T> const& __v )
    {
        throw std::logic_error( "[vetor eigen] ERROR dot function not yet implemented" );
        return 0;
    }
    //@}



protected:

private:

    /**
     * check vector consistency
     */
    void checkInvariant() const;

private:
    vector_type M_vec;
};

/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the element product of \p v1 and \p v2
 */
template <typename T>
VectorEigen<T>
element_product( VectorEigen<T> const& v1, VectorEigen<T> const& v2 )
{
    return v1.vec().array()*v2.vec().array();
}

/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the inner product of \p v1 and \p v2
 */
template <typename T>
VectorEigen<T>
element_product( boost::shared_ptr<VectorEigen<T> > const& v1,
                 boost::shared_ptr<VectorEigen<T> > const& v2 )
{
    return v1->vec().array()*v2->vec().array();
}

/**
 * FEELPP_INSTANTIATE_VECTOREIGEN is never defined except in vectoreigen.cpp
 * where we do the instantiate. This allows to reduce the VectorEigen
 * instantiation to the strict minimum
 */
#if !defined( FEELPP_INSTANTIATE_VECTOREIGEN )
extern template class VectorEigen<double>;
extern template class VectorEigen<std::complex<double>>;
#endif

} // Feel


#endif /* __VectorEigen_H */
