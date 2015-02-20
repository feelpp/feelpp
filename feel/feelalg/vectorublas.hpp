/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-13

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2010 Universit� Joseph Fourier (Grenoble I)

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
   \file vectorublas.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-11-13
 */
#ifndef __VectorUblas_H
#define __VectorUblas_H 1

#include <set>
#include <boost/operators.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/vector.hpp>



namespace Feel
{
namespace ublas = boost::numeric::ublas;
/*!
 * \class VectorUblas
 * \brief interface to vector
 *
 * \code
 * VectorUblas<T> m;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T, typename Storage = ublas::vector<T> >
class VectorUblas
    : public Vector<T>
    , boost::addable<VectorUblas<T,Storage> >
    , boost::subtractable<VectorUblas<T,Storage> >
    , boost::multipliable<VectorUblas<T,Storage>, T >
{
    typedef Vector<T> super1;
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename type_traits<value_type>::real_type real_type;

    typedef Storage vector_type;
    typedef typename vector_type::difference_type difference_type;
    typedef ublas::basic_range<size_type, difference_type> range_type;
    typedef ublas::basic_slice<size_type, difference_type> slice_type;
    typedef Vector<value_type> clone_type;
    typedef boost::shared_ptr<clone_type> clone_ptrtype;
    typedef VectorUblas<value_type, Storage> this_type;

    typedef typename vector_type::iterator iterator;
    typedef typename vector_type::const_iterator const_iterator;

    struct range
    {
        typedef ublas::vector_range<ublas::vector<value_type> > subtype;
        typedef VectorUblas<value_type,subtype> type;
    };

    struct slice
    {
        typedef ublas::vector_slice<ublas::vector<value_type> > subtype;
        typedef VectorUblas<value_type,subtype> type;
    };


    typedef typename super1::datamap_type datamap_type;
    typedef typename super1::datamap_ptrtype datamap_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    VectorUblas();

    VectorUblas( size_type __s );

    VectorUblas( datamap_ptrtype const& dm );

    VectorUblas( size_type __s, size_type __n_local );

    VectorUblas( VectorUblas const & m );

    VectorUblas( VectorUblas<value_type>& m, range_type const& range, datamap_ptrtype const& dm );

    VectorUblas( ublas::vector<value_type>& m, range_type const& range );

    VectorUblas( VectorUblas<value_type>& m, slice_type const& slice );

    VectorUblas( ublas::vector<value_type>& m, slice_type const& slice );


    ~VectorUblas();

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


    /**
     * Creates a copy of this vector and returns it in an
     * \p shared_ptr<>.
     */
    clone_ptrtype clone() const
    {
        return clone_ptrtype( new this_type( *this ) );
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     *  \f$U = V\f$: copy all components.
     */
    Vector<value_type>& operator= ( const Vector<value_type> &V );

    template<typename AE>
    VectorUblas<value_type, Storage>& operator=( ublas::vector_expression<AE> const& e )
    {
        this->outdateGlobalValues();
        M_vec.operator=( e );
        return *this;
    }

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
        this->outdateGlobalValues();
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
        this->outdateGlobalValues();
        return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    Vector<T>& operator+=( const Vector<T>& v )
    {
        checkInvariant();
        this->outdateGlobalValues();
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
        this->outdateGlobalValues();
        add( -1., v );

        return *this;
    }

    /**
     * multiplication by a scalar value
     */
    Vector<T>& operator*=( T const& v )
    {
        checkInvariant();
        this->outdateGlobalValues();
        this->scale( v );

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    iterator begin()
    {
        return M_vec.begin();
    }
    const_iterator begin() const
    {
        return M_vec.begin();
    }

    iterator end()
    {
        return M_vec.end();
    }
    const_iterator end() const
    {
        return M_vec.end();
    }

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
     * \c close the ublas vector, that will copy the content of write
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
     * Returns the read optimized ublas vector.
     */
    vector_type const& vec () const
    {
        return M_vec;
    }

    /**
     * Returns the read optimized ublas vector.
     */
    vector_type & vec ()
    {
        return M_vec;
    }

    /**
     * update global values array
     */
    bool areGlobalValuesUpdated() const
    {
        return M_global_values_updated;
    }

    /**
     * update global values
     */
    void updateGlobalValues() const
    {
        //this->localize( M_global_values );
        M_global_values_updated = true;
    }

    /**
     * get the \p i -th global value
     */
    value_type globalValue( size_type i ) const
    {
        return this->operator()( i );
    }
    //{ return M_global_values( i ); }

    //@


    /** @name  Mutators
     */
    //@{

    /**
     * outdate global values array e.g. they must be update
     */
    void outdateGlobalValues()
    {
        M_global_values_updated = false;
    }

    /**
     * set the \p i -th global value
     */
    void setGlobalValue( size_type i, value_type v ) const
    {
        //M_global_values( i ) = v;
    }

    /**
     * set the entries to the constant \p v
     */
    void setConstant( value_type v )
    {
        M_vec = ublas::scalar_vector<double>( M_vec.size(), v );
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
        //M_vec.clear();
        this->outdateGlobalValues();
        std::fill( this->begin(), this->end(), value_type( 0 ) );
    }

    void zero ( size_type /*start1*/, size_type /*stop1*/ )
    {
        //ublas::project( (*this), ublas::range( start1, stop1 ) ) = ublas::zero_vector<value_type>( stop1 );
    }

    /**
     * Add \p value to the value already accumulated
     */
    void add ( const size_type i, const value_type& value )
    {
        checkInvariant();
        this->outdateGlobalValues();
        ( *this )( i ) += value;
    }

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    void addVector ( int* i, int n, value_type* v )
    {
        for ( int j = 0; j < n; ++j )
            ( *this )( i[j] ) += v[j];
    }

    /**
     * set to \p value
     */
    void set ( size_type i, const value_type& value )
    {
        checkInvariant();
        this->outdateGlobalValues();
        ( *this )( i ) = value;
    }

    /**
     * v([i1,i2,...,in]) = [value1,...,valuen]
     */
    void setVector ( int* i, int n, value_type* v )
    {
        for ( int j = 0; j < n; ++j )
            ( *this )( i[j] ) = v[j];
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
        this->outdateGlobalValues();

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
        this->outdateGlobalValues();

        for ( size_type i=0; i<V.size(); i++ )
            this->add ( dof_indices[i], V( i ) );
    }


    /**
     * \f$ U+=A*V\f$, add the product of a \p MatrixSparse \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const Vector<value_type>& /*V_in*/,
                     const MatrixSparse<value_type>& /*A_in*/ )
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

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
        this->outdateGlobalValues();

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
    void insert ( const ublas::vector<T>& /*V*/,
                  const std::vector<size_type>& /*dof_indices*/ )
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * Scale each element of the
     * vector by the given factor.
     */
    void scale ( const T factor )
    {
        this->outdateGlobalValues();
        M_vec.operator *=( factor );
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
        return this->min( true );
    }
    real_type min( bool parallel ) const
    {
        checkInvariant();

        real_type local_min = (this->localSize()>0)?
            *std::min_element( M_vec.begin(), M_vec.end() ) :
            std::numeric_limits<real_type>::max() ;
        real_type global_min = local_min;

#ifdef FEELPP_HAS_MPI
        if ( parallel && this->comm().size() > 1 )
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
        return this->max( true );
    }
    real_type max( bool parallel ) const
    {
        checkInvariant();

        real_type local_max = (this->localSize()>0)?
            *std::max_element( M_vec.begin(), M_vec.end() ) :
            std::numeric_limits<real_type>::min() ;
        real_type global_max = local_max;

#ifdef FEELPP_HAS_MPI
        if ( parallel && this->comm().size() > 1 )
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
        double local_l1 = ublas::norm_1( M_vec );

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
        real_type local_norm2 = 0, global_norm2=0;

        if ( this->comm().size() == 1 )
            {
                local_norm2 = ublas::inner_prod( M_vec, M_vec );
                global_norm2 = local_norm2;
            }
        else
            {
                size_type s = this->localSize();
                size_type start = this->firstLocalIndex();
                for ( size_type i = 0; i < s; ++i )
                    {
                        if ( !this->localIndexIsGhost( start + i ) )
                            local_norm2 += std::pow(M_vec.operator()( start + i ),2);
                    }
#ifdef FEELPP_HAS_MPI
                mpi::all_reduce( this->comm(), local_norm2, global_norm2, std::plus<real_type>() );
#endif
            }

        return math::sqrt( global_norm2 );
    }

    /**
     * @return the maximum absolute value of the elements of this
     * vector, which is the \f$l_\infty\f$-norm of a vector.
     */
    real_type linftyNorm() const
    {
        checkInvariant();
        real_type local_norminf = ublas::norm_inf( M_vec );
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
        value_type local_sum = 0, global_sum=0;

        if ( this->comm().size() == 1 )
        {
            local_sum = ublas::sum( M_vec );
            global_sum = local_sum;
        }
        else
        {
            size_type s = this->localSize();
            size_type start = this->firstLocalIndex();
            for ( size_type i = 0; i < s; ++i )
            {
                if ( !this->localIndexIsGhost( start + i ) )
                    local_sum += M_vec.operator()( start + i );
            }
#ifdef FEELPP_HAS_MPI
            mpi::all_reduce( this->comm(), local_sum, global_sum, std::plus<value_type>() );
#endif
        }

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
        this->outdateGlobalValues();

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
        this->outdateGlobalValues();
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
        this->outdateGlobalValues();

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
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( ublas::vector<value_type>& v_local ) const;

    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( ublas::vector_range<ublas::vector<value_type> >& v_local ) const;

    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( ublas::vector_slice<ublas::vector<value_type> >& v_local ) const;

    /**
     * Same, but fills a \p NumericVector<T> instead of
     * a \p std::vector.
     */
    void localize ( Vector<T>& v_local ) const;

    /**
     * Creates a local vector \p v_local containing
     * only information relevant to this processor, as
     * defined by the \p send_list.
     */
    void localize ( Vector<T>& v_local,
                    const std::vector<size_type>& send_list ) const;

    /**
     * Updates a local vector with selected values from neighboring
     * processors, as defined by \p send_list.
     */
    void localize ( const size_type first_local_idx,
                    const size_type last_local_idx,
                    const std::vector<size_type>& send_list );

    /**
     * Creates a local copy of the global vector in
     * \p v_local only on processor \p proc_id.  By
     * default the data is sent to processor 0.  This method
     * is useful for outputting data from one processor.
     */
    void localizeToOneProcessor ( ublas::vector<T>& v_local,
                                  const size_type proc_id = 0 ) const;

    /**
     * Creates a local copy of the global vector in
     * \p v_local only on processor \p proc_id.  By
     * default the data is sent to processor 0.  This method
     * is useful for outputting data from one processor.
     */
    void localizeToOneProcessor ( std::vector<T>& v_local,
                                  const size_type proc_id = 0 ) const;


    value_type dot( Vector<T> const& __v );
    //@}



protected:

private:

    /**
     * check vector consistency
     */
    void checkInvariant() const;

private:
    vector_type M_vec;
    mutable bool M_global_values_updated;
    //mutable ublas::vector<value_type> M_global_values;
};

/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the element product of \p v1 and \p v2
 */
template <typename T>
VectorUblas<T>
element_product( VectorUblas<T> const& v1, VectorUblas<T> const& v2 )
{
    FEELPP_ASSERT( v1.localSize() == v2.localSize() &&
                   v1.size() == v2.size() )
    ( v1.localSize() )( v2.localSize() )
    ( v1.size() )( v2.size() ).error( "incompatible vector sizes" );

    typedef typename type_traits<T>::real_type real_type;

    VectorUblas<real_type> _t( v1.mapPtr() );
    size_type s = v1.localSize();
    size_type start = v1.firstLocalIndex();

    for ( size_type i = 0; i < s; ++i )
        _t.operator()( start+i ) = v1.operator()( start + i )* v2.operator()( start + i );

    return _t;
}

/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the inner product of \p v1 and \p v2
 */
template <typename T>
VectorUblas<T>
element_product( boost::shared_ptr<VectorUblas<T> > const& v1,
                 boost::shared_ptr<VectorUblas<T> > const& v2 )
{
    return element_product( *v1, *v2 );
}

/**
 * FEELPP_INSTANTIATE_VECTORUBLAS is never defined except in vectorublas.cpp
 * where we do the instantiate. This allows to reduce the VectorUblas
 * instantiation to the strict minimum
 */
#if !defined( FEELPP_INSTANTIATE_VECTORUBLAS )
extern template class VectorUblas<double,ublas::vector<double> >;
extern template class VectorUblas<double,ublas::vector_range<ublas::vector<double> > >;
extern template class VectorUblas<double,ublas::vector_slice<ublas::vector<double> > >;
#endif

} // Feel


#endif /* __VectorUblas_H */
