/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-13

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
#include <boost/make_shared.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/vector.hpp>



namespace Feel
{
template<typename T> class VectorPetsc;
template<typename T> class VectorPetscMPI;
template<typename T> class VectorPetscMPIRange;

namespace ublas = boost::numeric::ublas;

namespace detail
{
// patch for shallow_array_adaptor from
// http://stackoverflow.com/questions/1735841/initializing-a-ublas-vector-from-a-c-array
template<typename T>
class FEELPP_EXPORT shallow_array_adaptor
    :
        public boost::numeric::ublas::shallow_array_adaptor<T>
{
public:
    typedef boost::numeric::ublas::shallow_array_adaptor<T> base_type;
    typedef typename base_type::size_type                   size_type;
    typedef typename base_type::pointer                     pointer;

    shallow_array_adaptor(size_type n) : base_type(n) {}
    shallow_array_adaptor(size_type n, pointer data) : base_type(n,data) {}
    shallow_array_adaptor(const shallow_array_adaptor& c) : base_type(c) {}

    // This function must swap the values of the items, not the data pointers.
    void swap(shallow_array_adaptor& a) {
        if (base_type::begin() != a.begin())
            std::swap_ranges(base_type::begin(), base_type::end(), a.begin());
    }
};

} // namespace detail

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
class FEELPP_EXPORT VectorUblas
    : public Vector<T>
//    , boost::addable<VectorUblas<T,Storage> >
//    , boost::subtractable<VectorUblas<T,Storage> >
//    , boost::multipliable<VectorUblas<T,Storage>, T >
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
    typedef ublas::basic_range<typename vector_type::size_type, difference_type> range_type;
    typedef ublas::basic_slice<typename vector_type::size_type, difference_type> slice_type;
    typedef Vector<value_type> clone_type;
    typedef std::shared_ptr<clone_type> clone_ptrtype;
    typedef VectorUblas<value_type, Storage> this_type;
    using self_t = this_type;
    typedef typename vector_type::iterator iterator;
    typedef typename vector_type::const_iterator const_iterator;

    struct shallow_array_adaptor
    {
        typedef ublas::vector<value_type, Feel::detail::shallow_array_adaptor<value_type> > subtype;
        typedef ublas::vector_range<subtype> rangesubtype;
        typedef ublas::vector_slice<subtype> slicesubtype;
        typedef VectorUblas<value_type,subtype> type;
    };
    static const bool is_shallow_array_adaptor_vector =
        boost::is_same<vector_type,typename this_type::shallow_array_adaptor::subtype>::value ||
        boost::is_same<vector_type,typename this_type::shallow_array_adaptor::rangesubtype>::value ||
        boost::is_same<vector_type,typename this_type::shallow_array_adaptor::slicesubtype>::value;
    static const bool is_extarray_vector = is_shallow_array_adaptor_vector;

    struct range
    {
        typedef typename mpl::if_< mpl::bool_< is_shallow_array_adaptor_vector >,
                                   typename this_type::shallow_array_adaptor::rangesubtype,
                                   ublas::vector_range<ublas::vector<value_type> > >::type subtype;
        typedef VectorUblas<value_type,subtype> type;
    };

    struct slice
    {
        typedef typename mpl::if_< mpl::bool_< is_shallow_array_adaptor_vector >,
                                   typename this_type::shallow_array_adaptor::slicesubtype,
                                   ublas::vector_slice<ublas::vector<value_type> > >::type subtype;
        typedef VectorUblas<value_type,subtype> type;
    };

    static const bool is_range_vector = boost::is_same<vector_type,typename this_type::range::subtype>::value;
    static const bool is_slice_vector = boost::is_same<vector_type,typename this_type::slice::subtype>::value;
    static const bool has_non_contiguous_ghosts = is_range_vector || is_slice_vector || is_shallow_array_adaptor_vector;


    typedef typename super1::datamap_type datamap_type;
    typedef typename super1::datamap_ptrtype datamap_ptrtype;
    using size_type = typename datamap_type::size_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    VectorUblas();

    VectorUblas( size_type __s );

    VectorUblas( datamap_ptrtype const& dm );

    VectorUblas( size_type __s, size_type __n_local );

    VectorUblas( VectorUblas const & m );

    FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, range_type const& range, datamap_ptrtype const& dm );
    VectorUblas( VectorUblas<value_type>& m, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm );
    VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm );

    FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, slice_type const& range, datamap_ptrtype const& dm );
    VectorUblas( VectorUblas<value_type>& m, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm );
    VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm );

    FEELPP_DEPRECATED VectorUblas( ublas::vector<value_type>& m, range_type const& range );
    VectorUblas( ublas::vector<value_type>& m, range_type const& range, datamap_ptrtype const& dm );

    FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, slice_type const& slice );

    FEELPP_DEPRECATED VectorUblas( ublas::vector<value_type>& m, slice_type const& slice );
    FEELPP_DEPRECATED VectorUblas( ublas::vector<value_type>& m, slice_type const& slice, datamap_ptrtype const& dm );
    VectorUblas( ublas::vector<value_type>& mActive, slice_type const& sliceActive,
                 ublas::vector<value_type>& mGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm );
    VectorUblas( typename this_type::shallow_array_adaptor::subtype& mActive, slice_type const& sliceActive,
                 typename this_type::shallow_array_adaptor::subtype& mGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm );

    VectorUblas( size_type nActiveDof, value_type* arrayActiveDof,
                 size_type nGhostDof, value_type* arrayGhostDof,
                 datamap_ptrtype const& dm );

    ~VectorUblas() override;

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
                const bool      fast=false ) override;

    /**
     * call init with n_local = N,
     */
    void init ( const size_type n,
                const bool      fast=false ) override;

    /**
     * init from a \p DataMap
     */
    void init( datamap_ptrtype const& dm ) override;


    /**
     * Creates a copy of this vector and returns it in an
     * \p shared_ptr<>.
     */
    clone_ptrtype clone() const override
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
    Vector<value_type>& operator= ( const Vector<value_type> &V ) override;
    Vector<value_type>& operator= ( const this_type &V );

    template<typename AE>
    VectorUblas<value_type, Storage>& operator=( ublas::vector_expression<AE> const& e )
    {
        M_vec.operator=( e );
        return *this;
    }

    /**
     * Access components, returns \p u(i).
     */
    T operator()( size_type i ) const override
    {
        checkIndex(i);
        if ( has_non_contiguous_ghosts )
        {
            auto const& dm = this->map();
            const size_type nLocalActiveDof = dm.nLocalDofWithoutGhost();
            if ( i < nLocalActiveDof )
                return M_vec.operator()( i );
            else
                return M_vecNonContiguousGhosts.operator()( i-nLocalActiveDof );
        }
        else
            return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Access components, returns \p u(i).
     */
    T& operator()( size_type i ) override
    {
        checkIndex(i);
        if ( has_non_contiguous_ghosts )
        {
            auto const& dm = this->map();
            const size_type nLocalActiveDof = dm.nLocalDofWithoutGhost();
            if ( i < nLocalActiveDof )
                return M_vec.operator()( i );
            else
                return M_vecNonContiguousGhosts.operator()( i-nLocalActiveDof );
        }
        else
            return M_vec.operator()( i-this->firstLocalIndex() );
    }

    /**
     * Access components, returns \p u(i).
     */
    T operator[]( size_type i ) const
    {
        //checkIndex(i);
        //return M_vec.operator()( i-this->firstLocalIndex() );
        return this->operator()( i );
    }

    /**
     * Access components, returns \p u(i).
     */
    T& operator[]( size_type i )
    {
        //checkIndex(i);
        //return M_vec.operator()( i-this->firstLocalIndex() );
        return this->operator()( i );
    }

    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    Vector<T>& operator+=( const Vector<T>& v ) override
    {
        add( 1., v );
        return *this;
    }
    /**
     * Addition operator.
     */
    self_t operator+( const self_t& v ) const
    {
        self_t w( *this );
        w += v;
        return w;
    }
    /**
     * Addition operator with a scalar
     */
    self_t operator+( T const& f ) const
    {
        self_t v( *this );
        v.setConstant( f );
        v.add( 1., *this );
        return v;
    }
    /**
     * @brief operator f + v
     * 
     * @param f scalar
     * @param v VectorUblas
     * @return self_t new vector
     */
    friend self_t operator+(T f, const self_t &v) { self_t w( v ); w.setConstant( f ); w.add(1.,v); return w;  }

    /**
     * Subtraction operator.
     * Fast equivalent to \p U.add(-1, V).
     */
    Vector<T>& operator-=( const Vector<T>& v ) override { add( -1., v ); return *this; }
    /**
     * Addition operator.
     */
    self_t operator-( const self_t& v ) const { self_t w( *this ); w -= v; return w; }
    /**
     * substraction operator with a scalar
     */
    self_t operator-( T const& f ) const
    {
        self_t v( *this );
        v.setConstant( -f );
        v.add( 1., *this );
        return v;
    }
    /**
     * @brief operator f - v
     * 
     * @param f scalar
     * @param v VectorUblas
     * @return self_t new vector
     */
    friend self_t operator-(T f, const self_t &v) { self_t w( v ); w.setConstant( f ); w.add(-1.,v); return w;  }

    /**
     * multiplication by a scalar value
     */
    Vector<T>& operator*=( T const& v )
    {
        this->scale( v );
        return *this;
    }
    /**
     * @return negative of this vector
     */
    self_t operator-() const
    {
        self_t v( *this );
        v.scale( -1. );
        return v;
    }
    /**
     * @return v * f
     */
    self_t operator*( T f ) const
    {
        self_t v( *this );
        v.scale( f );
        return v;
    }
    /**
     * @brief operator f * v
     * 
     * @param f scalar
     * @param v VectorUblas
     * @return self_t scaled vector
     */
    friend self_t operator*(T f, const self_t &v) { self_t w( v ); w.scale( f ); return w;  }

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
    iterator beginGhost()
    {
        return M_vecNonContiguousGhosts.begin();
    }
    const_iterator beginGhost() const
    {
        return M_vecNonContiguousGhosts.begin();
    }

    iterator endGhost()
    {
        return M_vecNonContiguousGhosts.end();
    }
    const_iterator endGhost() const
    {
        return M_vecNonContiguousGhosts.end();
    }

    /**
     * if the vector is a range, return the first index of the range,
     * otherwise returns 0
     */
    size_type start() const;

    /**
     * if the vector is a range, return the first index of the range,
     * otherwise returns 0
     */
    size_type startNonContiguousGhosts() const;

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
    bool isInitialized() const override
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
    bool closed() const override
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
     * Return ghost ublas vector (can be not used)
     */
    vector_type const& vecNonContiguousGhosts() const
    {
        return M_vecNonContiguousGhosts;
    }

    /**
     * Return ghost ublas vector (can be not used)
     */
    vector_type & vecNonContiguousGhosts()
    {
        return M_vecNonContiguousGhosts;
    }

    /**
     * update global values array
     */
    bool areGlobalValuesUpdated() const
    {
        return true;
    }

    /**
     * update global values
     */
    void updateGlobalValues() const
    {
        //this->localize( M_global_values );
        //M_global_values_updated = true;
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
        //M_global_values_updated = false;
    }

    /**
     * set the entries to the constant \p v
     */
    void setConstant( value_type v ) override
    {
        M_vec = ublas::scalar_vector<double>( M_vec.size(), v );
        if ( has_non_contiguous_ghosts )
            M_vecNonContiguousGhosts = ublas::scalar_vector<double>( M_vecNonContiguousGhosts.size(), v );
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
    void clear () override;

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero () override
    {
        std::fill( this->begin(), this->end(), value_type( 0 ) );
        if ( has_non_contiguous_ghosts )
            std::fill( this->beginGhost(), this->endGhost(), value_type( 0 ) );
    }

    void zero ( size_type /*start1*/, size_type /*stop1*/ ) override
    {
        this->zero();
        //ublas::project( (*this), ublas::range( start1, stop1 ) ) = ublas::zero_vector<value_type>( stop1 );
    }

    /**
     * Add \p value to the value already accumulated
     */
    void add ( const size_type i, const value_type& value ) override
    {
#if !defined(NDEBUG)
        checkInvariant();
#endif
        ( *this )( i ) += value;
    }

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    void addVector ( int* i, int n, value_type* v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override
    {
        for ( int j = 0; j < n; ++j )
            ( *this )( i[j] ) += v[j];
    }

    /**
     * set to \p value
     */
    void set ( size_type i, const value_type& value ) override
    {
#if !defined(NDEBUG)
        checkInvariant();
#endif
        ( *this )( i ) = value;
    }

    /**
     * v([i1,i2,...,in]) = [value1,...,valuen]
     */
    void setVector ( int* i, int n, value_type* v ) override
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
    void addVector ( const Vector<value_type>& /*V_in*/,
                     const MatrixSparse<value_type>& /*A_in*/ ) override
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
    void insert ( const Vector<T>& /*V*/,
                  const std::vector<size_type>& /*dof_indices*/ ) override
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
                  const std::vector<size_type>& /*dof_indices*/ ) override
    {
        FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    }

    /**
     * Scale each element of the
     * vector by the given factor.
     */
    void scale ( const T factor ) override
    {
        M_vec.operator *=( factor );
        if ( has_non_contiguous_ghosts )
            M_vecNonContiguousGhosts.operator *=( factor );

    }

    /**
     * Print the contents of the vector in Matlab's
     * sparse vector forvec. Optionally prints the
     * vector to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL", bool renumber = false ) const override;

    void close() override {}

    /**
     * @return the minimum element in the vector.  In case of complex
     * numbers, this returns the minimum Real part.
     */
    real_type min() const override
    {
        return this->min( true );
    }
    real_type min( bool parallel ) const
    {
        checkInvariant();

        size_type nActiveDof = this->map().nLocalDofWithoutGhost();
        size_type nGhostDof = this->map().nLocalGhosts();
        real_type local_min = (nActiveDof>0)?
            *std::min_element( M_vec.begin(), ( has_non_contiguous_ghosts || ( nGhostDof == 0 ) )? M_vec.end() : M_vec.begin()+nActiveDof ) :
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
    real_type max() const override
    {
        return this->max( true );
    }
    real_type max( bool parallel ) const
    {
        checkInvariant();

        size_type nActiveDof = this->map().nLocalDofWithoutGhost();
        size_type nGhostDof = this->map().nLocalGhosts();
        real_type local_max = (nActiveDof>0)?
            *std::max_element( M_vec.begin(), ( has_non_contiguous_ghosts || ( nGhostDof == 0 ) )? M_vec.end() : M_vec.begin()+nActiveDof ) :
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
    real_type l1Norm() const override
    {
        checkInvariant();

        double local_l1 = 0;
        if ( has_non_contiguous_ghosts || this->comm().size() == 1 )
            local_l1 = ublas::norm_1( M_vec );
        else
            local_l1 = ublas::norm_1( ublas::project( M_vec, ublas::range( 0,this->map().nLocalDofWithoutGhost() ) ) );

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
    real_type l2Norm() const override
    {
        checkInvariant();
        real_type local_norm2 = 0;

        if ( has_non_contiguous_ghosts || this->comm().size() == 1 )
        {
            local_norm2 = ublas::inner_prod( M_vec, M_vec );
        }
        else
        {
            auto vecActiveDof = ublas::project( M_vec, ublas::range( 0,this->map().nLocalDofWithoutGhost() ) );
            local_norm2 = ublas::inner_prod( vecActiveDof,vecActiveDof );
        }
        real_type global_norm2 = local_norm2;

#ifdef FEELPP_HAS_MPI
        if ( this->comm().size() > 1 )
            mpi::all_reduce( this->comm(), local_norm2, global_norm2, std::plus<real_type>() );
#endif

        return math::sqrt( global_norm2 );
    }

    /**
     * @return the maximum absolute value of the elements of this
     * vector, which is the \f$l_\infty\f$-norm of a vector.
     */
    real_type linftyNorm() const override
    {
        checkInvariant();
        real_type local_norminf = 0;
        if ( has_non_contiguous_ghosts || this->comm().size() == 1 )
            local_norminf = ublas::norm_inf( M_vec );
        else
            local_norminf = ublas::norm_inf( ublas::project( M_vec, ublas::range( 0,this->map().nLocalDofWithoutGhost() ) ) );

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
    value_type sum() const override
    {
        checkInvariant();
        value_type local_sum = 0;

        if ( has_non_contiguous_ghosts || this->comm().size() == 1 )
        {
            local_sum = ublas::sum( M_vec );
        }
        else
        {
            local_sum = ublas::sum( ublas::project( M_vec, ublas::range( 0,this->map().nLocalDofWithoutGhost() ) ) );
        }
        value_type global_sum = local_sum;
#ifdef FEELPP_HAS_MPI
        if ( this->comm().size() > 1 )
            mpi::all_reduce( this->comm(), local_sum, global_sum, std::plus<value_type>() );
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
    void add( const T& a ) override
    {
        checkInvariant();
#if 1
        for ( size_type i = 0; i < this->localSize(); ++i )
            this->operator[]( i ) += a;
        //M_vec.operator[]( i ) += a;

        return;
#else
        M_vec += a*ublas::unit_vector<value_type>();
        if ( has_non_contiguous_ghosts )
            M_vecNonContiguousGhosts += a;
#endif
    }

    /**
     * \f$U+=V\f$.
     * Simple vector addition, equal to the \p operator+=.
     */
    void add( const Vector<T>& v ) override
    {
        add( 1., v );
        return;
    }

    /**
     * \f$U+=a*V\f$.
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add( const T& a, const Vector<T>& v ) override;

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
    void localize ( ublas::vector<value_type,Feel::detail::shallow_array_adaptor<T> >& v_local ) const {}

    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( ublas::vector_range<ublas::vector<value_type> >& v_local ) const;
    void localize ( ublas::vector_range<ublas::vector<value_type,Feel::detail::shallow_array_adaptor<T> > >& v_local ) const {}

    /**
     * Creates a copy of the global vector in the
     * local vector \p v_local.
     */
    void localize ( ublas::vector_slice<ublas::vector<value_type> >& v_local ) const;
    void localize ( ublas::vector_slice<ublas::vector<value_type,Feel::detail::shallow_array_adaptor<T> > >& v_local ) const {}

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


    value_type dot( Vector<T> const& __v ) const override;
    //@}

#ifdef FEELPP_HAS_HDF5
    void saveHDF5( std::string const& filename, std::string tableName = "element", bool appendMode = false ) const;
    void loadHDF5( std::string const& filename, std::string tableName = "element" );
#endif

    //! create a VectorUblas as a view of another Vector type (as VectorPetsc)
    static
    typename this_type::shallow_array_adaptor::type
    createView( Vector<T> const& vec );

protected:

private:

    template <typename OtherStorage>
    void assignWithUblasImpl( VectorUblas<T,OtherStorage> const& v );
    template <typename OtherStorage>
    void addWithUblasImpl( const T& a, VectorUblas<T,OtherStorage> const& v );
    template <typename OtherStorage>
    value_type dotWithUblasImpl( VectorUblas<T,OtherStorage> const& v ) const;

    /**
     * check vector consistency
     */
    void checkInvariant() const;

    inline void checkIndex( size_type i ) const
        {
            #if 0
            DCHECK(  this->isInitialized() ) << "vector not initialized";
            DCHECK (( i >= this->firstLocalIndex() ) &&
                    ( i <=  this->lastLocalIndex() ) )
                << "invalid index " << i << " min=" <<  this->firstLocalIndex() << " max=" << this->lastLocalIndex() ;
            #endif
        }


private:

    //! vector with in first active dofs and then ghost, all dofs is contiguous, if Storage is not a view
    //! else is Storage is a view, the vector represent only the active dofs contiguous
    vector_type M_vec;
    //! non contiguous ghost values : defined only with range view vector
    vector_type M_vecNonContiguousGhosts;

};


template<typename T, typename Storage>
template <typename OtherStorage>
void
VectorUblas<T,Storage>::assignWithUblasImpl( VectorUblas<T,OtherStorage> const& v )
{
    size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
    size_type nLocalGhosts = this->map().nLocalGhosts();
    size_type nLocalDofWithoutGhostV = v.map().nLocalDofWithoutGhost();
    size_type nLocalGhostsV = v.map().nLocalGhosts();

    static const bool otherstorage_has_non_contiguous_ghosts = VectorUblas<T,OtherStorage>::has_non_contiguous_ghosts;

    if ( this->comm().localSize() == 1 )
    {
        this->vec() = v.vec();
    }
    else if ( has_non_contiguous_ghosts )
    {
        if ( otherstorage_has_non_contiguous_ghosts )
        {
            if ( nLocalDofWithoutGhost > 0 )
                this->vec() = v.vec();
            if ( nLocalGhosts > 0 )
                this->vecNonContiguousGhosts() = v.vecNonContiguousGhosts();
        }
        else
        {
            if ( nLocalDofWithoutGhostV > 0 )
                this->vec() = ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostV ) );
            if ( nLocalGhostsV > 0 )
                this->vecNonContiguousGhosts() = ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostV, nLocalDofWithoutGhostV+nLocalGhostsV ) );
        }
    }
    else
    {
        if ( otherstorage_has_non_contiguous_ghosts )
        {
            if ( nLocalDofWithoutGhost > 0 )
                ublas::project( this->vec(), ublas::range( 0, nLocalDofWithoutGhost ) ) = v.vec();
            if ( nLocalGhosts > 0 )
                ublas::project( this->vec(), ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) = v.vecNonContiguousGhosts();
        }
        else
        {
            this->vec() = v.vec();
        }
    }

}


template<typename T, typename Storage>
template <typename OtherStorage>
void
VectorUblas<T,Storage>::addWithUblasImpl( const T& a, VectorUblas<T,OtherStorage> const& v )
{
    size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
    size_type nLocalGhosts = this->map().nLocalGhosts();
    size_type nLocalDofWithoutGhostV = v.map().nLocalDofWithoutGhost();
    size_type nLocalGhostsV = v.map().nLocalGhosts();

    static const bool otherstorage_has_non_contiguous_ghosts = VectorUblas<T,OtherStorage>::has_non_contiguous_ghosts;

    if ( this->comm().localSize() == 1 )
    {
        this->vec() += a*v.vec();
    }
    else if ( has_non_contiguous_ghosts )
    {
        if ( otherstorage_has_non_contiguous_ghosts )
        {
            if ( nLocalDofWithoutGhost > 0 )
                this->vec() += a*v.vec();
            if ( nLocalGhosts > 0 )
                this->vecNonContiguousGhosts() += a*v.vecNonContiguousGhosts();
        }
        else
        {
            if ( nLocalDofWithoutGhostV > 0 )
                this->vec() += a*ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostV ) );
            if ( nLocalGhostsV > 0 )
                this->vecNonContiguousGhosts() += a*ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostV, nLocalDofWithoutGhostV+nLocalGhostsV ) );
        }
    }
    else
    {
        if ( otherstorage_has_non_contiguous_ghosts )
        {
            if ( nLocalDofWithoutGhost > 0 )
                ublas::project( this->vec(), ublas::range( 0, nLocalDofWithoutGhost ) ) += a*v.vec();
            if ( nLocalGhosts > 0 )
                ublas::project( this->vec(), ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) += a*v.vecNonContiguousGhosts();
        }
        else
        {
            this->vec() += a*v.vec();
        }
    }
}

template<typename T, typename Storage>
template <typename OtherStorage>
typename VectorUblas<T,Storage>::value_type
VectorUblas<T,Storage>::dotWithUblasImpl( VectorUblas<T,OtherStorage> const& v ) const
{
    static const bool otherstorage_has_non_contiguous_ghosts = VectorUblas<T,OtherStorage>::has_non_contiguous_ghosts;
    value_type localResult = 0;
    size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
    size_type nLocalDofWithoutGhostV = v.map().nLocalDofWithoutGhost();

    if ( this->comm().localSize() == 1 )
    {
        localResult = ublas::inner_prod( this->vec(), v.vec() );
    }
    else if ( has_non_contiguous_ghosts )
    {
        if ( otherstorage_has_non_contiguous_ghosts )
        {
            localResult = ublas::inner_prod( this->vec(), v.vec() );
        }
        else
        {
            localResult = ublas::inner_prod( this->vec(),
                                             ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostV ) ) );
        }
    }
    else
    {
        if ( otherstorage_has_non_contiguous_ghosts )
        {
            localResult = ublas::inner_prod( ublas::project( this->vec(), ublas::range( 0, nLocalDofWithoutGhost ) ),
                                             v.vec() );
        }
        else
        {
            localResult = ublas::inner_prod( ublas::project( this->vec(), ublas::range( 0, nLocalDofWithoutGhost ) ),
                                             ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostV ) ) );
        }
    }

    value_type globalResult = localResult;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localResult, globalResult, std::plus<value_type>() );
#endif

    return globalResult;
}


/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the element product of \p v1 and \p v2
 */
template <typename T, typename Storage1, typename Storage2>
VectorUblas<T>
element_product( VectorUblas<T,Storage1> const& v1, VectorUblas<T,Storage2> const& v2 )
{
    FEELPP_ASSERT( v1.localSize() == v2.localSize() &&
                   v1.size() == v2.size() )
    ( v1.localSize() )( v2.localSize() )
    ( v1.size() )( v2.size() ).error( "incompatible vector sizes" );

    typedef typename type_traits<T>::real_type real_type;

    static const bool storage1_has_non_contiguous_ghosts = VectorUblas<T,Storage1>::has_non_contiguous_ghosts;
    static const bool storage2_has_non_contiguous_ghosts = VectorUblas<T,Storage2>::has_non_contiguous_ghosts;

    VectorUblas<real_type> _t( v1.mapPtr() );

    size_type nLocalDofWithoutGhost = _t.map().nLocalDofWithoutGhost();
    size_type nLocalGhosts = _t.map().nLocalGhosts();

    if ( _t.comm().localSize() == 1 )
    {
        _t.vec() = ublas::element_prod( v1.vec(),v2.vec() );
    }
    else if ( storage1_has_non_contiguous_ghosts == storage2_has_non_contiguous_ghosts )
    {
        if ( storage1_has_non_contiguous_ghosts )
        {
            if ( nLocalDofWithoutGhost > 0 )
                ublas::project( _t.vec(), ublas::range( 0, nLocalDofWithoutGhost ) ) =
                    ublas::element_prod( v1.vec(),v2.vec() );
            if ( nLocalGhosts > 0 )
                ublas::project( _t.vec(), ublas::range( nLocalDofWithoutGhost,nLocalDofWithoutGhost+nLocalGhosts ) ) =
                    ublas::element_prod( v1.vecNonContiguousGhosts(),v2.vecNonContiguousGhosts() );
        }
        else
        {
            _t.vec() = ublas::element_prod( v1.vec(),v2.vec() );
        }
    }
    else if ( storage1_has_non_contiguous_ghosts )
    {
        size_type nLocalDofWithoutGhostV2 = v1.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsV2 = v1.map().nLocalGhosts();

        if ( nLocalDofWithoutGhost > 0 )
            ublas::project( _t.vec(), ublas::range( 0,nLocalDofWithoutGhost ) ) =
                ublas::element_prod( v1.vec(),
                                     ublas::project( v2.vec(), ublas::range( 0, nLocalDofWithoutGhostV2 ) ) );
        if ( nLocalGhosts > 0 )
            ublas::project( _t.vec(), ublas::range( nLocalDofWithoutGhost,nLocalDofWithoutGhost+nLocalGhosts ) ) =
                ublas::element_prod( v1.vecNonContiguousGhosts(),
                                     ublas::project( v2.vec(), ublas::range( nLocalDofWithoutGhostV2,nLocalDofWithoutGhostV2+nLocalGhostsV2 ) ) );
    }
    else if ( storage2_has_non_contiguous_ghosts )
    {
        size_type nLocalDofWithoutGhostV1 = v1.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsV1 = v1.map().nLocalGhosts();

        if ( nLocalDofWithoutGhost > 0 )
            ublas::project( _t.vec(), ublas::range( 0,nLocalDofWithoutGhost ) ) =
                ublas::element_prod( ublas::project( v1.vec(), ublas::range( 0, nLocalDofWithoutGhostV1 ) ),
                                     v2.vec() );
        if ( nLocalGhosts > 0 )
            ublas::project( _t.vec(), ublas::range( nLocalDofWithoutGhost,nLocalDofWithoutGhost+nLocalGhosts ) ) =
                ublas::element_prod( ublas::project( v1.vec(), ublas::range( nLocalDofWithoutGhostV1,nLocalDofWithoutGhostV1+nLocalGhostsV1 ) ),
                                     v2.vecNonContiguousGhosts() );
    }
    else
    {
        size_type s = v1.localSize();
        size_type start = v1.firstLocalIndex();
        for ( size_type i = 0; i < s; ++i )
            _t.operator()( start+i ) = v1.operator()( start + i )* v2.operator()( start + i );
    }

    return _t;
}

/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the inner product of \p v1 and \p v2
 */
template <typename T, typename Storage1, typename Storage2>
VectorUblas<T>
element_product( std::shared_ptr<VectorUblas<T,Storage1> > const& v1,
                 std::shared_ptr<VectorUblas<T,Storage2> > const& v2 )
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

}

#if FEELPP_HAS_PETSC
#include <feel/feelalg/vectorpetsc.hpp>

namespace Feel
{
/**
 * returns a VectorPetsc from a VectorUblas
 *
 * here is a sample code:
 * @code
 * auto Xh = Pch<1>( mesh );
 * auto v = Xh->element();
 * auto b = backend(); // default backend is petsc type
 * auto vp = toPETSc( v ); // get the VectorPetsc
 * @endcode
 *
 * \warning one must be careful that the VectorUblas will provide contiguous
 * data access.
 */
template<typename T,typename Storage>
inline VectorPetsc<T>
toPETSc( VectorUblas<T,Storage> & v )
{
    if ( v.comm().size() > 1 )
    {
        if ( VectorUblas<T,Storage>::is_range_vector || VectorUblas<T,Storage>::is_extarray_vector )
            return VectorPetscMPIRange<T>( v );
        else
            return VectorPetscMPI<T>( v );
    }
    else
        return VectorPetsc<T>( v );
}

/**
 * returns a VectorPetsc shared_ptr from a VectorUblas
 *
 * here is a sample code:
 * @code
 * auto Xh = Pch<1>( mesh );
 * auto v = Xh->element();
 * auto b = backend(); // default backend is petsc type
 * auto vp = toPETScPtr( v ); // get the shared_ptr<VectorPetsc<>>
 * @endcode
 *
 * \warning one must be careful that the VectorUblas will provide contiguous
 * data access.
 */
template<typename T,typename Storage>
inline std::shared_ptr<VectorPetsc<T>>
toPETScPtr( VectorUblas<T,Storage> const& v )
{
    if ( v.comm().size() > 1 )
    {
        if ( VectorUblas<T,Storage>::is_range_vector || VectorUblas<T,Storage>::is_extarray_vector )
            return std::make_shared<VectorPetscMPIRange<T>>( v );
        else
            return std::make_shared<VectorPetscMPI<T>>( v );
    }
    else
        return std::make_shared<VectorPetsc<T>>( v );
}

/**
 * returns a VectorPetsc shared_ptr from a VectorUblas
 *
 * here is a sample code:
 * @code
 * auto Xh = Pch<1>( mesh );
 * auto v = Xh->elementPtr();
 * auto b = backend(); // default backend is petsc type
 * auto vp = toPETSc( v ); // get the shared_ptr<VectorPetsc<>>
 * @endcode
 *
 * \warning one must be careful that the VectorUblas will provide contiguous
 * data access.
 */
template<typename T,typename Storage>
inline std::shared_ptr<VectorPetsc<T>>
toPETSc( std::shared_ptr<VectorUblas<T,Storage>> & v )
{
    if ( v->comm().size() > 1 )
    {
        if ( VectorUblas<T,Storage>::is_range_vector || VectorUblas<T,Storage>::is_extarray_vector )
            return std::make_shared<VectorPetscMPIRange<T>>( *v );
        else
            return std::make_shared<VectorPetscMPI<T>>( *v );
    }
    else
        return std::make_shared<VectorPetsc<T>>( *v );
}



} // Feel
#endif


#endif /* __VectorUblas_H */
