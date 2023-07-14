/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
       Date: 2021-10-29

  Copyright (C) 2021 Feel++ Consortium
  Copyright (C) 2021 INRIA

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
#ifndef _FEELPP_VECTORUBLAS_HPP
#define _FEELPP_VECTORUBLAS_HPP 1

#include <variant>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
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
class FEELPP_EXPORT shallow_array_adaptor:
        public boost::numeric::ublas::shallow_array_adaptor<T>
{
public:
    typedef boost::numeric::ublas::shallow_array_adaptor<T> base_type;
    typedef typename base_type::size_type                   size_type;
    typedef typename base_type::pointer                     pointer;

    shallow_array_adaptor() : base_type() {}
    shallow_array_adaptor(size_type n) : base_type(n) {}
    shallow_array_adaptor(size_type n, pointer data) : base_type(n,data) {}
    shallow_array_adaptor(const shallow_array_adaptor& c) : base_type(c) {}

    // This function must swap the values of the items, not the data pointers.
    void swap(shallow_array_adaptor& a) {
        if (base_type::begin() != a.begin())
            std::swap_ranges(base_type::begin(), base_type::end(), a.begin());
    }
};

/*-----------------------------------------------------------------------------*/
// Forward declarations
template< typename T >
class VectorUblasBase;
template< typename T >
class VectorUblasContiguousGhostsBase;
template< typename T, typename Storage >
class VectorUblasContiguousGhosts;
template< typename T >
class VectorUblasNonContiguousGhostsBase;
template< typename T, typename Storage >
class VectorUblasNonContiguousGhosts;
template< typename T, typename Storage >
class VectorUblasRange;
template< typename T, typename Storage >
class VectorUblasSlice;

} // namespace detail

template< typename T > class VectorUblas;

template< typename T > 
VectorUblas<T> element_product( const VectorUblas<T> & v1, const VectorUblas<T> & v2 );

/*-----------------------------------------------------------------------------*/
template< typename T >
class FEELPP_EXPORT VectorUblas : public Vector<T>
{
    public:
        // Typedefs
        typedef VectorUblas<T> self_type;
        typedef Vector<T> super_type;

        using clone_ptrtype = typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;

        typedef Feel::detail::VectorUblasBase<value_type> vector_impl_type;
        typedef std::unique_ptr<vector_impl_type> vector_impl_ptrtype;
        typedef std::unique_ptr<const vector_impl_type> vector_impl_cstptrtype;

        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;
        
        // Ublas range and slice types
        typedef ublas::basic_range<typename ublas::vector<value_type>::size_type, typename ublas::vector<value_type>::difference_type> range_type;
        typedef ublas::basic_slice<typename ublas::vector<value_type>::size_type, typename ublas::vector<value_type>::difference_type> slice_type;

        // Iterators types
        typedef typename Feel::detail::VectorUblasBase<T>::iterator_type iterator;
        typedef typename Feel::detail::VectorUblasBase<T>::const_iterator_type const_iterator;

    public:
        // Constructors/Destructor
        VectorUblas();
        VectorUblas( size_type s );
        VectorUblas( datamap_ptrtype const& dm );
        VectorUblas( size_type s, size_type n_local );
        //FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, range_type const& range, datamap_ptrtype const& dm );
        VectorUblas( VectorUblas<value_type>& v, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm );
        //VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm );
        //FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, slice_type const& range, datamap_ptrtype const& dm );
        VectorUblas( VectorUblas<value_type>& v, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm );
        //VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm );
        //FEELPP_DEPRECATED VectorUblas( ublas::vector<value_type>& m, range_type const& range );
        VectorUblas( ublas::vector<value_type>& m, range_type const& range, datamap_ptrtype const& dm );
        //FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, slice_type const& slice );
        //FEELPP_DEPRECATED VectorUblas( ublas::vector<value_type>& m, slice_type const& slice );
        //FEELPP_DEPRECATED VectorUblas( ublas::vector<value_type>& m, slice_type const& slice, datamap_ptrtype const& dm );
        VectorUblas( ublas::vector<T>& vActive, slice_type const& sliceActive,
                     ublas::vector<T>& vGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm );
        //VectorUblas( typename this_type::shallow_array_adaptor::subtype& mActive, slice_type const& sliceActive,
                     //typename this_type::shallow_array_adaptor::subtype& mGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm );
        VectorUblas( size_type nActiveDof, value_type * arrayActiveDof,
                     size_type nGhostDof, value_type * arrayGhostDof,
                     datamap_ptrtype const& dm );

        VectorUblas( const VectorUblas<T> & other ): super_type( other ), M_vectorImpl( other.M_vectorImpl ? other.M_vectorImpl->clonePtr() : nullptr ) { }
        VectorUblas( VectorUblas<T> && other ): super_type( std::move( other ) ) { swap( *this, other ); }
        friend void swap( VectorUblas<T> & first, VectorUblas<T> & other ) { using std::swap; swap( first.M_vectorImpl, other.M_vectorImpl ); }

        ~VectorUblas() override = default;

        Vector<T> & operator=( const Vector<T> & v ) override;
        //VectorUblas<T> & operator=( VectorUblas<T> other ) { swap( *this, other ); return *this; }
        VectorUblas<T> & operator=( const VectorUblas<T> & other );

        clone_ptrtype clone() const override { return clone_ptrtype( new self_type( *this ) ); }

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override { return M_vectorImpl->init( n, n_local, fast ); }
        void init( const size_type n, const bool fast = false ) override { return M_vectorImpl->init( n, fast ); }
        void init( const datamap_ptrtype & dm ) override { return M_vectorImpl->init( dm ); }
        
        void resize( size_type n ) { return M_vectorImpl->resize( n ); }
        void clear() override { return M_vectorImpl->clear(); }

        void setMap( const datamap_ptrtype & dm ) { super_type::setMap( dm ); M_vectorImpl->setMap( dm ); }

        const vector_impl_type & vectorImpl() const { return *M_vectorImpl; }

        // Status API
        bool isInitialized() const override { return true; }
        void close() override { }
        bool closed() const override { return true; }
        bool areGlobalValuesUpdated() const { return true; }
        void updateGlobalValues() const { }
        void outdateGlobalValues() { }

        //void checkInvariants() const { DCHECK( M_vectorImpl ) << "vector impl not initialized"; if( M_vectorImpl ) M_vectorImpl->checkInvariants(); }

        // Operators API
        virtual value_type operator()( size_type i ) const override { return M_vectorImpl->operator()( i ); }
        virtual value_type& operator()( size_type i ) override { return M_vectorImpl->operator()( i ); }
        value_type operator[]( size_type i ) const { return this->operator()( i ); }
        value_type& operator[]( size_type i ) { return this->operator()( i ); }

        Vector<T>& operator+=( const Vector<T>& v ) override { this->add( v ); return *this; }
        Vector<T>& operator-=( const Vector<T>& v ) override { this->sub( v ); return *this; }
        Vector<T>& operator*=( const value_type v ) { this->scale( v ); return *this; }

        self_type operator+( const self_type & v ) const { return self_type( M_vectorImpl->operator+( *v.M_vectorImpl ) ); }
        self_type operator+( value_type a ) const { return self_type( M_vectorImpl->operator+( a ) ); }
        friend self_type operator+( value_type a, const self_type & v ) { return self_type( a + (*v.M_vectorImpl) ); }
        self_type operator-( const self_type & v ) const { return self_type( M_vectorImpl->operator-( *v.M_vectorImpl ) ); }
        self_type operator-( value_type a ) const { return self_type( M_vectorImpl->operator-( a ) ); }
        friend self_type operator-( value_type a, const self_type & v ) { return self_type( a - (*v.M_vectorImpl) ); }
        self_type operator-() const { return self_type( M_vectorImpl->operator-() ); }
        self_type operator*( value_type a ) const { return self_type( M_vectorImpl->operator*( a ) ); }
        friend self_type operator*( value_type a, const self_type & v ) { return self_type( a * (*v.M_vectorImpl) ); }

        // Iterators API
        auto begin() { return M_vectorImpl->begin(); }
        auto begin() const { return M_vectorImpl->begin(); }
        auto end() { return M_vectorImpl->end(); }
        auto end() const { return M_vectorImpl->end(); }
        
        auto beginActive() { return M_vectorImpl->beginActive(); }
        auto beginActive() const { return M_vectorImpl->beginActive(); }
        auto endActive() { return M_vectorImpl->endActive(); }
        auto endActive() const { return M_vectorImpl->endActive(); }

        auto beginGhost() { return M_vectorImpl->beginGhost(); }
        auto beginGhost() const { return M_vectorImpl->beginGhost(); }
        auto endGhost() { return M_vectorImpl->endGhost(); }
        auto endGhost() const { return M_vectorImpl->endGhost(); }

        //size_type rowStart() const { return M_vectorImpl->rowStart(); }
        //size_type rowStop() const { return M_vectorImpl->rowStop(); }

        size_type startActive() const { return M_vectorImpl->startActive(); }
        size_type startGhost() const { return M_vectorImpl->startGhost(); }
        size_type start() const { return this->startActive(); }

        // Setters API
        void setConstant( value_type a ) override { return M_vectorImpl->setConstant( a ); }
        void setZero() override { return M_vectorImpl->setZero(); }
        void zero() override { return this->setZero(); }
        void zero( size_type /*start*/, size_type /*stop*/ ) override { CHECK(false) << "unsupported"; }

        void set( const size_type i, const value_type & value ) override { return M_vectorImpl->set( i, value ); }
        void setVector( int * i, int n, value_type * v ) override { return M_vectorImpl->setVector( i, n, v ); }

        void add( const size_type i, const value_type & value ) override { return M_vectorImpl->add( i, value ); }
        void add( const value_type & value ) override { return M_vectorImpl->add( value ); }
        void add( const Vector<T> & v ) override { return M_vectorImpl->add( v ); }
        void add( const value_type & a, const Vector<T> & v ) override { return M_vectorImpl->add( a, v ); }
        void add( const eigen_vector_type<Eigen::Dynamic, value_type>& a, const std::vector<vector_ptrtype>& v ) override { return M_vectorImpl->add( a, v ); }
        void addVector( int * i, int n, value_type * v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override { return M_vectorImpl->addVector( i, n, v, K, K2 ); }
        void addVector( const std::vector<value_type> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->addVector( v, dof_ids ); }
        void addVector( const Vector<T> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->addVector( v, dof_ids ); }
        void addVector( const ublas::vector<value_type> & v, const std::vector<size_type> & dof_ids ) { return M_vectorImpl->addVector( v, dof_ids ); }
        void addVector( const Vector<T> & v, const MatrixSparse<value_type> & A ) override { return M_vectorImpl->addVector( v, A ); }
        
        void sub( const Vector<T> & v ) { return M_vectorImpl->sub( v ); }
        void sub( const value_type & a, const Vector<T> & v ) { return M_vectorImpl->sub( a, v ); }

        void insert( const std::vector<value_type> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->insert( v, dof_ids ); }
        void insert( const Vector<value_type> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->insert( v, dof_ids ); }
        void insert( const ublas::vector<value_type> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->insert( v, dof_ids ); }

        void scale( const value_type factor ) override { return M_vectorImpl->scale( factor ); }

        // Utilities
        real_type min() const override { return this->min( true ); }
        real_type min( bool parallel ) const { return M_vectorImpl->min( parallel ); }
        real_type max() const override { return this->max( true ); }
        real_type max( bool parallel ) const { return M_vectorImpl->max( parallel ); }

        real_type l1Norm() const override { return M_vectorImpl->l1Norm(); }
        real_type l2Norm() const override { return M_vectorImpl->l2Norm(); }
        real_type linftyNorm() const override { return M_vectorImpl->linftyNorm(); }

        value_type sum() const override { return M_vectorImpl->sum(); }

        value_type dot( const Vector<T> & v ) const override { return M_vectorImpl->dot( v ); }
        eigen_vector_type<Eigen::Dynamic,value_type> mDot( const std::vector<std::shared_ptr<Vector<T>>> & vs ) const override { return M_vectorImpl->mDot( vs ); }

        self_type sqrt() const { return self_type( M_vectorImpl->sqrt() ); }
        self_type pow( int n ) const { return self_type( M_vectorImpl->pow( n ) ); }
        
        // Exports
        void printMatlab( const std::string filename = "NULL", bool renumber = false ) const override { return M_vectorImpl->printMatlab( filename, renumber ); }
#ifdef FEELPP_HAS_HDF5
        FEELPP_DONT_INLINE
        void saveHDF5( const std::string & filename, const std::string & tableName = "element", bool appendMode = false ) const { return M_vectorImpl->saveHDF5( filename, tableName, appendMode ); }
        FEELPP_DONT_INLINE
        void loadHDF5( const std::string & filename, const std::string & tableName = "element" ) { return M_vectorImpl->loadHDF5( filename, tableName ); }
        FEELPP_DONT_INLINE
        void loadHDF5( const std::string & filename, std::optional<std::vector<index_type>> const& mappingFromInput, const std::string & tableName = "element" ) { return M_vectorImpl->loadHDF5( filename, tableName, mappingFromInput ); }
#endif

        // Range and slice API
        self_type range( const range_type & rangeActive, const range_type & rangeGhost, const datamap_ptrtype & dm ) { return self_type( *this, rangeActive, rangeGhost, dm ); }
        self_type slice( const slice_type & sliceActive, const slice_type & sliceGhost, const datamap_ptrtype & dm ) { return self_type( *this, sliceActive, sliceGhost, dm ); }

        //! VectorUblas view on another Vector (eg a VectorPetsc)
        static self_type createView( Vector<T> & vec );

        // Localization (parallel global to one proc local)
        void localizeToOneProcessor( ublas::vector<T> & v_local, const size_type proc_id = 0 ) const { return M_vectorImpl->localizeToOneProcessor( v_local, proc_id ); }
        void localizeToOneProcessor( std::vector<T> & v_local, const size_type proc_id = 0 ) const { return M_vectorImpl->localizeToOneProcessor( v_local, proc_id ); }

    protected:
        VectorUblas( vector_impl_ptrtype && vectorImpl ): super_type( vectorImpl->mapPtr() ), M_vectorImpl( std::move( vectorImpl ) ) { }

        friend self_type element_product<>( const self_type &, const self_type & );
        
    private:
        vector_impl_ptrtype M_vectorImpl;
};

namespace detail 
{

template< typename T >
class FEELPP_EXPORT VectorUblasBase: public Vector<T>
{
    public:
        // Typedefs
        typedef Vector<T> super_type;
        typedef VectorUblasBase<T> self_type;

        using clone_ptrtype = typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        // Actual storage variants
        typedef ublas::vector<value_type> vector_storage_type;
        typedef ublas::vector_range<vector_storage_type> vector_range_storage_type;
        typedef ublas::vector_slice<vector_storage_type> vector_slice_storage_type;
        typedef ublas::vector<value_type, Feel::detail::shallow_array_adaptor<value_type>> vector_map_storage_type;
        typedef ublas::vector_range<vector_map_storage_type> vector_range_map_storage_type;
        typedef ublas::vector_slice<vector_map_storage_type> vector_slice_map_storage_type;
        typedef std::variant<
            vector_storage_type,
            vector_range_storage_type,
            vector_slice_storage_type,
            vector_map_storage_type,
            vector_range_map_storage_type,
            vector_slice_map_storage_type
                > vector_variant_type;
        typedef std::variant<
            vector_storage_type *,
            vector_range_storage_type *,
            vector_slice_storage_type *,
            vector_map_storage_type *,
            vector_range_map_storage_type *,
            vector_slice_map_storage_type *
                > vector_ptr_variant_type;
        typedef std::variant<
            const vector_storage_type *,
            const vector_range_storage_type *,
            const vector_slice_storage_type *,
            const vector_map_storage_type *,
            const vector_range_map_storage_type *,
            const vector_slice_map_storage_type *
                > vector_cstptr_variant_type;

        // Ublas range and slice types
        typedef ublas::basic_range<typename vector_storage_type::size_type, typename vector_storage_type::difference_type> range_type;
        typedef ublas::basic_slice<typename vector_storage_type::size_type, typename vector_storage_type::difference_type> slice_type;

        // Iterator class
        template< typename ValueType >
        class iterator
        {
            public:
                using iterator_category = std::random_access_iterator_tag;
                using value_type = std::decay_t<ValueType>;
                using reference = ValueType &;
                using difference_type = std::ptrdiff_t;
                using pointer = ValueType *;

            public:
                template< typename It >
                iterator( const It & it ): M_iteratorImpl( new iterator_impl<It>( it ) ) { }
                iterator() = delete;
                
                iterator( const iterator & other ): M_iteratorImpl( other.M_iteratorImpl ? other.M_iteratorImpl->clone() : nullptr ) { }
                iterator & swap( iterator & other ) { std::swap( M_iteratorImpl, other.M_iteratorImpl ); return *this; }
                ~iterator() { if( M_iteratorImpl ) delete M_iteratorImpl; }

                iterator & operator=( const iterator & other ) { swap( iterator( other ) ); return *this; }

                reference operator*() const { return M_iteratorImpl->current(); }

                iterator & operator++() { M_iteratorImpl->next(); return *this; }
                iterator operator++(int) { iterator old = *this; operator++(); return old; }

                iterator & operator--() { M_iteratorImpl->previous(); return *this; }
                iterator operator--(int) { iterator old = *this; operator--(); return old; }

                iterator & operator+=( difference_type n ) { M_iteratorImpl->incr( n ); return *this; }
                iterator operator+( difference_type n ) { iterator tmp = *this; return tmp += n; }

                iterator & operator-=( difference_type n ) { M_iteratorImpl->decr( n ); return *this; }
                iterator operator-( difference_type n ) { iterator tmp = *this; return tmp -= n; }

                difference_type operator-( const iterator & other ) const { return M_iteratorImpl->sub( *other.M_iteratorImpl ); }

                bool operator==( const iterator & other ) const { return M_iteratorImpl->equal( *other.M_iteratorImpl ); }
                bool operator!=( const iterator & other ) const { return !( *this == other ); }
                bool operator<( const iterator & other ) const { return M_iteratorImpl->cmp( *other.M_iteratorImpl ); }
                bool operator>( const iterator & other ) const { return other < *this; }
                bool operator<=( const iterator & other ) const { return !( other < *this ); }
                bool operator>=( const iterator & other ) const { return !( *this < other ); }

            private:
                class iterator_impl_base
                {
                    public:
                        virtual ~iterator_impl_base() = default;

                        virtual iterator_impl_base * clone() const = 0;

                        virtual void assign( const iterator_impl_base & other ) = 0;

                        virtual reference current() const = 0;
                        virtual void next() = 0;
                        virtual void incr( difference_type n ) = 0;
                        virtual void previous() = 0;
                        virtual void decr( difference_type n ) = 0;
                        virtual difference_type sub( const iterator_impl_base & other ) const = 0;
                        virtual bool equal( const iterator_impl_base & other ) const = 0;
                        virtual bool cmp( const iterator_impl_base & other ) const = 0;
                };
                
                template< typename It >
                class iterator_impl: public iterator_impl_base
                {
                    public:
                        iterator_impl( const It & it ) { M_it = it; }
                        ~iterator_impl() override = default;
                        iterator_impl<It> * clone() const override { return new iterator_impl<It>( M_it ); }

                        void assign( const iterator_impl_base & other ) override { M_it = dynamic_cast<const iterator_impl<It>&>( other ).M_it; }

                        reference current() const override { return *M_it; }
                        void next() override { ++M_it; }
                        void incr( difference_type n ) override { M_it += n; }
                        void previous() override { --M_it; }
                        void decr( difference_type n ) override { M_it -= n; }
                        difference_type sub( const iterator_impl_base & other ) const override { return M_it.operator-( dynamic_cast<const iterator_impl<It>&>( other ).M_it ); }
                        bool equal( const iterator_impl_base & other ) const override { return M_it.operator==( dynamic_cast<const iterator_impl<It>&>( other ).M_it ); }
                        bool cmp( const iterator_impl_base & other ) const override { return M_it.operator<( dynamic_cast<const iterator_impl<It>&>( other ).M_it ); }

                    private:
                        It M_it;
                };

            private:
                iterator_impl_base * M_iteratorImpl;
        };

        typedef iterator<value_type> iterator_type;
        typedef iterator<const value_type> const_iterator_type;

    public:
        // Constructors/Destructor
        VectorUblasBase( ) = default;
        VectorUblasBase( const VectorUblasBase<T> & v ): super_type(v) { }
        VectorUblasBase( size_type s );
        VectorUblasBase( const datamap_ptrtype & dm );
        VectorUblasBase( size_type s, size_type n_local );

        ~VectorUblasBase() override = default;

        Vector<T> & operator=( const Vector<T> & v ) override;
        VectorUblasBase<T> & operator=( const VectorUblasBase<T> & v );

        virtual clone_ptrtype clone() const override { return clone_ptrtype( this->clonePtr() ); }
        virtual self_type * clonePtr() const = 0;
        virtual self_type * emptyPtr() const = 0;

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override;
        void init( const size_type n, const bool fast = false ) override;
        void init( const datamap_ptrtype & dm ) override;
        
        virtual void resize( size_type n ) = 0;
        virtual void clear() override = 0;

        // Status API
        bool isInitialized() const override { return true; }
        void close() override { }
        bool closed() const override { return true; }
        bool areGlobalValuesUpdated() const { return true; }
        void updateGlobalValues() const { }
        void outdateGlobalValues() { }

        // Operators API
        virtual value_type operator()( size_type i ) const override = 0;
        virtual value_type& operator()( size_type i ) override = 0;
        value_type operator[]( size_type i ) const { return this->operator()( i ); }
        value_type& operator[]( size_type i ) { return this->operator()( i ); }

        Vector<T>& operator+=( const Vector<T>& v ) override { this->add( v ); return *this; }
        //VectorUblasExpression<T> operator+( const Vector<T>& v ) const;
        Vector<T>& operator-=( const Vector<T>& v ) override { this->sub( v ); return *this; }
        Vector<T>& operator*=( const value_type & a ) { this->scale( a ); return *this; }

        std::unique_ptr<VectorUblasBase<T>> operator+( const self_type & v ) const;
        std::unique_ptr<VectorUblasBase<T>> operator+( value_type a ) const;
        friend std::unique_ptr<VectorUblasBase<T>> operator+( value_type a, const self_type & v ) { std::unique_ptr<VectorUblasBase<T>> res( v.emptyPtr() ); res->init( v.mapPtr() ); res->setConstant( a ); res->addVector( v ); return res; }
        std::unique_ptr<VectorUblasBase<T>> operator-( const self_type & v ) const;
        std::unique_ptr<VectorUblasBase<T>> operator-( value_type a ) const;
        friend std::unique_ptr<VectorUblasBase<T>> operator-( value_type a, const self_type & v ) { std::unique_ptr<VectorUblasBase<T>> res( v.emptyPtr() ); res->init( v.mapPtr() ); res->setConstant( a ); res->subVector( v ); return res; }
        std::unique_ptr<VectorUblasBase<T>> operator-() const;
        std::unique_ptr<VectorUblasBase<T>> operator*( value_type a ) const;
        friend std::unique_ptr<VectorUblasBase<T>> operator*( value_type a, const self_type & v ) { std::unique_ptr<VectorUblasBase<T>> res( v.clonePtr() ); res->scale( a ); return res; }

        // Iterators API
        virtual iterator_type begin() = 0;
        virtual const_iterator_type begin() const = 0;
        virtual iterator_type end() = 0;
        virtual const_iterator_type end() const = 0;
        
        virtual iterator_type beginActive() = 0;
        virtual const_iterator_type beginActive() const = 0;
        virtual iterator_type endActive() = 0;
        virtual const_iterator_type endActive() const = 0;

        virtual iterator_type beginGhost() = 0;
        virtual const_iterator_type beginGhost() const = 0;
        virtual iterator_type endGhost() = 0;
        virtual const_iterator_type endGhost() const = 0;
        
        virtual size_type startActive() const = 0;
        virtual size_type startGhost() const = 0;

        //virtual size_type start() const = 0;
        //virtual size_type startNonContiguousGhosts() const = 0;

        //size_type rowStart() const { checkInvariants(); return 0; }
        //size_type rowStop() const { checkInvariants(); return 0; }

        // Setters API
        virtual void setConstant( value_type v ) override = 0;
        virtual void setZero() override = 0;
        void zero() override { this->setZero(); }
        void zero( size_type /*start*/, size_type /*stop*/ ) override { CHECK(false) << "unsupported"; }

        void set( const size_type i, const value_type & value ) override;
        void set( const Vector<T> & v );
        void setVector( int * i, int n, value_type * v ) override;

        void add( const size_type i, const value_type & value ) override;
        void add( const value_type & value ) override;
        void add( const Vector<T> & v ) override;
        void add( const value_type & a, const Vector<T> & v ) override;
        void add( const eigen_vector_type<Eigen::Dynamic,value_type> & a, const std::vector<vector_ptrtype> & v ) override;
        void addVector( int * i, int n, value_type * v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override;
        void addVector( const std::vector<value_type> & v, const std::vector<size_type> & dof_ids ) override;
        void addVector( const Vector<T> & v, const std::vector<size_type> & dof_ids ) override;
        void addVector( const ublas::vector<value_type> & v, const std::vector<size_type> & dof_ids );
        void addVector( const Vector<T> & /*v*/, const MatrixSparse<value_type> & /*A*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }

        void sub( const Vector<T> & v );
        void sub( const value_type & a, const Vector<T> & v );
        virtual void sub( const value_type & a );
        
        virtual void scale( const value_type factor ) override = 0;

        void insert( const std::vector<value_type> & /*v*/, const std::vector<size_type> & /*dof_ids*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }
        void insert( const Vector<value_type> & /*v*/, const std::vector<size_type> & /*dof_ids*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }
        void insert( const ublas::vector<value_type> & /*v*/, const std::vector<size_type> & /*dof_ids*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }

        // Multiple dispatch operations
        void setVector( const VectorUblasBase<T> & v ) { return v.applySetVector( *this ); }
        virtual void setVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        void addVector( const VectorUblasBase<T> & v ) { return v.applyAddVector( *this ); }
        virtual void addVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        void maddVector( const value_type & a, const VectorUblasBase<T> & v ) { return v.applyMaddVector( a, *this ); }
        virtual void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        void subVector( const VectorUblasBase<T> & v ) { return v.applySubVector( *this ); }
        virtual void subVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        virtual void msubVector( const value_type & a, const VectorUblasBase<T> & v ) { return v.applyMsubVector( a, *this ); }
        virtual void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        void mulVector( const VectorUblasBase<T> & v ) { return v.applyMulVector( *this ); }
        virtual void mulVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void mulVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        value_type dotVector( const VectorUblasBase<T> & v ) const { return v.applyDotVector( *this ); }
        virtual value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const = 0;
        virtual value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const = 0;

#if FEELPP_HAS_PETSC
        virtual void setVector( const VectorPetsc<T> & v ) = 0;
        virtual void addVector( const VectorPetsc<T> & v ) = 0;
        virtual void maddVector( const value_type & a, const VectorPetsc<T> & v ) = 0;
        virtual void subVector( const VectorPetsc<T> & v ) = 0;
        virtual void msubVector( const value_type & a, const VectorPetsc<T> & v ) = 0;
        virtual value_type dotVector( const VectorPetsc<T> & v ) const = 0;
#endif

        // Utilities
        real_type min() const override { return this->min( true ); }
        virtual real_type min( bool parallel ) const = 0;
        real_type max() const override { return this->max( true ); }
        virtual real_type max( bool parallel ) const = 0;

        virtual real_type l1Norm() const override = 0;
        virtual real_type l2Norm() const override = 0;
        virtual real_type linftyNorm() const override = 0;

        virtual value_type sum() const override = 0;

        value_type dot( const Vector<T> & v ) const override;
        eigen_vector_type<Eigen::Dynamic,value_type> mDot( const std::vector<std::shared_ptr<Vector<T>>> & v ) const override;

        virtual std::unique_ptr<VectorUblasBase<T>> sqrt() const;
        virtual std::unique_ptr<VectorUblasBase<T>> pow( int n ) const;
        
        // Exports
        void printMatlab( const std::string filename = "NULL", bool renumber = false ) const override;
#ifdef FEELPP_HAS_HDF5
        void saveHDF5( const std::string & filename, const std::string & tableName = "element", bool appendMode = false ) const;
        void loadHDF5( const std::string & filename, const std::string & tableName = "element", std::optional<std::vector<index_type>> const& mappingFromInput = {} );
#endif

        // Localization (parallel global to one proc local)
        void localizeToOneProcessor( ublas::vector<T> & v_local, const size_type proc_id = 0 ) const;
        void localizeToOneProcessor( std::vector<T> & v_local, const size_type proc_id = 0 ) const;
        
        // Range and slice API
        std::unique_ptr<VectorUblasBase<T>> range( const range_type & rangeActive, const range_type & rangeGhost ) { return std::unique_ptr<VectorUblasBase<T>>( this->rangeImpl( rangeActive, rangeGhost ) ); }
        std::unique_ptr<VectorUblasBase<T>> slice( const slice_type & sliceActive, const slice_type & sliceGhost ) { return std::unique_ptr<VectorUblasBase<T>>( this->sliceImpl( sliceActive, sliceGhost ) ); }

        // Views
#if FEELPP_HAS_PETSC
        std::unique_ptr<VectorPetsc<T>> vectorPetsc() const { return std::unique_ptr<VectorPetsc<T>>( this->vectorPetscImpl() ); }
#endif

    protected:
        virtual void checkInvariants() const = 0;

        virtual void applySetVector( VectorUblasBase<T> & v ) const = 0;
        virtual void applyAddVector( VectorUblasBase<T> & v ) const = 0;
        virtual void applyMaddVector( const value_type & a, VectorUblasBase<T> & v ) const = 0;
        virtual void applySubVector( VectorUblasBase<T> & v ) const = 0;
        virtual void applyMsubVector( const value_type & a, VectorUblasBase<T> & v ) const = 0;
        virtual void applyMulVector( VectorUblasBase<T> & v ) const = 0;
        virtual value_type applyDotVector( const VectorUblasBase<T> & v ) const = 0;

        virtual VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) = 0;
        virtual VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) = 0;
#if FEELPP_HAS_PETSC
        virtual VectorPetsc<T> * vectorPetscImpl() const = 0;
#endif
};

template< template < typename > class V, typename T >
class SettableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applySetVector( VectorUblasBase<T> & b ) const override { return b.setVector( static_cast< const V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class AddableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applyAddVector( VectorUblasBase<T> & b ) const override { return b.addVector( static_cast< const V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class MaddableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        using typename VectorUblasBase<T>::value_type;
        void applyMaddVector( const value_type & a, VectorUblasBase<T> & b ) const override { return b.maddVector( a, static_cast< const V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class SubtractableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applySubVector( VectorUblasBase<T> & b ) const override { return b.subVector( static_cast< const V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class MsubtractableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        using typename VectorUblasBase<T>::value_type;
        void applyMsubVector( const value_type & a, VectorUblasBase<T> & b ) const override { return b.msubVector( a, static_cast< const V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class MultipliableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applyMulVector( VectorUblasBase<T> & b ) const override { return b.mulVector( static_cast< const V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class DottableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        typename VectorUblasBase<T>::value_type applyDotVector( const VectorUblasBase<T> & b ) const override { return b.dotVector( static_cast< const V<T>& >( *this ) ); }
};

/****************************************************************************/
template< typename T >
class VectorUblasContiguousGhostsBase: 
    public SettableVectorUblas<VectorUblasContiguousGhostsBase, T>,
    public AddableVectorUblas<VectorUblasContiguousGhostsBase, T>,
    public MaddableVectorUblas<VectorUblasContiguousGhostsBase, T>,
    public SubtractableVectorUblas<VectorUblasContiguousGhostsBase, T>,
    public MsubtractableVectorUblas<VectorUblasContiguousGhostsBase, T>,
    public MultipliableVectorUblas<VectorUblasContiguousGhostsBase, T>,
    public DottableVectorUblas<VectorUblasContiguousGhostsBase, T>
{
    public:
        // Typedefs
        typedef VectorUblasBase<T> super_type;

        using typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

    public:
        // Constructors/Destructor
        VectorUblasContiguousGhostsBase( ) = default;
        ~VectorUblasContiguousGhostsBase() override = default;
        using super_type::operator=;
        
        // Storage API
        virtual void resize( size_type n ) override = 0;
        virtual void clear() override = 0;

        virtual vector_cstptr_variant_type vec() const = 0;

        // Operators API
        virtual value_type operator()( size_type i ) const override = 0;
        virtual value_type& operator()( size_type i ) override = 0;

        // Setters API
        virtual void setConstant( value_type v ) override = 0;
        virtual void setZero() override = 0;

        virtual void scale( const value_type factor ) override = 0;

        // Multiple dispatch operations
        void setVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void addVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;
        
        void subVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void mulVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void mulVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const override = 0;
        value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const override = 0;
        
        // Utilities
        virtual real_type min( bool parallel ) const override = 0;
        virtual real_type max( bool parallel ) const override = 0;

        virtual real_type l1Norm() const override = 0;
        virtual real_type l2Norm() const override = 0;
        virtual real_type linftyNorm() const override = 0;

        virtual value_type sum() const override = 0;

    protected:
        virtual vector_ptr_variant_type vec() = 0;

        void checkInvariants() const override = 0;

        VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) override = 0;
        VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) override = 0;
};

template< typename T, typename Storage = ublas::vector<T> >
class VectorUblasContiguousGhosts: public VectorUblasContiguousGhostsBase<T>
{
    public:
        // Typedefs
        typedef VectorUblasContiguousGhostsBase<T> super_type;
        typedef VectorUblasBase<T> base_type;
        typedef VectorUblasContiguousGhosts<T, Storage> self_type;

        using typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        typedef Storage storage_type;

        using typename super_type::vector_storage_type;
        using typename super_type::vector_range_storage_type;
        using typename super_type::vector_slice_storage_type;
        using typename super_type::vector_map_storage_type;
        using typename super_type::vector_range_map_storage_type;
        using typename super_type::vector_slice_map_storage_type;
        // VectorUblasContiguousGhosts can only hold simple vector or vector view, no proxy
        static_assert(
                std::is_same_v< storage_type, vector_storage_type > ||
                //std::is_same_v< storage_type, vector_range_storage_type > ||
                //std::is_same_v< storage_type, vector_slice_storage_type > || 
                std::is_same_v< storage_type, vector_map_storage_type >,
                //std::is_same_v< storage_type, vector_range_map_storage_type > ||
                //std::is_same_v< storage_type, vector_slice_map_storage_type >,
                "unsupported storage type" );
        static constexpr bool is_vector_slice = false;
        static constexpr bool is_vector_range = false;
        static constexpr bool is_vector_proxy = false;

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

        using typename base_type::iterator_type;
        using typename base_type::const_iterator_type;

        friend VectorUblasRange<T, Storage>;
        friend VectorUblasSlice<T, Storage>;

    public:
        VectorUblasContiguousGhosts( ): base_type() { }
        VectorUblasContiguousGhosts( size_type s ): base_type(s) { }
        VectorUblasContiguousGhosts( const datamap_ptrtype & dm ): base_type( dm ) { }
        VectorUblasContiguousGhosts( size_type s, size_type n_local ): base_type( s, n_local ) { }

        VectorUblasContiguousGhosts( const Storage & s, const datamap_ptrtype & dm ): base_type( dm ), M_vec( s ) { }
        VectorUblasContiguousGhosts( Storage && s, const datamap_ptrtype & dm ): base_type( dm ), M_vec( std::move(s) ) { }
        VectorUblasContiguousGhosts( VectorUblasContiguousGhosts<T> const& v ): base_type(v), M_vec( v.M_vec ) { }

        ~VectorUblasContiguousGhosts() override = default;

        using base_type::operator=;

        using base_type::init;
        void init( const size_type n, const size_type n_local, const bool fast = false ) override;

        virtual self_type * clonePtr() const override;
        virtual base_type * emptyPtr() const override;

        // Storage API
        virtual void resize( size_type n ) override;
        virtual void clear() override;

        virtual vector_cstptr_variant_type vec() const override { return &M_vec; }

        // Operators API
        virtual value_type operator()( size_type i ) const override;
        virtual value_type& operator()( size_type i ) override;
        
        // Iterators API
        virtual iterator_type begin() override { return iterator_type( M_vec.begin() ); }
        virtual const_iterator_type begin() const override { return const_iterator_type( M_vec.begin() ); }
        virtual iterator_type end() override { return iterator_type( M_vec.end() ); }
        virtual const_iterator_type end() const override { return const_iterator_type( M_vec.end() ); }
        
        virtual iterator_type beginActive() override { return iterator_type( M_vec.begin() ); }
        virtual const_iterator_type beginActive() const override { return const_iterator_type( M_vec.begin() ); }
        virtual iterator_type endActive() override { return iterator_type( M_vec.find( this->map().nLocalDofWithoutGhost() ) ); }
        virtual const_iterator_type endActive() const override { return const_iterator_type( M_vec.find( this->map().nLocalDofWithoutGhost() ) ); }

        virtual iterator_type beginGhost() override { return iterator_type( M_vec.find( this->map().nLocalDofWithoutGhost() ) ); }
        virtual const_iterator_type beginGhost() const override { return const_iterator_type( M_vec.find( this->map().nLocalDofWithoutGhost() ) ); }
        virtual iterator_type endGhost() override { return iterator_type( M_vec.end() ); }
        virtual const_iterator_type endGhost() const override { return const_iterator_type( M_vec.end() ); }
        
        size_type startActive() const override;
        size_type startGhost() const override;

        // Setters API
        virtual void setConstant( value_type v ) override;
        virtual void setZero() override;

        virtual void add( const value_type & a ) override;
        virtual void sub( const value_type & a ) override;
        virtual void scale( const value_type factor ) override;

        // Multiple dispatch operations
        using base_type::setVector;
        void setVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::addVector;
        void addVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::maddVector;
        void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;
        
        using base_type::subVector;
        void subVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::msubVector;
        void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::mulVector;
        void mulVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void mulVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::dotVector;
        value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const override;
        value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const override;

#if FEELPP_HAS_PETSC
        virtual void setVector( const VectorPetsc<T> & v ) override;
        virtual void addVector( const VectorPetsc<T> & v ) override;
        virtual void maddVector( const value_type & a, const VectorPetsc<T> & v ) override;
        virtual void subVector( const VectorPetsc<T> & v ) override;
        virtual void msubVector( const value_type & a, const VectorPetsc<T> & v ) override;
        virtual value_type dotVector( const VectorPetsc<T> & v ) const override;
#endif
        
        // Utilities
        virtual real_type min( bool parallel ) const override;
        virtual real_type max( bool parallel ) const override;

        virtual real_type l1Norm() const override;
        virtual real_type l2Norm() const override;
        virtual real_type linftyNorm() const override;

        virtual value_type sum() const override;

        virtual std::unique_ptr<VectorUblasBase<T>> sqrt() const override;
        virtual std::unique_ptr<VectorUblasBase<T>> pow( int n ) const override;

    protected:
        virtual vector_ptr_variant_type vec() override { return &M_vec; }
        
        void checkInvariants() const override;
         
        VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) override;
        VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) override;
#if FEELPP_HAS_PETSC
        VectorPetsc<T> * vectorPetscImpl() const override;
#endif

    protected:
        storage_type M_vec;

};

/****************************************************************************/
template< typename T >
class VectorUblasNonContiguousGhostsBase: 
    public SettableVectorUblas<VectorUblasNonContiguousGhostsBase, T>,
    public AddableVectorUblas<VectorUblasNonContiguousGhostsBase, T>,
    public MaddableVectorUblas<VectorUblasNonContiguousGhostsBase, T>,
    public SubtractableVectorUblas<VectorUblasNonContiguousGhostsBase, T>,
    public MsubtractableVectorUblas<VectorUblasNonContiguousGhostsBase, T>,
    public MultipliableVectorUblas<VectorUblasNonContiguousGhostsBase, T>,
    public DottableVectorUblas<VectorUblasNonContiguousGhostsBase, T>
{
    public:
        // Typedefs
        typedef VectorUblasBase<T> super_type;

        using typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

    public:
        // Constructors/Destructor
        VectorUblasNonContiguousGhostsBase( ) = default;
        ~VectorUblasNonContiguousGhostsBase() override = default;
        using super_type::operator=;
        
        // Storage API
        virtual void resize( size_type n ) override = 0;
        virtual void clear() override = 0;

        virtual vector_cstptr_variant_type vec() const = 0;
        virtual vector_cstptr_variant_type vecNonContiguousGhosts() const = 0;

        // Operators API
        virtual value_type operator()( size_type i ) const override = 0;
        virtual value_type& operator()( size_type i ) override = 0;

        // Setters API
        virtual void setConstant( value_type v ) override = 0;
        virtual void setZero() override = 0;

        virtual void scale( const value_type factor ) override = 0;

        // Multiple dispatch operations
        void setVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void addVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;
        
        void subVector( const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override = 0;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override = 0;

        value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const override = 0;
        value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const override = 0;

        // Utilities
        virtual real_type min( bool parallel ) const override = 0;
        virtual real_type max( bool parallel ) const override = 0;

        virtual real_type l1Norm() const override = 0;
        virtual real_type l2Norm() const override = 0;
        virtual real_type linftyNorm() const override = 0;

        virtual value_type sum() const override = 0;

    protected:
        virtual vector_ptr_variant_type vec() = 0;
        virtual vector_ptr_variant_type vecNonContiguousGhosts() = 0;
        
        void checkInvariants() const override = 0;
 
        VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) override = 0;
        VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) override = 0;

};

template< typename T, typename Storage = ublas::vector<T> >
class VectorUblasNonContiguousGhosts: public VectorUblasNonContiguousGhostsBase<T>
{
    public:
        // Typedefs
        typedef VectorUblasNonContiguousGhostsBase<T> super_type;
        typedef VectorUblasBase<T> base_type;
        typedef VectorUblasNonContiguousGhosts<T, Storage> self_type;

        using typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        typedef Storage storage_type;
        
        using typename super_type::vector_storage_type;
        using typename super_type::vector_range_storage_type;
        using typename super_type::vector_slice_storage_type;
        using typename super_type::vector_map_storage_type;
        using typename super_type::vector_range_map_storage_type;
        using typename super_type::vector_slice_map_storage_type;
        static_assert( 
                std::is_same_v< storage_type, vector_storage_type > ||
                std::is_same_v< storage_type, vector_range_storage_type > ||
                std::is_same_v< storage_type, vector_slice_storage_type > || 
                std::is_same_v< storage_type, vector_map_storage_type > ||
                std::is_same_v< storage_type, vector_range_map_storage_type > ||
                std::is_same_v< storage_type, vector_slice_map_storage_type >,
                "unsupported storage type" );
        static constexpr bool is_vector_range = 
            std::is_same_v< storage_type, vector_range_storage_type > ||
            std::is_same_v< storage_type, vector_range_map_storage_type >;
        static constexpr bool is_vector_slice = 
            std::is_same_v< storage_type, vector_slice_storage_type > || 
            std::is_same_v< storage_type, vector_slice_map_storage_type >;
        static constexpr bool is_vector_proxy = is_vector_range || is_vector_slice;

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

        using typename base_type::iterator_type;
        using typename base_type::const_iterator_type;

        friend VectorUblasRange<T, Storage>;
        friend VectorUblasSlice<T, Storage>;

    public:
        VectorUblasNonContiguousGhosts( ): base_type() { }
        VectorUblasNonContiguousGhosts( size_type s ): base_type(s) { }
        VectorUblasNonContiguousGhosts( const datamap_ptrtype & dm ): base_type( dm ) { }
        VectorUblasNonContiguousGhosts( size_type s, size_type n_local ): base_type( s, n_local ) { }

        VectorUblasNonContiguousGhosts( const Storage & s, const Storage & sGhost, const datamap_ptrtype & dm ): base_type( dm ), M_vec( s ), M_vecNonContiguousGhosts( sGhost ) { }
        VectorUblasNonContiguousGhosts( Storage && s, Storage && sGhost, const datamap_ptrtype & dm ): base_type( dm ), M_vec( std::move(s) ), M_vecNonContiguousGhosts( std::move(sGhost) ) { }

        VectorUblasNonContiguousGhosts( VectorUblasNonContiguousGhosts<T, Storage> const& v ): base_type(v), M_vec( v.M_vec ), M_vecNonContiguousGhosts( v.M_vecNonContiguousGhosts ) { }

        ~VectorUblasNonContiguousGhosts() override = default;

        using base_type::operator=;
        
        using base_type::init;
        void init( const size_type n, const size_type n_local, const bool fast = false ) override;

        virtual self_type * clonePtr() const override;
        virtual base_type * emptyPtr() const override;

        // Storage API
        virtual void resize( size_type n ) override;
        virtual void clear() override;

        virtual vector_cstptr_variant_type vec() const override { return &M_vec; }
        virtual vector_cstptr_variant_type vecNonContiguousGhosts() const override { return &M_vecNonContiguousGhosts; }

        // Operators API
        virtual value_type operator()( size_type i ) const override;
        virtual value_type& operator()( size_type i ) override;
        
        // Iterators API
        virtual iterator_type begin() override { return iterator_type( M_vec.begin() ); }
        virtual const_iterator_type begin() const override { return const_iterator_type( M_vec.begin() ); }
        virtual iterator_type end() override { return iterator_type( M_vec.end() ); }
        virtual const_iterator_type end() const override { return const_iterator_type( M_vec.end() ); }

        virtual iterator_type beginActive() override { return iterator_type( M_vec.begin() ); }
        virtual const_iterator_type beginActive() const override { return const_iterator_type( M_vec.begin() ); }
        virtual iterator_type endActive() override { return iterator_type( M_vec.end() ); }
        virtual const_iterator_type endActive() const override { return const_iterator_type( M_vec.end() ); }

        virtual iterator_type beginGhost() override { return iterator_type( M_vecNonContiguousGhosts.begin() ); }
        virtual const_iterator_type beginGhost() const override { return const_iterator_type( M_vecNonContiguousGhosts.begin() ); }
        virtual iterator_type endGhost() override { return iterator_type( M_vecNonContiguousGhosts.end() ); }
        virtual const_iterator_type endGhost() const override { return const_iterator_type( M_vecNonContiguousGhosts.end() ); }

        size_type startActive() const override;
        size_type startGhost() const override;

        // Setters API
        virtual void setConstant( value_type v ) override;
        virtual void setZero() override;

        virtual void add( const value_type & a ) override;
        virtual void sub( const value_type & a ) override;
        virtual void scale( const value_type factor ) override;

        // Multiple dispatch operations
        using base_type::setVector;
        void setVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::addVector;
        void addVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::maddVector;
        void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;
        
        using base_type::subVector;
        void subVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::msubVector;
        void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::mulVector;
        void mulVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void mulVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        using base_type::dotVector;
        value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const override;
        value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const override;

#if FEELPP_HAS_PETSC
        virtual void setVector( const VectorPetsc<T> & v ) override;
        virtual void addVector( const VectorPetsc<T> & v ) override;
        virtual void maddVector( const value_type & a, const VectorPetsc<T> & v ) override;
        virtual void subVector( const VectorPetsc<T> & v ) override;
        virtual void msubVector( const value_type & a, const VectorPetsc<T> & v ) override;
        virtual value_type dotVector( const VectorPetsc<T> & v ) const override;
#endif

        // Utilities
        virtual real_type min( bool parallel ) const override;
        virtual real_type max( bool parallel ) const override;

        virtual real_type l1Norm() const override;
        virtual real_type l2Norm() const override;
        virtual real_type linftyNorm() const override;

        virtual value_type sum() const override;

        virtual std::unique_ptr<VectorUblasBase<T>> sqrt() const override;
        virtual std::unique_ptr<VectorUblasBase<T>> pow( int n ) const override;

    protected:
        virtual vector_ptr_variant_type vec() override { return &M_vec; }
        virtual vector_ptr_variant_type vecNonContiguousGhosts() override { return &M_vecNonContiguousGhosts; }
        
        void checkInvariants() const override;
         
        VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) override;
        VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) override;
#if FEELPP_HAS_PETSC
        VectorPetsc<T> * vectorPetscImpl() const override;
#endif

    protected:
        storage_type M_vec;
        storage_type M_vecNonContiguousGhosts;

};

template< typename T, typename Storage >
class VectorUblasRange: public VectorUblasNonContiguousGhosts<T, ublas::vector_range<Storage> >
{
    public:
        typedef VectorUblasNonContiguousGhosts<T, ublas::vector_range<Storage>> super_type;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;
        
        typedef Storage storage_type;
        using typename super_type::vector_storage_type;
        using typename super_type::vector_map_storage_type;
        static_assert( 
                std::is_same_v< storage_type, vector_storage_type > ||
                std::is_same_v< storage_type, vector_map_storage_type >,
                "range not supported for storage type" );

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

    public:
        //VectorUblasRange( VectorUblasContiguousGhosts<T, Storage> & v, const range_type & rangeActive, const range_type & rangeGhost );
        //VectorUblasRange( VectorUblasNonContiguousGhosts<T, Storage> & v, const range_type & rangeActive, const range_type & rangeGhost );
        VectorUblasRange( Storage & v, const range_type & range, const datamap_ptrtype & dm );
        VectorUblasRange( Storage & v, const range_type & rangeActive, const range_type & rangeGhost, const datamap_ptrtype & dm );
        VectorUblasRange( Storage & vActive, const range_type & rangeActive, 
                Storage & vGhost, const range_type & rangeGhost, const datamap_ptrtype & dm );

};

template< typename T, typename Storage >
class VectorUblasSlice: public VectorUblasNonContiguousGhosts<T, ublas::vector_slice<Storage> >
{
    public:
        typedef VectorUblasNonContiguousGhosts<T, ublas::vector_slice<Storage>> super_type;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;
        
        typedef Storage storage_type;
        using typename super_type::vector_storage_type;
        using typename super_type::vector_map_storage_type;
        static_assert( 
                std::is_same_v< storage_type, vector_storage_type > ||
                std::is_same_v< storage_type, vector_map_storage_type >,
                "slice not supported for storage type" );

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

    public:
        //VectorUblasSlice( VectorUblasContiguousGhosts<T, Storage> & v, const slice_type & sliceActive, const slice_type & sliceGhost );
        //VectorUblasSlice( VectorUblasNonContiguousGhosts<T, Storage> & v, const slice_type & sliceActive, const slice_type & sliceGhost );
        VectorUblasSlice( Storage & v, const slice_type & slice, const datamap_ptrtype & dm );
        VectorUblasSlice( Storage & v, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm );
        VectorUblasSlice( Storage & vActive, slice_type const& sliceActive,
                Storage & vGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm );

};

} // detail

/**
 * FEELPP_INSTANTIATE_VECTORUBLAS is never defined except in vectorublas.cpp
 * where we do the instantiate. This allows to reduce the VectorUblas
 * instantiation to the strict minimum
 */
#if !defined( FEELPP_INSTANTIATE_VECTORUBLAS )
extern template class VectorUblas<double>;
namespace detail {
extern template class VectorUblasBase< double >;
extern template class VectorUblasContiguousGhosts< double, ublas::vector<double> >;
extern template class VectorUblasNonContiguousGhosts< double, ublas::vector<double> >;
extern template class VectorUblasRange< double, ublas::vector<double> >;
extern template class VectorUblasSlice< double, ublas::vector<double> >;
}
#endif
} // Feel

#if FEELPP_HAS_PETSC
#include <feel/feelalg/vectorpetsc.hpp>

namespace Feel {
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
 * \warning VectorPetsc view is invalidated by VectorUblas re-assignment
 */
//template< typename T >
//inline VectorPetsc<T>
//toPETSc( const VectorUblas<T> & v )
//{
    //using namespace detail;
    //if ( v.comm().size() > 1 )
    //{
        //const VectorUblasContiguousGhosts<T> * vecContiguousGhosts = dynamic_cast<const VectorUblasContiguousGhosts<T> *>( & v.vectorImpl() );
        //if ( vecContiguousGhosts )
            //return VectorPetscMPI<T>( std::addressof( * vecContiguousGhosts->begin() ), vecContiguousGhosts->mapPtr() );

        //const VectorUblasNonContiguousGhosts<T> * vecNonContiguousGhosts = dynamic_cast<const VectorUblasNonContiguousGhosts<T> *>( & v.vectorImpl() );
        //if ( vecNonContiguousGhosts )
            //return VectorPetscMPIRange<T>( std::addressof( * vecNonContiguousGhosts->beginActive() ), std::addressof( * vecNonContiguousGhosts->beginGhost() ), vecNonContiguousGhosts->mapPtr() );
        
        //else
            //CHECK( false ) << "Unsupported VectorUblas for VectorPetsc view";
    //}
    //else
        //return VectorPetsc<T>( std::addressof( * v.vectorImpl().begin() ), v.mapPtr() );
//}

template< typename T >
inline std::shared_ptr<VectorPetsc<T>>
toPETScPtr( const VectorUblas<T> & v )
{
    return std::shared_ptr<VectorPetsc<T>>( v.vectorImpl().vectorPetsc() );
}

template< typename T >
inline VectorPetsc<T>
toPETSc( const VectorUblas<T> & v )
{
    return *toPETScPtr( v );
}

} // Feel

#endif // FEELPP_HAS_PETSC

namespace Feel {
/* Utility functions */
namespace detail {

template< typename T >
std::unique_ptr<VectorUblasBase<T>> 
element_product( const VectorUblasBase<T> & v1, const VectorUblasBase<T> & v2 )
{
    VectorUblasContiguousGhosts<T> * prod = new VectorUblasContiguousGhosts<T>();
    *prod = v1;
    prod->mulVector( v2 );
    return std::unique_ptr<VectorUblasBase<T>>( prod );
}

} //detail

template< typename T >
VectorUblas<T>
element_product( const VectorUblas<T> & v1, const VectorUblas<T> & v2 )
{
    return VectorUblas<T>( Feel::detail::element_product( v1.vectorImpl(), v2.vectorImpl() ) );
}

} // Feel

#endif /* _FEELPP_VECTORUBLAS_HPP */
