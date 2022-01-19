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

    shallow_array_adaptor(size_type n) : base_type(n) {}
    shallow_array_adaptor(size_type n, pointer data) : base_type(n,data) {}
    shallow_array_adaptor(const shallow_array_adaptor& c) : base_type(c) {}

    // This function must swap the values of the items, not the data pointers.
    void swap(shallow_array_adaptor& a) {
        if (base_type::begin() != a.begin())
            std::swap_ranges(base_type::begin(), base_type::end(), a.begin());
    }
};

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

template< typename T >
class VectorUblas : public Vector<T>
{
    public:
        // Typedefs
        typedef VectorUblas<T> self_type;
        typedef Vector<T> super_type;

        using clone_ptrtype = typename super_type::clone_ptrtype;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;

        typedef detail::VectorUblasBase<value_type> vector_impl_type;
        typedef std::unique_ptr<vector_impl_type> vector_impl_ptrtype;
        typedef std::unique_ptr<const vector_impl_type> vector_impl_cstptrtype;

        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;
        
        // Ublas range and slice types
        typedef ublas::basic_range<typename ublas::vector<value_type>::size_type, typename ublas::vector<value_type>::difference_type> range_type;
        typedef ublas::basic_slice<typename ublas::vector<value_type>::size_type, typename ublas::vector<value_type>::difference_type> slice_type;

    public:
        // Constructors/Destructor
        VectorUblas();
        VectorUblas( size_type s );
        VectorUblas( datamap_ptrtype const& dm );
        VectorUblas( size_type s, size_type n_local );
        //FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, range_type const& range, datamap_ptrtype const& dm );
        //TODO:VectorUblas( VectorUblas<value_type>& v, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm );
        //VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm );
        //FEELPP_DEPRECATED VectorUblas( VectorUblas<value_type>& m, slice_type const& range, datamap_ptrtype const& dm );
        //TODO:VectorUblas( VectorUblas<value_type>& v, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm );
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

        VectorUblas( const VectorUblas & other ): M_vectorImpl( other.M_vectorImpl ? other.M_vectorImpl->clone() : nullptr ) { }
        VectorUblas & swap( VectorUblas & other ) { std::swap( M_vectorImpl, other.M_vectorImpl ); return *this; }
        ~VectorUblas() override { delete M_vectorImpl; }

        clone_ptrtype clone() const override { return clone_ptrtype( new self_type( *this ) ); }

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override { return M_vectorImpl->init( n, n_local, fast ); }
        void init( const size_type n, const bool fast = false ) override { return M_vectorImpl->init( n, fast ); }
        void init( const datamap_ptrtype & dm ) override { return M_vectorImpl->init( dm ); }
        
        void resize( size_type n ) { return M_vectorImpl->resize( n ); }
        void clear() override { return M_vectorImpl->clear(); }

        const vector_impl_type & vectorImpl() const { return *M_vectorImpl; }

        // Status API
        bool isInitialized() const override { return true; }
        void close() override { }
        bool closed() const override { return true; }
        bool areGlobalValuesUpdated() const { return true; }
        void updateGlobalValues() const { }
        void outdateGlobalValues() { }

        // Operators API
        Vector<value_type>& operator=( const Vector<value_type> & V ) override;
        VectorUblas & operator=( const VectorUblas & other ) { swap( VectorUblas( other ) ); return *this; }

        virtual value_type operator()( size_type i ) const override { return M_vectorImpl->operator()( i ); }
        virtual value_type& operator()( size_type i ) override { return M_vectorImpl->operator()( i ); }
        value_type operator[]( size_type i ) const { return this->operator()( i ); }
        value_type& operator[]( size_type i ) { return this->operator()( i ); }

        Vector<T>& operator+=( const Vector<T>& v ) override { this->add( v ); return *this; }
        //VectorUblasExpression<T> operator+( const Vector<T>& v ) const;

        //// Iterators API
        //virtual iterator begin() = 0;
        //virtual iterator end() = 0;

        //virtual iterator beginGhost() = 0;
        //virtual iterator endGhost() = 0;

        //virtual size_type start() const = 0;
        //virtual size_type startNonContiguousGhosts() const = 0;

        //size_type rowStart() const { return M_vectorImpl->rowStart(); }
        //size_type rowStop() const { return M_vectorImpl->rowStop(); }

        // Setters API
        void setConstant( value_type a ) override { return M_vectorImpl->setConstant( a ); }
        void setZero() override { return M_vectorImpl->setZero(); }
        void zero() override { return this->setZero(); }

        void set( const size_type i, const value_type & value ) override { return M_vectorImpl->set( i, value ); }
        void setVector( int * i, int n, value_type * v ) override { return M_vectorImpl->setVector( i, v ); }

        void add( const size_type i, const value_type & value ) override { return M_vectorImpl->add( i, value ); }
        void add( const value_type & value ) override { return M_vectorImpl->add( value ); }
        void add( const Vector<T> & v ) override { return M_vectorImpl->add( v ); }
        void add( const value_type & a, const Vector<T> & v ) override { return M_vectorImpl->add( a, v ); }
        void addVector( int * i, int n, value_type * v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override { return M_vectorImpl->addVector( i, n, v, K, K2 ); }
        void addVector( const std::vector<value_type> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->addVector( v, dof_ids ); }
        void addVector( const Vector<T> & v, const std::vector<size_type> & dof_ids ) override { return M_vectorImpl->addVector( v, dof_ids ); }
        void addVector( const ublas::vector<value_type> & v, const std::vector<size_type> & dof_ids ) { return M_vectorImpl->addVector( v, dof_ids ); }
        void addVector( const Vector<T> & v, const MatrixSparse<value_type> & A ) override { return M_vectorImpl->addVector( v, A ); }

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
        
        // Exports
        void printMatlab( const std::string filename = "NULL", bool renumber = false ) const override { return M_vectorImpl->printMatlab( filename, renumber ); }
#ifdef FEELPP_HAS_HDF5
        void saveHDF5( const std::string & filename, const std::string & tableName = "element", bool appendMode = false ) const { return M_vectorImpl->saveHDF5( filename, tableName, appendMode ); }
        void loadHDF5( const std::string & filename, const std::string & tableName = "element" ) { return M_vectorImpl->loadHDF5( filename, tableName ); }
#endif

        // Localization (parallel global to one proc local)
        void localizeToOneProcessor( ublas::vector<T> & v_local, const size_type proc_id = 0 ) const { return M_vectorImpl->localizeToOneProcessor( v_local, proc_id ); }
        void localizeToOneProcessor( std::vector<T> & v_local, const size_type proc_id = 0 ) const { return M_vectorImpl->localizeToOneProcessor( v_local, proc_id ); }
        
    private:
        vector_impl_ptrtype M_vectorImpl;
};

namespace detail 
{

template< typename T >
class VectorUblasBase: public Vector<T>
{
    public:
        // Typedefs
        typedef Vector<T> super_type;

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
        class iterator
        {
            public:
                template< typename It >
                iterator( const It & it ): M_iteratorImpl( new iterator_impl<It>( it ) ) { }
                iterator() = delete;
                
                iterator( const iterator & other ): M_iteratorImpl( other.M_iteratorImpl ? other.M_iteratorImpl->clone() : nullptr ) { }
                iterator & swap( iterator & other ) { std::swap( M_iteratorImpl, other.M_iteratorImpl ); return *this; }
                ~iterator() { if( M_iteratorImpl ) delete M_iteratorImpl; }

                iterator & operator=( const iterator & other ) { swap( iterator( other ) ); return *this; }
                value_type & operator*() const { return M_iteratorImpl->current(); }
                iterator & operator++() { M_iteratorImpl->next(); return *this; }
                iterator & operator--() { M_iteratorImpl->previous(); return *this; }
                bool operator==( const iterator & other ) const { return M_iteratorImpl->equal( *other.M_iteratorImpl ); }

            private:
                class iterator_impl_base
                {
                    public:
                        iterator_impl_base * clone() const = 0;

                        virtual value_type & current() const = 0;
                        virtual void next() = 0;
                        virtual void previous() = 0;
                        virtual bool equal( const iterator_impl_base & other ) const = 0;
                        virtual void assign( const iterator_impl_base & other ) = 0;
                };
                
                template< typename It >
                class iterator_impl: public iterator_impl_base
                {
                    public:
                        iterator_impl( const It & it ) { M_it = it; }
                        iterator_impl<It> * clone() const override { return new iterator_impl<It>( M_it ); }

                        value_type & current() const override { return *M_it; }
                        void next() override { ++M_it; }
                        void previous() override { --M_it; }
                        bool equal( const iterator_impl_base & other ) const { return M_it.operator==( dynamic_cast<iterator_impl<It>&>( other ).M_it ); }
                        void assign( const iterator_impl_base & other ) { M_it = dynamic_cast<iterator_impl<It>&>( other ).M_it; }

                    private:
                        It M_it;
                };

            private:
                iterator_impl_base * M_iteratorImpl;
        };

    public:
        // Constructors/Destructor
        VectorUblasBase( ) = default;
        VectorUblasBase( VectorUblasBase<T> const & v ): super_type(v) { }
        VectorUblasBase( size_type s );
        VectorUblasBase( const datamap_ptrtype & dm );
        VectorUblasBase( size_type s, size_type n_local );

        ~VectorUblasBase() override = default;

        virtual clone_ptrtype clone() const override = 0;

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override = 0;
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
        Vector<value_type>& operator=( const Vector<value_type> & v ) override { this->set( v ); return *this; }

        virtual value_type operator()( size_type i ) const override = 0;
        virtual value_type& operator()( size_type i ) override = 0;
        value_type operator[]( size_type i ) const { return this->operator()( i ); }
        value_type& operator[]( size_type i ) { return this->operator()( i ); }

        Vector<T>& operator+=( const Vector<T>& v ) override { this->add( v ); return *this; }
        //VectorUblasExpression<T> operator+( const Vector<T>& v ) const;
        Vector<T>& operator-=( const Vector<T>& v ) override { this->sub( v ); return *this; }
        Vector<T>& operator*=( const value_type & a ) { this->scale( a ); return *this; }

        // Iterators API
        virtual iterator begin() = 0;
        virtual iterator end() = 0;
        
        virtual iterator beginActive() = 0;
        virtual iterator endActive() = 0;

        virtual iterator beginGhost() = 0;
        virtual iterator endGhost() = 0;

        //virtual size_type start() const = 0;
        //virtual size_type startNonContiguousGhosts() const = 0;

        //size_type rowStart() const { checkInvariants(); return 0; }
        //size_type rowStop() const { checkInvariants(); return 0; }

        // Setters API
        virtual void setConstant( value_type v ) override = 0;
        virtual void setZero() override = 0;
        void zero() override { this->setZero(); }

        void set( const size_type i, const value_type & value ) override;
        void set( const Vector<T> & v );
        void setVector( int * i, int n, value_type * v ) override;

        void add( const size_type i, const value_type & value ) override;
        void add( const value_type & value ) override;
        void add( const Vector<T> & v ) override;
        void add( const value_type & a, const Vector<T> & v ) override;
        void addVector( int * i, int n, value_type * v, size_type K = 0, size_type K2 = invalid_v<size_type> ) override;
        void addVector( const std::vector<value_type> & v, const std::vector<size_type> & dof_ids ) override;
        void addVector( const Vector<T> & v, const std::vector<size_type> & dof_ids ) override;
        void addVector( const ublas::vector<value_type> & v, const std::vector<size_type> & dof_ids );
        void addVector( const Vector<T> & /*v*/, const MatrixSparse<value_type> & /*A*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }

        void sub( const Vector<T> & v );
        void sub( const value_type & a, const Vector<T> & v );
        
        virtual void scale( const value_type factor ) override = 0;

        void insert( const std::vector<value_type> & /*v*/, const std::vector<size_type> & /*dof_ids*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }
        void insert( const Vector<value_type> & /*v*/, const std::vector<size_type> & /*dof_ids*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }
        void insert( const ublas::vector<value_type> & /*v*/, const std::vector<size_type> & /*dof_ids*/ ) override { FEELPP_ASSERT( 0 ).error( "not implemented" ); }

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
        
        // Exports
        void printMatlab( const std::string filename = "NULL", bool renumber = false ) const override;
#ifdef FEELPP_HAS_HDF5
        void saveHDF5( const std::string & filename, const std::string & tableName = "element", bool appendMode = false ) const;
        void loadHDF5( const std::string & filename, const std::string & tableName = "element" );
#endif

        // Localization (parallel global to one proc local)
        void localizeToOneProcessor( ublas::vector<T> & v_local, const size_type proc_id = 0 ) const;
        void localizeToOneProcessor( std::vector<T> & v_local, const size_type proc_id = 0 ) const;
        
        // Range and slice API
        std::unique_ptr<VectorUblasBase<T>> range( const range_type & rangeActive, const range_type & rangeGhost ) { return std::unique_ptr<VectorUblasBase<T>>( this->rangeImpl( rangeActive, rangeGhost ) ); }
        std::unique_ptr<VectorUblasBase<T>> slice( const slice_type & sliceActive, const slice_type & sliceGhost ) { return std::unique_ptr<VectorUblasBase<T>>( this->sliceImpl( sliceActive, sliceGhost ) ); }

    protected:
        virtual void checkInvariants() const = 0;

        void setVector( const VectorUblasBase<T> & v ) { return v.applySetVector( *this ); }
        virtual void applySetVector( VectorUblasBase<T> & v ) const = 0;
        virtual void setVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;
        
        void addVector( const VectorUblasBase<T> & v ) { return v.applyAddVector( *this ); }
        virtual void applyAddVector( VectorUblasBase<T> & v ) const = 0;
        virtual void addVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        void maddVector( const value_type & a, const VectorUblasBase<T> & v ) { return v.applyMaddVector( a, *this ); }
        virtual void applyMaddVector( const value_type & a, VectorUblasBase<T> & v ) const = 0;
        virtual void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;
        
        void subVector( const VectorUblasBase<T> & v ) { return v.applySubVector( *this ); }
        virtual void applySubVector( VectorUblasBase<T> & v ) const = 0;
        virtual void subVector( const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        virtual void msubVector( const value_type & a, const VectorUblasBase<T> & v ) { return v.applyMsubVector( a, *this ); }
        virtual void applyMsubVector( const value_type & a, VectorUblasBase<T> & v ) const = 0;
        virtual void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) = 0;
        virtual void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) = 0;

        virtual value_type dotVector( const VectorUblasBase<T> & v ) const { return v.applyDotVector( *this ); }
        virtual value_type applyDotVector( const VectorUblasBase<T> & v ) const = 0;
        virtual value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const = 0;
        virtual value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const = 0;

        virtual VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) = 0;
        virtual VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) = 0;
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
        VectorUblasContiguousGhostsBase( VectorUblasContiguousGhostsBase<T> const& v ): super_type(v) { }
        
        VectorUblasContiguousGhostsBase( size_type s ): super_type( s ) { }
        VectorUblasContiguousGhostsBase( const datamap_ptrtype & dm ): super_type( dm ) { }
        VectorUblasContiguousGhostsBase( size_type s, size_type n_local ): super_type( s, n_local ) { }

        ~VectorUblasContiguousGhostsBase() override = default;
        
        void init( const size_type n, const size_type n_local, const bool fast = false ) override = 0;
        
        virtual clone_ptrtype clone() const override = 0;

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

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

        using typename base_type::iterator;

        friend VectorUblasRange<T, Storage>;
        friend VectorUblasSlice<T, Storage>;

    public:
        using super_type::VectorUblasContiguousGhostsBase;
        VectorUblasContiguousGhosts( const Storage & s, const datamap_ptrtype & dm ): super_type( dm ), M_vec( s ) { }
        VectorUblasContiguousGhosts( Storage && s, const datamap_ptrtype & dm ): super_type( dm ), M_vec( std::move(s) ) { }
        VectorUblasContiguousGhosts( VectorUblasContiguousGhosts<T> const& v ): super_type(v), M_vec( v.M_vec ) { }

        void init( const size_type n, const size_type n_local, const bool fast = false ) override;

        virtual clone_ptrtype clone() const override;

        // Storage API
        virtual void resize( size_type n ) override;
        virtual void clear() override;

        virtual vector_cstptr_variant_type vec() const override { return &M_vec; }

        // Operators API
        virtual value_type operator()( size_type i ) const override;
        virtual value_type& operator()( size_type i ) override;
        
        // Iterators API
        virtual iterator begin() override { return iterator( M_vec.begin() ); }
        virtual iterator end() override { return iterator( M_vec.end() ); }
        
        virtual iterator beginActive() override { return iterator( M_vec.begin() ); }
        virtual iterator endActive() override { return iterator( M_vec.find( this->map()->nLocalDofWithoutGhost() ) ); }

        virtual iterator beginGhost() override { return iterator( M_vec.find( this->map()->nLocalDofWithoutGhost() ) ); }
        virtual iterator endGhost() override { return iterator( M_vec.end() ); }

        // Setters API
        virtual void setConstant( value_type v ) override;
        virtual void setZero() override;

        virtual void scale( const value_type factor ) override;
        
        // Utilities
        virtual real_type min( bool parallel ) const override;
        virtual real_type max( bool parallel ) const override;

        virtual real_type l1Norm() const override;
        virtual real_type l2Norm() const override;
        virtual real_type linftyNorm() const override;

        virtual value_type sum() const override;

    protected:
        virtual vector_ptr_variant_type vec() override { return &M_vec; }
        
        void checkInvariants() const override;
        
        void setVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        void addVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;
        
        void subVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const override;
        value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const override;
        
        VectorUblasRange<T, Storage> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) override;
        VectorUblasSlice<T, Storage> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) override;

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
        VectorUblasNonContiguousGhostsBase( VectorUblasNonContiguousGhostsBase<T> const& v ): super_type(v) { }
        
        VectorUblasNonContiguousGhostsBase( size_type s ): super_type( s ) { }
        VectorUblasNonContiguousGhostsBase( const datamap_ptrtype & dm ): super_type( dm ) { }
        VectorUblasNonContiguousGhostsBase( size_type s, size_type n_local ): super_type( s, n_local ) { }

        ~VectorUblasNonContiguousGhostsBase() override = default;
        
        void init( const size_type n, const size_type n_local, const bool fast = false ) override = 0;
        
        virtual clone_ptrtype clone() const override = 0;

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

        using typename super_type::vector_variant_type;
        using typename super_type::vector_ptr_variant_type;
        using typename super_type::vector_cstptr_variant_type;
        using typename super_type::range_type;
        using typename super_type::slice_type;

        using typename base_type::iterator;

        friend VectorUblasRange<T, Storage>;
        friend VectorUblasSlice<T, Storage>;

    public:
        using super_type::VectorUblasNonContiguousGhostsBase;
        VectorUblasNonContiguousGhosts( const Storage & s, const Storage & sGhost, const datamap_ptrtype & dm ): super_type( dm ), M_vec( s ), M_vecNonContiguousGhosts( sGhost ) { }
        VectorUblasNonContiguousGhosts( Storage && s, Storage && sGhost, const datamap_ptrtype & dm ): super_type( dm ), M_vec( std::move(s) ), M_vecNonContiguousGhosts( std::move(sGhost) ) { }
        VectorUblasNonContiguousGhosts( VectorUblasNonContiguousGhosts<T, Storage> const& v ): super_type(v), M_vec( v.M_vec ), M_vecNonContiguousGhosts( v.M_vecNonContiguousGhosts ) { }
        
        void init( const size_type n, const size_type n_local, const bool fast = false ) override;

        virtual clone_ptrtype clone() const override;

        // Storage API
        virtual void resize( size_type n ) override;
        virtual void clear() override;

        virtual vector_cstptr_variant_type vec() const override { return &M_vec; }
        virtual vector_cstptr_variant_type vecNonContiguousGhosts() const override { return &M_vecNonContiguousGhosts; }

        // Operators API
        virtual value_type operator()( size_type i ) const override;
        virtual value_type& operator()( size_type i ) override;
        
        // Iterators API
        virtual iterator begin() override { return iterator( M_vec.begin() ); }
        virtual iterator end() override { return iterator( M_vec.end() ); }

        virtual iterator beginActive() override { return iterator( M_vec.begin() ); }
        virtual iterator endActive() override { return iterator( M_vec.end() ); }

        virtual iterator beginGhost() override { return iterator( M_vecNonContiguousGhosts.begin() ); }
        virtual iterator endGhost() override { return iterator( M_vecNonContiguousGhosts.end() ); }

        // Setters API
        virtual void setConstant( value_type v ) override;
        virtual void setZero() override;

        virtual void scale( const value_type factor ) override;

        // Utilities
        virtual real_type min( bool parallel ) const override;
        virtual real_type max( bool parallel ) const override;

        virtual real_type l1Norm() const override;
        virtual real_type l2Norm() const override;
        virtual real_type linftyNorm() const override;

        virtual value_type sum() const override;

    protected:
        virtual vector_ptr_variant_type vec() override { return &M_vec; }
        virtual vector_ptr_variant_type vecNonContiguousGhosts() override { return &M_vecNonContiguousGhosts; }
        
        void checkInvariants() const override;
        
        void setVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void setVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        void addVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void addVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        void maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;
        
        void subVector( const VectorUblasContiguousGhostsBase<T> & v ) override;
        void subVector( const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        void msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v ) override;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v ) override;

        value_type dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const override;
        value_type dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const override;
        
        VectorUblasBase<T> * rangeImpl( const range_type & rangeActive, const range_type & rangeGhost ) override;
        VectorUblasBase<T> * sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost ) override;

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
        VectorUblasRange( VectorUblasContiguousGhosts<T, Storage> & v, const range_type & rangeActive, const range_type & rangeGhost );
        VectorUblasRange( VectorUblasNonContiguousGhosts<T, Storage> & v, const range_type & rangeActive, const range_type & rangeGhost );
        VectorUblasRange( ublas::vector<T> & v, const range_type & range, const datamap_ptrtype & dm );

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
        VectorUblasSlice( VectorUblasContiguousGhosts<T, Storage> & v, const slice_type & sliceActive, const slice_type & sliceGhost );
        VectorUblasSlice( VectorUblasNonContiguousGhosts<T, Storage> & v, const slice_type & sliceActive, const slice_type & sliceGhost );
        VectorUblasSlice( Storage & v, const slice_type & slice, const datamap_ptrtype & dm );
        VectorUblasSlice( Storage & vActive, slice_type const& sliceActive,
                Storage & vGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm );

};

} // detail
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
    using namespace detail;
    if ( v.comm().size() > 1 )
    {
        const VectorUblasContiguousGhosts<T> * vecContiguousGhosts = dynamic_cast<const VectorUblasContiguousGhosts<T> *>( & v.vectorImpl() );
        if ( vecContiguousGhosts )
            return std::make_shared<VectorPetscMPI<T>>( std::addressof( * vecContiguousGhosts->begin() ), vecContiguousGhosts->mapPtr() );

        const VectorUblasNonContiguousGhosts<T> * vecNonContiguousGhosts = dynamic_cast<const VectorUblasNonContiguousGhosts<T> *>( & v.vectorImpl() );
        if ( vecNonContiguousGhosts )
            return std::make_shared<VectorPetscMPIRange<T>>( std::addressof( * vecNonContiguousGhosts->beginActive() ), std::addressof( * vecNonContiguousGhosts->beginGhost() ), vecNonContiguousGhosts->mapPtr() );
        
        else
            CHECK( false ) << "Unsupported VectorUblas for VectorPetsc view";
    }
    else
        return std::make_shared<VectorPetsc<T>>( std::addressof( * v.vectorImpl().begin() ), v.mapPtr() );
}

template< typename T >
inline VectorPetsc<T>
toPETSc( const VectorUblas<T> & v )
{
    return *toPetscPtr( v );
}

} // Feel

#endif // FEELPP_HAS_PETSC


#endif /* _FEELPP_VECTORUBLAS_HPP */
