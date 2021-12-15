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

template< typename T >
class VectorUblasBase;
template< typename T >
class VectorUblasContiguousGhosts;
template< typename T >
class VectorUblasNonContiguousGhosts;

} // namespace detail

template< typename T >
class VectorUblas : public Vector<T>
{
    public:
        // Typedefs
        typedef VectorUblas<T> self_type;
        typedef Vector<T> super_type;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;

        //typedef ublas::vector<value_type> storage_type;

        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

    public:
        // Constructors/Destructor
        VectorUblas();

        ~VectorUblas() override;
        //TODO
        super_type::clone_ptrtype clone() const override = 0;

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override { return M_vectorImpl->init( n, n_local, fast ); }
        void init( const size_type n, const bool fast = false ) override { return M_vectorImpl->init( n, fast ); }
        void init( const datamap_ptrtype & dm ) override { return M_vectorImpl->init( dm ); }
        
        void resize( size_type n ) { return M_vectorImpl->resize( n ); }
        void clear() override { return M_vectorImpl->clear(); }

        // Status API
        bool isInitialized() const override { return true; }
        void close() const override { }
        bool closed() const override { return true; }
        bool areGlobalValuesUpdated() const { return true; }
        void updateGlobalValues() const { }
        void outdateGlobalValues() { }

        // Operators API
        Vector<value_type>& operator=( const Vector<value_type> & V ) override;
        Vector<value_type>& operator=( const this_type & V );

        virtual value_type operator()( size_type i ) const override = 0;
        virtual value_type& operator()( size_type i ) override = 0;
        value_type operator[]( size_type i ) const { return this->operator()( i ); }
        value_type& operator[]( size_type i ) { return this->operator()( i ); }

        Vector<T>& operator+=( const Vector<T>& v ) override { this->add( v ); return *this; }
        VectorUblasExpression<T> operator+( const Vector<T>& v ) const;

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
        void setConstant( value_type a ) { return M_vectorImpl->setConstant( a ); }
        void setZero() { return M_vectorImpl->setZero(); }
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
        void printMatlab( const std::string & filename = "NULL", bool renumber = false ) const override { return M_vectorImpl->printMatlab( filename, renumber ); }
#ifdef FEELPP_HAS_HDF5
        void saveHDF5( const std::string & filename, const std::string & tableName = "element", bool appendMode = false ) const { return M_vectorImpl->saveHDF5( filename, tableName, appendMode ); }
        void loadHDF5( const std::string & filename, const std::string & tableName = "element" ) { return M_vectorImpl->loadHDF5( filename, tableName ); }
#endif

        // Localization (parallel global to one proc local)
        void localizeToOneProcessor( ublas::vector<T> & v_local, const size_type proc_id = 0 ) const { return M_vectorImpl->localizeToOneProcessor( v_local, proc_id ); }
        void localizeToOneProcessor( std::vector<T> & v_local, const size_type proc_id = 0 ) const { return M_vectorImpl->localizeToOneProcessor( v_local, proc_id ); }
        
    private:
        std::unique_ptr<VectorUblasBase> M_vectorImpl;
};

namespace detail 
{

template< typename T >
class VectorUblasBase: public Vector<T>
{
    public:
        // Typedefs
        typedef Vector<T> super_type;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        // Actual storage variants
        typedef ublas::vector<value_type> vector_storage_type;
        typedef ublas::vector_range<ublas::vector<value_type>> vector_range_storage_type;
        typedef ublas::vector_slice<ublas::vector<value_type>> vector_slice_storage_type;
        typedef ublas::vector<value_type, Feel::detail::shallow_array_adaptor<value_type>> vector_map_storage_type;
        typedef std::variant<
            vector_storage_type *,
            vector_range_storage_type *,
            vector_slice_storage_type *,
            vector_map_storage_type * > vector_ptr_variant_type;
        typedef std::variant<
            const vector_storage_type *,
            const vector_range_storage_type *,
            const vector_slice_storage_type *,
            const vector_map_storage_type * > vector_cstptr_variant_type;

        // Iterator class
        class iterator
        {
            public:
                template< typename It >
                iterator( const It & it ): M_iteratorImpl( new iterator_impl<It>( it ) ) { }
                iterator() = delete;
                
                iterator( const iterator & other ): M_iteratorImpl( other.M_iteratorImpl ? other.M_iteratorImpl->clone() : nullptr ) { }
                iterator & swap( iterator & other ) { std::swap( M_iteratorImpl, other.M_iteratorImpl ); return *this; }
                iterator & operator=( const iterator & other ) { swap( iterator( other ) ); return *this; }
                ~iterator() { delete M_iteratorImpl; }

                value_type & operator*() const { return M_iteratorImpl->current(); }
                iterator & operator++() { M_iteratorImpl->next(); return *this; }
                iterator & operator--() { M_iteratorImpl->previous(); return *this; }
                bool operator==( const iterator & other ) const { return M_iteratorImpl->equal( *other.M_iteratorImpl ); }
                iterator & operator=( const iterator & other ) { assign( other ); return *this; }

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

        ~VectorUblasBase() override;

        virtual super_type::clone_ptrtype clone() const override = 0;

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override;
        void init( const size_type n, const bool fast = false ) override;
        void init( const datamap_ptrtype & dm ) override;
        
        virtual void resize( size_type n ) = 0;
        virtual void clear() override = 0;

        // Status API
        bool isInitialized() const override { return true; }
        void close() const { }
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
        VectorUblasExpression<T> operator+( const Vector<T>& v ) const;
        Vector<T>& operator-=( const Vector<T>& v ) override { this->sub( v ); return *this; }
        Vector<T>& operator*=( const value_type & a ) { this->scale( a ); return *this; }

        //// Iterators API
        //virtual iterator begin() = 0;
        //virtual iterator end() = 0;

        //virtual iterator beginGhost() = 0;
        //virtual iterator endGhost() = 0;

        //virtual size_type start() const = 0;
        //virtual size_type startNonContiguousGhosts() const = 0;

        //size_type rowStart() const { checkInvariants(); return 0; }
        //size_type rowStop() const { checkInvariants(); return 0; }

        // Setters API
        virtual void setConstant( value_type v ) = 0;
        virtual void setZero() = 0;
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
        void printMatlab( const std::string & filename = "NULL", bool renumber = false ) const override;
#ifdef FEELPP_HAS_HDF5
        void saveHDF5( const std::string & filename, const std::string & tableName = "element", bool appendMode = false ) const;
        void loadHDF5( const std::string & filename, const std::string & tableName = "element" );
#endif

        // Localization (parallel global to one proc local)
        void localizeToOneProcessor( ublas::vector<T> & v_local, const size_type proc_id = 0 ) const;
        void localizeToOneProcessor( std::vector<T> & v_local, const size_type proc_id = 0 ) const;

    protected:
        void setVector( const VectorUblasBase<T> & v ) { return v.applySetVector( *this ); }
        virtual void applySetVector( const VectorUblasBase<T> & v ) = 0;
        virtual void setVector( const VectorUblasContiguousGhosts<T> & v ) = 0;
        virtual void setVector( const VectorUblasNonContiguousGhosts<T> & v ) = 0;
        
        void addVector( const VectorUblasBase<T> & v ) { return v.applyAddVector( *this ); }
        virtual void applyAddVector( const VectorUblasBase<T> & v ) = 0;
        virtual void addVector( const VectorUblasContiguousGhosts<T> & v ) = 0;
        virtual void addVector( const VectorUblasNonContiguousGhosts<T> & v ) = 0;

        void maddVector( const value_type & a, const VectorUblasBase<T> & v ) { return v.applyMaddVector( a, *this ); }
        virtual void applyMaddVector( const value_type & a, const VectorUblasBase<T> & v ) = 0;
        virtual void maddVector( const value_type & a, const VectorUblasContiguousGhosts<T> & v ) = 0;
        virtual void maddVector( const value_type & a, const VectorUblasNonContiguousGhosts<T> & v ) = 0;
        
        void subVector( const VectorUblasBase<T> & v ) { return v.applySubVector( *this ); }
        virtual void applySubVector( const VectorUblasBase<T> & v ) = 0;
        virtual void subVector( const VectorUblasContiguousGhosts<T> & v ) = 0;
        virtual void subVector( const VectorUblasNonContiguousGhosts<T> & v ) = 0;

        virtual void msubVector( const value_type & a, const VectorUblasBase<T> & v ) { return v.applyMsubVector( a, *this ); }
        virtual void applyMsubVector( const value_type & a, const VectorUblasBase<T> & v ) = 0;
        virtual void msubVector( const value_type & a, const VectorUblasContiguousGhosts<T> & v ) = 0;
        virtual void msubVector( const value_type & a, const VectorUblasNonContiguousGhosts<T> & v ) = 0;

        virtual value_type dotVector( const VectorUblasBase<T> & v ) { return v.applyDotVector( *this ); }
        virtual value_type applyDotVector( const VectorUblasBase<T> & v ) = 0;
        virtual value_type dotVector( const VectorUblasContiguousGhosts<T> & v ) = 0;
        virtual value_type dotVector( const VectorUblasNonContiguousGhosts<T> & v ) = 0;

};

template< template < typename > class V, typename T >
class SettableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applySetVector( VectorUblasBase<T> & b ) override { return b.setVector( static_cast< V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class AddableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applyAddVector( VectorUblasBase<T> & b ) override { return b.addVector( static_cast< V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class MaddableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applyMaddVector( VectorUblasBase<T> & b ) override { return b.maddVector( static_cast< V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class SubtractableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applySubVector( VectorUblasBase<T> & b ) override { return b.subVector( static_cast< V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class MsubtractableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        void applyMsubVector( VectorUblasBase<T> & b ) override { return b.msubVector( static_cast< V<T>& >( *this ) ); }
};
template< template < typename > class V, typename T >
class DottableVectorUblas: public virtual VectorUblasBase<T>
{
    protected:
        typename VectorUblasBase<T>::value_type applyDotVector( VectorUblasBase<T> & b ) override { return b.dotVector( static_cast< V<T>& >( *this ) ); }
};

template< typename T >
class VectorUblasContiguousGhosts: 
    public SettableVectorUblas<VectorUblasContiguousGhosts, T>,
    public AddableVectorUblas<VectorUblasContiguousGhosts, T>,
    public MaddableVectorUblas<VectorUblasContiguousGhosts, T>,
    public SubtractableVectorUblas<VectorUblasContiguousGhosts, T>,
    public MsubtractableVectorUblas<VectorUblasContiguousGhosts, T>,
    public DottableVectorUblas<VectorUblasContiguousGhosts, T>
{
    public:
        // Typedefs
        typedef VectorUblasBase<T> super_type;

        typedef T value_type;
        typedef typename type_traits<value_type>::real_type real_type;
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;

        using super_type::vector_ptr_variant_type;
        using super_type::vector_cstptr_variant_type;

    public:
        // Constructors/Destructor
        VectorUblasContiguousGhosts( ) = default;
        VectorUblasContiguousGhosts( VectorUblasContiguousGhosts<T> const& v ): super_type(v), M_vec( v.M_vec ) { }

        VectorUblasContiguousGhosts( size_type s );
        VectorUblasContiguousGhosts( const datamap_ptrtype & dm );
        VectorUblasContiguousGhosts( size_type s, size_type n_local );

        ~VectorUblasContiguousGhosts() override;
        
        virtual super_type::clone_ptrtype clone() const override;

        // Storage API
        virtual void resize( size_type n ) override;
        virtual void clear() override;

        vector_cstptr_variant_type vec() const { return &M_vec; }

        // Operators API
        virtual value_type operator()( size_type i ) const override;
        virtual value_type& operator()( size_type i ) override;

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
        void setVector( const VectorUblasContiguousGhosts<T> & v ) override;
        void setVector( const VectorUblasNonContiguousGhosts<T> & v ) override;

        void addVector( const VectorUblasContiguousGhosts<T> & v ) override;
        void addVector( const VectorUblasNonContiguousGhosts<T> & v ) override;

        void maddVector( const value_type & a, const VectorUblasContiguousGhosts<T> & v ) override;
        void maddVector( const value_type & a, const VectorUblasNonContiguousGhosts<T> & v ) override;
        
        void subVector( const VectorUblasContiguousGhosts<T> & v ) override;
        void subVector( const VectorUblasNonContiguousGhosts<T> & v ) override;

        void msubVector( const value_type & a, const VectorUblasContiguousGhosts<T> & v ) override;
        void msubVector( const value_type & a, const VectorUblasNonContiguousGhosts<T> & v ) override;

        value_type dotVector( const VectorUblasContiguousGhosts<T> & v ) override;
        value_type dotVector( const VectorUblasNonContiguousGhosts<T> & v ) override;

    protected:
        ublas::vector<value_type> M_vec;

};

}

} // Feel
#endif


#endif /* _FEELPP_VECTORUBLAS_HPP */
