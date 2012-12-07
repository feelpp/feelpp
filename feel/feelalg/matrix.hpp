/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-05-22

  Copyright (C) 2006 EPFL

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
   \file matrix.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-05-22
 */
#if !defined(FEELPP_MATRIX_HPP)
#define FEELPP_MATRIX_HPP 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>

namespace boost
{
namespace numeric
{
namespace ublas
{

/// \cond detail
// Identity matrix class
template<class T>
class anti_identity_matrix:
    public matrix_container<anti_identity_matrix<T> >
{

    typedef const T *const_pointer;
    typedef anti_identity_matrix<T> self_type;
public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
    using matrix_container<self_type>::operator ();
#endif
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef T value_type;
    typedef const T &const_reference;
    typedef T &reference;
    typedef const matrix_reference<const self_type> const_closure_type;
    typedef matrix_reference<self_type> closure_type;
    typedef sparse_tag storage_category;
    typedef unknown_orientation_tag orientation_category;

    // Construction and destruction
    BOOST_UBLAS_INLINE
    anti_identity_matrix ():
        matrix_container<self_type> (),
        size1_ ( 0 ), size2_ ( 0 ), size_common_ ( 0 ) {}
    BOOST_UBLAS_INLINE
    anti_identity_matrix ( size_type size ):
        matrix_container<self_type> (),
        size1_ ( size ), size2_ ( size ), size_common_ ( ( std::min ) ( size1_, size2_ ) ) {}
    BOOST_UBLAS_INLINE
    anti_identity_matrix ( size_type size1, size_type size2 ):
        matrix_container<self_type> (),
        size1_ ( size1 ), size2_ ( size2 ), size_common_ ( ( std::min ) ( size1_, size2_ ) ) {}
    BOOST_UBLAS_INLINE
    anti_identity_matrix ( const anti_identity_matrix &m ):
        matrix_container<self_type> (),
        size1_ ( m.size1_ ), size2_ ( m.size2_ ), size_common_ ( ( std::min ) ( size1_, size2_ ) ) {}

    // Accessors
    BOOST_UBLAS_INLINE
    size_type size1 () const
    {
        return size1_;
    }
    BOOST_UBLAS_INLINE
    size_type size2 () const
    {
        return size2_;
    }

    // Resizing
    BOOST_UBLAS_INLINE
    void resize ( size_type size, bool preserve = true )
    {
        size1_ = size;
        size2_ = size;
    }
    BOOST_UBLAS_INLINE
    void resize ( size_type size1, size_type size2, bool /*preserve*/ = true )
    {
        size1_ = size1;
        size2_ = size2;
    }

    // Element access
    BOOST_UBLAS_INLINE
    const_reference operator () ( size_type i, size_type j ) const
    {
        if ( i == size2_-( j+1 ) )
            return one_;

        else
            return zero_;
    }

    // Assignment
    BOOST_UBLAS_INLINE
    anti_identity_matrix &operator = ( const anti_identity_matrix &m )
    {
        size1_ = m.size1_;
        size2_ = m.size2_;
        return *this;
    }
    BOOST_UBLAS_INLINE
    anti_identity_matrix &assign_temporary ( anti_identity_matrix &m )
    {
        swap ( m );
        return *this;
    }

    // Swapping
    BOOST_UBLAS_INLINE
    void swap ( anti_identity_matrix &m )
    {
        if ( this != &m )
        {
            std::swap ( size1_, m.size1_ );
            std::swap ( size2_, m.size2_ );
        }
    }
    BOOST_UBLAS_INLINE
    friend void swap ( anti_identity_matrix &m1, anti_identity_matrix &m2 )
    {
        m1.swap ( m2 );
    }

    // Iterator types
private:
    // Use an index
    typedef size_type const_subiterator_type;

public:
    class const_iterator1;
    class const_iterator2;
    typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
    typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;

    // Element lookup
    BOOST_UBLAS_INLINE
    const_iterator1 find1 ( int rank, size_type i, size_type j ) const
    {
        if ( rank == 1 )
        {
            i = ( std::max ) ( i, j );
            i = ( std::min ) ( i, j + 1 );
        }

        return const_iterator1 ( *this, i );
    }
    BOOST_UBLAS_INLINE
    const_iterator2 find2 ( int rank, size_type i, size_type j ) const
    {
        if ( rank == 1 )
        {
            j = ( std::max ) ( j, i );
            j = ( std::min ) ( j, i + 1 );
        }

        return const_iterator2 ( *this, j );
    }


    class const_iterator1:
        public container_const_reference<anti_identity_matrix>,
        public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
        const_iterator1, value_type>
    {
    public:
        typedef typename anti_identity_matrix::value_type value_type;
        typedef typename anti_identity_matrix::difference_type difference_type;
        typedef typename anti_identity_matrix::const_reference reference;
        typedef typename anti_identity_matrix::const_pointer pointer;

        typedef const_iterator2 dual_iterator_type;
        typedef const_reverse_iterator2 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator1 ():
            container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator1 ( const self_type &m, const const_subiterator_type &it ):
            container_const_reference<self_type> ( m ), it_ ( it ) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator1 &operator ++ ()
        {
            BOOST_UBLAS_CHECK ( it_ < ( *this ) ().size1 (), bad_index () );
            ++it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator1 &operator -- ()
        {
            BOOST_UBLAS_CHECK ( it_ > 0, bad_index () );
            --it_;
            return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const
        {
            return one_;
        }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 begin () const
        {
            return const_iterator2 ( ( *this ) (), it_ );
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator2 end () const
        {
            return const_iterator2 ( ( *this ) (), it_ + 1 );
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 rbegin () const
        {
            return const_reverse_iterator2 ( end () );
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator2 rend () const
        {
            return const_reverse_iterator2 ( begin () );
        }
#endif

        // Indices
        BOOST_UBLAS_INLINE
        size_type index1 () const
        {
            return it_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const
        {
            return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator1 &operator = ( const const_iterator1 &it )
        {
            container_const_reference<self_type>::assign ( &it () );
            it_ = it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == ( const const_iterator1 &it ) const
        {
            BOOST_UBLAS_CHECK ( &( *this ) () == &it (), external_logic () );
            return it_ == it.it_;
        }

    private:
        const_subiterator_type it_;
    };

    typedef const_iterator1 iterator1;

    BOOST_UBLAS_INLINE
    const_iterator1 begin1 () const
    {
        return const_iterator1 ( *this, 0 );
    }
    BOOST_UBLAS_INLINE
    const_iterator1 end1 () const
    {
        return const_iterator1 ( *this, size_common_ );
    }

    class const_iterator2:
        public container_const_reference<anti_identity_matrix>,
        public bidirectional_iterator_base<sparse_bidirectional_iterator_tag,
        const_iterator2, value_type>
    {
    public:
        typedef typename anti_identity_matrix::value_type value_type;
        typedef typename anti_identity_matrix::difference_type difference_type;
        typedef typename anti_identity_matrix::const_reference reference;
        typedef typename anti_identity_matrix::const_pointer pointer;

        typedef const_iterator1 dual_iterator_type;
        typedef const_reverse_iterator1 dual_reverse_iterator_type;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        const_iterator2 ():
            container_const_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        const_iterator2 ( const self_type &m, const const_subiterator_type &it ):
            container_const_reference<self_type> ( m ), it_ ( it ) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator2 &operator ++ ()
        {
            BOOST_UBLAS_CHECK ( it_ < ( *this ) ().size_common_, bad_index () );
            ++it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator2 &operator -- ()
        {
            BOOST_UBLAS_CHECK ( it_ > 0, bad_index () );
            --it_;
            return *this;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference operator * () const
        {
            return one_;
        }

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 begin () const
        {
            return const_iterator1 ( ( *this ) (), it_ );
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_iterator1 end () const
        {
            return const_iterator1 ( ( *this ) (), it_ + 1 );
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 rbegin () const
        {
            return const_reverse_iterator1 ( end () );
        }
        BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
        typename self_type::
#endif
        const_reverse_iterator1 rend () const
        {
            return const_reverse_iterator1 ( begin () );
        }
#endif

        // Indices
        BOOST_UBLAS_INLINE
        size_type index1 () const
        {
            return it_;
        }
        BOOST_UBLAS_INLINE
        size_type index2 () const
        {
            return it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator2 &operator = ( const const_iterator2 &it )
        {
            container_const_reference<self_type>::assign ( &it () );
            it_ = it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == ( const const_iterator2 &it ) const
        {
            BOOST_UBLAS_CHECK ( &( *this ) () == &it (), external_logic () );
            return it_ == it.it_;
        }

    private:
        const_subiterator_type it_;
    };

    typedef const_iterator2 iterator2;

    BOOST_UBLAS_INLINE
    const_iterator2 begin2 () const
    {
        return const_iterator2 ( *this, 0 );
    }
    BOOST_UBLAS_INLINE
    const_iterator2 end2 () const
    {
        return const_iterator2 ( *this, size_common_ );
    }

    // Reverse iterators

    BOOST_UBLAS_INLINE
    const_reverse_iterator1 rbegin1 () const
    {
        return const_reverse_iterator1 ( end1 () );
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator1 rend1 () const
    {
        return const_reverse_iterator1 ( begin1 () );
    }

    BOOST_UBLAS_INLINE
    const_reverse_iterator2 rbegin2 () const
    {
        return const_reverse_iterator2 ( end2 () );
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator2 rend2 () const
    {
        return const_reverse_iterator2 ( begin2 () );
    }

private:
    size_type size1_;
    size_type size2_;
    size_type size_common_;
    static const value_type zero_;
    static const value_type one_;
};

template<class T>
const typename anti_identity_matrix<T>::value_type anti_identity_matrix<T>::zero_ ( 0 );
template<class T>
const typename anti_identity_matrix<T>::value_type anti_identity_matrix<T>::one_ ( 1 );
/// \endcond detail

}
}
}

#endif /* FEELPP_MATRIX_HPP */

