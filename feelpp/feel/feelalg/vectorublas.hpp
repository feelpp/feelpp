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

    protected:
        // Implementation classes
        class VectorUblasImplBase

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
        typedef typename super_type::datamap_type datamap_type;
        typedef typename super_type::datamap_ptrtype datamap_ptrtype;
        using size_type = typename datamap_type::size_type;


    public:
        // Constructors/Destructor
        VectorUblasBase( size_type s );
        VectorUblasBase( const datamap_ptrtype & dm );
        VectorUblasBase( size_type s, size_type n_local );

        ~VectorUblasBase() override;

        // Storage API
        void init( const size_type n, const size_type n_local, const bool fast = false ) override;
        void init( const size_type n, const bool fast = false ) override;
        void init( const datamap_ptrtype & dm ) override;

        // Operators API
        Vector<value_type>& operator=( const Vector<value_type> & V ) override;
        Vector<value_type>& operator=( const this_type & V );

        T operator()( size_type i ) const override = 0;
        T& operator()( size_type i ) override = 0;
        T operator[]( size_type i ) const { return this->operator()( i ); }
        T& operator[]( size_type i ) { return this->operator()( i ); }

        Vector<T>& operator+=( const Vector<T>& v ) override { this->add( v ); return *this; }
        //TODO

        // Iterators API

};

}

} // Feel
#endif


#endif /* _FEELPP_VECTORUBLAS_HPP */
