//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 28 Jul 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_VECTORCONDENSED_HPP
#define FEELPP_VECTORCONDENSED_HPP 1

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelalg/staticcondensation.hpp>


namespace Feel {

//!
//! Vector to represent statically condensed matrices 
//!
template<typename T, bool Condense = false>
class VectorCondensed  : public VectorBlockBase<T>
{
public:
    
    using super = VectorBlockBase<T>;
    using value_type = T;
    using vector_ptrtype = typename super::vector_ptrtype;
    using backend_type = typename super::backend_type;
    using backend_ptrtype = typename super::backend_ptrtype;
    using sc_type = StaticCondensation<value_type>;
    using sc_ptrtype = boost::shared_ptr<sc_type>;
    using this_vector_ptrtype = boost::shared_ptr<VectorCondensed<value_type>>;
    static const bool do_condense=Condense;
    
    VectorCondensed()
        :
        super(Environment::worldComm()),
        M_sc( new sc_type )
        {}
    VectorCondensed( WorldComm const& wc )
        :
        super( wc ),
        M_sc( new sc_type )
        {}
    VectorCondensed( vf::BlocksBase<vector_ptrtype> const & b, backend_type& ba, bool copy_values=true )
        :
        super( b, ba, copy_values ),
        M_sc( new sc_type )
        {}
    template<typename... Args>
    VectorCondensed( Args&&... params )
        :
         super( std::forward<Args>( params )... ),
         M_sc( new sc_type )
        {}
    VectorCondensed( VectorCondensed const& ) = default;
    VectorCondensed& operator=( VectorCondensed const& ) = default;
    ~VectorCondensed() = default;


    virtual void addVector ( int* rows, int nrows,
                             value_type* data,
                             size_type K, size_type K2
                             ) override
        {
            return addVectorImpl( rows, nrows, data, K, K2, mpl::bool_<do_condense>() );
        }
    void addVectorImpl ( int* rows, int nrows,
                         value_type* data,
                         size_type K, size_type K2,
                         mpl::bool_<true>
                         )
        {
            tic();
            //if ( K != invalid_size_type_value && K2 != invalid_size_type_value )
            M_sc->addLocalVector( rows, nrows, data, K, K2 );
            //M_sc->addLocalVector( rows, nrows, cols, ncols, data, K, K2 );
            toc("sc.addlocalvector",FLAGS_v>0);
        }
    void addVectorImpl ( int* rows, int nrows,
                         value_type* data,
                         size_type K, size_type K2,
                         mpl::bool_<false>
                         )
        {
            tic();
            super::addVector( rows, nrows, data, K, K2 );
            toc("mono.addlocalvector",FLAGS_v>0);
        }
    sc_ptrtype sc() { return M_sc; }
    sc_ptrtype const& sc() const { return M_sc; }
    sc_ptrtype const& sc( int row ) const { M_sc->block( row );return M_sc; }
    vector_ptrtype block( int row )
        {
            return block( row, mpl::bool_<Condense>() );
        }
    vector_ptrtype block( int row, mpl::bool_<true> )
        {
            //sc->setVector( this->shared_from_this() );
            M_sc->block( row );
            return this->shared_from_this();
        }
    vector_ptrtype block( int row, mpl::bool_<false> )
        {
            return this->shared_from_this();
        }
    
private:
    sc_ptrtype M_sc;
};

template<typename T, bool Condense=true>
using condensed_vector_t = VectorCondensed<T,Condense>;
template<typename T,bool Condense=true>
using condensed_vector_ptr_t = boost::shared_ptr<condensed_vector_t<T,Condense>>;

//!
//! Create a shared pointer \p VectorCondensed<T> from \p Args
//! @code
//! // create a VectorCondensed boost::shared_ptr
//! auto mc = makeSharedVectorCondensed<double>();
//! @endcode
//!
template< class T, class... Args >
condensed_vector_ptr_t<T>
makeSharedVectorCondensed( Args&&... args )
{
    return boost::make_shared<VectorCondensed<T>>( args... );
}

#if !defined(FEELPP_VECTORCONDENSED_NOEXTERN)
extern template class VectorCondensed<double>;
//extern template class Backend<std::complex<double>>;
#endif



}
#endif
