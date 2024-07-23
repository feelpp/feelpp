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
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/staticcondensation.hpp>


namespace Feel {

//!
//! Vector to represent statically condensed matrices 
//!
template<typename T>
class VectorCondensed  : public VectorBlockBase<T>
{
public:
    
    using super = VectorBlockBase<T>;
    using value_type = T;
    using size_type = typename super::size_type;
    using vector_ptrtype = typename super::vector_ptrtype;
    using backend_type = typename super::backend_type;
    using backend_ptrtype = typename super::backend_ptrtype;
    using sc_type = StaticCondensation<value_type>;
    using sc_ptrtype = std::shared_ptr<sc_type>;
    using this_vector_ptrtype = std::shared_ptr<VectorCondensed<value_type>>;
    
    VectorCondensed()
        :
        super(),
        M_sc( new sc_type ),
        M_strategy( solve::strategy::static_condensation )
        {}

    template<typename BackendT>
    VectorCondensed( vf::BlocksBase<vector_ptrtype> const & b, BackendT&& ba, bool copy_values=true )
        :
        super( b, std::forward<BackendT>(ba), copy_values ),
        M_sc( new sc_type ),
        M_strategy( solve::strategy::static_condensation )
        {}
    //!
    //! construct a VectorCondensed according to strategy \p s
    //! @param b blocks base data structure holding block vectors
    //! @param ba backend to use
    //!
    template<typename BackendT>
    VectorCondensed( solve::strategy s, vf::BlocksBase<vector_ptrtype> const & b, BackendT&& ba, bool copy_values=true )
        :
        super( b, std::forward<BackendT>(ba), copy_values ),
        M_sc( new sc_type ),
        M_strategy( s )
        {}
#if 0    
    template<typename... Args>
    VectorCondensed( Args&&... params )
        :
         super( std::forward<Args>( params )... ),
         M_sc( new sc_type ),
         M_strategy( solve::strategy::static_condensation )
        {}
#endif
    VectorCondensed( VectorCondensed const& ) = default;
    VectorCondensed& operator=( VectorCondensed const& ) = default;
    ~VectorCondensed() override = default;

    //!
    //! @return true if strategy is monolithic, false otherwise
    //!
    bool monolithic() const { return M_strategy == solve::strategy::monolithic; }
    //!
    //! @return true if strategy is static condensation, false otherwise
    //!
    bool staticCondensation() const { return M_strategy == solve::strategy::static_condensation; }

    //!
    //! @return true if strategy is local, false otherwise
    //!
    bool localSolve() const { return M_strategy == solve::strategy::local; }

    //!
    //! get the strategy 
    //!
    solve::strategy solveStrategy() const { return M_strategy; }
    
    //!
    //! set the strategy \p s
    //!
    void setStrategy( solve::strategy s ) { M_strategy = s; }

    //!
    //! add local vector according to strategy
    //!
    void addVector ( int* rows, int nrows,
                             value_type* data,
                             size_type K, size_type K2
                             ) override
        {
            tic();
            if ( staticCondensation() || localSolve() )
            {
                auto add_v = [&]()
                    {
                        M_sc->addLocalVector( rows, nrows, data, K, K2 );
                    };
                //M_futs.push_back(std::async( std::launch::async, add_v ));
                add_v();
            }
            else
                super::addVector( rows, nrows, data, K, K2 );
            toc("Vector::addVector",FLAGS_v>2);
        }

    sc_ptrtype sc() { getFuture();return M_sc; }
    sc_ptrtype const& sc() const { getFuture();return M_sc; }
    sc_ptrtype const& sc( int row ) const { M_sc->block( row );return M_sc; }
    vector_ptrtype block( int row )
        {
            if ( staticCondensation() || localSolve() )
                M_sc->block( row );
            return this->shared_from_this();
        }
    void zero() override
        {
            if ( staticCondensation() )
                M_sc->zeroVector();
            else
                super::zero();
        }
    //!
    //! zero block @param n1
    //!
    void zeroBlock( int n1 )
        {
            if ( staticCondensation() )
                M_sc->zero( n1 );
        }
protected:
    void getFuture() const
        {
            for( auto& f : M_futs )
            {
                f.get();
            }
        }
private:
    sc_ptrtype M_sc;
    solve::strategy M_strategy;
    mutable std::vector<std::future<void>> M_futs;
};

template<typename T>
using condensed_vector_t = VectorCondensed<T>;
template<typename T>
using condensed_vector_ptr_t = std::shared_ptr<condensed_vector_t<T>>;

template<typename V>
struct is_vectorblock: std::is_base_of<VectorBlockBase<typename V::value_type>,V> {};

template<typename V>
constexpr bool is_vectorblock_v = is_vectorblock<V>::value;

template<typename V>
struct is_vectorcondensed: std::is_base_of<VectorCondensed<typename V::value_type>,V> {};

template<typename V>
constexpr bool is_vectorcondensed_v = is_vectorcondensed<V>::value;




//!
//! Create a shared pointer \p VectorCondensed<T> from \p Args
//! @code
//! // create a VectorCondensed std::shared_ptr
//! auto mc = makeSharedVectorCondensed<double>();
//! @endcode
//!
template< class T, class... Args >
condensed_vector_ptr_t<T>
makeSharedVectorCondensed( Args&&... args )
{
    return std::make_shared<VectorCondensed<T>>( args... );
}

#if !defined(FEELPP_VECTORCONDENSED_NOEXTERN)
extern template class VectorCondensed<double>;
#endif



}
#endif
