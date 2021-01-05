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
//! @date 11 Apr 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_MATRIXCONDENSED_HPP
#define FEELPP_MATRIXCONDENSED_HPP 1

#include <feel/feelalg/matrixblock.hpp>


#include <feel/feelalg/staticcondensation.hpp>



namespace Feel {

//!
//! Matrix to represent statically condensed matrices 
//!
template<typename T>
class MatrixCondensed  : public MatrixBlockBase<T>
{
public:
    
    using super = MatrixBlockBase<T>;
    using size_type = typename super::size_type;
    using graph_ptrtype = typename super::graph_ptrtype;
    using value_type = T;
    using sc_type = StaticCondensation<value_type>;
    using sc_ptrtype = std::shared_ptr<sc_type>;
    using this_matrix_ptrtype = std::shared_ptr<MatrixCondensed<value_type>>;
        
    MatrixCondensed()
        :
        super(Environment::worldCommPtr()),
        M_sc( new sc_type ),
        M_strategy( solve::strategy::static_condensation )
        {}
    explicit MatrixCondensed( worldcomm_ptr_t const& wc )
        :
        super( wc ),
        M_sc( new sc_type ),
        M_strategy( solve::strategy::static_condensation )
        {}
    template<typename BackendT>
    MatrixCondensed( vf::BlocksBase<graph_ptrtype> const& graph,
                     BackendT&& b,
                     bool diag_is_nonzero = true )
        :
        super( graph, std::forward<BackendT>(b), true ),
        M_sc( new sc_type ),
        M_strategy( solve::strategy::static_condensation )
        {}
    template<typename BackendT>
    MatrixCondensed( solve::strategy s,
                     vf::BlocksBase<graph_ptrtype> const& graph,
                     BackendT&& b,
                     bool diag_is_nonzero = true )
        :
        super( graph, std::forward<BackendT>(b), true ),
        M_sc( new sc_type ),
        M_strategy( s )
        {}
#if 0
    template<typename... Args>
    MatrixCondensed( Args&&... params )
        :
         super( std::forward<Args>( params )... ),
         M_sc( new sc_type ),
         M_strategy( solve::strategy::static_condensation )
        {}
#endif
    MatrixCondensed( MatrixCondensed const& ) = default;
    MatrixCondensed& operator=( MatrixCondensed const& ) = default;
    ~MatrixCondensed() override = default;

    //!
    //! @return true if strategy is monolithic, false otherwise
    //!
    bool monolithic() const { return M_strategy == solve::strategy::monolithic; }
    //!
    //! @return true if strategy is static condensation, false otherwise
    //!
    bool staticCondensation() const { return M_strategy == solve::strategy::static_condensation; }

    //!
    //! @return true if strategy is static condensation, false otherwise
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


    void addMatrix ( int* rows, int nrows,
                             int* cols, int ncols,
                             value_type* data,
                             size_type K, size_type K2
                             ) override
        {
            tic();
            if ( staticCondensation() || localSolve() )
                M_sc->addLocalMatrix( rows, nrows, cols, ncols, data, K, K2 );
            else
                super::addMatrix( rows, nrows, cols, ncols, data );
            toc("addMatrix",FLAGS_v>2);
        }
    void addMatrix( const value_type a, MatrixSparse<value_type> const& M, Feel::MatrixStructure matStruc = Feel::SAME_NONZERO_PATTERN ) override
        {
            super::addMatrix( a, M, matStruc );
        }
    sc_ptrtype sc() { return M_sc; }
    sc_ptrtype const& sc() const { return M_sc; }
    sc_ptrtype const& sc( int row, int col ) const { M_sc->block( row, col );return M_sc; }
    sparse_matrix_ptrtype block( int row, int col )
        {
            if ( staticCondensation() || localSolve() )
                M_sc->block( row, col );
            return this->shared_from_this();
        }
    //!
    //! @return the number of non-zero entries
    //!
    size_type nnz() const override
        {
            if ( staticCondensation() )
                return M_sc->nnz();
            else
                return super::nnz();
        }
    void zero() override
        {
            if ( staticCondensation() )
                M_sc->zeroMatrix();
            else
                super::zero();
        }
    //!
    //! zero block with coordinates @p n1,n2
    //!
    void zeroBlock( int n1, int n2 )
        {
            if ( staticCondensation() )
                M_sc->zero( n1, n2 );
        }
    /**
     * transpose block \p n1,n2  and store it in \p n2,n1
     *
     * \param n1 row of the block to be transposed
     * \param n2 column of the block to be transposed
     */
    void transposeBlock( int n1, int n2 )
        {
            if ( staticCondensation() )
                M_sc->transpose( n1, n2 );
        }
private:
    sc_ptrtype M_sc;
    solve::strategy M_strategy;
};

template<typename T>
using condensed_matrix_t = MatrixCondensed<T>;
template<typename T>
using condensed_matrix_ptr_t = std::shared_ptr<condensed_matrix_t<T>>;

//!
//! Create a shared pointer \p MatrixCondensed<T> from \p Args
//! @code
//! // create a MatrixCondensed std::shared_ptr
//! auto mc = makeSharedMatrixCondensed<double>();
//! @endcode
//!
template< class T, class... Args >
condensed_matrix_ptr_t<T>
makeSharedMatrixCondensed( Args&&... args )
{
    return std::make_shared<MatrixCondensed<T>>( args... );
}

#if !defined(FEELPP_MATRIXCONDENSED_NOEXTERN)
extern template class MatrixCondensed<double>;
//extern template class Backend<std::complex<double>>;
#endif

}


#endif
