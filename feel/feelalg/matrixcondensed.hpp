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
template<typename T, bool Condense = false>
class MatrixCondensed  : public MatrixBlockBase<T>
{
public:
    
    using super = MatrixBlockBase<T>;
    using graph_ptrtype = typename super::graph_ptrtype;
    using value_type = T;
    using sc_type = StaticCondensation<value_type>;
    using sc_ptrtype = boost::shared_ptr<sc_type>;
    using this_matrix_ptrtype = boost::shared_ptr<MatrixCondensed<value_type>>;
    static const bool do_condense=Condense;
    
    MatrixCondensed()
        :
        super(Environment::worldComm()),
        M_sc( new sc_type )
        {}
    MatrixCondensed( WorldComm const& wc )
        :
        super( wc ),
        M_sc( new sc_type )
        {}
    MatrixCondensed( vf::BlocksBase<graph_ptrtype> const& graph,
                     backend_ptrtype backend,
                     bool diag_is_nonzero = true )
        :
        super( graph, backend, true ),
        M_sc( new sc_type )
        {}
#if 0    
    MatrixCondensed ( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, WorldComm const& worldComm )
        :
        super( dmRow, dmCol, worldComm ),
        M_sc( new sc_type )
        {}
#endif
    template<typename... Args>
    MatrixCondensed( Args&&... params )
        :
         super( std::forward<Args>( params )... ),
        M_sc( new sc_type )
        {}
    MatrixCondensed( MatrixCondensed const& ) = default;
    MatrixCondensed& operator=( MatrixCondensed const& ) = default;
    ~MatrixCondensed() = default;


    virtual void addMatrix ( int* rows, int nrows,
                             int* cols, int ncols,
                             value_type* data,
                             size_type K, size_type K2
                             ) override
        {
            return addMatrixImpl( rows, nrows, cols, ncols, data, K, K2, mpl::bool_<do_condense>() );
        }
    void addMatrixImpl ( int* rows, int nrows,
                         int* cols, int ncols,
                         value_type* data,
                         size_type K, size_type K2,
                         mpl::bool_<true>
                         )
        {
            tic();
            //if ( K != invalid_size_type_value && K2 != invalid_size_type_value )
            M_sc->addLocalMatrix( rows, nrows, cols, ncols, data, K, K2 );
            //M_sc->addLocalMatrix( rows, nrows, cols, ncols, data, K, K2 );
            toc("sc.addlocalmatrix",FLAGS_v>0);
        }
    void addMatrixImpl ( int* rows, int nrows,
                         int* cols, int ncols,
                         value_type* data,
                         size_type K, size_type K2,
                         mpl::bool_<false>
                         )
        {
            tic();
            super::addMatrix( rows, nrows, cols, ncols, data );
            toc("mono.addlocalmatrix",FLAGS_v>0);
        }
    sc_ptrtype sc() { return M_sc; }
    sc_ptrtype const& sc() const { return M_sc; }
    sc_ptrtype const& sc( int row, int col ) const { M_sc->block( row, col );return M_sc; }
    sparse_matrix_ptrtype block( int row, int col )
        {
            return block( row, col, mpl::bool_<Condense>() );
        }
    sparse_matrix_ptrtype block( int row, int col, mpl::bool_<true> )
        {
            //sc->setMatrix( this->shared_from_this() );
            M_sc->block( row, col );
            return this->shared_from_this();
        }
    sparse_matrix_ptrtype block( int row, int col, mpl::bool_<false> )
        {
            return this->shared_from_this();
        }
private:
    sc_ptrtype M_sc;
};

template<typename T, bool Condense=true>
using condensed_matrix_t = MatrixCondensed<T,Condense>;
template<typename T,bool Condense=true>
using condensed_matrix_ptr_t = boost::shared_ptr<condensed_matrix_t<T,Condense>>;

//!
//! Create a shared pointer \p MatrixCondensed<T> from \p Args
//! @code
//! // create a MatrixCondensed boost::shared_ptr
//! auto mc = makeSharedMatrixCondensed<double>();
//! @endcode
//!
template< class T, bool C, class... Args >
condensed_matrix_ptr_t<T,C>
makeSharedMatrixCondensed( Args&&... args )
{
    return boost::make_shared<MatrixCondensed<T,C>>( args... );
}

#if !defined(FEELPP_MATRIXCONDENSED_NOEXTERN)
extern template class MatrixCondensed<double,true>;
extern template class MatrixCondensed<double,false>;
//extern template class Backend<std::complex<double>>;
#endif

}


#endif
