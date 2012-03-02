/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-15

  Copyright (C) 2007, 2009 Université Joseph Fourier (Grenoble I)

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
   \file matrixgmm.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-15
 */
#include <feel/feelalg/matrixgmm.hpp>

namespace Feel
{

template <typename T, typename LayoutType>
MatrixGmm<T,LayoutType>::MatrixGmm()
    :
    super(),
    _M_is_initialized( false ),
    _M_is_closed( false ),
    _M_mat(),
    _M_wmat()
{}
template <typename T, typename LayoutType>
MatrixGmm<T,LayoutType>::MatrixGmm( size_type r, size_type c )
    :
    super(),
    _M_is_initialized( false ),
    _M_is_closed( false ),
    _M_mat( r, c ),
    _M_wmat( r, c )
{}

template <typename T, typename LayoutType>
MatrixGmm<T,LayoutType>::MatrixGmm( MatrixGmm const & m )
    :
    super( m ),
    _M_is_initialized( m._M_is_initialized ),
    _M_is_closed( m._M_is_closed ),
    _M_mat( m._M_mat ),
    _M_wmat( m._M_wmat )

{}

template <typename T, typename LayoutType>
MatrixGmm<T,LayoutType>::~MatrixGmm()
{}

template <typename T, typename LayoutType>
void
MatrixGmm<T,LayoutType>::init (const size_type m,
                               const size_type n,
                               const size_type /*m_l*/,
                               const size_type /*n_l*/,
                               const size_type /*nnz*/,
                               const size_type /*noz*/)
{
    if ((m==0) || (n==0))
        return;

    _M_wmat.resize(m,n);
    this->zero ();
}
template <typename T, typename LayoutType>
void
MatrixGmm<T,LayoutType>::init ( const size_type m,
                                const size_type n,
                                const size_type m_l,
                                const size_type n_l,
                                graph_ptrtype const& graph )
{
    Feel::detail::ignore_unused_variable_warning(m);
    Feel::detail::ignore_unused_variable_warning(n);
    Feel::detail::ignore_unused_variable_warning(m_l);
    Feel::detail::ignore_unused_variable_warning(n_l);
    Feel::detail::ignore_unused_variable_warning(graph);

    if ((m==0) || (n==0))
        return;

    _M_wmat.resize(m,n);
    this->zero ();
}

template<typename T, typename LayoutType>
void
MatrixGmm<T, LayoutType>::fill( pattern_type const& /*__pattern*/ )
{
}

template<typename T, typename LayoutType>
void
MatrixGmm<T, LayoutType>::resize( size_type nr, size_type nc, bool /*preserve*/ )
{
    gmm::resize( _M_wmat, nr, nc );
}

template<typename T, typename LayoutType>
void
MatrixGmm<T, LayoutType>::close() const
{
    Debug(7015) << "[MatrixGmm<T, LayoutType>::close()] nr = " << _M_mat.nr << "\n";
    Debug(7015) << "[MatrixGmm<T, LayoutType>::close()] nc = " << _M_mat.nc << "\n";
#if 1
    //if (_M_mat.pr) { delete[] pr; delete[] ir; delete[] jc; }
    //if (_M_mat.pr) delete[] _M_mat.pr;
    _M_mat.pr = 0;
    //delete[] _M_mat.ir;
    _M_mat.ir = 0;
    //delete[] _M_mat.jc;
    _M_mat.jc = 0;
#endif
    _M_mat.init_with_good_format( _M_wmat );
    _M_is_closed = true;
    // release memory of write optimized matrix
    gmm::resize( _M_wmat, 0, 0 );
}

template<typename T, typename LayoutType>
void
MatrixGmm<T, LayoutType>::transpose( MatrixSparse<value_type>& Mt ) const
{
    FEELPP_ASSERT( 0 ).warn( "not implemented yet" );
}



template<typename T, typename LayoutType>
typename MatrixGmm<T, LayoutType>::value_type
MatrixGmm<T, LayoutType>::energy( Vector<value_type> const& __v,
                                  Vector<value_type> const& __u,
                                  bool tranpose ) const
{
    VectorUblas<T> const& v( dynamic_cast<VectorUblas<T> const&>( __v ) );
    VectorUblas<T> const& u( dynamic_cast<VectorUblas<T> const&>( __u ) );

    std::vector<value_type> __v1( __u.size() ), __u1( __u.size() ), __t( __u.size() );
    std::copy( v.vec().begin(), v.vec().end(), __v1.begin() );
    std::copy( u.vec().begin(), u.vec().end(), __u1.begin() );
    if ( !tranpose )
        gmm::mult( _M_mat, __u1, __t );
    else
        gmm::mult( gmm::transposed( _M_mat ), __u1, __t );

    return gmm::vect_sp( __v1, __t );
}

template<typename T, typename LayoutType>
void
MatrixGmm<T, LayoutType>::addMatrix(value_type v, MatrixSparse<value_type>& _m )
{
    MatrixGmm<value_type, LayoutType>* m = dynamic_cast<MatrixGmm<value_type,LayoutType>*>(&_m);
    FEELPP_ASSERT( m != 0 ).error( "invalid sparse matrix type, should be MatrixGmm" );
    FEELPP_ASSERT( m->closed() ).error( "invalid sparse matrix type, should be closed" );
    if ( !m )
        throw std::invalid_argument("m");
    if ( !this->closed() )
        // matrices can have a different pattern
        gmm::add( gmm::scaled(m->mat(),v), this->wmat(), this->wmat() );
    //else
    // be careful : matrices should have the _same_ pattern here
    //gmm::add( gmm::scaled(m->mat(),v), this->mat(), this->mat() );
}


template<typename T, typename LayoutType>
void
MatrixGmm<T, LayoutType>::printMatlab(const std::string filename ) const
{
    std::string name = filename;
    std::string separator = " , ";

    // check on the file name
    int i = filename.find( "." );

    if ( i <= 0 )
        name = filename + ".m";
    else
        {
            if ( ( size_type ) i != filename.size() - 2 ||
                 filename[ i + 1 ] != 'm' )
                {
                    std::cerr << "Wrong file name ";
                    name = filename + ".m";
                }
        }

    std::ofstream file_out( name.c_str() );

    FEELPP_ASSERT( file_out)( filename ).error("[Feel::spy] ERROR: File cannot be opened for writing.");

    file_out << "S = [ ";
    file_out.precision( 16 );
    file_out.setf( std::ios::scientific );
    for (size_type i = 0; i < gmm::mat_nrows(_M_mat); ++i)
        {
            for (size_type j = 0; j < gmm::mat_ncols(_M_mat); ++j)
                {
                    value_type v = _M_mat(i,j);
                    if ( v != typename gmm::linalg_traits<matrix_type>::value_type(0))
                        {
                            file_out << i + 1 << separator
                                     << j + 1 << separator
                                     << v  << std::endl;
                        }
                }
        }
    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3);" << std::endl;
    file_out << "A=sparse(I,J,S); spy(A);" << std::endl;
}


//
// Explicit instantiations
//
template class MatrixGmm<double,gmm::row_major>;
//template class MatrixGmm<double,gmm::col_major>;
}
namespace gmm
{
template class linalg_traits<ublas::vector<double> >;
//template class linalg_traits<ublas::vector<long double> >;
//template class linalg_traits<ublas::vector<dd_real> >;
}
