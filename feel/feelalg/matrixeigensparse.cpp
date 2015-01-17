/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-15

  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-15
 */
#define FEELPP_INSTANTIATE_MATRIXEIGENSPARSE
#include <boost/unordered_map.hpp>

#include <feel/feelalg/matrixeigensparse.hpp>

namespace Feel
{

template <typename T>
MatrixEigenSparse<T>::MatrixEigenSparse()
    :
    super(),
    M_is_initialized( false ),
    M_is_closed( false ),
    M_mat()
{}
template <typename T>
MatrixEigenSparse<T>::MatrixEigenSparse( size_type r, size_type c, WorldComm const& worldComm )
    :
    super(worldComm),
    M_is_initialized( false ),
    M_is_closed( false ),
    M_mat( r, c )
{}

template <typename T>
MatrixEigenSparse<T>::MatrixEigenSparse( MatrixEigenSparse const & m )
    :
    super( m ),
    M_is_initialized( m.M_is_initialized ),
    M_is_closed( m.M_is_closed ),
    M_mat( m.M_mat )

{}

template <typename T>
MatrixEigenSparse<T>::~MatrixEigenSparse()
{}

template <typename T>
void
MatrixEigenSparse<T>::init ( const size_type m,
                                const size_type n,
                                const size_type /*m_l*/,
                                const size_type /*n_l*/,
                                const size_type /*nnz*/,
                                const size_type /*noz*/ )
{
    if ( ( m==0 ) || ( n==0 ) )
        return;

    M_mat.resize( m,n );
    this->zero ();
}
template <typename T>
void
MatrixEigenSparse<T>::init ( const size_type m,
                                const size_type n,
                                const size_type m_l,
                                const size_type n_l,
                                graph_ptrtype const& graph )
{
    Feel::detail::ignore_unused_variable_warning( m );
    Feel::detail::ignore_unused_variable_warning( n );
    Feel::detail::ignore_unused_variable_warning( m_l );
    Feel::detail::ignore_unused_variable_warning( n_l );
    Feel::detail::ignore_unused_variable_warning( graph );

    if ( ( m==0 ) || ( n==0 ) )
        return;

    M_mat.resize( m,n );
    M_tripletList.reserve(graph->ja().size());

}


template<typename T>
void
MatrixEigenSparse<T>::addMatrix ( int* rows, int nrows,
                                 int* cols, int ncols,
                                 value_type* data )
{

    for( int i=0; i < nrows; ++i )
        for( int j=0; j < ncols; ++j )
        {
            M_tripletList.push_back(triplet(rows[i], cols[j], data[i*ncols+j]) );
        }

}
template<typename T>
void
MatrixEigenSparse<T>::scale( const T a )
{
    M_mat *= a;
}

template<typename T>
void
MatrixEigenSparse<T>::resize( size_type nr, size_type nc, bool /*preserve*/ )
{
    M_mat.resize( nr, nc );
}

template<typename T>
void
MatrixEigenSparse<T>::close() const
{

    if ( !M_is_closed )
    {
        LOG(INFO) << "Closing matrix";
        M_mat.setFromTriplets(M_tripletList.begin(), M_tripletList.end());
        M_mat.makeCompressed();
        std::vector<triplet>().swap(M_tripletList);
        M_is_closed = true;
    }
}

template<typename T>
void
MatrixEigenSparse<T>::transpose( MatrixSparse<value_type>& Mt, size_type options ) const
{
    if(M_is_closed) {
        MatrixEigenSparse<T>* Atrans = dynamic_cast<MatrixEigenSparse<T>*>(&Mt);
        Atrans->M_mat = M_mat.transpose().eval();
        Atrans->M_is_closed = Atrans->M_is_initialized = true;
    }
}


template<typename T>
void
MatrixEigenSparse<T>::diagonal ( Vector<T>& dest ) const
{
#if 0
    VectorEigen<T>& _dest( dynamic_cast<VectorEigen<T>&>( dest ) );
    _dest = M_mat.diagonal();
#endif
}
template<typename T>
typename MatrixEigenSparse<T>::real_type
MatrixEigenSparse<T>::energy( Vector<value_type> const& __v,
                             Vector<value_type> const& __u,
                             bool tranpose ) const
{
#if 0
    VectorUblas<T> const& v( dynamic_cast<VectorUblas<T> const&>( __v ) );
    VectorUblas<T> const& u( dynamic_cast<VectorUblas<T> const&>( __u ) );

    std::vector<value_type> __v1( __u.size() ), __u1( __u.size() ), __t( __u.size() );
    std::copy( v.vec().begin(), v.vec().end(), __v1.begin() );
    std::copy( u.vec().begin(), u.vec().end(), __u1.begin() );

    if ( !tranpose )
        gmm::mult( M_mat, __u1, __t );

    else
        gmm::mult( gmm::transposed( M_mat ), __u1, __t );

    return gmm::vect_sp( __v1, __t );
#endif
    return 0;
}

template<typename T>
void
MatrixEigenSparse<T>::addMatrix( value_type v, MatrixSparse<value_type>& _m )
{
    MatrixEigenSparse<value_type>* m = dynamic_cast<MatrixEigenSparse<value_type>*>( &_m );
    FEELPP_ASSERT( m != 0 ).error( "invalid sparse matrix type, should be MatrixEigenSparse" );
    FEELPP_ASSERT( m->closed() ).error( "invalid sparse matrix type, should be closed" );

    if ( !m )
        throw std::invalid_argument( "m" );

    if ( !this->closed() )
    {
        M_mat += v * m->M_mat;
    }
}
template<typename T>
void
MatrixEigenSparse<T>::updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m,
                                      std::vector<size_type> start_i,
                                      std::vector<size_type> start_j )
{
    LOG(ERROR) << "Invalid call to updateBlockMat, not yet implemented\n";
}

template<typename T>
void
MatrixEigenSparse<T>::zeroRows( std::vector<int> const& rows,
                                Vector<value_type> const& vals,
                                Vector<value_type>& rhs,
                                Context const& on_context )
{
    LOG(INFO) << "zero out " << rows.size() << " rows except diagonal is row major: " << M_mat.IsRowMajor;
    Feel::detail::ignore_unused_variable_warning( rhs );
    Feel::detail::ignore_unused_variable_warning( vals );
    boost::unordered_map<int,std::set<int>> m;
    for (int k=0; k<rows.size(); ++k)
    {
        for (typename matrix_type::InnerIterator it(M_mat,rows[k]); it; ++it)
        {
            m[it.row()].insert(it.col());
            value_type value = 1.0;
            if ( on_context.test( ContextOn::ELIMINATION|ContextOn::KEEP_DIAGONAL ) )
                value = it.value();
            rhs.add( it.row(), -it.value() * vals(rows[k]) );
            it.valueRef() = 0;

            if ( it.row() == it.col() )
            {
                it.valueRef() = value;
                rhs.set( it.row(), value * vals(rows[k]) );
            }
        }
    }
    for(auto rit = m.begin();rit != m.end();++rit )
    {
        for (typename matrix_type::InnerIterator it(M_mat,rit->first); it; ++it)
        {
            double value = 1.0;
            if( rit->second.find( it.row() ) != rit->second.end() )
            {
                // don't change diagonal, it was done in the first pass
                if ( it.row() != it.col() )
                    it.valueRef() = 0;
            }
        }
    }
}

template<typename T>
void
MatrixEigenSparse<T>::printMatlab( const std::string filename ) const
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

    FEELPP_ASSERT( file_out )( filename ).error( "[Feel::spy] ERROR: File cannot be opened for writing." );

		std::string varName = "var_" + filename.substr(0,filename.find("."));
    file_out << varName << " = [ ";
    file_out.precision( 16 );
    file_out.setf( std::ios::scientific );

    for (int k=0; k<M_mat.outerSize(); ++k)
    {
        for (typename matrix_type::InnerIterator it(M_mat,k); it; ++it)
        {
            value_type v = it.value();

            file_out << it.row() + 1 << separator
                     << it.col() + 1 << separator
                     << v  << std::endl;
        }
    }

    file_out << "];" << std::endl;
    file_out << "I="<<varName<<"(:,1); J="<<varName<<"(:,2); "<<varName<<"="<<varName<<"(:,3);" << std::endl;
    file_out << "spy("<<varName<<");" << std::endl;
}


//
// Explicit instantiations
//
template class MatrixEigenSparse<double>;
template class MatrixEigenSparse<std::complex<double>>;
}
