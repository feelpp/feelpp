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
#define FEELPP_INSTANTIATE_MATRIXEIGENDENSE 1
#include <Eigen/Dense>
#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feelalg/vectoreigen.hpp>

namespace Feel
{

template <typename T>
MatrixEigenDense<T>::MatrixEigenDense()
    :
    super(),
    M_is_initialized( false ),
    M_is_closed( false ),
    M_mat()
{}
template <typename T>
MatrixEigenDense<T>::MatrixEigenDense( size_type r, size_type c, WorldComm const& worldComm )
    :
    super(worldComm),
    M_is_initialized( false ),
    M_is_closed( false ),
    M_mat( r, c )
{}

template <typename T>
MatrixEigenDense<T>::MatrixEigenDense( MatrixEigenDense const & m )
    :
    super( m ),
    M_is_initialized( m.M_is_initialized ),
    M_is_closed( m.M_is_closed ),
    M_mat( m.M_mat )

{}

template <typename T>
MatrixEigenDense<T>::~MatrixEigenDense()
{}

template <typename T>
void
MatrixEigenDense<T>::init ( const size_type m,
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
MatrixEigenDense<T>::init ( const size_type m,
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
    this->zero ();
}


template<typename T>
void
MatrixEigenDense<T>::addMatrix ( int* rows, int nrows,
                                 int* cols, int ncols,
                                 value_type* data )
{
    for( int i=0; i < nrows; ++i )
        for( int j=0; j < ncols; ++j )
        {
            M_mat( rows[i], cols[j] ) += data[i*ncols+j];
        }

}
template<typename T>
void
MatrixEigenDense<T>::scale( const T a )
{
    M_mat *= a;
}

template<typename T>
void
MatrixEigenDense<T>::resize( size_type nr, size_type nc, bool /*preserve*/ )
{
    M_mat.resize( nr, nc );
}

template<typename T>
void
MatrixEigenDense<T>::close() const
{
    M_is_closed = true;
}

template<typename T>
void
MatrixEigenDense<T>::transpose( MatrixSparse<value_type>& Mt, size_type options ) const
{
    FEELPP_ASSERT( 0 ).warn( "not implemented yet" );
}


template<typename T>
void
MatrixEigenDense<T>::diagonal ( Vector<T>& dest ) const
{
#if 0
    VectorEigen<T>& _dest( dynamic_cast<VectorEigen<T>&>( dest ) );
    _dest = M_mat.diagonal();
#endif
}
template<typename T>
typename MatrixEigenDense<T>::real_type
MatrixEigenDense<T>::energy( Vector<value_type> const& __v,
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
MatrixEigenDense<T>::addMatrix( value_type v, MatrixSparse<value_type>& _m )
{
    MatrixEigenDense<value_type>* m = dynamic_cast<MatrixEigenDense<value_type>*>( &_m );
    FEELPP_ASSERT( m != 0 ).error( "invalid sparse matrix type, should be MatrixEigenDense" );
    //FEELPP_ASSERT( m->closed() ).error( "invalid sparse matrix type, should be closed" );

    if ( !m )
        throw std::invalid_argument( "m" );

    if ( !this->closed() )
    {
        M_mat += v * m->M_mat;
    }
}

template<typename T>
void
MatrixEigenDense<T>::updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m,
                                     std::vector<size_type> start_i,
                                     std::vector<size_type> start_j )
{
    LOG(ERROR) << "invalid call to updateBlockMat(), not yet implemented\n";

}

template<typename T>
void
MatrixEigenDense<T>::printMatlab( const std::string filename ) const
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
    file_out << varName  << " = [ ";
    file_out.precision( 16 );
    file_out.setf( std::ios::scientific );

    for ( size_type i = 0; i < M_mat.rows(); ++i )
    {
        for ( size_type j = 0; j < M_mat.cols(); ++j )
        {
            value_type v = M_mat( i,j );

            file_out << i + 1 << separator
                     << j + 1 << separator
                     << v  << std::endl;
        }
    }

    file_out << "];" << std::endl;
    //file_out << "I=S(:,1); J=S(:,2); S=S(:,3);" << std::endl;
    file_out << "I="<<varName<<"(:,1); J="<<varName<<"(:,2);" << std::endl;
    file_out << "spy("<<varName<<");" << std::endl;
}

template<typename T>
void
MatrixEigenDense<T>::zeroRows( std::vector<int> const& rows,
                               Vector<value_type> const& vals,
                               Vector<value_type>& rhs,
                               Context const& on_context )
{
    Feel::detail::ignore_unused_variable_warning( rhs );
    Feel::detail::ignore_unused_variable_warning( vals );

    VectorEigen<T>* prhs = dynamic_cast<VectorEigen<T>*> ( &rhs );
    VectorEigen<T> const* pvals = dynamic_cast<const VectorEigen<T>*> ( &vals );

    for ( size_type i = 0; i < rows.size(); ++i )
    {
        value_type value = 1.0;

        if ( on_context.test( ContextOn::ELIMINATION|ContextOn::KEEP_DIAGONAL ) )
            value = M_mat( rows[i], rows[i] );
        M_mat.row( rows[i] ).setZero();

        prhs->vec() -= M_mat.col(rows[i])*vals(rows[i]);
        M_mat.col( rows[i] ).setZero();

        // set diagonal
        M_mat( rows[i], rows[i] ) = value;

        // multiply rhs by value of the diagonal entry value
        rhs.set( rows[i], value * vals(rows[i]) );
    }
}

template<typename T>
void
MatrixEigenDense<T>::sqrt( MatrixSparse<value_type>& _m ) const
{
    MatrixEigenDense<value_type>* _md = dynamic_cast<MatrixEigenDense<value_type>*>(&_m);
    _md->mat() = this->M_mat.sqrt();
}

template <typename T>
void
MatrixEigenDense<T>::eigenValues ( std::vector<std::complex<double>> &Eingvs )
{
    auto eigen_vals = M_mat.eigenvalues();
    for (size_type i=0; i < eigen_vals.size(); ++i )
    {
        Eingvs.push_back(eigen_vals[i]);
    }
}


template <typename T>
void
MatrixEigenDense<T>::matMatMult ( MatrixSparse<T> const& In, MatrixSparse<T> &Res )
{
    //FEELPP_ASSERT ( this->isInitialized() ).error( "eigendense matrix not initialized" );

    FEELPP_ASSERT( this->size2() == In.size1() )( this->size2() )( In.size1() ).error( "incompatible dimension" );

    MatrixEigenDense<T> const* X = dynamic_cast<MatrixEigenDense<T> const*> ( &In );
    MatrixEigenDense<T>* Y = dynamic_cast<MatrixEigenDense<T>*> ( &Res );

    FEELPP_ASSERT ( X != 0 ).error( "invalid eigendense matrix" );

    Y->mat() = this->M_mat*X->mat();

    // int ierr=0;
    // ierr = MatMatMult(this->M_mat, X->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y->mat());
    // CHKERRABORT( this->comm(),ierr );

}

template <typename T>
void
MatrixEigenDense<T>::matInverse ( MatrixSparse<T> &Inv )
{
    //FEELPP_ASSERT ( this->isInitialized() ).error( "eigendense matrix not initialized" );

    FEELPP_ASSERT( this->size1() == Inv.size1() )( this->size1() )( Inv.size1() ).error( "incompatible dimension" );

    //MatrixEigenDense<T> const* X = dynamic_cast<MatrixEigenDense<T> const*> ( &In );
    MatrixEigenDense<T>* X = dynamic_cast<MatrixEigenDense<T>*> ( &Inv );

    FEELPP_ASSERT ( X != 0 ).error( "invalid eigendense matrix" );

    X->mat() = this->M_mat.inverse();

    // int ierr=0;
    // ierr = MatMatMult(this->M_mat, X->mat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y->mat());
    // CHKERRABORT( this->comm(),ierr );

}


template<typename T>
void
MatrixEigenDense<T>::applyInverseSqrt( Vector<T>& vec_in, Vector<T>& vec_out )
{

    VectorEigen<T>* _rhs = dynamic_cast<VectorEigen<T>*> ( &vec_in );
    VectorEigen<T>* res = dynamic_cast<VectorEigen<T>*> ( &vec_out );


    std::vector<std::complex<double>> eigen_vals;
    this->eigenValues(eigen_vals);

    std::sort(eigen_vals.begin(), eigen_vals.end(),
              [&]( std::complex<double> const& x, std::complex<double> const& y )->bool{
                  return real(x) < real(y); } );

    double m = real(eigen_vals.front());
    double M = real(eigen_vals.back());

    double k = (std::pow(M/m,1./4)-1)/(std::pow(M/m,1./4)+1);
    double L = -std::log(k)/pi;

    double K, Kp;
    math::ellipkkp(L, K, Kp);

    size_type N = 15;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> AA;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> Amb(N);
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> Am1b(N);
    std::vector<std::complex<double>> t, sn ,cn, dn;

    for (double dt = 0.5; dt < N; ++dt)
    {
        t.push_back(std::complex<double>(0,0.5*Kp) - std::complex<double>(K,0) + std::complex<double>(dt*2*K/N,0));
    }

    math::ellipjc(t,L,sn,cn,dn);

    for (size_type it = 0; it < N; ++it)
    {
        std::complex<double> w = math::pow(m*M,0.25)*((std::complex<double>(1./k,0)+sn[it])/(std::complex<double>(1./k,0)-sn[it]));
        std::complex<double> dzdt = cn[it]*dn[it]/(std::pow(std::complex<double>(1./k,0)-sn[it],2.));

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> A = this->M_mat.template cast<std::complex<double>>();
        AA = A;
        A *= std::complex<double>(-1.,0);

        for (size_type iti = 0; iti < A.rows(); ++iti)
        {
            A(iti,iti) = std::pow(w,2.)-(this->M_mat).operator()(iti,iti);
        }

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> rhs = _rhs->vec().template cast<std::complex<double>>();

        Am1b = (A.inverse())*rhs;
        auto fw = std::pow(w,-0.5);
        Amb += Am1b*dzdt*(std::pow(fw,2.)/w);
    }

    Amb = -8*K*std::pow(m*M,1./4)*AA*Amb.imag()/(k*pi*N);

    for (size_type it = 0; it < Amb.size(); ++it)
    {
        res->set(it, Amb.real().operator()(it));
    }
}


template <typename T>
void
MatrixEigenDense<T>::createSubmatrix( MatrixSparse<T>& submatrix,
                                 const std::vector<size_type>& rows,
                                 const std::vector<size_type>& cols ) const
{
    MatrixEigenDense<T>* A = dynamic_cast<MatrixEigenDense<T>*> ( &submatrix );

    // eigen can not extract non-contiguous submatrices, hence we cannot use block
    for (size_type i=0; i<rows.size(); ++i)
    {
        for (size_type j=0; j<rows.size(); ++j)
        {
            A->mat()(i,j) = this->M_mat(rows[i],cols[j]);
        }
    }
}

//
// Explicit instantiations
//
template class MatrixEigenDense<double>;
template class MatrixEigenDense<std::complex<double>>;
//template class MatrixEigenDense<double,gmm::col_major>;
}
