/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-04-26

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_SERIALIZATION_HPP
#define FEELPP_SERIALIZATION_HPP 1

#include <boost/multi_array.hpp>
#include <boost/detail/identifier.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include <feel/feelcore/disablewarnings.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include <feel/feelcore/reenablewarnings.hpp>

namespace boost
{
namespace serialization
{
template<class Archive>
void load( Archive & ar,
           boost::multi_array<double,2> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<double,2> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    size_ n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( boost::extents[n0][n1] );
    ar >> make_array( t.data(), t.num_elements() );
}
template<typename Archive>
void save( Archive & ar,
           const boost::multi_array<double,2> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<double,2> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0 = ( t.shape()[0] );
    ar << BOOST_SERIALIZATION_NVP( n0 );
    size_ n1 = ( t.shape()[1] );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(),
                                            t.num_elements() );
}
template<class Archive>
void serialize( Archive & ar,
                boost::multi_array<double,2>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


// 3
template<class Archive>
void load( Archive & ar,
           boost::multi_array<double,3> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<double,3> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    size_ n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    size_ n2;
    ar >> BOOST_SERIALIZATION_NVP( n2 );
    t.resize( boost::extents[n0][n1][n2] );
    ar >> make_array( t.data(), t.num_elements() );
}
template<typename Archive>
void save( Archive & ar,
           const boost::multi_array<double,3> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<double,3> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0 = ( t.shape()[0] );
    ar << BOOST_SERIALIZATION_NVP( n0 );
    size_ n1 = ( t.shape()[1] );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    size_ n2 = ( t.shape()[2] );
    ar << BOOST_SERIALIZATION_NVP( n2 );
    ar << boost::serialization::make_array( t.data(),
                                            t.num_elements() );
}
template<class Archive>
void serialize( Archive & ar,
                boost::multi_array<double,3>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

// 4
template<class Archive>
void load( Archive & ar,
           boost::multi_array<double,4> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<double,4> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    size_ n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    size_ n2;
    ar >> BOOST_SERIALIZATION_NVP( n2 );
    size_ n3;
    ar >> BOOST_SERIALIZATION_NVP( n3 );
    t.resize( boost::extents[n0][n1][n2][n3] );
    ar >> make_array( t.data(), t.num_elements() );
}
template<typename Archive>
void save( Archive & ar,
           const boost::multi_array<double,4> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<double,4> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0 = ( t.shape()[0] );
    ar << BOOST_SERIALIZATION_NVP( n0 );
    size_ n1 = ( t.shape()[1] );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    size_ n2 = ( t.shape()[2] );
    ar << BOOST_SERIALIZATION_NVP( n2 );
    size_ n3 = ( t.shape()[3] );
    ar << BOOST_SERIALIZATION_NVP( n3 );
    ar << boost::serialization::make_array( t.data(),
                                            t.num_elements() );
}
template<class Archive>
void serialize( Archive & ar,
                boost::multi_array<double,4>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


//
// MatrixXd
//
template<class Archive>
void load( Archive & ar,
           Eigen::MatrixXd & t,
           const unsigned int file_version )
{
    int n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    int n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( n0, n1 );
    ar >> make_array( t.data(), t.rows()*t.cols() );
}
template<typename Archive>
void save( Archive & ar,
           const Eigen::MatrixXd & t,
           const unsigned int file_version )
{
    int n0 = t.rows();
    ar << BOOST_SERIALIZATION_NVP( n0 );
    int n1 = t.cols();
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(),
                                            t.rows()*t.cols() );
}
template<class Archive>
void serialize( Archive & ar,
                Eigen::MatrixXd& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}
template<class Archive>
void load( Archive & ar,
           Eigen::MatrixXf & t,
           const unsigned int file_version )
{
    int n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    int n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( n0, n1 );
    ar >> make_array( t.data(), t.rows()*t.cols() );
}
template<typename Archive>
void save( Archive & ar,
           const Eigen::MatrixXf & t,
           const unsigned int file_version )
{
    int n0 = t.rows();
    ar << BOOST_SERIALIZATION_NVP( n0 );
    int n1 = t.cols();
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(),
                                            t.rows()*t.cols() );
}
template<class Archive>
void serialize( Archive & ar,
                Eigen::MatrixXf& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

//
// VectorXd
//
template<class Archive>
void load( Archive & ar,
           Eigen::VectorXd & t,
           const unsigned int file_version )
{
    int n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    t.resize( n0 );
    ar >> make_array( t.data(), t.size() );
}
template<typename Archive>
void save( Archive & ar,
           const Eigen::VectorXd & t,
           const unsigned int file_version )
{
    int n0 = t.size();
    ar << BOOST_SERIALIZATION_NVP( n0 );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<class Archive>
void serialize( Archive & ar,
                Eigen::VectorXd& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


//
// Matrix<N,M>
//
template<typename T, int N, int M, class Archive>
void load( Archive & ar,
           Eigen::Matrix<T,N,M> & t,
           const unsigned int file_version )
{
    int n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    int n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    ar >> make_array( t.data(), n0*n1 );
}
template<typename T, int N, int M, typename Archive>
void save( Archive & ar,
           const Eigen::Matrix<T,N,M> & t,
           const unsigned int file_version )
{
    int n0 = t.rows();
    int n1 = t.cols();
    ar << BOOST_SERIALIZATION_NVP( n0 );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(), n0*n1 );
}
template<typename T, int N, int M, class Archive>
void serialize( Archive & ar,
                Eigen::Matrix<T,N,M>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

//
// Matrix<RowsAtCompileTime,ColsAtCompileTime,Options,MaxRowsAtCompileTime,MaxColsAtCompileTime>
//
template<typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime, class Archive>
void load( Archive & ar,
           Eigen::Matrix<T,RowsAtCompileTime,ColsAtCompileTime,Options,MaxRowsAtCompileTime,MaxColsAtCompileTime> & t,
           const unsigned int file_version )
{
    int n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    int n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( n0,n1 );
    ar >> make_array( t.data(), n0*n1 );
}
template<typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime, typename Archive>
void save( Archive & ar,
           const Eigen::Matrix<T,RowsAtCompileTime,ColsAtCompileTime,Options,MaxRowsAtCompileTime,MaxColsAtCompileTime> & t,
           const unsigned int file_version )
{
    int n0 = t.rows();
    int n1 = t.cols();
    ar << BOOST_SERIALIZATION_NVP( n0 );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(), n0*n1 );
}
template<typename T, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime, class Archive>
void serialize( Archive & ar,
                Eigen::Matrix<T,RowsAtCompileTime,ColsAtCompileTime,Options,MaxRowsAtCompileTime,MaxColsAtCompileTime>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

template<class Archive>
void load( Archive & ar,
           boost::multi_array<Eigen::MatrixXd,2> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<Eigen::MatrixXd,2> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    size_ n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( boost::extents[n0][n1] );
    ar >> make_array( t.data(), t.num_elements() );
}
template<typename Archive>
void save( Archive & ar,
           const boost::multi_array<Eigen::MatrixXd,2> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<Eigen::MatrixXd,2> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0 = ( t.shape()[0] );
    ar << BOOST_SERIALIZATION_NVP( n0 );
    size_ n1 = ( t.shape()[1] );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(),
                                            t.num_elements() );
}
template<class Archive>
void serialize( Archive & ar,
                boost::multi_array<Eigen::MatrixXd,2>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

template<class Archive>
void load( Archive & ar,
           boost::multi_array<Eigen::VectorXd,2> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<Eigen::VectorXd,2> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    size_ n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( boost::extents[n0][n1] );
    ar >> make_array( t.data(), t.num_elements() );
}
template<typename Archive>
void save( Archive & ar,
           const boost::multi_array<Eigen::VectorXd,2> & t,
           const unsigned int file_version )
{
    typedef boost::multi_array<Eigen::VectorXd,2> multi_array_;
    typedef typename multi_array_::size_type size_;
    size_ n0 = ( t.shape()[0] );
    ar << BOOST_SERIALIZATION_NVP( n0 );
    size_ n1 = ( t.shape()[1] );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(),
                                            t.num_elements() );
}
template<class Archive>
void serialize( Archive & ar,
                boost::multi_array<Eigen::VectorXd,2>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


//
// Eigen::Tensor
//
template<typename T, class Archive>
void load( Archive & ar,
           Eigen::Tensor<T,1> & t,
           const unsigned int file_version )
{
    int n0,n1=1,n2=1;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    t.resize( n0 );
    ar >> make_array( t.data(), n0 );
}
template<typename T, class Archive>
void load( Archive & ar,
           Eigen::Tensor<T,2> & t,
           const unsigned int file_version )
{
    int n0,n1=1,n2=1;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    t.resize( n0, n1 );
    ar >> make_array( t.data(), n1*n0 );
}
template<typename T, class Archive>
void load( Archive & ar,
           Eigen::Tensor<T,3> & t,
           const unsigned int file_version )
{
    int n0,n1=1,n2=1;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    ar >> BOOST_SERIALIZATION_NVP( n2 );
    t.resize( n0, n1, n2 );
    ar >> make_array( t.data(), n1*n2*n0 );
}
template<typename T, int N, typename Archive>
void save( Archive & ar,
           const Eigen::Tensor<T,N> & t,
           const unsigned int file_version )
{
    int n0 = t.dimension(0);
    ar << BOOST_SERIALIZATION_NVP( n0 );
    if ( N >= 2 )
    {
        int n1 = t.dimension(1);
        ar << BOOST_SERIALIZATION_NVP( n1 );
    }
    if ( N >= 3 )
    {
        int n2 = t.dimension(2);
        ar << BOOST_SERIALIZATION_NVP( n2 );
    }

    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, int N,  class Archive>
void serialize( Archive & ar,
                Eigen::Tensor<T,N>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


//
// Eigen::TensorFixedSize
//
template<typename T, typename Archive>
void save( Archive & ar,
           const Eigen::TensorFixedSize<T,Eigen::Sizes<1,1>> & t,
           const unsigned int file_version )
{
    int m = t.dimension(0);
    int n = t.dimension(1);
    ar << BOOST_SERIALIZATION_NVP( m );
    ar << BOOST_SERIALIZATION_NVP( n );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, typename Archive>
void load( Archive & ar,
           Eigen::TensorFixedSize<T,Eigen::Sizes<1,1>> & t,
           const unsigned int file_version )
{
    int m;
    int n;
    ar >> BOOST_SERIALIZATION_NVP( m );
    ar >> BOOST_SERIALIZATION_NVP( n );
    DCHECK( m == t.dimension(0) ) << "invalid nulmber of rows: " << m << " should be " << t.dimension(0);
    DCHECK( n == t.dimension(1) ) << "invalid nulmber of cols: " << n << " should be " << t.dimension(1);
    ar >> boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, class Archive>
void serialize( Archive & ar,
                Eigen::TensorFixedSize<T,Eigen::Sizes<1,1>>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

template<typename T, typename Archive>
void save( Archive & ar,
           const Eigen::TensorFixedSize<T,Eigen::Sizes<1,2>> & t,
           const unsigned int file_version )
{
    int m = t.dimension(0);
    int n = t.dimension(1);
    ar << BOOST_SERIALIZATION_NVP( m );
    ar << BOOST_SERIALIZATION_NVP( n );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, typename Archive>
void load( Archive & ar,
           Eigen::TensorFixedSize<T,Eigen::Sizes<1,2>> & t,
           const unsigned int file_version )
{
    int m;
    int n;
    ar >> BOOST_SERIALIZATION_NVP( m );
    ar >> BOOST_SERIALIZATION_NVP( n );
    DCHECK( m == t.dimension(0) ) << "invalid nulmber of rows: " << m << " should be " << t.dimension(0);
    DCHECK( n == t.dimension(1) ) << "invalid nulmber of cols: " << n << " should be " << t.dimension(1);
    ar >> boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, class Archive>
void serialize( Archive & ar,
                Eigen::TensorFixedSize<T,Eigen::Sizes<1,2>>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


template<typename T, typename Archive>
void save( Archive & ar,
           const Eigen::TensorFixedSize<T,Eigen::Sizes<1,3>> & t,
           const unsigned int file_version )
{
    int m = t.dimension(0);
    int n = t.dimension(1);
    ar << BOOST_SERIALIZATION_NVP( m );
    ar << BOOST_SERIALIZATION_NVP( n );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, typename Archive>
void load( Archive & ar,
           Eigen::TensorFixedSize<T,Eigen::Sizes<1,3>> & t,
           const unsigned int file_version )
{
    int m;
    int n;
    ar >> BOOST_SERIALIZATION_NVP( m );
    ar >> BOOST_SERIALIZATION_NVP( n );
    DCHECK( m == t.dimension(0) ) << "invalid nulmber of rows: " << m << " should be " << t.dimension(0);
    DCHECK( n == t.dimension(1) ) << "invalid nulmber of cols: " << n << " should be " << t.dimension(1);
    ar >> boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, class Archive>
void serialize( Archive & ar,
                Eigen::TensorFixedSize<T,Eigen::Sizes<1,3>>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

template<typename T, typename Archive>
void save( Archive & ar,
           const Eigen::TensorFixedSize<T,Eigen::Sizes<2,1>> & t,
           const unsigned int file_version )
{
    int m = t.dimension(0);
    int n = t.dimension(1);
    ar << BOOST_SERIALIZATION_NVP( m );
    ar << BOOST_SERIALIZATION_NVP( n );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, typename Archive>
void load( Archive & ar,
           Eigen::TensorFixedSize<T,Eigen::Sizes<2,1>> & t,
           const unsigned int file_version )
{
    int m;
    int n;
    ar >> BOOST_SERIALIZATION_NVP( m );
    ar >> BOOST_SERIALIZATION_NVP( n );
    DCHECK( m == t.dimension(0) ) << "invalid nulmber of rows: " << m << " should be " << t.dimension(0);
    DCHECK( n == t.dimension(1) ) << "invalid nulmber of cols: " << n << " should be " << t.dimension(1);
    ar >> boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, class Archive>
void serialize( Archive & ar,
                Eigen::TensorFixedSize<T,Eigen::Sizes<2,1>>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

template<typename T, typename Archive>
void save( Archive & ar,
           const Eigen::TensorFixedSize<T,Eigen::Sizes<3,1>> & t,
           const unsigned int file_version )
{
    int m = t.dimension(0);
    int n = t.dimension(1);
    ar << BOOST_SERIALIZATION_NVP( m );
    ar << BOOST_SERIALIZATION_NVP( n );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, typename Archive>
void load( Archive & ar,
           Eigen::TensorFixedSize<T,Eigen::Sizes<3,1>> & t,
           const unsigned int file_version )
{
    int m;
    int n;
    ar >> BOOST_SERIALIZATION_NVP( m );
    ar >> BOOST_SERIALIZATION_NVP( n );
    DCHECK( m == t.dimension(0) ) << "invalid nulmber of rows: " << m << " should be " << t.dimension(0);
    DCHECK( n == t.dimension(1) ) << "invalid nulmber of cols: " << n << " should be " << t.dimension(1);
    ar >> boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, class Archive>
void serialize( Archive & ar,
                Eigen::TensorFixedSize<T,Eigen::Sizes<3,1>>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

template<typename T, int M, int N, typename Archive>
void save( Archive & ar,
           const Eigen::TensorFixedSize<T,Eigen::Sizes<M,N>> & t,
           const unsigned int file_version )
{
    int m = t.dimension(0);
    int n = t.dimension(1);
    ar << BOOST_SERIALIZATION_NVP( m );
    ar << BOOST_SERIALIZATION_NVP( n );
    ar << boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, int M, int N, typename Archive>
void load( Archive & ar,
           Eigen::TensorFixedSize<T,Eigen::Sizes<M,N>> & t,
           const unsigned int file_version )
{
    int m;
    int n;
    ar >> BOOST_SERIALIZATION_NVP( m );
    ar >> BOOST_SERIALIZATION_NVP( n );
    DCHECK( m == t.dimension(0) ) << "invalid nulmber of rows: " << m << " should be " << t.dimension(0);
    DCHECK( n == t.dimension(1) ) << "invalid nulmber of cols: " << n << " should be " << t.dimension(1);
    ar >> boost::serialization::make_array( t.data(),
                                            t.size() );
}
template<typename T, int M,  int N, class Archive>
void serialize( Archive & ar,
                Eigen::TensorFixedSize<T,Eigen::Sizes<M,N>>& t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


//
// boost::tuple<T1,T2>
//

template<typename T1, typename T2, class Archive>
void load( Archive & ar,
           boost::tuple<T1,T2> & t,
           const unsigned int file_version )
{
    ar >> BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
}
template<typename T1, typename T2, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2> const& t,
           const unsigned int file_version )
{
    ar << BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
}
template<typename T1, typename T2, class Archive>
void serialize( Archive & ar,
                boost::tuple<T1,T2> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}



//
// boost::tuple<T1,T2,T3>
//

template<typename T1, typename T2, typename T3, class Archive>
void load( Archive & ar,
           boost::tuple<T1,T2,T3> & t,
           const unsigned int file_version )
{
    ar >> BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
}
template<typename T1, typename T2, typename T3, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3> const& t,
           const unsigned int file_version )
{
    ar << BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
}
template<typename T1, typename T2, typename T3, class Archive>
void serialize( Archive & ar,
                boost::tuple<T1,T2,T3> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

//
// boost::tuple<T1,T2,T3,T4>
//

template<typename T1, typename T2, typename T3, typename T4, class Archive>
void load( Archive & ar,
           boost::tuple<T1,T2,T3,T4> & t,
           const unsigned int file_version )
{
    ar >> BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<3>( t ) );
}
template<typename T1, typename T2, typename T3, typename T4, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3,T4> const& t,
           const unsigned int file_version )
{
    ar << BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<3>( t ) );
}
template<typename T1, typename T2, typename T3, typename T4, class Archive>
void serialize( Archive & ar,
                boost::tuple<T1,T2,T3,T4> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

//
// boost::tuple<T1,T2,T3,T4,T5>
//

template<typename T1, typename T2, typename T3, typename T4, typename T5, class Archive>
void load( Archive & ar,
           boost::tuple<T1,T2,T3,T4,T5> & t,
           const unsigned int file_version )
{
    ar >> BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<3>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<4>( t ) );
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3,T4,T5> const& t,
           const unsigned int file_version )
{
    ar << BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<3>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<4>( t ) );
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, class Archive>
void serialize( Archive & ar,
                boost::tuple<T1,T2,T3,T4,T5> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

//
// boost::tuple<T1,T2,T3,T4,T5,T6>
//

template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, class Archive>
void load( Archive & ar,
           boost::tuple<T1,T2,T3,T4,T5,T6> & t,
           const unsigned int file_version )
{
    ar >> BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<3>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<4>( t ) );
    ar >> BOOST_SERIALIZATION_NVP( boost::get<5>( t ) );
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3,T4,T5,T6> const& t,
           const unsigned int file_version )
{
    ar << BOOST_SERIALIZATION_NVP( boost::get<0>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<1>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<2>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<3>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<4>( t ) );
    ar << BOOST_SERIALIZATION_NVP( boost::get<5>( t ) );
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, class Archive>
void serialize( Archive & ar,
                boost::tuple<T1,T2,T3,T4,T5,T6> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}


//
// boost::detail::identifier<T, D>
//
template<typename T, typename D, class Archive>
void load( Archive & ar,
           boost::detail::identifier<T, D> & t,
           const unsigned int file_version )
{
    T value;
    ar >> BOOST_SERIALIZATION_NVP( value );
    t.assign( value );
}
template<typename T, typename D, class Archive>
void save( Archive & ar,
           boost::detail::identifier<T, D> const& t,
           const unsigned int file_version )
{
    T value = t.value();
    ar << BOOST_SERIALIZATION_NVP( value );
}
template<typename T, typename D, class Archive>
void serialize( Archive & ar,
                boost::detail::identifier<T, D> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}

} // serialization
} //boost

#endif /* __Serialization_H */
