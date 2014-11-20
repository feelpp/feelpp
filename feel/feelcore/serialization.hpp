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
/**
   \file serialization.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-04-26
 */
#ifndef __Serialization_H
#define __Serialization_H 1

#include <boost/multi_array.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <Eigen/Core>

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
template<int N, int M, class Archive>
void load( Archive & ar,
           Eigen::Matrix<double,N,M> & t,
           const unsigned int file_version )
{
    int n0;
    ar >> BOOST_SERIALIZATION_NVP( n0 );
    int n1;
    ar >> BOOST_SERIALIZATION_NVP( n1 );
    ar >> make_array( t.data(), n0*n1 );
}
template<int N, int M, typename Archive>
void save( Archive & ar,
           const Eigen::Matrix<double,N,M> & t,
           const unsigned int file_version )
{
    int n0 = t.rows();
    int n1 = t.cols();
    ar << BOOST_SERIALIZATION_NVP( n0 );
    ar << BOOST_SERIALIZATION_NVP( n1 );
    ar << boost::serialization::make_array( t.data(), n0*n1 );
}
template<int N, int M, class Archive>
void serialize( Archive & ar,
                Eigen::Matrix<double,N,M>& t,
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
// boost::tuple<T1,T2>
//

template<typename T1, typename T2, class Archive>
void load( Archive & ar,
           boost::tuple<T1,T2> & t,
           const unsigned int file_version )
{
    T1 x;
    ar >> BOOST_SERIALIZATION_NVP( x );
    T2 y;
    ar >> BOOST_SERIALIZATION_NVP( y );
    t = boost::make_tuple( x,y );
}
template<typename T1, typename T2, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2> const& t,
           const unsigned int file_version )
{
    T1 x;
    T2 y;
    boost::tie( x,y ) = t;
    ar << BOOST_SERIALIZATION_NVP( x );
    ar << BOOST_SERIALIZATION_NVP( y );
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
    T1 x;
    ar >> BOOST_SERIALIZATION_NVP( x );
    T2 y;
    ar >> BOOST_SERIALIZATION_NVP( y );
    T3 z;
    ar >> BOOST_SERIALIZATION_NVP( z );
    t = boost::make_tuple( x,y,z );
}
template<typename T1, typename T2, typename T3, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3> const& t,
           const unsigned int file_version )
{
    T1 x;
    T2 y;
    T3 z;
    boost::tie( x,y,z ) = t;
    ar << BOOST_SERIALIZATION_NVP( x );
    ar << BOOST_SERIALIZATION_NVP( y );
    ar << BOOST_SERIALIZATION_NVP( z );

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
    T1 x1;
    ar >> BOOST_SERIALIZATION_NVP( x1 );
    T2 x2;
    ar >> BOOST_SERIALIZATION_NVP( x2 );
    T3 x3;
    ar >> BOOST_SERIALIZATION_NVP( x3 );
    T4 x4;
    ar >> BOOST_SERIALIZATION_NVP( x4 );
    t = boost::make_tuple( x1,x2,x3,x4 );
}
template<typename T1, typename T2, typename T3, typename T4, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3,T4> const& t,
           const unsigned int file_version )
{
    T1 x1;
    T2 x2;
    T3 x3;
    T4 x4;
    boost::tie( x1,x2,x3,x4 ) = t;
    ar << BOOST_SERIALIZATION_NVP( x1 );
    ar << BOOST_SERIALIZATION_NVP( x2 );
    ar << BOOST_SERIALIZATION_NVP( x3 );
    ar << BOOST_SERIALIZATION_NVP( x4 );

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
    T1 x1;
    ar >> BOOST_SERIALIZATION_NVP( x1 );
    T2 x2;
    ar >> BOOST_SERIALIZATION_NVP( x2 );
    T3 x3;
    ar >> BOOST_SERIALIZATION_NVP( x3 );
    T4 x4;
    ar >> BOOST_SERIALIZATION_NVP( x4 );
    T5 x5;
    ar >> BOOST_SERIALIZATION_NVP( x5 );
    t = boost::make_tuple( x1,x2,x3,x4,x5 );
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3,T4,T5> const& t,
           const unsigned int file_version )
{
    T1 x1;
    T2 x2;
    T3 x3;
    T4 x4;
    T5 x5;
    boost::tie( x1,x2,x3,x4,x5 ) = t;
    ar << BOOST_SERIALIZATION_NVP( x1 );
    ar << BOOST_SERIALIZATION_NVP( x2 );
    ar << BOOST_SERIALIZATION_NVP( x3 );
    ar << BOOST_SERIALIZATION_NVP( x4 );
    ar << BOOST_SERIALIZATION_NVP( x5 );

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
    T1 x1;
    ar >> BOOST_SERIALIZATION_NVP( x1 );
    T2 x2;
    ar >> BOOST_SERIALIZATION_NVP( x2 );
    T3 x3;
    ar >> BOOST_SERIALIZATION_NVP( x3 );
    T4 x4;
    ar >> BOOST_SERIALIZATION_NVP( x4 );
    T5 x5;
    ar >> BOOST_SERIALIZATION_NVP( x5 );
    T6 x6;
    ar >> BOOST_SERIALIZATION_NVP( x6 );
    t = boost::make_tuple( x1,x2,x3,x4,x5,x6 );
}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename Archive>
void save( Archive & ar,
           boost::tuple<T1,T2,T3,T4,T5,T6> const& t,
           const unsigned int file_version )
{
    T1 x1;
    T2 x2;
    T3 x3;
    T4 x4;
    T5 x5;
    T6 x6;
    boost::tie( x1,x2,x3,x4,x5,x6 ) = t;
    ar << BOOST_SERIALIZATION_NVP( x1 );
    ar << BOOST_SERIALIZATION_NVP( x2 );
    ar << BOOST_SERIALIZATION_NVP( x3 );
    ar << BOOST_SERIALIZATION_NVP( x4 );
    ar << BOOST_SERIALIZATION_NVP( x5 );
    ar << BOOST_SERIALIZATION_NVP( x6 );

}
template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, class Archive>
void serialize( Archive & ar,
                boost::tuple<T1,T2,T3,T4,T5,T6> & t,
                const unsigned int file_version )
{
    split_free( ar, t, file_version );
}



} // serialization
} //boost

#endif /* __Serialization_H */
