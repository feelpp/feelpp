/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 May 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_GAUSSHYPERCUBE_HPP
#define FEELPP_GAUSSHYPERCUBE_HPP 1

#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>

#include <feel/feelpoly/gauss.hpp>
#include <feel/feelpoly/imfactory.hpp>



namespace Feel {

namespace detail{

template<int D, typename T>
class GaussHypercube
{
public :
    typedef T value_type;

    typedef IMBase<T> super;
    
    //static inline const uint16_type this->numberOfPoints() = ( Integration_this->numberOfPoints()+1 )/2+1;
    //static inline const uint32_type Npoints = this->numberOfPoints();

    GaussHypercube() = default;
    GaussHypercube( GaussHypercube const& ) = default;
    GaussHypercube( GaussHypercube && ) = default;
    GaussHypercube& operator=( GaussHypercube const& ) = default;
    GaussHypercube& operator=( GaussHypercube && ) = default;

    //! @return the degree of the quadrature given the order of the
    //! polynomial degree to integrate
    constexpr uint16_type gaussdegree( uint16_type o ) noexcept { return (o+1)/2+1; }
    
    //!
    //! create a Gauss quadrature on hypercube integrating up to order o
    //!
    explicit GaussHypercube( uint16_type o )
                    
        :
        M_order( o ) 
        {}
    ~GaussHypercube() = default;

    std::unique_ptr<IMBase<T>> operator()()
    {
        return operator()( mpl::int_<D>() );
    }
    std::unique_ptr<IMBase<T>> operator()( mpl::int_<1> ) 
        {
            std::unique_ptr<IMBase<T>> p( std::make_unique<IMBase<T>>(D,gaussdegree(M_order), ipow(gaussdegree(M_order), D) ) );
            // build rules in x and y direction
            std::vector<T> wx( p->numberOfPoints() );
            std::vector<T> px( p->numberOfPoints() );
            Feel::details::dyna::gaussjacobi( p->numberOfPoints(), wx, px, 0.0, 0.0 );
#if 0
            VLOG(1) << "[gauss<SP<2,1>] jacobi p = " << px << "\n";
            VLOG(1) << "[gauss<SP<2,1>] jacobi w = " << wx << "\n";
#endif
            p->q.resize( (D+1)*p->numberOfPoints() );
            for ( int i = 0; i < p->degree(); ++i )
            {
                // computes the weight of the k-th node
                p->q[(D+1)*i] = px[i];
                p->q[(D+1)*i+1] = wx[i];
            }


#if 0
            VLOG(1) << "[gauss<SP<2,1>] p = " << p->M_points << "\n";
            VLOG(1) << "[gauss<SP<2,1>] w = " << p->M_w << "\n";
#endif
            p->setDefined();
            return p;
        }

    std::unique_ptr<IMBase<T>> operator()( mpl::int_<2> ) 
    {
        std::unique_ptr<IMBase<T>> p( std::make_unique<IMBase<T>>(D,gaussdegree(M_order), ipow(gaussdegree(M_order), D) ) );
        // build rules in x and y direction
        std::vector<T> wx( p->degree() );
        std::vector<T> px( p->degree() );
        Feel::details::dyna::gaussjacobi( p->degree(), wx, px, 0.0, 0.0 );
        
        p->q.resize( (D+1)*p->numberOfPoints() );
        for ( int i = 0,  k = 0; i < p->degree(); ++i )
        {
            for ( int j = 0; j < p->degree(); ++j, ++k )
            {
                // computes the weight of the k-th node
                p->q[(D+1)*k] = px[i];
                p->q[(D+1)*k+1] = px[j];
                p->q[(D+1)*k+2] = wx[ i ] * wx[ j ];
            }
        }

        p->setDefined();
        return p;
    }
    std::unique_ptr<IMBase<T>> operator()( mpl::int_<3> ) 
        {
            std::unique_ptr<IMBase<T>> p( std::make_unique<IMBase<T>>(D,gaussdegree(M_order), ipow(gaussdegree(M_order), D) ) );
            // build rules in x and y direction
            std::vector<T> wx( p->degree() );
            std::vector<T> px( p->degree() );
            Feel::details::dyna::gaussjacobi( p->degree(), wx, px, 0.0, 0.0 );
            
            p->q.resize( (D+1)*p->numberOfPoints() );
            for ( int i = 0,  k = 0; i < p->degree(); ++i )
            {
                for ( int j = 0; j < p->degree(); ++j )
                {
                    for ( int l = 0; l < p->degree() ; ++l, ++k )
                    {
                        // computes the weight of the k-th node
                        p->q[(D+1)*k] = px[i];
                        p->q[(D+1)*k+1] = px[j];
                        p->q[(D+1)*k+2] = px[l];
                        p->q[(D+1)*k+3] = wx[ i ] * wx[ j ] * wx[ l ];
                    }
                }
            }
            return p;
            p->setDefined();
        }
    

#if 0
            GaussHypercube( uint16_type o, mpl::int_<4> )
                :
                super( 4, o, (o+1)/2+1 )
            {
                // build rules in x and y direction
                std::vector<T> wx( p->numberOfPoints() );
                std::vector<T> px( p->numberOfPoints() );
                Feel::details::dyna::gaussjacobi( p->numberOfPoints(), wx, px, 0.0, 0.0 );

                p->q.resize( (D+1)*p->numberOfPoints() );
                for ( int i = 0,  k = 0; i < p->numberOfPoints(); ++i )
                {
                    for ( int j = 0; j < p->numberOfPoints(); ++j )
                    {
                        for ( int l = 0; l < p->numberOfPoints() ; ++l )
                        {
                            for ( int r = 0; r < p->numberOfPoints() ; ++r, ++k )
                            {
                                // computes the weight of the k-th node
                                p->q[(D+1)*k] = px[i];
                                p->q[(D+1)*k+1] = px[j];
                                p->q[(D+1)*k+2] = px[l];
                                p->q[(D+1)*k+3] = px[r];
                                p->q[(D+1)*k+4] = wx[ i ] * wx[ j ] * wx[ l ]* wx[ r ];
                            }
                        }
                    }
                }
            }

            GaussHypercube( uint16_type o, mpl::int_<5> )
                :
                super( 5, o, (o+1)/2+1)
            {
                // build rules in x and y direction
                std::vector<T> wx( p->numberOfPoints() );
                std::vector<T> px( p->numberOfPoints() );
                Feel::details::dyna::gaussjacobi( p->numberOfPoints(), wx, px, 0.0, 0.0 );

                p->q.resize( (D+1)*p->numberOfPoints() );
                for ( int i = 0,  k = 0; i < p->numberOfPoints(); ++i )
                {
                    for ( int j = 0; j < p->numberOfPoints(); ++j )
                    {
                        for ( int l = 0; l < p->numberOfPoints() ; ++l )
                        {
                            for ( int r = 0; r < p->numberOfPoints() ; ++r )
                            {
                                for ( int s = 0; s < p->numberOfPoints() ; ++s, ++k )
                                {
                                    p->q[(D+1)*k] = px[i];
                                    p->q[(D+1)*k+1] = px[j];
                                    p->q[(D+1)*k+2] = px[l];
                                    p->q[(D+1)*k+3] = px[r];
                                    p->q[(D+1)*k+4] = px[s];
                                    p->q[(D+1)*k+5] = wx[ i ] * wx[ j ] * wx[ l ]* wx[ r ] * wx[ s ];
                                }
                            }
                        }
                    }
                }
            }
#endif
    int M_order = 0;
};

} // DETAILS
/// \endcond

const bool im10gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,0,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(0) );
const bool im11gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,1,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(1) );
const bool im12gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,2,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(2) );
const bool im13gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,3,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(3) );
const bool im14gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,4,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(4) );
const bool im15gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,5,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(5) );
const bool im16gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,6,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(6) );
const bool im17gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,7,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(7) );
const bool im18gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,8,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(8) );
const bool im19gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,9,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(9) );
const bool im110gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,10,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(10) );
const bool im111gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,11,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(11) );
const bool im112gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,12,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(12) );
const bool im113gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,13,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(13) );
const bool im114gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,14,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(14) );
const bool im115gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,15,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(15) );
const bool im116gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,16,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(16) );
const bool im117gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,17,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(17) );
const bool im118gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,18,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(18) );
const bool im119gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,19,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(19) );
const bool im120gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,20,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(20) );
const bool im121gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,21,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(21) );
const bool im122gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,22,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(22) );
const bool im123gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,23,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(23) );
const bool im124gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,24,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(24) );
const bool im125gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,25,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(25) );
const bool im126gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,26,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(26) );
const bool im127gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,27,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(27) );
const bool im128gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,28,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(28) );
const bool im129gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,29,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(29) );
const bool im130gausshypercube = IMFactory<double>::instance().registerProduct( "im(1,30,hypercube)", Feel::detail::GaussHypercube< 1 ,double>(30) );
const bool imMaxOrder1gausshypercube = IMMaxOrderFactory<double>::instance().insert( std::pair< std::string, IMMaxOrder<double> >( "im(1,hypercube)",30 ) ).second;

const bool im20gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,0,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(0) );
const bool im21gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,1,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(1) );
const bool im22gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,2,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(2) );
const bool im23gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,3,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(3) );
const bool im24gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,4,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(4) );
const bool im25gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,5,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(5) );
const bool im26gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,6,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(6) );
const bool im27gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,7,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(7) );
const bool im28gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,8,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(8) );
const bool im29gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,9,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(9) );
const bool im210gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,10,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(10) );
const bool im211gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,11,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(11) );
const bool im212gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,12,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(12) );
const bool im213gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,13,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(13) );
const bool im214gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,14,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(14) );
const bool im215gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,15,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(15) );
const bool im216gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,16,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(16) );
const bool im217gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,17,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(17) );
const bool im218gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,18,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(18) );
const bool im219gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,19,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(19) );
const bool im220gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,20,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(20) );
const bool im221gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,21,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(21) );
const bool im222gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,22,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(22) );
const bool im223gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,23,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(23) );
const bool im224gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,24,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(24) );
const bool im225gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,25,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(25) );
const bool im226gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,26,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(26) );
const bool im227gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,27,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(27) );
const bool im228gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,28,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(28) );
const bool im229gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,29,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(29) );
const bool im230gausshypercube = IMFactory<double>::instance().registerProduct( "im(2,30,hypercube)", Feel::detail::GaussHypercube< 2 ,double>(30) );
const bool imMaxOrder2gausshypercube = IMMaxOrderFactory<double>::instance().insert( std::pair< std::string, IMMaxOrder<double> >( "im(2,hypercube)",30 ) ).second;

const bool im30gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,0,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(0) );
const bool im31gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,1,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(1) );
const bool im32gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,2,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(2) );
const bool im33gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,3,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(3) );
const bool im34gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,4,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(4) );
const bool im35gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,5,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(5) );
const bool im36gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,6,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(6) );
const bool im37gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,7,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(7) );
const bool im38gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,8,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(8) );
const bool im39gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,9,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(9) );
const bool im310gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,10,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(10) );
const bool im311gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,11,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(11) );
const bool im312gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,12,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(12) );
const bool im313gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,13,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(13) );
const bool im314gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,14,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(14) );
const bool im315gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,15,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(15) );
const bool im316gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,16,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(16) );
const bool im317gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,17,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(17) );
const bool im318gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,18,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(18) );
const bool im319gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,19,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(19) );
const bool im320gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,20,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(20) );
const bool im321gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,21,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(21) );
const bool im322gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,22,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(22) );
const bool im323gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,23,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(23) );
const bool im324gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,24,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(24) );
const bool im325gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,25,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(25) );
const bool im326gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,26,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(26) );
const bool im327gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,27,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(27) );
const bool im328gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,28,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(28) );
const bool im329gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,29,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(29) );
const bool im330gausshypercube = IMFactory<double>::instance().registerProduct( "im(3,30,hypercube)", Feel::detail::GaussHypercube< 3 ,double>(30) );
const bool imMaxOrder3gausshypercube = IMMaxOrderFactory<double>::instance().insert( std::pair< std::string, IMMaxOrder<double> >( "im(3,hypercube)",30 ) ).second;

const bool im10gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,0,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(0) );
const bool im11gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,1,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(1) );
const bool im12gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,2,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(2) );
const bool im13gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,3,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(3) );
const bool im14gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,4,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(4) );
const bool im15gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,5,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(5) );
const bool im16gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,6,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(6) );
const bool im17gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,7,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(7) );
const bool im18gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,8,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(8) );
const bool im19gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,9,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(9) );
const bool im110gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,10,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(10) );
const bool im111gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,11,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(11) );
const bool im112gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,12,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(12) );
const bool im113gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,13,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(13) );
const bool im114gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,14,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(14) );
const bool im115gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,15,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(15) );
const bool im116gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,16,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(16) );
const bool im117gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,17,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(17) );
const bool im118gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,18,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(18) );
const bool im119gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,19,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(19) );
const bool im120gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,20,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(20) );
const bool im121gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,21,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(21) );
const bool im122gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,22,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(22) );
const bool im123gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,23,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(23) );
const bool im124gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,24,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(24) );
const bool im125gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,25,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(25) );
const bool im126gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,26,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(26) );
const bool im127gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,27,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(27) );
const bool im128gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,28,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(28) );
const bool im129gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,29,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(29) );
const bool im130gausshypercubef = IMFactory<float>::instance().registerProduct( "im(1,30,hypercube)", Feel::detail::GaussHypercube< 1 ,float>(30) );
const bool imMaxOrder1gausshypercubef = IMMaxOrderFactory<float>::instance().insert( std::pair< std::string, IMMaxOrder<float> >( "im(1,hypercube)",30 ) ).second;

const bool im20gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,0,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(0) );
const bool im21gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,1,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(1) );
const bool im22gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,2,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(2) );
const bool im23gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,3,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(3) );
const bool im24gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,4,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(4) );
const bool im25gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,5,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(5) );
const bool im26gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,6,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(6) );
const bool im27gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,7,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(7) );
const bool im28gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,8,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(8) );
const bool im29gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,9,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(9) );
const bool im210gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,10,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(10) );
const bool im211gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,11,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(11) );
const bool im212gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,12,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(12) );
const bool im213gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,13,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(13) );
const bool im214gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,14,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(14) );
const bool im215gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,15,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(15) );
const bool im216gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,16,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(16) );
const bool im217gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,17,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(17) );
const bool im218gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,18,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(18) );
const bool im219gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,19,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(19) );
const bool im220gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,20,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(20) );
const bool im221gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,21,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(21) );
const bool im222gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,22,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(22) );
const bool im223gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,23,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(23) );
const bool im224gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,24,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(24) );
const bool im225gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,25,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(25) );
const bool im226gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,26,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(26) );
const bool im227gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,27,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(27) );
const bool im228gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,28,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(28) );
const bool im229gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,29,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(29) );
const bool im230gausshypercubef = IMFactory<float>::instance().registerProduct( "im(2,30,hypercube)", Feel::detail::GaussHypercube< 2 ,float>(30) );
const bool imMaxOrder2gausshypercubef = IMMaxOrderFactory<float>::instance().insert( std::pair< std::string, IMMaxOrder<float> >( "im(2,hypercube)",30 ) ).second;

const bool im30gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,0,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(0) );
const bool im31gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,1,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(1) );
const bool im32gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,2,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(2) );
const bool im33gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,3,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(3) );
const bool im34gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,4,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(4) );
const bool im35gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,5,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(5) );
const bool im36gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,6,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(6) );
const bool im37gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,7,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(7) );
const bool im38gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,8,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(8) );
const bool im39gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,9,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(9) );
const bool im310gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,10,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(10) );
const bool im311gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,11,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(11) );
const bool im312gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,12,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(12) );
const bool im313gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,13,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(13) );
const bool im314gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,14,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(14) );
const bool im315gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,15,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(15) );
const bool im316gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,16,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(16) );
const bool im317gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,17,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(17) );
const bool im318gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,18,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(18) );
const bool im319gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,19,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(19) );
const bool im320gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,20,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(20) );
const bool im321gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,21,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(21) );
const bool im322gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,22,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(22) );
const bool im323gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,23,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(23) );
const bool im324gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,24,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(24) );
const bool im325gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,25,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(25) );
const bool im326gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,26,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(26) );
const bool im327gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,27,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(27) );
const bool im328gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,28,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(28) );
const bool im329gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,29,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(29) );
const bool im330gausshypercubef = IMFactory<float>::instance().registerProduct( "im(3,30,hypercube)", Feel::detail::GaussHypercube< 3 ,float>(30) );
const bool imMaxOrder3gausshypercubef = IMMaxOrderFactory<float>::instance().insert( std::pair< std::string, IMMaxOrder<float> >( "im(3,hypercube)",30 ) ).second;

}
#endif
