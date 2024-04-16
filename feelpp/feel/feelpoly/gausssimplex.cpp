/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 10 May 2015
 
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
#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feelpoly/imfactory.hpp>

namespace Feel {

namespace detail{
template<int D, typename T>
class GaussSimplex
{
public :
    typedef T value_type;
    using super = IMBase<T>;
    //static inline const uint16_type Degree = ( Integration_Degree+1 )/2+1;
    //static inline const uint32_type Npoints = Degree;

    typedef Gauss<Simplex<D-1,1>,invalid_uint16_type_value, T> face_quad_type;

    GaussSimplex() = default;

    explicit GaussSimplex( uint16_type o )
        :
        M_order( o )
        {}

    std::unique_ptr<IMBase<T>> operator()() 
        {
            std::unique_ptr<IMBase<T>> p( std::make_unique<IMBase<T>>(D,M_order,(M_order+1)/2+1) );

            std::vector<T> px( p->numberOfPoints() );
            std::vector<T> w( p->numberOfPoints() );
            
            Feel::details::dyna::gaussjacobi<T, std::vector<T>, std::vector<T> >( p->numberOfPoints(), w, px );
            
            p->q.resize( (D+1)*p->numberOfPoints() );
            for( int i = 0; i < p->numberOfPoints(); ++i )
            {
                p->q[(D+1)*i] = px[i];
                p->q[(D+1)*i+1] = w[i];
            }
            p->setDefined();
            return p;
        }
    ~GaussSimplex() = default;
    int M_order = 0;
    
    
};
}
const bool im10gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,0,simplex)", Feel::detail::GaussSimplex< 1 ,double>(0) );
const bool im11gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,1,simplex)", Feel::detail::GaussSimplex< 1 ,double>(1) );
const bool im12gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,2,simplex)", Feel::detail::GaussSimplex< 1 ,double>(2) );
const bool im13gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,3,simplex)", Feel::detail::GaussSimplex< 1 ,double>(3) );
const bool im14gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,4,simplex)", Feel::detail::GaussSimplex< 1 ,double>(4) );
const bool im15gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,5,simplex)", Feel::detail::GaussSimplex< 1 ,double>(5) );
const bool im16gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,6,simplex)", Feel::detail::GaussSimplex< 1 ,double>(6) );
const bool im17gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,7,simplex)", Feel::detail::GaussSimplex< 1 ,double>(7) );
const bool im18gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,8,simplex)", Feel::detail::GaussSimplex< 1 ,double>(8) );
const bool im19gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,9,simplex)", Feel::detail::GaussSimplex< 1 ,double>(9) );
const bool im110gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,10,simplex)", Feel::detail::GaussSimplex< 1 ,double>(10) );
const bool im111gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,11,simplex)", Feel::detail::GaussSimplex< 1 ,double>(11) );
const bool im112gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,12,simplex)", Feel::detail::GaussSimplex< 1 ,double>(12) );
const bool im113gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,13,simplex)", Feel::detail::GaussSimplex< 1 ,double>(13) );
const bool im114gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,14,simplex)", Feel::detail::GaussSimplex< 1 ,double>(14) );
const bool im115gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,15,simplex)", Feel::detail::GaussSimplex< 1 ,double>(15) );
const bool im116gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,16,simplex)", Feel::detail::GaussSimplex< 1 ,double>(16) );
const bool im117gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,17,simplex)", Feel::detail::GaussSimplex< 1 ,double>(17) );
const bool im118gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,18,simplex)", Feel::detail::GaussSimplex< 1 ,double>(18) );
const bool im119gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,19,simplex)", Feel::detail::GaussSimplex< 1 ,double>(19) );
const bool im120gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,20,simplex)", Feel::detail::GaussSimplex< 1 ,double>(20) );
#if 1
const bool im121gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,21,simplex)", Feel::detail::GaussSimplex< 1 ,double>(21) );
const bool im122gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,22,simplex)", Feel::detail::GaussSimplex< 1 ,double>(22) );
const bool im123gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,23,simplex)", Feel::detail::GaussSimplex< 1 ,double>(23) );
const bool im124gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,24,simplex)", Feel::detail::GaussSimplex< 1 ,double>(24) );
const bool im125gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,25,simplex)", Feel::detail::GaussSimplex< 1 ,double>(25) );
const bool im126gausssimplex = IMFactory<double>::instance().registerProduct( "im(1,26,simplex)", Feel::detail::GaussSimplex< 1 ,double>(26) );
#endif
const bool imMaxOrder1gausssimplex = IMMaxOrderFactory<double>::instance().insert( std::pair< std::string, IMMaxOrder<double> >( "im(1,simplex)",26 ) ).second;

const bool im10gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,0,simplex)", Feel::detail::GaussSimplex< 1 ,float>(0) );
const bool im11gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,1,simplex)", Feel::detail::GaussSimplex< 1 ,float>(1) );
const bool im12gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,2,simplex)", Feel::detail::GaussSimplex< 1 ,float>(2) );
const bool im13gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,3,simplex)", Feel::detail::GaussSimplex< 1 ,float>(3) );
const bool im14gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,4,simplex)", Feel::detail::GaussSimplex< 1 ,float>(4) );
const bool im15gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,5,simplex)", Feel::detail::GaussSimplex< 1 ,float>(5) );
const bool im16gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,6,simplex)", Feel::detail::GaussSimplex< 1 ,float>(6) );
const bool im17gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,7,simplex)", Feel::detail::GaussSimplex< 1 ,float>(7) );
const bool im18gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,8,simplex)", Feel::detail::GaussSimplex< 1 ,float>(8) );
const bool im19gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,9,simplex)", Feel::detail::GaussSimplex< 1 ,float>(9) );
const bool im110gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,10,simplex)", Feel::detail::GaussSimplex< 1 ,float>(10) );
const bool im111gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,11,simplex)", Feel::detail::GaussSimplex< 1 ,float>(11) );
const bool im112gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,12,simplex)", Feel::detail::GaussSimplex< 1 ,float>(12) );
const bool im113gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,13,simplex)", Feel::detail::GaussSimplex< 1 ,float>(13) );
const bool im114gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,14,simplex)", Feel::detail::GaussSimplex< 1 ,float>(14) );
const bool im115gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,15,simplex)", Feel::detail::GaussSimplex< 1 ,float>(15) );
const bool im116gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,16,simplex)", Feel::detail::GaussSimplex< 1 ,float>(16) );
const bool im117gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,17,simplex)", Feel::detail::GaussSimplex< 1 ,float>(17) );
const bool im118gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,18,simplex)", Feel::detail::GaussSimplex< 1 ,float>(18) );
const bool im119gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,19,simplex)", Feel::detail::GaussSimplex< 1 ,float>(19) );
const bool im120gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,20,simplex)", Feel::detail::GaussSimplex< 1 ,float>(20) );
#if 1
const bool im121gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,21,simplex)", Feel::detail::GaussSimplex< 1 ,float>(21) );
const bool im122gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,22,simplex)", Feel::detail::GaussSimplex< 1 ,float>(22) );
const bool im123gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,23,simplex)", Feel::detail::GaussSimplex< 1 ,float>(23) );
const bool im124gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,24,simplex)", Feel::detail::GaussSimplex< 1 ,float>(24) );
const bool im125gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,25,simplex)", Feel::detail::GaussSimplex< 1 ,float>(25) );
const bool im126gausssimplexf = IMFactory<float>::instance().registerProduct( "im(1,26,simplex)", Feel::detail::GaussSimplex< 1 ,float>(26) );
#endif
const bool imMaxOrder1gausssimplexf = IMMaxOrderFactory<float>::instance().insert( std::pair< std::string, IMMaxOrder<float> >( "im(1,simplex)",26 ) ).second;

}
