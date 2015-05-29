/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s) : Thomas Lantz
       Date: 2015-04-27

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
   \file multiscaleimage.hpp
   \author Thomas Lantz 
   \date 2015-04-27
 */

//#ifndef _MULTISCALEIMAGE_HPP_
//#define _MULTISCALEIMAGE_HPP_
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/special_functions/round.hpp>
using namespace boost::numeric;

namespace Feel
{
template <typename T = float>
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

template<typename T>
class MultiScaleImage
{
public :
    using value_type = T;
    
    MultiScaleImage(holo3_image<value_type> const& im, float L)
        :
        image(im),level(L)
    {
    }

    value_type 
    operator()(ublas::vector<double> const& c) const
        {
            double x = c[0];
            double y = c[1];
    
            std::cout <<" x =" << x <<", y =" << y << std::endl;
 
            double i = boost::math::iround(x/*pow(2,level)*//dx);
            double j = boost::math::iround(y/*pow(2,level)*//dy);
            
            std::cout <<" i =" << i <<", j =" << j << std::endl;


            return image(j,i);
        }
/*
    value_type operator()(ublas::vector<double> const& c, int L)
        {     
            double x = c[0];
            double y = c[1];
            
            int i = x/dx;
            int j = y/dy;
            
            return image(L*j,L*i);
        }
  */   
private :
    double dx =8.9e-3;
    double dy =8.9e-3;
    holo3_image<value_type> image;
    int level;
};

} // Feel
