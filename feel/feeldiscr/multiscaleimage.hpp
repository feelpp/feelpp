/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
enum { ComputeGradient = 1 << 0  };

template <typename T = float>
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

template<typename T, int _Options = 0>
class MultiScaleImage
{
public :
    using value_type = T;
    static const int Options = _Options;
    using needs_gradient_t = mpl::bool_<Options&ComputeGradient>;
    using do_compute_gradient_t = mpl::bool_<true>;
    using no_compute_gradient_t = mpl::bool_<false>;    
    // true if must compute gradient, false otherwise.
    static const bool needs_gradient = needs_gradient_t::value;

    MultiScaleImage(holo3_image<value_type> const& im, float L)
        :
        dx(doption("msi.pixelsize")),dy(doption("msi.pixelsize")),image(im),level(L)
    {
    }
    /**
     * @return the component \c c of the gradient of the image at point \c real
     * in the coarse grid
     */
    value_type 
    operator()(int c, ublas::vector<double> const& real,ublas::vector<double> const& ref ) const
        {
            double x = real[0];
            double y = real[1];
             
            int i = boost::math::iround(x/dx);
            //int j = image.cols()-1-boost::math::iround(y/dy);
            int j = boost::math::iround(y/dy);
            
            double v=  image(j,i);
            // v has the value if the image, now must compute the basis functions
            // note that it would be differently handled if we use the fft and ifft
            return v;
        }
    /**
     * @return the value of the image at point \c real in the coarse grid
     */
    value_type 
    operator()(ublas::vector<double> const& real,ublas::vector<double> const& ref ) const
        {
            double x = real[0];
            double y = real[1];
             
            int i = boost::math::iround(x/dx);
            //int j = image.cols()-1-boost::math::iround(y/dy);
            int j = boost::math::iround(y/dy);
            
            double v=  image(j,i);
#if 0
            std::cout << "Value " << v << " Coarse real (" << x <<"," << y 
                      << ") Ref : ("<< ref[0] << "," << ref[1]  
                      << ") Fine image coord. i =" << i <<", j =" << j << std::endl;
#endif

            return v;
        }
private :
    double dx;
    double dy;
    holo3_image<value_type> image;
    int level;

};

} // Feel
