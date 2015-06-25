/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2014-06-04

Copyright (C) 2014 Feel++ Consortium

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
#include <Eigen/Core>
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <math.h>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/expr.hpp>
#include <fftw3.h>

using namespace boost::numeric;

namespace Feel
{
template <typename T = float>
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

class FFTFeel
{
public :    

fftw_complex* hbfToFFT ( holo3_image<float> const& im )
{
    int nLi=im.rows();
    int nCo=im.cols();
    fftw_complex* imfft;
    imfft=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nLi*nCo);

    for (int i=0;i<nLi;i++)
    {
        for(int j=0;j<nCo;j++)
        {
            imfft[nCo*i+j][0]=im(i,j);
            imfft[nCo*i+j][1]=0.0f;
        }
    }
    return imfft;
}

holo3_image<float> FFTtoHbf (fftw_complex* const& imfft, int nLi, int nCo )
{
    holo3_image<float> im (nLi,nCo);
    for (int i=0;i<nLi;i++)
    {
        for(int j=0;j<nCo;j++)
        {
            im(i,j)=imfft[nCo*i+j][0];
        }
    }
    return im;
}

void fftshift(fftw_complex* const& data)
{
    int k = 0;
    int count = sizeof(data)/sizeof(data[0]);
    int c = (int) math::floor((float)count/2);
    fftw_complex tmp;
    if (count % 2 == 0)
    {
        for (k = 0; k < c; k++)
        {
            tmp[0] = data[k][0];
            tmp[1] = data[k][1];
            data[k][0] = data[k+c][0];
            data[k][1] = data[k+c][1];
            data[k+c][0]= tmp[0];
            data[k+c][1]= tmp[1];
        }
    }
    else
    {
        tmp[0] = data[0][0];
        tmp[1] = data[0][1];
        for (k = 0; k < c; k++)
        {
            data[k][0] = data[c + k + 1][0];
            data[k][1] = data[c + k + 1][1];
            data[c + k + 1][0] = data[k + 1][0];
            data[c + k + 1][1] = data[k + 1][1];
        }
        data[c][0] = tmp[0];
        data[c][1] = tmp[0];
    }
}

//holo3_image<float> gradImage (holo3_image<float> const& im )
void gradImage (holo3_image<float> const& im )
{

    // define useful parameters 
    int meshp=std::pow(2,ioption("msi.level")); 
    int nLi=im.rows(); 
    int nCo=im.cols(); 

    // define wave number 
    fftw_complex* nbOndes;
    nbOndes=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*meshp);
    for (int i=0;i<meshp;i++)
    {
        nbOndes[i][0]=-meshp/2.+i;
        nbOndes[i][1]=0.0;
    }
    fftshift(nbOndes);
    for (int i=0;i<meshp;i++)
        nbOndes[i][0]*=2*(4./math::atan(1.));

    // define image fft
    fftw_complex* image=hbfToFFT(im);
    fftw_complex* imageFFTtmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nLi*nCo);
    fftw_complex* imageFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nLi*nCo);

    fftw_plan forward= fftw_plan_dft_2d(nLi,nCo,image,imageFFTtmp, FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // multi by i 
    for (int i=0;i<nLi*nCo;i++)
    {
        imageFFT[i][0]=-1*imageFFTtmp[i][1];
        imageFFT[i][1]=imageFFTtmp[i][0];
    }

    // estimation ifft entry
    fftw_complex* FFTx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nLi*nCo);
    fftw_complex* FFTy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nLi*nCo);

    for (int i=0;i<nLi;i++)
    {
        for (int j=0;j<nCo;j++)
        {
            FFTx[nCo*i+j][0]=imageFFT[nCo*i+j][0]*nbOndes[i%(meshp+1)][0];
            FFTx[nCo*i+j][1]=imageFFT[nCo*i+j][1]*nbOndes[i%(meshp+1)][0];

            FFTy[nCo*i+j][0]=imageFFT[nCo*i+j][0]*nbOndes[j%(meshp+1)][0];
            FFTy[nCo*i+j][1]=imageFFT[nCo*i+j][1]*nbOndes[j%(meshp+1)][0];
        }
    }

    fftw_plan backwardX = fftw_plan_dft_2d(nLi,nCo,FFTx,FFTx,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_plan backwardY = fftw_plan_dft_2d(nLi,nCo,FFTy,FFTy,FFTW_BACKWARD,FFTW_ESTIMATE);

    fftw_execute(backwardX);
    fftw_execute(backwardY);

    fftw_destroy_plan(backwardX);
    fftw_destroy_plan(backwardY);

    imageX =  holo3_image<float> (nLi,nCo);
    imageY =  holo3_image<float> (nLi,nCo);
    

    for (int i=0;i<nLi;i++)
    {
        for (int j=0;j<nCo;j++)
        {
            imageX(i,j)=FFTx[nCo*i+j][0]/(meshp*meshp);
            imageY(i,j)=FFTy[nCo*i+j][0]/(meshp*meshp);
        }
    }


    fftw_free(nbOndes);
    fftw_free(image);
    fftw_free(imageFFTtmp);
    fftw_free(imageFFT);
    fftw_free(FFTx);
    fftw_free(FFTy);



    //return imY;
}

holo3_image<float> getX()
{
    return imageX;
}

holo3_image<float> getY()
{
    return imageY;
}

private :

holo3_image<float> imageX;
holo3_image<float> imageY;
};


}// Feel
