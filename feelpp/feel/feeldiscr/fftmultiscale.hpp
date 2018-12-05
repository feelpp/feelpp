/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Thomas Lantz <lantz.thomas0@gmail.com>
   Date: 2015-07-31

Copyright (C) 2014-2016 Feel++ Consortium

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
#ifndef FEELPP_FFTMULTISCALE_HPP
#define FEELPP_FFTMULTISCALE_HPP 1

#include <Eigen/Core>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wredeclared-class-member"
#endif
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <math.h>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeltiming/tic.hpp>

#if defined( FEELPP_HAS_FFTW )
#include <fftw3.h>


// define a class to estimate gradient by FFT 
using namespace boost::numeric;

namespace Feel
{
template <typename T = float>
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

class FFTFeel
{
public :    

// Modify real matrix structure to complex matrix array 
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

// Rebuild the image structure from an complex array 
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

// exchange left and right part of an complex array
void fftshift(fftw_complex* data, int s)
{
    int k = 0;
    int count = s;//sizeof(data)/sizeof(data[0]);
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


// main function : Esimate X and Y gradient field by the FFT method
void gradImage (holo3_image<float> const& im )
{
    // define useful parameters 
    int meshp=std::pow(2,ioption("msi.level")); 
    int nLi=im.rows()/meshp; 
    int nCo=im.cols()/meshp;

    // define wave number 
    fftw_complex* nbOndes;
    nbOndes=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(meshp+1));
    if (meshp==1)
    {
        nbOndes[0][0]=-1.0;
        nbOndes[0][1]=0.0;
        nbOndes[1][0]=1.0;
        nbOndes[1][1]=0.0;
    }
    else
    {   
        for (int i=0;i<=meshp/2;i++)
        {
            nbOndes[i][0]=i;
        }

        for ( int i=0;i<meshp/2;i++)
        {
            nbOndes[meshp/2+1+i][0]=-nbOndes[meshp/2-i][0];
        }
    }   
    for (int i=0;i<meshp+1;i++)
    {
        nbOndes[i][0]*=2*(4.*math::atan(1.))/(doption("msi.pixelsize")*doption("msi.pixelsize"));
        
        std::ofstream fichierOndes ("Ondes.txt",std::ios::out|std::ios::trunc);
        if (fichierOndes)
        {
            for (int i=0;i<std::pow(2,ioption("msi.level"))+1;i++)
            {
                fichierOndes << nbOndes[i][0] << " " ;
            }
            fichierOndes.close();
        }


        //std::cout << "nb2 :\t" << nbOndes[i][0] << std::endl;
    }
    
    // define arrays which will store gradient
    tabFFTX=new holo3_image<float>* [nLi];
    tabFFTY=new holo3_image<float>* [nLi];

    for (int i=0;i<nLi;i++)
    {
        tabFFTX[i]=new holo3_image<float> [nCo];
        tabFFTY[i]=new holo3_image<float> [nCo]; 
        for (int j=0;j<nCo;j++)
        {
            tabFFTX[i][j]= holo3_image<float>(meshp+1,meshp+1);
            tabFFTY[i][j]= holo3_image<float>(meshp+1,meshp+1);
            fftw_complex* image = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(meshp+1)*(meshp+1));
            
            for (int k=0;k<meshp+1;k++)
            {
                for (int l=0;l<meshp+1;l++)
                {
                    image[(meshp+1)*k+l][0]=im((i)*(meshp)+k,(j)*(meshp)+l);
                    image[(meshp+1)*k+l][1]=0;
                   //std::cout << "image \t:" <<tabFFTX[i][j](k,l) << std::endl;
                }
            }

            // define image fft
            std::ofstream fichierImage ("vImage.txt",std::ios::out|std::ios::trunc);
            if (fichierImage)
            {
                for (int i=0;i<std::pow(2,ioption("msi.level"))+1;i++)
                {
                    for (int j=0;j<std::pow(2,ioption("msi.level"))+1;j++)
                    {
                        fichierImage << image[(meshp+1)*i+j][0] << "+i" <<  image[(meshp+1)*i+j][1] << " " ;
                    }
                    fichierImage << "\n" ;
                }
                fichierImage.close();
            }

            fftw_complex* imageFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(meshp+1)*(meshp+1));
            fftw_complex* imageFFTtmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(meshp+1)*(meshp+1));

            fftw_plan forward= fftw_plan_dft_2d(meshp+1,meshp+1,image,imageFFTtmp,FFTW_FORWARD,FFTW_ESTIMATE);
            fftw_execute(forward);
            fftw_destroy_plan(forward);

            std::ofstream fichierFFT ("vPreFFT.txt",std::ios::out|std::ios::trunc);
            if (fichierFFT)
            {
                for (int i=0;i<std::pow(2,ioption("msi.level"))+1;i++)
                {
                    for (int j=0;j<std::pow(2,ioption("msi.level"))+1;j++)
                    {
                        fichierFFT << imageFFTtmp[(meshp+1)*i+j][0] << "+i" <<  imageFFTtmp[(meshp+1)*i+j][1] << " " ;
                    }
                    fichierFFT << "\n" ;
                }
                fichierFFT.close();
            }
             
            // multiply by i 
            for (int k=0;k<(meshp+1)*(meshp+1);k++)
            {
                imageFFT[k][0]=-1*imageFFTtmp[k][1];
                imageFFT[k][1]=imageFFTtmp[k][0];
                //std::cout << "\t FFT :\t" << imageFFT[k][1] << std::endl;
            }

            std::ofstream fichierFFT2 ("vFFT.txt",std::ios::out|std::ios::trunc);
            if (fichierFFT)
            {
                for (int i=0;i<std::pow(2,ioption("msi.level"))+1;i++)
                {
                    for (int j=0;j<std::pow(2,ioption("msi.level"))+1;j++)
                    {
                        fichierFFT2 << imageFFT[(meshp+1)*i+j][0] << "+i" <<  imageFFT[(meshp+1)*i+j][1] << " " ;
                    }
                    fichierFFT2 << "\n" ;
                }
                fichierFFT2.close();
            }   


            // add waves number contributions 
            fftw_complex* FFTx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(meshp+1)*(meshp+1));
            fftw_complex* FFTy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(meshp+1)*(meshp+1));

            for (int k=0;k<meshp+1;k++)
            {
                for (int l=0;l<meshp+1;l++)
                {
                    FFTy[(meshp+1)*k+l][0]=imageFFT[(meshp+1)*k+l][0]*nbOndes[k][0];
                    FFTy[(meshp+1)*k+l][1]=imageFFT[(meshp+1)*k+l][1]*nbOndes[k][0];

                    //std::cout << "\t FFT2 :\t" << FFTy[(meshp+1)*k+l][1] << std::endl;

                    FFTx[(meshp+1)*k+l][0]=imageFFT[(meshp+1)*k+l][0]*nbOndes[l][0];
                    FFTx[(meshp+1)*k+l][1]=imageFFT[(meshp+1)*k+l][1]*nbOndes[l][0];
                }
            }

            std::ofstream fichierFFT3 ("vOFFT.txt",std::ios::out|std::ios::trunc);
            if (fichierFFT3)
            {
                for (int i=0;i<std::pow(2,ioption("msi.level"))+1;i++)
                {
                    for (int j=0;j<std::pow(2,ioption("msi.level"))+1;j++)
                    {
                        fichierFFT3 << FFTx[(meshp+1)*i+j][0] << "+i" <<  FFTx[(meshp+1)*i+j][1] << " " ;
                    }
                    fichierFFT3 << "\n" ;
                }
                fichierFFT3.close();
            }

            //estimation ifft entry  
            fftw_plan backwardX = fftw_plan_dft_2d(meshp+1,meshp+1,FFTx,FFTx,FFTW_BACKWARD,FFTW_ESTIMATE);
            fftw_plan backwardY = fftw_plan_dft_2d(meshp+1,meshp+1,FFTy,FFTy,FFTW_BACKWARD,FFTW_ESTIMATE);
            fftw_execute(backwardX);
            fftw_execute(backwardY);

            fftw_destroy_plan(backwardX);
            fftw_destroy_plan(backwardY);
            
            std::ofstream fichierFFT4 ("vBFFT.txt",std::ios::out|std::ios::trunc);
            if (fichierFFT4)
            {
                for (int i=0;i<std::pow(2,ioption("msi.level"))+1;i++)
                {
                    for (int j=0;j<std::pow(2,ioption("msi.level"))+1;j++)
                    {
                        fichierFFT4 << FFTx[(meshp+1)*i+j][0] << "+i" <<  FFTx[(meshp+1)*i+j][1] << " " ;
                    }
                    fichierFFT4 << "\n" ;
                }
                fichierFFT4.close();
            }

            // Retrieve the real part and store it
            for (int k=0;k<meshp+1;k++)
            {
                for (int l=0;l<meshp+1;l++)
                {
                    tabFFTX[i][j](k,l)=FFTx[(meshp+1)*k+l][0]/((meshp+1));
                    
                    tabFFTY[i][j](k,l)=FFTy[(meshp+1)*k+l][0]/((meshp+1));
                    //std::cout << "\t final :\t"<<tabFFTX[i][j] (k,l) << std::endl;
                }
            }

            // free space
            fftw_free(image);
            fftw_free(imageFFTtmp);
            fftw_free(imageFFT);
            fftw_free(FFTx);
            fftw_free(FFTy);
        }
    }
}

// return the X gradient field stored
holo3_image<float>** getFFTX()
{
    return tabFFTX;
}

// return the Y gradient field stored
holo3_image<float>** getFFTY()
{
    return tabFFTY;
}


// Test an alternative method with extand element to obtain periodic "images"
/*
void gradImage2 (holo3_image<float> const& im, ublas::vector<double> const& begin )
{
    // define useful parameters 
    int meshp=std::pow(2,ioption("msi.level"))+1; 
    //int nLi=im.rows(); 
    //int nCo=im.cols();
    
    holo3_image<float> work = holo3_image<float> (meshp,meshp);
    for (int i=0;i<meshp;i++)
    {
        for (int j=0;j<meshp;j++)
        {
            work(i,j)=im(begin[0]+i,begin[1]+j);
            std::cout << work(i,j) << std::endl;
        }
    }

    // define wave number 
    fftw_complex* nbOndes;
    nbOndes=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*meshp);
    for (int i=0;i<meshp;i++)
    {
        nbOndes[i][0]=-(meshp+1)/2.+i;
        nbOndes[i][1]=0.0;
        //std::cout << nbOndes[i][0] << std::endl;
    }
    fftshift(nbOndes);
    for (int i=0;i<meshp;i++)
    {
        nbOndes[i][0]*=2*(4.*math::atan(1.));
        
    }
    // define image fft
    fftw_complex* image=hbfToFFT(work);
    fftw_complex* imageFFTtmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*meshp*meshp);
    fftw_complex* imageFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*meshp*meshp);

    fftw_plan forward= fftw_plan_dft_2d(meshp,meshp,image,imageFFTtmp, FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(forward);
    fftw_destroy_plan(forward);

    // multi by i 
    for (int i=0;i<meshp*meshp;i++)
    {
        imageFFT[i][0]=-1*imageFFTtmp[i][1];
        imageFFT[i][1]=imageFFTtmp[i][0];
    }

    // estimation ifft entry
    fftw_complex* FFTx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*meshp*meshp);
    fftw_complex* FFTy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*meshp*meshp);

    for (int i=0;i<meshp;i++)
    {
        for (int j=0;j<meshp;j++)
        {
            FFTy[meshp*i+j][0]=imageFFT[meshp*i+j][0]*nbOndes[i][0];
            FFTy[meshp*i+j][1]=imageFFT[meshp*i+j][1]*nbOndes[i][0];

            FFTx[meshp*i+j][0]=imageFFT[meshp*i+j][0]*nbOndes[j][0];
            FFTx[meshp*i+j][1]=imageFFT[meshp*i+j][1]*nbOndes[j][0];
        }
    }

    fftw_plan backwardX = fftw_plan_dft_2d(meshp,meshp,FFTx,FFTx,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_plan backwardY = fftw_plan_dft_2d(meshp,meshp,FFTy,FFTy,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(backwardX);
    fftw_execute(backwardY);

    fftw_destroy_plan(backwardX);
    fftw_destroy_plan(backwardY);

    imageX =  holo3_image<float> (meshp,meshp);
    imageY =  holo3_image<float> (meshp,meshp);
    

    for (int i=0;i<meshp;i++)
    {
        for (int j=0;j<meshp;j++)
        {
            imageX(i,j)=FFTx[meshp*i+j][0]/(meshp*meshp);
            imageY(i,j)=FFTy[meshp*i+j][0]/(meshp*meshp);
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
*/

private :

// store the X gradient field
holo3_image<float>** tabFFTX;

// store the Y gradient field
holo3_image<float>** tabFFTY;


};



}// Feel

#endif // FEELPP_HAS_FFTW


#endif // FEELPP_FFTMULTISCALE_HPP
