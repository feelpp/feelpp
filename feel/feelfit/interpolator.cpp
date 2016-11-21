/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Vincent Huber <vincent.huber@cemosis.Fr>
Date: 01-04-2016

Copyright (C) 2007-2008 Universite Joseph Fourier (Grenoble I)

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

#include <feel/feelfit/interpolator.hpp>

std::unique_ptr<Interpolator> Interpolator::New( InterpolationType type, std::vector<pair_type> data)
{
    switch(type)
    {
    case P0:
    {
        return std::make_unique<InterpolatorP0>(data, static_cast<InterpolationType_P0>(Feel::ioption("fit.P0")));
        break;
    }
    case P1:
    {
        return std::make_unique<InterpolatorP1>(data,
                                                static_cast<ExtrapolationType_P1>(Feel::ioption("fit.P1_right")), 
                                                static_cast<ExtrapolationType_P1>(Feel::ioption("fit.P1_left"))
                                                );
        break;
    }
    case Spline:
    {
        return std::make_unique<InterpolatorSpline>(data, 
                                                    static_cast<ExtrapolationType_spline>(Feel::ioption("fit.Spline_right")),
                                                    static_cast<ExtrapolationType_spline>(Feel::ioption("fit.Spline_left"))
                                                    );
        break;
    }
    case Akima:
    {
        return std::make_unique<InterpolatorAkima>(data);
        break;
    }
    }
}

std::unique_ptr<Interpolator> Interpolator::New( InterpolationType type, std::string dataFile)
{
    std::vector<std::pair<double, double>> data;
    std::ifstream infile(Feel::Environment::expand(dataFile));
    std::string line;
    double a, b;
    if (infile.is_open()) 
    {
        while (getline(infile,line)) 
        {
            if(line[0] != '#')
            {
                std::istringstream iss(line);
                iss >> a >> b;
                data.push_back({a,b});
            }
        }
        infile.close();
    }
    else
    {
        Feel::cout << "Error opening file" << Feel::Environment::expand(dataFile) << "\n";
    }
    return Interpolator::New( type, data );
}

Interpolator::value_type InterpolatorP0::diff(double _x) { return 0.; }
Interpolator::value_type InterpolatorP0::operator()(double _x)
{
    pair_type lower(_x,0.);
    double v1 = 0.;
    double v2 = 0.;
    double x1 = 0.;
    double x2 = 0.;
    auto it = std::lower_bound(this->M_data.begin(), this->M_data.end(), lower, [] (pair_type a, pair_type b){ return a.first < b.first ; } );
    if(it == M_data.end())
    {
        x1 = x2 = (it-1)->first;
        v2 = v1 = (it-1)->second;
    }
    else if(it == M_data.begin()) 
    {
        x1 = x2 =  (it)->first;
        v1 = v2 = it->second;
    }
    else
    {
        x1 = (it-1)->first;
        x2 = it->first;
        v1 = (it-1)->second;
        v2 = it->second;
    }
    switch(M_iType)
    {
        case 0: return v2; //left
        case 1: return v1; //right
        case 2: {
                    if(_x < 0.5*(x1+x2)) 
                        return v1;
                    else
                        return v2;
                };
        default: return v2;
    }

    return 0;
}

Interpolator::value_type InterpolatorP1::operator()(double _x)
{
    double y1 = 0.;
    double y2 = 0.;
    double x1 = 0.;
    double x2 = 0.;

    Coefficients(_x, y1, y2, x1, x2);
    return ((y1-y2)*_x - (x2*y1-x1*y2))/(x1-x2);
}
Interpolator::value_type InterpolatorP1::diff(double _x) 
{ 
    double y1 = 0.;
    double y2 = 0.;
    double x1 = 0.;
    double x2 = 0.;
    Coefficients(_x, y1, y2, x1, x2);
    return ((y1-y2))/(x1-x2);
}
        void InterpolatorP1::Coefficients(value_type _x, 
                value_type &y1,
                value_type &y2,
                value_type &x1,
                value_type &x2)
        {
            pair_type lower(_x,0.);
            // (it-1)->first < _x < it->first
            auto it = std::upper_bound(M_data.begin(), M_data.end(), lower, [] (pair_type a, pair_type b){ return a.first < b.first ; } );
            if(it == M_data.end()) 
            {
                switch(M_l_type)
                {
                    case zero:
                        {
                            y1 = y2 = 0.;
                            x1 = (it-1)->first;
                            x2 = 2*x1+1; // to avoid division by zero
                            break;
                        }
                    case constant:
                        {
                            y1 = y2 = (it-1)->second;
                            x1 = (it-1)->first;
                            x2 = 2*x1+1; // to avoid division by zero
                            break;
                        }
                    case extrapol:
                        {
                            y1 = (it-2)->second;
                            y2 = (it-1)->second;
                            x1 = (it-2)->first;
                            x2 = (it-1)->first;
                            break;
                        }
                }
            }
            else if(it == M_data.begin()) 
            {
                switch(M_r_type)
                {
                    case zero:
                        {
                            y1 = y2 = 0.;
                            x1 = it->first;
                            x2 = 2*x1+1; // to avoid division by zero
                            break;
                        }
                    case constant:
                        {
                            y1 = y2 = it->second;
                            x1 = it->first;
                            x2 = 2*x1+1; // to avoid division by zero
                            break;
                        }
                    case extrapol:
                        {
                            y1 = (it)->second;
                            y2 = (it+1)->second;
                            x1 = (it)->first;
                            x2 = (it+1)->first;
                            break;
                        }
                }
            }
            else
            {
                y1 = (it-1)->second;
                y2 = (it-1+1)->second;
                x1 = (it-1)->first;
                x2 = (it-1+1)->first;
            }
        }

/*
 * Interpolate with a*x**3 + b*x**2 + c*x + d
 * TODO: interpolate a(x-xi)**3 + ... instead of ax**3 (only to generate a easiest to inverse linear system)
 */
InterpolatorSpline::InterpolatorSpline( std::vector<pair_type> data, ExtrapolationType_spline r = natural, ExtrapolationType_spline l = natural ) : super(data), M_r_type(r), M_l_type(l)
{
    SpMat A((this->M_data.size()-1)*4,(this->M_data.size()-1)*4);
    std::vector<T> coef;
    Eigen::VectorXd b = Eigen::VectorXd::Zero((this->M_data.size()-1)*4);
    for(int i = 0; i<this->M_data.size()-1; i++)
    {
        double x  = this->M_data.at(i).first;
        double xp = this->M_data.at(i+1).first;
        //si(xi) = yi
        coef.push_back(T(1+i*4,i*4+0, x*x*x));
        coef.push_back(T(1+i*4,i*4+1, x*x  ));
        coef.push_back(T(1+i*4,i*4+2, x    ));
        coef.push_back(T(1+i*4,i*4+3, 1.   ));
        //si(xi+1) = yi+1
        coef.push_back(T(1+i*4+1,i*4+0,xp*xp*xp));
        coef.push_back(T(1+i*4+1,i*4+1,xp*xp   ));
        coef.push_back(T(1+i*4+1,i*4+2,xp      ));
        coef.push_back(T(1+i*4+1,i*4+3,1.      ));
        if(i < M_data.size()-2)
        {
            // si'(xi+1) - si+1'(xi+1) = 0
            coef.push_back(T(1+i*4+2,i*4+0, 3.*xp*xp));
            coef.push_back(T(1+i*4+2,i*4+1, 2.*xp   ));
            coef.push_back(T(1+i*4+2,i*4+2, 1       ));
            coef.push_back(T(1+i*4+2,i*4+0+4, -3.*xp*xp));
            coef.push_back(T(1+i*4+2,i*4+1+4, -2.*xp   ));
            coef.push_back(T(1+i*4+2,i*4+2+4, -1       ));

            // si''(xi+1) - si+1''(xi+1) = 0
            coef.push_back(T(1+i*4+3,i*4+0, 6.*xp));
            coef.push_back(T(1+i*4+3,i*4+1, 2.   ));
            coef.push_back(T(1+i*4+3,i*4+0+4, -6.*xp));
            coef.push_back(T(1+i*4+3,i*4+1+4, -2.   ));
        }
        b(1+i*4)   = M_data.at(i).second;
        b(1+i*4+1) = M_data.at(i+1).second;
    }
    double x0 = M_data.at(0).first;
    switch(M_l_type)
    {
        case natural:
            {
                coef.push_back(T(0,0, 6.*x0));
                coef.push_back(T(0,1, 2.   ));
                b(0) = 0.;
                break;
            }
        case clamped:
            {
                //s0'(x0) = 0
                coef.push_back(T(0,0, 3.*x0*x0));
                coef.push_back(T(0,1, 2.*x0   ));
                coef.push_back(T(0,2, 1       ));
                b(0) = 1.; //diff_l;
                break;

            }
    }
    double xn = M_data.at(M_data.size()-1).first;
    switch(M_r_type)
    {
        case natural:
            {
                coef.push_back(T(4*(M_data.size()-1)-1,(M_data.size()-2)*4+0, 6.*xn));
                coef.push_back(T(4*(M_data.size()-1)-1,(M_data.size()-2)*4+1, 2    ));
                b(4*(M_data.size()-1)-1) = 0.;
                break;
            }
        case clamped:
            {
                //sn'(xn) = 0
                coef.push_back(T(4*(M_data.size()-1)-1,(M_data.size()-2)*4+0, 3.*xn*xn));
                coef.push_back(T(4*(M_data.size()-1)-1,(M_data.size()-2)*4+1, 2.*xn   ));
                coef.push_back(T(4*(M_data.size()-1)-1,(M_data.size()-2)*4+2, 1       ));
                b(4*(M_data.size()-1)-1) = 0.;//diff_r;
                break;
            }
    }
    A.setFromTriplets(coef.begin(), coef.end());
    Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int>> solver(A);
    sol = solver.solve(b);
}
Interpolator::value_type InterpolatorSpline::operator()(double _x)
{
    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    Coefficients(_x, a, b, c, d);
    return a*_x*_x*_x+b*_x*_x+c*_x+d;
}
Interpolator::value_type InterpolatorSpline::diff(double _x)
{
    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    Coefficients(_x, a, b, c, d);
    return 3.*a*_x*_x+2.*b*_x+c;
}
void InterpolatorSpline::Coefficients(value_type _x,
        value_type &a,
        value_type &b,
        value_type &c,
        value_type &d)
{
    pair_type lower(_x,0.);
    auto it = std::lower_bound(M_data.begin(), M_data.end(), lower, [] (pair_type a, pair_type b){ return a.first < b.first ; } );
    if(it == M_data.end()) // _x > M_data.end()->first 
    {
        a = sol((M_data.size()-2)*4+0);
        b = sol((M_data.size()-2)*4+1);
        c = sol((M_data.size()-2)*4+2);
        d = sol((M_data.size()-2)*4+3);
    }
    else if(it == M_data.begin()) // _x < M_data.begin()->first 
    {
        a = sol(0);
        b = sol(1);
        c = sol(2);
        d = sol(3);
    }
    else
    {
        size_t index = std::distance(M_data.begin(), it)-1;
        a = sol(index*4+0);
        b = sol(index*4+1);
        c = sol(index*4+2);
        d = sol(index*4+3);
    }
}

// http://www.iue.tuwien.ac.at/phd/rottinger/node60.html
/*
 * The Akima's interpolation differs with the Cspline by the first order derivative imposed at nodes
 * TODO: interpolate a(x-xi)**3 + ... instead of ax**3
 */
InterpolatorAkima::InterpolatorAkima( std::vector<pair_type> data) : super(data)
{
    SpMat A((this->M_data.size()-1)*4,(this->M_data.size()-1)*4);
    std::vector<T> coef;
    Eigen::VectorXd b = Eigen::VectorXd::Zero((this->M_data.size()-1)*4);
    double x, xp;
    for(int i = 0; i<this->M_data.size()-1; i++)
    {
        x  = this->M_data.at(i).first;
        xp = this->M_data.at(i+1).first;
        //si(xi) = yi
        coef.push_back(T(i*4,i*4+0, x*x*x));
        coef.push_back(T(i*4,i*4+1, x*x  ));
        coef.push_back(T(i*4,i*4+2, x    ));
        coef.push_back(T(i*4,i*4+3, 1.   ));
        b(i*4)   = M_data.at(i).second;
        //si(xi+1) = yi+1
        coef.push_back(T(i*4+1,i*4+0,xp*xp*xp));
        coef.push_back(T(i*4+1,i*4+1,xp*xp   ));
        coef.push_back(T(i*4+1,i*4+2,xp      ));
        coef.push_back(T(i*4+1,i*4+3,1.      ));
        b(i*4+1) = M_data.at(i+1).second;
        // si'(xi) = 
        coef.push_back(T(i*4+2,i*4+0, 3*x*x));
        coef.push_back(T(i*4+2,i*4+1, 2*x  ));
        coef.push_back(T(i*4+2,i*4+2, 1    ));
        b(i*4+2)   = sp(i);
        // si'(xi+1) = 
        coef.push_back(T(i*4+3,i*4+0,3*xp*xp));
        coef.push_back(T(i*4+3,i*4+1,2*xp   ));
        coef.push_back(T(i*4+3,i*4+2,1      ));
        b(i*4+3) = sp(i+1);

    }
    A.setFromTriplets(coef.begin(), coef.end());
    Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int>> solver(A);
    sol = solver.solve(b);
}
Interpolator::value_type InterpolatorAkima::operator()(double _x)
{
    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    Coefficients(_x, a, b, c, d);
    return a*_x*_x*_x+b*_x*_x+c*_x+d;
}
Interpolator::value_type InterpolatorAkima::diff(double _x)
{
    double a = 0.;
    double b = 0.;
    double c = 0.;
    double d = 0.;
    Coefficients(_x, a, b, c, d);
    return 3.*a*_x*_x+2.*b*_x+c;
}
Interpolator::value_type InterpolatorAkima::sp(int i)
{
    value_type dm2;
    value_type dm1;
    value_type di;
    value_type dp;

    value_type wm, w;

    dm2 = d(i-2);
    dm1 = d(i-1);
    di  = d(i  );
    dp  = d(i+1);
    wm = std::abs(dp-di);
    w  = std::abs(dm1-dm2);

    if(dm2 == dm1 && di != dp)
        return dm1;
    else if(di == dp && dm2 != dm1)
        return di;
    else if(dm1 == di)
        return di;
    else if (dm2==dm1 && di == dp && dp != dm2)
        return 0.5*(dm1+di);
    else
        return (wm*dm1+w*di)/std::max((wm+w),1e-12);
}
Interpolator::value_type InterpolatorAkima::d(int j)
{
    if(j < 0)
        return 2*d(j+1)-d(j+2);
    else if(j > M_data.size()-2 )
    return 2*d(j-1)-d(j-2);
    else
    {
        return (M_data.at(j+1).second - M_data.at(j).second) /  (M_data.at(j+1).first - M_data.at(j).first);
    }
}

void InterpolatorAkima::Coefficients(value_type _x,
        value_type &a,
        value_type &b,
        value_type &c,
        value_type &d)
{
    pair_type lower(_x,0.);
    auto it = std::lower_bound(M_data.begin(), M_data.end(), lower, [] (pair_type a, pair_type b){ return a.first < b.first ; } );
    if(it == M_data.end()) // _x > M_data.end()->first 
    {
        a = sol((M_data.size()-2)*4+0);
        b = sol((M_data.size()-2)*4+1);
        c = sol((M_data.size()-2)*4+2);
        d = sol((M_data.size()-2)*4+3);
    }
    else if(it == M_data.begin()) // _x < M_data.begin()->first 
    {
        a = sol(0);
        b = sol(1);
        c = sol(2);
        d = sol(3);
    }
    else
    {
        size_t index = std::distance(M_data.begin(), it)-1;
        a = sol(index*4+0);
        b = sol(index*4+1);
        c = sol(index*4+2);
        d = sol(index*4+3);
    }
}


