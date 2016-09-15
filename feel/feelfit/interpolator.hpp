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

#include <feel/feel.hpp>
#include <feel/feelalg/backendeigen.hpp>
#include <feel/feelfit/enums.hpp>

#ifndef __FEELPP_INTERPOLATOR_H
#define __FEELPP_INTERPOLATOR_H 1

/**
 * \class Interpolator
 * \brief  Given an un-ordered dataset, provide a representation for the interpolant
 * Options are available for the extrapolation
 * Given a two column file 'x - y', we compute the interpolation of that cloud.
 * We provide the first order derivative of the interpolant as well.
 */
class Interpolator
{
    public:
        using value_type = double;
        using pair_type = std::pair<value_type, value_type>;

        static std::unique_ptr<Interpolator> New( InterpolationType type, std::vector<pair_type> data);
        static std::unique_ptr<Interpolator> New( InterpolationType type, std::string dataFile);

        /**
         * Constructor
         * Sort the data
         * \param std::vector<pair<double,double>> data to interpolate
         */
        Interpolator(std::vector<pair_type> data) : M_data(data)
        {
            std::sort(M_data.begin(), M_data.end(), []( pair_type a, pair_type b){ return a.first < b.first; });
        }
        Interpolator(const Interpolator &) = default;
        Interpolator(Interpolator &&) = default;
        ~Interpolator() = default;
        /**
         *  Evaluate the interpolant
         */
        virtual value_type operator()(double _x)
        {
            return 0.;
        }
        /**
         *  Evaluate the interpolant derivative
         */
        virtual value_type diff(double _x)
        {
            return 0.;
        }
        /**
         * Return the type of the interpolant (P0, P1, Akima, Spline)
         */
        virtual int type(void) { return -1; }
        /**
         * Access to the data
         */
        std::vector<std::pair<value_type, value_type>> data() { return M_data;}
    protected:
        std::vector<std::pair<value_type, value_type>> M_data;
};// class interpolBase

class InterpolatorP0 : public Interpolator
{
    public:
        typedef Interpolator super;
        using value_type = double;
        using pair_type = std::pair<value_type, value_type>;
        /**
         * P0-style interpolant
         * \param data the Dataset
         * \param iType left, right or center
         */
        InterpolatorP0( std::vector<pair_type> data, InterpolationType_P0 iType = left ) : super(data), M_iType(iType) {}
        InterpolatorP0(const InterpolatorP0 &) = default;
        /**
         *  Evaluate the interpolant derivative
         */
        value_type diff(double _x); 
        /**
         *  Evaluate the interpolant
         */
        value_type operator()(double _x);
        int type(void) { return 0; }
    private:
        InterpolationType_P0 M_iType;

}; // InterpolatorP0

class InterpolatorP1 : public Interpolator
{
    public:
        typedef Interpolator super;
        using value_type = double;
        using pair_type = std::pair<value_type, value_type>;
        /**
         * Interpolant P1 style
         * Extrapolation method at left or right can be different
         * \param r right kind of extrapolation
         * \param l left kind of extrapolation
         */
        InterpolatorP1(const InterpolatorP1 &) = default;
        InterpolatorP1( std::vector<pair_type> data, ExtrapolationType_P1 r = zero, ExtrapolationType_P1 l = zero ) : super(data), M_r_type(r), M_l_type(l){}
        /**
         *  Evaluate the interpolant
         */
        value_type operator()(double _x);
        /**
         *  Evaluate the interpolant derivative
         */
        value_type diff(double _x) ;
        int type(void) { return 1; }
    private:
        void Coefficients(value_type _x, 
                value_type &y1,
                value_type &y2,
                value_type &x1,
                value_type &x2);
        ExtrapolationType_P1 M_r_type;
        ExtrapolationType_P1 M_l_type;

}; // InterpolatorP0

/**
 * \class InterpolatorSpline 
 * \brief Interpolate with a*x**3 + b*x**2 + c*x + d
 * TODO: interpolate a(x-xi)**3 + ... instead of ax**3 (only to generate a easiest to inverse linear system)
 */
class InterpolatorSpline : public Interpolator
{
    public:
        typedef Interpolator super;
        using value_type = double;
        using pair_type = std::pair<value_type, value_type>;
        using SpMat = Eigen::SparseMatrix<value_type>; 
        using T = Eigen::Triplet<value_type>;
        InterpolatorSpline( std::vector<pair_type> data, ExtrapolationType_spline r , ExtrapolationType_spline l );
        InterpolatorSpline(const InterpolatorSpline &) = default;
        ~InterpolatorSpline() = default; // the vector is cleared by Eigen
        value_type operator()(double _x);
        value_type diff(double _x);
        int type(void) { return 2; }

    private:
        ExtrapolationType_spline M_r_type;
        ExtrapolationType_spline M_l_type;

        Eigen::VectorXd sol;

        void Coefficients(value_type _x,
                value_type &a,
                value_type &b,
                value_type &c,
                value_type &d);

}; // InterpolatorSpline

/**
 * \class InterpolatorAkima
 * \brief The Akima's interpolation differs with the Cspline by the first order derivative imposed at nodes
 *  http://www.iue.tuwien.ac.at/phd/rottinger/node60.html
 * TODO: interpolate a(x-xi)**3 + ... instead of ax**3
 */
class InterpolatorAkima : public Interpolator
{
    public:
        typedef Interpolator super;
        using value_type = double;
        using pair_type = std::pair<value_type, value_type>;
        using SpMat = Eigen::SparseMatrix<value_type>; 
        using T = Eigen::Triplet<value_type>;
        InterpolatorAkima( std::vector<pair_type> data) ;
        InterpolatorAkima(const InterpolatorAkima &) = default;
        ~InterpolatorAkima() = default; // the vector is cleared by Eigen
        value_type operator()(double _x);
        value_type diff(double _x);
        int type(void) { return 3; }
    private:
        /**
         * See provided links for explanation
         */
        value_type sp(int i);
        /**
         * See provided links for explanation
         */
        value_type d(int j);
        Eigen::VectorXd sol;

        void Coefficients(value_type _x,
                value_type &a,
                value_type &b,
                value_type &c,
                value_type &d);

}; // InterpolatorAkima

std::unique_ptr<Interpolator> Interpolator::New( InterpolationType type, std::vector<pair_type> data)
{
    switch(type)
    {
        case P0:
            {
                return std::unique_ptr<InterpolatorP0>(new InterpolatorP0(data, static_cast<InterpolationType_P0>(Feel::ioption("fit.P0"))));
                break;
            }
        case P1:
            {
                return std::unique_ptr<InterpolatorP1>(new InterpolatorP1(data,
                            static_cast<ExtrapolationType_P1>(Feel::ioption("fit.P1_right")), 
                            static_cast<ExtrapolationType_P1>(Feel::ioption("fit.P1_left"))
                            ));
                break;
            }
        case Spline:
            {
                return std::unique_ptr<InterpolatorSpline>(new InterpolatorSpline(data, 
                            static_cast<ExtrapolationType_spline>(Feel::ioption("fit.Spline_right")),
                            static_cast<ExtrapolationType_spline>(Feel::ioption("fit.Spline_left"))
                            ));
                break;
            }
        case Akima:
            {
                return std::unique_ptr<InterpolatorAkima>(new InterpolatorAkima(data));
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
    if (infile.is_open()) {
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
    else{
        Feel::cout << "Error opening file" << Feel::Environment::expand(dataFile) << "\n";
    }
    return std::move(Interpolator::New( type, data));
}

#endif //__FEELPP_INTERPOLATOR_H
