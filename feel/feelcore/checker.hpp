//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 23 Jul 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_CHECKER_HPP
#define FEELPP_CHECKER_HPP 1

#include <feel/feelcore/environment.hpp>


namespace Feel {

namespace rate {

class FEELPP_EXPORT hp
{
public:
    //!
    //! 
    //!
    hp( double h, int p );

    //!
    //! 
    //!
    bool operator()( std::string const& solution, std::pair<std::string,double> const& r, double otol = 1e-1, double etol=1e-15 );
private:
    struct FEELPP_NO_EXPORT data
    {
        struct errors_t
        {
            double order;
            std::vector<double> values;
        };
        data() = default;
        data( std::vector<double> const& _hs, std::map<std::string,errors_t> const& _errors, bool _exact = false )
            :
            hs( _hs ),
            errors( _errors ),
            exact(_exact)
            {}
        data( std::vector<double> && _hs, std::map<std::string,errors_t> && _errors, bool _exact = false )
                :
                hs( _hs ),
                errors( _errors ),
                exact(_exact)
            {}
        ~data() = default;
            
        std::vector<double> hs;
        std::map<std::string,errors_t> errors;
        bool exact = false;
    };
private:
    double M_h;
    int M_p;
    int M_expected;
    pt::ptree M_ptree;
    // [fonction][order]{.hs,.errors}
    std::map<std::string,std::map<int, data>> M_data;
};

}
//!
//! Checker class
//!
//! for PDE solve When solving PDEs, we want to check that we get
//! proper convergence rates using a priori error estimates.
//!
//! The \c Checker class 
//!
class FEELPP_EXPORT Checker
{
public:
    Checker();
    Checker( Checker const& c ) = default;
    Checker( Checker && c ) = default;
    Checker & operator=( Checker const& c ) = default;
    Checker & operator=( Checker && c ) = default;
    ~Checker() = default;

    bool check() const { return M_check; }
    void setCheck( bool c ) { M_check = c; }
    
    std::string const& solution() const { return M_solution; }
    void setSolution( std::string const& s ) { M_solution = s; }
    
    template<typename ErrorFn, typename ErrorLaw>
    int
    runOnce( ErrorFn fn, ErrorLaw law );

private:
    bool M_check;
    std::string M_solution;
    double M_etol, M_otol;
};


template<typename ErrorFn, typename ErrorRate>
int
Checker::runOnce( ErrorFn fn, ErrorRate rate )
{
    if ( M_check )
    {
        
        auto err = fn( M_solution );
        bool status = true;
        for( auto const& e : err )
        {
            cout << "||u-u_h||_" << e.first << "=" << e.second  << std::endl;
            if ( rate(M_solution, e, M_otol, M_etol) )
            {
                cout << tc::green << e.first << " error norm check successful: " << e.second << tc::reset << std::endl;
            }
            else
            {
                cout << tc::red << e.first <<  " error norm check failed: " << e.second << tc::reset << std::endl;
                status = false;
            }
        }
        return status;
    }

    return 0;
}

//!
//! check create function
//!
Checker checker( std::string const& p  = "" );

} // Feel

#endif
