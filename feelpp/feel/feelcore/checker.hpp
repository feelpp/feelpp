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

enum class Checks
{
    NONE=0,
    EXACT,
    CONVERGENCE_ORDER
};

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
    Checks operator()( std::string const& solution, std::pair<std::string,double> const& r, double otol = 1e-1, double etol=1e-13 );
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

} // rate

struct CheckerConvergenceFailed : public std::logic_error
{
    CheckerConvergenceFailed() = delete;

    CheckerConvergenceFailed( double expected, double got, double tol )
        :
        std::logic_error( "Checker convergence order  verification failed" ),
        M_expected_order( expected ),
        M_got_order( got ),
        M_tol( tol )
        {}
    double expectedOrder() const { return M_expected_order; }
    double computedOrder() const { return M_got_order; }
    double tolerance() const { return M_tol; }
private:
    double M_expected_order, M_got_order, M_tol;
};
struct CheckerExactFailed : public std::logic_error
{
    CheckerExactFailed() = delete;
    CheckerExactFailed( double got, double tol )
        :
        std::logic_error( "Checker exact verification failed" ),
        M_got_error( got ),
        M_tol( tol )
        {}
    double computedError() const { return M_got_error; }
    double tolerance() const { return M_tol; }
private:
    double M_got_error, M_tol;
};

//!
//! Checker class
//!
//! for PDE solve When solving PDEs, we want to check that we get
//! proper convergence rates using a priori error estimates.
//!
//! The \c Checker class
//!
class FEELPP_EXPORT Checker : public JournalWatcher
{
    using super = JournalWatcher;
public:
    explicit Checker( std::string const& name = "" );
    Checker( Checker const& c ) = default;
    Checker( Checker && c ) = default;

    //!
    //! use @param s as solution to check the numerical results
    //!
    Checker & operator=( Checker const& c ) = default;
    Checker & operator=( Checker && c ) = default;
    ~Checker() = default;

    bool check() const { return M_check; }
    void setCheck( bool c ) { M_check = c; }
    bool verbose() const { return M_verbose; }
    void setVerbose( bool v ) { M_verbose = v; }

    std::string const& solution() const { return M_solution; }
    void setSolution( std::string const& s ) { M_solution = s; }

    template<typename ErrorFn, typename ErrorLaw>
    int
    runOnce( ErrorFn fn, ErrorLaw law, std::string metric = "||u-u_h||_" );

private:
    bool M_check;
    bool M_verbose;
    std::string M_solution;
    double M_etol, M_otol;
};


template<typename ErrorFn, typename ErrorRate>
int
Checker::runOnce( ErrorFn fn, ErrorRate rate, std::string metric )
{
    if ( !M_check )
        return true;
    cout << "================================================================================\n"
         << "[Checker] " << this->journalWatcherInstanceName() << "\n";
    auto err = fn( M_solution );
    std::vector<bool> checkSuccess;
    pt::ptree pt;
    for( auto const& e : err )
    {
        pt.put( e.first, e.second );
        //cout << "||u-u_h||_" << e.first << "=" << e.second  << std::endl;
        
        //cout << "||u-u_h||_" << e.first << "=" << e.second  << std::endl;
        try
        {
            Checks c = rate(M_solution, e, M_otol, M_etol);
            
            switch( c ) 
            {
            case Checks::NONE:
                cout << tc::yellow << "[no checks]" << metric << e.first << "=" << e.second << tc::reset << std::endl;
                break;
            case Checks::EXACT:
                cout << tc::green << "[exact verification success]" << metric << e.first <<  "=" << e.second << tc::reset << std::endl;
                break;
            case Checks::CONVERGENCE_ORDER:
                cout << tc::green << "[convergence order verification success]" << metric << e.first <<  "=" << e.second << tc::reset << std::endl;
                break;
            }
            checkSuccess.push_back( true );
        }
        catch( CheckerConvergenceFailed const& ex )
        {
            cout << tc::red
                 << "Checker convergence order verification failed for " << metric << e.first << std::endl
                 << " Solution " << M_solution << std::endl	
                 << " Computed order " << ex.computedOrder() << std::endl
                 << " Expected order " << ex.expectedOrder() << std::endl
                 << " Tolerance " << ex.tolerance()  << tc::reset << std::endl;
            checkSuccess.push_back( false );
        }
        catch( CheckerExactFailed const& ex )
        {
            cout << tc::red
                 << "Checker exact verification failed for " << metric << e.first << std::endl
                 << " Solution " << M_solution << std::endl	
                 << " Computed error " << ex.computedError() << std::endl
                 << " Tolerance " << ex.tolerance()  << tc::reset << std::endl;
            checkSuccess.push_back( false );
        }
        catch( std::exception const& ex )
        {
            cout << tc::red << "Checker Caught exception " << ex.what() << tc::reset << std::endl
                 << " Solution " << M_solution << std::endl;
            checkSuccess.push_back( false );
        }
    }
    this->putInformationObject( pt );
    return ( std::find(checkSuccess.begin(),checkSuccess.end(), false ) == checkSuccess.end() );
}

//!
//! check create function
//! @param s solution to pass to checker and check the numerical results
//!
FEELPP_EXPORT Checker checker( std::string const& s  = "", std::string const& prefix = "" );

} // Feel

#endif
