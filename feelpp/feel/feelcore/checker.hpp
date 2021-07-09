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

#include <optional>
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
    using variables_t = std::map<std::string,std::string>;

    explicit Checker( std::string const& name = "" );
    Checker( Checker const& c ) = default;
    Checker( Checker && c ) = default;

    //!
    //! use @param s as solution to check the numerical results
    //!
    Checker & operator=( Checker const& c ) = default;
    Checker & operator=( Checker && c ) = default;
    ~Checker() override = default;

    bool check() const { return M_check; }
    void setCheck( bool c ) { M_check = c; }
    bool verbose() const { return M_verbose; }
    void setVerbose( bool v ) { M_verbose = v; }


    /*
     * setter and getter for solution expression
     */
    std::string const& solution() const { return M_solution; }
    void setSolution( std::string const& s, std::string const& key = "solution" ) { M_solution = s; M_solution_key = key; }
    
    /*
     * setter and getter for gradient expression
     */
    bool hasGradient() const { return M_gradient.has_value(); }
    std::optional<std::string> const& gradient() const { return M_gradient; }
    void setGradientKey( std::string const& key ) { M_gradient_key = key; }
    std::string const& gradientKey() const { return M_gradient_key; } 
    void setGradient( std::string const& s, std::string const& key = "gradient" ) { M_gradient = s; M_gradient_key = key; }

    //! set the script to configure a PDE from a set of coefficient using a script
    void setScript( std::string const& s, variables_t const& inputs = {} , std::map<std::string,double> const& params = {}, bool use = true );

    //! @return true if use the script, false otherwise
    void setUseScript( bool u ) { M_use_script = u; }
    //! @return true if use the script, false otherwise
    bool useScript() const { return M_use_script; }

    void setParameterValues( std::map<std::string,double> const& p ) { M_param_values = p; }
    std::map<std::string,double> const& parameterValues() const { return M_param_values; }
    /*
     * compute manufactured solution from a python script
     * @param vm variables map
     * @return an optional variables map
     */
    variables_t runScript();

    template<typename ErrorFn, typename ErrorLaw>
    int
    runOnce( ErrorFn fn, ErrorLaw law, std::string metric = "||u-u_h||_" );

private:
    bool M_check;
    bool M_verbose;
    std::string M_solution;
    std::optional<std::string> M_gradient;
    std::string M_solution_key{"solution"}, M_gradient_key{"gradient"};
    double M_etol, M_otol;
    std::string M_script;
    std::map<std::string,std::string> M_script_in;
    std::map<std::string,double> M_param_values;
    bool M_use_script;
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
    nl::json pt;
    for( auto const& e : err )
    {
        pt.emplace( e.first, e.second );
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
//! @param c_name name of the checker
//! @param solution to pass to checker and check the numerical results, if empty use checker.solution
//! @param gradient optional expression for the gradient, if empty use checker.gradient
//! if @p gradient and option checker.gradient are empty then we do not set the gradient expression
//!
BOOST_PARAMETER_FUNCTION(
    ( Checker ), // return type
    checker,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( name, (std::string))
      ( solution_key, (std::string))
      ) // 4. one required parameter, and

    ( optional
      ( inputs,(std::map<std::string,std::string>), (std::map<std::string,std::string>{})  )
      ( parameter_values,(std::map<std::string,double>), (std::map<std::string,double>{})  )
      ( solution, (std::string), inputs.count(solution_key)?inputs.at(solution_key):soption("checker.solution"))
      ( gradient_key, (std::string), std::string{"grad_"}+solution_key)
      ( gradient, (std::string), inputs.count(gradient_key)?inputs.at(gradient_key):soption("checker.gradient"))
      ( script,(std::string), soption("checker.script"))    
      ( compute_pde_coefficients, (bool), boption("checker.compute-pde-coefficients") )
      ( prefix, (std::string), "" )
      ( verbose,   (bool), boption(_prefix=prefix,_name="checker.verbose") )
      )
    )
{
    using Feel::cout;
    
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsequenced"
#endif
    //if ( solution.empty() && soption("checker.solution" ).empty() )
    //    throw std::logic_error("Invalid setup of Checker system, no solution provided");
    Checker c{name};
    c.setVerbose( verbose );
    c.setSolution( solution, solution_key );
    if ( !gradient.empty() )
        c.setGradient( gradient, gradient_key );
    else
        c.setGradientKey( gradient_key );
    c.setScript( script, inputs, parameter_values, compute_pde_coefficients ); 
    return c;
}

} // Feel

#endif
