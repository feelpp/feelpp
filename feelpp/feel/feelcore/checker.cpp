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
//! @date 24 Jul 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/checker.hpp>
#include <feel/feelmath/polyfit.hpp>
#include <feel/feelmath/vector.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

namespace Feel {

namespace rate {

hp::hp( double h, int p )
    :
    M_h(h),
    M_p(p)

{
    bool verbose = boption(_name="checker.verbose");
    fs::path path ( Environment::expand( soption(_name="checker.filename" ) ) );
    if ( fs::exists( path ) )
    {
        std::ifstream f ( path.string().c_str() );
        pt::read_json( f, M_ptree);
        if ( verbose )
            pt::write_json( std::cout,  M_ptree );
        // build data structure
        for( auto itfn = M_ptree.begin(); itfn != M_ptree.end(); ++itfn )
        {
            auto const& fn = *itfn;
            for( auto it = fn.second.begin(); it != fn.second.end(); ++it )
            {
                auto const& v = *it;

                int order = boost::lexical_cast<int>(v.first);

                M_data[fn.first][order].exact = v.second.get("exact",false);
                if ( M_data[fn.first][order].exact )
                    for( auto o : irange( order, 20 ) )
                        M_data[fn.first][o].exact = true;
                if ( M_data[fn.first][order].exact == false )
                {
                    std::vector<double> hs;
                    for (pt::ptree::value_type const &h : v.second.get_child("h"))
                    {
                        // hs stores an array of h
                        hs.push_back(boost::lexical_cast<double>(h.second.data()));
                    }
                    LOG(INFO) << "hs: " << hs;
                    std::map<std::string,data::errors_t> errors;
                    for (pt::ptree::value_type const &error : v.second.get_child("errors"))
                    {
                        errors[error.first].order = error.second.get<double>( "order" );
                        // error values stores a set of arrays associated to one or more error norms

                        for ( pt::ptree::value_type const &e : error.second.get_child("values" ) )
                        {
                            errors[error.first].values.push_back(boost::lexical_cast<double>(e.second.data()) );
                        }
                        LOG(INFO) << error.first << " error norms(" << errors[error.first].order << "): " << errors[error.first].values;
                    }


                    M_data[fn.first][order] = data( std::move(hs), std::move( errors ) );
                }
                if ( M_data[fn.first][order].hs.size() >= 3 )
                {
                    for( auto const& error: M_data[fn.first][order].errors )
                    {
                        auto c = polyfit( log(M_data[fn.first][order].hs), log(error.second.values), 1 );
                        LOG(INFO) << "order = " << error.second.order << " c[1]=" << c[1] << std::endl;
                        if ( c[1] >= error.second.order-1e-1 )
                            LOG(INFO) << "db for " << error.first << " error norm is consistent";
                        else
                            throw std::logic_error( std::string("inconsistent convergence database, expected ")+std::to_string(error.second.order)+" got " + std::to_string(c[1]) );
                    }
                }
            }
        }
    }


}

Checks
hp::operator()( std::string const& sol, std::pair<std::string,double> const& r, double otol, double etol )
{
    std::string solution = sol;
    boost::trim_right_if(solution,boost::is_any_of(":"));
    boost::erase_all(solution, " ");   
    // check that we have some data on solution for the specific order
    if ( M_data.count(solution) && M_data.at(solution).count(M_p) && ( M_data.at(solution).at(M_p).exact || M_data.at(solution).at(M_p).errors.count(r.first) ) )
    {
        hp::data d = M_data.at(solution).at(M_p);

        if ( d.exact )
        {

            if ( r.second < etol )
                return Checks::EXACT;
            throw CheckerExactFailed( r.second, etol );

        }
        else
        {
            if ( d.hs.size() > 1 )
            {
                d.hs.push_back(M_h);
                d.errors.at(r.first).values.push_back(r.second);
                auto c = polyfit( log(d.hs), log(d.errors.at(r.first).values), 1 );
                LOG(INFO) << "order = " << d.errors.at(r.first).order << " c[1]=" << c[1] << std::endl;
                //std::cout << "order = " << d.errors.at(r.first).order << " c[1]=" << c[1] << std::endl;
                if (  c[1]>=d.errors.at(r.first).order - otol )
                    return Checks::CONVERGENCE_ORDER;
                throw CheckerConvergenceFailed( d.errors.at(r.first).order, c[1], otol );
            }
            return Checks::NONE;
        }

    }
#if 0
    pt::ptree p_node;
    pt::ptree h_node;
    h_node.put( "", std::to_string(M_h));
    p_node.add_child( "h", h_node );
    pt::ptree error_node;
    error_node.put("", r );
    p_node.add_child( "error", error_node );
    pt::ptree node;
    node.add_child( std::to_string(M_p), p_node );
    pt::ptree fnode;
    fnode.add_child( solution, node );
    //M_ptree.add_child( std::to_string(M_p), node );
    pt::write_json(soption(_name="checker.name")+".json", fnode);
    //pt::write_json(std::cout, node);
#endif
    return Checks::NONE;
}


} // rate
#if 0
Checker
checker( std::string const& c_name, std::string const& solution, std::string const& gradient )
{
    if ( solution.empty() && soption(_name="checker.solution" ).empty() )
        throw std::logic_error("Invalid setup of Checker system, no solution provided");
    auto sol = solution.empty()?soption(_name="checker.solution" ):solution;
    Checker c{c_name};
    c.setSolution( solution );
    if ( !gradient.empty() || !soption(_name="checker.gradient" ).empty() )
        c.setGradient( gradient.empty()?soption(_name="checker.gradient" ):gradient );
    return c;
}
#endif
Checker::Checker( std::string const& name )
    :
    super( "Checker", name ),
    M_check( boption(_name="checker.check") ),
    M_verbose( boption(_name="checker.verbose" ) ),
    M_solution( soption(_name="checker.solution" ) ),
    M_etol( doption(_name="checker.tolerance.exact" ) ),
    M_otol( doption(_name="checker.tolerance.order" ) )
{
}

void
Checker::setScript( std::string const& s, variables_t const& in, std::map<std::string,double> const& p, bool u )
{
    M_use_script = u;
    M_script_in = in;
    M_script = s;
    std::cout << fmt::format( "script: {}", s ) << std::endl;
    M_param_values= p;
    std::cout << fmt::format( "param_values: {}", M_param_values ) << std::endl;
}
Checker::variables_t
Checker::runScript()
{
    variables_t locals{ M_script_in };
    locals[M_solution_key]=M_solution;
    if ( M_gradient )
        locals[M_gradient_key]=*M_gradient;
    else
        locals[M_gradient_key]="";
    std::cout << "gradient(" << M_gradient_key << "):" << locals[M_gradient_key] << std::endl;
    locals["compute_pde_coefficients"]=boption(_name="checker.compute-pde-coefficients")?"true":"false";
    Feel::pyexprFromFile( Environment::expand(M_script), locals );
    std::cout << "gradient 2:" << locals[M_gradient_key] << std::endl;
    M_solution=locals[M_solution_key];
    M_gradient=locals[M_gradient_key];
    return locals;
}

}
