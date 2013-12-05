/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   Date: 2012-06-02

   Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file error.hpp
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \date 2012-06-02
*/

#ifndef __ERROR_HPP
#define __ERROR_HPP 1

#include <string>
#include <sstream>

#include <boost/range/algorithm/for_each.hpp>
#include <boost/algorithm/string/split.hpp>

#include <boost/parameter.hpp>
#include <feel/feelcore/parameter.hpp>

#include <feel/feel.hpp>

// cf bdf class for instance :

// should add default option
// params
// exact
// rhs
// convergence_study

// should
// compute error(s) (both norm and "error")
// run cvg study
// display result (gnuplot_iostream?)

// should
// run external code
// read output from code
// display some error indicator

namespace Feel
{
    using  GiNaC::symbol;
    using  GiNaC::lst;
    using  GiNaC::ex;

    BOOST_PARAMETER_NAME(Xh)
    BOOST_PARAMETER_NAME(T)
    BOOST_PARAMETER_NAME(exact_solution)
    BOOST_PARAMETER_NAME(parameters)
    BOOST_PARAMETER_NAME(convergence)
    BOOST_PARAMETER_NAME(convergence_steps)

    /**
     * \class ErrorBase
     *
     * Provides abstract interface to error calculations
     *
     * \tparam Dim dimension of space
     */
    template<int Dim, int Order>
    class ErrorBase
    {
    public:

        //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
        typedef Simplex<Dim> convex_type;
        //! mesh type
        typedef Mesh<convex_type> mesh_type;
        //! mesh shared_ptr<> type
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

        //! the basis type of our approximation space
        typedef bases<Lagrange<Order,Scalar> > basis_type;
        //! the approximation function space type
        typedef FunctionSpace<mesh_type, basis_type> space_type;
        //! the approximation function space type (shared_ptr<> type)
        typedef boost::shared_ptr<space_type> space_ptrtype;
        //! an element type of the approximation function space
        typedef typename space_type::element_type element_type;

        /** @name Constructors, destructor
         */
        //@{
        ErrorBase( po::variables_map const& vm, std::string const& prefix)
            :
            M_exact( vm[prefixvm( prefix, "error.exact" )].as<std::string>() ),
            M_exact_params( vm[prefixvm( prefix, "error.params" )].as<std::string>() ),
            M_rhs( vm[prefixvm( prefix, "error.rhs" )].as<std::string>() ),
            M_rhs_computed(vm[prefixvm( prefix, "error.rhs.computed" )].as<bool>() ),
            M_convergence(vm[prefixvm( prefix, "error.convergence" )].as<bool>() ),
            M_convergence_max(vm[prefixvm( prefix, "error.convergence.steps" )].as<int>() )
        {
        }

        ErrorBase( std::string const& name, std::string const& params)
            :
            M_exact( name ),
            M_exact_params(params),
            M_rhs(""),
            M_rhs_computed( false ),
            M_convergence(0),
            M_convergence_max(0)
        {
        }

        ErrorBase( ErrorBase const &e )
        :
        M_exact( e.M_exact ),
        M_exact_params( e.M_exact_params ),
        M_rhs( e.M_rhs ),
        M_rhs_computed( e.M_rhs_computed ),
        M_convergence( e.M_convergence ),
        M_convergence_max( e.M_convergence_max )
        {
        }

        ErrorBase& operator=( ErrorBase const& e )
        {
            if ( this != &e )
                {
                    M_exact = e.M_exact ;
                    M_exact_params = e.M_exact_params;
                    M_rhs = e.M_rhs;
                    M_rhs_computed = e.M_rhs_computed;
                    M_convergence = e.M_convergence;
                    M_convergence_max = e.M_convergence_max;
                }

            return *this;
        }

        virtual ~ErrorBase() {}
        //@}


        /** \name  Mutators
         */
        //@{
        //! set the parameters string
        void setParams ( std::string params )
        {
            M_exact_params=params;
        }

        //! set the rhs_computed flag
        void setComputedRhs ( bool doComputeRhs )
        {
            M_rhs_computed=doComputeRhs;
        }
        //! set the convergence flag
        void setConvergence ( bool doConvergence )
        {
            M_convergence=doConvergence;
        }

        //! set the convergence interations
        void setnumberOfConvergenceSteps( int n)
        {
            M_convergence_max=n;
        }

        void setSolution(std::string const& expression ="", std::string const& p ="")
        {
            if ( expression.empty() )
                FEELPP_ASSERT( M_rhs.size() ).error( "undefined rhs exact" );
            else
                M_exact = expression;

            if ( !p.empty() )
                M_exact_params = p;

            vars = symbols<Dim>();
            std::vector<std::string> lst_params;
            boost::split(lst_params, M_exact_params, boost::is_any_of(";"));
            LOG(INFO) << "Loading additionnal parameters : ";
            boost::for_each(lst_params, [](const std::string &s){LOG(INFO) << s << ", ";});
            LOG(INFO) << std::endl;
            parameters = symbols(lst_params);

            LOG(INFO) << "Loading function : " << M_exact << std::endl;
            boost::for_each( parameters, [](symbol const& s ) {LOG(INFO) << "Symbol " << s.get_name() << " found\n";});
            exact = parse(M_exact, vars, parameters);
            LOG(INFO) << "exact solution is : " << exact << "\n";
        }

        void setRhs(std::string const& expression ="")
        {
            if ( expression.empty() )
                FEELPP_ASSERT( M_rhs.size() ).error( "undefined rhs exact" );
            else
                M_rhs = expression;

            LOG(INFO) << "Loading rhs function : " << M_rhs << std::endl;
            rhs = parse(M_rhs, vars, parameters);
            LOG(INFO) << "rhs is : " << rhs << "\n";
        }

        typedef typename boost::function<ex (ex u, std::vector<symbol> vars, std::vector<symbol> p)> t_edp_type;
        typedef typename boost::function<ex (ex u, std::vector<symbol> vars)> edp_type;

        void setRhs(t_edp_type * edp)
        {
            FEELPP_ASSERT( parameters.size() ).error( "setRhs: inconsistant numbers of parameters" );
            rhs = (*edp)(exact, vars, parameters);
            LOG(INFO) << "computed rhs is : " << rhs << "\n";

            std::ostringstream rhs_expression;
            rhs_expression << rhs;
            M_rhs = rhs_expression.str();
        }

        void setRhs(edp_type * edp)
        {
            rhs = (*edp)(exact, vars);
            LOG(INFO) << "computed rhs is : " << rhs << "\n";

            std::ostringstream rhs_expression;
            rhs_expression << rhs;
            M_rhs = rhs_expression.str();
        }
        //@}

        /** \name Accessors
         */
        //@{

        /**
         * \return the computedrhs flag
         */
        bool computedrhs() const
        {
            return M_rhs_computed;
        }

        /**
         * \return the convergence flag
         */
        bool convergence() const
        {
            return M_convergence;
        }

        /**
         * \return the number of convergence steps
         */
        int numberOfConvergenceSteps() const
        {
            return M_convergence_max;
        }

        /**
         * \return the real space variables  GiNaC symbols
         */
        std::vector<symbol> getVars() const
        {
            return vars;
        }

        /**
         * \return the real spaceadditional parameters GiNaC symbols
         */
        std::vector<symbol> getParams() const
        {
            return parameters;
        }


        /**
         * \return GiNaC exact solution
         */
        std::string getExactSolution() const
        {
            return M_exact;
        }

        /**
         * \return GiNaC rhs
         */
        std::string getExactRhs() const
        {
            return M_rhs;
        }

        /**
         * \return GiNaC rhs
         */
        ex getRhs() const
        {
            // check if exact is defined
            FEELPP_ASSERT( M_rhs.size() ).error( "undefined rhs" );
            return rhs;
        }

        /**
         * \return GiNaC rhs
         */
        ex getRhs(const double & value) const
        {
            // check if exact is defined
            FEELPP_ASSERT( M_rhs.size() ).error( "undefined rhs" );

            // check if values is consistant with params
            FEELPP_ASSERT( parameters.size() == 1 ).error( "inconsistant values and parameters size" );

            ex tmp_sol = rhs;
            tmp_sol = substitute(tmp_sol, parameters[0], value);

            return tmp_sol;
        }

        /**
         * \return GiNaC rhs
         */
        ex getRhs(const std::vector<double> & values) const
        {
            // check if exact is defined
            FEELPP_ASSERT( M_rhs.size() ).error( "undefined exact solution" );

            // check if values is consistant with params
            FEELPP_ASSERT( values.size() == parameters.size() ).error( "inconsistant values and parameters size" );

            ex tmp_sol = rhs;
            for (unsigned int i=0; i<values.size(); i++) // may use boost zip iterator
                tmp_sol = substitute(tmp_sol, parameters[i], values[i]);

            return tmp_sol;
        }

        /**
         * \return GiNaC exact solution
         */
        ex getSolution() const
        {
            // check if exact is defined
            FEELPP_ASSERT( M_exact.size() ).error( "undefined exact solution" );
            return exact;
        }

        /**
         * \return GiNaC exact solution
         */
        ex getSolution(const double & value) const
        {
            // check if exact is defined
            FEELPP_ASSERT( M_exact.size() ).error( "undefined exact solution" );

            // check if values is consistant with params
            FEELPP_ASSERT( parameters.size() == 1 ).error( "inconsistant values and parameters size" );

            ex tmp_sol = exact;
            tmp_sol = substitute(tmp_sol, parameters[0], value);

            return tmp_sol;
        }

        /**
         * \return GiNaC exact solution
         */
        ex getSolution(const std::vector<double> & values) const
        {
            // check if exact is defined
            FEELPP_ASSERT( M_exact.size() ).error( "undefined exact solution" );

            // check if values is consistant with params
            FEELPP_ASSERT( values.size() == parameters.size() ).error( "inconsistant values and parameters size" );

            ex tmp_sol = exact;
            for (unsigned int i=0; i<values.size(); i++) // may use boost zip iterator
                tmp_sol = substitute(tmp_sol, parameters[i], values[i]);

            return tmp_sol;
        }
        //@}

        void print() const
        {
            LOG(INFO) << "============================================================\n";
            LOG(INFO) << "Error Information\n";
            LOG(INFO) << "   exact : " << getSolution() << "\n";
            LOG(INFO) << "   params : ";
            boost::for_each( parameters, [](symbol const& s ) {LOG(INFO) << s.get_name() << " ";});
            LOG(INFO) << "\n";
            if ( !M_rhs.empty() )
                {
                    LOG(INFO) << "   rhs : " << getRhs() << "\n";
                }
            LOG(INFO) << "   convergence : " << convergence() << "\n";
            LOG(INFO) << "   convergence steps : " << numberOfConvergenceSteps() << "\n";
        }

        /****************/
        element_type rhs_project(const space_ptrtype & Xh)
        {
            ex solution = getRhs();

            auto mesh = Xh->mesh();
            auto gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=expr(solution,vars) );
            return gproj;
        }

        element_type exact_project(const space_ptrtype & Xh)
        {
            ex solution = getSolution();

            auto mesh = Xh->mesh();
            auto gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=expr(solution,vars) );
            return gproj;
        }

        element_type grad_exact_project(const space_ptrtype & Xh)
        {
            ex solution = getSolution();
            auto gradg = expr<1,Dim,2>(grad(solution,vars), vars );

            auto mesh = Xh->mesh();
            auto gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=expr(gradg,vars) );
            return gproj;
        }

        element_type error_project(const space_ptrtype & Xh, const element_type & T)
        {
            auto gproj = exact_project(Xh);
            return (gproj - T);
        }

        element_type grad_error_project(const space_ptrtype & Xh, const element_type & T)
        {
            auto gproj = grad_exact_project(Xh);
            return (gproj - gradv(T));
        }

        element_type rhs_project(const space_ptrtype & Xh, const double & val)
        {
            ex solution = getRhs(val);

            auto mesh = Xh->mesh();
            auto gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=expr(solution,vars) );
            return gproj;
        }

        element_type exact_project(const space_ptrtype & Xh, const double & val)
        {
            ex solution = getSolution(val);

            auto mesh = Xh->mesh();
            auto gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=expr(solution,vars) );
            return gproj;
        }

        element_type grad_exact_project(const space_ptrtype & Xh, const double & val)
        {
            ex solution = getSolution(val);
            auto gradg = expr<1,Dim,2>(grad(solution,vars), vars );

            auto mesh = Xh->mesh();
            auto gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=expr(gradg,vars) );
            return gproj;
        }

        element_type error_project(const space_ptrtype & Xh, const element_type & T, const double & val)
        {
            auto gproj = exact_project(Xh, val);
            return (gproj - T);
        }

        element_type grad_error_project(const space_ptrtype & Xh, const element_type & T, const double & val)
        {
            auto gproj = grad_exact_project(Xh, val);
            return (gproj - gradv(T));
        }
        /****************/

        // may add optionnal weight for AXI???
        double L2_error(const space_ptrtype & Xh, const element_type & T) const
        {
            // replace params by specified values
            ex solution = getSolution();

            auto g = expr(solution,vars);
            auto mesh = Xh->mesh();
            return  math::sqrt( integrate( elements(mesh), Px()*(idv(T) - g)*(idv(T) - g) ).evaluate()(0,0) );
        }

        double H1_error(const space_ptrtype & Xh, const element_type & T) const
        {
            // replace params by specified values
            ex solution = getSolution();

            auto gradg = expr<1,Dim,2>(grad(solution,vars), vars );
            auto mesh = Xh->mesh();

            double L2error = L2_error(Xh, T);
            double H1seminorm = math::sqrt( integrate( elements(mesh), Px()*(gradv(T) - gradg)*trans(gradv(T) - gradg) ).evaluate()(0,0) );
            return math::sqrt( L2error*L2error + H1seminorm*H1seminorm);
        }

        double L2_error(const space_ptrtype & Xh, const element_type & T, const double & values) const
        {
            // replace params by specified values
            ex solution = getSolution(values);

            auto g = expr(solution,vars);
            auto mesh = Xh->mesh();
            return  math::sqrt( integrate( elements(mesh), Px()*(idv(T) - g)*(idv(T) - g) ).evaluate()(0,0) );
        }

        double H1_error(const space_ptrtype & Xh, const element_type & T, const double & values) const
        {
            // replace params by specified values
            ex solution = getSolution(values);

            auto gradg = expr<1,Dim,2>(grad(solution,vars), vars );
            auto mesh = Xh->mesh();

            double L2error = L2_error(Xh, T, values);
            double H1seminorm = math::sqrt( integrate( elements(mesh), Px()*(gradv(T) - gradg)*trans(gradv(T) - gradg) ).evaluate()(0,0) );
            return math::sqrt( L2error*L2error + H1seminorm*H1seminorm);
        }


    protected:

        //! name of the exact solution
        std::string M_exact;

        //! name of the exact solution parameters
        std::string M_exact_params;

        //! name of the rhs
        std::string M_rhs;
        bool M_rhs_computed;

        // convergence study
        bool M_convergence;
        int M_convergence_max;

        ex exact;
        ex rhs;
        std::vector<symbol> vars;
        std::vector<symbol> parameters;


    };


}
#endif /* __ERRROR_HPP */
