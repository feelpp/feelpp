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
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2012-06-02
*/

#ifndef __ERROR_HPP
#define __ERROR_HPP 1

#include <string>
#include <sstream>

#include <boost/algorithm/string/split.hpp>

#include <boost/parameter.hpp>
#include <feel/feelcore/parameter.hpp>

#include <feel/feelvf/inner.hpp>

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
    using  GiNaC::matrix;

    BOOST_PARAMETER_NAME(sol)
    BOOST_PARAMETER_NAME(val)
    BOOST_PARAMETER_NAME(axi)
    BOOST_PARAMETER_NAME(component)

    /**
     * \class ErrorBase
     *
     * Provides abstract interface to error calculations
     *
     * \tparam Dim dimension of space
     */
        template<int Dim, int Order, template<uint16_type PsetDim> class PsetType>
        class ErrorBase
        {
        public:
            //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
            typedef Simplex<Dim> convex_type;
            //! mesh type
            typedef Mesh<convex_type> mesh_type;
            typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

            //! the basis type of our approximation space
            //typedef bases<Lagrange<Order,Scalar> > basis_type;
            typedef bases<Lagrange<Order, PsetType> > basis_type;

            //! the approximation function space type
            typedef FunctionSpace<mesh_type, basis_type > space_type;
            typedef boost::shared_ptr<space_type> space_ptrtype;
            typedef typename space_type::element_type element_type;

            static const uint16_type nComp = PsetType<Dim>::nComponents;
            typedef typename mpl::if_<mpl::bool_<space_type::is_scalar>, ex, matrix>::type ex_type;

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
                M_convergence_max(vm[prefixvm( prefix, "error.convergence.steps" )].as<int>() ),
                vars( symbols<Dim>() )
            {
                //resize
                if( nComp != 1 ) //boost::is_same<ex_type, matrix>::value )
                    {
                        exact = matrix(1,nComp);
                        grad_exact = matrix(nComp,Dim);
                        curl_exact = matrix(nComp, Dim); // not valid in 2D
                        div_exact = matrix(nComp, Dim); //???
                        rhs = matrix(1,nComp);
                    }
                else
                    {
                        grad_exact = matrix(1, Dim);
                    }
            }

            ErrorBase( std::string const& name, std::string const& params)
                :
                M_exact( name ),
                M_exact_params(params),
                M_rhs(""),
                M_rhs_computed( false ),
                M_convergence(0),
                M_convergence_max(0),
                vars( symbols<Dim>() )
            {
                //resize
                if( nComp != 1 ) //boost::is_same<ex_type, matrix>::value )
                    {
                        exact = matrix(1,nComp);
                        grad_exact = matrix(nComp, Dim);
                        curl_exact = matrix(nComp, Dim); // not valid in 2D
                        div_exact = matrix(nComp, Dim); //???
                        rhs = matrix(1,nComp);
                    }
                else
                    {
                        grad_exact = matrix(1, Dim);
                    }
            }

            ErrorBase( ErrorBase const &e )
            :
            M_exact( e.M_exact ),
            M_exact_params( e.M_exact_params ),
            M_rhs( e.M_rhs ),
            M_rhs_computed( e.M_rhs_computed ),
            M_convergence( e.M_convergence ),
            M_convergence_max( e.M_convergence_max ),
            exact_file( e.exact_file ),
            grad_exact_file( e.grad_exact_file ),
            curl_exact_file( e.curl_exact_file ),
            div_exact_file( e.div_exact_file ),
            rhs_file( e.rhs_file ),
            vars( e.vars )
            {
                //resize
                if( nComp != 1 ) //boost::is_same<ex_type, matrix>::value )
                    {
                        exact = matrix(1,nComp);
                        grad_exact = matrix(nComp,Dim);
                        curl_exact = matrix(nComp, Dim); // not valid in 2D
                        div_exact = matrix(nComp, Dim); //???
                        rhs = matrix(1,nComp);
                    }
                else
                    {
                        grad_exact = matrix(1, Dim);
                    }
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

                        exact_file = e.exact_file;
                        grad_exact_file = e.grad_exact_file;
                        curl_exact_file = e.curl_exact_file;
                        div_exact_file = e.div_exact_file;
                        rhs_file = e.rhs_file;
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
                space_vars = vars;

                std::vector<std::string> lst_params;
                boost::split(lst_params, M_exact_params, boost::is_any_of(",;"));
                LOG(INFO) << "Loading additionnal parameters : ";
                for(const std::string &s: lst_params){LOG(INFO) << s << ", ";};
                LOG(INFO) << std::endl;
                for(const std::string &s: lst_params){vars.push_back(symbol(s));};

                LOG(INFO) << "Loading function : " << M_exact << "\n";
                for(symbol const& s : parameters) {LOG(INFO) << "Symbol " << s.get_name() << " found\n";};
                LOG(INFO) << "exact solution is : " << exact << "\n";

                string_to_ex(M_exact, exact);

                // Define grad_exact
                grad_exact = grad(exact, space_vars);

                // Set names by default for Ginac shared libs
                exact_file = "_exact";
                grad_exact_file = "_grad_exact";
            }

            void setRhs(std::string const& expression ="")
            {
                if ( expression.empty() )
                    FEELPP_ASSERT( M_rhs.size() ).error( "undefined rhs exact" );
                else
                    M_rhs = expression;

                LOG(INFO) << "Loading rhs function : " << M_rhs << std::endl;

                string_to_ex(M_rhs, rhs);
                rhs_file = "_rhs";
            }

            typedef typename boost::function<ex_type (ex_type & u, std::vector<symbol> & vars, std::vector<symbol> & space_vars)> edp_type;

            void setRhs(edp_type * edp)
            {
                rhs = (*edp)(exact, vars, space_vars);
                LOG(INFO) << "computed rhs is : " << rhs << "\n";
                std::cout << "rhs=" << rhs << "\n" << std::flush;

                std::ostringstream rhs_expression;
                rhs_expression << rhs;
                M_rhs = rhs_expression.str();

                rhs_file = "_computedrhs";
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
             * \return GiNaC solution file
             */
            std::string getSolutionFile() const
            {
                return exact_file;
            }

            /**
             * \return GiNaC solution file
             */
            std::string getGradSolutionFile() const
            {
                return grad_exact_file;
            }

            /**
             * \return GiNaC solution file
             */
            std::string getCurlSolutionFile() const
            {
                return curl_exact_file;
            }

            /**
             * \return GiNaC solution file
             */
            std::string getDivSolutionFile() const
            {
                return div_exact_file;
            }

            /**
             * \return GiNaC rhs file
             */
            std::string getRhsFile() const
            {
                return rhs_file;
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
            ex_type getRhs() const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_rhs.size() ).error( "undefined rhs" );
                return rhs;
            }

            /**
             * \return GiNaC exact solution
             */
            ex_type getSolution() const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_exact.size() ).error( "undefined exact solution" );
                return exact;
            }

            /**
             * \return GiNaC grad exact solution
             */
            ex_type getGradSolution() const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_grad_exact.size() ).error( "undefined grad(exact)" );
                return grad_exact;
            }
            /**
             * \return GiNaC curl exact solution
             */
            ex_type getCurlSolution() const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_curl_exact.size() ).error( "undefined curl(exact)" );
                return curl_exact;
            }
            /**
             * \return GiNaC exact solution
             */
            ex_type getDivSolution() const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_div_exact.size() ).error( "undefined div(exact)" );
                return div_exact;
            }

            /**
             * \return GiNaC exact gradient solution
             */
            ex_type getGradSolution(const int & i) const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_exact.size() ).error( "undefined exact solution" );

                ex_type grad_i;
                if ( nComp == 1)
                    grad_i = grad_exact(0,i);
                else
                    {
                        lst l;
                        for(int n=0; n<Dim; n++)
                            l.append(grad_exact(i,n));
                        grad_i = matrix(1,nComp, l);
                    }
                return grad_i;
            }

            /**
             * \return GiNaC exact curl solution
             */
            ex_type getCurlSolution(const int & i) const
            {
                // check if exact is defined
                FEELPP_ASSERT( M_exact.size() ).error( "undefined exact solution" );

                ex_type curl_i;
                if ( nComp == 1) // ?? what to do in 2D
                    curl_i = curl_exact(0,i);
                else
                    {
                        lst l;
                        for(int n=0; n<Dim; n++)
                            l.append(curl_exact(i,n));
                        curl_i = matrix(1,nComp, l);
                    }
                return curl_i;
            }

            //@}

            void print() const
            {
                LOG(INFO) << "============================================================\n";
                LOG(INFO) << "Error Information\n";
                LOG(INFO) << "   exact : " << getSolution() << "\n";
                LOG(INFO) << "   grad(exact) : " << grad_exact << "\n";
                LOG(INFO) << "\n";
                if ( !M_rhs.empty() )
                    {
                        LOG(INFO) << "   rhs : " << getRhs() << "\n";
                    }
                LOG(INFO) << "   convergence : " << convergence() << "\n";
                LOG(INFO) << "   convergence steps : " << numberOfConvergenceSteps() << "\n";
            }

            /* The call of projection depend on ex_type (matrix or ex) */
            // Projection of a Vectorial element
            template<typename T>
            typename boost::enable_if< boost::is_same<T, matrix>, Expr< GinacMatrix<nComp,1,2> > >::type
            ex_to_project(T& expression, const std::string & filename ="")
            {
                return expr<nComp,1,2>(expression, vars, filename);
            }
            // Projection of a Scalar element
            template<typename T>
            typename boost::disable_if< boost::is_same<T, matrix>, Expr< GinacEx<2> > >::type
            ex_to_project(T& expression, const std::string & filename ="")
            {
                return expr(expression, vars, filename);
            }

            /* Transformation from string into ginac expression depend on ex_type (matrix or ex) */
            // Ex of a Vectorial element (GiNac::matrix)
            template<typename T>
            typename boost::enable_if< boost::is_same<T, matrix>, void>::type
            string_to_ex(std::string exname, T& expression)
            {
                std::vector<std::string> lst_comp;
                boost::split(lst_comp, exname, boost::is_any_of(";,"));
                for(int i=0; i<lst_comp.size(); i++)
                    expression.set(0, i, parse(lst_comp[i], vars));
            }
            // Ex of a Scalar element (GiNac::ex)
            template<typename T>
            typename boost::disable_if< boost::is_same<T, matrix>, void>::type
            string_to_ex(std::string exname, T& expression)
            {
                expression = parse(exname, vars);
            }

            typedef std::map<std::string,double> val_type;
            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( element_type ), // return type
                                            rhs_project,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( space, * ))
                                            (optional
                                             ( range, *, elements(space->mesh()) )
                                             ( val, *, val_type() ))
                                            )
            {
                auto rhs_expr = ex_to_project( rhs, rhs_file);
                rhs_expr.expression().setParameterValues(val);
                return vf::project( _space=space, _range=range, _expr=rhs_expr);
            }

            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( element_type ), // return type
                                            exact_project,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( space, * ))
                                            (optional
                                             ( range, *, elements(space->mesh()))
                                             ( val, *, val_type() ) )
                                            )
            {
                auto solution_expr = ex_to_project( exact, exact_file);
                solution_expr.expression().setParameterValues(val);
                return vf::project( _space=space, _range=range, _expr=solution_expr);
            }

            // wrong return type
            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( element_type ), // return type
                                            grad_exact_project,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( space, * )
                                             ( component, *))
                                            (optional
                                             ( range, *, elements(space->mesh()))
                                             ( val, *, val_type() ))
                                            )
            {
                FEELPP_ASSERT( component >= 0 && component < Dim ).error( "wrong component index" );

                std::string tmpfile = grad_exact_file;
                std::ostringstream cmp; cmp << "_" << component;
                tmpfile += cmp.str();

                // actual return type depend on exact
                // should distinct is exact is ex or matrix
                // if exact is matrix, gradg[i] is matrix
                // ........... ex, grad[i] is ex

                ex_type grad_i = getGradSolution(component);
                auto gradg = ex_to_project(grad_i, tmpfile);
                gradg.expression().setParameterValues(val);

                return vf::project( _space=space, _range=range, _expr=gradg );
            }

            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( element_type ), // return type
                                            error_project,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( space, * )
                                             ( sol, (typename ErrorBase<Dim, Order, PsetType>::element_type) ) )
                                            (optional
                                             ( range, *, elements(space->mesh()))
                                             ( val, *, val_type() ))
                                            )
            {
                auto gproj = exact_project(_space=space, _range=range, _val=val);
                return (gproj - sol);
            }

            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( element_type ), // return type
                                            grad_error_project,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( space, * )
                                             ( sol, (typename ErrorBase<Dim, Order, PsetType>::element_type) )
                                             ( component, *) )
                                            (optional
                                             ( range, *, elements(space->mesh()))
                                             ( val, *, val_type() ) )
                                            )
            {
                FEELPP_ASSERT( component >= 0 && component < Dim ).error( "wrong component index" );

                auto gradg = grad_exact_project(_space=space, _range=range, _component=component, _val=val);

                switch (component)
                    {
                    case 0:
                        {
                            auto gsol = vf::project( _space=space, _range=range, _expr=dxv(sol));
                            return (gradg - gsol);
                        }
                    case 1:
                        {
                            auto gsol = vf::project( _space=space, _range=range, _expr=dyv(sol));
                            return (gradg - gsol);
                        }
                    case 2:
                        {
                            auto gsol = vf::project( _space=space, _range=range, _expr=dzv(sol));
                            return (gradg - gsol);
                        }
                    default:
                        {
                            element_type e(space);
                            return (e);
                        }
                    }

            }

            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( double ), // return type
                                            L2_error,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( space, * )
                                             ( sol, (typename ErrorBase<Dim, Order, PsetType>::element_type) ) )
                                            (optional
                                             ( range, *, elements(space->mesh()) )
                                             ( axi, (bool) , false)
                                             ( val, *, val_type()  ))
                                            )
            {
                auto g = ex_to_project(exact, exact_file);
                g.expression().setParameterValues(val);

                double L2;
                if( axi )
                    L2 = integrate( _range=range, _expr=Px()*inner( idv(sol) - g, idv(sol) - g ) ).evaluate()(0,0);
                else
                    L2 = integrate( _range=range, _expr=inner( idv(sol) - g, idv(sol) - g ) ).evaluate()(0,0);

                return math::sqrt(L2);
            }

            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( double ), // return type
                                            H1_seminorm,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( sol, * )
                                             ( space, * ))
                                            (optional
                                             ( range, *, elements(space->mesh()) )
                                             ( axi, (bool), false)
                                             ( val, *, val_type() ))
                                            )
            {
                //auto gradg = ex_to_project(grad_exact, grad_exact_file);
                auto gradg = expr<nComp,Dim,2>(grad_exact, vars, grad_exact_file);
                gradg.expression().setParameterValues(val);

                double H1seminorm;
                if( axi )
                    H1seminorm = integrate( _range=range, _expr=Px()*(gradv(sol) - gradg)*trans(gradv(sol) - gradg) ).evaluate()(0,0);
                else
                    H1seminorm = integrate( _range=range, _expr=(gradv(sol) - gradg)*trans(gradv(sol) - gradg) ).evaluate()(0,0);

                return math::sqrt(H1seminorm);

            }

            BOOST_PARAMETER_MEMBER_FUNCTION(
                                            ( double ), // return type
                                            H1_error,    // 2. function name
                                            tag,           // 3. namespace of tag types
                                            (required
                                             ( sol, * )
                                             ( space, * ))
                                            (optional
                                             ( range, *, elements(space->mesh()) )
                                             ( axi, (bool), false)
                                             ( val, *, val_type() ))
                                            )
            {
                double L2 = L2_error(_range=range, _space=space, _sol=sol, _axi=axi, _val=val);
                double H1seminorm = H1_seminorm(_range=range, _space=space, _sol=sol, _axi=axi, _val=val);

                return math::sqrt( L2*L2 + H1seminorm*H1seminorm );

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

            ex_type exact;
            matrix grad_exact;
            matrix curl_exact;
            matrix div_exact;
            ex_type rhs;

            std::string exact_file;
            std::string grad_exact_file;
            std::string curl_exact_file;
            std::string div_exact_file;
            std::string rhs_file;

            std::vector<symbol> vars;
            std::vector<symbol> space_vars;
            std::vector<symbol> parameters;


        };

}
#endif /* __ERRROR_HPP */
