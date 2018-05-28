/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2012-01-19

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
/**
 \file modelalgebraicfactory.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2012-01-19
 */

#ifndef FEELPP_MODELSALGEBRAICFACTORY_HPP
#define FEELPP_MODELSALGEBRAICFACTORY_HPP 1

#include <feel/feelmodels/modelcore/modelalgebraic.hpp>


namespace Feel
{
namespace FeelModels
{

    class ModelAlgebraicFactory
    {
    public :
        typedef ModelAlgebraicFactory self_type;

        typedef FeelModels::ModelAlgebraic model_type;
        typedef boost::shared_ptr<model_type> model_ptrtype;
        typedef boost::weak_ptr<model_type> model_weakptrtype;

        typedef typename model_type::value_type value_type;
        typedef typename model_type::backend_type backend_type;
        typedef typename model_type::backend_ptrtype backend_ptrtype;

        typedef GraphCSR graph_type;
        typedef boost::shared_ptr<graph_type> graph_ptrtype;

        typedef typename model_type::vector_ptrtype vector_ptrtype;
        typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
        typedef typename sparse_matrix_ptrtype::element_type sparse_matrix_type;

        typedef Preconditioner<value_type> preconditioner_type;
        typedef boost::shared_ptr<preconditioner_type> preconditioner_ptrtype;

        typedef typename backend_type::solve_return_type solve_return_type;
        typedef typename backend_type::nl_solve_return_type nl_solve_return_type;

        typedef model_type::indexsplit_type indexsplit_type;
        typedef model_type::indexsplit_ptrtype indexsplit_ptrtype;


        typedef boost::function<void ( ModelAlgebraic::DataUpdateLinear& )> function_assembly_linear_type;
        typedef boost::function<void ( ModelAlgebraic::DataUpdateJacobian& )> function_assembly_jacobian_type;
        typedef boost::function<void ( ModelAlgebraic::DataUpdateResidual& )> function_assembly_residual_type;
        typedef boost::function<void ( vector_ptrtype& )> function_newton_initial_guess_type;

        typedef typename backend_type::pre_solve_type pre_solve_type;
        typedef typename backend_type::post_solve_type post_solve_type;

        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        ModelAlgebraicFactory(model_ptrtype const& __app, backend_ptrtype const& __backend);

        ModelAlgebraicFactory(model_ptrtype const&__app,
                              backend_ptrtype const& __backend,
                              graph_ptrtype const& graph,
                              indexsplit_ptrtype const& indexSplit );

        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        template <typename SpaceType>
        void
        initFromFunctionSpace(boost::shared_ptr<SpaceType> const& space )
        {

            if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".MethodNum","initFromFunctionSpace", "start",
                                                                this->model()->worldComm(),this->model()->verboseAllProc());

            this->model()->timerTool("Constructor").start();
            auto graph=stencil(_test=space,_trial=space,
                               _pattern_block=this->model()->blockPattern(),
                               _diag_is_nonzero=true)->graph();
            auto indexSplit= space->dofIndexSplit();
            this->model()->timerTool("Constructor").elapsed("graph");

            this->model()->timerTool("Constructor").restart();
            this->buildMatrixVector(graph,indexSplit);
            this->model()->timerTool("Constructor").elapsed("matrixVector");

            this->model()->timerTool("Constructor").restart();
            this->buildOthers();
            this->model()->timerTool("Constructor").elapsed("algebraicOthers");

            if (this->model()->verbose()) Feel::FeelModels::Log(this->model()->prefix()+".MethodNum","initFromFunctionSpace", "finish",
                                                                this->model()->worldComm(),this->model()->verboseAllProc());
        }

        //---------------------------------------------------------------------------------------------------------------//

        void
        rebuildMatrixVector( graph_ptrtype const& graph,
                             indexsplit_ptrtype const& indexSplit);

        void reset(backend_ptrtype __backend,
                   graph_ptrtype const& graph,
                   indexsplit_ptrtype const& indexSplit);

        boost::shared_ptr<std::ostringstream> getInfo() const;
        void printInfo() const;

        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        //! return the model used by the algebraic factory
        model_ptrtype
        model() const { return M_model.lock(); }

        //! return the backend
        backend_ptrtype const&
        backend() const { return M_backend; }

        //! return matrix used for assembly linear system or the jacobian
        sparse_matrix_ptrtype matrix() const { return M_J; }

        //! return vector used for assembly rhs in linear system or the residual
        vector_ptrtype rhs() const { return M_R; }

        //! return the preconditioner
        preconditioner_ptrtype const&
        preconditionerTool() const { return M_PrecondManage; }

        //! return sparsity matrix graph
        graph_ptrtype const&
        sparsityMatrixGraph() const { CHECK(M_J) << "no matrix register"; return M_J->graph(); }

        void
        attachNullSpace( NullSpace<value_type> const& nullSpace );
        void
        attachNearNullSpace( NullSpace<value_type> const& nearNullSpace );
        void
        attachNearNullSpace( int k, NullSpace<value_type> const& nearNullSpace );


        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        void solve( std::string const& type,vector_ptrtype& sol );
        FEELPP_DEPRECATED void linearSolver( vector_ptrtype &U ) { this->solveLinear( U ); }
        FEELPP_DEPRECATED void AlgoNewton2( vector_ptrtype &U ) { this->solveNewton( U ); }
        FEELPP_DEPRECATED void AlgoPicard(vector_ptrtype &U) { this->solvePicard( U ); }
        void solveLinear( vector_ptrtype &U );
        void solveNewton( vector_ptrtype &U );
        void solvePicard( vector_ptrtype &U );

        void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J );
        void updateResidual( const vector_ptrtype& X, vector_ptrtype& R);

        void rebuildCstJacobian( vector_ptrtype U );
        void rebuildCstLinearPDE( vector_ptrtype U );
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//

        void addFunctionLinearAssembly( function_assembly_linear_type const& func, std::string const& key = "" );
        void addFunctionLinearPostAssembly( function_assembly_linear_type const& func, std::string const& key = "" );
        void addFunctionNewtonInitialGuess( function_newton_initial_guess_type const& func, std::string const& key = "" );
        void addFunctionJacobianAssembly( function_assembly_jacobian_type const& func, std::string const& key = "" );
        void addFunctionResidualAssembly( function_assembly_residual_type const& func, std::string const& key = "" );
        void addFunctionJacobianPostAssembly( function_assembly_jacobian_type const& func, std::string const& key = "" );
        void addFunctionResidualPostAssembly( function_assembly_residual_type const& func, std::string const& key = "" );

    private :

        void
        init(graph_ptrtype const& graph,
             indexsplit_ptrtype const& indexSplit);

        void
        buildMatrixVector(graph_ptrtype const& graph,
                          indexsplit_ptrtype const& indexSplit);
        void
        buildOthers();

    private :

        model_weakptrtype M_model;

        backend_ptrtype M_backend;

        preconditioner_ptrtype M_PrecondManage;

        vector_ptrtype M_R;
        vector_ptrtype M_CstR;
        sparse_matrix_ptrtype M_J;
        sparse_matrix_ptrtype M_CstJ;
        sparse_matrix_ptrtype M_Prec;
        sparse_matrix_ptrtype M_Extended;

        bool M_hasBuildLinearJacobian;
        bool M_hasBuildResidualCst;
        bool M_hasBuildLinearSystemCst;

        std::map<std::string,function_assembly_linear_type> M_addFunctionLinearAssembly;
        std::map<std::string,function_assembly_linear_type> M_addFunctionLinearPostAssembly;
        std::map<std::string,function_newton_initial_guess_type> M_addFunctionNewtonInitialGuess;
        std::map<std::string,function_assembly_jacobian_type> M_addFunctionJacobianAssembly;
        std::map<std::string,function_assembly_residual_type> M_addFunctionResidualAssembly;
        std::map<std::string,function_assembly_jacobian_type> M_addFunctionJacobianPostAssembly;
        std::map<std::string,function_assembly_residual_type> M_addFunctionResidualPostAssembly;
    };


} // namespace FeelModels
} // namespace Feel


#endif //FEELPP_MODELSALGEBRAICFACTORY_HPP
